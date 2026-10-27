"""Program S for the plugin: the specification as data, a sample-by-sample reference implementation, and reference
IRs to check a port against (docs/model_s_results.md §17, §38.1; docs/NEXT_STEPS.md §10a).

    .venv/bin/python tools/s_spec.py

Writes analysis/s/s_plugin_spec.json (every delay, coefficient, table and tap of the adopted baseline, each with its
label and source section) and analysis/s/ref/s_<LF>_<RT>_<HF>[_amp<N>].wav (device rate 26,373.6 Hz, float32, (L, R),
the input impulse at sample 0, tap 0 at sample 0, direct-tap units: the AMPL-90 cluster tap = 1).

The model is the adopted baseline: EMT_KERNEL=chain3, EMT_PMAX=0.8475, EMT_HFRT_MIN=0.12, the AP230 table
ap230_mltable.json, the loop_cv fold set hfcells_chain3_cap0.8475_apml and pulse_check's default compact fit. The
resolved loop filters are checked against q2.final_G's PB / R_i for every setting, and simulate() is checked against
q1.Joint.render with the measurement band applied (the band is library-side and not part of the spec).

    simulate(n, setting) -> the unit at device rate, (n, 2): direct taps + [tap injection -> four loops -> AP20 -> AP72
    -> AP230(damping b)] -> sum. No analysis band (q1.LFBAND), no capture-chain kernel, no sub-100 Hz capture transfer,
    no output kernel: those belong to the 96 kHz library, not to the unit. The plugin's own analogue chain (converter
    filters, sample-and-hold, 20 Hz high-pass, L half a sample late) is separate and follows this spec.
"""
import copy
import json
import os
import sys

os.environ["EMT_KERNEL"] = "chain3"
os.environ["EMT_PMAX"] = "0.8475"
os.environ["EMT_HFRT_MIN"] = "0.12"
os.environ.pop("EMT_AP230", None)
os.environ.pop("EMT_COMPACT", None)
os.environ.pop("PULSE_CV", None)
os.environ.pop("EMT_DCFLOOR", None)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import soundfile as sf  # noqa: E402
from scipy.signal import lfilter, sosfreqz  # noqa: E402

import loop_fir as lf_  # noqa: E402
import model_s as ms  # noqa: E402
import pulse_check as pc  # noqa: E402
import q1_generalise as q1  # noqa: E402
import q2_injection as q2  # noqa: E402

FS = q1.FS
P = [int(v) for v in q1.P]
S0 = q1.S0
PRE = q1.PRE
OUT = os.path.join(ROOT, "analysis", "s")
REF = os.path.join(OUT, "ref")
AP_G = -0.70710678
REF_SETTINGS = [("0.1_0.3_0.2", 90), ("0.4_0.6_0.4", 90), ("1.0_1.0_1.0", 90), ("1.0_0.8_0.2", 90),
                ("1.0_1.8_1.0", 90), ("2.0_1.4_1.0", 90), ("2.0_2.5_0.6", 90), ("0.1_2.5_0.2", 90),
                ("1.0_1.0_1.0", 30)]

INFO = None
LAW0 = None
FOLDS = None
CELLS = None
AP230 = None
V = None
TAPS = None
ALPHA = None


def _load():
    global INFO, LAW0, FOLDS, CELLS, AP230, V, TAPS, ALPHA
    INFO, _, _, _, _ = q2.q1_tables()
    LAW0 = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    AP230 = json.load(open(os.path.join(q2.OUT, lf_.AP230)))
    CELLS = json.load(open(os.path.join(pc.CV, "cells.json")))
    FOLDS = {rt: json.load(open(os.path.join(pc.CV, f"fold_{rt}.json"))) for rt in sorted({g["rt"] for g in INFO.values()})}
    V, _, ALPHA = q2.final_injection(compact=pc.COMPACT)
    TAPS = q2.tap_table()


def _lab(key):
    g = INFO[key]
    return dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])


def _b230(key):
    return float(AP230["steps"][AP230["settings"][key]["step"]])


def shelf_b(dc_gain, nyquist_gain, pole):
    """First-order shelf (b0 + b1 z^-1) / (1 - pole z^-1), DC and Nyquist gains as model_s.shelf."""
    b0 = (dc_gain * (1 - pole) + nyquist_gain * (1 + pole)) / 2
    b1 = (dc_gain * (1 - pole) - nyquist_gain * (1 + pole)) / 2
    return b0, b1, -pole


def loop_filter(lab, b230, fold, cell):
    """The per-loop filter G_i of the adopted baseline, as (scalar gain, sections), from the same law and tables that
    pulse_check.setting_G resolves: the G_capped / final_Gs form at the baseline switches (loopk = bbtilt = tilt = 0)."""
    law = copy.deepcopy(LAW0)
    law["lf"] = {float(k): v for k, v in law["lf"].items()}
    for lf_k, dc in fold.get("lfdc", {}).items():
        law["lf"][float(lf_k)]["dc"] = float(dc)
    gamma, fl, hfdp = fold["gamma"], {}, fold["hfdp"]
    dc_min, beta, gamma_dc = q2.DC_MIN, q2.FINAL_BETA, True
    if cell and "dc_min" in cell:
        law["lf"][0.1]["dc"] = cell["dc"]
        gamma, dc_min, beta, gamma_dc = cell["gamma"], cell["dc_min"], cell["beta"], cell.get("gamma_dc", True)
    elif cell and "floor_form" in cell:
        fl = dict(floor_form=cell["floor_form"])
    if cell and "hf_poles" in cell:
        hfdp = {"poles": cell["hf_poles"]}
    loopcell = dict(gain=cell["lf_gain"], dc=cell["lf_dc"]) if cell and "lf_gain" in cell else None
    rt, lf, hf = lab["rt"], lab["lf"], lab["hf"]
    L, pramp, lfrt = law["lf"][lf], tuple(fold["pramp"]), fold["lfrt"]
    out = []
    for i, Pi in enumerate(P):
        kk = -3 * Pi / FS
        g, sections, floored = 10 ** (kk / (L["hi"] * rt)), [], False
        if "floor_form" in fl and L["dc"] < L["hi"]:
            ff = fl["floor_form"]
            r0 = 10 ** (kk * (1 / (L["dc"] * rt) - 1 / (L["hi"] * rt)))
            x, fd, c = 20 * np.log10(g * r0), 20 * np.log10(ff["F"]), ff["c"]
            dd = max(x, fd - c) + c * np.exp(-max(x - fd + c, 0.0) / ff["tau"])
            u = 20 * np.log10(ff["F0"]) - x
            if u > 0:
                g *= 10 ** ((-ff["beta"] * u + ff["gamma"] * 60 / (L["hi"] * rt) * Pi / FS) / 20)
            sections.append(dict(kind="shelf", dc_gain=10 ** (dd / 20) / g, nyquist_gain=1.0, pole=law["a"]))
        elif L["dc"] != L["hi"]:
            r0 = 10 ** (kk * (1 / (L["dc"] * rt) - 1 / (L["hi"] * rt)))
            if g * r0 < dc_min:
                floored = True
                lift = 20 * np.log10(dc_min / (g * r0))
                g *= 10 ** (-beta * lift / 20)
            if floored and gamma and not gamma_dc:
                g *= 10 ** (gamma * 60 / (L["hi"] * rt) * Pi / (20 * FS))
            sections.append(dict(kind="shelf", dc_gain=max(r0, dc_min / g), nyquist_gain=1.0, pole=law["a"]))
            if floored and gamma and gamma_dc:
                sections.append(dict(kind="gain", db=gamma * 60 / (L["hi"] * rt) * Pi / FS))
        if hf < 1.0:
            hf_e = max(hf, q2.HFRT_MIN / rt) if q2.HFRT_MIN else hf
            p_law = float(q1.onepole_for(np.array(10 ** (kk / rt * L["k"] * (1 / hf_e - 1))), law["f_ref"]))
            p = min(p_law + pramp[0] * max(0.0, p_law - pramp[1]), q2.P_MAX) if pramp else min(p_law, q2.P_MAX)
            if hfdp and "poles" in hfdp:
                p = hfdp["poles"][i]
            elif (hfdp and str(hf) in hfdp and hfdp.get("p_above", 0.0) <= p < hfdp.get("p_below", 1.0)
                  and lf not in hfdp.get("exclude_lf", [])):
                p += hfdp[str(hf)] + hfdp.get("slope", 0.0) * (Pi - np.mean(P)) / np.mean(P)
            sections.append(dict(kind="onepole", pole=float(p)))
        if str(lf) in lfrt["coefs"]:
            a_, b_ = lfrt["coefs"][str(lf)]
            rel = a_ + b_ * rt
            if str(lf) in lfrt.get("clip", []) and rel * a_ <= 0:
                rel = 0.0
            g *= 10 ** (rel * (60 / (L["hi"] * rt)) * Pi / (20 * FS))
        if loopcell:
            g *= 10 ** (loopcell["gain"][i] / 20)
            sections.append(dict(kind="shelf", dc_gain=10 ** (loopcell["dc"][i] / 20), nyquist_gain=1.0, pole=law["a"]))
        out.append(dict(P=Pi, gain=float(g), sections=sections))
    return out


def loop_response(loop, w):
    z = np.exp(-1j * np.asarray(w))
    h = loop["gain"] * np.ones(len(z), complex)
    for s in loop["sections"]:
        if s["kind"] == "gain":
            h = h * 10 ** (s["db"] / 20)
        elif s["kind"] == "shelf":
            b0, b1, a1 = shelf_b(s["dc_gain"], s["nyquist_gain"], s["pole"])
            h = h * (b0 + b1 * z) / (1 + a1 * z)
        else:
            h = h * (1 - s["pole"]) / (1 - s["pole"] * z)
    return h


def loop_poly(loop):
    """G_i as a direct-form polynomial (b length m+1, a length n+1, a[0] = 1)."""
    b, a = np.array([1.0]), np.array([1.0])
    for s in loop["sections"]:
        if s["kind"] == "gain":
            b = b * 10 ** (s["db"] / 20)
        elif s["kind"] == "shelf":
            b0, b1, a1 = shelf_b(s["dc_gain"], s["nyquist_gain"], s["pole"])
            b, a = np.convolve(b, [b0, b1]), np.convolve(a, [1.0, a1])
        else:
            b, a = np.convolve(b, [1.0 - s["pole"]]), np.convolve(a, [1.0, -s["pole"]])
    return loop["gain"] * b / a[0], a / a[0]


def resolved(key):
    lab = _lab(key)
    loops = loop_filter(lab, _b230(key), FOLDS[lab["rt"]], CELLS.get(key))
    for loop in loops:
        b, a = loop_poly(loop)
        loop["b"], loop["a"] = [float(v) for v in b], [float(v) for v in a]
        w1 = 2 * np.pi * 1000 / FS
        loop["at_1khz"] = [float(v) for v in (np.abs(loop_response(loop, [w1])[0]), np.angle(loop_response(loop, [w1])[0]))]
        loop["dc_gain"] = float(np.real(loop_response(loop, [0.0])[0]))
    return lab, _b230(key), loops


def simulate(n, loops, cluster_gains, b230, ampl=90):
    """The reference implementation: the unit at device rate, (n, 2), direct-tap units, the input impulse at sample 0.
    loops: the resolved per-loop filters; cluster_gains: (2, 4) in the model's taps order; ampl: the AMPL step."""
    y = np.zeros((n, 2, len(P)))
    for i, Pi in enumerate(P):
        b, a = np.asarray(loops[i]["b"]), np.asarray(loops[i]["a"])
        for c in (0, 1):
            u = np.zeros(n)
            for (loop, ch, delay, sign, _) in TAPS:
                if loop == i and ch == c:
                    u[delay + S0] += sign * ALPHA
            zi, yy, vv = np.zeros(len(a) - 1), y[:, c, i], np.zeros(n)
            for start in range(0, n, Pi):
                stop = min(n, start + Pi)
                w, zi = lfilter(b, a, vv[start:stop], zi=zi)
                yy[start:stop] = u[start:stop] + w
                if stop + Pi <= n:
                    vv[start + Pi:stop + Pi] = yy[start:stop]
                else:
                    vv[start + Pi:n] = yy[start:n - Pi]
    s = y.sum(axis=2)

    def ap(x, M):
        bb = np.zeros(M + 1); bb[0], bb[M] = -AP_G, 1.0
        aa = np.zeros(M + 1); aa[0], aa[M] = 1.0, -AP_G
        return lfilter(bb, aa, x, axis=0)

    z = ap(ap(s, 20), 72)
    f0 = 1 - b230
    bb = np.zeros(232); bb[0], bb[230], bb[231] = -AP_G, f0, b230
    aa = np.zeros(232); aa[0], aa[230], aa[231] = 1.0, -AP_G * f0, -AP_G * b230
    z = lfilter(bb, aa, z, axis=0)
    for c in (0, 1):
        for t, gain in zip(ms.CLUSTER[c], cluster_gains[c]):
            if t < n:
                z[t, c] += (ampl / 90) * gain
    return z


def _check(key, loops, b230):
    """Resolved G_i against the model's PB / R_i, over 25 Hz-11 kHz (PB has zeros at the band edges)."""
    G = pc.setting_G(_lab(key), b230, copy.deepcopy(LAW0), FOLDS[_lab(key)["rt"]], CELLS.get(key))
    w = 2 * np.pi * np.fft.rfftfreq(G["N"])
    f = w * FS / (2 * np.pi)
    m = (f > 25) & (f < 11000)
    worst = 0.0
    for i, Pi in enumerate(P):
        true = np.zeros_like(G["PB"])
        true[m] = (1 - G["PB"][m] / G["R"][i][m]) * np.exp(1j * w[m] * Pi)
        got = loop_response(loops[i], w)
        worst = max(worst, float(np.max(np.abs(got[m] - true[m])) / np.max(np.abs(true[m]))))
    return G, worst


def _check_sim(key, loops, b230, G, n):
    """simulate() against Joint.render with the measurement band (q1.LFBAND) applied the way the model applies it."""
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((G["N"], 2)), w=np.ones(G["N"]), amp=90)])])
    model = J.render(V, G, G["N"], 90)
    gains = np.array([100.0 * V[:4], 100.0 * V[4:8]])
    sim = simulate(n, loops, gains, b230)
    w = 2 * np.pi * np.fft.rfftfreq(n)
    band = np.abs(sosfreqz(q1.LFBAND, worN=w)[1]) ** 2
    plain = np.zeros((n, 2))
    for c in (0, 1):
        for t, gain in zip(ms.CLUSTER[c], gains[c]):
            if t < n:
                plain[t, c] += gain
    reverb = np.fft.irfft(np.fft.rfft(sim - plain, axis=0) * band[:, None], n, axis=0)
    check = reverb + J.cluster(V[:8])[PRE:PRE + n]
    ref = model[PRE:PRE + n]
    return {f"{a}-{b}": [float(10 * np.log10(np.sum((check[int(a / 1000 * FS):int(b / 1000 * FS), c]
                                                    - ref[int(a / 1000 * FS):int(b / 1000 * FS), c]) ** 2)
                                           / np.sum(ref[int(a / 1000 * FS):int(b / 1000 * FS), c] ** 2)))
                         for c in (0, 1)]
            for a, b in ((0, 50), (50, 150), (150, 300), (300, 600), (600, 1000))}


def spec():
    lab = lambda v, c, src: dict(value=v, label=c, source=src)
    gains = np.array([100.0 * V[:4], 100.0 * V[4:8]])
    settings, worst_filter, worst_sim = {}, 0.0, -999.0
    for key in sorted(INFO):
        lab_, b230, loops = resolved(key)
        _, err = _check(key, loops, b230)
        worst_filter = max(worst_filter, err)
        settings[key] = dict(lf=lab_["lf"], rt=lab_["rt"], hf=lab_["hf"], ap230_b=b230, loops=loops)
    refs = {}
    for key, ampl in REF_SETTINGS:
        b230 = _b230(key)
        loops = [dict(loop) for loop in settings[key]["loops"]]
        G = pc.setting_G(_lab(key), b230, copy.deepcopy(LAW0), FOLDS[_lab(key)["rt"]], CELLS.get(key))
        n = G["N"] - PRE
        nulls = _check_sim(key, loops, b230, G, n)
        worst_sim = max([worst_sim] + [v for cc in nulls.values() for v in cc])
        name = f"s_{key}" + (f"_amp{ampl}" if ampl != 90 else "")
        refs[name] = dict(file=f"analysis/s/ref/{name}.wav", setting=key, ampl=ampl, null_vs_model_db=nulls)
    law = LAW0
    return {
        "program": "S (EMT 246 program S), docs/model_s_results.md §17, §38.1",
        "device_rate_hz": FS,
        "units": "direct-tap units: the AMPL-90 cluster tap = 1; delays in integer device samples from tap 0; mono in, "
                 "(L, R) out; the input impulse at sample 0, tap 0 at sample 0",
        "signal_path": "input -> tap injection -> four parallel loops (filter G_i in the feedback, delay P_i) -> sum -> "
                       "AP(20) -> AP(72) -> AP(230, damping b) -> output; the direct path (four taps per channel) adds in "
                       "parallel. The reverb level does not move with AMPL; the direct path is x AMPL/90",
        "measurement_side": "not in this spec and not in the reference IRs: the 20 Hz-11 kHz analysis band (q1.LFBAND), "
                            "the chain3 capture-chain kernel, the shared sub-100 Hz capture transfer "
                            "(tools/lowband_joint.py), and the deconvolution's band-limit window on the direct taps. The "
                            "plugin's own analogue chain (converter filters, sample-and-hold, 20 Hz high-pass, L half a "
                            "sample late) is separate and follows this spec",
        "baseline": {"kernel": "chain3", "p_max": q2.P_MAX, "hrt_min": q2.HFRT_MIN, "dc_min": q2.DC_MIN,
                     "beta": q2.FINAL_BETA, "fold_set": os.path.relpath(pc.CV, ROOT),
                     "compact_fit": pc.COMPACT, "ap230_table": lf_.AP230, "ap230_g": AP230["g230"]},
        "direct_path": lab(dict(taps={"L": list(ms.CLUSTER[0]), "R": list(ms.CLUSTER[1])},
                                gains={"L": gains[0].tolist(), "R": gains[1].tolist()},
                                ampl="direct taps x AMPL / 90; the reverb does not move with AMPL (S's linear cluster law)",
                                onset="tap 0 at sample 0; the model's band-limited atom is measurement-side"),
                           "S", "§17, §21 (chain2 refit: L = R within 0.1 %), §38.1"),
        "injection": lab(dict(taps=[dict(loop=i, channel="LR"[c], delay=int(t) + S0, sign=int(sg), tap=j) for i, c, t, sg, j in TAPS],
                              magnitude=float(ALPHA), format="a sum of signed impulses of magnitude `magnitude` at `delay` samples "
                                                             "from tap 0, one per loop and channel, added to the loop input",
                              note="fitted, held out by RT; the same shape in every setting, scaled by one magnitude", ),
                         "S (fitted, held out) / H (the physical injection)", "§16-17, §22, §38.1"),
        "loops": lab(dict(P=P, form="y_i[n] = x_i[n] + G_i * (z^-P_i y_i)[n]; G_i is the loop filter in the feedback, "
                                   "stable, 0-3 first-order sections (shelf, gain, one-pole)",
                          loop_filter="settings.<key>.loops[i]: gain, sections (exact form) and b / a (direct form, a[0] = 1); "
                                      "shelf(x) = (b0 + b1 z^-1) / (1 - pole z^-1) with b0, b1 from (dc_gain, nyquist_gain, pole) "
                                      "as model_s.shelf; onepole = (1 - pole) / (1 - pole z^-1)"),
                     "S", "§17, §22, §32"),
        "diffuser": lab(dict(allpasses=[dict(M=20, g=AP_G), dict(M=72, g=AP_G)],
                             ap230=dict(M=230, g=float(AP230["g230"]), b_steps=[float(v) for v in AP230["steps"]],
                                        form="H = (-g + z^-M F) / (1 - g z^-M F), F = (1 - b) + b z^-1",
                                        b_per_setting="settings.<key>.ap230_b"),
                             form="all-passes in series, g = -1/sqrt(2); order immaterial (LTI)"),
                        "S", "§10, §17, §25.14, §38.1"),
        "loop_filter_law": lab(dict(law={str(k): v for k, v in law["lf"].items()}, a=law["a"], f_ref=law["f_ref"],
                                    g="10^(-3 P / (fs hi RT))", shelf_dc="10^(-3 P / fs (1/(dc RT) - 1/(hi RT)))",
                                    hf_pole="one-pole at f_ref, T = 10^(-3 P / (fs RT) k (1/HF - 1)), p = min(p_law + c max(0, p_law - p0), P_MAX)",
                                    dc_floor="g r0 >= DC_MIN; the cells' dc_min / beta / gamma and the folds' lfdc / hfdp / pramp / lfrt",
                                    cells="analysis/loop_cv/hfcells_chain3_cap0.8475_apml/cells.json (dc, dc_min, beta, gamma, "
                                          "gamma_dc, hf_poles, lf_gain, lf_dc, or floor_form)",
                                    folds="analysis/loop_cv/hfcells_chain3_cap0.8475_apml/fold_<RT>.json (gamma, lfdc, hfdp, "
                                          "pramp, lfrt)"),
                               "S at the library's steps / H between them (nearest RT fold, nearest sampled cell; the "
                               "cells do not interpolate, §32.8)", "§22, §25, §29, §32, §34-35, §38.1"),
        "settings": settings,
        "controls": {
            "RT": lab(dict(values_s=sorted({g["rt"] for g in INFO.values()})), "C (library steps) / H (between)", "§1"),
            "LF": lab(dict(values_s=sorted({g["lf"] for g in INFO.values()})), "C (library steps) / H (between)", "§1"),
            "HF": lab(dict(values_s=sorted({g["hf"] for g in INFO.values()})), "C (library steps) / H (between)", "§1"),
            "AMPL": lab("direct taps x AMPL / 90 (AMPL 90 = the reference); the reverb is fixed", "S", "§17, §21"),
            "DELAY": lab("a pure delay of the whole output; the library is all DELAY 0", "H", "operation card"),
            "LOW_PASS": lab("the panel's LOW PASS has no varying files in the library; response unknown", "H", "§7"),
        },
        "not_resolved": ["the LOW PASS response and its steps (no library data)", "DELAY steps (always 0 in the library)",
                         "the law between the library's LF / RT / HF steps (H: nearest fold and sampled cell)",
                         "the direct taps' and injection's physical cause (fitted; S's AMPL is a linear cluster law)"],
        "references": dict(files=refs,
                           note="device rate, float32, (L, R), input impulse at sample 0, tap 0 at sample 0, direct-tap "
                                "units. A port must null these sample for sample (float precision). They are the unit "
                                "structure without the measurement band; nulling against the library model instead needs "
                                "the band applied, where simulate() is at -71 dB or better (2.0_1.4_1.0: -63.8 dB, the "
                                "model's 2.48 s FFT grid wraps its tail)"),
        "validation": dict(loop_filter_max_rel_error=worst_filter,
                           simulate_vs_model_worst_null_db=worst_sim,
                           simulate_note="time-domain simulate() against q1.Joint.render with the measurement band "
                                         "applied to both, reference settings, 0-1000 ms; the band is library-side"),
    }


def main():
    _load()
    os.makedirs(REF, exist_ok=True)
    sp = spec()
    for name, r in sp["references"]["files"].items():
        key, ampl = r["setting"], r["ampl"]
        b230 = _b230(key)
        loops = [dict(loop) for loop in sp["settings"][key]["loops"]]
        n = pc.setting_G(_lab(key), b230, copy.deepcopy(LAW0), FOLDS[_lab(key)["rt"]], CELLS.get(key))["N"] - PRE
        gains = np.array([100.0 * V[:4], 100.0 * V[4:8]])
        y = simulate(n, loops, gains, b230, ampl=ampl).astype(np.float32)
        sf.write(os.path.join(ROOT, r["file"]), y, int(round(FS)), subtype="FLOAT")
        print(f"{name}: {y.shape} at {int(round(FS))} Hz")
    json.dump(sp, open(os.path.join(OUT, "s_plugin_spec.json"), "w"), indent=1)
    print(f"loop filters vs model: {sp['validation']['loop_filter_max_rel_error']:.2e} relative")
    print(f"simulate vs model: worst {sp['validation']['simulate_vs_model_worst_null_db']:.1f} dB (0-1000 ms)")
    print("wrote analysis/s/s_plugin_spec.json")


if __name__ == "__main__":
    main()
