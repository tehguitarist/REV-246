"""Programs M and L for the plugin: the specification as data, a sample-by-sample reference implementation, and
reference IRs to check a port against (docs/ml_structure.md §7-9, §14-26; docs/NEXT_STEPS.md §10a).

    .venv/bin/python tools/ml_spec.py

Writes analysis/ml/ml_plugin_spec.json (every delay, coefficient, table and tap of the adopted baseline, each with its
label and source section) and analysis/ml/ref/{m,l}_<LF>_<RT>_<HF>.wav (device rate 26,373.6 Hz, float32, (L, R), the
input impulse at sample 0, the output kernel applied causally: cluster tap 0 at sample k_h = 10 (M) / 30 (L); direct-tap
units).

The model is the adopted baseline: ML_KFIX refit kernels (kernel_fix_full.json), ML_BREFIT=steps (one b table for M and
L), the HF pole curve (ML_POLE_CURVE=32), HFRT_MIN=0.12 and L's RT 2.2 floor 0.40, the chain3 library. The resolved
line filters are checked against the model's line_coefs, and simulate() is checked against ml_lines.Model.parts(band=False)
(the 100 Hz-10 kHz analysis band is library-side and not part of the spec). The time variance (theta) is described but not
in the static references: it is always on in the plugin and is checked statistically (ml_structure.md §9).

    simulate(prog, key, n) -> the static network at device rate, (n, 2): direct taps + [tap injection -> four coupled
    lines (filter G_j in the feedback, delay P_j, coupling F) -> sum -> diffuser] -> output kernel. No analysis band, no
    capture clock / shift / gains: those belong to the 96 kHz library and its captures.
"""
import copy
import json
import os
import sys

for _v in ("ML_KFIX", "ML_B0"):
    os.environ.pop(_v, None)
for _v in ("ML_FP", "ML_HK", "ML_KFIX_FILE", "ML_BSTEP", "ML_HFRT", "ML_L22"):
    os.environ.pop(_v, None)
os.environ["ML_BREFIT"] = "steps"
os.environ["ML_POLE_CURVE"] = "32"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "harness"))

import numpy as np  # noqa: E402
import soundfile as sf  # noqa: E402
from scipy.signal import lfilter  # noqa: E402

import ml_carryover as cc  # noqa: E402
import ml_lines as mll  # noqa: E402
import ml_structure as ml  # noqa: E402

FS, K0 = mll.FS, mll.K0
AP_G = -1 / np.sqrt(2)
OUT = os.path.join(ROOT, "analysis", "ml")
REF = os.path.join(OUT, "ref")
PROGS = ("M", "L")
REF_SETTINGS = {
    "M": ["0.1_0.6_0.1", "0.7_1.2_0.4", "1.2_1.8_0.9", "2.0_3.0_0.2", "2.0_5.0_0.9"],
    "L": ["0.1_2.2_0.1", "0.8_4.5_0.5", "1.2_6.0_0.8", "2.0_10.0_0.3", "2.0_20.0_1.0"],
}


def line_poly(coef):
    """G_j as a causal filter: the S shelf (DC gain r0, Nyquist 1, pole LAW['a']) then the HF one-pole p."""
    g, r0, p = coef
    a = mll.LAW["a"]
    b0, b1 = (r0 * (1 - a) + (1 + a)) / 2, (r0 * (1 - a) - (1 + a)) / 2
    b, aa = np.array([b0, b1], float), np.array([1.0, -a])
    if p:
        b, aa = np.convolve(b, [1 - p]), np.convolve(aa, [1.0, -p])
    return g * b / aa[0], aa / aa[0]


def simulate(M, key, n):
    """The static network at device rate, (n, 2), direct-tap units, the input impulse at sample 0, the output kernel
    causal (cluster tap 0 at sample k_h)."""
    lf, rt, hf = (float(x) for x in key.split("_"))
    coefs = mll.line_coefs(M.P, rt, lf, hf)
    t = np.zeros((n, 2, len(M.P)))
    for ch in (0, 1):
        for p, sg, i in zip(np.round(M.pos[ch]), M.sgn[ch], M.loop[ch]):
            if int(p) < n:
                t[int(p), ch, int(i)] += sg
    x = np.zeros((n, 2, len(M.P)))
    filt = [line_poly(c) for c in coefs]
    zi = [[np.zeros(len(filt[j][1]) - 1) for _ in (0, 1)] for j in range(len(M.P))]
    block = min(M.P)
    for s0 in range(0, n, block):
        e = min(n, s0 + block)
        vs = []
        for j, Pj in enumerate(M.P):
            v = np.zeros((e - s0, 2))
            for k in range(len(M.P)):
                idx = np.arange(s0, e) - Pj
                ok = idx >= 0
                v[ok] += M.F[j, k] * x[idx[ok], :, k]
            vs.append(v)
        for j in range(len(M.P)):
            b, a = filt[j]
            for ch in (0, 1):
                out, zi[j][ch] = lfilter(b, a, vs[j][:, ch], zi=zi[j][ch])
                x[s0:e, ch, j] = t[s0:e, ch, j] + out
    rev = x.sum(axis=2)
    b = M.b.get(key, 0.0)
    for k, Md in enumerate(cc.DIFF[M.prog]):
        bb = b if k == len(cc.DIFF[M.prog]) - 1 else 0.0
        f0 = 1 - bb
        bn = np.zeros(Md + 2); bn[0], bn[Md], bn[Md + 1] = -AP_G, f0, bb
        an = np.zeros(Md + 2); an[0], an[Md], an[Md + 1] = 1.0, -AP_G * f0, -AP_G * bb
        rev = lfilter(bn, an, rev, axis=0)
    cl = np.zeros((n, 2))
    for ch in (0, 1):
        for p, sg in M.cl[ch]:
            if int(p) < n:
                cl[int(p), ch] += sg
    s = cl + rev
    return np.stack([np.convolve(s[:, ch], M.h)[:n] for ch in (0, 1)], axis=1)


def check(M, key, n):
    """simulate() against ml_lines.Model.parts(band=False), aligned for the causal output kernel."""
    extra = K0 - M.kh
    parts = M.parts(n + extra, key, band=False)
    ref = (parts[0] + parts[1])[K0 - M.kh:K0 - M.kh + n]
    got = simulate(M, key, n)
    d = got - ref
    return {f"{a}-{b}": [float(10 * np.log10(np.sum(d[int(a / 1000 * FS):int(b / 1000 * FS), c] ** 2)
                                             / np.sum(ref[int(a / 1000 * FS):int(b / 1000 * FS), c] ** 2)))
                         for c in (0, 1)]
            for a, b in ((0, 50), (50, 150), (150, 300), (300, 600), (600, 1000))}


def spec(M):
    prog = M.prog
    lab = lambda v, c, src: dict(value=v, label=c, source=src)
    settings = {}
    for key in sorted(M.b):
        lf, rt, hf = (float(x) for x in key.split("_"))
        coefs = mll.line_coefs(M.P, rt, lf, hf)
        settings[key] = dict(lf=lf, rt=rt, hf=hf, b=float(M.b[key]),
                             lines=[dict(P=int(Pi), gain=float(g), shelf_dc=float(r0), hf_pole=float(p))
                                    for Pi, (g, r0, p) in zip(M.P, coefs)])
    refs, worst = {}, -999.0
    for key in REF_SETTINGS[prog]:
        lf, rt, hf = (float(x) for x in key.split("_"))
        n = min(int(3 * rt + 0.5), 12) * int(round(FS))
        nulls = check(M, key, n)
        worst = max([worst] + [v for cc_ in nulls.values() for v in cc_])
        refs[key] = dict(file=f"analysis/ml/ref/{prog.lower()}_{key}.wav", n=n, null_vs_model_db=nulls)
    return {
        "program": f"{prog} ('{prog}' on the 246), docs/ml_structure.md §7-9, §14-26",
        "device_rate_hz": FS,
        "units": "direct-tap units (unit cluster gain); delays in integer device samples from tap 0; mono in, (L, R) out; "
                 "the input impulse at sample 0",
        "signal_path": "static/reference path: input -> tap injection -> four coupled lines (filter G_j in the feedback, "
                       "delay P_j, coupling F) -> sum -> diffuser (all-passes, the last damped by b) -> output kernel; the "
                       "direct path (four taps per channel) adds in parallel. The plugin's time-varying path is the "
                       "read-out form in modulation.application: the diffused input feeds the lines and each channel's "
                       "taps read u_i; with theta = 0 the two forms are equivalent (ml_variance.core, about -174 dB)",
        "measurement_side": "not in this spec and not in the reference IRs: the 100 Hz-10 kHz analysis band (ml.BAND), "
                            "the capture clock ratio and shift (ml_lines.warp), the per-capture cluster and reverb gains, "
                            "the file fades and usable ends, and the chain3 capture-chain kernel. The plugin's own "
                            "analogue chain and time variance are separate and follow this spec",
        "baseline": {"kernel": "chain3", "output_kernel": "analysis/ml_lines/kernel_fix_full.json (ML_KFIX)",
                     "damping": "ML_BREFIT=steps, analysis/ml_lines/b_steps_ml.json (one table for M and L)",
                     "pole_curve": "ML_POLE_CURVE=32, analysis/lowband/chain3_pole_cells_allf.json",
                     "hrt_min": mll.HFRT_MIN, "l22_hf": mll.L22_HF, "p_max": mll.P_MAX, "dc_min": mll.DC_MIN, "beta": mll.BETA},
        "lines": lab(dict(P=list(M.P), form="x_j[n] = t_j[n] + G_j * (sum_k F[j,k] x_k)[n - P_j]; the four filters are in "
                                          "settings.<key>.lines (gain, shelf DC gain, HF pole), the shelf pole is the "
                                          "S law's a = 0.916883 and the shelf's Nyquist gain is 1 (the HF one-pole's is "
                                          "(1 - p) / (1 + p))",
                           filter_form="G_j = gain * (b0 + b1 z^-1)/(1 - a z^-1) * (1 - p)/(1 - p z^-1), "
                                       "b0 = (r0 (1 - a) + (1 + a))/2, b1 = (r0 (1 - a) - (1 + a))/2, r0 = shelf_dc; "
                                       "p = hf_pole and p = 0 bypasses the HF one-pole",
                           law="S's line law at the line's own length: g = 10^(-3 P / (fs hi RT)), shelf DC = 10^(-3 P / fs "
                               "(1/(dc RT) - 1/(hi RT))), HF one-pole at f_ref with T = 10^(-3 P / (fs RT) k (1/HF_e - 1)), "
                               "HF_e = max(HF, HFRT_MIN x program scale / RT); the DC floor lifts g by 10^(-beta lift / 20) "
                               "and sets the shelf to dc_min / g when g r0 < DC_MIN; the p_max cap is used only when the "
                               "pole curve is off (ML_POLE_CURVE=0) or a non-default p_max is passed (the adopted curve "
                               "maps and caps at its own knots / DROP_BACK). L's RT 2.2 floor 0.40 and HFRT_MIN 0.12 x the "
                               "program scale (x2 M, x4 L)",
                           law_tables={str(k): v for k, v in mll.LAW["lf"].items()}, a=mll.LAW["a"], f_ref=mll.LAW["f_ref"],
                           pole_curve=dict(knots=json.load(open(os.path.join(ROOT, "analysis", "lowband",
                                                                              "chain3_pole_cells_allf.json")))["knots"],
                                           values=json.load(open(os.path.join(ROOT, "analysis", "lowband",
                                                                             "chain3_pole_cells_allf.json")))["all"],
                                           fold=json.load(open(os.path.join(ROOT, "analysis", "lowband",
                                                                           "chain3_pole_cells_allf.json")))["fold"],
                                           drop_back=mll.DROP_BACK)),
                       "S at the library's steps / H between them (LF log-interpolated between S's steps; RT and HF "
                       "continuous; b from the nearest (RT, HF) cell)", "§8, §8.3, §8.7, §16, §25.20"),
        "coupling": lab(dict(F=M.F.tolist(), form="y = h * [C + D * 1^T (I - Lambda F)^-1 T], Lambda = diag(G_j z^-P_j); "
                                                  "F is orthogonal (the bipartite matrix, entries +-1/sqrt(2))"),
                        "S", "§7.3, §9"),
        "direct_path": lab(dict(taps={c: [[int(p), int(sg)] for p, sg in M.cl[i]] for i, c in enumerate("LR")},
                                gains="unit; AMPL scales the cluster only (the pairs measure that the reverb tails are "
                                      "identical after one gain); the absolute AMPL scale is the plugin's calibration (H)",
                                onset="tap 0 at sample k_h after the output kernel (causal), not at sample 0"),
                           "S (taps) / H (AMPL scale)", "§2, §7.5, §9"),
        "injection": lab(dict(taps=[dict(channel="LR"[ch], delay=int(p), sign=int(sg), line=int(i))
                                    for ch in (0, 1) for p, sg, i in zip(np.round(M.pos[ch]), M.sgn[ch], M.loop[ch])],
                              note="the first-pass/read taps, one sum per channel and line; fitted on 6 training "
                                   "settings, held out by setting and AMPL"),
                         "S", "§7.2-7.4"),
        "diffuser": lab(dict(delays=list(cc.DIFF[prog]), g=float(cc.G_AP),
                             damping="the last all-pass's two-tap loop filter (1 - b) + b z^-1, b per setting in "
                                     "settings.<key>.b; the step comes from the setting's (RT, HF) cell in "
                                     "analysis/ml_lines/b_steps.json, and the value from the shared table in "
                                     "analysis/ml_lines/b_steps_ml.json",
                             form="all-passes H = (-g + z^-M F)/(1 - g z^-M F) in series"),
                        "S", "§7.4, §19, §22"),
        "output_kernel": lab(dict(taps=[float(v) for v in M.h], centre=int(M.kh),
                                  note="causal: the kernel's centre tap k_h lands at the cluster tap 0's time, so the "
                                       "reference IRs place cluster tap 0 at sample k_h"),
                             "S", "§7.5, §18, §25"),
        "modulation": lab(dict(rate_hz=0.83, amplitude_rad=0.18, rms_rad=0.13,
                               measured_rms="0.13 rad (0.11-0.144, bootstrap over groups), against Tomkin's AMPL pairs with 12 evenly spaced sim phases (analysis/plugin_check/theta_depth.json, ml_structure.md 9.8); the earlier 0.15 (9.7) came from 6 random phases per group that happened to understate the decorrelation",
                               form="theta(t) = amplitude sin(2 pi rate t + phase); phase random per capture; a sine "
                                    "shape is H, the depth is S",
                               coupling="one angle of both 2x2 blocks: F(theta) = F R(theta) per block, R = [[cos, -sin], "
                                        "[sin, cos]]; F stays orthogonal",
                               application="read-out form: u_i = x + sum_j F(theta)_ji v_j, v_j = G_j u_j delayed P_j; "
                                           "each channel's taps read u_i. With theta = 0 it reproduces the static network "
                                           "to about -100 dB (ml_variance.py)",
                               status="always on in the plugin; the static reference IRs are theta = 0. Check the moving "
                                      "part statistically: Tomkin's AMPL pairs decorrelate like theta renders from the "
                                      "first recirculation to 2-4 s (ml_structure.md §9.5-9.7)"),
                         "S (rate, depth, one angle) / H (sine shape and phase)", "§9.5-9.7"),
        "settings": settings,
        "controls": {
            "RT": lab(dict(values_s=sorted({rt for _, rt, _ in (tuple(map(float, k.split("_"))) for k in M.b)})),
                      "C (library steps) / H (between)", "§1"),
            "LF": lab(dict(values_s=sorted({lf for lf, _, _ in (tuple(map(float, k.split("_"))) for k in M.b)})),
                      "C (library steps) / H (between)", "§1"),
            "HF": lab(dict(values_s=sorted({hf for _, _, hf in (tuple(map(float, k.split("_"))) for k in M.b)})),
                      "C (library steps) / H (between)", "§1"),
            "AMPL": lab(f"the cluster only; the reverb is fixed (the AMPL pairs are repeat captures after one gain). The "
                        f"library steps: {sorted({r['amp'] for r in json.load(open(os.path.join(cc.OUT, 'firstpass.json')))[prog]['rows']})}",
                        "S (cluster-only) / H (scale)", "§2"),
            "DELAY": lab("a pure delay of the whole output; the library is all DELAY 0", "H", "operation card"),
            "LOW_PASS": lab("the panel's LOW PASS has no varying files in the library; response unknown", "H", "§7"),
        },
        "not_resolved": ["theta's shape and phase (the sine is H; the 0.13 rad rms depth is S)",
                         "the LOW PASS response and DELAY steps (no library data)",
                         "the law between the library's LF / RT / HF steps (H: LF log-interpolated, RT / HF continuous; "
                         "b nearest (RT, HF) cell)",
                         "the absolute AMPL scale (per-capture gains are fitted, so the model does not fix it)"],
        "references": dict(files=refs, wav_sample_rate_hz=int(round(FS)),
                           note="device rate (the exact device_rate_hz; the WAV header is the integer "
                                "wav_sample_rate_hz), float32, (L, R), input impulse at sample 0, direct-tap units, the "
                                "output kernel causal so cluster tap 0 is at sample k_h. A port must null these sample for "
                                "sample (float precision). They are the static (theta = 0) network without the measurement "
                                "band, where simulate() is at -130 dB or better against ml_lines.Model.parts"),
        "validation": dict(simulate_vs_model_worst_null_db=worst,
                           simulate_note="time-domain simulate() against ml_lines.Model.parts(band=False) at the reference "
                                         "settings, 0-1000 ms, aligned for the causal output kernel"),
    }


def main():
    os.makedirs(REF, exist_ok=True)
    out, worst = {}, -999.0
    for prog in PROGS:
        M = mll.Model(prog)
        sp = spec(M)
        out[prog] = sp
        for key, r in sp["references"]["files"].items():
            y = simulate(M, key, r["n"]).astype(np.float32)
            sf.write(os.path.join(ROOT, r["file"]), y, int(round(FS)), subtype="FLOAT")
            print(f"{prog} {key}: {y.shape} at {int(round(FS))} Hz")
        worst = max(worst, sp["validation"]["simulate_vs_model_worst_null_db"])
        print(f"{prog}: simulate vs model worst {sp['validation']['simulate_vs_model_worst_null_db']:.1f} dB")
    out["validation"] = dict(simulate_vs_model_worst_null_db=worst,
                             simulate_note="time-domain simulate() against ml_lines.Model.parts(band=False) at the "
                                           "reference settings, 0-1000 ms, aligned for the causal output kernel")
    json.dump(out, open(os.path.join(OUT, "ml_plugin_spec.json"), "w"), indent=1)
    print("wrote analysis/ml/ml_plugin_spec.json")


if __name__ == "__main__":
    main()
