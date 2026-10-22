"""Altiverb 246 S presets against the consolidated S model (model_s_results.md §17), held out, and the differences
from the CLA and Tomkin captures. Results: model_s_results.md §18, docs/EMT246_unit_facts.md §5c.

    .venv/bin/python tools/altiverb_predict.py [--stage all|predict|direct|diag] [--out analysis/altiverb]

Needs analysis/altiverb/cache/ (tools/altiverb_compare.py --stage deconv) and analysis/q1, analysis/q2, analysis/shortRT.
Writes <out>/predict.json.

  predict   The §17 S model (trained on the Tomkin library, nothing refitted) at LF 1.0, HF 0.9, AMPL 80 and RT from the
            preset name (settings inferred in altiverb_compare.py --stage settings; HF 0.9 is not a Tomkin step, as on CLA).
            Per preset: clock ratio, one common shift and one gain fitted on the Altiverb IR (3 numbers), as for CLA in §17,
            then a 7-tap FIR per channel fitted on the cluster-only region (14 numbers; without it the null is -12/-15 dB because
            the Altiverb's R is 0.74 x CLA's and its L has a different chain shape). Nulls 100 Hz-10 kHz per window and channel.
  direct    The Altiverb RT 1.0 IR against the CLA IR (the same panel setting: LF 1.0, RT 1.0, HF 0.9, AMPL 80): clock,
            shift and gain, then with the per-channel FIR. Shows how far the two units and captures are apart without a model.
  diag      Polarity, L/R offset, device-domain clock, and the first-arrival structure of the Altiverb, CLA and Tomkin S files.

Floor: the sweep deconvolution agrees with an exact frequency-domain deconvolution to -74 dB (100 Hz-10 kHz, first 2 s),
and two renders at different levels agree to -104 dB (altiverb_compare.py), so the data are good to about -74 dB.
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import q1_generalise as q1  # noqa: E402
import q2_injection as q2  # noqa: E402
import analyze_ir as ai  # noqa: E402

FS, PRE = q1.FS, q1.PRE            # device rate, samples before tap 0
FS48 = 48000
WIN = q2.W_COMP[:6]                # 0-50 ... 1000-1500 ms
FLOOR_DB = -74.0
LAB = dict(lf=1.0, hf=0.9, amp=80)
PRESETS = {"0.6": 0.6, "1.0": 1.0, "1.6": 1.6, "2.0": 2.0, "2.2": 2.2}
E0 = {}
ONLY = os.environ.get("ALT_ONLY")  # smoke test: one preset, e.g. ALT_ONLY=1.0


def device_ir(x48, fs=FS48):
    """Deconvolved capture (2 x N or N x 2) -> device-rate IR with tap 0 at index PRE; returns (h, info)."""
    x = x48 if x48.shape[1] == 2 else x48.T
    h, t0, offs = ai.recover(x, fs, ai.FS_DEV)
    cl = ai.cluster_taps(h)
    k0 = cl[0]["first_tap_index"]
    sign = float(np.sign(h[k0, :].sum()))
    return h[k0 - PRE: k0 + int(1.9 * FS)], dict(k0=int(k0), channel_offsets_us=[float(o * 1e6) for o in offs],
                                                  tap0_sign=sign, tap0_sign_LR=[float(np.sign(h[k0, 0])), float(np.sign(h[k0, 1]))])


def altiverb_h(rt):
    ir = np.load(os.path.join(ROOT, "analysis", "altiverb", "cache", f"NON-LIN SHORT RT {rt}s.npy")).T
    fa = int(np.argmax(np.abs(ir[:, 0]) > np.abs(ir[:, 0]).max() * 10 ** (-40 / 20)))
    h, info = device_ir(ir[fa - 240: fa + int(3.5 * FS48)])
    return h, info


def cla_h():
    return np.load(q1.CLA_DEV)[__import__("model_s").K0 - PRE: __import__("model_s").K0 + int(1.9 * FS)]


def fit_clock(yb, mdl, starts=range(-300, 301, 100)):
    """Clock, shift and gain fitted against the model; the starting clock is scanned (the Altiverb's is not known)."""
    best = None
    for e0 in starts:
        _, (eg, shg, gain) = q1._clock_job((yb, mdl, 1.0, float(e0), True))
        y = q1.apply_clock(yb, eg, shg)[:len(mdl)] / gain
        r = q1.nulls(y, mdl, 1.5, [(50, 300)], narrow=True)["50-300"]
        if best is None or np.mean(r) < best[0]:
            best = (float(np.mean(r)), eg, shg, gain, y)
    return best[1:4]


def model_for(rt, n, amp, v, law, b230=0.0):
    w = 2 * np.pi * np.fft.rfftfreq(q1.grid_n(rt))
    lab = dict(LAB, rt=rt, amp=amp)
    G = q2.make_G(lab, law, q2.G_capped(lab, law, w, beta=q2.FINAL_BETA), b230, q2.FINAL_DIFF)
    return q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=amp)])]).render(v, G, n, amp)


def channel_fir(y, mdl, K=3):
    """Per-channel (2K+1)-tap FIR taking the model to the data, fitted on the cluster-only region (before each channel's
    reverb starts: L 225, R 555 device samples after tap 0), so no reverb data are used. Counted: 2 x (2K+1) numbers."""
    ends = {0: PRE + 225, 1: PRE + 555}
    out, hs = np.zeros_like(mdl), []
    for c in (0, 1):
        A = np.stack([np.roll(mdl[:, c], k) for k in range(-K, K + 1)], 1)
        hc, *_ = np.linalg.lstsq(A[PRE - 40: ends[c]], y[PRE - 40: ends[c], c], rcond=None)
        out[:, c] = A @ hc
        hs.append(hc)
    return out, hs


def stage_predict(out):
    q1.EXT = 1700
    v, law, a = q2.final_injection()
    res = {}
    print("S model (9 numbers + library laws), nothing fitted on the Altiverb except as stated")
    print("window (ms): " + " ".join(f"{a_}-{b_}".rjust(11) for a_, b_ in WIN) + "   (L / R, dB re the Altiverb IR, 100 Hz-10 kHz)")
    for name, rt in PRESETS.items():
        if ONLY and name != ONLY:
            continue
        h, info = altiverb_h(name)
        yb = sosfiltfilt(q1.LFBAND, h * (info["tap0_sign"] or 1.0), axis=0)
        n = len(yb) - 600
        r = {"info": info}
        for amp in (80,):  # AMPL cannot be scanned here: the cluster FIR absorbs the cluster gain (AMPL comes from the tap/reverb ratio)
            mdl = model_for(rt, n, amp, v, law)
            eg, shg, gain = fit_clock(yb, mdl, starts=[E0["eps"]] if "eps" in E0 else range(-300, 301, 100))
            E0["eps"] = eg
            y = q1.apply_clock(yb, eg, shg)[:n] / gain
            r[f"amp{amp}_plain"] = dict(eps_ppm=eg, shift=shg, gain_db=float(20 * np.log10(gain)), null=q1.nulls(y, mdl, 1.5, WIN, narrow=True))
            mdl2, hs = channel_fir(y, mdl)
            eg2, shg2, gain2 = fit_clock(yb, mdl2, starts=[eg])
            y2 = q1.apply_clock(yb, eg2, shg2)[:n] / gain2
            r[f"amp{amp}_fir"] = dict(eps_ppm=eg2, shift=shg2, gain_db=float(20 * np.log10(gain2)), fir=[np.round(x, 4).tolist() for x in hs],
                                      null=q1.nulls(y2, mdl2, 1.5, WIN, narrow=True))
        res[name] = r
        for tag in ("plain", "fir"):
            a80 = r[f"amp80_{tag}"]
            print(f"RT {name} {tag:5s}: clock {a80['eps_ppm']:+.1f} ppm shift {a80['shift']:+.3f} | "
                  + " ".join((f"{x[0]:.1f}/{x[1]:.1f}" if x else "--").rjust(11) for x in a80["null"].values()))
        print("   channel FIR (L, R), 7 taps each:", r["amp80_fir"]["fir"])
    return res


def stage_direct(out):
    ha, _ = altiverb_h("1.0")
    hc = cla_h()
    n = min(len(ha), len(hc)) - 600
    ya = sosfiltfilt(q1.LFBAND, ha * np.sign(ha[PRE].sum()), axis=0)[:n + 600]
    yc = sosfiltfilt(q1.LFBAND, hc * np.sign(hc[PRE].sum()), axis=0)[:n + 600]
    res = {}
    for who, (data, ref) in {"Altiverb vs CLA": (ya, yc[:n]), "CLA vs Altiverb": (yc, ya[:n])}.items():
        eg, shg, gain = fit_clock(data, ref)
        y = q1.apply_clock(data, eg, shg)[:n] / gain
        res[who] = dict(eps_ppm=eg, shift=shg, gain_db=float(20 * np.log10(gain)), null=q1.nulls(y, ref, 1.5, WIN, narrow=True))
        print(f"{who} (gain only): clock {eg:+.1f} ppm shift {shg:+.3f} gain {res[who]['gain_db']:+.2f} dB | "
              + " ".join((f"{x[0]:.1f}/{x[1]:.1f}" if x else "--").rjust(11) for x in res[who]["null"].values()))
        ref2, hs = channel_fir(y, ref)       # per-channel 7-tap FIR on the cluster-only region, as in predict
        eg2, shg2, gain2 = fit_clock(data, ref2, starts=[eg])
        y2 = q1.apply_clock(data, eg2, shg2)[:n] / gain2
        res[who + " + channel FIR"] = dict(eps_ppm=eg2, shift=shg2, null=q1.nulls(y2, ref2, 1.5, WIN, narrow=True))
        print(f"{who} + channel FIR (14): clock {eg2:+.1f} ppm shift {shg2:+.3f} | "
              + " ".join((f"{x[0]:.1f}/{x[1]:.1f}" if x else "--").rjust(11) for x in res[who + ' + channel FIR']["null"].values()))
    return res


def stage_diag(out):
    import soundfile as sf
    from scipy.signal import resample_poly
    res = {}
    for nm, path in (("Tomkin S RT 1.0 AMPL 90", os.path.join(ROOT, "irs", "EMT 246 captures", "96kHz", "Short RT", "12_0_90_10_10_10_dc.wav")),
                     ("CLA", os.path.join(ROOT, "irs", "eugtone_EMT246_CLA_shortRT.wav"))):
        x, fs = sf.read(path, always_2d=True)
        h, info = device_ir(x, fs)
        res[nm] = info
    h, info = altiverb_h("1.0")
    res["Altiverb S RT 1.0"] = info
    for k, v in res.items():
        print(f"{k:26s} tap-0 sign L/R {v['tap0_sign_LR']}  channel offsets {np.round(v['channel_offsets_us'], 2)} us")
    return res


STAGES = {"predict": stage_predict, "direct": stage_direct, "diag": stage_diag}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=["all", *STAGES])
    ap.add_argument("--out", default=os.path.join(ROOT, "analysis", "altiverb"))
    a = ap.parse_args()
    path = os.path.join(a.out, "predict.json")
    summary = json.load(open(path)) if os.path.exists(path) else {}
    for name, fn in STAGES.items():
        if a.stage in ("all", name):
            print(f"== {name}")
            summary[name] = fn(a.out)
            json.dump(summary, open(path, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
