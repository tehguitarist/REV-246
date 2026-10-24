"""The Altiverb MED / LONG presets (programs M and L of another unit) against the M / L model (ml_structure §16b): the pulse
check's M / L metrics, nulls and G6, per preset. Nothing of the model is fitted on the Altiverb; per preset and channel only its
chain: an 11-tap FIR on the cluster-only region (as altiverb_predict --stage ml), a cluster gain and a reverb gain, and the
clock (the plugin's +42.8 ppm from S, refined on a +-6 ppm grid by the recirculation-1 null).

    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 8 -- .venv/bin/python tools/altiverb_ml.py
                                                       # -> analysis/altiverb/ml_pulse.{json,log}

Settings: RT from the preset name; LF / HF as inferred from these files' decay by altiverb_compare --stage settings and snapped
to one value per family (MED LF 1.1 / HF 0.7, LONG LF 1.2 / HF 0.4; library_notes §6). So the decay metrics are not held out
(the settings were read from the same decay); AP damping step b is chosen per preset from S's four values on the first
pass (the presets are not library steps), so the first pass is not held out where b > 0; the recirculation nulls and the structural G6 metrics (EDR,
coloration, echo density, IACC) were not fitted. Each preset is one frozen realisation of a time-variant program (θ ≠ 0), so
the reverb nulls are expected near the library's model-vs-capture level, not deeper. Usable end: where the R channel's 50 ms
energy falls 60 dB below its 50-300 ms level, at most 5 s (the files are good to about -74 dB).
"""
import glob, json, os, sys

import numpy as np
from scipy.signal import sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import analyze_ir as ai  # noqa: E402
import ml_carryover as cc  # noqa: E402
import ml_lines as mll  # noqa: E402
import ml_structure as mls  # noqa: E402
import pulse_check as pc  # noqa: E402
import q1_generalise as q1  # noqa: E402

FS, K0 = mll.FS, mll.K0
SET = {"M": ("MED RT", 1.1, 0.7), "L": ("LONG RT", 1.2, 0.4)}
KH = 5                                                       # 11-tap chain FIR per channel
EPS = 42.8                                                   # the plugin's clock, ppm (altiverb_predict, from S)
B_STEPS = (0.0, 0.17, 0.33, 0.50)                           # AP230-type damping steps (model_s_results §14.3)
CACHE = os.path.join(ROOT, "analysis", "altiverb", "cache")
OUT = os.path.join(ROOT, "analysis", "altiverb", "ml_pulse")


def device_ir(name, n):
    """Preset (cache name) -> device-rate IR (n x 2), cluster tap 0 at K0, the plugin clock removed, band-limited."""
    ir = np.load(os.path.join(CACHE, name + ".npy")).T
    fa = int(np.argmax(np.abs(ir[:, 0]) > np.abs(ir[:, 0]).max() * 10 ** (-40 / 20)))
    x = ir[fa - 240: fa + int((n / FS + 0.2) * 48000)]
    h, _, _ = ai.recover(x, 48000, ai.FS_DEV)
    k0 = ai.cluster_taps(h)[0]["first_tap_index"]
    h = q1.apply_clock(h * (np.sign(h[k0].sum()) or 1.0), EPS, 0.0)
    h = h[k0 - K0: k0 - K0 + n]
    if len(h) < n:
        h = np.vstack([h, np.zeros((n - len(h), 2))])
    return sosfiltfilt(mls.BAND, h, axis=0)


def usable(y):
    e = np.convolve(y[:, 1] ** 2, np.ones(int(0.05 * FS)), "same")
    ref = np.mean(e[K0 + int(0.05 * FS): K0 + int(0.3 * FS)])
    below = np.nonzero(e[K0 + int(0.3 * FS):] < ref * 1e-6)[0]
    return K0 + int(0.3 * FS) + (below[0] if len(below) else len(e) - K0 - int(0.3 * FS) - 1)


def fit(prog, y, C, R, i1):
    """Per channel: chain FIR on the cluster-only region, then cluster and reverb gains on [K0 - 30, i1)."""
    out = np.zeros_like(y)
    for ch in (0, 1):
        sl = slice(K0 - 30, K0 + cc.ONSET[prog][ch] - 3 - KH)
        X = np.stack([np.roll(C[:, ch], k)[sl] for k in range(-KH, KH + 1)], 1)
        h = np.linalg.lstsq(X, y[sl, ch], rcond=None)[0]
        f = lambda a: sum(h[k + KH] * np.roll(a, k) for k in range(-KH, KH + 1))
        Cf, Rf = f(C[:, ch]), f(R[:, ch])
        s = slice(K0 - 30, i1)
        g = np.linalg.lstsq(np.stack([Cf[s], Rf[s]], 1), y[s, ch], rcond=None)[0]
        out[:, ch] = g[0] * Cf + g[1] * Rf
    return out


def run_one(prog, name):
    pre, lf, hf = SET[prog]
    rt = float(name.split(" ")[-1][:-1])
    n = K0 + int(min(5.0, 1.2 * rt + 0.5) * FS)
    y = device_ir(name, n)
    i1 = min(usable(y), n - 1)
    key = f"{lf}_{rt}_{hf}"
    M = mll.Model(prog)
    fb = None                                               # AP damping step: S's four values, chosen on the first pass
    for b in B_STEPS:                                       # (the presets are not library steps; short RT uses b > 0)
        M.b[key] = b
        C, R = M.parts(n, key)
        fpn = pc._ml_nulls(prog, y, fit(prog, y, C, R, i1), i1)["first_pass"]
        if fb is None or fpn < fb[1]:
            fb = (b, fpn, C, R)
    b, _, C0, R0 = fb
    best = None
    for d in np.arange(-6.0, 6.01, 1.0):                   # clock refinement around the plugin's
        C, R = (mll.warp(C0, 0.0, d), mll.warp(R0, 0.0, d)) if d else (C0, R0)
        mdl = fit(prog, y, C, R, i1)
        nl = pc._ml_nulls(prog, y, mdl, i1)
        if best is None or (nl["recirc1"] or 0) < (best[1]["recirc1"] or 0):
            best = (d, nl, mdl)
    d, nl, mdl = best
    v = []                                                  # cluster and gap (in-sample: the chain FIR's region)
    for ch in (0, 1):
        s = slice(K0 - 30, K0 + cc.ONSET[prog][ch] - 3)
        ya, ma = y[s, ch], mdl[s, ch]
        v.append(10 * np.log10(np.sum((ya - np.dot(ya, ma) / np.dot(ma, ma) * ma) ** 2) / np.sum(ya ** 2)))
    nl = dict(cluster_gap=float(np.mean(v)), **nl)
    t0 =(cc.ONSET[prog][1] + 20) / FS
    m = pc.compare(pc.prep(mdl[:i1], dtype=np.float32), pc.prep(y[:i1], dtype=np.float32), i1, K0, t0, (t0, t0 + 0.25))
    print(f"{prog} {name}: b {b:.2f}, clock {EPS + d:+.1f} ppm, end {(i1 - K0) / FS:.2f} s | " + " ".join(f"{w} {x:.1f}" for w, x in nl.items() if x is not None), flush=True)
    return dict(prog=prog, preset=name, rt=rt, lf=lf, hf=hf, b=b, ppm=EPS + d, t_end=(i1 - K0) / FS, null=nl, m=m)


def main():
    res = []
    for prog, (pre, _, _) in SET.items():
        names = sorted((os.path.basename(p)[:-4] for p in glob.glob(os.path.join(CACHE, pre + " *.npy"))), key=lambda t: float(t.split(" ")[-1][:-1]))
        res += [run_one(prog, nm) for nm in names]
    json.dump(res, open(OUT + ".json", "w"))
    res = json.loads(json.dumps(res))                       # summarise expects the JSON form (string keys)
    lib = {p: [json.load(open(f)) for f in glob.glob(os.path.join(pc.REF, "ml", f"{p}_*.json"))] for p in SET}
    out = ["## Altiverb MED / LONG against the M / L model (static, θ = 0): per preset; library model vs capture "
           "[library capture pair] medians for reference (analysis/pulse/ml)"]
    for prog in SET:
        R_ = [r for r in res if r["prog"] == prog]
        L_ = lib[prog]
        out.append(f"\n### {prog}: {len(R_)} Altiverb presets")
        ws = ["cluster_gap", "first_pass", "recirc1", "recirc2", "recirc3"]
        lm = {w: [r["null"][w] for d in L_ for r in d["rows"] if r["null"].get(w) is not None] for w in ws}
        lp = {w: [p["null"][w] for d in L_ for p in d["pairs"] if p["null"].get(w) is not None] for w in ws}
        for r in R_:
            out.append(f"- {r['preset']} (b {r['b']:.2f}, end {r['t_end']:.1f} s): " + "; ".join(f"{w} {r['null'][w]:.1f}" for w in ws if r["null"].get(w) is not None))
        out.append("- library: " + "; ".join(f"{w} {np.median(lm[w]):.1f} [{np.median(lp[w]):.1f}]" for w in ws if lm[w] and lp[w]))
        rows = [pc.summarise(r["m"]) for r in R_]
        lrows = [pc.summarise(r["m"]) for d in L_ for r in d["rows"]]
        lreps = [pc.summarise(p["m"]) for d in L_ for p in d["pairs"]]
        md = lambda rr, k, f=abs: np.median([f(r[k]) for r in rr if r.get(k) is not None and np.isfinite(r[k])] or [np.nan])
        out.append("- octave decay |ratio - 1| % (not held out: LF / HF read from these decays), Altiverb {library model} [library pair]: "
                   + " ".join(f"{f}: {100 * md(rows, f'oct_t{f}'):.1f} {{{100 * md(lrows, f'oct_t{f}'):.1f}}} [{100 * md(lreps, f'oct_t{f}'):.1f}]" for f in pc.OCT))
        out.append("- largest EDC difference (dB): " + " ".join(f"{f}: {md(rows, f'oct_dev{f}'):.2f} {{{md(lrows, f'oct_dev{f}'):.2f}}} [{md(lreps, f'oct_dev{f}'):.2f}]" for f in pc.OCT))
        out.append(f"- echo density |diff| {md(rows, 'ned'):.3f} {{{md(lrows, 'ned'):.3f}}} [{md(lreps, 'ned'):.3f}]; fine structure R "
                   + " ".join(f"{lo}-{hi}: {md(rows, f'fine_{lo}-{hi}', lambda x: x):.3f} {{{md(lrows, f'fine_{lo}-{hi}', lambda x: x):.3f}}}" for lo, hi in pc.FS_BANDS))
        out.append("- G6, Altiverb (model - Altiverb) [library capture pair]:")
        out += pc.g6_lines(rows, lreps, indent="  - ")
        out.append("- G6, library (model - capture) [library capture pair], for comparison:")
        out += pc.g6_lines(lrows, lreps, indent="  - ")
    open(OUT + ".log", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
