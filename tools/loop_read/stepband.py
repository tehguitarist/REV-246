"""The lag step of the 16 stepped S captures (model_s_results §25.16), per octave band and model-free: each stepped
capture against a step-free capture of the same setting (the other take), relative lag in 40 ms windows (hop 20 ms,
cross-correlation peak, R) per octave band 250 Hz - 4 kHz, from 5 ms after tap 0. A unit that changes its timing
during the reverb steps in time, the same in every band; a slip in the sweep capture shifts the frequencies swept
after it, at every lag, so it steps in frequency and only looks like a time step where the high bands die early
(short RT, strong HF). Then, per setting: each capture's late shift against the setting's other captures (fractional,
fitted at 130-250 ms, 300 Hz-6 kHz) and the step's time to the sample, in reverb time and in file time (from the start
of the raw 96 kHz file); a capture is stepped when its shift against the reference is >= 0.5 device
samples (the reference: the capture that follows the model best). The 9 settings rescored from the pulse check's per-capture nulls (analysis/pulse/s/) with the stepped
captures and the bad AMPL-0 captures (library_notes §2) left out. Then the baseline's HF cells (HF 0.2 / RT 0.3,
HF 0.4 / RT 0.3, HF 0.2 / RT 0.6) re-chosen on step-free captures, by the baseline's rule (per setting from its cell's
other LF steps), from the per-capture hfcell grid (lowband_shelf.py --stage hfcell with SHELF_TAG=_percap), and scored
on the step-free captures with the old and the new value. model_s_results §28.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds SHELF_TAG=_percap \
        .venv/bin/python tools/lowband_shelf.py --stage hfcell --workers 6 --only <the 12 settings of the three cells>

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/stepband.py [keys]     # -> analysis/lowband/chain3_stepband.{json,log}"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.signal import butter, sosfiltfilt
import dense_mid as dm
import warp as W

BANDS = [250, 500, 1000, 2000, 4000]
STEPPED = ["0.1_0.3_0.2", "0.4_0.3_0.2", "1.0_0.3_0.2", "2.0_0.3_0.2", "0.1_0.3_0.4", "0.4_0.3_0.4", "1.0_0.3_0.4",
           "1.0_0.6_0.2", "2.0_0.6_0.2"]
CONTROL = ["1.0_0.3_1.0", "1.0_0.6_1.0"]
WIN, HOP, MAXLAG = 0.040, 0.020, 12
OUT = os.path.join(os.path.dirname(dm.OUT.rstrip("/")).replace("dense_mid", "lowband"), f"{dm.q1.ms.KTAG}_stepband")


def track(a, b, fs, f, i0, i1):
    """Lag of a against b (samples; + = a late) per window in one octave band: (t ms after i0's origin, lag, energy)."""
    sos = butter(3, [f / np.sqrt(2), min(f * np.sqrt(2), 11000)], btype="band", fs=fs, output="sos")
    ab, bb = sosfiltfilt(sos, a), sosfiltfilt(sos, b)
    Wn, H = int(WIN * fs), int(HOP * fs)
    out = []
    for s in range(i0, i1 - Wn, H):
        xs, ms = ab[s:s + Wn], bb[s - MAXLAG:s + Wn + MAXLAG]
        c = np.array([xs @ ms[MAXLAG + k:MAXLAG + k + Wn] for k in range(-MAXLAG, MAXLAG + 1)])
        j = int(np.argmax(c))
        if 0 < j < len(c) - 1:
            d = 0.5 * (c[j - 1] - c[j + 1]) / (c[j - 1] - 2 * c[j] + c[j + 1])
            out.append(((s + Wn / 2 - dm.PRE) / fs * 1000, j - MAXLAG + d, float(xs @ xs)))
    return out


def step_of(z, j):
    """The broadband step of capture j against its model (steps.py's measure)."""
    y, m = z[f"y{j}"][:, 1].astype(float), z[f"m{j}"][:, 1].astype(float)
    i0, i1 = dm.PRE + int(0.06 * dm.FS), min(dm.PRE + int(min(float(z[f"te{j}"]), 0.3) * dm.FS), len(y) - W.MAXLAG)
    T, L, _ = W.lag_track(y, m, dm.FS, i0, i1)
    t = T / dm.FS * 1000 + 60
    e, l = L[(t >= 60) & (t <= 100)], L[(t >= 120) & (t <= 300)]
    return float(np.median(l) - np.median(e)) if len(e) and len(l) else np.nan


def fshift(x, d):
    N = 1 << int(np.ceil(np.log2(len(x) + 64)))
    f = np.fft.rfftfreq(N)
    return np.fft.irfft(np.fft.rfft(x, N) * np.exp(-2j * np.pi * f * d), N)[:len(x)]


def onset_ms(path):
    """Tap 0 in the raw file (ms from its first sample), as ml_structure._make finds the onset."""
    import soundfile as sf
    x, fs = sf.read(path, always_2d=True)
    a = np.abs(x).max(1)
    main = int(np.argmax(a > a.max() * 0.1))
    lo = max(main - int(0.002 * fs), 0)
    return (lo + int(np.argmax(a[lo:main + 1] > a.max() * 1e-2))) / fs * 1000


from defects import BAD  # noqa: E402  library_notes §2


def classify(k, z, caps, amp):
    """Per capture: late shift (device samples) against the setting's reference capture, stepped flag, step time."""
    sos = butter(4, [300, 6000], btype="band", fs=dm.FS, output="sos")
    y = {j: sosfiltfilt(sos, z[f"y{j}"][:, 1].astype(float)) for j in caps}
    good = [j for j in caps if (k, amp[j]) not in BAD]
    i0, i1 = dm.PRE + int(0.13 * dm.FS), dm.PRE + int(0.25 * dm.FS)
    ds = np.linspace(-2.5, 2.5, 501)

    def shift(a, b):
        r = [np.sum((y[a][i0:i1] - np.dot(y[a][i0:i1], v := fshift(y[b], d)[i0:i1]) / np.dot(v, v) * v) ** 2) for d in ds]
        return float(-ds[int(np.argmin(r))])
    ref = min(good, key=lambda j: abs(step_of(z, j)))        # the capture that follows the model best
    sh = {a: (0.0 if a == ref else shift(a, ref)) for a in caps}
    med = 0.0
    out = {}
    names = dm.q2.q1_tables()[0][k]["files"]
    paths = {os.path.basename(p): p for p in dm.q1.ms.files("S")}
    for j in caps:
        d = sh[j] - med
        row = dict(shift=d, stepped=bool(abs(d) >= 0.5 and (k, amp[j]) not in BAD))
        if row["stepped"]:
            ys = fshift(y[ref], -d)
            lo, hi = dm.PRE + int(0.06 * dm.FS), dm.PRE + int(0.16 * dm.FS)
            res = [np.sum((y[j][lo:hi] - np.r_[y[ref][lo:t], ys[t:hi]]) ** 2) for t in range(lo + 50, hi - 50)]
            t = lo + 50 + int(np.argmin(res))
            on = onset_ms(paths[names[str(amp[j])]["file"]])
            row.update(t_reverb_ms=(t - dm.PRE) / dm.FS * 1000, t_file_ms=on + (t - dm.PRE) / dm.FS * 1000, onset_ms=on,
                       shift_96k=d * 96000 / dm.FS)
        out[amp[j]] = row
    return out


def rescore(k, stepped):
    """Median R null per window from the pulse check, all good captures against step-free ones; repeat pairs likewise."""
    P = json.load(open(os.path.join(dm.ROOT, "analysis", "pulse", "s", k + ".json")))
    good = [r for r in P["rows"] if (k, r["amp"]) not in BAD]
    free = [r for r in good if r["amp"] not in stepped]
    ok = lambda p: all((k, a) not in BAD for a in p["pair"])
    pairs = [p for p in P["pairs"] if ok(p)]
    fpairs = [p for p in pairs if not set(p["pair"]) & set(stepped)]
    med = lambda rr, w: (float(np.median([r["null"][w][1] for r in rr if r["null"].get(w)])) if any(r["null"].get(w) for r in rr) else None)
    W = ["50-150", "150-300", "300-600"]
    return dict(model_all={w: med(good, w) for w in W}, model_free={w: med(free, w) for w in W},
                pairs_all={w: med([p for p in pairs if p["cross"]] or pairs, w) for w in W},
                pairs_free={w: med(fpairs, w) for w in W}, free_pairs=[p["pair"] + [p["cross"]] for p in fpairs], n_free=len(free))


GRID = os.path.join(dm.ROOT, "analysis", "lowband", f"{dm.q1.ms.KTAG}_shelf_hfcell_cap0.8475_percap")
CELLS = os.path.join(dm.ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475", "cells.json")


def hf_cells(stepped, res):
    """The three cells' HF pole shifts re-chosen without the stepped and bad captures, and scored on step-free captures."""
    import glob
    D = [json.load(open(p)) for p in glob.glob(os.path.join(GRID, "*.json"))]
    if not D:
        return []
    old = json.load(open(CELLS))
    keep = lambda d: [i for i, a in enumerate(d["amps"]) if (d["key"], a) not in stepped and (d["key"], a) not in BAD]
    rel = {d["key"]: np.array([10 * np.log10(sum(g["costs"][i] for i in keep(d)) / sum(d["base_costs"][i] for i in keep(d))) for g in d["grid"]]) for d in D}
    dps = [g["dp"] for g in D[0]["grid"]]
    out = ["\n# The HF cells re-chosen on step-free captures (held out across LF in each cell); R nulls of the step-free captures, "
           "median, old -> new dp (50-150 / 150-300 / 300-600 ms)"]
    for d in sorted(D, key=lambda d: (d["hf"], d["rt"], d["lf"])):
        oth = [rel[e["key"]] for e in D if e["hf"] == d["hf"] and e["rt"] == d["rt"] and e["key"] != d["key"]]
        new = dps[int(np.argmin(np.mean(oth, 0)))]
        was = old.get(d["key"], {}).get("hf_dp", 0.0)
        iw, inew = int(np.argmin(np.abs(np.array(dps) - was))), dps.index(new)
        nl = lambda j, w: np.median([d["grid"][j]["nulls"][i][w][1] for i in keep(d) if d["grid"][j]["nulls"][i].get(w)] or [np.nan])
        res.setdefault("cells", {})[d["key"]] = dict(old=was, new=new, nulls_old={w: nl(iw, w) for w in ("50-150", "150-300", "300-600")},
                                                    nulls_new={w: nl(inew, w) for w in ("50-150", "150-300", "300-600")}, n=len(keep(d)))
        c = res["cells"][d["key"]]
        out.append(f"- {d['key']}: dp {was:+.4f} -> {new:+.4f} ({c['n']} captures): " + " / ".join(f"{c['nulls_old'][w]:.1f} -> {c['nulls_new'][w]:.1f}" for w in c["nulls_old"]))
    return out


def main(keys):
    res, lines = {}, []
    for k in keys:
        z = np.load(os.path.join(dm.OUT, "prep", k + ".npz"))
        caps = sorted(int(n[3:]) for n in z.files if n.startswith("res"))
        st = {j: step_of(z, j) for j in caps}
        amp = {j: int(z[f"amp{j}"]) for j in caps}
        free = [j for j in caps if abs(st[j]) < 0.5]
        stepped = [j for j in caps if abs(st[j]) >= 0.5]
        pairs = [(a, b) for a in stepped for b in free] if stepped else [(free[0], free[-1])]
        res[k] = {"step": {str(amp[j]): st[j] for j in caps}, "pairs": {}}
        lines.append(f"\n## {k}: broadband step per AMPL " + " ".join(f"{amp[j]}:{st[j]:+.2f}" for j in caps))
        for a, b in pairs:
            ya, yb = z[f"y{a}"][:, 1].astype(float), z[f"y{b}"][:, 1].astype(float)
            i1 = min(dm.PRE + int(min(float(z[f"te{a}"]), float(z[f"te{b}"]), 0.6) * dm.FS), len(ya) - MAXLAG)
            tr = {f: track(ya, yb, dm.FS, f, dm.PRE + int(0.005 * dm.FS), i1) for f in BANDS}
            res[k]["pairs"][f"{amp[a]}-{amp[b]}"] = {str(f): v for f, v in tr.items()}
            lines.append(f"- AMPL {amp[a]} against {amp[b]}: lag (samples) per band in 5-45 / 45-100 / 100-150 / 150-250 / 250-400 ms "
                         "(energy-weighted median; level of the band re its first window)")
            for f, v in tr.items():
                v = np.array(v)
                if not len(v):
                    continue
                e0 = v[0, 2]
                cells = []
                for t0, t1 in ((5, 45), (45, 100), (100, 150), (150, 250), (250, 400)):
                    s = (v[:, 0] >= t0) & (v[:, 0] < t1)
                    cells.append(f"{np.median(v[s, 1]):+5.2f} ({10 * np.log10(v[s, 2].mean() / e0):+4.0f} dB)" if s.sum() else "      --       ")
                lines.append(f"  {f:5d} Hz  " + "  ".join(cells))
    lines.append("\n# Step per capture against the setting's other captures, and the rescore (pulse check, R, median)")
    for k in [k for k in keys if k in STEPPED]:
        z = np.load(os.path.join(dm.OUT, "prep", k + ".npz"))
        caps = sorted(int(n[3:]) for n in z.files if n.startswith("res"))
        amp = {j: int(z[f"amp{j}"]) for j in caps}
        cl = classify(k, z, caps, amp)
        st = [a for a, r in cl.items() if r["stepped"]]
        sc = rescore(k, st)
        res[k].update(classify={str(a): r for a, r in cl.items()}, rescore=sc)
        lines.append(f"\n## {k}: shift per AMPL " + " ".join(f"{a}:{r['shift']:+.2f}" for a, r in cl.items()))
        for a, r in cl.items():
            if r["stepped"]:
                lines.append(f"- AMPL {a}: {r['shift']:+.2f} device samples ({r['shift_96k']:+.2f} at 96 kHz) from "
                             f"{r['t_reverb_ms']:.2f} ms after tap 0 = {r['t_file_ms']:.2f} ms file time (tap 0 at {r['onset_ms']:.2f} ms)")
        f_ = lambda d: " / ".join("--" if d[w] is None else f"{d[w]:.1f}" for w in d)
        lines.append(f"- model null all good captures {f_(sc['model_all'])}; step-free ({sc['n_free']}) {f_(sc['model_free'])}  (50-150 / 150-300 / 300-600 ms)")
        lines.append(f"- repeat null cross-take pairs {f_(sc['pairs_all'])}; step-free pairs {f_(sc['pairs_free'])} {sc['free_pairs']} [pair, cross-take]")
    stepped = {(k, int(a)) for k in res if "classify" in res[k] for a, r in res[k]["classify"].items() if r["stepped"]}
    lines += hf_cells(stepped, res)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT + ".json", "w"))
    open(OUT + ".log", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1].split(",") if len(sys.argv) > 1 else STEPPED + CONTROL)
