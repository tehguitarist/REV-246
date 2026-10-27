"""With NL_WARP2=1 and the exact resampler gain in the model: is the small complex residual above 11.5 kHz worth keeping?
Pooled over captures (train = even index of a seeded sample, test = odd), complex gain on knots fitted on train, scored on test.
The bench_warp2 fits (o, ppm, gain, front shift) are used read-only. usage: resid_eq.py [n_captures]"""
import os, sys, json, random
os.environ["NL_WARP2"] = "1"
os.environ.pop("NL_BK", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter, sosfiltfilt
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
import defects
random.seed(11)
ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
sample = random.sample(ks, int(sys.argv[1]) if len(sys.argv) > 1 else 160)
knots = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
SOS = {nm: butter(4, [lo, hi], btype="band", fs=FS, output="sos") for nm, (lo, hi) in dict(MID=(5000, 10000), HI1=(9000, 11000), S1=(11000, 11600), S2=(11600, 12300), S3=(12300, 12800), C=(10000, 12000)).items()}
rec = []
for i, key in enumerate(sample):
    b = common.bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
    a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if b2 - a < int(0.05 * FS): continue
    n = b2 + 200
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = nlt.unwarp(rawp, b["ppm"])[a:b2]
    m = (M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"]))[a:b2]
    rec.append(dict(key=key, y=y, m=m))
print("captures", len(rec), "warp2", nlt.WARP2, "EQ file used:", False)
def fit(rs, cplx):
    num = np.zeros(len(knots), complex); den = np.zeros(len(knots))
    for r in rs:
        L = len(r["m"]); w = np.hanning(L); f = np.fft.rfftfreq(L) * FS
        Y = np.fft.rfft(r["y"] * w[:, None], axis=0); Mm = np.fft.rfft(r["m"] * w[:, None], axis=0)
        for ci in (0, 1):
            x = Y[:, ci] * np.conj(Mm[:, ci]); d = np.abs(Mm[:, ci]) ** 2
            for kk in range(len(knots)):
                lo = knots[kk - 1] if kk else knots[0] - 1; hi = knots[kk + 1] if kk < len(knots) - 1 else knots[-1] + 1
                wg = np.clip(np.minimum((f - lo) / (knots[kk] - lo), (hi - f) / (hi - knots[kk])), 0, 1)
                num[kk] += np.sum(wg * x); den[kk] += np.sum(wg * d)
    g = num / den
    return g if cplx else g.real + 0j
def apply(r, g):
    L = len(r["m"]); f = np.fft.rfftfreq(L) * FS
    kx = np.r_[9500.0, knots]; E = np.interp(f, kx, np.r_[1.0, g.real]) + 1j * np.interp(f, kx, np.r_[0.0, g.imag])   # unity below 9.5 kHz
    return np.fft.irfft(np.fft.rfft(r["m"], axis=0) * E[:, None], L, axis=0)
def null(y, m, nm):
    yb, mb = sosfiltfilt(SOS[nm], y, axis=0), sosfiltfilt(SOS[nm], m, axis=0)
    return np.array([10 * np.log10(np.sum((yb[:, c] - mb[:, c]) ** 2) / np.sum(yb[:, c] ** 2)) for c in (0, 1)])
out = {}
for tag, tr_i, te_i in (("fold A: train even / test odd", 0, 1), ("fold B: train odd / test even", 1, 0)):
    train = rec[tr_i::2]; test = rec[te_i::2]
    gr, gc = fit(train, False), fit(train, True)
    print(tag, "| real residual dB at knots", np.round(20 * np.log10(np.abs(gr)), 2), "| complex phase mrad", np.round(np.angle(gc) * 1e3, 1))
    for nm in ("MID", "HI1", "S1", "S2", "S3", "C"):
        row = {}
        for vn, g in (("exact gain only", None), ("+ real residual", gr), ("+ complex residual", gc)):
            v = np.array([null(r["y"], r["m"] if g is None else apply(r, g), nm) for r in test])
            row[vn] = v
        out[(tag, nm)] = row
        print("   %-4s test median L/R: " % nm + " | ".join("%s %.1f/%.1f" % (vn, np.median(v[:, 0]), np.median(v[:, 1])) for vn, v in row.items())
              + " | worse than exact by >0.5 dB (real, cplx): %d, %d of %d" % tuple([int(((row[k] - row["exact gain only"]) > 0.5).sum()) for k in ("+ real residual", "+ complex residual")] + [row["exact gain only"].size]))
# who gets worse with the complex residual? (fold B, band C)
train = rec[1::2]; test = rec[0::2]; gc = fit(train, True)
rows = []
for r in test:
    e, c = null(r["y"], r["m"], "C"), null(r["y"], apply(r, gc), "C")
    rows.append((r["key"], e.mean(), c.mean(), float(10 * np.log10(np.mean(r["y"] ** 2) + 1e-30))))
worse = [x for x in rows if x[2] - x[1] > 0.5]
print("worse by >0.5 dB:", len(worse), "of", len(rows))
def grp(k): return ("EARLY<0" if k[1] < 0 else "EARLY10" if k[1] >= 10 else "EARLY>0") + (" LATE0" if k[2] == 0 else " LATE>0")
import collections
tot = collections.Counter(grp(x[0]) for x in rows); wc = collections.Counter(grp(x[0]) for x in worse)
print({g: f"{wc[g]}/{tot[g]}" for g in tot})
print("median exact null of the worse group %.1f dB, of the rest %.1f dB (band C); mean level dB: worse %.1f rest %.1f" % (
    np.median([x[1] for x in worse]), np.median([x[1] for x in rows if x not in worse]), np.mean([x[3] for x in worse]), np.mean([x[3] for x in rows if x not in worse])))
print("worst examples:", [(x[0], round(x[1], 1), round(x[2], 1)) for x in sorted(worse, key=lambda x: x[1] - x[2])[:5]])
