"""Repeat-capture ceiling per sub-band, early section only (SPAN 0.5, same EARLY, LATE differs; window 0.05-0.45 s).
Pair (i, i+1) by LATE. Fit one fractional delay + gain on band A (8-11 kHz), score all bands (so HF bands are held out of the fit)."""
import sys, json, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
from scipy.signal import butter, sosfiltfilt
from scipy.optimize import minimize_scalar
SOS = {nm: butter(4, [lo, hi], btype="band", fs=FS, output="sos") for nm, (lo, hi) in dict(LO=(1000, 5000), MID=(5000, 10000), A=(8000, 11000), HI2=(11000, 12000), HI3=(11800, 12800), C=(10000, 12000), B=(11000, 12800)).items()}
def load(key, unw):
    b = bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; a = M.K0 + int(fr) + int(0.05 * FS); e = M.K0 + int(fr) + int(0.45 * FS)
    rawp = raw[:e + 64]
    return unw(rawp, b["ppm"])[a:e], b
def shift(x, tau):
    f = np.fft.rfftfreq(len(x))
    return np.fft.irfft(np.fft.rfft(x, axis=0) * np.exp(-2j * np.pi * f * tau)[:, None], len(x), axis=0)
def nullb(y, m, nm):
    yb, mb = sosfiltfilt(SOS[nm], y, axis=0), sosfiltfilt(SOS[nm], m, axis=0)
    return np.array([10 * np.log10(np.sum((yb[:, c] - mb[:, c]) ** 2) / np.sum(yb[:, c] ** 2)) for c in (0, 1)])
res = {"old": [], "new": []}
for early in (10.0, 0.5, 1.0, -1.0):
    lates = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    for tag, unw in (("old", nlt.unwarp), ("new", unwarp2)):
        ys = [load((0.5, early, l), unw)[0] for l in lates]
        for i in range(len(ys) - 1):
            y1, y2 = ys[i], ys[i + 1]
            def cost(t):
                y2s = shift(y2, t); g = (y1 * y2s).sum() / (y2s * y2s).sum()
                return float(np.mean(nullb(y1, g * y2s, "A")))
            r = minimize_scalar(cost, bounds=(-0.3, 0.3), method="bounded", options={"xatol": 1e-4})
            y2s = shift(y2, r.x); g = (y1 * y2s).sum() / (y2s * y2s).sum()
            res[tag].append(dict(early=early, pair=(lates[i], lates[i + 1]), tau=float(r.x), **{nm: nullb(y1, g * y2s, nm).tolist() for nm in SOS}))
for tag in ("old", "new"):
    print(tag, "pairs", len(res[tag]), " median repeat null dB (L, R):")
    for nm in SOS:
        v = np.array([r[nm] for r in res[tag]]); print("   %-4s %s" % (nm, np.round(np.median(v, 0), 1)))
    print("   median |tau| between repeats %.4f samples" % np.median([abs(r["tau"]) for r in res[tag]]))
json.dump(res, open(OUT + "/repeat_ceiling.json", "w"), default=float)
