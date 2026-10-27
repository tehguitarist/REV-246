"""Fork C, step 2: per-channel complex residual on the REVERB of captures disjoint from the front-tap set (seed 23, c_front.py),
and the held-out test of the front-tap-derived correction (applied with unity below 9.5 kHz, nothing fitted on these captures).
usage: c_reverb.py N [front_pkl] [variant]"""
import os, sys, json, random, pickle
os.environ["NL_WARP2"] = "1"; os.environ.pop("NL_BK", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter, sosfiltfilt
import c_front as CF
import defects
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 150
pkl = sys.argv[2] if len(sys.argv) > 2 else "c_front_150.pkl"
variant = sys.argv[3] if len(sys.argv) > 3 else "all"
KNOTS = CF.KNOTS
front = pickle.load(open(OUT + "/" + pkl, "rb"))[variant]
fkeys = {tuple(r["key"]) for r in front}
corr = {ch: CF.pooled(front, ch) for ch in (0, 1)}
ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and k not in fkeys and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
random.seed(29); sample = random.sample(ks, N)
SOS = {nm: butter(4, [lo, hi], btype="band", fs=FS, output="sos") for nm, (lo, hi) in dict(
    B20=(20, 100), B100=(100, 1000), B1k=(1000, 5000), B5k=(5000, 10000), HI1=(9000, 11000), S1=(11000, 11600), S2=(11600, 12300), S3=(12300, 12800), C=(10000, 12000)).items()}
def null(y, m, nm, a=None, b2=None):
    """a, b2: y and m are the FULL arrays, filtered whole and cut afterwards (the low bands see no segment-edge transient)"""
    yb, mb = sosfiltfilt(SOS[nm], y, axis=0), sosfiltfilt(SOS[nm], m, axis=0)
    if a is not None: yb, mb = yb[a:b2], mb[a:b2]
    return np.array([10 * np.log10(np.sum((yb[:, c] - mb[:, c]) ** 2) / np.sum(yb[:, c] ** 2)) for c in (0, 1)])
def gain(f, g, ch):
    kx = np.r_[9500.0, KNOTS[:-1]]; gg = g[:-1]          # knots to 12.8 kHz; unity below 9.5 kHz and above 12.8 kHz (the analysis band ends)
    return np.interp(f, kx, np.r_[1.0, gg.real], right=1.0) + 1j * np.interp(f, kx, np.r_[0.0, gg.imag], right=0.0)
def apply(m, gs):
    """the correction on the whole model, zero-padded (no circular wrap), then the scored segment is cut by the caller"""
    n = len(m); Np = 1 << int(np.ceil(np.log2(2 * n))); f = np.fft.rfftfreq(Np) * FS; Mf = np.fft.rfft(m, Np, axis=0)
    return np.fft.irfft(np.stack([Mf[:, c] * gain(f, gs[c], c) for c in (0, 1)], 1), Np, axis=0)[:n]
rec = []
for key in sample:
    b = common.bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
    a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if b2 - a < int(0.05 * FS): continue
    n = b2 + 200
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    yfull = nlt.unwarp(rawp, b["ppm"])[:n]; y = yfull[a:b2]
    mf = M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"])
    rec.append(dict(key=key, y=y, yf=yfull, m=mf[a:b2], mf=mf, a=a, b2=b2, ppm=b["ppm"]))
print("reverb captures", len(rec), "front set", len(fkeys), "disjoint", not (fkeys & {r["key"] for r in rec}))
# per-channel complex residual on the reverb
num = {c: np.zeros(len(KNOTS), complex) for c in (0, 1)}; den = {c: np.zeros(len(KNOTS)) for c in (0, 1)}
for r in rec:
    L = len(r["m"]); w = np.hanning(L); f = np.fft.rfftfreq(L) * FS
    Y = np.fft.rfft(r["y"] * w[:, None], axis=0); Mm = np.fft.rfft(r["m"] * w[:, None], axis=0)
    W = CF.hats(f)
    for c in (0, 1):
        num[c] += (W * (Y[:, c] * np.conj(Mm[:, c]))[None, :]).sum(1); den[c] += (W * (np.abs(Mm[:, c]) ** 2)[None, :]).sum(1)
rv = {c: num[c] / den[c] for c in (0, 1)}
print("knots kHz", np.round(KNOTS / 1000, 2))
for c in (0, 1):
    print(f"ch{c} REVERB level dB ", np.round(20 * np.log10(np.abs(rv[c])), 3)); print(f"ch{c} FRONT  level dB ", np.round(20 * np.log10(np.abs(corr[c])), 3))
    print(f"ch{c} REVERB phase mrad", np.round(np.angle(rv[c]) * 1e3, 1)); print(f"ch{c} FRONT  phase mrad", np.round(np.angle(corr[c]) * 1e3, 1))
# held-out test of the front correction (and of an R-only version applied to both channels, and phase-only)
cands = {"front per channel": [corr[0], corr[1]], "front R on both": [corr[1], corr[1]],
         "front hybrid (L own to 12.25 kHz, R's above)": [np.r_[corr[0][:8], corr[1][8:]], corr[1]],
         "front phase only": [np.exp(1j * np.angle(corr[0])), np.exp(1j * np.angle(corr[1]))]}
bands = ["B20", "B100", "B1k", "B5k", "HI1", "S1", "S2", "S3", "C"]
base = {nm: np.array([null(r["yf"], r["mf"], nm, r["a"], r["b2"]) for r in rec]) for nm in bands}
res = {}
for cn, gs in cands.items():
    cm = [apply(r["mf"], gs) for r in rec]
    new = {nm: np.array([null(r["yf"], cm[i], nm, r["a"], r["b2"]) for i, r in enumerate(rec)]) for nm in bands}
    res[cn] = new
    print(f"--- correction: {cn} (held out: derived on {len(fkeys)} other captures' front taps), test captures {len(rec)}")
    for nm in bands:
        d = new[nm] - base[nm]
        print(f"   {nm:4s} median L/R before {np.median(base[nm][:,0]):6.1f}/{np.median(base[nm][:,1]):6.1f} after {np.median(new[nm][:,0]):6.1f}/{np.median(new[nm][:,1]):6.1f} | diff median {np.median(d):+.2f} dB; worse >+0.5 dB: {int((d > 0.5).sum())}/{d.size}; worst {d.max():+.2f}; better >0.5: {int((d < -0.5).sum())}")
# who gets worse in band C (front R on both)?
import collections
d = res["front R on both"]["C"] - base["C"]
def grp(k): return ("EARLY<0" if k[1] < 0 else "EARLY10" if k[1] >= 10 else "EARLY>0") + (" LATE0" if k[2] == 0 else " LATE>0")
tot = collections.Counter(); wc = collections.Counter()
for i, r in enumerate(rec):
    for c in (0, 1):
        tot[grp(r["key"])] += 1; wc[grp(r["key"])] += int(d[i, c] > 0.5)
print("worse >0.5 dB in C by group (front R on both):", {g: f"{wc[g]}/{tot[g]}" for g in tot})
wb = [(i, c) for i in range(len(rec)) for c in (0, 1) if d[i, c] > 0.5]
print("exact null before of the worse ones: median %.1f dB; of the rest %.1f dB" % (np.median([base["C"][i, c] for i, c in wb]), np.median([base["C"][i, c] for i in range(len(rec)) for c in (0, 1) if (i, c) not in wb])))
# what separates the captures that get worse? (front R on both, band C, mean over L / R)
dm = d.mean(1); span = np.array([r["key"][0] for r in rec]); late = np.array([r["key"][2] for r in rec]); dppm = np.array([r["ppm"] for r in rec]) - 11.6
from scipy.stats import spearmanr
print("improvement (band C, dB, negative = better) vs SPAN rho %.2f, vs LATE %.2f, vs ppm-11.6 %.2f, vs |ppm-11.6| %.2f, vs exact null before %.2f" % (
    spearmanr(dm, span)[0], spearmanr(dm, late)[0], spearmanr(dm, dppm)[0], spearmanr(dm, np.abs(dppm))[0], spearmanr(dm, base["C"].mean(1))[0]))
for sp in sorted(set(span)): print("   SPAN", sp, "median diff %+.2f dB, worse >0.5: %d/%d" % (np.median(dm[span == sp]), int((dm[span == sp] > 0.5).sum()), int((span == sp).sum())))
pickle.dump(dict(rv=rv, corr=corr), open(OUT + "/c_reverb_" + variant + ".pkl", "wb"))
