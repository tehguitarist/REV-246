"""Program N's low-band transfer T_N(f) tested on Program S (read-only on analysis/fullband/s/*.npz: per setting and channel, per 2 Hz bin, sum Y M*, |M|^2, |Y|^2
over S's reverb 150 ms .. compensated end, S's bench model with its own per-capture clock/shift/gain fit, held out by RT fold).
(a) T_S (pooled over LF >= 1.0, as sml_fullband.fold_T, smoothed +-2 bins) against T_N: complex correlation of (T-1) per band.
(b) the S error energy per band with T = 1, T = T_N, T = T_S(other RT fold) [held out], from the stored sums."""
import glob, os, json, re
import numpy as np
D = "/Users/leighpierce/Git/246/analysis/fullband/s"
TBIN, TNB = 2.0, int(12300.0 / 2.0)
FOLDS = {"A": (0.3, 0.8, 1.4, 2.2), "B": (0.6, 1.0, 1.8, 2.5)}
files = sorted(glob.glob(os.path.join(D, "*.npz")))
S = {}
for p in files:
    k = os.path.basename(p)[:-4]; lf, rt, hf = (float(x) for x in k.split("_")); z = np.load(p)
    S[k] = dict(lf=lf, rt=rt, hf=hf, sym=z["sym"].astype(complex), smm=z["smm"].astype(float), syy=z["syy"].astype(float))
print(len(S), "S settings; LF values", sorted({v["lf"] for v in S.values()}), "RT", sorted({v["rt"] for v in S.values()}))
fc = (np.arange(TNB) + 0.5) * TBIN
z = np.load("/Users/leighpierce/Git/246/analysis/nl_lowband30/fineT_all_sm6.npz"); frN = z["fr"]
TN = [np.interp(fc, frN, z[k].real) + 1j * np.interp(fc, frN, z[k].imag) for k in ("T0", "T1")]
def pooled(ks):
    sym = sum(S[k]["sym"] for k in ks); smm = sum(S[k]["smm"] for k in ks)
    ker = np.ones(5)
    num = np.stack([np.convolve(sym[c], ker, "same") for c in (0, 1)]); den = np.stack([np.convolve(smm[c], ker, "same") for c in (0, 1)])
    return num / np.maximum(den, 1e-30)
fold_T = {f: pooled([k for k, v in S.items() if v["rt"] in FOLDS[f] and v["lf"] >= 1.0]) for f in FOLDS}
TS_all = pooled([k for k, v in S.items() if v["lf"] >= 1.0])
BD = [(20, 40), (40, 80), (80, 160), (160, 270), (270, 400), (400, 800), (800, 1600), (1600, 3200)]
print("\n(a) T_S (S, LF>=1.0 pooled) vs T_N (N, this work): per band, rms|T-1| in % and complex correlation of (T-1); channel L / R")
for lo, hi in BD:
    m = (fc >= lo) & (fc < hi)
    row = []
    for c in (0, 1):
        a, b = TS_all[c][m] - 1, TN[c][m] - 1
        rho = abs(np.vdot(a, b)) / np.sqrt(np.vdot(a, a).real * np.vdot(b, b).real)
        phase = np.degrees(np.angle(np.vdot(a, b)))
        row.append(f"ch{c}: S {100*np.sqrt(np.mean(abs(a)**2)):.2f}% N {100*np.sqrt(np.mean(abs(b)**2)):.2f}% corr {rho:.2f} (phase {phase:.0f} deg)")
    print(f"  {lo:5d}-{hi:<5d} " + " | ".join(row))
def err(v, T, m, c):
    t = T[c][m]
    return np.sum(v["syy"][c][m] + abs(t) ** 2 * v["smm"][c][m] - 2 * np.real(t * np.conj(v["sym"][c][m])))
print("\n(b) S, all settings: median band error energy (dB re |Y|^2) over settings: T=1 | T=T_N | T=T_S of the other RT fold (held out) ; and the per-setting count with T_N worse by > 0.5 dB")
for c in (0, 1):
    print(f" channel {'LR'[c]}")
    for lo, hi in BD:
        m = (fc >= lo) & (fc < hi); r0 = []; rN = []; rS = []
        for k, v in S.items():
            other = "B" if v["rt"] in FOLDS["A"] else ("A" if v["rt"] in FOLDS["B"] else None)
            den = np.sum(v["syy"][c][m])
            if den <= 0: continue
            e0 = err(v, [np.ones(TNB)] * 2, m, c); eN = err(v, TN, m, c)
            r0.append(10 * np.log10(e0 / den)); rN.append(10 * np.log10(max(eN, 1e-30) / den))
            if other: rS.append(10 * np.log10(max(err(v, fold_T[other], m, c), 1e-30) / den))
        r0, rN = np.array(r0), np.array(rN)
        print(f"   {lo:5d}-{hi:<5d} T=1 {np.median(r0):6.1f} | T_N {np.median(rN):6.1f} (worse by >0.5 dB: {(rN - r0 > 0.5).sum():3d}/{len(r0)}, better by >0.5: {(rN - r0 < -0.5).sum():3d}) | T_S held out {np.median(rS):6.1f}")
