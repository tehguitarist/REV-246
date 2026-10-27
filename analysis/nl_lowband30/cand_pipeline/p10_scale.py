"""Optimal scale s of the DC step response in e = y - m (20-100 Hz band, LB model per env), per time range; s = 1 means the sim fully explains it."""
import sys, json; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=220, precision=2, suppress=True)
B20 = butter(4, [20, 100], btype="band", fs=FSD, output="sos")
RNG = ((0.0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.6))
S = []; keys = [tuple(float(v) for v in a.split(",")) for a in sys.argv[1:]]
for key in keys:
    r = align(key); fr = r["front_shift"]; x, _ = sf.read(file_path(key), always_2d=True); nf = len(x)
    raw = nl.ir(key); h = sosfiltfilt(ne.BAND, raw, axis=0); n = min(M.usable_end(h), K0 + int(fr) + int(2 * FSD)); i0, _ = M.window_start(h, fr)
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))]); yf = nlt.unwarp(rawp, r["ppm"])[:n]
    mf = M.render(key, n, r["o"], fr, band=False) * np.array(r["gain"])
    hdc, _ = pipeline(np.ones((nf, 2)) * x[:30].mean(0)); hp = np.vstack([hdc[:n + 64], np.zeros((max(0, n + 64 - len(hdc)), 2))]); corr = nlt.unwarp(hp, r["ppm"])[:n]
    eb, cb = sosfiltfilt(B20, yf - mf, axis=0), sosfiltfilt(B20, corr, axis=0)
    row = []
    for t0, t1 in RNG:
        a, b = max(i0, K0 + int(fr) + int(t0 * FSD)), min(n, K0 + int(fr) + int(t1 * FSD))
        row.append([(eb[a:b, c] @ cb[a:b, c]) / (cb[a:b, c] @ cb[a:b, c]) if b - a > 100 else np.nan for c in (0, 1)])
    S.append(row)
S = np.array(S)   # keys, ranges, ch
for i, (t0, t1) in enumerate(RNG):
    s = S[:, i, :].ravel(); s = s[np.isfinite(s)]
    print(f"{t0}-{t1} s: optimal scale median {np.median(s):.2f}  IQR {np.percentile(s,25):.2f}..{np.percentile(s,75):.2f}  n {len(s)}")
