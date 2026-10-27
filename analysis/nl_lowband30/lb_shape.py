import sys
from lb_common import *
B = bench_all()
np.set_printoptions(linewidth=220, precision=1, suppress=True)
CL = {}
for k, r in sorted(B.items()):
    if k[1] not in (10.0, 3.0) or k[2] not in (0.0, 0.2, 0.4, 0.6) or r["trimmed"] or abs(r["o"] - 0.5) > 1.2: continue
    y, m, _ = load(k, r, n=3000)
    yf, mf = frame(y, r["front_shift"]), frame(m, r["front_shift"])
    pk = np.abs(mf[K0:K0 + 800, 0]).max()
    CL.setdefault((round(r["front_shift"]), k[1]), []).append((yf - mf) / pk * 1e4)
idx = np.array([-50, -40, -30, -20, -10, 0, 10, 20, 30, 38, 42, 46, 60, 100, 200, 300, 400, 408, 412, 460, 466, 500, 600, 680]) + K0
print("columns: front-frame sample t (front taps at 41, 410, 463 (R: 41, 278, 463)); values: cluster-mean full-band e, 1e-4 of the 410-tap peak (cluster-mean over captures)")
print("t:", (idx - K0))
for key in sorted(CL):
    E = np.array(CL[key]).mean(0)
    for c in (0, 1):
        print(f"fs {key[0]:3d} EARLY {key[1]:4.1f} n{len(CL[key])} ch{c}: ", E[idx, c])
