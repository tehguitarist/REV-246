import sys
from lb_common import *
B = bench_all()
np.set_printoptions(linewidth=220, precision=1, suppress=True)
CL = {}
for k, r in sorted(B.items()):
    if k[1] not in (10.0, 3.0) or k[2] not in (0.0, 0.2, 0.4, 0.6) or r["trimmed"] or abs(r["o"] - 0.5) > 1.2: continue
    y, m, _ = load(k, r, n=1500)
    z = int(round(K0 + r["front_shift"])); pk = np.abs(m[z:z + 800, 0]).max()
    CL.setdefault(round(r["front_shift"]), []).append((y - m) / pk * 1e4)
idx = np.arange(34, 70, 2)
print("absolute cache index (file sample 0 = index 50):", idx)
for key in sorted(CL):
    E = np.array(CL[key]).mean(0); n = len(CL[key])
    for c in (0, 1):
        print(f"fs {key:3d} n{n:2d} ch{c} mean e 1e-4:", E[idx, c], "  | std over captures at idx 40:", np.round(np.array(CL[key])[:, 40, c].std(), 1))
