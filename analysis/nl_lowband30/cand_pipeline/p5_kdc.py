import sys; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=220, precision=3, suppress=True)
idx = np.array([10, 30, 40, 44, 45, 46, 47, 48, 50, 52, 56, 60, 80, 120, 200, 500, 1000, 3000, 6000, 10000, 16000, 20000])
for nf in (85044, 160000, 502655):
    h, d = pipeline(np.ones((nf, 2)))
    print("nf", nf, "len h", len(h)); print(" idx", idx[idx < len(h)])
    for c in (0, 1): print("  ch", c, "response to unit DC (x1e4)", h[idx[idx < len(h)], c] * 1e4)
