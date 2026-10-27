import sys, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
for key, shift in [((0.5, -0.2, 0.0), 9008), ((0.4, -0.2, 0.0), 5243), ((0.5, -0.5, 0.0), 5244), ((0.5, 10.0, 0.0), 0), ((0.5, 0.1, 0.0), 64)]:
    y, m, info = pair(key)
    print(key, "o %.1f ppm %.2f i0 %d n %d (response time of file t=0: %.3f s)" % (info["b"]["o"], info["b"]["ppm"], info["i0"], info["n"], -info["b"]["o"] / FS))
    for lo, hi in [(10000, 11000), (11000, 11600), (11600, 12300), (12300, 12800)]:
        r = local_delay(y[:, 0], m[:, 0], lo, hi, win=2048, hop=2048, start=info["i0"] + 40)
        sel = r[r[:, 4] > r[:, 4].max() - 40][:7]
        print("  %5d-%5d " % (lo, hi) + " ".join("%.2f:%+.3f" % (t, tp) for t, tp, tg, c, e in sel))
