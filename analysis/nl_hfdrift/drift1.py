import sys, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
keys = [(0.0, 10.0, 2.0), (0.5, 10.0, 2.0), (0.3, 1.0, 1.0)]
bands = [(6000, 8000), (8000, 9500), (9500, 10500), (10500, 11000), (11000, 11600), (11600, 12300), (12300, 12800)]
for key in keys:
    if tuple(key) not in nl.table():
        print("missing", key); continue
    y, m, info = pair(key)
    print(key, "ppm %.2f o %.3f n %d" % (info["b"]["ppm"], info["b"]["o"], info["n"]))
    for ch in (0, 1):
        for lo, hi in bands:
            r = local_delay(y[:, ch], m[:, ch], lo, hi, win=4096, hop=4096, start=info["i0"])
            sel = r[(r[:, 4] > r[:, 4].max() - 45)]
            s = " ".join("%.2fs:%+.3f(%.2f)" % (t, tp, c) for t, tp, tg, c, e in sel[:8])
            print(" ch%d %5d-%5d %s" % (ch, lo, hi, s))
