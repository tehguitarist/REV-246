"""20 Hz-1 kHz error of the model in the repeat window: relative error spectrum E/M (coherent mean, rms) and the repeat floor in the same bins,
time course, step-removal effect.  SPAN 0.5, EARLY >= 0.7, LATE <= 1.0 (late starts after the window), bench_bk_a0.65 alignment.
env: NL_BK=1 NL_BK_ALPHA=0.65 LB_BENCH=bench_bk_a0.65"""
import sys, pickle, itertools
from lb_common import *
B = bench_all()
span = 0.5
keys = sorted(k for k in B if k[0] == span and k[1] >= 0.7 and k[2] <= 1.0 and not B[k]["trimmed"])
n = 14000
sel = {}
for k in keys:
    r = B[k]
    y, m, _ = load(k, r, n=n, band=False)
    z = K0 + r["front_shift"]
    sel[k] = dict(y=y, m=m, z=z, o=r["o"], fs=r["front_shift"], pk=np.abs(m[int(z):int(z) + 800, 0]).max())
pickle.dump(sel, open("low_sel.pkl", "wb"))
print(len(sel), "captures; fs:", sorted({round(v["fs"], 1) for v in sel.values()}))
