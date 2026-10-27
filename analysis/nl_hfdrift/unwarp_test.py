import sys, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
key = (0.0, 10.0, 2.0)
y, m, info = pair(key)
ppm = 12.0
z = nlt.unwarp(np.vstack([m, np.zeros((64, 2))]), ppm)[:len(m)]
print("unwarp(model) vs model: delay per band (the stretch itself is common to all bands)")
for lo, hi in [(6000, 8000), (9500, 10500), (11000, 11600), (11600, 12300), (12300, 12800)]:
    r = local_delay(z[:, 0], m[:, 0], lo, hi, win=4096, hop=4096, start=info["i0"])
    sel = r[r[:, 4] > r[:, 4].max() - 40][:6]
    print("  %5d-%5d " % (lo, hi) + " ".join("%.2f:%+.3f" % (t, tp) for t, tp, tg, c, e in sel))
print("expected common delay = -ppm*1e-6*t*FS:", " ".join("%.2f:%+.3f" % (t, -ppm * 1e-6 * t * FS) for t in (0.08, 0.24, 0.39, 0.55)))
# stretch via exact FFT zero-padded sinc interpolation, chunk-wise, as reference (chunks of 8192 with margin)
