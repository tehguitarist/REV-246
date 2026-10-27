import sys, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
key = (0.0, 10.0, 2.0)
y, m, info = pair(key)
ch = 0
edges = np.arange(10500, 12801, 250)
print("phase (rad) vs f at several times; also model/capture amplitude ratio")
for t0 in (0.1, 0.4, 0.8, 1.5, 2.5, 4.0):
    a = int(t0 * FS)
    win = 4096
    w = np.hanning(win)
    Y = np.fft.rfft(y[a:a + win, ch] * w); Mm = np.fft.rfft(m[a:a + win, ch] * w)
    f = np.fft.rfftfreq(win) * FS
    s = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (f >= lo) & (f < hi)
        x = np.sum(Y[k] * np.conj(Mm[k]))
        s.append("%5.0f: ph %+.3f |Y/M| %.3f" % ((lo + hi) / 2, np.angle(x), np.sqrt(np.sum(abs(Y[k]) ** 2) / np.sum(abs(Mm[k]) ** 2))))
    print("t=%.1fs" % t0); print("  " + "\n  ".join(s))
# long-range tau in the 11.6-12.3 band, each 8192 window
for lo, hi in [(11600, 12300), (12300, 12800)]:
    r = local_delay(y[:, ch], m[:, ch], lo, hi, win=8192, hop=8192, start=int(0.02 * FS))
    print(lo, hi, " ".join("%.2f:%+.3f(%.0fdB)" % (t, tp, e) for t, tp, tg, c, e in r))
