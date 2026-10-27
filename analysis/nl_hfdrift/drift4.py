import sys, json, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
M.EQ = False; M._proc.cache_clear()
H = json.load(open(OUT + "/heldout.json")); knots = np.array(H["knots"]); g = np.array(H["g_new_tr"])
key = (0.0, 10.0, 2.0)
b = bench(key); raw = nl.ir(nl.table()[key]); fr = b["front_shift"]
n = 60000
rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
yo = nlt.unwarp(rawp, b["ppm"])[:n]; yn = unwarp2(rawp, b["ppm"])[:n]
m = M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"])
f = np.fft.rfftfreq(n) * FS
m = np.fft.irfft(np.fft.rfft(m, axis=0) * np.interp(f, knots, g)[:, None], n, axis=0)
edges = np.arange(11000, 12801, 300)
for nm, y in (("old unwarp", yo), ("new unwarp", yn)):
    print(nm)
    for t0 in (0.1, 0.4, 0.8, 1.2):
        a = int(t0 * FS); win = 4096; w = np.hanning(win); ff = np.fft.rfftfreq(win) * FS
        Y = np.fft.rfft(y[a:a + win, 0] * w); Mm = np.fft.rfft(m[a:a + win, 0] * w)
        s = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (ff >= lo) & (ff < hi); x = np.sum(Y[k] * np.conj(Mm[k]))
            s.append("%+.3f/%.3f" % (np.angle(x), np.sqrt(np.sum(abs(Y[k]) ** 2) / np.sum(abs(Mm[k]) ** 2))))
        print("  t=%.1f  ph/amp at %s: %s" % (t0, [int((lo + hi) / 2) for lo, hi in zip(edges[:-1], edges[1:])], " ".join(s)))
