"""Time course of e = y - m (100 Hz-1 kHz and full band), z-relative, normalised by the model's peak, per (SPAN, session) subgroup."""
import sys, collections, pickle
from lb_common import *
B = bench_all()
sos = bands(100, 1000)
SUB = [("0.0", 19.6), ("0.1", 16.4), ("0.2", 22.9), ("0.3", -0.1), ("0.4", 19.6), ("0.5", -0.1), ("0.5", 22.9)]
rng = np.random.default_rng(2)
OUT = {}
for sp, fsc in SUB:
    ks = sorted(k for k, r in B.items() if not r["trimmed"] and k[1] >= 0.3 and k[0] == float(sp) and abs(r["front_shift"] - fsc) < 0.2)
    ks = [ks[i] for i in sorted(rng.choice(len(ks), min(24, len(ks)), replace=False))]
    for k in ks:
        r = B[k]
        y, m, _ = load(k, r, n=2200)
        z = int(round(K0 + r["front_shift"]))
        rec = {}
        for c in (0, 1):
            yb, mb = sosfiltfilt(sos, y[:, c]), sosfiltfilt(sos, m[:, c])
            sc = np.abs(m[z:z + 800, c]).max()
            rec[c] = dict(e=(y[:, c] - m[:, c]) / sc, eb=(yb - mb) / sc, mb=mb / sc, m=m[:, c] / sc, y=y[:, c] / sc, z=z)
        OUT[(sp, fsc, k)] = rec
    print(sp, fsc, len(ks), flush=True)
pickle.dump(OUT, open("time.pkl", "wb"))
