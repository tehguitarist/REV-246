"""Scan: per capture, null in 100 Hz-1 kHz over [content start, 30 ms] and [30, 60 ms] from the front frame zero."""
import sys, json, time
from lb_common import *
B = bench_all()
sos = bands(100, 1000)
out = {}
t0 = time.time()
keys = [k for k in B if not B[k]["trimmed"] and (len(sys.argv) < 2 or k[0] in [float(s) for s in sys.argv[1:]])]
for i, k in enumerate(keys):
    r = B[k]
    y, m, _ = load(k, r)
    fr = r["front_shift"]
    z = K0 + fr
    res = {}
    for c in (0, 1):
        yb, mb = sosfiltfilt(sos, y[:, c]), sosfiltfilt(sos, m[:, c])
        for name, (a, b) in {"w30": (0, 0.03), "w60": (0.03, 0.06)}.items():
            ia = max(int(z + a * FS), K0 + 20 if a == 0 else 0)
            ib = int(z + b * FS)
            e = np.sum((yb[ia:ib] - mb[ia:ib]) ** 2)
            res[f"{name}_{c}"] = float(10 * np.log10(e / np.sum(yb[ia:ib] ** 2)))
    out["_".join(map(str, k))] = dict(fs=fr, **res)
    if i % 100 == 0:
        print(i, len(keys), time.time() - t0, flush=True)
json.dump(out, open(os.path.join(OUTD, "scan.json"), "w"))
