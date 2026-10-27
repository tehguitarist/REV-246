"""Per capture (untrimmed, EARLY>=0.3): baseline offset (pre-tap), nulls in w30 for several bands. Saves scan2.json."""
import json, time
from lb_common import *
B = bench_all()
BANDS = {"20-100": bands(20, 100), "100-1k": bands(100, 1000), "150-1k": bands(150, 1000, 6), "300-1k": bands(300, 1000, 6), "1k-5k": bands(1000, 5000)}
keys = [k for k in B if not B[k]["trimmed"] and k[1] >= 0.3]
out = {}; t0 = time.time()
for i, k in enumerate(keys):
    r = B[k]
    y, m, _ = load(k, r, n=2200)
    z = int(round(K0 + r["front_shift"]))
    res = dict(fs=r["front_shift"])
    for c in (0, 1):
        pk = np.abs(m[z:z + 800, c]).max()
        e = y[:, c] - m[:, c]
        ia = K0 + 20; ib = int(round(z + 0.03 * FS))
        # baseline over pre-tap file samples 20..(z+30)
        pa, pb = K0 + 10, min(z + 30, z + 36)
        res[f"base_y{c}"] = float(y[pa:pb, c].mean() / pk); res[f"base_m{c}"] = float(m[pa:pb, c].mean() / pk)
        res[f"base_e{c}"] = float(e[pa:pb].mean() / pk)
        for nm, sos in BANDS.items():
            yb, mb = sosfiltfilt(sos, y[:, c]), sosfiltfilt(sos, m[:, c])
            res[f"w30_{nm}_{c}"] = float(10 * np.log10(np.sum((yb[ia:ib] - mb[ia:ib]) ** 2) / np.sum(yb[ia:ib] ** 2)))
    out["_".join(map(str, k))] = res
    if i % 100 == 0: print(i, len(keys), time.time() - t0, flush=True)
json.dump(out, open("scan2.json", "w"))
