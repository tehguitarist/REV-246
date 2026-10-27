"""Oracle for the file-start step: c_ch = mean(e[8:40]) (pre-roll, outside every scored window, 1 number per channel and capture); the data's pre-roll
baseline step c*u(t < 45.5) is removed from y (smoothed edge) and the first-30 ms nulls are recomputed. Saves oracle.json."""
import json, time
from lb_common import *
B = bench_all()
BD = {"20-100": bands(20, 100), "100-1k": bands(100, 1000), "150-1k": bands(150, 1000, 6), "300-1k": bands(300, 1000, 6)}
keys = [k for k in B if not B[k]["trimmed"] and k[1] >= 0.3]
n = 2200; t = np.arange(n)
step = 0.5 * (1 - np.tanh((t - 45.5) / 1.5))                  # 1 before the file start (-4 samples), 0 after
out = {}; t0 = time.time()
for i, k in enumerate(keys):
    r = B[k]
    y, m, _ = load(k, r, n=n)
    z = int(round(K0 + r["front_shift"]))
    res = dict(fs=r["front_shift"])
    for c in (0, 1):
        e = y[:, c] - m[:, c]
        cc = e[8:40].mean()
        yc = y[:, c] - cc * step
        res[f"c{c}"] = float(cc / np.abs(m[z:z + 800, c]).max())
        ia, ib = K0 + 20, int(round(z + 0.03 * FS))
        for nm, sos in BD.items():
            mb = sosfiltfilt(sos, m[:, c])
            for tag, yy in (("raw", y[:, c]), ("cor", yc)):
                yb = sosfiltfilt(sos, yy)
                res[f"{tag}_{nm}_{c}"] = float(10 * np.log10(np.sum((yb[ia:ib] - mb[ia:ib]) ** 2) / np.sum(yb[ia:ib] ** 2)))
    out["_".join(map(str, k))] = res
    if i % 100 == 0: print(i, len(keys), time.time() - t0, flush=True)
json.dump(out, open("oracle.json", "w"))
