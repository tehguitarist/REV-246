"""Held-out-by-setting version: the step c(SPAN-session, LATE) predicted from the OTHER EARLY values of the same SPAN / session
(leave-one-EARLY-out; L: constant per group, R: a + b*LATE per group; 3 numbers per group) instead of measured on the capture. Saves oracle2.json."""
import json, time, collections
from lb_common import *
B = bench_all()
O = json.load(open("oracle.json"))
BD = {"20-100": bands(20, 100), "100-1k": bands(100, 1000), "150-1k": bands(150, 1000, 6), "300-1k": bands(300, 1000, 6)}
grp = collections.defaultdict(list)
for k, v in O.items():
    s, e, l = (float(x) for x in k.split("_"))
    grp[(s, round(v["fs"]) if abs(v["fs"]) > 0.5 else 0)].append((e, l, v["c0"], v["c1"]))
def predict(ses, e0, l0):
    A = np.array([(l, c0, c1) for e, l, c0, c1 in grp[ses] if e != e0])
    if len(A) < 8: return None
    c0 = np.median(A[:, 1])
    p = np.polyfit(A[:, 0], A[:, 2], 1)
    return c0, p[1] + p[0] * l0
n = 2200; t = np.arange(n)
step = 0.5 * (1 - np.tanh((t - 45.5) / 1.5))
out = {}; t0 = time.time()
keys = [k for k in B if "_".join(map(str, k)) in O]
for i, k in enumerate(keys):
    r = B[k]; v0 = O["_".join(map(str, k))]
    ses = (k[0], round(r["front_shift"]) if abs(r["front_shift"]) > 0.5 else 0)
    pr = predict(ses, k[1], k[2])
    if pr is None: continue
    y, m, _ = load(k, r, n=n)
    z = int(round(K0 + r["front_shift"]))
    res = dict(fs=r["front_shift"])
    for c in (0, 1):
        pk = np.abs(m[z:z + 800, c]).max()
        yc = y[:, c] - pr[c] * pk * step
        ia, ib = K0 + 20, int(round(z + 0.03 * FS))
        for nm, sos in BD.items():
            mb = sosfiltfilt(sos, m[:, c]); yb = sosfiltfilt(sos, yc)
            res[f"pred_{nm}_{c}"] = float(10 * np.log10(np.sum((yb[ia:ib] - mb[ia:ib]) ** 2) / np.sum(yb[ia:ib] ** 2)))
    out["_".join(map(str, k))] = res
    if i % 100 == 0: print(i, len(keys), time.time() - t0, flush=True)
json.dump(out, open("oracle2.json", "w"))
