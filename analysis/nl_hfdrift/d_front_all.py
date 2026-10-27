"""Fork D: c_front's estimator on every usable capture (variant 'all'), keeping per-capture num/den (11 knots), tap offsets,
gain, the bench's ppm / o / front_shift and the capture timestamp. usage: d_front_all.py N|all [out.pkl]"""
import os, sys, json, pickle, time
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import c_front as C
from c_front import *
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
times = json.load(open(OUT + "/d_times.json"))
n = sys.argv[1] if len(sys.argv) > 1 else "10"
ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
ks = [k for k in ks if usable(k, common.bench(k))]
if n != "all": ks = ks[::max(1, len(ks) // int(n))][:int(n)]
print(len(ks), "captures", flush=True)
rows = []; t0 = time.time()
for i, key in enumerate(ks):
    b = common.bench(key); h = nl.ir(nl.table()[key])
    row = dict(key=key, ppm=b["ppm"], o=b["o"], fr=b["front_shift"], t=times.get("_".join(str(v) for v in key)))
    for ch in (0, 1):
        row[ch] = C.fit_channel(key, ch, b, h, C.TAPS_USED[ch])
    rows.append(row)
    if i % 50 == 0: print(i, "%.0f s" % (time.time() - t0), flush=True)
pickle.dump(rows, open(OUT + "/" + (sys.argv[2] if len(sys.argv) > 2 else "d_front_smoke.pkl"), "wb"))
