"""Fork H: build the scoring segments (capture y and uncorrected model m, 0.05-0.6 s window, new warp, NL_BK=1, lb2 default) for a
stratified seeded sample, once; saved outside the repo. usage: h_seg.py N_PER_CLUSTER OUT.pkl"""
import os, sys, json, random, pickle
os.environ["NL_BK"] = "1"; os.environ.pop("NL_HF", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import sosfiltfilt
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7_lb2"
import defects
npc = int(sys.argv[1]); outp = sys.argv[2]
times = json.load(open(OUT + "/d_times.json"))
keys = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and "_".join(str(v) for v in k) in times
              and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
def cluster(k):
    b = common.bench(k)
    if b["trimmed"]: return "trimmed"
    return "S0.1" if k[0] == 0.1 else "S0.3" if k[0] == 0.3 else "rest"
cl = {}
for k in keys: cl.setdefault(cluster(k), []).append(k)
print({c: len(v) for c, v in cl.items()}, flush=True)
random.seed(21)
sample = []
for c, v in cl.items(): sample += [(c, k) for k in random.sample(v, min(npc, len(v)))]
out = []
for i, (c, key) in enumerate(sample):
    b = common.bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
    a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if b2 - a < int(0.05 * FS): continue
    n = b2 + 200
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = nlt.unwarp(rawp, b["ppm"])[a:b2]
    m = (M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"]))[a:b2]
    out.append(dict(key=key, cluster=c, t=times["_".join(str(v) for v in key)], y=y.astype(np.float32), m=m.astype(np.float32)))
    if i % 25 == 0: print(i, len(sample), flush=True)
pickle.dump(out, open(outp, "wb")); print("saved", len(out))
