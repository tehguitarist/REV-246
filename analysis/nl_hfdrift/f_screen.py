"""Bench the HF correction at one strength (env NL_HF_ALPHA) on a subset: 100 captures whose 10-12 kHz band got worse by > 0.5 dB at alpha 1, plus 50 evenly sampled."""
import os, sys, json, glob
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import nl_model as M
import defects
if __name__ == "__main__":
    new, old = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7_lb2_hf", "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7_lb2"
    keys = sorted(k for k in M.nl.table() if k not in defects.NL_DEFECTS)
    worse = []
    for k in keys:
        n = "_".join(str(v) for v in k) + ".json"
        a, b = json.load(open(f"{new}/{n}")), json.load(open(f"{old}/{n}"))
        if "10000-12000" in a["bands"] and max(np.array(a["bands"]["10000-12000"]) - np.array(b["bands"]["10000-12000"])) > 0.5: worse.append(k)
    rng = np.random.default_rng(3); sel = set(rng.choice(len(worse), 100, replace=False).tolist()); sel = {worse[i] for i in sel} | set(keys[::21][:50])
    sel = sorted(sel); print(len(sel), "captures; folder", M.BENCH, flush=True)
    json.dump([list(k) for k in sel], open("/Users/leighpierce/Git/246/analysis/nl_hfdrift/f_screen_keys.json", "w"))
    with ProcessPoolExecutor(4) as ex:
        for _ in ex.map(M.bench_file, sel, chunksize=2): pass
