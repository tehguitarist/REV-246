"""Run the bench (warp2 + refit kernel) at one strength on a subset: the captures whose kernel-vs-none loss exceeded 0.3 dB at
alpha 0.7 in any of the first four windows, plus every 10th capture. Env NL_WARP2=1 NL_BK=1 NL_BK_FILE=... NL_BK_ALPHA=a set by the caller."""
import os, sys, json, glob
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
from concurrent.futures import ProcessPoolExecutor
import nl_model as M
import defects
if __name__ == "__main__":
    new, base = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7", "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
    keys = sorted(k for k in M.nl.table() if k not in defects.NL_DEFECTS)
    sel = set(keys[::10])
    for k in keys:
        n = "_".join(str(v) for v in k) + ".json"
        a, b = json.load(open(f"{new}/{n}")), json.load(open(f"{base}/{n}"))
        for w in ("0.0-0.05", "0.05-0.15", "0.15-0.3", "0.3-0.6"):
            if w in a["windows"] and w in b["windows"] and max(a["windows"][w][f"null{c}"] - b["windows"][w][f"null{c}"] for c in (0, 1)) > 0.3:
                sel.add(k)
    sel = sorted(sel); print(len(sel), "captures; bench folder", M.BENCH, flush=True)
    with ProcessPoolExecutor(4) as ex:
        for _ in ex.map(M.bench_file, sel, chunksize=2):
            pass
