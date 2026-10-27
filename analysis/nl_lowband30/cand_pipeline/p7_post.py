"""Is the file's pre-roll DC (first 30 samples) also the DC that explains the 20-100 Hz residual AFTER the step (idx 46..T)?
LS fit of band-limited (cache - LTI) onto the pipeline's unit-DC response k(t), over several index ranges; beta = effective DC (1e-6)."""
import sys, json; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=220, precision=1, suppress=True)
B = butter(4, [20, 100], btype="band", fs=FSD, output="sos")
keys = [tuple(float(v) for v in a.split(",")) for a in sys.argv[2:]]
out = {}
for key in keys:
    r = align(key); x, _ = sf.read(file_path(key), always_2d=True); nf = len(x); cache = nl.ir(key); n = len(cache)
    mp = M.render(key, n, r["o"], r["front_shift"], band=False, proc=True) * np.array(r["gain"])
    e = cache - forward_warp(mp, r["ppm"])
    k, _ = pipeline(np.ones((nf, 2))); k = k[:n]
    eb, kb = sosfiltfilt(B, e, axis=0), sosfiltfilt(B, k, axis=0)
    dc = x[:30].mean(0)
    row = {}
    for nm, a, b in (("post 46-400", 46, 400), ("post 46-1500", 46, 1500), ("post 400-3000", 400, 3000), ("post 3000-16000", 3000, 16000)):
        if b > n: continue
        row[nm] = [float((eb[a:b, c] @ kb[a:b, c]) / (kb[a:b, c] @ kb[a:b, c])) * 1e6 for c in (0, 1)]
    row["pre dc"] = (dc * 1e6).tolist()
    out["_".join(map(str, key))] = row
    print(key, " | ".join(f"{nm} {np.array(v)}" for nm, v in row.items()), flush=True)
json.dump(out, open(sys.argv[1], "w"))
