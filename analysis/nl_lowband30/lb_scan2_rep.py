import json, numpy as np, collections
d = json.load(open("scan2.json"))
g = collections.defaultdict(list)
for k, v in d.items():
    s, e, l = k.split("_")
    g[(s, round(v["fs"] / 0.5) * 0.5 if abs(v["fs"]) < 0.5 else round(v["fs"], 0))].append((float(e), v))
bn = ["20-100", "100-1k", "150-1k", "300-1k", "1k-5k"]
print("SPAN  fs    n |  median w30 null dB (mean of L,R) per band:", bn, "| pre-tap baseline (1e-4 of 410-peak): data  model  err(med)  |err| ")
for key in sorted(g):
    L = [v for _, v in g[key]]
    if len(L) < 6: continue
    f = lambda nm: np.median([0.5 * (v[f"w30_{nm}_0"] + v[f"w30_{nm}_1"]) for v in L])
    by = np.median([abs(v["base_y0"]) for v in L]) * 1e4; bm = np.median([abs(v["base_m0"]) for v in L]) * 1e4
    be = np.median([v["base_e0"] for v in L]) * 1e4; ba = np.median([abs(v["base_e0"]) for v in L]) * 1e4
    print(f"{key[0]:4s} {key[1]:5.1f} {len(L):3d} | " + " ".join(f"{f(b):6.1f}" for b in bn) + f" | {by:5.1f} {bm:5.1f} {be:5.1f} {ba:5.1f}")
