import json, numpy as np, collections
d = json.load(open("scan.json"))
g = collections.defaultdict(list)
for k, v in d.items():
    g[(k.split("_")[0], round(v["fs"], 1))].append(v)
print("SPAN  fs     n   w30 L/R   w60 L/R")
for key in sorted(g, key=lambda s: (s[0], -len(g[s]))):
    L = g[key]
    if len(L) < 5: continue
    f = lambda n: np.median([v[n] for v in L])
    print(f"{key[0]:5s} {key[1]:6.1f} {len(L):4d}  {f('w30_0'):6.1f} {f('w30_1'):6.1f}  {f('w60_0'):6.1f} {f('w60_1'):6.1f}")
