import json, numpy as np, collections
d = json.load(open("scan.json"))
g = collections.defaultdict(list)
for k, v in d.items():
    s, e, l = k.split("_")
    if float(e) < 0.3: continue
    g[(s, round(v["fs"], 1))].append(v)
print("EARLY>=0.3 only.  SPAN  fs     n   w30 L/R   w60 L/R  (median, 100 Hz-1 kHz)")
for key in sorted(g):
    L = g[key]
    if len(L) < 4: continue
    f = lambda n: np.median([v[n] for v in L])
    print(f"{key[0]:5s} {key[1]:6.1f} {len(L):4d}  {f('w30_0'):6.1f} {f('w30_1'):6.1f}  {f('w60_0'):6.1f} {f('w60_1'):6.1f}")
