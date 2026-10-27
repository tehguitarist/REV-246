import json, numpy as np, collections
d = json.load(open("scan.json"))
g = collections.defaultdict(list)
for k, v in d.items():
    g[round(v["fs"], 1)].append((k, v))
print("front shift | n | spans | w30 L/R med | w60 L/R med  (100 Hz-1 kHz null dB, untrimmed)")
for fs in sorted(g, key=lambda s: -len(g[s])):
    L = g[fs]
    if len(L) < 5: continue
    sp = collections.Counter(k.split("_")[0] for k, _ in L)
    f = lambda n: np.median([v[n] for _, v in L])
    print(f"{fs:7.2f} {len(L):4d} {dict(sp)}  {f('w30_0'):6.1f} {f('w30_1'):6.1f} | {f('w60_0'):6.1f} {f('w60_1'):6.1f}")
