import json, numpy as np, collections
d = json.load(open("scan.json"))
g = collections.defaultdict(list)
for k, v in d.items():
    s, e, l = k.split("_")
    cl = "fs0" if abs(v["fs"]) < 0.2 else "other"
    g[(float(e), cl)].append(0.5 * (v["w30_0"] + v["w30_1"]))
print("EARLY   fs0: n med | other: n med")
for e in sorted({k[0] for k in g}):
    a, b = g.get((e, "fs0"), []), g.get((e, "other"), [])
    print(f"{e:5.1f}  {len(a):4d} {np.median(a) if a else float('nan'):6.1f} | {len(b):4d} {np.median(b) if b else float('nan'):6.1f}")
