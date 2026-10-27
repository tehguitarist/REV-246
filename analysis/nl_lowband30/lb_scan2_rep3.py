import json, numpy as np, collections
d = json.load(open("scan2.json"))
def tab(filt, title):
    g = collections.defaultdict(list)
    for k, v in d.items():
        s, e, l = (float(x) for x in k.split("_"))
        if filt(s, e, l): g[(s, round(v["fs"]) if abs(v["fs"]) > 0.5 else 0)].append(v)
    print(title)
    for key in sorted(g):
        L = g[key]
        if len(L) < 3: continue
        f = lambda nm: np.median([0.5 * (v[f"w30_{nm}_0"] + v[f"w30_{nm}_1"]) for v in L])
        print(f"  SPAN {key[0]} fs {key[1]:3d} n {len(L):3d} | 20-100 {f('20-100'):6.1f} 100-1k {f('100-1k'):6.1f} 150-1k {f('150-1k'):6.1f} 300-1k {f('300-1k'):6.1f} 1k-5k {f('1k-5k'):6.1f}")
tab(lambda s, e, l: l == 0.0, "LATE 0 (no late section), EARLY>=0.3")
tab(lambda s, e, l: l == 0.0 and e in (3.0, 10.0), "LATE 0, EARLY 3 or 10")
tab(lambda s, e, l: l in (0.2, 0.4), "LATE 0.2-0.4")
tab(lambda s, e, l: l >= 1.4, "LATE 1.4-2.0")
tab(lambda s, e, l: l >= 1.4 and e in (3.0, 10.0), "LATE 1.4-2.0, EARLY 3 or 10")
