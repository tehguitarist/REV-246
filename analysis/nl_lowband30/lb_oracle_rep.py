import json, numpy as np, collections
d = json.load(open("oracle.json"))
g = collections.defaultdict(list)
for k, v in d.items():
    s, e, l = (float(x) for x in k.split("_"))
    g[(s, round(v["fs"]) if abs(v["fs"]) > 0.5 else 0)].append((e, l, v))
def mean2(v, tag, nm): return 0.5 * (v[f"{tag}_{nm}_0"] + v[f"{tag}_{nm}_1"])
print("first-30ms null, median over captures (EARLY>=0.3, untrimmed): raw -> pre-roll-step removed  [dB]")
print("SPAN  fs   n |  20-100          | 100-1k           | 150-1k           | 300-1k           | median |c|(1e-4 re peak) L / R")
tot = collections.defaultdict(list)
for key in sorted(g):
    L = g[key]
    if len(L) < 5: continue
    row = f"{key[0]:4.1f} {key[1]:3d} {len(L):3d} |"
    for nm in ("20-100", "100-1k", "150-1k", "300-1k"):
        a = np.median([mean2(v, "raw", nm) for _, _, v in L]); b = np.median([mean2(v, "cor", nm) for _, _, v in L])
        row += f" {a:6.1f} -> {b:6.1f} |"
    row += f" {np.median([abs(v['c0']) for *_, v in L])*1e4:4.1f} / {np.median([abs(v['c1']) for *_, v in L])*1e4:4.1f}"
    print(row)
    for _, _, v in L:
        for nm in ("20-100", "100-1k"):
            tot[("fs0" if key[1] == 0 else "other", nm, "raw")].append(mean2(v, "raw", nm)); tot[("fs0" if key[1] == 0 else "other", nm, "cor")].append(mean2(v, "cor", nm))
for s in ("fs0", "other"):
    print(s, " ".join(f"{nm}: {np.median(tot[(s,nm,'raw')]):.1f} -> {np.median(tot[(s,nm,'cor')]):.1f}" for nm in ("20-100", "100-1k")))
