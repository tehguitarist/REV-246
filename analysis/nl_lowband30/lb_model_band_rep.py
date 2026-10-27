import json, sys, collections, numpy as np
d = json.load(open(sys.argv[1])); rp = json.load(open(sys.argv[2]))
bs = ["20-100", "100-1000", "1000-5000", "5000-10000", "10000-12000"]
g = collections.defaultdict(list)
for k, v in d.items(): g[float(k.split("_")[1])].append(v)
print("SPAN 0.5, window z+0.05..0.45 s.  Model null (fixed bench alignment | alignment re-fitted per band)  vs  repeat level (pairs, raw), median L/R mean")
print("EARLY  n |" + "".join(f" {b:>11s}" for b in bs))
tot = collections.defaultdict(list)
for e in sorted(g, reverse=True):
    L = g[e]
    f = lambda b, j: np.median([0.5 * (v[b][0][j] + v[b][1][j]) for v in L])
    print(f"{e:5.1f} {len(L):3d} fixed |" + "".join(f" {f(b,0):11.1f}" for b in bs))
    print(f"          fit   |" + "".join(f" {f(b,1):11.1f}" for b in bs))
    if str(e) in rp: print(f"          repeat|" + "".join(f" {0.5*(rp[str(e)][b][0]+rp[str(e)][b][1]):11.1f}" for b in bs))
    for b in bs:
        for v in L: tot[(b, "fixed")].append(0.5 * (v[b][0][0] + v[b][1][0])); tot[(b, "fit")].append(0.5 * (v[b][0][1] + v[b][1][1]))
print("pooled model fixed:", " ".join(f"{np.median(tot[(b,'fixed')]):7.1f}" for b in bs))
print("pooled model fit  :", " ".join(f"{np.median(tot[(b,'fit')]):7.1f}" for b in bs))
