import json, numpy as np, collections
O = json.load(open("oracle.json")); P = json.load(open("oracle2.json"))
g = collections.defaultdict(list)
for k, p in P.items():
    v = O[k]; s = float(k.split("_")[0]); fs = p["fs"]
    g[(s, round(fs) if abs(fs) > 0.5 else 0)].append((v, p))
m2 = lambda v, t, nm: 0.5 * (v[f"{t}_{nm}_0"] + v[f"{t}_{nm}_1"])
print("first-30ms median null [dB]: raw | measured step (2 numbers/capture, pre-roll) | predicted step (leave-one-EARLY-out, 3 numbers/group)")
print("SPAN fs   n |      20-100              |      100-1k              |      150-1k")
allr = collections.defaultdict(list)
for key in sorted(g):
    L = g[key]
    if len(L) < 5: continue
    row = f"{key[0]:3.1f} {key[1]:3d} {len(L):3d} |"
    for nm in ("20-100", "100-1k", "150-1k"):
        a = [m2(v, "raw", nm) for v, p in L]; b = [m2(v, "cor", nm) for v, p in L]; c = [m2(p, "pred", nm) for v, p in L]
        row += f" {np.median(a):6.1f} {np.median(b):6.1f} {np.median(c):6.1f} |"
        for t, x in (("raw", a), ("cor", b), ("pred", c)): allr[(nm, t, key[1] == 0)] += x
    print(row)
for fs0 in (True, False):
    print("fs0 sessions" if fs0 else "other sessions", " ".join(f"{nm}: {np.median(allr[(nm,'raw',fs0)]):.1f}/{np.median(allr[(nm,'cor',fs0)]):.1f}/{np.median(allr[(nm,'pred',fs0)]):.1f}" for nm in ("20-100", "100-1k", "150-1k")))
# fraction of captures made worse by the predicted correction
w = [m2(p, "pred", "100-1k") - m2(v, "raw", "100-1k") for v, p in sum(g.values(), [])]
print("100-1k: predicted correction changes null by median %.2f dB; worse by >0.5 dB in %d of %d captures; better by >0.5: %d" % (np.median(w), sum(x > 0.5 for x in w), len(w), sum(x < -0.5 for x in w)))
