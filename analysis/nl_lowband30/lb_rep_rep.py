import json, sys, numpy as np
d = json.load(open(sys.argv[1]))
bs = ["20-100", "100-1000", "1000-5000", "5000-10000", "10000-12000"]
print("EARLY   n |", " ".join(f"{b:>11s}" for b in bs), "   (median pair null, L/R mean, raw; minus 3 dB for one capture is -3)")
pool = {b: [] for b in bs}; pool2 = {b: [] for b in bs}
for k, v in d.items():
    if k.endswith("_pairs"): continue
    n = v[bs[0]][2]
    print(f"{k:>6s} {n:3d} |", " ".join(f"{0.5*(v[b][0]+v[b][1]):11.1f}" for b in bs))
    for b in bs:
        P = d[k + "_pairs"][b]; pool[b] += P[0] + P[1]
        if float(k) >= 0.7: pool2[b] += P[0] + P[1]
print("pooled median, all EARLY      :", " ".join(f"{np.median(pool[b]):11.1f}" for b in bs))
print("pooled median, EARLY >= 0.7   :", " ".join(f"{np.median(pool2[b]):11.1f}" for b in bs))
print("(doc 11.1: 100-1k -66.5, 1-5k -60.6, 5-10k -59.4, 10-12k -58.9)")
