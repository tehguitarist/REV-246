import json, glob, os, numpy as np
rows = []
for f in sorted(glob.glob("/Users/leighpierce/Git/246/analysis/nl/bench/*.json")):
    b = json.load(open(f)); k = os.path.basename(f)[:-5].split("_")
    rows.append((float(k[0]), float(k[1]), float(k[2]), b["o"], b["ppm"], b["front_shift"], b["trimmed"]))
a = np.array([r[:6] for r in rows]); tr = np.array([r[6] for r in rows])
print(len(rows), "trimmed", tr.sum(), "ppm mean %.2f sd %.2f" % (a[:, 4].mean(), a[:, 4].std()))
for r in rows:
    if r[6] and r[0] >= 0.3 and r[2] == 0.0: print(r)
