"""M / L (sml_fullband --stage ml) with and without ML_LOWBAND=1 on the common groups: median model null per band and window (R / L channel pooled), and the repeat level."""
import json, glob, os, sys
import numpy as np
R = "/Users/leighpierce/Git/246/analysis/fullband/"
EDGES = [20, 40, 80, 160, 270, 400, 800, 1600, 3200]; BN = [f"{a}-{b}" for a, b in zip(EDGES[:-1], EDGES[1:])]
new = {os.path.basename(f): f for f in glob.glob(R + "ml_lowband/*.json")}
def load(path): return json.load(open(path))
for prog in "ML":
    names = sorted(n for n in new if n.startswith(prog + "_") and os.path.exists(R + "ml/" + n))
    print(f"\n=== {prog}: {len(names)} common groups.  median null dB (model | repeat) base -> lowband")
    for w in ("cluster_gap", "first_pass", "recirc1"):
        print(f"  window {w}")
        for b in BN:
            v = {"base": [], "new": [], "rep": []}
            for n in names:
                for tag, root in (("base", R + "ml/"), ("new", R + "ml_lowband/")):
                    for row in load(root + n):
                        x = row["bands"].get(b, {}).get(w)
                        if not x: continue
                        for ch in (0, 1):
                            if x[ch] and x[ch][0] > 0 and x[ch][1] > 0:
                                v[tag].append(10 * np.log10(x[ch][1] / x[ch][0]))
                                if tag == "base" and x[ch][2] > 0: v["rep"].append(10 * np.log10(x[ch][2] / x[ch][0]))
            if v["base"]:
                print(f"    {b:>9s}: {np.median(v['base']):6.1f} -> {np.median(v['new']):6.1f}   (repeat {np.median(v['rep']):6.1f}; n {len(v['base'])}; worse by >0.5 dB median-shift {np.median(v['new'])-np.median(v['base']):+.2f})")
