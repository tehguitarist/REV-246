import json, numpy as np
d = json.load(open("scan2.json")); B = {}
import glob, os
for f in glob.glob("../nl/bench/*.json"):
    B[os.path.basename(f)[:-5]] = json.load(open(f))
for sp, e in (("0.4", "10.0"), ("0.5", "10.0"), ("0.5", "3.0"), ("0.3", "10.0")):
    print("SPAN", sp, "EARLY", e, " late: o  fs  null100-1k(L,R)  null20-100(L,R)  base_e(L,R)e-4")
    for l in [f"{x/5:.1f}" for x in range(0, 11)]:
        k = f"{sp}_{e}_{l}"
        if k in d:
            v = d[k]; b = B[k]
            print(f"  {l}: o {b['o']:8.2f} fs {v['fs']:6.2f}  {v['w30_100-1k_0']:6.1f} {v['w30_100-1k_1']:6.1f}   {v['w30_20-100_0']:6.1f} {v['w30_20-100_1']:6.1f}   {v['base_e0']*1e4:5.1f} {v['base_e1']*1e4:5.1f}")
