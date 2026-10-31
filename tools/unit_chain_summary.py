"""Condensed view of analysis/unit_chain/fullband13.json and joint_run/fullband_joint.json: per band, mean over the held-out
windows 50-150 ... 600-1000 ms of (variant - plain) null R in dB, per unit set (CLA; the five Altiverb presets)."""
import json, os, sys
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis", "unit_chain")
files = [os.path.join(OUT, "fullband13.json"), os.path.join(OUT, "joint_run", "fullband_joint.json")]
D = {}
for p in files:
    if os.path.exists(p):
        r = json.load(open(p))
        for u, v in r["irs"].items():
            D.setdefault(u, {}).update({k: x for k, x in v.items() if k != "rt"})
bands = json.load(open(files[0]))["bands"][:11]
W = ["50-150", "150-300", "300-600", "600-1000"]
want = sys.argv[1].split(",")
for sn, us in (("CLA", ["CLA"]), ("Altiverb (5 presets, median)", [u for u in D if u.startswith("Altiverb")])):
    print(f"## {sn}: mean over 50-1000 ms windows of null R minus plain (dB; negative = better)")
    print(f"{'band':12s}" + "".join(f"{w:>9s}" for w in want))
    for b in bands:
        row = []
        for vn in want:
            ds = []
            for u in us:
                if vn in D[u]:
                    ds.append(np.mean([D[u][vn]["null"][b][w][1] - D[u]["plain"]["null"][b][w][1] for w in W if D[u][vn]["null"][b][w]]))
            row.append(np.median(ds) if ds else np.nan)
        print(f"{b:12s}" + "".join(f"{x:+9.1f}" for x in row))
    print()
