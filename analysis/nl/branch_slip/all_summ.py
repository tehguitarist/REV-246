"""Whole sampled N set (every 4th capture, 261): plugin R null per window, before (phase 5, library slip left in) and after
(slip undone for the 165 slipped captures), medians, excess over the noise floor, number <= -50 dB. Untrimmed and trimmed apart."""
import json, glob, os
import numpy as np
S = os.path.dirname(os.path.abspath(__file__)); R = os.path.dirname(os.path.dirname(os.path.dirname(S)))
STEP = 96000 / 26373.6
orig = {os.path.basename(f): json.load(open(f)) for f in glob.glob(R + "/analysis/plugin_check/n_fix/*.json")}
des = {os.path.basename(f): json.load(open(f)) for f in glob.glob(S + "/deslip/n_fix/*.json")}
W = ["0.0-0.05", "0.05-0.15", "0.15-0.3", "0.3-0.6", "0.6-1.0", "1.0-2.0"]
for trimmed in (False, True):
    print("trimmed" if trimmed else "untrimmed")
    for t in W:
        row = {}
        for tag in ("before", "after"):
            v = []
            for n, d in orig.items():
                if d["trimmed"] != trimmed or t not in d["chain3"]["windows"]: continue
                dd = des.get(n, d) if tag == "after" else d
                if t not in dd["chain3"]["windows"]: continue
                w = dd["chain3"]["windows"][t]; v.append((w["null1"], w["null0"], w["null1"] - w["floor1"]))
            row[tag] = np.array(v)
        s = []
        for tag in ("before", "after"):
            a = row[tag]
            s.append("%s N=%d R %.1f (L %.1f) excess %+.1f <=-50: %d" % (tag, len(a), np.median(a[:, 0]), np.median(a[:, 1]), np.median(a[:, 2]), (a[:, 0] <= -50).sum()))
        print(t, " | ".join(s))
