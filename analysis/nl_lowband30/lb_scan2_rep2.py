import json, numpy as np, collections
d = json.load(open("scan2.json"))
# SPAN 0.5: which settings are in which session
for sp in ("0.5", "0.3"):
    ss = collections.defaultdict(list)
    for k, v in d.items():
        s, e, l = k.split("_")
        if s != sp: continue
        ss["fs0" if abs(v["fs"]) < 0.5 else "other"].append((float(e), float(l), 0.5 * (v["w30_20-100_0"] + v["w30_20-100_1"]), 0.5 * (v["w30_100-1k_0"] + v["w30_100-1k_1"])))
    for sn, L in ss.items():
        E = collections.Counter(e for e, *_ in L); Lt = sorted({l for _, l, *_ in L})
        print(f"SPAN {sp} {sn}: n={len(L)} EARLY counts {dict(sorted(E.items()))} LATE values {Lt}")
# by LATE and by EARLY, SPAN 0.5 and 0.3
for sp in ("0.5", "0.3", "0.2", "0.4"):
    print("SPAN", sp, " 20-100/100-1k median null by EARLY:")
    by = collections.defaultdict(list)
    for k, v in d.items():
        s, e, l = k.split("_")
        if s == sp: by[(float(e), "fs0" if abs(v["fs"]) < 0.5 else "oth")].append((0.5 * (v["w30_20-100_0"] + v["w30_20-100_1"]), 0.5 * (v["w30_100-1k_0"] + v["w30_100-1k_1"])))
    print("  " + "  ".join(f"{e}/{c}: {np.median([a for a,_ in L]):.0f}/{np.median([b for _,b in L]):.0f}(n{len(L)})" for (e, c), L in sorted(by.items())))
