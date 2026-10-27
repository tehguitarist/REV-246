"""Compare two bench folders per capture, per channel-window and per band (strict rule: no channel-window worse than +0.5 dB).
usage: compare_bench.py NEW_DIR OLD_DIR"""
import sys, json, os, glob, numpy as np
new, old = sys.argv[1], sys.argv[2]
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
import defects
keys = [os.path.basename(f)[:-5] for f in sorted(glob.glob(new + "/*.json"))]
keys = [k for k in keys if tuple(float(v) for v in k.split("_")) not in defects.NL_DEFECTS and os.path.exists(f"{old}/{k}.json")]
print(len(keys), "captures common to both")
W = ["0.0-0.05", "0.05-0.15", "0.15-0.3", "0.3-0.6", "0.6-1.0", "1.0-2.0"]
worst = {}
for w in W:
    d, nn, oo, fl = [], [], [], []
    for k in keys:
        a, b = json.load(open(f"{new}/{k}.json")), json.load(open(f"{old}/{k}.json"))
        if w in a["windows"] and w in b["windows"]:
            for c in (0, 1):
                x, y = a["windows"][w][f"null{c}"], b["windows"][w][f"null{c}"]
                d.append((x - y, k, c)); nn.append(x); oo.append(y); fl.append(a["windows"][w][f"floor{c}"])
    dd = np.array([v[0] for v in d]); ex = np.array(nn) - np.array(fl)
    wrs = max(d)
    print(f"{w:10s} n {len(d)//2:4d}  null L/R median new {np.median(nn):6.1f}  old {np.median(oo):6.1f} | diff median {np.median(dd):+.2f} dB, worse than +0.5: {int((dd > 0.5).sum())}, worst {wrs[0]:+.2f} ({wrs[1]} ch{wrs[2]}) | excess over floor median new {np.median(np.maximum(ex, 0)):.1f} old {np.median(np.maximum(np.array(oo) - np.array(fl), 0)):.1f}")
print("per band, 0.05-0.6 s (median L / R, new | old; worse than +0.5 dB count, worst):")
for bd in ["20-100", "100-1000", "1000-5000", "5000-10000", "10000-12000"]:
    v = [(json.load(open(f"{new}/{k}.json"))["bands"].get(bd), json.load(open(f"{old}/{k}.json"))["bands"].get(bd)) for k in keys]
    v = [(a, b) for a, b in v if a and b]
    A, B = np.array([a for a, b in v]), np.array([b for a, b in v])
    D = (A - B).ravel()
    print(f"  {bd:12s} new {np.median(A[:,0]):6.1f} / {np.median(A[:,1]):6.1f} | old {np.median(B[:,0]):6.1f} / {np.median(B[:,1]):6.1f} | diff median {np.median(D):+.2f}, >+0.5: {int((D>0.5).sum())}/{len(D)}, worst {D.max():+.2f}")
ppm = [(json.load(open(f"{new}/{k}.json"))["ppm"], json.load(open(f"{old}/{k}.json"))["ppm"], json.load(open(f"{new}/{k}.json"))["o"], json.load(open(f"{old}/{k}.json"))["o"]) for k in keys]
p = np.array(ppm); print("ppm new-old: mean %.3f sd %.3f; o new-old: median %.4f, max |.| %.3f" % ((p[:, 0] - p[:, 1]).mean(), (p[:, 0] - p[:, 1]).std(), np.median(p[:, 2] - p[:, 3]), np.abs(p[:, 2] - p[:, 3]).max()))
