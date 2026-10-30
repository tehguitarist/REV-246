"""Grid bench (NL_OGRID=1, o = 1 + n x 26373.6/96000) against the free-o bench: o agreement and nulls per window, untrimmed."""
import json, os, numpy as np
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = R + "/bench_warp2_bk_kf_a0.7_lb2_hf_g5/"; B = R + "/bench_warp2_bk_kf_a0.7_lb2_og_hf_g5/"
St = 96000 / 26373.6
W = ["0.0-0.05", "0.05-0.15", "0.15-0.3", "0.3-0.6", "0.6-1.0", "1.0-2.0"]
rows = {}; nd = 0; ns = 0
for f in sorted(os.listdir(B)):
    a, b = json.load(open(A + f)), json.load(open(B + f))
    if a["trimmed"]: continue
    if b["o"] < -150 or a["o"] < -150: continue            # SPAN 0 / LATE 0, no branch to fit
    ns += 1
    na = round((a["o"] - 1) * St); nb = round((b["o"] - 1) * St); nd += na != nb
    for t in W:
        if t in a["windows"] and t in b["windows"]:
            rows.setdefault(t, []).append((a["windows"][t]["null1"], b["windows"][t]["null1"], b["windows"][t]["null1"] - b["windows"][t]["floor1"]))
print("untrimmed captures with a branch:", ns, " n differs between free and grid:", nd)
for t in W:
    v = np.array(rows[t]); print(t, len(v), "free %.1f  grid %.1f  excess(grid) %+.1f  <=-50: %d -> %d  worse by >1 dB: %d" % (np.median(v[:,0]), np.median(v[:,1]), np.median(v[:,2]), (v[:,0]<=-50).sum(), (v[:,1]<=-50).sum(), (v[:,1]-v[:,0]>1).sum()))
