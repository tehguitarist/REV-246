"""Compare the N bench with the low-band switch (NL_LOWBAND 1 / 2) against the committed default (bench_warp2_bk_kf_a0.7): per window and per band medians,
per-capture changes (strict rule: channel-windows worse than +0.5 dB), and the split into the 51 captures the table was estimated from (SPAN >= 0.2, LATE <= 0.2, EARLY >= 0.3)
and the rest."""
import json, os, sys, glob
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
ROOT = "/Users/leighpierce/Git/246/analysis/nl"
DIRS = {"base": "bench_warp2_bk_kf_a0.7", "lb1": "bench_warp2_bk_kf_a0.7_lb1", "lb2": "bench_warp2_bk_kf_a0.7_lb2"}
import defects
def load(d):
    out = {}
    for f in glob.glob(os.path.join(ROOT, d, "*.json")):
        k = tuple(float(v) for v in os.path.basename(f)[:-5].split("_"))
        if k not in defects.NL_DEFECTS: out[k] = json.load(open(f))
    return out
R = {n: load(d) for n, d in DIRS.items() if os.path.isdir(os.path.join(ROOT, d))}
keys = sorted(set.intersection(*[set(r) for r in R.values()]))
print({n: len(r) for n, r in R.items()}, "common", len(keys))
est = lambda k: k[0] >= 0.2 and k[2] <= 0.2 and k[1] >= 0.3 and not R["base"][k].get("trimmed")
WIN = ["0.0-0.05", "0.05-0.15", "0.15-0.3", "0.3-0.6", "0.6-1.0", "1.0-2.0"]
BAND = ["20-100", "100-1000", "1000-5000", "5000-10000", "10000-12000"]
def table(sel, label):
    print(f"\n--- {label}: {sum(1 for k in keys if sel(k))} captures")
    for n in R:
        row = []
        for w in WIN:
            v = [R[n][k]["windows"][w][f"null{c}"] for k in keys if sel(k) and w in R[n][k]["windows"] for c in (0, 1)]
            row.append(f"{np.median(v):6.1f}" if v else "   -  ")
        print(f"  {n:5s} windows (null, 100 Hz-10 kHz) {' '.join(WIN)}:  " + " ".join(row))
    for n in R:
        row = []
        for b in BAND:
            v = [R[n][k]["bands"][b][c] for k in keys if sel(k) and b in R[n][k].get("bands", {}) for c in (0, 1)]
            row.append(f"{np.median(v):6.1f}" if v else "   -  ")
        print(f"  {n:5s} bands 0.05-0.6 s ({' | '.join(BAND)}):  " + " ".join(row))
    for n in R:
        if n == "base": continue
        out = []
        for w in WIN[:5]:
            d = [R[n][k]["windows"][w][f"null{c}"] - R["base"][k]["windows"][w][f"null{c}"] for k in keys if sel(k) and w in R[n][k]["windows"] and w in R["base"][k]["windows"] for c in (0, 1)]
            out.append(f"{w}: {np.sum(np.array(d) > 0.5)}/{len(d)} worse>0.5 (max {max(d):+.2f})" if d else "")
        db = [R[n][k]["bands"]["20-100"][c] - R["base"][k]["bands"]["20-100"][c] for k in keys if sel(k) and "20-100" in R[n][k].get("bands", {}) and "20-100" in R["base"][k].get("bands", {}) for c in (0, 1)]
        print(f"  {n}: " + "; ".join(out) + f"; 20-100 band worse>0.5: {np.sum(np.array(db) > 0.5)}/{len(db)}")
table(lambda k: True, "all captures")
table(est, "captures the table was estimated from (in-sample)")
table(lambda k: not est(k), "all other captures (not used in the estimation)")
