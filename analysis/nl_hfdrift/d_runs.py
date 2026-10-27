import sys, pickle, numpy as np, datetime
FS = 26373.626373626373
KN = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
rows = [r for r in pickle.load(open("/Users/leighpierce/Git/246/analysis/nl_hfdrift/d_front_all.pkl", "rb")) if r["t"] is not None]
rows.sort(key=lambda r: r["t"])
t = np.array([r["t"] for r in rows]); gaps = np.where(np.diff(t) > 20 * 60)[0]
runs = np.split(np.arange(len(rows)), gaps + 1)
print(len(rows), "captures,", len(runs), "runs (gap > 20 min); capture interval median %.0f s" % np.median(np.diff(t)))
def pooled(idx, ch, k):
    num = sum(rows[i][ch]["num"][k] for i in idx); den = sum(rows[i][ch]["den"][k] for i in idx); c = num / den
    return np.angle(c) * 1e3, 20 * np.log10(abs(c))
out = []
print(" run  start            n  spans   ppm(med)  R ph@11.5 @12.0 @12.5 @12.8 mrad | R gain@12.0 dB | L ph@12.0")
for j, idx in enumerate(runs):
    if len(idx) < 8: continue
    sp = sorted(set(rows[i]["key"][0] for i in idx)); ppm = np.median([rows[i]["ppm"] for i in idx])
    ph = [pooled(idx, 1, k)[0] for k in (4, 6, 8, 9)]; g = pooled(idx, 1, 6)[1]; lph = pooled(idx, 0, 6)[0]
    out.append((j, t[idx].mean(), ppm, ph, g, len(idx), lph))
    print("%3d  %s %3d  %-12s %5.2f   %6.1f %6.1f %6.1f %6.1f | %6.3f | %6.1f" % (j, datetime.datetime.fromtimestamp(t[idx[0]]).strftime("%m-%d %H:%M"), len(idx), ",".join("%.1f" % s for s in sp), ppm, *ph, g, lph))
P = np.array([[o[2], *o[3], o[4], o[6]] for o in out])
names = ["ppm", "R ph11.5", "R ph12.0", "R ph12.5", "R ph12.8", "R gain12", "L ph12"]
print("across-run correlation matrix (n=%d runs)" % len(out)); C = np.corrcoef(P.T)
print("            " + " ".join("%9s" % n for n in names))
for i, n in enumerate(names): print("%-11s " % n + " ".join("%9.2f" % c for c in C[i]))
