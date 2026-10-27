import sys, pickle, numpy as np, datetime
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
FS = 26373.626373626373
KN = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
rows = pickle.load(open("/Users/leighpierce/Git/246/analysis/nl_hfdrift/d_front_all.pkl", "rb"))
rows = [r for r in rows if r["t"] is not None]
print(len(rows), "usable captures with timestamps")
sel = [4, 6, 7, 8, 9]       # 11.5, 12, 12.25, 12.5, 12.8 kHz
def prof(rs, ch):
    num = sum(r[ch]["num"] for r in rs); den = sum(r[ch]["den"] for r in rs); c = num / den
    return np.angle(c), 20 * np.log10(np.abs(c))
print("knots kHz", KN[sel] / 1000)
spans = sorted(set(r["key"][0] for r in rows))
for ch in (0, 1):
    print("== channel", ("L", "R")[ch])
    P, G = prof(rows, ch); print("  pooled  phase mrad", np.round(P[sel] * 1e3, 1), " level dB", np.round(G[sel], 3), " delay-equiv samples", np.round(-P[sel] / (2 * np.pi * KN[sel] / FS), 4))
    for s in spans:
        rs = [r for r in rows if r["key"][0] == s]
        P, G = prof(rs, ch)
        print("  SPAN %.1f n=%3d phase mrad %s level dB %s  delay-eq %s" % (s, len(rs), np.round(P[sel] * 1e3, 0), np.round(G[sel], 3), np.round(-P[sel] / (2 * np.pi * KN[sel] / FS), 4)))
# per-capture scalar: phase at 12.5 kHz (mrad) and 12.0
def pc(r, ch, k):
    c = r[ch]["num"][k] / r[ch]["den"][k]; return np.angle(c)
ph125 = np.array([[pc(r, 0, 8), pc(r, 1, 8)] for r in rows]) * 1e3
ph120 = np.array([[pc(r, 0, 6), pc(r, 1, 6)] for r in rows]) * 1e3
ppm = np.array([r["ppm"] for r in rows]); o = np.array([r["o"] for r in rows]); fr = np.array([r["fr"] for r in rows]); t = np.array([r["t"] for r in rows]) / 3600
span = np.array([r["key"][0] for r in rows]); early = np.array([r["key"][1] for r in rows]); late = np.array([r["key"][2] for r in rows])
d0 = np.array([np.median(r[0]["d"]) for r in rows]); d1 = np.array([np.median(r[1]["d"]) for r in rows])
print("per-capture phase@12.5 kHz mrad: median L %.1f R %.1f, sd (robust) L %.1f R %.1f" % (np.median(ph125[:, 0]), np.median(ph125[:, 1]), 1.48 * np.median(abs(ph125[:, 0] - np.median(ph125[:, 0]))), 1.48 * np.median(abs(ph125[:, 1] - np.median(ph125[:, 1])))))
def rho(a, b, m=None):
    m = np.ones(len(a), bool) if m is None else m
    return np.corrcoef(a[m], b[m])[0, 1]
x = ph125[:, 1]
print("corr of R phase@12.5 with: ppm %.2f  o %.2f  front_shift %.2f  time %.2f  early %.2f  late %.2f  tapoffset(L) %.2f" % (rho(x, ppm), rho(x, np.clip(o, -2, 2)), rho(x, fr), rho(x, t), rho(x, early), rho(x, late), rho(x, d0)))
x = ph125[:, 0]
print("corr of L phase@12.5 with: ppm %.2f  o %.2f  front_shift %.2f  time %.2f" % (rho(x, ppm), rho(x, np.clip(o, -2, 2)), rho(x, fr), rho(x, t)))
print("within-session corr (R phase@12.5, ppm | time):", {("%.1f" % s): (round(rho(ph125[:, 1], ppm, span == s), 2), round(rho(ph125[:, 1], t, span == s), 2)) for s in spans})
print("session mean ppm:", {("%.1f" % s): round(ppm[span == s].mean(), 2) for s in spans}, " sd within:", {("%.1f" % s): round(ppm[span == s].std(), 2) for s in spans})
# time course: bins of 2 h across all sessions
order = np.argsort(t); tb = t[order]
nb = 12; edges = np.quantile(tb, np.linspace(0, 1, nb + 1))
print("time bins (hours from first): R phase@12.5 median mrad | ppm mean | n | sessions")
for i in range(nb):
    m = (t >= edges[i]) & (t <= edges[i + 1])
    print("  %6.1f-%6.1f h  %6.1f | %5.2f | %3d | %s" % (edges[i] - t.min(), edges[i + 1] - t.min(), np.median(ph125[m, 1]), ppm[m].mean(), m.sum(), sorted(set(span[m]))))
