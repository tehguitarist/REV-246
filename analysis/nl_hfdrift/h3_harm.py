"""How big is the harm, really? On the 300 stratified segments (h_seg.pkl): per-capture change of the 10-12 kHz null AND of the full-band
(20 Hz - 12.8 kHz, raw) null for: no correction, current best (Gaussian 5 min + hard skip p90), and soft variants (scalar shrunk smoothly by the
neighbour scatter, optional hard skip at a higher percentile), plus the oracle. Distribution of the harm, absolute null after the harm."""
import sys, json, pickle, collections
import numpy as np
from scipy.signal import butter, sosfiltfilt
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import h_est as H
seg = pickle.load(open(sys.argv[1], "rb"))
FS, KN, KI = H.FS, H.KN, H.KI
z0a, ua = H.fit_shape(H.Z, H.nrun)
def scatter(idx, w, z0, u):
    si = np.array([H.row_s(i, z0, u) for i in idx]); m = np.average(si, weights=w); return float(np.sqrt(np.average((si - m) ** 2, weights=w)))
sc_all = []
for i in range(0, len(H.rows), 2):
    s_n, inf_n, idx, w = H.s_nb(H.KEYS[i], H.T[i], [], z0a, ua); sc_all.append(scatter(idx, w, z0a, ua))
P50, P75, P90, P95 = np.percentile(sc_all, [50, 75, 90, 95])
SOS = butter(4, [10000, 12000], btype="band", fs=FS, output="sos")
def own_runs(tk): return [j for j, (a, b) in enumerate(H.span) if a <= tk <= b]
prep = []
for sg in seg:
    y = sg["y"].astype(float); m = sg["m"].astype(float); L = len(m)
    prep.append(dict(sg=sg, y=y, m=m, L=L, Mf=np.fft.rfft(m, axis=0), f=np.fft.rfftfreq(L) * FS, yb=sosfiltfilt(SOS, y, axis=0)))
def corrected(P, z):
    if z is None: return P["m"]
    g = np.ones(len(KN), complex); g[KI] = np.exp(z); g[10] = 1.0
    E = np.interp(P["f"], np.r_[9500.0, KN], np.r_[1.0, g.real]) + 1j * np.interp(P["f"], np.r_[9500.0, KN], np.r_[0.0, g.imag])
    return np.fft.irfft(P["Mf"] * E[:, None], P["L"], axis=0)
def nulls(P, z):
    mm = corrected(P, z)
    hb = 10 * np.log10(np.sum((P["yb"] - sosfiltfilt(SOS, mm, axis=0)) ** 2, 0) / np.sum(P["yb"] ** 2, 0))
    fb = 10 * np.log10(np.sum((P["y"] - mm) ** 2, 0) / np.sum(P["y"] ** 2, 0))
    return hb, fb
def est(P, mode, c=None, hard=None):
    key, tk = P["sg"]["key"], P["sg"]["t"]; own = own_runs(tk); z0, u = H.shape_for(own)
    s_n, inf, idx, w = H.s_nb(key, tk, own, z0, u, k=15); sc = scatter(idx, w, z0, u)
    if hard is not None and sc > hard: return None
    z, idx5 = H.est("gauss", dict(tau=5), tk, key, own, z0, u)
    if mode == "hard": return z
    s = H.proj(np.log(z), z0, u) if False else None
    # soft: scale the deviation of the correction from unity by 1/(1+(sc/c)^2)
    wgt = 1.0 / (1.0 + (sc / c) ** 2)
    return wgt * z
cands = collections.OrderedDict()
cands["none"] = lambda P: None
cands["current: gauss5, hard skip p90"] = lambda P: est(P, "hard", hard=P90)
cands["gauss5, no skip"] = lambda P: est(P, "hard")
for c in (0.005, 0.01, 0.02):
    cands[f"soft c={c} (no hard skip)"] = (lambda c: lambda P: est(P, "soft", c=c))(c)
    cands[f"soft c={c} + hard skip p95"] = (lambda c: lambda P: est(P, "soft", c=c, hard=P95))(c)
base = [nulls(P, None) for P in prep]
print("n", len(prep), " scatter p50/p75/p90/p95 %.4f %.4f %.4f %.4f" % (P50, P75, P90, P95))
print("%-34s | 10-12k: med gain, >0.5 / >1 / >2 / >3 dB worse (%%), worst | full band 20-12.8k: med gain, >0.5 / >1 dB worse (%%), worst | harmed: median abs null after (10-12k)" % "")
res = {}
for name, fn in cands.items():
    D = []; F = []; A = []; Bn = []
    for P, b in zip(prep, base):
        z = fn(P); hb, fb = nulls(P, z); D.append(hb - b[0]); F.append(fb - b[1]); A.append(hb); Bn.append(b[0])
    D, F, A, Bn = map(np.array, (D, F, A, Bn))
    h = D > 0.5
    print("%-34s | %+.1f, %4.1f / %4.1f / %4.1f / %4.1f, %+.1f | %+.2f, %4.1f / %4.1f, %+.2f | %.1f dB (was %.1f)" % (
        name, -np.median(D), 100 * (D > 0.5).mean(), 100 * (D > 1).mean(), 100 * (D > 2).mean(), 100 * (D > 3).mean(), D.max(),
        -np.median(F), 100 * (F > 0.5).mean(), 100 * (F > 1).mean(), F.max(), np.median(A[h]) if h.any() else np.nan, np.median(Bn[h]) if h.any() else np.nan), flush=True)
    res[name] = dict(D=D.tolist(), F=F.tolist())
json.dump(res, open("/Users/leighpierce/Git/246/analysis/nl_hfdrift/h3_harm.json", "w"))
