"""Fork D: per-run (capture-run) HF correction. Pooled shape + ONE scalar per run estimated on the run's FRONT TAPS (R channel), applied to
both channels, held out on the run's reverb windows (leave-one-run-out for the shape).
usage: d_perrun.py CAPS_PER_RUN"""
import os, sys, json, pickle, random, datetime
os.environ["NL_WARP2"] = "1"; os.environ.pop("NL_BK", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter, sosfiltfilt
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
import defects
KN = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
KI = [1, 2, 3, 4, 5, 6, 7, 8, 9]                            # 10.5 ... 12.8 kHz (the small level step starts below 11.5 kHz)
rows = [r for r in pickle.load(open(OUT + "/d_front_all.pkl", "rb")) if r["t"] is not None]
rows.sort(key=lambda r: r["t"]); t = np.array([r["t"] for r in rows])
runs = [np.array(ix) for ix in np.split(np.arange(len(rows)), np.where(np.diff(t) > 20 * 60)[0] + 1)]
runs = [ix for ix in runs if len(ix) >= 8]
times = json.load(open(OUT + "/d_times.json"))
def zvec(idx, ch=1):
    c = sum(rows[i][ch]["num"] for i in idx) / sum(rows[i][ch]["den"] for i in idx)
    return np.log(c[KI])                                      # complex log: ln|c| + i phase
Z = np.array([zvec(ix) for ix in runs]); nrun = np.array([len(ix) for ix in runs])
def split(z): return np.r_[z.real, z.imag]
def unsplit(v): return v[:len(KI)] + 1j * v[len(KI):]
def fit_shape(Zs, w):
    V = np.array([split(z) for z in Zs]); z0 = (w[:, None] * V).sum(0) / w.sum()
    U, S, Vt = np.linalg.svd((V - z0) * np.sqrt(w)[:, None], full_matrices=False)
    return z0, Vt[0], S
z0_all, u_all, S = fit_shape(Z, nrun.astype(float))
print("runs", len(runs), "n front", nrun.tolist(), " singular values", np.round(S, 3), " rank-1 share %.3f" % (S[0] ** 2 / (S ** 2).sum()))
print("shape u (gain part, phase part mrad-scaled):", np.round(u_all[:len(KI)], 3), np.round(u_all[len(KI):] * 1e3, 1))
def E_from(z):                                               # knots complex gain, unity below 9.5 kHz and above 12.8 kHz
    g = np.ones(len(KN), complex); g[KI] = np.exp(z); g[10] = 1.0
    return g
def apply(m, g):
    L = len(m); f = np.fft.rfftfreq(L) * FS
    kx = np.r_[9500.0, KN]; gg = np.r_[1.0, g]
    E = np.interp(f, kx, gg.real) + 1j * np.interp(f, kx, gg.imag)
    return np.fft.irfft(np.fft.rfft(m, axis=0) * E[:, None], L, axis=0)
SOS = {nm: butter(4, [lo, hi], btype="band", fs=FS, output="sos") for nm, (lo, hi) in dict(B5k=(5000, 10000), HI1=(9000, 11000), S1=(11000, 11600), S2=(11600, 12300), S3=(12300, 12800), C=(10000, 12000)).items()}
def null(y, m, nm):
    yb, mb = sosfiltfilt(SOS[nm], y, axis=0), sosfiltfilt(SOS[nm], m, axis=0)
    return np.array([10 * np.log10(np.sum((yb[:, c] - mb[:, c]) ** 2) / np.sum(yb[:, c] ** 2)) for c in (0, 1)])
ncap = int(sys.argv[1]) if len(sys.argv) > 1 else 5
random.seed(5)
allkeys = {tuple(float(v) for v in k.split("_")): tt for k, tt in times.items()}
res = []
for j, ix in enumerate(runs):
    t0, t1 = t[ix].min() - 600, t[ix].max() + 600
    ks = [k for k, tt in allkeys.items() if t0 <= tt <= t1 and k not in defects.NL_DEFECTS and k in nl.table()
          and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json")]
    ks = random.sample(sorted(ks), min(ncap, len(ks)))
    others = [i for i in range(len(runs)) if i != j]
    z0, u, _ = fit_shape(Z[others], nrun[others].astype(float))
    z0_pool = (nrun[others][:, None] * Z[others]).sum(0) / nrun[others].sum()
    s = float(np.dot(split(Z[j]) - z0, u) / np.dot(u, u))
    z_run = unsplit(z0 + s * u); z_pool = z0_pool
    z_oracle = Z[j]
    for key in ks:
        b = common.bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
        fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
        a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
        if b2 - a < int(0.05 * FS): continue
        n = b2 + 200
        rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
        y = nlt.unwarp(rawp, b["ppm"])[a:b2]
        m = (M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"]))[a:b2]
        zs = zvec([i for i in range(len(rows)) if rows[i]["key"][0] == key[0]])           # SPAN-session pooled (own front taps)
        tk = allkeys[key]; order = sorted((i for i in range(len(rows)) if rows[i]["key"] != key), key=lambda i: abs(rows[i]["t"] - tk))[:15]
        zl = zvec(order); s_loc = float(np.dot(split(zl) - z0, u) / np.dot(u, u)); z_loc = unsplit(z0 + s_loc * u)   # 15 nearest-in-time front-tap captures, shape from the other runs
        r = dict(run=j, key=key, s=s, s_loc=s_loc)
        for nm in SOS:
            r[nm] = {v: null(y, m if g is None else apply(m, E_from(g)), nm) for v, g in (("none", None), ("pooled", z_pool), ("span", zs), ("run", z_run), ("local", z_loc), ("oracle", z_oracle))}
        res.append(r)
    print("run %2d (%s) n_front %3d s=%+.2f  captures scored %d" % (j, datetime.datetime.fromtimestamp(t[ix].min()).strftime("%m-%d %H:%M"), len(ix), s, len([r for r in res if r["run"] == j])), flush=True)
pickle.dump(dict(res=res, runs=[ix.tolist() for ix in runs], Z=Z, z0=z0_all, u=u_all), open(OUT + "/d_perrun_%s.pkl" % ncap, "wb"))
print("\nheld out (shape from the other runs, scalar from the run's own front taps); median over captures, L / R null dB; (worse by >0.5 dB vs none)")
for nm in SOS:
    line = "%-4s" % nm
    for v in ("none", "pooled", "span", "run", "local", "oracle"):
        A = np.array([r[nm][v] for r in res]); B = np.array([r[nm]["none"] for r in res])
        line += " | %s %.1f/%.1f (%d/%d)" % (v, np.median(A[:, 0]), np.median(A[:, 1]), int(((A - B) > 0.5).sum()), A.size)
    print(line)
