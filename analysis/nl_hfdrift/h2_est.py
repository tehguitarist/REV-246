"""Fork H: front-tap-only estimators of the per-capture HF scalar, scored on the capture's reverb window (never used to choose s).
usage: h_est.py SEG.pkl [OUT.json]"""
import os, sys, json, pickle, collections
import numpy as np
from scipy.signal import butter, sosfiltfilt
ROOT = "/Users/leighpierce/Git/246"; OUT = ROOT + "/analysis/nl_hfdrift"
FS = 2400000 / 91
KI = [1, 2, 3, 4, 5, 6, 7, 8, 9]
KN = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
KH = KN[KI]
rows = [r for r in pickle.load(open(OUT + "/d_front_all.pkl", "rb")) if r["t"] is not None]
rows.sort(key=lambda r: r["t"]); T = np.array([r["t"] for r in rows])
runs = [np.array(ix) for ix in np.split(np.arange(len(rows)), np.where(np.diff(T) > 20 * 60)[0] + 1)]
runs = [ix for ix in runs if len(ix) >= 8]
span = [(T[ix].min() - 600, T[ix].max() + 600) for ix in runs]
NUM = {ch: np.array([r[ch]["num"][KI] for r in rows]) for ch in (0, 1)}
DEN = {ch: np.array([r[ch]["den"][KI] for r in rows]) for ch in (0, 1)}
KEYS = [r["key"] for r in rows]; KEYIDX = {k: i for i, k in enumerate(KEYS)}
SPANOF = np.array([r["key"][0] for r in rows])
def split(z): return np.r_[z.real, z.imag]
def zsum(idx, w, ch):
    return np.log((w[:, None] * NUM[ch][idx]).sum(0) / (w[:, None] * DEN[ch][idx]).sum(0))
Z = np.array([zsum(ix, np.ones(len(ix)), 1) for ix in runs]); nrun = np.array([len(ix) for ix in runs], float)
def fit_shape(Zs, w):
    V = np.array([split(z) for z in Zs]); z0 = (w[:, None] * V).sum(0) / w.sum()
    U, S, Vt = np.linalg.svd((V - z0) * np.sqrt(w)[:, None], full_matrices=False)
    u = Vt[0]
    if u[4] > 0: u = -u
    return z0, u
_shape = {}
def shape_for(own):
    k = tuple(own)
    if k not in _shape:
        others = [j for j in range(len(runs)) if j not in own]
        _shape[k] = fit_shape(Z[others], nrun[others])
    return _shape[k]
def proj(z, z0, u): return float(np.dot(split(z) - z0, u) / np.dot(u, u))
def zfrom(z0, u, s): v = z0 + s * u; return v[:9] + 1j * v[9:]

# ---- per-row single-capture estimates (shape-dependent) and information
def row_s(i, z0, u, ch=1):
    return proj(np.log(NUM[ch][i] / DEN[ch][i]), z0, u)
def info(idx_den, u):
    return float(np.sum(idx_den * (u[:9] ** 2 + u[9:] ** 2)))

# ---- estimators: return (z or None) ----
def est(kind, p, tk, key, own, z0, u):
    ok = np.array([KEYS[i] != key for i in range(len(rows))])
    dt = np.abs(T - tk)
    cand = np.where(ok)[0]
    if kind in ("knn", "knn_span", "knn_run"):
        pool = cand
        if kind == "knn_span": pool = cand[SPANOF[cand] == key[0]]
        if kind == "knn_run" and own: pool = np.array([i for i in cand if any(i in runs[j] for j in own)])
        if len(pool) < 3: pool = cand
        idx = pool[np.argsort(dt[pool])[:p["k"]]]
        return zfrom(z0, u, proj(zsum(idx, np.ones(len(idx)), 1), z0, u)), idx
    if kind == "gauss":
        w = np.exp(-0.5 * (dt[cand] / (p["tau"] * 60)) ** 2)
        sel = w > 1e-3
        if sel.sum() < 3:
            idx = cand[np.argsort(dt[cand])[:5]]; return zfrom(z0, u, proj(zsum(idx, np.ones(len(idx)), 1), z0, u)), idx
        idx = cand[sel]; return zfrom(z0, u, proj(zsum(idx, w[sel], 1), z0, u)), idx
    raise ValueError(kind)

def s_nb(key, tk, own, z0, u, k=15, ch=1, tau=None):
    ok = np.array([KEYS[i] != key for i in range(len(rows))]); dt = np.abs(T - tk); cand = np.where(ok)[0]
    if tau is None:
        idx = cand[np.argsort(dt[cand])[:k]]; w = np.ones(len(idx))
    else:
        w = np.exp(-0.5 * (dt[cand] / (tau * 60)) ** 2); sel = w > 1e-3; idx, w = cand[sel], w[sel]
        if len(idx) < 3: idx = cand[np.argsort(dt[cand])[:5]]; w = np.ones(len(idx))
    z = zsum(idx, w, ch)
    inf = info((w[:, None] * DEN[ch][idx]).sum(0), u)
    return proj(z, z0, u), inf, idx, w

if __name__ == "__main__":
    seg = pickle.load(open(sys.argv[1], "rb")); outj = sys.argv[2] if len(sys.argv) > 2 else None
    # ---------------- noise model from front-tap data only: LOO residual of per-row s vs 15-nearest neighbour s, by information
    z0a, ua = fit_shape(Z, nrun)
    res = []
    for i in range(len(rows)):
        s_i = row_s(i, z0a, ua); s_n, inf_n, idx, w = s_nb(KEYS[i], T[i], [], z0a, ua)
        inf_i = info(DEN[1][i], ua); res.append((s_i - s_n, inf_i, inf_n, s_i, s_n))
    res = np.array(res)
    # var(resid) ~ a + b (1/inf_i + 1/inf_n):  least squares on r^2
    X = np.c_[np.ones(len(res)), 1 / res[:, 1] + 1 / res[:, 2]]
    a_, b_ = np.linalg.lstsq(X, res[:, 0] ** 2, rcond=None)[0]; a_ = max(a_, 0); b_ = max(b_, 1e-12)
    print("front-tap noise model: var(s_i - s_nb) = %.4g + %.4g (1/inf_i + 1/inf_nb); median inf_i %.3g, resid sd %.3f" % (a_, b_, np.median(res[:, 1]), res[:, 0].std()))
    var_s = max(np.var(res[:, 3]) - b_ * np.mean(1 / res[:, 1]), 1e-4)       # between-capture variance of the true s
    print("var of s over front-tap captures %.4g (noise-corrected true %.4g)" % (np.var(res[:, 3]), var_s))
    sc_all = []
    # scatter of the neighbours' per-row s around their weighted mean (front-tap data only), for the skip rule
    def scatter(idx, w, z0, u):
        si = np.array([row_s(i, z0, u) for i in idx]); m = np.average(si, weights=w); return float(np.sqrt(np.average((si - m) ** 2, weights=w)))
    for i in range(0, len(rows), 2):
        s_n, inf_n, idx, w = s_nb(KEYS[i], T[i], [], z0a, ua); sc_all.append(scatter(idx, w, z0a, ua))
    sc_all = np.array(sc_all); print("neighbour scatter percentiles 50/75/90:", np.round(np.percentile(sc_all, [50, 75, 90]), 3))
    X75, X90 = np.percentile(sc_all, [75, 90])

    # ---------------- scoring machinery
    SOS = {nm: butter(4, [lo, hi], btype="band", fs=FS, output="sos") for nm, (lo, hi) in dict(C=(10000, 12000), S1=(11000, 11600), S2=(11600, 12300), S3=(12300, 12800)).items()}
    BANDS = list(SOS)
    prep = []
    for sg in seg:
        y = sg["y"].astype(float); m = sg["m"].astype(float); L = len(m)
        yb = {nm: sosfiltfilt(SOS[nm], y, axis=0) for nm in BANDS}
        prep.append(dict(sg=sg, y=y, m=m, L=L, Mf=np.fft.rfft(m, axis=0), f=np.fft.rfftfreq(L) * FS, yb=yb, yn={nm: np.sum(yb[nm] ** 2, 0) for nm in BANDS}))
    def nulls(P, z):
        if z is None: mm = P["m"]
        else:
            g = np.ones(len(KN), complex); g[KI] = np.exp(z); g[10] = 1.0
            kx = np.r_[9500.0, KN]; gg = np.r_[1.0, g]
            E = np.interp(P["f"], kx, gg.real) + 1j * np.interp(P["f"], kx, gg.imag)
            mm = np.fft.irfft(P["Mf"] * E[:, None], P["L"], axis=0)
        out = {}
        for nm in BANDS:
            mb = sosfiltfilt(SOS[nm], mm, axis=0); out[nm] = 10 * np.log10(np.sum((P["yb"][nm] - mb) ** 2, 0) / P["yn"][nm])
        return out
    base = [nulls(P, None) for P in prep]
    def own_runs(tk): return [j for j, (a, b) in enumerate(span) if a <= tk <= b]
    # ---------------- candidate estimators
    def run_est(name, fn):
        zs = []
        for P in prep:
            key = P["sg"]["key"]; tk = P["sg"]["t"]; own = own_runs(tk); z0, u = shape_for(own)
            zs.append(fn(key, tk, own, z0, u))
        return zs
    def f_knn(k, kind="knn"):
        return lambda key, tk, own, z0, u: est(kind, dict(k=k), tk, key, own, z0, u)[0]
    def f_gauss(tau): return lambda key, tk, own, z0, u: est("gauss", dict(tau=tau), tk, key, own, z0, u)[0]
    def f_own(mode, k=15, tau=None):
        def fn(key, tk, own, z0, u):
            s_n, inf_n, idx, w = s_nb(key, tk, own, z0, u, k=k, tau=tau)
            if key not in KEYIDX:
                return zfrom(z0, u, s_n)                                     # no front taps: neighbours
            i = KEYIDX[key]; s_i = row_s(i, z0, u); inf_i = info(DEN[1][i], u)
            if mode == "own": return zfrom(z0, u, s_i)
            # precision weighting: own variance b/inf_i, neighbour variance a + b/inf_n
            v_i = b_ / inf_i; v_n = a_ + b_ / inf_n
            wl = v_n / (v_i + v_n); return zfrom(z0, u, wl * s_i + (1 - wl) * s_n)
        return fn
    def f_lr(k=15):
        def fn(key, tk, own, z0, u):
            ok = np.array([KEYS[i] != key for i in range(len(rows))]); dt = np.abs(T - tk); cand = np.where(ok)[0]
            idx = cand[np.argsort(dt[cand])[:k]]
            zR = zsum(idx, np.ones(len(idx)), 1); infR = info((DEN[1][idx]).sum(0), u)
            sR = proj(zR, z0, u)
            keep = list(range(8)); kk = keep + [9 + j for j in keep]                 # L usable to 12.25 kHz = knots 0..7
            zL = zsum(idx, np.ones(len(idx)), 0); uL = u[kk]; z0L = z0[kk]
            sL = float(np.dot(split(zL)[kk] - z0L, uL) / np.dot(uL, uL))
            dL = DEN[0][idx].sum(0); infL = float(np.sum(dL[keep] * (u[keep] ** 2 + u[9 + np.array(keep)] ** 2)))
            s = (infR * sR + infL * sL) / (infR + infL); return zfrom(z0, u, s)
        return fn
    def f_shrink(target, k=15, tau=None, strength=1.0):
        def fn(key, tk, own, z0, u):
            s_n, inf_n, idx, w = s_nb(key, tk, own, z0, u, k=k, tau=tau)
            v = a_ + b_ / inf_n
            B = var_s / (var_s + strength * v)
            if target == "zero": s = B * s_n
            else:
                runm = 0.0
                if own:
                    ri = [i for j in own for i in runs[j] if KEYS[i] != key]
                    runm = proj(zsum(np.array(ri), np.ones(len(ri)), 1), z0, u)
                s = runm + B * (s_n - runm)
            return zfrom(z0, u, s)
        return fn
    def f_skip(X, k=15, mode="unity", trimmed_skip=False, nofront_skip=False):
        def fn(key, tk, own, z0, u):
            s_n, inf_n, idx, w = s_nb(key, tk, own, z0, u, k=k)
            sc = scatter(idx, w, z0, u)
            trimmed = P_trim[key]
            if sc > X or (trimmed_skip and trimmed) or (nofront_skip and key not in KEYIDX):
                return None if mode == "unity" else zfrom(z0, u, 0.0)
            return zfrom(z0, u, s_n)
        return fn
    P_trim = {P["sg"]["key"]: (P["sg"]["cluster"] == "trimmed") for P in prep}
    cands = collections.OrderedDict()
    cands["none"] = lambda key, tk, own, z0, u: None
    cands["mean z0 only"] = lambda key, tk, own, z0, u: zfrom(z0, u, 0.0)
    for k in (5, 15, 40): cands[f"knn k={k}"] = f_knn(k)
    cands["knn15 same-SPAN"] = f_knn(15, "knn_span"); cands["knn15 same-run"] = f_knn(15, "knn_run")
    for tau in (5, 20, 60): cands[f"gauss tau={tau}min"] = f_gauss(tau)
    cands["own front taps alone"] = f_own("own")
    cands["own+knn15 precision-weighted"] = f_own("mix")
    cands["own+gauss20 precision-weighted"] = f_own("mix", tau=20)
    cands["knn15 + L front taps"] = f_lr()
    cands["shrink knn15 -> 0"] = f_shrink("zero")
    cands["shrink knn15 -> run mean"] = f_shrink("run")
    cands["shrink x2 -> 0"] = f_shrink("zero", strength=2.0)
    cands["shrink x2 -> run mean"] = f_shrink("run", strength=2.0)
    cands["skip if scatter>p75 (unity)"] = f_skip(X75); cands["skip if scatter>p90 (unity)"] = f_skip(X90)
    cands["skip p75 + trimmed + no-front"] = f_skip(X75, trimmed_skip=True, nofront_skip=True)
    cands["skip p90 + trimmed"] = f_skip(X90, trimmed_skip=True)
    cands["skip p75 -> mean z0"] = f_skip(X75, mode="mean")
    def f_skipx(X, base_fn, mode="unity", k=15):
        def fn(key, tk, own, z0, u):
            s_n, inf_n, idx, w = s_nb(key, tk, own, z0, u, k=k); sc = scatter(idx, w, z0, u)
            if sc > X: return None if mode == "unity" else zfrom(z0, u, 0.0)
            return base_fn(key, tk, own, z0, u)
        return fn
    def f_half(base_fn, a):
        def fn(key, tk, own, z0, u):
            z = base_fn(key, tk, own, z0, u)
            return None if z is None else a * z
        return fn
    g5 = f_gauss(5); k5 = f_knn(5)
    cands["gauss5 + skip p90"] = f_skipx(X90, g5); cands["gauss5 + skip p75"] = f_skipx(X75, g5)
    cands["knn5 + skip p90"] = f_skipx(X90, k5); cands["knn5 + skip p75"] = f_skipx(X75, k5)
    cands["gauss5 + skip p90 -> mean z0"] = f_skipx(X90, g5, mode="mean")
    cands["gauss5 + skip p90, x0.8"] = f_half(f_skipx(X90, g5), 0.8)
    cands["gauss5 + skip p90, x0.6"] = f_half(f_skipx(X90, g5), 0.6)
    cands["gauss5 + skip p95"] = f_skipx(float(np.percentile(sc_all, 95)), g5)
    cands["gauss5 + skip p85"] = f_skipx(float(np.percentile(sc_all, 85)), g5)
    # ---------------- oracle (diagnostic, uses the reverb): best single s on a grid, 'none' included
    results = {}
    def score(name, zs):
        D = collections.defaultdict(list); N = collections.defaultdict(list)
        for P, z, b0 in zip(prep, zs, base):
            n1 = nulls(P, z)
            for nm in BANDS:
                N[nm].append(n1[nm]); D[nm].append(n1[nm] - b0[nm])
        row = {}
        for nm in BANDS:
            A = np.array(N[nm]); Dd = np.array(D[nm])
            row[nm] = dict(med=float(np.median(A)), medL=float(np.median(A[:, 0])), medR=float(np.median(A[:, 1])), harm=float((Dd > 0.5).mean()), worst=float(Dd.max()),
                           gain=float(-np.median(Dd)))
        row["cluster"] = {}
        for c in ("rest", "S0.1", "S0.3", "trimmed"):
            ii = [j for j, P in enumerate(prep) if P["sg"]["cluster"] == c]
            Dd = np.array([D["C"][j] for j in ii]); A = np.array([N["C"][j] for j in ii])
            row["cluster"][c] = dict(n=len(ii), med=float(np.median(A)), harm=float((Dd > 0.5).mean()), worst=float(Dd.max()))
        results[name] = row
        print("%-34s C med %.1f (L %.1f / R %.1f) gain %+.1f harm %4.1f%% worst %+.1f | S1 %.1f S2 %.1f S3 %.1f | " % (
            name, row["C"]["med"], row["C"]["medL"], row["C"]["medR"], row["C"]["gain"], 100 * row["C"]["harm"], row["C"]["worst"], row["S1"]["med"], row["S2"]["med"], row["S3"]["med"])
            + " ".join("%s %.1f/%.0f%%" % (c, v["med"], 100 * v["harm"]) for c, v in row["cluster"].items()), flush=True)
    for name, fn in cands.items():
        score(name, run_est(name, fn))
    # oracle
    grid = np.arange(-0.25, 0.51, 0.05)
    zs = []
    for P in prep:
        key = P["sg"]["key"]; tk = P["sg"]["t"]; own = own_runs(tk); z0, u = shape_for(own)
        best = (np.mean(base[prep.index(P)]["C"]), None)
        for s in grid:
            v = np.mean(nulls(P, zfrom(z0, u, s))["C"])
            if v < best[0]: best = (v, s)
        zs.append(None if best[1] is None else zfrom(z0, u, best[1]))
    score("ORACLE best s (uses reverb)", zs)
    # the chosen s of the current estimator vs oracle (diagnostic)
    if outj: json.dump(results, open(outj, "w"), indent=1)
