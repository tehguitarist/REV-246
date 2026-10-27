"""The UAD EMT 250 plugin's moving reads, measured from its 960 impulse snapshots at Time 0.4 (prior for Program E).

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/uad_moving.py --stage t1
    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/uad_moving.py --stage t2 [--sub ...]

Data (read-only, analysis/e250/uad/): dense_T0.4_24k.npy (960, 16800, 2) = (L, R) snapshots on the plugin's 24 kHz grid, R's
first tap at index 263; dense_T0.4_t.npy the snapshot times (0.25 s steps). Outputs: analysis/e250/uad/uad_moving_*.json
(plugin-derived: not committed).

Stage t1: L's moving part against R's (residual = snapshot - mean): the plain correlation lag pattern around the static
alignment (285 samples), a least-squares FIR from R's residual to L's (fitted on one half of the snapshots, tested on the
other), its structure, the transfer function (coherence, gain, phase) and the same on the static mean.
Stage t2: the moving reads themselves (strong taps of the base residual: bands, lockstep displacement g(t), loops).
"""
import argparse, json, os
import numpy as np
from scipy.signal import fftconvolve

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
UAD = os.path.join(ROOT, "analysis", "e250", "uad")
R0, LAG, FS = 263, 285, 24000.0
NFFT = 1 << 15


def load():
    X = np.load(os.path.join(UAD, "dense_T0.4_24k.npy"))
    t = np.load(os.path.join(UAD, "dense_T0.4_t.npy"))
    m = X.astype(np.float64).mean(0)
    return X, m, t


def make_base(m):
    """Deconvolver: chain kernel (tap 1's own response, from the mean) and the three all-passes divided out, regularised."""
    w = np.fft.rfftfreq(NFFT)
    z = lambda d: np.exp(-2j * np.pi * w * d)
    ap = np.ones_like(w, dtype=complex)
    for D, g in ((41, .70), (132, .69), (356, .67)): ap *= (-g + z(D)) / (1 - g * z(D))
    Kf = {}
    for ch, off in ((1, 0), (0, LAG)):
        k = m[R0 + off - 8:R0 + off + 370, ch] / m[R0 + off, ch]
        Kf[ch] = np.fft.rfft(k, NFFT) * z(-8)

    def base(x, ch, nout=9000):
        off = 0 if ch == 1 else LAG
        xx = x.astype(np.float64).copy(); xx[:R0 + off + 680] = 0; xx = xx[:NFFT]
        K = Kf[ch]
        y = np.fft.irfft(np.fft.rfft(xx, NFFT) / ap * np.conj(K) / (np.abs(K) ** 2 + 1e-2 * np.max(np.abs(K)) ** 2), NFFT)
        return y[R0 + off:R0 + off + nout]
    return base


def base_residuals(X, m, nout=9000, cache=True):
    """(n, nout) base residuals of R and L: (snapshot - mean) deconvolved; index = samples after that channel's tap 1."""
    pr = os.path.join(UAD, "uad_moving_base_R.npy"); pl = os.path.join(UAD, "uad_moving_base_L.npy")
    if cache and os.path.exists(pr) and os.path.exists(pl):
        return np.load(pr), np.load(pl)
    base = make_base(m)
    Rb = np.array([base(X[i, :, 1] - m[:, 1], 1, nout) for i in range(len(X))]).astype(np.float32)
    Lb = np.array([base(X[i, :, 0] - m[:, 0], 0, nout) for i in range(len(X))]).astype(np.float32)
    if cache: np.save(pr, Rb); np.save(pl, Lb)
    return Rb, Lb


def toeplitz_solve(rxx, rxy, K, ridge=1e-6):
    G = np.array([[rxx[abs(k - l)] for k in range(-K, K + 1)] for l in range(-K, K + 1)])
    G += ridge * np.trace(G) / G.shape[0] * np.eye(G.shape[0])
    return np.linalg.solve(G, rxy)


def fir_fit(x, y, K, sel, a, b):
    """y[n] ~ sum_k h[k] x[n-k], k=-K..K, n in [a, b) (R-relative; x at R0+n, y at R0+LAG+n); sums over snapshots sel."""
    xs = x[sel][:, R0 + a - K:R0 + b + K]
    ys = y[sel][:, R0 + LAG + a:R0 + LAG + b]
    # zero-pad y to x's span so the correlation lags are consistent
    Y = np.zeros_like(xs); Y[:, K:K + b - a] = ys
    Fx = np.fft.rfft(xs, NFFT, axis=1); Fy = np.fft.rfft(Y, NFFT, axis=1)
    rxx = np.fft.irfft((np.abs(Fx) ** 2).sum(0), NFFT)
    rxy = np.fft.irfft((np.conj(Fx) * Fy).sum(0), NFFT)   # sum_n x[n] y[n+tau] (tau relative to the padded frame)
    lags = np.arange(-K, K + 1)
    rhs = np.array([rxy[l] if l >= 0 else rxy[NFFT + l] for l in lags])
    h = toeplitz_solve(rxx, rhs, K)
    return h


def fir_apply(x, h, K, sel, a, b):
    xs = x[sel][:, R0 + a - K:R0 + b + K]
    yh = fftconvolve(xs, h[None, :], axes=1)        # full conv; y_hat[n] = sum_k h[k] x[n-k]
    # index: xs[j] <-> n = a-K+j ; h index i <-> k = i-K ; conv[j+i] <-> n = a-K+j+i-K ... align by position
    off = 2 * K                                      # conv index of n=a  (j=K, i=K)
    return yh[:, off:off + (b - a)]


def score(x, y, h, K, sel, a, b):
    yt = y[sel][:, R0 + LAG + a:R0 + LAG + b]
    yh = fir_apply(x, h, K, sel, a, b)
    num = ((yt - yh) ** 2).sum(); den = (yt ** 2).sum()
    return float(1 - num / den), float(den / sel.sum())


def corr_pattern(x, y, a, b, lags):
    xs = x[:, R0 + a:R0 + b]
    out = []
    for d in lags:
        ys = y[:, R0 + LAG + a + d:R0 + LAG + b + d]
        out.append(float((xs * ys).sum() / np.sqrt((xs ** 2).sum() * (ys ** 2).sum())))
    return out


def stage_t1(args):
    X, m, t = load()
    Rr = (X[:, :, 1] - m[:, 1]).astype(np.float64); Lr = (X[:, :, 0] - m[:, 0]).astype(np.float64)
    del X
    n = len(Rr); half = n // 2
    A = np.arange(n) < half; B = ~A
    out = {"snapshots": n}
    wins = {"405-805": (405, 805), "805-1600": (805, 1600), "1600-2400": (1600, 2400), "2400-4000": (2400, 4000),
            "4000-8000": (4000, 8000), "805-8000": (805, 8000), "405-8000": (405, 8000)}
    lags = list(range(-3, 4))
    out["corr_lags"] = lags
    out["corr"] = {k: [round(v, 3) for v in corr_pattern(Rr, Lr, a, b, lags)] for k, (a, b) in wins.items()}
    # energy of L's residual vs R's
    out["energy_L_over_R"] = {k: round(float((Lr[:, R0 + LAG + a:R0 + LAG + b] ** 2).sum() / (Rr[:, R0 + a:R0 + b] ** 2).sum()), 3)
                              for k, (a, b) in wins.items()}
    K = args.K
    res = {}
    for k, (a, b) in wins.items():
        r = {}
        for nm, tr, te in (("A->B", A, B), ("B->A", B, A)):
            h = fir_fit(Rr, Lr, K, tr, a, b)
            ex, _ = score(Rr, Lr, h, K, te, a, b)
            exin, _ = score(Rr, Lr, h, K, tr, a, b)
            # plain opposite-sign copy at the best integer delay (trained on tr, tested on te), and a 3-tap fit around it
            best = None
            for d in range(-8, 9):
                hh = np.zeros(2 * K + 1); hh[K + d] = 1.0
                # LS gain on train
                yt = Lr[tr][:, R0 + LAG + a:R0 + LAG + b]; yh = fir_apply(Rr, hh, K, tr, a, b)
                gain = float((yt * yh).sum() / (yh ** 2).sum())
                e, _ = score(Rr, Lr, gain * hh, K, te, a, b)
                if best is None or e > best[0]: best = (e, d, gain)
            top = np.argsort(-np.abs(h))[:8]
            r[nm] = {"fir_heldout_explained": round(ex, 4), "fir_insample_explained": round(exin, 4),
                     "single_tap_heldout": [round(best[0], 4), best[1], round(best[2], 3)],
                     "fir_top_taps": [[int(i - K), round(float(h[i]), 3)] for i in top]}
        res[k] = r
    out["fir"] = res
    # FIR on the 805-8000 window, trained on all, structure + transfer function
    a, b = 805, 8000
    hall = fir_fit(Rr, Lr, K, np.ones(n, bool), a, b)
    out["fir_all_805_8000"] = [round(float(v), 4) for v in hall]
    # transfer function: cross-spectra over snapshots, window a..b, Hann-free (rectangular window; many snapshots)
    xs = Rr[:, R0 + a:R0 + b]; ys = Lr[:, R0 + LAG + a:R0 + LAG + b]
    Fx = np.fft.rfft(xs, 16384, axis=1); Fy = np.fft.rfft(ys, 16384, axis=1)
    Sxx = (np.abs(Fx) ** 2).sum(0); Syy = (np.abs(Fy) ** 2).sum(0); Sxy = (np.conj(Fx) * Fy).sum(0)
    fr = np.fft.rfftfreq(16384, 1 / FS)
    bands = [(100, 250), (250, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 8000), (8000, 11000)]
    tf = {}
    for lo, hi in bands:
        s = (fr >= lo) & (fr < hi)
        sxy = Sxy[s].sum(); coh = abs(sxy) ** 2 / (Sxx[s].sum() * Syy[s].sum())
        H = Sxy[s] / Sxx[s]
        # gain and phase of -H averaged (weighted by Sxx)
        hw = (H * Sxx[s]).sum() / Sxx[s].sum()
        # group delay (samples) from the unwrapped phase of -H over the band, slope
        ph = np.unwrap(np.angle(-H)); w = Sxx[s]
        sl = np.polyfit(2 * np.pi * fr[s] / FS, ph, 1, w=np.sqrt(w))[0]
        tf[f"{lo}-{hi}"] = {"coherence": round(float(coh), 3), "gain_minusH": round(float(abs(hw)), 3),
                             "phase_minusH_deg": round(float(np.degrees(np.angle(-hw))), 1), "group_delay_samples": round(float(-sl), 2)}
    out["transfer_805_8000"] = tf
    # static part: FIR from R's mean to L's mean, direct region and reverb region (no splits; one realisation)
    Xs, ms, _ = load(); del Xs
    xm = ms[None, :, 1]; ym = ms[None, :, 0]
    st = {}
    for nm, (a2, b2) in {"direct_-5-700": (-5, 700), "reverb_805-4000": (805, 4000)}.items():
        h = fir_fit(xm, ym, K, np.ones(1, bool), a2, b2)
        ex, _ = score(xm, ym, h, K, np.ones(1, bool), a2, b2)
        top = np.argsort(-np.abs(h))[:6]
        st[nm] = {"explained": round(ex, 5), "top_taps": [[int(i - K), round(float(h[i]), 3)] for i in top]}
        # the chain kernel is long (410 samples): a plain +-delta copy for comparison
        for sg in (+1, -1):
            hh = np.zeros(2 * K + 1); hh[K] = sg
            e, _ = score(xm, ym, hh, K, np.ones(1, bool), a2, b2)
            st[nm][f"plain_{'+' if sg > 0 else '-'}delta_at_285_explained"] = round(e, 5)
    out["static"] = st
    p = os.path.join(UAD, "uad_moving_t1.json")
    json.dump(out, open(p, "w"), indent=1)
    print(json.dumps(out, indent=1)[:6000])


PL = (2009, 2280, 2523)           # loop lengths A / B / C (plugin, 24 kHz)


def strong_peaks(B, thr, lo=300):
    a = np.abs(B)
    out = []
    for i in range(len(B)):
        x = a[i]
        p = np.where((x[1:-1] > thr) & (x[1:-1] >= x[:-2]) & (x[1:-1] >= x[2:]))[0] + 1
        out.append(p[p >= lo])
    return out


def find_bands(Rb, thr=0.0075):
    """Bands (centre, sign) of strong taps in the first pass: runs of the smoothed position histogram; a run about twice
    the 66-sample travel is two reads."""
    cent = []
    for sg in (1, -1):
        h = np.zeros(Rb.shape[1])
        for i in range(len(Rb)):
            x = sg * Rb[i]
            p = np.where((x[1:-1] > thr) & (x[1:-1] >= x[:-2]) & (x[1:-1] >= x[2:]))[0] + 1
            h[p] += 1
        hs = np.convolve(h, np.ones(5), "same"); on = hs > 2; i = 400
        while i < 2600:
            if on[i]:
                j = i
                while j < 2600 and on[j]: j += 1
                if h[i:j].sum() >= 100:
                    w = j - i
                    if w > 90: cent += [((i + 33 + 0.5), sg), ((j - 33 - 0.5), sg)]
                    else: cent.append(((i + j) / 2, sg))
                i = j
            else: i += 1
    return sorted(cent)


def track_g(Rb, t, C, SG, iters=4, thr=0.0075, v0=70.0):
    """Common displacement g_i (and its sign-velocity s_i) per snapshot, read centres C refined: p = C_j + g(t_i + p/fs)."""
    n = len(Rb); g = np.zeros(n); s = np.zeros(n); C = np.array(C, float)
    gg = np.arange(-40, 40.01, 0.5)
    def taps(i, gi, si, tol):
        out = []
        for j in range(len(C)):
            pr = C[j] + gi + si * (C[j] / FS - 0.058); lo = int(round(pr - tol)); hi = int(round(pr + tol)) + 1
            seg = SG[j] * Rb[i, lo:hi]; k = int(np.argmax(seg))
            if seg[k] > thr: out.append((j, lo + k, seg[k]))
        return out
    for i in range(n):                                   # coarse g by matching all reads
        sc = np.zeros(len(gg))
        for j in range(len(C)):
            pos = np.round(C[j] + gg).astype(int)
            pj = np.max([np.clip(SG[j] * Rb[i, np.clip(pos + d, 0, Rb.shape[1] - 1)], 0, 0.02) for d in range(-2, 3)], 0)
            sc += pj
        g[i] = gg[np.argmax(sc)]
    for it in range(iters):
        for i in range(1, n - 1):
            d = g[i + 1] - g[i - 1]; s[i] = v0 * np.sign(d) if abs(d) > 12 else 0
        for i in range(n):
            tp = taps(i, g[i], s[i], 7)
            if tp: g[i] = np.median([p - C[j] - s[i] * (p / FS - 0.058) for j, p, a in tp])
        nc = []
        for j in range(len(C)):
            r = [p - g[i] - s[i] * (p / FS - 0.058) for i in range(n) for jj, p, a in taps(i, g[i], s[i], 4) if jj == j]
            nc.append(np.median(r))
        C = np.array(nc)
    return g, s, C, taps


def ginterp_factory(t, g, s):
    tn = t + 0.058; n = len(t)
    def gi(tau):
        k = np.clip(np.searchsorted(tn, tau) - 1, 0, n - 2); h = tn[k + 1] - tn[k]; u = (tau - tn[k]) / h
        return (2 * u ** 3 - 3 * u ** 2 + 1) * g[k] + (u ** 3 - 2 * u ** 2 + u) * h * s[k] + (-2 * u ** 3 + 3 * u ** 2) * g[k + 1] + (u ** 3 - u ** 2) * h * s[k + 1]
    return gi


def stage_t2(args):
    X, m, t = load()
    Rb, Lb = base_residuals(X, m)
    Rb = Rb.astype(np.float64); Lb = Lb.astype(np.float64); n = len(Rb)
    base = make_base(m)
    sb = base(m[:, 1], 1, 9000)                           # base of the static mean (R)
    stat = {p: float(sb[p - 2:p + 3][np.argmax(np.abs(sb[p - 2:p + 3]))]) for p in (724, 1631, 2008)}
    out = {"static_first_reads_base_units": stat}
    bands = find_bands(Rb)
    C0 = [c for c, sg in bands]; SG = np.array([sg for c, sg in bands], float)
    g, s, C, taps = track_g(Rb, t, C0, SG)
    out["read_centres"] = [round(float(c), 2) for c in C]; out["read_sign_R"] = [int(v) for v in SG]
    # per-read: detection, amplitude, L sign / ratio / position offset
    rd = []
    for j in range(len(C)):
        amps, detn, la, dl = [], 0, [], []
        for i in range(n):
            tp = [(p, a) for jj, p, a in taps(i, g[i], s[i], 5) if jj == j]
            if not tp: continue
            p, a = tp[0]; detn += 1; amps.append(a)
            seg = Lb[i, p - 3:p + 4]; k = int(np.argmax(np.abs(seg))); la.append(seg[k] / (SG[j] * a) * SG[j]); dl.append(k - 3)
        la = np.array(la); amps = np.array(amps)
        rd.append({"c": round(float(C[j]), 1), "signR": int(SG[j]), "detect_frac": round(detn / n, 3), "amp_median": round(float(np.median(amps)), 4),
                   "L_over_R_signed_median": round(float(np.median(la)), 3), "L_pos_offset_mean": round(float(np.mean(dl)), 2)})
    out["reads"] = rd
    # held out: first 120 s vs last 120 s (centres, detection, amplitude re-estimated on each half from the same g)
    half = {}
    for nm, sel in (("first120", t < 120), ("last120", t >= 120)):
        hh = []
        for j in range(len(C)):
            cc, am = [], []
            for i in np.where(sel)[0]:
                for jj, p, a in taps(i, g[i], s[i], 4):
                    if jj == j: cc.append(p - g[i] - s[i] * (p / FS - 0.058)); am.append(a)
            hh.append([round(float(np.median(cc)), 2), round(float(np.median(am)), 4), round(len(cc) / sel.sum(), 3)])
        half[nm] = hh
    out["halves_centre_amp_detect"] = half
    # displacement g(t): range, speed, runs, phase wander at fixed f
    d = np.diff(g); inner = [abs(d[i]) for i in range(1, len(d) - 1) if np.sign(d[i - 1]) == np.sign(d[i]) == np.sign(d[i + 1]) and abs(d[i]) > 8]
    out["g"] = {"percentiles_1_5_50_95_99": [round(float(v), 1) for v in np.percentile(g, [1, 5, 50, 95, 99])],
                "interior_step_per_0.25s_mean_sd": [round(float(np.mean(inner)), 2), round(float(np.std(inner)), 2)],
                "speed_samples_per_s": round(float(np.mean(inner)) / 0.25, 1)}
    for nm, sel in (("first120", t < 120), ("last120", t >= 120)):
        dd = np.diff(g[sel]); inn = [abs(dd[i]) for i in range(1, len(dd) - 1) if np.sign(dd[i - 1]) == np.sign(dd[i]) == np.sign(dd[i + 1]) and abs(dd[i]) > 8]
        out["g"][f"speed_{nm}"] = round(float(np.mean(inn)) / 0.25, 1)
        out["g"][f"range_{nm}"] = [round(float(np.percentile(g[sel], 1)), 1), round(float(np.percentile(g[sel], 99)), 1)]
    tri = lambda th: 4 * np.abs((th % 1) - 0.5) - 1
    ph = []
    for seg in range(8):
        sel = (t >= 30 * seg) & (t < 30 * seg + 30); best = (-9, 0)
        for p0 in np.arange(0, 1, 0.01):
            b = tri(0.5535 * (t[sel] + 0.058) + p0); c = np.corrcoef(b, g[sel])[0, 1]
            if c > best[0]: best = (c, p0)
        ph.append([round(best[1], 2), round(best[0], 2)])
    out["g"]["phase_per_30s_at_0.5535Hz_[phase,corr]"] = ph
    gi = ginterp_factory(t, g, s)
    # c-domain stack and comb test (loops)
    cg = np.arange(300, 9000, 0.5)
    def stack(B, sel):
        acc = np.zeros(len(cg)); pp = np.arange(300, 9000)
        for i in np.where(sel)[0]:
            c = pp - gi(t[i] + pp / FS); idx = np.round((c - 300) / 0.5).astype(int); ok = (idx >= 0) & (idx < len(cg))
            np.add.at(acc, idx[ok], B[i, pp[ok]])
        return acc
    allsel = np.ones(n, bool)
    SR = stack(Rb, allsel); SL = stack(Lb, allsel)
    sm = lambda x: np.convolve(x, np.ones(3), "same")
    SRs = sm(SR)
    def val(c, w=2.0):
        k = (cg >= c - w) & (cg <= c + w); v = SRs[k]; return float(v[np.argmax(np.abs(v))])
    sigma = float(np.median(np.abs(SRs[(cg > 2700) & (cg < 3400)])) * 1.4826)
    # presence of each read in each snapshot (tap >= 0.006 within +-2 of the predicted place), its statistics
    pres = np.zeros((len(C), n), bool)
    for j in range(len(C)):
        for i in range(n):
            pr = C[j] + g[i] + s[i] * (C[j] / FS - 0.058); lo = int(round(pr)) - 2
            pres[j, i] = (SG[j] * Rb[i, lo:lo + 5]).max() > 0.006
    out["presence_fraction"] = [round(float(v), 3) for v in pres.mean(1)]
    out["presence_corr"] = np.round(np.corrcoef(pres.astype(float)), 2).tolist()
    out["presence_lag1_autocorr"] = [round(float(np.corrcoef(pres[j, :-1], pres[j, 1:])[0, 1]), 3) for j in range(len(C))]
    out["n_present_hist"] = np.bincount(pres.sum(0)).tolist()
    comb = []
    for j in range(len(C)):
        sel = pres[j]; Sj = sm(stack(Rb, sel)); SRs_keep = SRs
        def valj(c, w=2.0):
            k = (cg >= c - w) & (cg <= c + w); v = Sj[k]; return float(v[np.argmax(np.abs(v))]) if len(v) else float('nan')
        v0 = valj(C[j]); row = {"c": round(float(C[j]), 1), "n_snapshots": int(sel.sum()), "stack0": round(v0, 2)}
        sg_ = np.sign(v0)
        for k, P in enumerate(PL):
            row["loop" + "ABC"[k]] = [round(valj(C[j] + m * P) / v0, 3) for m in (1, 2, 3)]
        row["control_P+-120"] = [round(valj(C[j] + P + 120) / v0, 3) for P in PL]
        row["control_P-3"] = [round(valj(C[j] + P - 3.0 + 60) / v0, 3) for P in PL]
        comb.append(row)
    out["comb_conditional_on_presence"] = comb; out["stack_sigma_3bin_2700-3400"] = round(sigma, 3)
    # L/R sign map in the c domain, trained on one half, tested on the other
    def signmap(selA, selB):
        num = np.zeros(len(cg)); den = np.zeros(len(cg)); pp = np.arange(405, 8000)
        for i in np.where(selA)[0]:
            c = pp - gi(t[i] + pp / FS); idx = np.round((c - 300) / 0.5).astype(int)
            np.add.at(num, idx, Lb[i, pp] * Rb[i, pp]); np.add.at(den, idx, Rb[i, pp] ** 2)
        rho = np.where(den > 0, num / np.maximum(den, 1e-30), 0)
        rho = np.where(np.abs(rho) > 0.5, np.clip(rho, -1, 1), 0)         # keep only clear +-1 relations
        eL = eE = 0.0; eC = 0.0
        for i in np.where(selB)[0]:
            c = pp - gi(t[i] + pp / FS); idx = np.round((c - 300) / 0.5).astype(int)
            Lh = rho[idx] * Rb[i, pp]; eL += (Lb[i, pp] ** 2).sum(); eE += ((Lb[i, pp] - Lh) ** 2).sum()
            eC += ((Lb[i, pp] + Rb[i, pp]) ** 2).sum()
        return 1 - eE / eL, 1 - eC / eL, int((rho != 0).sum())
    f, l = t < 120, t >= 120
    r1 = signmap(f, l); r2 = signmap(l, f)
    out["LR_signmap_heldout_explained"] = {"train_first->test_last": [round(r1[0], 3), "plain -R copy:", round(r1[1], 3), "nonzero bins", r1[2]],
                                           "train_last->test_first": [round(r2[0], 3), "plain -R copy:", round(r2[1], 3), "nonzero bins", r2[2]]}
    out["stack_spikes_R"] = None
    SLs = sm(SL)
    y = SRs; a = np.abs(y); thr = 0.4
    pk = [i for i in range(3, len(y) - 3) if a[i] > thr and a[i] == a[i - 3:i + 4].max()]
    sp = []
    for i in pk:
        if sp and cg[i] - sp[-1][0] < 3: continue
        sp.append((float(cg[i]), round(float(y[i]), 2), round(float(SLs[i]), 2)))
    out["stack_spikes_c_R_L"] = sp
    np.save(os.path.join(UAD, "uad_moving_g.npy"), np.c_[t, g, s])
    json.dump(out, open(os.path.join(UAD, "uad_moving_t2.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "stack_spikes_c_R_L"}, indent=0)[:9000])
    print("spikes", len(sp))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--K", type=int, default=60)
    args = ap.parse_args()
    {"t1": stage_t1, "t2": stage_t2}[args.stage](args)
