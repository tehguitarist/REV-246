"""Derivative columns of candidate model perturbations (frequency-domain component builder, equal to nl_model.render to <1e-9) and the per-capture
sufficient statistics (X'X, X'e, e'e, y'y) of the regression e = y - m ~ X beta, per channel and band.  X includes the capture's channel gain."""
from cr_common import *
from scipy.signal import lfilter
BK_RAW = {int(c): np.array(v) for c, v in json.load(open(M.BK_FILE)).items()}
BANDS3 = [(20, 100), (100, 1000), (1000, 5000)]
SOSB = {b: band_sos(*b) for b in BANDS3}
NT = 27


def comps(key, ch, N, r):
    """dict of complex spectra (N//2+1), each already multiplied by the capture's time shift, processing response and gain."""
    span, early, late = key
    assert late == 0 or span == 0.5
    ns = ne.N_STAGES[span]
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    sh = np.exp(ph * (K0 + r["front_shift"])) * M._proc(N)[:, ch] * r["gain"][ch]
    fdb, sdb = ne.law_db(span, early, ch)
    fr_ = [ne.FRONT_AMP[ch][i] * 10 ** (fdb[i] / 20) * np.exp(ph * f) for i, f in enumerate(ne.FRONT[ch])]
    d = ne.diffuser(ch, N // 2)
    Df = np.fft.rfft(d, N)
    kd = M.BK[ch]; h = (len(kd) - 1) // 2
    K = np.fft.rfft(np.roll(np.r_[kd, np.zeros(N - len(kd))], -h))
    k1 = BK_RAW[ch]; dl = np.zeros(len(k1)); dl[(len(k1) - 1) // 2] = 1.0
    K1m = np.fft.rfft(np.roll(np.r_[k1 - dl, np.zeros(N - len(k1))], -h))          # d K / d alpha
    taps = {}
    rdir = np.zeros(len(ph), complex)
    for t, s, i in ne.taps(ch, ns):
        g = s * 10 ** (sdb[[j for j, st in enumerate(ne.STAGES) if i in st][0]] / 20)
        if ch == 1 and i == 0:
            taps[i] = (g * (Df - d[0]) * np.exp(ph * t), t)
            rdir = g * d[0] * np.exp(ph * ne.R_FIRST_DIRECT)
        else:
            taps[i] = (g * Df * np.exp(ph * t), t)
    ob = np.exp(ph * r["o"])
    return dict(ph=ph, sh=sh, fr=fr_, taps=taps, rdir=rdir, K=K, K1m=K1m, ob=ob, Df=Df, d=d, ns=ns, span=span, early=early, ch=ch)


def ap_deriv(ch, N, which):
    """d Df / d g for the AP(124) (which=0) or the output AP (180 / 188) (which=1) coefficient, central difference."""
    out = []
    for dg in (+1e-5, -1e-5):
        g1 = ne.G_AP + (dg if which == 0 else 0); g2 = ne.G_AP + (dg if which == 1 else 0)
        dd = np.zeros(N // 2); dd[0] = 1.0
        dd = ne.ap(ne.ap(dd, ne.AP_IN, g1), ne.AP_OUT[ch], g2)
        out.append(np.fft.rfft(dd, N))
    return (out[0] - out[1]) / 2e-5


def model_spec(C):
    """Baseline model spectrum from the components."""
    br = sum(v[0] for v in C["taps"].values()) if C["taps"] else 0
    return (sum(C["fr"]) + C["rdir"] * C["K"] + br * C["K"] * C["ob"])


def build(key, r, n, a, b):
    """Return per channel: names, X (n, P) full band time series (before the band filter, segment not yet taken), baseline model (n,)."""
    N = 1 << int(np.ceil(np.log2(n + 2 * ne.TAPS[-1] + 4096)))
    out = []
    ramp = np.zeros(n); ramp[a:b] = np.linspace(-1, 1, b - a)          # window ramp (the segment is what gets scored)
    for ch in (0, 1):
        C = comps(key, ch, N, r); ph = C["ph"]; sh = C["sh"]; K = C["K"]; ob = C["ob"]
        irf = lambda S: np.fft.irfft(S * sh, N)[:n]
        mod = irf(model_spec(C))
        cols, names = [], []
        def add(nm, S, rampit=False):
            x = irf(S)
            if rampit: x = x * ramp
            cols.append(x); names.append(nm)
        taps = C["taps"]
        bsum = sum(v[0] for v in taps.values()) if taps else np.zeros(len(ph), complex)
        Btot = bsum * K * ob                                   # branch part (taps through K, shifted by o)
        Ftot = sum(C["fr"]) + C["rdir"] * K                    # front part (+ R direct term through K)
        # nuisance
        add("nuis_gain", model_spec(C))
        add("nuis_shift", ph * model_spec(C))
        add("nuis_ppm", ph * model_spec(C), True)
        # branch / front
        add("branch_dt", ph * Btot); add("branch_gain", Btot); add("front_dt", ph * Ftot); add("front_gain", Ftot)
        add("rdir_dt", ph * C["rdir"] * K if ch else 0 * ph); 
        add("branch_dt_ramp", ph * Btot, True); add("branch_gain_ramp", Btot, True)
        for j in range(3):
            add(f"fr{j}_gain", C["fr"][j]); add(f"fr{j}_dt", ph * C["fr"][j])
        for i in range(NT):
            if i in taps: add(f"tap{i}_gain", taps[i][0] * K * ob); add(f"tap{i}_dt", ph * taps[i][0] * K * ob)
            else: add(f"tap{i}_gain", 0 * ph); add(f"tap{i}_dt", 0 * ph)
        add("bk_alpha_all", bsum * C["K1m"] * ob)
        for i in range(NT):
            add(f"bk_tap{i}", (taps[i][0] * C["K1m"] * ob) if i in taps else 0 * ph)
        for j in range(-14, 15):
            add(f"shape{j:+d}", Btot * np.exp(ph * j))
        for w_, nm in ((0, "ap124"), (1, "apout")):
            dD = ap_deriv(ch, N, w_)
            S = 0
            for i, (c, t) in taps.items():
                S = S + c * (dD / C["Df"]) if ch == 0 or i != 0 else S + c * (dD - (dD[0] if False else 0)) / C["Df"]
            add(nm, S * K * ob)
        out.append((names, np.array(cols).T, mod))
    return out


def stats_for(key, r):
    y, m = cached(key, r)
    a, b = window(key, r)
    n = NCACHE
    res = build(key, r, n, a, b)
    st = {}
    for ch, (names, X, mod) in enumerate(res):
        if ch == 0: st["names"] = np.array(names)
        st[f"base_err{ch}"] = float(np.max(np.abs(mod - m[:, ch])[a:b]) / np.max(np.abs(m[:, ch])))
        for bd in BANDS3:
            sos = SOSB[bd]
            Xb = sosfiltfilt(sos, X, axis=0)[a:b]
            eb = sosfiltfilt(sos, y[:, ch] - m[:, ch], axis=0)[a:b]
            yb = sosfiltfilt(sos, y[:, ch], axis=0)[a:b]
            st[f"XtX_{ch}_{bd[0]}"] = Xb.T @ Xb; st[f"Xte_{ch}_{bd[0]}"] = Xb.T @ eb
            st[f"ete_{ch}_{bd[0]}"] = float(eb @ eb); st[f"yty_{ch}_{bd[0]}"] = float(yb @ yb)
    return st


if __name__ == "__main__":
    B = bench_all()
    EAR = (0.7, 1.0, 1.5, 2.0, 3.0, 10.0)
    sel = [k for k in sorted(B) if k[1] in EAR and not B[k]["trimmed"] and (k[0] == 0.5 or (k[0] in (0.2, 0.3, 0.4) and k[2] == 0.0))]
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else len(sel)
    sel = sel[::max(1, len(sel) // lim)][:lim] if lim < len(sel) else sel
    os.makedirs(os.path.join(HERE, "stats2"), exist_ok=True)
    t = time.time()
    for k in sel:
        f = os.path.join(HERE, "stats2", tag(k) + ".npz")
        if os.path.exists(f): continue
        st = stats_for(k, B[k])
        np.savez(f, **st)
        print(k, f"base err {st['base_err0']:.1e} {st['base_err1']:.1e}", f"{time.time()-t:.0f}s", flush=True)
