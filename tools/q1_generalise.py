"""NEXT_STEPS Q1: does the program-S grey-box model (Path B) generalise across the Tomkin 96 kHz library?
Results: docs/model_s_results.md §14.

    .venv/bin/python tools/q1_generalise.py [--stage all|data|loops|laws|ap230|joint|ceiling|cla|report]
        [--splits rt,hf,lf] [--variant x0+Gx1|x] [--rounds 1] [--workers 10] [--sample 8] [--out analysis/q1]

Structure kept fixed (model_s_results.md §8-13): loops 1321/1389/1461/1587, the diffuser
comb460 . AP20 . AP72 . AP230 in series, the 8-tap cluster (gain x AMPL/90), per-channel injection
(excitations of P_i + 300 samples and a 1,400-tap first-pass FIR, shared by every setting).
What changes with the controls: the loop responses (control laws) and a stepped damping in AP230.

Stages:
  data     device-rate IRs (analysis/ml_structure/cache; with EMT_KERNEL=<name> cache_<name> -> cache_<name>/ and
           data_<name>.json here) re-aligned per channel so cluster tap 0 sits on
           index PRE; AMPL-0 captures (no cluster) are aligned to their AMPL-90 partner by cross-correlation.
  loops    per capture: loop RT, LF shelf, HF lowpass and clock ratio from the annihilator (loops x the
           diffuser's denominator), window from 0.30 s to where the natural decay ends.
  laws     control laws for the loops from the per-capture fits (Q4-Q7); refitted on training settings in `joint`.
  ap230    per setting: the AP230 loop-filter step b, from the first pass (before any loop recirculates).
  joint    per split (settings held out by RT, HF or LF step): one injection fitted on the training settings,
           loops from the laws, a clock ratio per capture; held-out nulls per window up to each capture's fade.
  ceiling  repeat-capture null: AMPL partners against each other with the same clock and gain fits.
  cla      cross-unit check: the CLA IR predicted from the library-trained injection and laws.
  report   held-out tables (stdout) and summary.json.
  check    spot check: re-derive --sample held-out captures per split from the saved injection and compare
           with the stored nulls (a few minutes; use instead of a full re-run to verify committed results).
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import ml_structure as ms  # noqa: E402
import model_s as m  # noqa: E402   (PRE, W, shelf, bl_kernel, CLUSTER: this module follows its conventions)
from core import subsample_peak  # noqa: E402

FS = ms.FSD
PRE = 96                       # device index of cluster tap 0 in the aligned arrays (as model_s.PRE)
CACHE = ms.tagged("cache")     # per-file arrays in <out>/<CACHE>/ and their index <out>/<DATA> (EMT_KERNEL, ml_structure)
DATA = ms.tagged("data") + ".json"
P = np.array([1321, 1389, 1461, 1587])


def fshift(x, s, n_out=None):
    """x delayed by s samples (band-limited FFT shift), along axis 0."""
    n = len(x) if n_out is None else n_out
    N = 1 << int(np.ceil(np.log2(max(len(x), n) + abs(s) + 256)))
    X = np.fft.rfft(x, N, axis=0)
    ph = np.exp(-2j * np.pi * np.fft.rfftfreq(N) * s)
    return np.fft.irfft(X * (ph[:, None] if x.ndim == 2 else ph), N, axis=0)[:n]


# ---------------------------------------------------------------- data
def _tap0(h):
    """Sub-sample position of cluster tap 0 per channel: the first peak above half the early maximum
    (taps 143 (R) and 54 (L) are as large as tap 0). Returns (R, L)."""
    out = []
    for c in (1, 0):
        e = np.abs(h[:300, c])
        i = int(np.argmax(e > 0.5 * e.max()))
        i = i + int(np.argmax(e[i:i + 5]))
        out.append(subsample_peak(h[:, c], i)[0])
    return tuple(out)


def _xalign(ref, x, a=0.03, b=0.4, maxlag=600):
    """Shift s (samples) such that fshift(x, s) matches ref over [a, b] s after tap 0 (one channel, 100 Hz-10 kHz)."""
    from scipy.optimize import minimize_scalar
    i0, i1 = PRE + int(a * FS), PRE + int(b * FS)
    r = sosfiltfilt(ms.BAND, ref[:i1 + maxlag])
    y = sosfiltfilt(ms.BAND, x[:i1 + maxlag])
    rw = np.zeros_like(r); rw[i0:i1] = r[i0:i1]
    N = 1 << int(np.ceil(np.log2(2 * len(r) + 2 * maxlag)))
    c = np.fft.irfft(np.fft.rfft(rw, N) * np.conj(np.fft.rfft(y, N)), N)       # c[k] = sum r[n] y[n-k]
    k = np.r_[np.arange(maxlag + 1), np.arange(-maxlag, 0)]
    best = int(k[np.argmax(np.r_[c[:maxlag + 1], c[-maxlag:]])])

    def cost(s):
        z = fshift(y, s, i1)
        return -np.sum(r[i0:i1] * z[i0:i1]) / np.sqrt(np.sum(z[i0:i1] ** 2))
    return float(minimize_scalar(cost, bounds=(best - 1, best + 1), method="bounded", options={"xatol": 1e-5}).x)


def data_stage(ir, out):
    os.makedirs(os.path.join(out, CACHE), exist_ok=True)
    groups = ms.ampl_groups("S")
    info = {}
    for key, d in groups.items():
        ref = ir(d[90])
        t0R, t0L = _tap0(ref)
        r = np.stack([fshift(ref[:, 0], PRE - t0L), fshift(ref[:, 1], PRE - t0R)], 1)
        rows = {}
        for amp, f in sorted(d.items()):
            if amp == 90:
                x, s = r, [PRE - t0L, PRE - t0R]
            else:
                h = np.vstack([np.zeros((PRE, 2)), ir(f)])
                s = [_xalign(r[:, c], h[:, c]) for c in range(2)]
                x = np.stack([fshift(h[:, c], s[c]) for c in range(2)], 1)
                if amp:                                          # check: the cluster tap lands on PRE
                    rows.setdefault("_tap0_check", {})[amp] = [round(float(subsample_peak(x[:, c], PRE)[0] - PRE), 4) for c in range(2)]
                s = [v + PRE for v in s]
            np.save(os.path.join(out, CACHE, os.path.basename(f)[:-4] + ".npy"), x.astype(np.float32))
            rows[amp] = {"file": os.path.basename(f), "shift_LR": [round(float(v), 4) for v in s], "n": len(x)}
        info["_".join(map(str, key))] = {"lf": key[0], "rt": key[1], "hf": key[2], "tap0_LminusR": round(t0L - t0R, 4), "files": rows}
    json.dump(info, open(os.path.join(out, DATA), "w"), indent=1)
    lr = [v["tap0_LminusR"] for v in info.values()]
    print(f"{len(info)} groups; L-R tap-0 offset after the cache's 19.19 us shift: median {np.median(lr):.4f}, range {min(lr):.4f}..{max(lr):.4f}")
    return info


def load(out, name, band=False):
    h = np.load(os.path.join(out, CACHE, name[:-4] + ".npy")).astype(np.float64)
    return sosfiltfilt(ms.BAND, h, axis=0) if band else h


# ---------------------------------------------------------------- per-file loop fits (annihilator)

T_IN = 0.30                    # the input (with the diffuser poles in the annihilator) has ended by here
LFBAND = butter(2, [20, 11000], btype="band", fs=FS, output="sos")
DIFF = dict(g460=0.445, g_ap=-0.7, ap=(20, 72, 230))


def diff_den(g460=DIFF["g460"], g=DIFF["g_ap"], w=m.W):
    """Denominator of comb460 . AP20 . AP72 . AP230: its poles are part of the characteristic polynomial."""
    z = lambda M: np.exp(-1j * w * M)
    d = 1 - g460 * z(460)
    for M in DIFF["ap"]:
        d = d * (1 - g * z(M))
    return d


F_HF = 7000.0                  # HF reference frequency of the parametrisation (the data sheet's 7 kHz)
PM = 1440.0                    # pole scaling reference length (about the mean loop)


def onepole_for(T, f=F_HF):
    """Pole p of (1 - p)/(1 - p z^-1) with |H| = T at f (T <= 1 gives 0 <= p < 1)."""
    c = np.cos(2 * np.pi * f / FS)
    T2 = np.minimum(T ** 2, 0.999999)
    return ((1 - T2 * c) - np.sqrt((1 - T2 * c) ** 2 - (1 - T2) ** 2)) / (1 - T2)


def G_loops(v, w=m.W, scaled_pole=False, hf_shared=False):
    """v = [rt0, ldc, lny, a0, hfr, eps_ppm]: per-pass response of each loop, every gain Jot-scaled.
    gain:  10^(-3 P_i / (fs rt0))                         (rt0 = the loops' RT at Nyquist before the HF lowpass)
    LF:    first-order shelf, DC gain giving RT ldc*rt0, Nyquist gain giving RT lny*rt0 (lny = 1 in the fits),
           pole a0 shared by the loops (scaled_pole: a_i = a0^(P_i/PM) instead; fits worse)
    HF:    one-pole lowpass, pole set so the RT at F_HF is hfr * (RT the loop has there without it)
    clock: delay P_i (1 + eps)"""
    rt0, ldc, lny, a0, hfr, eps = v[:6]
    out = []
    for Pi in P:
        k = -3 * Pi / (FS * rt0)
        g = 10 ** k
        a = np.sign(a0) * abs(a0) ** (Pi / PM) if scaled_pole else a0
        H = m.shelf(w, 10 ** (k * (1 / ldc - 1)), 10 ** (k * (1 / lny - 1)), a)
        kp = -3 * PM / (FS * rt0) if hf_shared else k   # hf_shared: one pole for every loop (set at the mean length)
        T = 10 ** (kp * (1 / hfr - 1)) if hfr <= 1 else 1.0
        p = onepole_for(np.array(T))
        out.append(g * H * (1 - p) / (1 - p * np.exp(-1j * w)) * np.exp(-1j * w * eps * 1e-6 * Pi))
    return out


def _annihilate(Y, Gs, w, N):
    A = np.ones(len(w), complex)
    for G, Pi in zip(Gs, P):
        A = A * (1 - G * np.exp(-1j * w * Pi))
    return np.fft.irfft(Y * A[:, None], N, axis=0)


def _blocks(e, y, i0, i1, blk):
    return np.array([10 * np.log10(np.sum(e[j:j + blk] ** 2) / np.sum(y[j:j + blk] ** 2)) for j in range(i0, i1 - blk + 1, blk)])


def fit_file(args):
    """Annihilator fit of the loops on one file (diffuser poles fixed), window [T_IN, t_end]; t_end where
    the block residual rises 10 dB above its early level (Tomkin's fade or the format floor).
    Returns the parameters, rough 1-sigma, the window and the in-window residual per channel."""
    from scipy.optimize import least_squares
    out, name, lab, opts = args
    h = sosfiltfilt(LFBAND, load(out, name), axis=0)
    N = 1 << int(np.ceil(np.log2(len(h) + 8192)))
    w = 2 * np.pi * np.fft.rfftfreq(N)
    Y = np.fft.rfft(h, N, axis=0) * diff_den(w=w)[:, None]
    yd = np.fft.irfft(Y, N, axis=0)
    n = np.arange(N) - PRE
    nf = int(0.025 * FS)
    E = np.array([np.mean(h[PRE + j:PRE + j + nf] ** 2) for j in range(0, len(h) - PRE - nf, nf)])
    E = 10 * np.log10(E / E[2:12].mean() + 1e-30)
    rt = lab["rt"]
    t_end = min(max(np.argmax(E < -30) * 0.025, T_IN + 0.15), len(h) / FS - 0.05)
    free = [0, 4, 5] + ([] if lab["lf"] == 1.0 and not opts.get("lf_free") else [1, 3])     # shelf Nyquist gain fixed at 1
    v = np.array([rt, lab["lf"], 1.0, 0.98, lab["hf"], 13.0])
    lo = np.array([0.05, 0.02, 0.2, -0.5, 0.05, -300])
    hi = np.array([20, 20, 5, 0.99999, 1.5, 300])
    xs = np.array([0.001 * rt, 0.01, 0.01, 0.001, 0.01, 1.0])
    sp = opts.get("scaled_pole", False)
    hs = opts.get("hf_shared", False)
    for it in range(2):
        sel = (n >= T_IN * FS) & (n < t_end * FS)
        wt = 10 ** (3 * (n[sel] / FS) / rt)[:, None]

        def res(u):
            vv = v.copy(); vv[free] = u
            return (_annihilate(Y, G_loops(vv, w, sp, hs), w, N)[sel] * wt).ravel()
        r = least_squares(res, v[free], bounds=(lo[free], hi[free]), x_scale=xs[free], diff_step=1e-7)
        v[free] = r.x
        if it == 0 and opts.get("t_end"):                    # fixed window (to compare filter forms)
            t_end = opts["t_end"]
        elif it == 0:                                        # window end: where the natural decay ends
            e = _annihilate(Y, G_loops(v, w, sp, hs), w, N)
            blk = int(0.05 * FS)
            b = _blocks(e[:, 1] + e[:, 0], yd[:, 1] + yd[:, 0], PRE + int(T_IN * FS), len(h) - blk, blk)
            base = np.median(b[:3])
            bad = np.nonzero(b > base + 10)[0]
            t_end = max(T_IN + (bad[0] if len(bad) else len(b)) * 0.05, T_IN + 0.15)
    e = _annihilate(Y, G_loops(v, w, sp, hs), w, N)[sel]
    rel = [float(10 * np.log10(np.sum(e[:, c] ** 2) / np.sum(yd[sel, c] ** 2))) for c in range(2)]
    J = r.jac
    s2 = np.sum(r.fun ** 2) / max(len(r.fun) - len(free), 1)
    sd = np.full(len(v), np.nan)
    try:
        sd[free] = np.sqrt(np.diag(s2 * 50 * np.linalg.inv(J.T @ J)))   # x50: residual correlation time (cf. model_s §9.1)
    except np.linalg.LinAlgError:
        pass
    return dict(file=name, labels=lab, v=v.tolist(), sd=sd.tolist(), free=free, t_end=float(t_end), rel_db=rel)


def response(v, f_hz, scaled_pole=False):
    """Per-loop |G| at f_hz and the implied RT (s) per loop."""
    G = np.abs(np.array(G_loops(v, 2 * np.pi * np.asarray(f_hz, float) / FS, scaled_pole)))
    return G, -60 * P[:, None] / (FS * 20 * np.log10(G))


def loops_stage(out, info):
    from concurrent.futures import ProcessPoolExecutor
    jobs = [(out, f["file"], dict(lf=g["lf"], rt=g["rt"], hf=g["hf"], amp=int(k)), {})
            for g in info.values() for k, f in g["files"].items() if not k.startswith("_")]
    with ProcessPoolExecutor(10) as ex:
        res = list(ex.map(fit_file, jobs, chunksize=2))
    json.dump(res, open(os.path.join(out, "loops.json"), "w"), indent=0)
    return res


# ---------------------------------------------------------------- control laws (from training files only)
LAW_FIT = dict(rt_min=1.4, rel_max=-35.0)        # per-file fits used for the laws: long tails, clean annihilator


def fit_laws(fits, train_keys=None, f_grid=np.arange(6000, 7201, 100)):
    """Laws from per-file annihilator fits: per LF step the shelf (Nyquist-end and DC RT over the label,
    one shared pole), one HF reference frequency, and per LF step the HF factor k
    (one-pole excess loss at f_ref = k (1/HF - 1) x the label's loss rate)."""
    use = [r for r in fits if r["labels"]["rt"] >= LAW_FIT["rt_min"] and np.mean(r["rel_db"]) < LAW_FIT["rel_max"]
           and (train_keys is None or _key(r["labels"]) in train_keys)]
    lfs = sorted({r["labels"]["lf"] for r in use})
    law = {"n_files": len(use), "lf": {}}
    for lf in lfs:
        u = [r for r in use if r["labels"]["lf"] == lf]
        hi = np.median([r["v"][0] / r["labels"]["rt"] for r in u])
        dc = np.median([r["v"][1] * r["v"][0] / r["labels"]["rt"] for r in u]) if lf != 1.0 else 1.0
        law["lf"][lf] = {"hi": float(hi), "dc": float(dc), "n": len(u)}
    law["a"] = float(np.median([r["v"][3] for r in use if r["labels"]["lf"] != 1.0]))

    def kvals(r, f):
        v, l = np.array(r["v"]), r["labels"]
        out = []
        for Pi in P:
            kk = -3 * Pi / (FS * v[0])
            p = onepole_for(np.array(10 ** (kk * (1 / v[4] - 1))))
            H = abs((1 - p) / (1 - p * np.exp(-2j * np.pi * f / FS)))
            out.append(-np.log10(H) * FS * l["rt"] / (3 * Pi * (1 / l["hf"] - 1)))
        return float(np.mean(out))
    hfu = [r for r in use if r["labels"]["hf"] < 1.0]
    spread = []
    for f in f_grid:                                  # f_ref: where k does not depend on the HF step (LF 1.0 files)
        med = [np.median([kvals(r, f) for r in hfu if r["labels"]["lf"] == 1.0 and r["labels"]["hf"] == hf])
               for hf in sorted({r["labels"]["hf"] for r in hfu})]
        spread.append(np.std(np.log(med)))
    law["f_ref"] = float(f_grid[int(np.argmin(spread))])
    for lf in lfs:
        law["lf"][lf]["k"] = float(np.median([kvals(r, law["f_ref"]) for r in hfu if r["labels"]["lf"] == lf]))
    law["n_params"] = 3 * len(lfs) - 1 + 2          # (hi, dc, k) per LF step (dc at LF 1.0 fixed), shared pole, f_ref
    return law


def _key(l):
    return f"{l['lf']}_{l['rt']}_{l['hf']}"


def G_law(lab, law, w, lf_interp=False):
    """Loop responses at a setting from the laws (no clock: data are re-gridded to the nominal device rate)."""
    rt, lf, hf = lab["rt"], lab["lf"], lab["hf"]
    if lf in law["lf"] and not lf_interp:
        L = law["lf"][lf]
    else:                                             # LF step not in the training set: interpolate in log LF
        ks = sorted(law["lf"]); x = np.log(ks)
        L = {n: float(np.exp(np.interp(np.log(lf), x, np.log([law["lf"][k][n] for k in ks])))) for n in ("hi", "dc", "k")}
    out = []
    for Pi in P:
        kk = -3 * Pi / FS
        g = 10 ** (kk / (L["hi"] * rt))
        H = m.shelf(w, 10 ** (kk * (1 / (L["dc"] * rt) - 1 / (L["hi"] * rt))), 1.0, law["a"]) if L["dc"] != L["hi"] else 1.0
        if hf < 1.0:
            p = onepole_for(np.array(10 ** (kk / rt * L["k"] * (1 / hf - 1))), law["f_ref"])
            H = H * (1 - p) / (1 - p * np.exp(-1j * w))
        out.append(g * H)
    return out


def laws_stage(out):
    fits = json.load(open(os.path.join(out, "loops.json")))
    law = fit_laws(fits)
    json.dump(law, open(os.path.join(out, "laws.json"), "w"), indent=1)
    print("laws (all files):", json.dumps(law))
    return law


# ---------------------------------------------------------------- the AP230 step (switched damping in the diffuser)
AP230_REF = "1.0_1.8_1.0"              # a setting whose AP230 has no damping (b = 0)


def ap230_resp(b, w, g=None):
    """AP230 with a two-tap loop filter (1 - b) + b z^-1 (b = 0: the plain all-pass of model_s §11)."""
    g = DIFF["g_ap"] if g is None else g
    z = np.exp(-1j * w * 230)
    F = (1 - b) + b * np.exp(-1j * w)
    return (-g + z * F) / (1 - g * z * F)


def ap230_stage(out, info):
    """Per setting: the AP230 loop-filter step b, from the first pass (240-1550 samples after tap 0, before any
    loop recirculates) of the AMPL-0 capture against the reference setting's: x = ref * AP230(b) / AP230(0)."""
    from scipy.optimize import minimize_scalar
    N = 1 << 14
    w = 2 * np.pi * np.fft.rfftfreq(N)
    ref = sosfiltfilt(LFBAND, load(out, info[AP230_REF]["files"]["0"]["file"]), axis=0)
    R = np.fft.rfft(ref[:N], axis=0)
    i0, i1 = PRE + 240, PRE + 1550
    res = {}
    for key, g in info.items():
        x = sosfiltfilt(LFBAND, load(out, g["files"]["0"]["file"]), axis=0)

        def cost(b):
            mdl = np.fft.irfft(R * (ap230_resp(b, w) / ap230_resp(0.0, w))[:, None], N, axis=0)
            return np.sum((x[i0:i1] - mdl[i0:i1]) ** 2) / np.sum(x[i0:i1] ** 2)
        c0 = cost(0.0)
        r = minimize_scalar(cost, bounds=(-0.1, 0.9), method="bounded", options={"xatol": 1e-5})
        res[key] = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"], b=float(r.x), null_db=float(10 * np.log10(r.fun)),
                        null_b0_db=float(10 * np.log10(c0)))
    bs = np.array([v["b"] for v in res.values()])
    levels = []
    for b in np.sort(bs):                               # the steps: cluster the fitted values (gaps > 0.05)
        if not levels or b - levels[-1][-1] > 0.05:
            levels.append([b])
        else:
            levels[-1].append(b)
    steps = [float(np.median(l)) for l in levels]
    for v in res.values():
        v["step"] = int(np.argmin([abs(v["b"] - s) for s in steps]))
    out_ = dict(steps=steps, settings=res)
    json.dump(out_, open(os.path.join(out, "ap230.json"), "w"), indent=0)
    print("AP230 steps:", np.round(steps, 4).tolist())
    for rt in sorted({v["rt"] for v in res.values()}):
        print(f"  RT {rt}: " + "  ".join(f"HF {hf}: " + "/".join(str(v["step"]) for v in sorted(
            [v for v in res.values() if v["rt"] == rt and v["hf"] == hf], key=lambda v: v["lf"])) for hf in sorted({v["hf"] for v in res.values()})))
    return out_


# ---------------------------------------------------------------- shared injection, held-out prediction
S0, EXT, LFIR = 230, 300, 1400          # excitation start, excitation length P_i + EXT, first-pass FIR length (model_s Cs)
T_CAP = 1.0                             # half-flattening fit weight 10^(1.5 t / RT), capped here
FADE_GUARD = 0.70                       # Tomkin's fades start about 0.65 s before the last non-zero sample
FLOOR_DB = 60.0                         # windows also end where the level is this far below the 50-300 ms level
WINDOWS = [(0, 50), (50, 300), (300, 600), (600, 1000), (1000, 1500), (1500, 2500)]


def warp(x, ppm, origin=PRE, half=128, beta=12.0):
    """x read at origin + (n - origin)(1 + ppm 1e-6), Kaiser-windowed sinc (re-grids a capture to the nominal device clock)."""
    n = np.arange(len(x))
    t = origin + (n - origin) * (1 + ppm * 1e-6)
    i = np.floor(t).astype(int)
    fr = t - i
    k = np.arange(-half + 1, half + 1)
    xp = np.concatenate([np.zeros((half,) + x.shape[1:]), x, np.zeros((half + 2,) + x.shape[1:])])
    out = np.zeros(x.shape)
    for s0 in range(0, len(x), 4096):
        d = fr[s0:s0 + 4096, None] - k[None, :]
        wgt = np.sinc(d) * np.i0(beta * np.sqrt(np.clip(1 - (d / (half + 0.5)) ** 2, 0, None))) / np.i0(beta)
        idx = np.clip(i[s0:s0 + 4096, None] + k[None, :] + half, 0, len(xp) - 1)
        out[s0:s0 + 4096] = np.einsum("nk,nk...->n...", wgt, xp[idx])
    return out


def fade_times(out, info):
    """Per file: end of the usable natural decay (s after tap 0) = min(start of the fade seen against the
    AMPL partners, last non-zero sample - FADE_GUARD (or 45 % of the file for short files), level FLOOR_DB down)."""
    nf = int(0.025 * FS)
    res = {}
    for g in info.values():
        H = {f["file"]: load(out, f["file"], band=True) for k, f in g["files"].items() if not k.startswith("_")}
        E, last = {}, {}
        for name, h in H.items():
            nz = np.nonzero(np.any(h != 0, 1))[0][-1]
            last[name] = (nz - PRE) / FS
            E[name] = 10 * np.log10(np.array([np.mean(h[PRE + j * nf:PRE + (j + 1) * nf] ** 2) for j in range((nz - PRE) // nf)]) + 1e-30)
        L = max(len(e) for e in E.values())
        M = np.full((len(E), L), -300.0)
        for i, e in enumerate(E.values()):
            M[i, :len(e)] = e
        for i, (name, e) in enumerate(E.items()):
            d = e - np.max(np.delete(M, i, 0), 0)[:len(e)]
            j = next((j for j in range(4, len(e) - 2) if np.all(d[j:j + 3] < -0.5)), None)
            lv = e - np.mean(e[2:12])                    # level re 50-300 ms; stop 60 dB down (short files reach the format floor)
            k = next((k for k in range(2, len(e)) if lv[k] < -FLOOR_DB), None)
            res[name] = float(min(last[name] - min(FADE_GUARD, 0.45 * last[name]), j * 0.025 if j is not None else np.inf,
                                  k * 0.025 if k is not None else np.inf))
    return res


def grid_n(rt):
    return 1 << 17 if rt > 1.4 else 1 << 16     # loop tail at the FFT length below about -100 dB


class Joint:
    """y_f = (AMPL_f / 90) cluster  +  B post (F + sum_i X_i R_i^s), per channel, for every file f of every
    setting s; X (excitations), F (first-pass FIR) and the 8 cluster gains are shared by all settings.
    R_i^s = 1 / (1 - G_i^s z^-P_i) from the laws; post = comb460 . AP20 . AP72 . AP230; B = analysis band.
    Data and model are compared on the nominal device grid (each capture re-gridded by its clock ratio)."""

    def __init__(self, groups, ncl=96):
        self.groups = groups                         # [{'R': [K spectra incl. post and band], 'PB': post*band, 'N': n, 'files': [{'y','w','amp'}]}]
        K = len(groups[0]["R"])                      # 4, or 8 when each loop also has an excitation through one pass of G_i
        self.Lx = [int(P[i % 4]) + EXT for i in range(K)]
        self.nx = sum(self.Lx)
        self.nmax = max(len(f["y"]) for g in groups for f in g["files"])
        kern = m.bl_kernel(ncl)[1]
        self.cl = []                                 # 8 cluster atoms in time (band-limited, band-filtered)
        for c in (0, 1):
            for t in m.CLUSTER[c]:
                a = np.zeros(self.nmax + 4 * ncl)
                a[PRE + t - ncl:PRE + t + ncl + 1] = kern
                self.cl.append(sosfiltfilt(LFBAND, a)[:self.nmax])
        self.cs = 100.0

    def size(self):
        return 8 + 2 * LFIR + 2 * self.nx

    def split(self, v):
        cl = v[:8]
        f = v[8:8 + 2 * LFIR].reshape(2, LFIR)
        xs, o = [], 8 + 2 * LFIR
        for c in (0, 1):
            row = []
            for L in self.Lx:
                row.append(v[o:o + L]); o += L
            xs.append(row)
        return cl, f, xs

    def _spectra(self, f, xs, N):
        F = np.zeros((2, N // 2 + 1), complex)
        X = np.zeros((2, len(self.Lx), N // 2 + 1), complex)
        for c in (0, 1):
            b = np.zeros(N); b[PRE + S0:PRE + S0 + LFIR] = f[c]; F[c] = np.fft.rfft(b)
            for i, (x, L) in enumerate(zip(xs[c], self.Lx)):
                b = np.zeros(N); b[PRE + S0:PRE + S0 + L] = x; X[c, i] = np.fft.rfft(b)
        return F, X

    def reverb(self, g, F, X):
        S = F * g["PB"][None] + np.einsum("cik,ik->ck", X, np.array(g["R"]))
        return np.fft.irfft(S, g["N"], axis=1).T

    def cluster(self, cl):
        out = np.zeros((self.nmax, 2))
        for j in range(8):
            out[:, j // 4] += self.cs * cl[j] * self.cl[j]
        return out

    def forward(self, v):
        cl, f, xs = self.split(v)
        C = self.cluster(cl)
        spec = {}
        res = []
        for g in self.groups:
            if g["N"] not in spec:
                spec[g["N"]] = self._spectra(f, xs, g["N"])
            r = self.reverb(g, *spec[g["N"]])
            for fl in g["files"]:
                n = len(fl["y"])
                res.append(((fl["amp"] / 90) * C[:n] + r[:n]) * fl["w"][:, None])
        return np.concatenate([x.T.ravel() for x in res])

    def adjoint(self, z):
        o = 0
        gcl = np.zeros(8)
        acc = {}
        for g in self.groups:
            N = g["N"]
            rs = np.zeros((N, 2))
            for fl in g["files"]:
                n = len(fl["y"])
                r = z[o:o + 2 * n].reshape(2, n).T * fl["w"][:, None]; o += 2 * n
                rs[:n] += r
                for j in range(8):
                    gcl[j] += (fl["amp"] / 90) * self.cs * (r[:, j // 4] @ self.cl[j][:n])
            Z = np.fft.rfft(rs, axis=0).T                 # (2, N/2+1)
            if N not in acc:
                acc[N] = [np.zeros((2, N // 2 + 1), complex), np.zeros((2, len(self.Lx), N // 2 + 1), complex)]
            acc[N][0] += Z * np.conj(g["PB"])[None]
            acc[N][1] += Z[:, None, :] * np.conj(np.array(g["R"]))[None]
        f = np.zeros((2, LFIR))
        xs = [[np.zeros(L) for L in self.Lx] for _ in (0, 1)]
        for N, (AF, AX) in acc.items():
            for c in (0, 1):
                f[c] += np.fft.irfft(AF[c], N)[PRE + S0:PRE + S0 + LFIR]
                for i, L in enumerate(self.Lx):
                    xs[c][i] += np.fft.irfft(AX[c, i], N)[PRE + S0:PRE + S0 + L]
        return np.concatenate([gcl, f.ravel()] + [x for row in xs for x in row])

    def target(self):
        return np.concatenate([(fl["y"] * fl["w"][:, None]).T.ravel() for g in self.groups for fl in g["files"]])

    def solve(self, iters=300, x0=None):
        from scipy.sparse.linalg import LinearOperator, lsqr
        op = LinearOperator((len(self.target()), self.size()), matvec=self.forward, rmatvec=self.adjoint, dtype=float)
        return lsqr(op, self.target(), iter_lim=iters, atol=1e-12, btol=1e-12, x0=x0)[0]

    def render(self, v, g, n, amp):
        """Model for one file of setting g (unweighted), n samples."""
        cl, f, xs = self.split(v)
        r = self.reverb(g, *self._spectra(f, xs, g["N"]))
        if n > self.nmax:                                # longer than every training capture: extend the cluster atoms
            self.cl = [np.r_[a, np.zeros(n - self.nmax)] for a in self.cl]
            self.nmax = n
        return (amp / 90) * self.cluster(cl)[:n] + r[:n]


def group_spectra(lab, law, N, Gs=None, through_g=False, b230=0.0):
    w = 2 * np.pi * np.fft.rfftfreq(N)
    z = lambda M: np.exp(-1j * w * M)
    post = 1 / (1 - DIFF["g460"] * z(460))
    for M in DIFF["ap"]:
        post = post * (ap230_resp(b230, w) if M == 230 else (-DIFF["g_ap"] + z(M)) / (1 - DIFF["g_ap"] * z(M)))
    from scipy.signal import sosfreqz
    B = np.abs(sosfreqz(LFBAND, worN=w)[1]) ** 2       # filtfilt = |H|^2, zero phase
    PB = post * B
    Gs = G_law(lab, law, w) if Gs is None else Gs
    R = [PB / (1 - G * z(Pi)) for G, Pi in zip(Gs, P)]
    if through_g == "only":                           # every excitation passes the loop gain/filter once
        R = [r * G for r, G in zip(R, Gs)]
    elif through_g == "gain":                         # second excitation through the loop gain only (DC of G_i)
        R += [r * np.abs(G[0]) for r, G in zip(R, Gs)]
    elif through_g == "gain+G":                       # second through the gain, third through gain and filter
        R += [r * np.abs(G[0]) for r, G in zip(R, Gs)] + [r * G for r, G in zip(R, Gs)]
    elif through_g == "filter":                       # second through the loop filter only (G_i / its DC gain)
        R += [r * G / np.abs(G[0]) for r, G in zip(R, Gs)]
    elif through_g == "all":                          # 1, gain, filter, gain and filter
        R += [r * np.abs(G[0]) for r, G in zip(R, Gs)] + [r * G / np.abs(G[0]) for r, G in zip(R, Gs)] + [r * G for r, G in zip(R, Gs)]
    elif through_g:                                   # second excitation per loop that passes the loop gain/filter once
        R += [r * G for r, G in zip(R, Gs)]
    return dict(N=N, PB=PB, R=R)


def file_eps(fits):
    """Clock ratio (ppm) per file from its annihilator fit where the tail is long and clean, else the median."""
    ok = {r["file"]: r["v"][5] for r in fits if r["labels"]["rt"] >= 1.0 and np.mean(r["rel_db"]) < -30}
    med = float(np.median(list(ok.values())))
    return {r["file"]: ok.get(r["file"], med) for r in fits}, med


def prep(out, name, rt, eps, t_end, weight=True, shift=0.0):
    """Re-gridded, band-limited capture up to t_end, and its fit weight (decay-flattened, capped at T_CAP, unit energy)."""
    y = warp(fshift(load(out, name), shift) if shift else load(out, name), eps)
    n = min(PRE + int(t_end * FS), len(y))
    y = sosfiltfilt(LFBAND, y, axis=0)[:n]
    t = (np.arange(n) - PRE) / FS
    w = 10 ** (3 * np.clip(t, 0, T_CAP) / rt / 2) if weight else np.ones(n)
    w = w / np.sqrt(np.sum((y * w[:, None]) ** 2))
    return y, w


def nulls(y, mdl, t_end, windows=WINDOWS, narrow=False):
    """Null (dB) per window and channel, windows clipped at t_end; None where less than 50 ms remain.
    narrow: compare in 100 Hz-10 kHz (the band of the repeat-capture ceiling, ml_structure.md §4.1) instead of 20 Hz-11 kHz."""
    if narrow:
        y, mdl = sosfiltfilt(ms.BAND, y, axis=0), sosfiltfilt(ms.BAND, mdl, axis=0)
    out = {}
    for a_ms, b_ms in windows:
        i0, i1 = PRE + int(a_ms / 1000 * FS), min(PRE + int(b_ms / 1000 * FS), PRE + int(t_end * FS), len(y))
        if i1 - i0 < int(0.05 * FS):
            out[f"{a_ms}-{b_ms}"] = None; continue
        out[f"{a_ms}-{b_ms}"] = [float(10 * np.log10(np.sum((y[i0:i1, c] - mdl[i0:i1, c]) ** 2) / np.sum(y[i0:i1, c] ** 2))) for c in (0, 1)]
    return out


def _clock_job(args):
    """Per-capture clock ratio: minimise the in-band null of the re-gridded capture against a fixed model
    over [50 ms, t_end] (one parameter per capture; the per-capture clock rule, NEXT_STEPS §2), without and
    with one gain per capture. AMPL-0 captures have no cluster tap to anchor them, so they also get a shift.
    Returns (eps, shift), (eps, shift, gain) for the gain fit."""
    from scipy.optimize import minimize
    yb, mdl, t_end, eps0, shift_free = args
    i0, i1 = PRE + int(0.05 * FS), min(PRE + int(t_end * FS), len(mdl), len(yb))
    x = mdl[i0:i1]

    def cost(p, gain=False):
        y = warp(fshift(yb[:i1 + 400], p[1]) if shift_free else yb[:i1 + 400], p[0], half=32, beta=8.0)[i0:i1]
        g = np.sum(y * x) / np.sum(x * x) if gain else 1.0
        return np.sum((y - g * x) ** 2) / np.sum(y ** 2), g
    res = []
    for gain in (False, True):
        best = min(np.arange(eps0 - 10, eps0 + 10.1, 2.0), key=lambda e: cost([e, 0.0], gain)[0])
        p0 = [best, 0.0]
        if shift_free:
            r = minimize(lambda p: cost(p, gain)[0], p0, method="Nelder-Mead",
                         options=dict(initial_simplex=[p0, [p0[0] + 1, 0], [p0[0], 0.01]], xatol=0.01, fatol=1e-9))
            p = r.x
        else:
            from scipy.optimize import minimize_scalar
            p = [minimize_scalar(lambda e: cost([e, 0.0], gain)[0], bounds=(best - 2, best + 2), method="bounded", options={"xatol": 0.02}).x, 0.0]
        res.append((float(p[0]), float(p[1])) + ((float(cost(p, True)[1]),) if gain else ()))
    return res


def apply_clock(yb, eps, shift):
    return warp(fshift(yb, shift) if shift else yb, eps)


def _prep_job(a):
    return prep(*a[:5], weight=a[5], shift=a[6])


TRAIN_AMPL = (30, 90)                   # captures used for training (the reverb is the same in all four; AMPL 0 has no cluster)


def run_split(out, info, fits, fade, name, train, test, variant="x0+Gx1", rounds=2, iters=300, workers=10, ap230=None):
    from concurrent.futures import ProcessPoolExecutor
    law = fit_laws(fits, set(train))
    b230 = {k: (ap230["steps"][v["step"]] if ap230 else 0.0) for k, v in (ap230 or {"settings": {}})["settings"].items()}
    eps0, med = file_eps(fits)
    eps = dict(eps0)
    shift = {}

    def spectra(key, own=None):
        g = info[key]
        lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
        N = grid_n(g["rt"])
        Gs = None
        if own is not None:
            w = 2 * np.pi * np.fft.rfftfreq(N)
            v = np.array(own); v[5] = 0.0
            Gs = G_loops(v, w)
        return group_spectra(lab, law, N, Gs=Gs, through_g={"x": False, "x0+Gx1": True}[variant], b230=b230.get(key, 0.0))

    def files_of(key, train_only=False):
        return [(int(k), f["file"]) for k, f in sorted(info[key]["files"].items(), key=lambda kv: kv[0])
                if not k.startswith("_") and (not train_only or int(k) in TRAIN_AMPL)]

    def band(name_):
        return sosfiltfilt(LFBAND, load(out, name_), axis=0)

    def build(keys):
        jobs = [(out, nm, info[k]["rt"], eps[nm], fade[nm], True, shift.get(nm, 0.0)) for k in keys for _, nm in files_of(k, True)]
        with ProcessPoolExecutor(workers) as ex:
            ys = dict(zip([j[1] for j in jobs], ex.map(_prep_job, jobs, chunksize=4)))
        gs = []
        for k in keys:
            G = spectra(k)
            G["key"] = k
            G["files"] = [dict(y=ys[nm][0], w=ys[nm][1], amp=amp, name=nm) for amp, nm in files_of(k, True)]
            gs.append(G)
        return gs

    hist = []
    v = None
    for rd in range(rounds):
        J = Joint(build(train))
        v = J.solve(iters, x0=v)
        jobs, names = [], []
        for G in J.groups:
            for fl in G["files"]:
                mdl = J.render(v, G, len(fl["y"]), fl["amp"])
                jobs.append((band(fl["name"]), mdl, fade[fl["name"]], eps[fl["name"]], fl["amp"] == 0))
                names.append(fl["name"])
        with ProcessPoolExecutor(workers) as ex:
            cl = list(ex.map(_clock_job, jobs, chunksize=4))
        new = [c[0][0] for c in cl]
        shift.update({nm: c[0][1] for nm, c in zip(names, cl)})
        d = np.array([new[i] - eps[nm] for i, nm in enumerate(names)])
        hist.append({"round": rd, "eps_change_ppm_median_abs": float(np.median(np.abs(d)))})
        eps.update(dict(zip(names, new)))
        print(f"[{name}] round {rd}: clock change median |d| {np.median(np.abs(d)):.2f} ppm", flush=True)
        del J
    J = Joint(build(train))
    v = J.solve(iters, x0=v)

    own = {r["file"]: r for r in fits}
    rows = []
    for k in train[::max(len(train) // 16, 1)] + test:
        held = k in test
        g = info[k]
        G = spectra(k)
        fl_ = files_of(k)
        renders = {}
        jobs = []
        for amp, nm in fl_:
            yb = band(nm)
            n = min(PRE + int(fade[nm] * FS), len(yb))
            mdl = J.render(v, G, n, amp)
            renders[nm] = mdl
            jobs.append((yb, mdl, fade[nm], eps[nm] if not held else eps0[nm], amp == 0))
        with ProcessPoolExecutor(workers) as ex:
            ce = list(ex.map(_clock_job, jobs))
        cap, capg = {}, {}
        for (amp, nm), ((e, sh), (eg, shg, gain)) in zip(fl_, ce):
            mdl = renders[nm]
            y = apply_clock(band(nm), e, sh)[:len(mdl)]
            yg = apply_clock(band(nm), eg, shg)[:len(mdl)] / gain        # one gain per capture (as the repeat-capture ceiling)
            cap[nm], capg[nm] = y, yg
            row = dict(key=k, held_out=held, lf=g["lf"], rt=g["rt"], hf=g["hf"], amp=amp, file=nm, eps_ppm=e, t_end=fade[nm],
                       shift=sh, eps_ppm_with_gain=eg, gain_db=float(20 * np.log10(gain)),
                       null=nulls(y, mdl, fade[nm]), null_100_10k=nulls(y, mdl, fade[nm], narrow=True),
                       null_gain=nulls(yg, mdl, fade[nm]), null_gain_100_10k=nulls(yg, mdl, fade[nm], narrow=True))
            o = own.get(nm)
            if held and o and g["rt"] >= 1.0 and np.mean(o["rel_db"]) < -30:      # loops from the file's own tail (no injection refit)
                mo = J.render(v, spectra(k, own=o["v"]), len(mdl), amp)
                row["null_own_loops"] = nulls(y, mo, fade[nm])
                row["null_own_loops_gain_100_10k"] = nulls(yg, mo, fade[nm], narrow=True)
            rows.append(row)
        names = [nm for _, nm in fl_]                  # repeat-capture ceiling: AMPL partners against each other
        for i, nm in enumerate(names):
            pn = names[(i + 1) % len(names)]
            te = min(fade[nm], fade[pn])
            r = next(r for r in rows if r["file"] == nm)
            pw = [w for w in WINDOWS if w[0] >= 50]
            for key_, C in (("", cap), ("_gain", capg)):
                n = min(len(C[nm]), len(C[pn]))
                r["null_partner" + key_] = nulls(C[nm][:n], C[pn][:n], te, windows=pw)
                r["null_partner" + key_ + "_100_10k"] = nulls(C[nm][:n], C[pn][:n], te, windows=pw, narrow=True)
        print(f"[{name}] {'HELD ' if held else 'train'} {k}: " + " ".join(
            f"{w}:{'--' if x is None else '%.1f/%.1f' % tuple(x)}" for w, x in rows[-1]["null"].items()), flush=True)
    npar = J.size()
    res = dict(split=name, variant=variant, train=train, test=test, law=law, rounds=hist, params_injection=npar, ap230=bool(ap230),
               params_law=law["n_params"], params_per_capture="clock ratio (1); per-channel sub-sample alignment (2)",
               rows=rows)
    np.save(os.path.join(out, f"v_{name}.npy"), v)
    json.dump(res, open(os.path.join(out, f"split_{name}.json"), "w"), indent=0)
    return res


SPLITS = {
    "rt": lambda l: l["rt"] in (0.8, 1.8),          # RT steps held out (interpolation)
    "hf": lambda l: l["hf"] == 0.4,                 # HF step held out: the HF law interpolates it
    "lf": lambda l: l["lf"] == 0.4,                 # LF step held out: the shelf tables interpolate it in log LF
}


def joint_stage(out, info, which, variant, use_ap230=True, rounds=1, workers=10):
    fits = json.load(open(os.path.join(out, "loops.json")))
    ap230 = json.load(open(os.path.join(out, "ap230.json"))) if use_ap230 else None
    fade = fade_times(out, info)
    if len(which) == 1 and os.path.exists(os.path.join(out, "fade.json")):
        fade = json.load(open(os.path.join(out, "fade.json")))   # parallel runs share one table
    else:
        json.dump(fade, open(os.path.join(out, "fade.json"), "w"), indent=0)
    for sp in which:
        test = [k for k, g in info.items() if SPLITS[sp](g)]
        train = [k for k in info if k not in test]
        run_split(out, info, fits, fade, f"{sp}_{variant}" + ("" if use_ap230 else "_noap230"), train, test, variant,
                  rounds=rounds, workers=workers, ap230=ap230)


# ---------------------------------------------------------------- cross-unit check: the CLA IR
CLA_DEV = os.path.join(ROOT, "analysis", "shortRT", "h_dev.npy")
CLA_LAB = dict(lf=1.0, rt=1.0, hf=0.9, amp=80)          # model_s_results.md §2 (panel display)


def cla_stage(out, split="rt_x0+Gx1", eps0=98.4):
    """Predict the CLA IR (another unit; HF 0.9 is not in the library) from the library-trained injection
    and laws. Fitted on CLA: clock ratio, a common shift and one gain (3 numbers), nothing digital."""
    R = json.load(open(os.path.join(out, f"split_{split}.json")))
    v = np.load(os.path.join(out, f"v_{split}.npy"))
    law = R["law"]
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    h = np.load(CLA_DEV)[m.K0 - PRE:m.K0 + int(1.9 * FS)]
    yb = sosfiltfilt(LFBAND, h, axis=0)
    n = len(yb) - 600
    G = group_spectra(CLA_LAB, law, grid_n(CLA_LAB["rt"]), through_g=R["variant"] == "x0+Gx1")
    J = Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=CLA_LAB["amp"])])])
    mdl = J.render(v, G, n, CLA_LAB["amp"])
    (e, sh), (eg, shg, gain) = _clock_job((yb, mdl, 1.0, eps0, True))
    y = apply_clock(yb, eg, shg)[:n] / gain
    res = {"split": split, "eps_ppm": eg, "shift": shg, "gain_db": float(20 * np.log10(gain)),
           "null_gain_20_11k": nulls(y, mdl, 1.5, windows=WINDOWS[:5] + [(400, 600)]),
           "null_gain_100_10k": nulls(y, mdl, 1.5, windows=WINDOWS[:5] + [(400, 600)], narrow=True)}
    json.dump(res, open(os.path.join(out, "cla.json"), "w"), indent=1)
    print("CLA cross-unit:", json.dumps(res))
    return res


# ---------------------------------------------------------------- spot check (sampled re-evaluation)
def check_stage(out, n=8, seed=0, workers=8, tol_db=0.2):
    """Re-derive a random sample of held-out rows of every split from the saved injection (v_*.npy), the split's
    laws and the AP230 table, and compare with the stored nulls. Minutes instead of a full re-run; it checks the
    evaluation path (rendering, clock fit, windows), not the training."""
    import glob
    from concurrent.futures import ProcessPoolExecutor
    info = json.load(open(os.path.join(out, DATA)))
    fade = json.load(open(os.path.join(out, "fade.json")))
    eps0, _ = file_eps(json.load(open(os.path.join(out, "loops.json"))))
    ap = json.load(open(os.path.join(out, "ap230.json")))
    rng = np.random.default_rng(seed)
    worst = 0.0
    for path in sorted(glob.glob(os.path.join(out, "split_*.json"))):
        R = json.load(open(path))
        name = os.path.basename(path)[6:-5]
        v = np.load(os.path.join(out, f"v_{name}.npy"))
        law = R["law"]
        law["lf"] = {float(k): x for k, x in law["lf"].items()}
        held = [r for r in R["rows"] if r["held_out"]]
        sample = [held[i] for i in rng.choice(len(held), min(n, len(held)), replace=False)]
        jobs, mdls = [], []
        for r in sample:
            g = info[r["key"]]
            lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
            G = group_spectra(lab, law, grid_n(g["rt"]), through_g={"x": False, "x0+Gx1": True}[R["variant"]],
                              b230=ap["steps"][ap["settings"][r["key"]]["step"]] if R.get("ap230") else 0.0)
            yb = sosfiltfilt(LFBAND, load(out, r["file"]), axis=0)
            nn = min(PRE + int(fade[r["file"]] * FS), len(yb))
            J = Joint([dict(G, files=[dict(y=np.zeros((nn, 2)), w=np.ones(nn), amp=r["amp"])])])
            mdl = J.render(v, G, nn, r["amp"])
            mdls.append((yb, mdl))
            jobs.append((yb, mdl, fade[r["file"]], eps0[r["file"]], r["amp"] == 0))
        with ProcessPoolExecutor(workers) as ex:
            ce = list(ex.map(_clock_job, jobs))
        for r, (yb, mdl), (_, (eg, shg, gain)) in zip(sample, mdls, ce):
            y = apply_clock(yb, eg, shg)[:len(mdl)] / gain
            new = nulls(y, mdl, fade[r["file"]], narrow=True)
            d = [abs(a_ - b_) for w in new if new[w] and r["null_gain_100_10k"].get(w)
                 for a_, b_ in zip(new[w], r["null_gain_100_10k"][w])]
            worst = max(worst, max(d))
            print(f"[{name}] {r['file']}: 50-300 ms {new['50-300'][0]:.1f}/{new['50-300'][1]:.1f} "
                  f"(stored {r['null_gain_100_10k']['50-300'][0]:.1f}/{r['null_gain_100_10k']['50-300'][1]:.1f}); "
                  f"max |diff| over windows {max(d):.2f} dB", flush=True)
    print(f"spot check: worst |diff| {worst:.2f} dB ({'PASS' if worst <= tol_db else 'FAIL'} at {tol_db} dB)")
    return worst


# ---------------------------------------------------------------- repeat-capture ceiling (AMPL partners)
def _pair_job(a):
    out, fa, fb, ea, t_end = a
    ya = warp(sosfiltfilt(LFBAND, load(out, fa), axis=0), ea)
    yb = sosfiltfilt(LFBAND, load(out, fb), axis=0)
    n = min(len(ya), len(yb))
    ya = ya[:n]
    _, (eg, _, g) = _clock_job((yb[:n], ya, t_end, ea, False))
    y = apply_clock(yb[:n], eg, 0.0) / g
    return dict(a=fa, b=fb, eps_rel_ppm=eg - ea, null_100_10k=nulls(ya, y, t_end, windows=WINDOWS[1:], narrow=True),
                null_20_11k=nulls(ya, y, t_end, windows=WINDOWS[1:]))


def ceiling_stage(out, info, workers=8):
    """Repeat-capture null: AMPL 30/60/90 captures of each setting against each other, the second re-gridded to the
    first with its own clock ratio and one gain (the same two numbers a model null gets). Windows as for the model."""
    from concurrent.futures import ProcessPoolExecutor
    fade = json.load(open(os.path.join(out, "fade.json")))
    eps, _ = file_eps(json.load(open(os.path.join(out, "loops.json"))))
    jobs = []
    for g in info.values():
        f = {int(k): v["file"] for k, v in g["files"].items() if not k.startswith("_")}
        for a_, b_ in ((30, 60), (60, 90), (90, 30)):
            jobs.append((out, f[a_], f[b_], eps[f[a_]], min(fade[f[a_]], fade[f[b_]])))
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_pair_job, jobs, chunksize=4))
    json.dump(res, open(os.path.join(out, "repeat_ceiling.json"), "w"), indent=0)
    return res


# ---------------------------------------------------------------- report
def _med(rows, field, w, c=None):
    v = [r[field][w] for r in rows if field in r and r[field].get(w) is not None]
    if not v:
        return None
    v = np.array(v)
    return np.median(v, 0) if c is None else float(np.median(v[:, c]))


def report(out):
    """Held-out tables per split (markdown to stdout, numbers to summary.json)."""
    import glob
    summ = {}
    cpath = os.path.join(out, "repeat_ceiling.json")
    ceil = {r["a"]: r for r in json.load(open(cpath))} if os.path.exists(cpath) else {}
    for path in sorted(glob.glob(os.path.join(out, "split_*.json"))):
        R = json.load(open(path))
        rows = [r for r in R["rows"] if r["held_out"]]
        tr = [r for r in R["rows"] if not r["held_out"]]
        ws = [f"{a}-{b}" for a, b in WINDOWS]
        S = {"variant": R["variant"], "params_injection": R["params_injection"], "params_law": R["params_law"],
             "held_out_files": len(rows), "train_groups": len(R["train"]), "test_groups": len(R["test"]), "law": R["law"], "windows": {}}
        print(f"\n## {R['split']}  ({R['variant']}; injection {R['params_injection']} + laws {R['params_law']} parameters;"
              f" {len(R['train'])} training settings, {len(R['test'])} held out = {len(rows)} captures)")
        print("| window (ms) | n | held out, 100 Hz-10 kHz, +gain | no gain | 20 Hz-11 kHz, +gain | IQR (R, first col.) | 16 training settings | own loops | repeat-capture | excess over repeat (R) |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for w in ws:
            n = sum(1 for r in rows if r["null_gain_100_10k"].get(w) is not None)
            if not n:
                continue
            f = lambda x: "--" if x is None else f"{x[0]:.1f} / {x[1]:.1f}"
            a1 = _med(rows, "null_gain_100_10k", w)
            a2 = _med(rows, "null_100_10k", w)
            a3 = _med(rows, "null_gain", w)
            iq = np.percentile([r["null_gain_100_10k"][w][1] for r in rows if r["null_gain_100_10k"].get(w)], [25, 75])
            ins = _med(tr, "null_gain_100_10k", w)
            own = _med(rows, "null_own_loops_gain_100_10k", w)
            cr = [ceil[r["file"]] for r in rows if r["file"] in ceil and ceil[r["file"]]["null_100_10k"].get(w)]
            pr = np.median([c["null_100_10k"][w] for c in cr], 0) if cr else None
            ex = [r["null_gain_100_10k"][w][1] - ceil[r["file"]]["null_100_10k"][w][1] for r in rows
                  if r["file"] in ceil and r["null_gain_100_10k"].get(w) and ceil[r["file"]]["null_100_10k"].get(w)]
            exs = f"{np.median(ex):+.1f}" if ex else "--"
            print(f"| {w} | {n} | {f(a1)} | {f(a2)} | {f(a3)} | {iq[0]:.1f}…{iq[1]:.1f} | {f(ins)} | {f(own)} | {f(pr)} | {exs} |")
            S["windows"][w] = {"n": n, "held_out_gain_100_10k": a1.tolist(), "held_out_100_10k": a2.tolist(), "held_out_gain_20_11k": a3.tolist(),
                               "iqr_R": iq.tolist(), "in_sample": None if ins is None else ins.tolist(),
                               "own_loops": None if own is None else own.tolist(), "repeat": None if pr is None else pr.tolist(),
                               "excess_R": float(np.median(ex)) if ex else None}
        for ctrl in ("rt", "hf", "lf", "amp"):
            vals = sorted({r[ctrl] for r in rows})
            if len(vals) < 2:
                continue
            line = []
            for v_ in vals:
                rr = [r for r in rows if r[ctrl] == v_]
                x = _med(rr, "null_gain_100_10k", "50-300")
                y = _med(rr, "null_gain_100_10k", "300-600")
                line.append(f"{ctrl} {v_}: {x[0]:.1f}/{x[1]:.1f}" + (f", {y[0]:.1f}/{y[1]:.1f}" if y is not None else ""))
            print(f"- by {ctrl} (50-300 ms, 300-600 ms; L/R): " + "; ".join(line))
        g = np.array([r["gain_db"] for r in rows])
        print(f"- per-capture gain: median {np.median(g):+.3f} dB, IQR {np.percentile(g, 25):+.3f}…{np.percentile(g, 75):+.3f}; by AMPL " +
              ", ".join(f"{a_}: {np.median([r['gain_db'] for r in rows if r['amp'] == a_]):+.3f}" for a_ in sorted({r['amp'] for r in rows})))
        S["gain_db_by_amp"] = {a_: float(np.median([r["gain_db"] for r in rows if r["amp"] == a_])) for a_ in sorted({r["amp"] for r in rows})}
        summ[R["split"]] = S
    # the same capture held out in one split and in the training set of another
    import itertools
    S = {path: json.load(open(path)) for path in sorted(glob.glob(os.path.join(out, "split_*.json")))}
    rows = {k: {r["file"]: r for r in v["rows"]} for k, v in S.items()}
    pairs = [(rows[a_][f], rows[b_][f]) for a_, b_ in itertools.permutations(S, 2)
             for f, r in rows[a_].items() if r["held_out"] and f in rows[b_] and not rows[b_][f]["held_out"]]
    summ["paired"] = {"captures": len(pairs), "windows": {}}
    print(f"\n## Paired: {len(pairs)} captures held out in one split and in the training set of another (100 Hz-10 kHz, +gain)")
    print("| window (ms) | n | held out | in-sample | held out − in-sample (median of pairs) |\n|---|---|---|---|---|")
    for w in [f"{a}-{b}" for a, b in WINDOWS]:
        d = [(p[0]["null_gain_100_10k"][w], p[1]["null_gain_100_10k"][w]) for p in pairs
             if p[0]["null_gain_100_10k"].get(w) and p[1]["null_gain_100_10k"].get(w)]
        if d:
            h, t = np.array([x[0] for x in d]), np.array([x[1] for x in d])
            print(f"| {w} | {len(d)} | {np.median(h[:, 0]):.1f} / {np.median(h[:, 1]):.1f} | {np.median(t[:, 0]):.1f} / {np.median(t[:, 1]):.1f} |"
                  f" {np.median(h[:, 0] - t[:, 0]):+.1f} / {np.median(h[:, 1] - t[:, 1]):+.1f} |")
            summ["paired"]["windows"][w] = {"n": len(d), "held_out": np.median(h, 0).tolist(), "in_sample": np.median(t, 0).tolist(),
                                           "diff": np.median(h - t, 0).tolist()}
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    return summ


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=["all", "data", "loops", "laws", "ap230", "joint", "ceiling", "report", "cla", "check"])
    ap.add_argument("--no-ap230", action="store_true", help="joint: leave the AP230 step out (b = 0 everywhere)")
    ap.add_argument("--rounds", type=int, default=1, help="joint: clock-refit rounds before the final fit")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--sample", type=int, default=8, help="check: held-out captures re-derived per split")
    ap.add_argument("--out", default=os.path.join(ROOT, "analysis", "q1"))
    ap.add_argument("--splits", default="rt,hf,lf")
    ap.add_argument("--variant", default="x0+Gx1", choices=["x", "x0+Gx1"])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ir = ms.IRs(ms.CACHE)
    if a.stage in ("all", "data"):
        data_stage(ir, a.out)
    info = json.load(open(os.path.join(a.out, DATA)))
    if a.stage in ("all", "loops"):
        loops_stage(a.out, info)
    if a.stage in ("all", "laws"):
        laws_stage(a.out)
    if a.stage in ("all", "ap230"):
        ap230_stage(a.out, info)
    if a.stage in ("all", "joint"):
        joint_stage(a.out, info, a.splits.split(","), a.variant, not a.no_ap230, a.rounds, a.workers)
    if a.stage == "check":                            # not part of "all": a quick spot check of committed results
        check_stage(a.out, a.sample)
    if a.stage in ("all", "ceiling"):
        ceiling_stage(a.out, info)
    if a.stage in ("all", "cla"):
        cla_stage(a.out)
    if a.stage in ("all", "report"):
        report(a.out)


if __name__ == "__main__":
    main()
