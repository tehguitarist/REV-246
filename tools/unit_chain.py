"""NEXT_STEPS §0 accuracy audit item 4: are the cross-unit floors (S -35...-38 dB, M / L first pass -35...-36 dB against CLA
and the Altiverb) the units' capture chains or the algorithm? A compact per-unit chain model, counted and fitted on one
window / preset only, with the digital model fixed at the Tomkin-fitted baseline. Results: docs/model_s_results.md §37.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 .venv/bin/python tools/memguard.py --limit-gb 9 -- \\
        .venv/bin/python tools/unit_chain.py --stage s|ml|diag|report

  s     CLA (LF 1.0 / RT 1.0 / HF 0.9, AMPL 80) and the five Altiverb S presets (LF 1.0 / HF 0.9, RT from the name) against the
        S baseline (`_pshoulder`: fold_1.0's constants, the (LF 1.0, RT) shelf cells where RT is a library step, every
        loop's HF pole from §32.6's pole curve at the law's HF 0.9 pole; law elsewhere), and against `units` (the model
        pulse_check --stage units renders: fold_all, the ramp, no cells) for reference. Chain variants (per-unit
        numbers in brackets):
          fir        §26.2: clock, shift, gain on 50 ms-1.5 s, 7-tap FIR per channel on the cluster-only region [17]
          fir+bq     the same after the other unit's 10 low-band biquad numbers (bq.json; held out by unit) [17 (+10 other)]
          chain@W    the compact chain fitted on window W only: clock, shift, a 7-tap FIR per channel (gain, delay, HF tilt)
                     and the 10-number low-band cascade (§27's form) [26]; W = 0-50 (scored on 50-600 ms), 150-300 (scored
                     on 0-150 and 300-600), all (0-600, in-sample)
          np@W       the LTI ceiling: clock / shift of chain@W, then a free complex transfer per channel and 1/12 octave
                     (40 Hz-11 kHz) measured on W [~4 x 100]
          rt1.6      (Altiverb only) chain@all fitted on the RT 1.6 preset, applied to the other presets with only a
                     shift and gain per preset fitted on 0-50 ms [2 per preset; 26 shared]
        Nulls 100 Hz-10 kHz per window (0-50 / 50-150 / 150-300 / 300-600 ms) and per octave (125 Hz-8 kHz), L and R; the
        noise-limited null per window and octave from each file's tail (CLA 1.6-2.0 s; the Altiverb: the RT 0.6 / 1.0
        presets' 1.5-1.9 s tails, the same absolute noise for every preset) -> analysis/unit_chain/s.json.
  ml    The Altiverb MED / LONG presets (altiverb_ml.py's fit: per preset and channel an 11-tap FIR on the cluster-only
        region, cluster and reverb gains, clock refined on recirculation 1; AP damping step from S's four values on the
        first pass) with the model parts first through: none (§17), the CLA biquads (held out by unit), the Altiverb RT
        1.6 S preset's chain@all low-band cascade (held out by program), or its whole chain (FIR + cascade, a shift and the
        two gains per preset and channel only) -> analysis/unit_chain/ml.json.
  diag  The chain@0-50 null at HF 0.80-1.00 (CLA, the Altiverb RT 1.6; fine steps 0.88-0.92) and RT 0.97-1.03 (CLA): a
        sensitivity read of the setting mismatch, not a fit (hf_fine.json: the fine HF steps, first run as a separate scan), and the correlation of CLA's residual with the model's residual on Tomkin's nearest setting
        (1.0_1.0_1.0) -> analysis/unit_chain/diag.json.
  report tables -> analysis/unit_chain/report.log
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import pulse_check as pc  # noqa: E402
import unit_transfer as ut  # noqa: E402

q1, q2, lf_ = pc.q1, pc.q2, pc.lf_
FS, PRE = pc.FS, pc.PRE
OUT = os.path.join(ROOT, "analysis", "unit_chain")
CV = pc.CV
WIN = [(0, 50), (50, 150), (150, 300), (300, 600)]
OCT = [125, 250, 500, 1000, 2000, 4000, 8000]
OCT_SOS = {f: butter(3, [f / np.sqrt(2), min(f * np.sqrt(2), 11000)], btype="band", fs=FS, output="sos") for f in OCT}
BAND = q1.ms.BAND                                               # 100 Hz-10 kHz, as q1.nulls(narrow=True)
TWELFTH = 40 * 2 ** (np.arange(0, 101) / 12)                    # 40 Hz-11.3 kHz
KF = 3                                                          # 7-tap FIR per channel
FIT_W = {"0-50": [(0, 50)], "150-300": [(150, 300)], "all": [(0, 600)]}
LIB_RT = (0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5)


# ---------------------------------------------------------------- the digital model (fixed: Tomkin-fitted)
def _law():
    q1.EXT = 1700
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    ap = json.load(open(os.path.join(q2.OUT, lf_.AP230)))
    return v, law, ap


def pole_curve():
    p = json.load(open(os.path.join(ROOT, "analysis", "lowband", "chain3_pole_cells_allf.json")))
    kn, th, fold = np.array(p["knots"]), np.array(p["all"]), float(p["fold"])
    return lambda x: float(np.interp(x, kn, th)) if x < fold else 0.8465


def law_poles(lab, law):
    """Each loop's law HF pole at lab, with the RT x HF floor (as q2.G_capped)."""
    L = {float(k): x for k, x in law["lf"].items()}[lab["lf"]]
    hf_e = max(lab["hf"], q2.HFRT_MIN / lab["rt"]) if q2.HFRT_MIN else lab["hf"]
    return [float(q1.onepole_for(np.array(10 ** (-3 * Pi / FS / lab["rt"] * L["k"] * (1 / hf_e - 1))), law["f_ref"])) for Pi in q1.P]


def model_G(lab, mode, law, ap):
    """mode base: the _pshoulder baseline at lab (fold_1.0 constants, which don't act at LF 1.0 / HF 0.9; the (LF, RT) shelf
    cells of the HF 1.0 setting where RT is a library step; the pole curve at the law's pole); units: pulse_check's units
    model (fold_all, the ramp's median, no cells)."""
    if mode == "units":
        fold = json.load(open(os.path.join(ROOT, "analysis", "loop_cv", "fold_all.json")))
        fold["pramp"] = [float(np.median([json.load(open(p))["pramp"][0] for p in glob.glob(os.path.join(
            ROOT, "analysis", "loop_cv", f"cap{os.environ['EMT_PMAX']}_folds", "fold_*.json"))])), 0.62]
        return pc.setting_G(lab, ap["steps"][0], law, fold, None)
    rt = min(LIB_RT, key=lambda r: abs(r - lab["rt"]))
    fold = json.load(open(os.path.join(CV, f"fold_{rt}.json")))
    cells = json.load(open(os.path.join(CV, "cells.json")))
    fc = {}
    key = f"{lab['lf']}_{lab['rt']}_1.0"
    if abs(rt - lab["rt"]) < 1e-9 and key in cells:
        fc.update({k: cells[key][k] for k in ("lf_gain", "lf_dc") if k in cells[key]})
    if lab["hf"] < 1.0:
        pm = pole_curve()
        fc["hf_poles"] = [pm(p) for p in law_poles(lab, law)]
    return pc.setting_G(lab, ap["steps"][0], law, fold, fc or None)


def render(G, v, n):
    return q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=80)])]).render(v, G, n, 80)


# ---------------------------------------------------------------- the units
def unit_irs(which="all"):
    import altiverb_predict as av
    out = []
    if which in ("all", "CLA"):
        out.append(("CLA", 1.0, av.cla_h(), 98.4))
    if which in ("all", "Altiverb"):
        out += [(f"Altiverb RT {k}", rt, av.altiverb_h(k)[0], 42.8) for k, rt in av.PRESETS.items()]
    return [(nm, rt, sosfiltfilt(q1.LFBAND, h * (np.sign(h[PRE].sum()) or 1.0), axis=0), e0) for nm, rt, h, e0 in out]


def noise_power():
    """Mean in-band (100 Hz-10 kHz) and per-octave noise power per channel from tails where the reverb has decayed
    below the noise (CLA RT 1.0: 1.6-2.0 s after tap 0, the reverb 96 dB down; the Altiverb RT 0.6 / 1.0 presets'
    1.5-1.9 s), in the units of the band-limited IRs."""
    import altiverb_predict as av
    import model_s as m
    h = np.load(q1.CLA_DEV)
    h = sosfiltfilt(q1.LFBAND, h * (np.sign(h[m.K0].sum()) or 1.0), axis=0)
    tails = {"CLA": [h[m.K0 + int(1.6 * FS): len(h) - 200]]}
    tails["Altiverb"] = []
    for k in ("0.6", "1.0"):
        a = av.altiverb_h(k)[0]
        a = sosfiltfilt(q1.LFBAND, a * (np.sign(a[PRE].sum()) or 1.0), axis=0)
        tails["Altiverb"].append(a[PRE + int(1.5 * FS): len(a) - 200])
    out = {}
    for u, segs in tails.items():
        d = {"bb": np.mean([np.mean(sosfiltfilt(BAND, s, axis=0)[500:-500] ** 2, 0) for s in segs], 0).tolist()}
        for f, sos in OCT_SOS.items():
            d[str(f)] = np.mean([np.mean(sosfiltfilt(sos, s, axis=0)[500:-500] ** 2, 0) for s in segs], 0).tolist()
        out[u] = d
    return out


def floors(y, npow):
    """Noise-limited null per window (100 Hz-10 kHz) and per octave and window, per channel: noise energy over the window
    / the unit's in-band energy there (y already aligned and scaled to the model)."""
    yb = sosfiltfilt(BAND, y, axis=0)
    out = {"bb": {}, "oct": {}}
    for a, b in WIN:
        i0, i1 = PRE + int(a / 1000 * FS), PRE + int(b / 1000 * FS)
        out["bb"][f"{a}-{b}"] = [float(10 * np.log10(npow["bb"][c] * (i1 - i0) / np.sum(yb[i0:i1, c] ** 2))) for c in (0, 1)]
    for f, sos in OCT_SOS.items():
        ybf = sosfiltfilt(sos, y, axis=0)
        out["oct"][str(f)] = [[float(10 * np.log10(npow[str(f)][c] * (PRE + int(b / 1000 * FS) - PRE - int(a / 1000 * FS))
                                                   / np.sum(ybf[PRE + int(a / 1000 * FS): PRE + int(b / 1000 * FS), c] ** 2))) for c in (0, 1)] for a, b in WIN]
    return out


def score(y, m, te):
    """Broadband (100 Hz-10 kHz) and per-octave nulls per window and channel; None past te or under 25 ms."""
    out = {"bb": {}, "oct": {}}
    yb, mb = sosfiltfilt(BAND, y, axis=0), sosfiltfilt(BAND, m, axis=0)

    def nl(a_, b_, ya, ma):
        i0, i1 = PRE + int(a_ / 1000 * FS), min(PRE + int(b_ / 1000 * FS), PRE + int(te * FS), len(ya))
        if i1 - i0 < int(0.025 * FS):
            return None
        return [float(10 * np.log10(np.sum((ya[i0:i1, c] - ma[i0:i1, c]) ** 2) / np.sum(ya[i0:i1, c] ** 2))) for c in (0, 1)]
    for a, b in WIN:
        out["bb"][f"{a}-{b}"] = nl(a, b, yb, mb)
    for f, sos in OCT_SOS.items():
        yf, mf = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, m, axis=0)
        out["oct"][str(f)] = [nl(a, b, yf, mf) for a, b in WIN]
    return out


# ---------------------------------------------------------------- the compact chain
def rows(wins, n):
    """Sample indices of the fit windows (0 ms starts 40 samples before tap 0, as the cluster-only FIR)."""
    idx = [np.arange(PRE + int(a / 1000 * FS) - (40 if a == 0 else 0), min(PRE + int(b / 1000 * FS), n)) for a, b in wins]
    return np.concatenate(idx)


X3 = [420.0, -0.2, 2.0]                                         # the third peak (500 Hz: §27's 10 numbers leave -35 dB there)
LO3, HI3 = ut.LO + [350.0, -3.0, 0.5], ut.HI + [700.0, 3.0, 8.0]


def lb_resp(x, f):
    """§27's cascade (10 numbers: high-pass, low shelf, two peaks) and, with 13 numbers, a third peaking section."""
    H = ut.bq_resp(x[:10], f)
    if len(x) > 10:
        s_ = ut.peak_sos(*x[10:13])
        z = np.exp(-2j * np.pi * np.asarray(f) / FS)
        H = H * (s_[0] + s_[1] * z + s_[2] * z ** 2) / (s_[3] + s_[4] * z + s_[5] * z ** 2)
    return H


def lowband(m, x):
    return ut.through(m, lambda f: np.repeat(lb_resp(x, f)[:, None], 2, 1)) if x is not None else m


def fir_fit(mb, yb, idx, K=KF):
    """Per-channel (2K+1)-tap FIR on rows idx (least squares; band-limited inputs). Returns taps (2, 2K+1)."""
    hs = []
    for c in (0, 1):
        A = np.stack([np.roll(mb[:, c], k)[idx] for k in range(-K, K + 1)], 1)
        hs.append(np.linalg.lstsq(A, yb[idx, c], rcond=None)[0])
    return np.array(hs)


def fir_apply(m, hs):
    K = (hs.shape[1] - 1) // 2
    return np.stack([sum(hs[c, k + K] * np.roll(m[:, c], k) for k in range(-K, K + 1)) for c in (0, 1)], 1)


def chain_fit(yb_raw, m0, wins, eps0, x0, iters=3, fit_bq=True):
    """The compact chain on windows wins: clock and shift of the data (2), per-channel 7-tap FIR (14) and the low-band
    cascade (10), alternating (clock / shift | cascade with the FIR projected out). Returns dict(eps, shift, bq, fir)."""
    n = len(m0)
    idx = rows(wins, n)
    i1 = int(idx.max()) + 1
    m0b = sosfiltfilt(BAND, m0, axis=0)
    ybb = sosfiltfilt(BAND, yb_raw, axis=0)
    eps, sh, x = eps0, 0.0, (list(x0) if x0 is not None else None)
    for it in range(iters):
        mb = lowband(m0b, x)
        hs = fir_fit(mb, q1.apply_clock(ybb, eps, sh)[:n], idx)
        mc = fir_apply(mb, hs)

        def cost(p):
            y = q1.warp(q1.fshift(ybb[:i1 + 400], p[1]), p[0], half=32, beta=8.0)[:i1]
            return np.sum((y[idx] - mc[idx]) ** 2) / np.sum(y[idx] ** 2)
        best = min(np.arange(eps - 6, eps + 6.1, 1.0), key=lambda e: cost([e, sh])) if it == 0 else eps
        r = minimize(cost, [best, sh], method="Nelder-Mead",
                     options=dict(initial_simplex=[[best, sh], [best + 1, sh], [best, sh + 0.02]], xatol=0.005, fatol=1e-10))
        eps, sh = float(r.x[0]), float(r.x[1])
        y = q1.apply_clock(ybb, eps, sh)[:n]
        if fit_bq and x is not None:
            def res(v):
                mb_ = lowband(m0b, v)
                h_ = fir_fit(mb_, y, idx)
                return (y[idx] - fir_apply(mb_, h_)[idx]).ravel() / np.sqrt(np.mean(y[idx] ** 2))
            lo, hi = (ut.LO, ut.HI) if len(x) == 10 else (LO3, HI3)
            s_ = least_squares(res, np.clip(x, lo, hi), bounds=(lo, hi), x_scale="jac", max_nfev=400)
            x = s_.x.tolist()
    mb = lowband(m0b, x)
    hs = fir_fit(mb, q1.apply_clock(ybb, eps, sh)[:n], idx)
    return dict(eps=eps, shift=sh, bq=x, fir=hs.tolist())


def chain_model(m0, ch):
    return fir_apply(lowband(m0, ch["bq"]), np.array(ch["fir"]))


def np_transfer(y, m, idx):
    """Free complex transfer per channel and 1/12 octave measured on rows idx (Tukey-tapered), interpolated in log f,
    1 / 0 deg held outside 40 Hz-11.3 kHz. Returns H(f) callable (len, 2)."""
    seg = np.zeros(len(y))
    seg[idx] = 1.0
    w = np.convolve(seg, np.hanning(64) / np.hanning(64).sum(), "same")
    N = 1 << int(np.ceil(np.log2(len(y) + 4096)))
    f = np.fft.rfftfreq(N, 1 / FS)
    Y, M = (np.fft.rfft(x * w[:, None], N, axis=0) for x in (y, m))
    T = np.ones((len(TWELFTH), 2), complex)
    for c in (0, 1):
        for i, fc in enumerate(TWELFTH):
            s = (f >= fc * 2 ** (-1 / 24)) & (f < fc * 2 ** (1 / 24))
            if s.sum():
                T[i, c] = np.sum(Y[s, c] * np.conj(M[s, c])) / np.sum(np.abs(M[s, c]) ** 2)
    lf = np.log(TWELFTH)

    def H(ff):
        u = np.log(np.maximum(ff, 1.0))
        return np.stack([np.interp(u, lf, np.abs(T[:, c]), left=abs(T[0, c]), right=abs(T[-1, c]))
                         * np.exp(1j * np.interp(u, lf, np.unwrap(np.angle(T[:, c])), left=np.angle(T[0, c]), right=np.unwrap(np.angle(T[:, c]))[-1]))
                         for c in (0, 1)], 1)
    return H


def shift_gain_fit(ybb, mb, eps, idx):
    """Per-preset shift and one gain on rows idx (clock fixed). Returns (shift, gain)."""
    n = len(mb)

    def cost(s):
        y = q1.apply_clock(ybb, eps, s)[:n]
        g = np.sum(y[idx] * mb[idx]) / np.sum(mb[idx] ** 2)
        return np.sum((y[idx] - g * mb[idx]) ** 2) / np.sum(y[idx] ** 2), g
    from scipy.optimize import minimize_scalar
    s = minimize_scalar(lambda s: cost(s)[0], bounds=(-1.5, 1.5), method="bounded", options={"xatol": 1e-3}).x
    return float(s), float(cost(s)[1])


# ---------------------------------------------------------------- stage s
def s_stage():
    import altiverb_predict as av
    os.makedirs(OUT, exist_ok=True)
    v, law, ap = _law()
    bq = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))
    npow = noise_power()
    res = {"noise_power": npow, "irs": {}}
    shared = None
    irs = [r for r in unit_irs() if not os.environ.get("UC_ONLY") or r[0] in os.environ["UC_ONLY"].split(",")]   # smoke test
    order = sorted(irs, key=lambda r: (r[0] != "CLA", r[0] != "Altiverb RT 1.6"))      # CLA, then RT 1.6 (the shared chain)
    for name, rt, yb, e0 in order:
        unit = "CLA" if name == "CLA" else "Altiverb"
        lab = dict(lf=1.0, rt=rt, hf=0.9)
        R = res["irs"][name] = {"rt": rt}
        for mode in ("units", "base"):
            G = model_G(lab, mode, law, ap)
            n = min(len(yb) - 600, G["N"])
            m0 = render(G, v, n)
            if mode == "base":
                R["hf_poles"] = law_poles(lab, law), [pole_curve()(p) for p in law_poles(lab, law)]
            # fir / fir+bq: the §26.2 / §27 procedure
            for vn, x in (("fir", None), ("fir+bq", bq["Altiverb RT 1.6" if unit == "CLA" else "CLA"]["x"])):
                mdl = lowband(m0, x)
                eg, shg, gain = av.fit_clock(yb, mdl, starts=[e0])
                y = q1.apply_clock(yb, eg, shg)[:n] / gain
                mdl2, hs = av.channel_fir(y, mdl)
                eg2, shg2, gain2 = av.fit_clock(yb, mdl2, starts=[eg])
                y2 = q1.apply_clock(yb, eg2, shg2)[:n] / gain2
                te = min(1.5, q2.comp_end(y2, np.ones_like(y2)))
                R[f"{mode}/{vn}"] = dict(n_par=17, eps=eg2, t_end=te, s=score(y2, mdl2, te))
                if mode == "base" and vn == "fir":
                    R["floor"] = floors(y2, npow[unit])
                    R["te"] = te
                    gref = gain2
                print(f"{name} {mode}/{vn}: " + " ".join(f"{w}:{x_[1]:.1f}/{x_[0]:.1f}" for w, x_ in R[f"{mode}/{vn}"]["s"]["bb"].items() if x_), flush=True)
            if mode != "base":
                continue
            te = R["te"]
            xo = bq["Altiverb RT 1.6" if unit == "CLA" else "CLA"]["x"]
            for wn, wins in FIT_W.items():
                for tag, x0 in (("chain10", xo), ("chain", xo + X3)):
                    ch = chain_fit(yb / gref, m0, wins, e0, x0)
                    y = q1.apply_clock(yb / gref, ch["eps"], ch["shift"])[:n]
                    R[f"{tag}@{wn}"] = dict(n_par=16 + len(x0), chain=ch, s=score(y, chain_model(m0, ch), te))
                Hn = np_transfer(sosfiltfilt(BAND, y, axis=0), sosfiltfilt(BAND, m0, axis=0), rows(wins, n))
                R[f"np@{wn}"] = dict(n_par=2 + 4 * len(TWELFTH), s=score(y, ut.through(m0, Hn), te))
                for vn in (f"chain10@{wn}", f"chain@{wn}", f"np@{wn}"):
                    print(f"{name} {vn}: " + " ".join(f"{w}:{x_[1]:.1f}/{x_[0]:.1f}" for w, x_ in R[vn]["s"]["bb"].items() if x_), flush=True)
            if name == "Altiverb RT 1.6":
                shared = R["chain@all"]["chain"]
            if unit == "Altiverb" and shared is not None and name != "Altiverb RT 1.6":
                mc = sosfiltfilt(BAND, chain_model(m0, shared), axis=0)
                ybb = sosfiltfilt(BAND, yb / gref, axis=0)
                s_, g_ = shift_gain_fit(ybb, mc, shared["eps"], rows([(0, 50)], n))
                y = q1.apply_clock(yb / gref, shared["eps"], s_)[:n] / g_
                R["rt1.6"] = dict(n_par=2, shift=s_, gain=g_, s=score(y, chain_model(m0, shared), te))
                print(f"{name} rt1.6: " + " ".join(f"{w}:{x_[1]:.1f}/{x_[0]:.1f}" for w, x_ in R["rt1.6"]["s"]["bb"].items() if x_), flush=True)
        json.dump(res, open(os.path.join(OUT, "s.json"), "w"))


# ---------------------------------------------------------------- stage ml
def ml_stage():
    import altiverb_ml as am
    import ml_carryover as cc
    import ml_lines as mll
    S = json.load(open(os.path.join(OUT, "s.json")))
    bq = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))
    shared = S["irs"]["Altiverb RT 1.6"]["chain@all"]["chain"]
    K0 = am.K0
    variants = {"none": None, "bq_CLA": bq["CLA"]["x"], "bq_altS": shared["bq"]}
    res = []
    for prog, (pre, lf, hf) in am.SET.items():
        names = sorted((os.path.basename(p)[:-4] for p in glob.glob(os.path.join(am.CACHE, pre + " *.npy"))), key=lambda t: float(t.split(" ")[-1][:-1]))
        for name in names:
            rt = float(name.split(" ")[-1][:-1])
            n = K0 + int(min(1.5, 1.2 * rt + 0.5) * FS)
            y = am.device_ir(name, n)
            i1 = min(am.usable(y), n - 1)
            key = f"{lf}_{rt}_{hf}"
            M = mll.Model(prog)
            row = dict(prog=prog, preset=name, rt=rt)
            # b and the clock as altiverb_ml (chosen once, on the plain model), then each variant through the same fit
            fb = None
            for b in am.B_STEPS:
                M.b[key] = b
                C, Rr = M.parts(n, key)
                fpn = pc._ml_nulls(prog, y, am.fit(prog, y, C, Rr, i1), i1)["first_pass"]
                if fb is None or fpn < fb[1]:
                    fb = (b, fpn, C, Rr)
            b, _, C0, R0 = fb
            best = None
            for d in np.arange(-6.0, 6.01, 1.0):
                C, Rr = (mll.warp(C0, 0.0, d), mll.warp(R0, 0.0, d)) if d else (C0, R0)
                nl = pc._ml_nulls(prog, y, am.fit(prog, y, C, Rr, i1), i1)
                if best is None or (nl["recirc1"] or 0) < (best[1]["recirc1"] or 0):
                    best = (d, nl, C, Rr)
            d, _, C, Rr = best
            row.update(b=b, ppm=am.EPS + d)
            W = mll.windows(prog)[1]
            for vn, x in variants.items():
                Cv, Rv = lowband(C, x), lowband(Rr, x)
                mdl = am.fit(prog, y, Cv, Rv, i1)
                row[vn] = dict(n_par=22 + 4 + 1, null=dict(cluster_gap=cluster_gap(prog, y, mdl), **pc._ml_nulls(prog, y, mdl, i1)), fp_oct=fp_oct(y, mdl, K0, W))
            # the whole S chain (FIR + cascade, Altiverb RT 1.6): a shift per preset and the two gains per channel
            hs = np.array(shared["fir"])
            Cs, Rs = fir_apply(lowband(C, shared["bq"]), hs), fir_apply(lowband(Rr, shared["bq"]), hs)
            sl = slice(K0 - 30, K0 + cc.ONSET[prog][1] - 3)

            def cost(s):
                Cx = q1.fshift(Cs, s)
                return -np.dot(y[sl, 1], Cx[sl, 1]) ** 2 / np.dot(Cx[sl, 1], Cx[sl, 1])
            from scipy.optimize import minimize_scalar
            s_ = float(minimize_scalar(cost, bounds=(-1.5, 1.5), method="bounded", options={"xatol": 1e-3}).x)
            Cx, Rx = q1.fshift(Cs, s_), q1.fshift(Rs, s_)
            mdl = np.zeros_like(y)
            for c in (0, 1):
                s2 = slice(K0 - 30, i1)
                g = np.linalg.lstsq(np.stack([Cx[s2, c], Rx[s2, c]], 1), y[s2, c], rcond=None)[0]
                mdl[:, c] = g[0] * Cx[:, c] + g[1] * Rx[:, c]
            row["chain_altS"] = dict(n_par=1 + 4, null=dict(cluster_gap=cluster_gap(prog, y, mdl), **pc._ml_nulls(prog, y, mdl, i1)), fp_oct=fp_oct(y, mdl, K0, W))
            res.append(row)
            print(f"{prog} {name}: b {b:.2f} | " + " | ".join(f"{vn} fp {row[vn]['null']['first_pass']:.1f} r1 {row[vn]['null']['recirc1']:.1f}"
                                                         for vn in list(variants) + ["chain_altS"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "ml.json"), "w"))


def cluster_gap(prog, y, m):
    """Cluster-and-gap null (both channels' means, one gain each), as altiverb_ml.run_one."""
    import altiverb_ml as am
    import ml_carryover as cc
    v = []
    for ch in (0, 1):
        s = slice(am.K0 - 30, am.K0 + cc.ONSET[prog][ch] - 3)
        ya, ma = y[s, ch], m[s, ch]
        v.append(10 * np.log10(np.sum((ya - np.dot(ya, ma) / np.dot(ma, ma) * ma) ** 2) / np.sum(ya ** 2)))
    return float(np.mean(v))


def fp_oct(y, m, K0, W):
    """R first-pass null per octave (one gain per band)."""
    a0, a1 = K0 + W["first_pass"][0], K0 + W["first_pass"][1]
    out = {}
    for f, sos in OCT_SOS.items():
        ya, ma = sosfiltfilt(sos, y[:, 1])[a0:a1], sosfiltfilt(sos, m[:, 1])[a0:a1]
        out[str(f)] = float(10 * np.log10(np.sum((ya - np.dot(ya, ma) / np.dot(ma, ma) * ma) ** 2) / np.sum(ya ** 2)))
    return out


# ---------------------------------------------------------------- stage diag (CLA)
def diag_stage():
    S = json.load(open(os.path.join(OUT, "s.json")))
    v, law, ap = _law()
    bq = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))
    (name, _, yb, e0), = unit_irs("CLA")
    R = S["irs"]["CLA"]
    te = R["te"]
    out = {"scan": [], "scan_alt": []}
    (_, _, ya16, ea), = [r for r in unit_irs("Altiverb") if r[0] == "Altiverb RT 1.6"]
    for hf in (0.6, 0.7, 0.8, 0.88, 0.89, 0.895, 0.9, 0.905, 0.91, 0.92, 1.0):   # the Altiverb's inferred HF 0.9 (sensitivity, not a fit)
        lab = dict(lf=1.0, rt=1.6, hf=hf)
        G = model_G(lab, "base", law, ap)
        n = min(len(ya16) - 600, G["N"])
        m0 = render(G, v, n)
        ch = chain_fit(ya16, m0, FIT_W["0-50"], ea, bq["CLA"]["x"] + X3)
        y = q1.apply_clock(ya16, ch["eps"], ch["shift"])[:n]
        sc = score(y, chain_model(m0, ch), S["irs"]["Altiverb RT 1.6"]["te"])
        out["scan_alt"].append(dict(hf=hf, rt=1.6, s=sc))
        print(f"Altiverb RT 1.6 HF {hf}: " + " ".join(f"{w}:{x_[1]:.1f}/{x_[0]:.1f}" for w, x_ in sc["bb"].items() if x_)
              + " | 8k R " + "/".join("--" if x_ is None else f"{x_[1]:.1f}" for x_ in sc["oct"]["8000"]), flush=True)
    for hf, rt in [(h, 1.0) for h in (0.8, 0.85, 0.88, 0.89, 0.895, 0.9, 0.905, 0.91, 0.92, 0.95, 1.0)] + [(0.9, r) for r in (0.97, 0.985, 1.015, 1.03)]:
        lab = dict(lf=1.0, rt=rt, hf=hf)
        G = model_G(lab, "base", law, ap)
        n = min(len(yb) - 600, G["N"])
        m0 = render(G, v, n)
        ch = chain_fit(yb, m0, FIT_W["0-50"], e0, bq["Altiverb RT 1.6"]["x"] + X3)
        y = q1.apply_clock(yb, ch["eps"], ch["shift"])[:n]
        sc = score(y, chain_model(m0, ch), te)
        out["scan"].append(dict(hf=hf, rt=rt, s=sc))
        print(f"HF {hf} RT {rt}: " + " ".join(f"{w}:{x_[1]:.1f}/{x_[0]:.1f}" for w, x_ in sc["bb"].items() if x_), flush=True)
    # residual correlation: CLA (chain@0-50) against the model's residual on Tomkin 1.0_1.0_1.0
    lab = dict(lf=1.0, rt=1.0, hf=0.9)
    G = model_G(lab, "base", law, ap)
    n = min(len(yb) - 600, G["N"])
    m0 = render(G, v, n)
    ch = R["chain@0-50"]["chain"]
    yc = sosfiltfilt(BAND, q1.apply_clock(yb, ch["eps"], ch["shift"])[:n], axis=0)
    rc = yc - sosfiltfilt(BAND, chain_model(m0, ch), axis=0)
    info, fade, _, _, eps0 = q2.q1_tables()
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    k = "1.0_1.0_1.0"
    g_ = info[k]
    cells = json.load(open(os.path.join(CV, "cells.json")))
    Gt = pc.setting_G(dict(lf=1.0, rt=1.0, hf=1.0), ap["steps"][ap["settings"][k]["step"]], law, json.load(open(os.path.join(CV, "fold_1.0.json"))), cells.get(k))
    out["corr"] = []
    for amp, nm in q2.files_of(g_):
        yt, mt, gt, tet = q2.fit_capture(Gt, v, amp, nm, fade[nm], eps0[nm], ends[nm], 1700)
        rt_ = sosfiltfilt(BAND, yt - mt * gt, axis=0)
        row = dict(amp=amp, file=nm)
        for a, b in WIN[1:]:
            i0, i1 = PRE + int(a / 1000 * FS), PRE + int(b / 1000 * FS)
            i1 = min(i1, len(rt_), len(rc), PRE + int(tet * FS))
            row[f"{a}-{b}"] = [float(np.dot(rc[i0:i1, c], rt_[i0:i1, c]) / np.sqrt(np.dot(rc[i0:i1, c], rc[i0:i1, c]) * np.dot(rt_[i0:i1, c], rt_[i0:i1, c]))) for c in (0, 1)]
            row[f"null_{a}-{b}"] = [float(10 * np.log10(np.sum(rt_[i0:i1, c] ** 2) / np.sum(sosfiltfilt(BAND, yt, axis=0)[i0:i1, c] ** 2))) for c in (0, 1)]
        out["corr"].append(row)
        print(f"Tomkin {k} AMPL {amp}: corr with CLA residual " + " ".join(f"{w}:{row[w][1]:+.2f}/{row[w][0]:+.2f}" for w in ("50-150", "150-300", "300-600")), flush=True)
    # captures against each other (the repeat's residual correlation with CLA's, a reference for chance)
    json.dump(out, open(os.path.join(OUT, "diag.json"), "w"))


# ---------------------------------------------------------------- report
def report_stage():
    S = json.load(open(os.path.join(OUT, "s.json")))
    f2 = lambda x: "--" if x is None else f"{x[1]:.1f} ({x[0]:.1f})"
    out = ["# Per-unit chain model (tools/unit_chain.py); nulls 100 Hz-10 kHz, R (L), dB; [n] per-unit numbers", ""]
    for name, R in S["irs"].items():
        fl = R["floor"]["bb"]
        out.append(f"## {name} (t_end {R['te']:.2f} s); noise-limited null " + " ".join(f"{w}: {f2(x)}" for w, x in fl.items()))
        if "hf_poles" in R:
            out.append(f"   HF 0.9 law poles {np.round(R['hf_poles'][0], 4).tolist()} -> curve {np.round(R['hf_poles'][1], 4).tolist()}")
        for vn in [k for k in R if "/" in k or k.startswith(("chain", "np", "rt1.6"))]:
            r = R[vn]
            line = f"- {vn:16s} [{r['n_par']}] " + " | ".join(f"{w}: {f2(x)}" for w, x in r["s"]["bb"].items())
            ex = [f"{x[1] - fl[w][1]:+.1f}" if x else "--" for w, x in r["s"]["bb"].items()]
            out.append(line + "   excess R " + " / ".join(ex))
            if vn.startswith(("chain", "np")) or vn.endswith("fir+bq") or vn == "base/fir":
                out.append("    oct R " + "; ".join(f"{f}: " + "/".join("--" if x is None else f"{x[1]:.1f}" for x in r["s"]["oct"][str(f)]) for f in OCT))
        out.append("    floor oct R " + "; ".join(f"{f}: " + "/".join(f"{x[1]:.0f}" for x in R["floor"]["oct"][str(f)]) for f in OCT))
        out.append("")
    if os.path.exists(os.path.join(OUT, "ml.json")):
        M = json.load(open(os.path.join(OUT, "ml.json")))
        out.append("## M / L first pass, R (recirc 1), Altiverb MED / LONG")
        for r in M:
            out.append(f"- {r['prog']} {r['preset']} (b {r['b']:.2f}): " + " | ".join(
                f"{vn} {r[vn]['null']['first_pass']:.1f} ({r[vn]['null']['recirc1']:.1f})" for vn in ("none", "bq_CLA", "bq_altS", "chain_altS")))
        for prog in ("M", "L"):
            rr = [r for r in M if r["prog"] == prog and r["b"] == 0.0]
            out.append(f"- {prog} median, b = 0 presets ({len(rr)}): " + " | ".join(
                f"{vn} fp {np.median([r[vn]['null']['first_pass'] for r in rr]):.1f} cl {np.median([r[vn]['null']['cluster_gap'] for r in rr]) if 'cluster_gap' in rr[0][vn]['null'] else float('nan'):.1f} "
                f"r1 {np.median([r[vn]['null']['recirc1'] for r in rr]):.1f}" for vn in ("none", "bq_CLA", "bq_altS", "chain_altS")))
            for vn in ("none", "bq_CLA", "chain_altS"):
                out.append(f"    {vn} fp oct R " + " ".join(f"{f}: {np.median([r[vn]['fp_oct'][str(f)] for r in rr]):.1f}" for f in OCT))
    if os.path.exists(os.path.join(OUT, "diag.json")):
        D = json.load(open(os.path.join(OUT, "diag.json")))
        out.append("\n## CLA diag: chain@0-50 null vs HF / RT (sensitivity, not a fit)")
        for r in D["scan"]:
            out.append(f"- HF {r['hf']} RT {r['rt']}: " + " | ".join(f"{w}: {f2(x)}" for w, x in r["s"]["bb"].items()))
        for r in D.get("scan_alt", []):
            out.append(f"- Altiverb RT 1.6 HF {r['hf']}: " + " | ".join(f"{w}: {f2(x)}" for w, x in r["s"]["bb"].items())
                       + " | oct R 50-150 " + " ".join(f"{f}: {r['s']['oct'][str(f)][1][1]:.1f}" for f in OCT))
        out.append("- residual correlation CLA (chain@0-50) vs model-capture on Tomkin 1.0_1.0_1.0, R (L):")
        for r in D["corr"]:
            out.append(f"  - AMPL {r['amp']}: " + " ".join(f"{w}: {r[w][1]:+.2f} ({r[w][0]:+.2f}) [Tomkin null {r['null_' + w][1]:.1f}]" for w in ("50-150", "150-300", "300-600")))
    print("\n".join(out))
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["s", "ml", "diag", "report"], required=True)
    a = ap.parse_args()
    {"s": s_stage, "ml": ml_stage, "diag": diag_stage, "report": lambda: None}[a.stage]()
    report_stage()


if __name__ == "__main__":
    main()
