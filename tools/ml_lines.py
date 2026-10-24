"""NEXT_STEPS §5 item 5, after the carry-over: the line gains and filters of programs M and L, per setting, through the
coupled model of docs/ml_structure.md §7 (Tomkin 96 kHz library, device-rate IRs of tools/ml_structure.py).
Results: docs/ml_structure.md §8.

    .venv/bin/python tools/ml_lines.py [--stage law|report|fit|all] [--prog ML] [--workers 4] [--only N] [--chunk i/n]

The model is the carry-over's (ml_carryover.py: cluster, diffuser with the per-setting damping b, 20 equal sign-only
taps per channel, output kernel; positions, signs and kernel trained on 6 settings, firstpass.json) closed through
the four lines:

    y = h * [ g_c C + g_t D 1^T (I - Lambda F)^-1 T ],   Lambda = diag(G_j(z) z^-P_j)

T: per S loop i the sum of its read taps; F: the coupling (M and L: the orthogonal bipartite matrix of §7.3; S: I);
G_j: line j's per-pass response from S's control laws (model_s_results §14.2) at line j's own length (Jot gain,
LF shelf with the fixed pole, HF one-pole), with S's coefficient limits (HF pole <= 0.84, per-pass DC gain >= 0.18,
beta = 0.019; §16.15-16.16). LF steps that S has no law for (M 0.7 / 1.2, L 0.8 / 1.2) are interpolated in log LF.
No loop number is fitted on M or L. Rendered in the frequency domain on an exponentially weighted grid (z -> r z, no
time aliasing at long RT) with integer delays; the capture clock and shift are applied afterwards (Kaiser-sinc warp).

Per capture: one sub-sample shift, one clock ratio and two gains per channel (cluster, reverb), fitted on the
deterministic part [0, onset + shortest line) only, which no line reaches. Everything after it is predicted:

  law     per AMPL group (one checkpoint file each in analysis/ml_lines/law/; finished groups are skipped):
          nulls (100 Hz-10 kHz) per window, cluster and gap / first pass / recirculations 1-3 (each one shortest line
          long), up to the capture's fade, for the law, the lines uncoupled (F = I, control), the law without the
          limits, and scans of each limit where it can act; the repeat-capture null with the same freedoms (another
          capture of the group as the model). Per octave band, decay and level against the model: 10 log10 of the
          band energy ratio measured / model per frame, from the onset to the fade; its slope (dB/s) is the decay
          error (energy-based, so the M/L time variance does not bias it) and gives RT(f) = 60 / (60 / RT_law(f) - slope).
  report  breakdowns of the law files by control, the capture-to-capture / model error split, the limit scans.
  ML_HK=<chain_kernel json> in the environment replaces the stored output kernels (firstpass.json) by the ones
  refitted with the chain kernel (tools/chain_kernel.py), and S's by a unit tap.
  ML_FP=<firstpass json> reads the first-pass parameters (positions, signs, kernel, b per setting) from that file
  instead of analysis/ml_carryover/firstpass.json, e.g. firstpass_chain2.json refitted on the chain-kernel IRs.
  Tools that honour either suffix their outputs with SUFFIX (_fp, _hk).
  fitted(prog, key)  (a function, no stage) rebuilds one group's captures with their models from the law files.
  fit     per setting, the line law's three numbers (Nyquist-end RT, DC RT ratio, HF RT factor at f_ref; shelf pole
          and f_ref fixed, length-scaled over the four lines) fitted to the band decay of all captures of the group:
          against S's laws per LF and HF step.
"""
import argparse, json, os, sys
for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                          # one thread per worker: the pool supplies the parallelism
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.optimize import minimize
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import ml_carryover as cc  # noqa: E402
import ml_structure as ml  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_lines")
FS, K0 = cc.FS, cc.K0
LAW = json.load(open(os.path.join(ROOT, "analysis", "q1", "laws.json")))
P_MAX, DC_MIN, BETA = 0.84, 0.18, 0.019                    # S's coefficient limits (model_s_results §16.15-16.16)
DROP_BACK = 0.8465                                          # S's HF 0.2 / RT 0.6 pole past the drop-back point (model_s_results 32.4)


def _pole_curve():
    """The HF pole as one curve of the law's pole: S's curve of model_s_results 32.6 (all-data fit, piecewise linear in the
    law's pole) with its drop-back past the fold point to S's HF 0.2 / RT 0.6 cell value, the default for M / L's lines
    since ml_structure 8.7. ML_POLE_CURVE=25: §25.17's curve (law, ramp above ~0.62, ceiling ~0.8475; the default from
    §25.20); ML_POLE_CURVE=0: the capped law (min(p, P_MAX))."""
    mode = os.environ.get("ML_POLE_CURVE", "32")
    if mode == "0":
        return None
    if mode == "25":
        pc = json.load(open(os.path.join(ROOT, "analysis", "lowband", "chain3_polecurve.json")))
        kn, th = np.array(pc["knots"]), np.array(pc["all"])
        return lambda p: min(p + float(np.interp(p, kn, th[:len(kn)])), float(th[len(kn)]))
    pc = json.load(open(os.path.join(ROOT, "analysis", "lowband", "chain3_pole_cells_allf.json")))
    kn, th, fold = np.array(pc["knots"]), np.array(pc["all"]), float(pc["fold"])
    return lambda p: float(np.interp(p, kn, th)) if p < fold else DROP_BACK


POLE_MAP = _pole_curve()
HFRT_MIN = 0.0 if os.environ.get("ML_HFRT") == "0" else 0.12   # the HF control's floor, RT x HF >= 0.12 s in S's time,
                                                            # x 2 / x 4 in M / L (model_s_results 29); ML_HFRT=0: off
L22_HF = 0.0 if os.environ.get("ML_L22") == "0" else 0.40     # L at RT 2.2: HF acts as max(HF, 0.40), not 0.48 / RT = 0.218
                                                            # (poles 0.55-0.61 instead of the 0.8475 ceiling); ML_L22=0: off
F_ML = np.array([[0, 1, 0, 1], [1, 0, -1, 0], [0, 1, 0, -1], [-1, 0, -1, 0]]) / np.sqrt(2)   # rows: lines; columns: loops
OCT = [125, 250, 500, 1000, 2000, 4000, 8000]
OCT_SOS = {f: butter(3, [f / np.sqrt(2), min(f * np.sqrt(2), 10500)], btype="band", fs=FS, output="sos") for f in OCT}
FRAME = 1024                                                # band-energy frame (39 ms)
FADE = 0.70 + 0.02                                          # Tomkin's fade ramp plus a guard (s before the last sample)
HK_FILE = os.environ.get("ML_HK", "")                       # output kernels from a chain_kernel.py write (S: none)
FP_FILE = os.environ.get("ML_FP", "")                       # first-pass parameters other than firstpass.json
SUFFIX = ("_fp" if FP_FILE else "") + ("_hk" if HK_FILE else "")


def fp_path():
    """The first-pass parameter file (ML_FP, else ml_carryover's firstpass.json)."""
    if not FP_FILE:
        return os.path.join(cc.OUT, "firstpass.json")
    return FP_FILE if os.path.isabs(FP_FILE) else os.path.join(ROOT, FP_FILE)
ALIAS_DB = 160.0                                            # exponential weighting: wrap-around suppressed this far


# ---------------------------------------------------------------- the line law
def lf_law(lf):
    """S's LF law numbers (hi, dc, k) at an LF value, interpolated in log LF between S's steps where needed."""
    if f"{lf}" in LAW["lf"]:
        return LAW["lf"][f"{lf}"]
    ks = sorted(LAW["lf"], key=float)
    x = np.log([float(k) for k in ks])
    return {n: float(np.exp(np.interp(np.log(lf), x, np.log([LAW["lf"][k][n] for k in ks])))) for n in ("hi", "dc", "k")}


def onepole_for(T, f):
    c = np.cos(2 * np.pi * f / FS)
    T2 = min(T ** 2, 0.999999)
    return ((1 - T2 * c) - np.sqrt((1 - T2 * c) ** 2 - (1 - T2) ** 2)) / (1 - T2)


def line_coefs(P, rt, lf, hf, caps=True, v=None, p_max=P_MAX, dc_min=DC_MIN, beta=BETA):
    """Per line (g, shelf r0, HF pole p): S's laws at line length P, or v = (rt_hi, dc ratio, HF factor) in place of
    the law's (L['hi'] rt, L['dc'] / L['hi'], hf ** (1 / k)); the limits applied when caps (p_max / dc_min None: that
    limit off)."""
    L = lf_law(lf)
    if caps and HFRT_MIN and v is None:                     # the HF control saturates (model_s_results 29), in program time
        hf = max(hf, HFRT_MIN * round(min(P) / 1321) / rt)
        if L22_HF and round(min(P) / 1321) == 4 and abs(rt - 2.2) < 1e-6:     # L's shortest RT step: the floor is 0.40 (ml_structure 8.6)
            hf = max(hf, L22_HF)
    if v is None:
        rt_hi, rt_dc, k = L["hi"] * rt, L["dc"] * rt, L["k"]
    else:
        rt_hi, rt_dc, k = v[0], v[1] * v[0], 1.0
    out = []
    for Pi in P:
        kk = -3 * Pi / FS
        g = 10 ** (kk / rt_hi)
        r0 = 10 ** (kk * (1 / rt_dc - 1 / rt_hi))
        if caps and dc_min is not None and g * r0 < dc_min:
            lift = 20 * np.log10(dc_min / (g * r0))
            g *= 10 ** (-beta * lift / 20)
            r0 = dc_min / g
        if v is None:
            T = 10 ** (kk / rt * k * (1 / hf - 1)) if hf < 1.0 else 1.0
        else:
            T = 10 ** (kk / rt_hi * (1 / v[2] - 1)) if v[2] < 1.0 else 1.0
        p = onepole_for(T, LAW["f_ref"]) if T < 1.0 else 0.0
        if caps and p_max is not None and p > 0.0:
            p = POLE_MAP(p) if POLE_MAP and p_max == P_MAX else min(p, p_max)   # §25.20: the pole curve unless another limit is asked
        out.append((g, r0, p))
    return out


def line_resp(coefs, zi):
    """Per-pass responses G_j(zi), zi = z^-1 (any complex grid)."""
    a = LAW["a"]
    out = []
    for g, r0, p in coefs:
        b0, b1 = (r0 * (1 - a) + (1 + a)) / 2, (r0 * (1 - a) - (1 + a)) / 2
        out.append(g * (b0 + b1 * zi) / (1 - a * zi) * (1 - p) / (1 - p * zi))
    return np.stack(out, 1)


def rt_of(coefs, P, f):
    """RT (s) at frequency f from the per-pass response (mean decay rate over the lines)."""
    zi = np.exp(-2j * np.pi * np.atleast_1d(f) / FS)
    G = np.abs(line_resp(coefs, zi))
    rate = np.mean(20 * np.log10(np.maximum(G, 1e-12)) / np.array(P)[None] * FS, 1)    # dB/s
    return -60 / rate


# ---------------------------------------------------------------- the model
class Model:
    """The carry-over model of one program with its lines closed; parts rendered per capture clock."""

    def __init__(self, prog):
        fp = json.load(open(fp_path()))[prog]
        self.prog = prog
        self.P = cc.LINES[prog]
        self.pos = np.array(fp["positions"])
        self.sgn = np.array(fp["signs"])
        self.h = np.array(fp["kernel"])
        if HK_FILE:                                         # ML_HK: S a unit tap, M and L the kernels refitted there
            hk = json.load(open(HK_FILE if os.path.isabs(HK_FILE) else os.path.join(ROOT, HK_FILE)))["output_kernels"]
            self.h = np.array(hk[prog]) if prog in hk else np.array([1.0])
        self.kh = (len(self.h) - 1) // 2
        self.b = {}
        for r in fp["rows"]:
            self.b[r["key"]] = r["b"]
        Ls, Rs = cc.s_table()
        self.loop = [np.array([i for _, _, i in Ls]), np.array([i for _, _, i in Rs])]
        self.cl = cc.cluster_taps(prog)
        self.F = np.eye(4) if prog == "S" else F_ML

    def parts(self, n, key, coefs=None, F=None, band=True):
        """(cluster, reverb), (n, 2) each, unit gain, band-limited, on the nominal device clock with cluster tap 0 at K0.
        All delays are integers (the taps' fractional parts are the capture clock, ml_structure §7.2), so the weighted
        grid stays causal; clock and shift come after (warp)."""
        lf, rt, hf = (float(x) for x in key.split("_"))
        coefs = coefs if coefs is not None else line_coefs(self.P, rt, lf, hf)
        F = self.F if F is None else F
        N = 1 << int(np.ceil(np.log2(n + 4096)))
        r = 10 ** (ALIAS_DB / 20 / N)
        w = 2 * np.pi * np.fft.rfftfreq(N)
        zi = np.exp(-1j * w) / r
        zd = lambda d: np.exp(-1j * w * d) / r ** d            # z^-d (integer d)
        D = np.ones(len(w), complex)
        Ms = cc.DIFF[self.prog]
        for k, M in enumerate(Ms):
            bb = self.b.get(key, 0.0) if k == len(Ms) - 1 else 0.0
            Fz = zd(M) * ((1 - bb) + bb * zi)
            D *= (-cc.G_AP + Fz) / (1 - cc.G_AP * Fz)
        Hk = sum(hk * zd(K0 + k - self.kh) for k, hk in enumerate(self.h))
        G = line_resp(coefs, zi) * np.stack([zd(Pj) for Pj in self.P], 1)          # (nb, 4)
        A = np.eye(4)[None] - G[:, :, None] * F[None]
        C = np.zeros((len(w), 2), complex)
        T = np.zeros((len(w), 4, 2), complex)
        for ch in (0, 1):
            for p_, sg in self.cl[ch]:
                C[:, ch] += sg * zd(p_)
            for p_, sg, i in zip(np.round(self.pos[ch]), self.sgn[ch], self.loop[ch]):
                T[:, i, ch] += sg * zd(p_)
        R = np.linalg.solve(A, T).sum(1) * D[:, None]
        grow = r ** np.arange(n)[:, None]
        raw = [np.fft.irfft(X * Hk[:, None], N, axis=0)[:n] * grow for X in (C, R)]
        return [sosfiltfilt(ml.BAND, x, axis=0) for x in raw] if band else raw     # band=False: full band (outofband.py)


HALF = 48                                                   # warp: Kaiser-windowed sinc half-length
PHASES = 1024                                               # warp: table phases (linear between them, error < -110 dB)
BLOCK = 8192                                                # warp: rows per block (bounds the temporaries)


def _warp_table(beta=10.0):
    k = np.arange(-HALF + 1, HALF + 1)
    d = np.arange(PHASES + 1)[:, None] / PHASES - k[None]
    return np.sinc(d) * np.i0(beta * np.sqrt(np.clip(1 - (d / (HALF + 0.5)) ** 2, 0, None))) / np.i0(beta)


WTAB = _warp_table()


def warp(x, s, ppm):
    """x read at K0 + (n - K0) / (1 + ppm 1e-6) - s (the nominal-clock signal as a capture with clock ratio ppm and
    delay s sees it), per column; Kaiser-windowed sinc from a phase table, in blocks."""
    k = np.arange(-HALF + 1, HALF + 1)
    xp = np.vstack([np.zeros((HALF + 1, x.shape[1])), x, np.zeros((HALF + 2, x.shape[1]))])
    out = np.empty(x.shape)
    for b0 in range(0, len(x), BLOCK):
        n = np.arange(b0, min(b0 + BLOCK, len(x)))
        t = K0 + (n - K0) / (1 + ppm * 1e-6) - s
        i = np.floor(t).astype(int)
        ph = (t - i) * PHASES
        j = np.minimum(ph.astype(int), PHASES - 1)
        fr = (ph - j)[:, None]
        wgt = WTAB[j] * (1 - fr) + WTAB[j + 1] * fr
        idx = np.clip(i[:, None] + k[None] + HALF + 1, 0, len(xp) - 1)
        out[n] = np.einsum("nk,nkc->nc", wgt, xp[idx])
    return out


def windows(prog):
    """Per channel, lag windows (start, end): cluster and gap, first pass, recirculations 1-3 (one shortest line each)."""
    Pm = min(cc.LINES[prog])
    out = []
    for ch in (0, 1):
        o = cc.ONSET[prog][ch]
        out.append({"cluster_gap": (0, o - 3), "first_pass": (o - 3, o + Pm - 10), "recirc1": (o + Pm - 10, o + 2 * Pm - 10),
                    "recirc2": (o + 2 * Pm - 10, o + 3 * Pm - 10), "recirc3": (o + 3 * Pm - 10, o + 4 * Pm - 10)})
    return out


_FADES = None


def fade_start(f, prog):
    """Seconds before a file's end where Tomkin's fade starts, plus a guard: M / L from tools/ml_fade.py's per-file table
    (ml_structure §15: M a linear 1.0 s ramp, ~0.3 s on files under 1.5 s; L linear, ~0.3 x the file length), its
    program fallback where the reverb meets the floor before the fade; S (or no file) S's 0.70 s ramp (FADE)."""
    global _FADES
    if f is None or prog not in ("M", "L"):
        return FADE, None
    if _FADES is None:
        _FADES = json.load(open(os.path.join(ROOT, "analysis", "ml_fade", "fade.json")))
    e = _FADES.get(os.path.basename(f))
    if e is not None:
        return e["start"] + _FADES["_guard"], None
    return (None, _FADES["_fallback"]["L"]) if prog == "L" else (_FADES["_fallback"]["M"] + _FADES["_guard"], None)


def usable_end(y, f=None, prog=None):
    """Last index before Tomkin's fade (per channel): the file's end (the last sample above 1e-6 of the peak) minus the
    fade and a guard; f (the capture's file) and prog select the measured M / L fade, else S's 0.72 s (FADE)."""
    a = np.abs(y)
    last = [int(np.nonzero(a[:, ch] > a[:, ch].max() * 1e-6)[0][-1]) for ch in (0, 1)]
    t, frac = fade_start(f, prog)
    if frac is not None:                            # L fallback: a fraction of the file's length (from tap 0)
        return [int(e - frac * (e - K0) - _FADES["_guard"] * FS) for e in last]
    return [int(e - t * FS) for e in last]


def gains(y, parts, sl):
    """Per channel least-squares gains of the parts over the slices sl[ch]; returns the model and the gains."""
    mdl = np.zeros_like(y)
    gs = []
    for ch in (0, 1):
        X = np.stack([p[sl[ch], ch] for p in parts], 1)
        c = np.linalg.lstsq(X, y[sl[ch], ch], rcond=None)[0]
        mdl[:, ch] = sum(ci * p[:, ch] for ci, p in zip(c, parts))
        gs.append(c.tolist())
    return mdl, gs


def null_db(y, mdl, a, b):
    return float(10 * np.log10(np.sum((y[a:b] - mdl[a:b]) ** 2) / np.sum(y[a:b] ** 2))) if b - a > 200 else None


def fit_window(prog, has_cluster):
    W = windows(prog)
    sl = [slice(K0 + W[ch]["cluster_gap" if has_cluster else "first_pass"][0], K0 + W[ch]["first_pass"][1]) for ch in (0, 1)]
    return sl, K0 + max(W[ch]["first_pass"][1] for ch in (0, 1)) + 2 * HALF


def fit_warp(y, cols, sl, nf, x0=(0.0, 14.0)):
    """Shift s and clock (ppm) of the columns (each (n, 2)), with per-channel gains, on the slices sl; returns the
    warped full-length model, (s, ppm), gains."""
    def cost(u):
        mdl, _ = gains(y[:nf], [warp(c[:nf], *u) for c in cols], sl)
        return sum(np.sum((y[sl[ch], ch] - mdl[sl[ch], ch]) ** 2) for ch in (0, 1))
    u = minimize(cost, x0, method="Nelder-Mead",
                 options=dict(xatol=1e-3, fatol=1e-14, initial_simplex=[x0, (x0[0] + 0.2, x0[1]), (x0[0], x0[1] + 10)])).x
    wc = [warp(c, *u) for c in cols]
    _, gs = gains(y[:nf], [c[:nf] for c in wc], sl)
    mdl = np.stack([sum(g * c[:, ch] for g, c in zip(gs[ch], wc)) for ch in (0, 1)], 1)
    return mdl, (float(u[0]), float(u[1])), gs, wc


def fit_repeat(y, other, cl_part, M, has_cluster):
    """The repeat-capture model: another capture of the group (clock and shift relative to this one) plus a cluster
    column, with the same freedoms and window as the model."""
    sl, nf = fit_window(M.prog, has_cluster)
    def cost(u):
        mdl, _ = gains(y[:nf], [cl_part[:nf], warp(other[:nf], *u)], sl)
        return sum(np.sum((y[sl[ch], ch] - mdl[sl[ch], ch]) ** 2) for ch in (0, 1))
    u = minimize(cost, (0.0, 0.0), method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-14, initial_simplex=[(0, 0), (0.2, 0), (0, 10)])).x
    o = warp(other, *u)
    _, gs = gains(y[:nf], [cl_part[:nf], o[:nf]], sl)
    return np.stack([gs[ch][0] * cl_part[:, ch] + gs[ch][1] * o[:, ch] for ch in (0, 1)], 1)


def band_decay(y, mdl, t0, t1):
    """Per octave band and channel: slope (dB/s) and level (dB at t0) of 10 log10(E_meas / E_model) per frame."""
    out = {}
    for f, sos in OCT_SOS.items():
        yb, mb = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, mdl, axis=0)
        row = []
        for ch in (0, 1):
            a, b = t0[ch], t1[ch]
            k = (b - a) // FRAME
            if k < 4:
                row.append(None)
                continue
            ey = (yb[a:a + k * FRAME, ch] ** 2).reshape(k, FRAME).mean(1)
            em = (mb[a:a + k * FRAME, ch] ** 2).reshape(k, FRAME).mean(1)
            ok = (ey > ey.max() * 1e-7) & (em > em.max() * 1e-7)
            t = (np.arange(k) + 0.5) * FRAME / FS
            d = 10 * np.log10(ey[ok] / em[ok])
            wgt = np.sqrt(em[ok] / em[ok].max())                  # louder frames count more (fewer of them fluctuate)
            A = np.stack([np.ones(ok.sum()), t[ok]], 1) * wgt[:, None]
            c = np.linalg.lstsq(A, d * wgt, rcond=None)[0]
            row.append([float(c[1]), float(c[0]), float(t[ok][-1])])
        out[f] = row
    return out


# ---------------------------------------------------------------- law stage
P_SCAN = [0.80, 0.82, 0.84, 0.86, 0.88, 0.90, None]
DC_SCAN = [0.10, 0.14, 0.18, 0.22, 0.26, None]
BETA_SCAN = [0.0, 0.019, 0.04]


def variants(prog, key):
    """The line variants of one setting: the law, its two controls, and the limit scans where a limit can act."""
    lf, rt, hf = (float(x) for x in key.split("_"))
    P = cc.LINES[prog]
    raw = line_coefs(P, rt, lf, hf, caps=False)
    active = dict(hf=bool(max(p for _, _, p in raw) > min(x for x in P_SCAN if x)),
                  dc=bool(min(g * r0 for g, r0, _ in raw) < max(x for x in DC_SCAN if x)))
    v = {"law": {}, "nocaps": dict(coefs=raw), "uncoupled": dict(F=np.eye(4))}
    if active["hf"]:
        v.update({f"p{x}": dict(coefs=line_coefs(P, rt, lf, hf, p_max=x)) for x in P_SCAN})
    if active["dc"]:
        v.update({f"dc{x}": dict(coefs=line_coefs(P, rt, lf, hf, dc_min=x)) for x in DC_SCAN})
        v.update({f"beta{x}": dict(coefs=line_coefs(P, rt, lf, hf, beta=x)) for x in BETA_SCAN})
    return v, active


def _law_job(a):
    """One AMPL group. Per capture, shift, clock and two gains per channel are fitted once, on [0, onset + shortest
    line), which no line reaches; every line variant then only swaps the rendered reverb (one render per group)."""
    prog, key, path = a
    M = Model(prog)
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, _, amps = cc.captures(prog, key, n)
    n = len(caps[amps[-1]])
    W = windows(prog)
    var, active = variants(prog, key)
    C, R0 = M.parts(n, key)
    rev = {vn: (R0 if vn == "law" else M.parts(n, key, **kw)[1]) for vn, kw in var.items()}
    rows = []
    for amp in amps:
        y = caps[amp]
        end = usable_end(y, g[amp], prog)
        sl, nf = fit_window(prog, amp > 0)
        _, u, gs, (Cw, _) = fit_warp(y, [C, R0], sl, nf)
        row = dict(key=key, amp=amp, end=[(e - K0) / FS for e in end], clock=u, gains=gs, active=active)
        nulls = lambda mdl: {w: [null_db(y[:, ch], mdl[:, ch], K0 + W[ch][w][0], min(K0 + W[ch][w][1], end[ch])) for ch in (0, 1)]
                             for w in W[0] if amp > 0 or w != "cluster_gap"}
        for vn, Rv in rev.items():
            Rw = warp(Rv, *u)
            mdl = np.stack([gs[ch][0] * Cw[:, ch] + gs[ch][1] * Rw[:, ch] for ch in (0, 1)], 1)
            row[vn] = nulls(mdl)
            if vn == "law":
                row["bands"] = band_decay(y, mdl, [K0 + cc.ONSET[prog][ch] for ch in (0, 1)], end)
        cl_part = Cw * np.array([gs[ch][0] for ch in (0, 1)])[None]                     # the fitted model's cluster
        other = max((x for x in amps if x != amp), key=lambda x: min(usable_end(caps[x], g[x], prog)))   # the longest other capture
        rep = fit_repeat(y, caps[other], cl_part, M, amp > 0)
        oe = usable_end(caps[other], g[other], prog)
        row["repeat"] = {w: [null_db(y[:, ch], rep[:, ch], K0 + W[ch][w][0], min(K0 + W[ch][w][1], end[ch], oe[ch])) for ch in (0, 1)]
                         for w in W[0] if amp > 0 or w != "cluster_gap"}
        row["repeat_with"] = other
        rows.append(row)
    json.dump(rows, open(path, "w"))
    return prog, key


def group_dir(tag=""):
    return os.path.join(OUT, f"law{tag}")


def group_path(prog, key, tag=""):
    return os.path.join(group_dir(tag), f"{prog}_{key}.json")


def law_stage(progs, workers, only=None, chunk=None):
    """Every AMPL group of the programs (smoke: `only` per program), one checkpoint file each; groups already done are
    skipped, so an interrupted run resumes. chunk = (i, n): only every n-th group from the i-th."""
    tag = "_smoke" if only else ""
    os.makedirs(group_dir(tag), exist_ok=True)
    jobs = []
    for p in progs:
        ks = sorted(cc.groups(p))
        if only:
            ks = ks[::max(1, len(ks) // only)][:only]
        jobs += [(p, k, group_path(p, k, tag)) for k in ks]
    if chunk:
        jobs = jobs[chunk[0]::chunk[1]]                      # chunks of the full list, so they don't shift as groups finish
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    print(f"law: {len(jobs)} groups to do", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:     # a fresh process per group returns its memory
        fut = [ex.submit(_law_job, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            p, k = f.result()
            print(f"  {i + 1}/{len(jobs)} {p} {k}", flush=True)


def load_rows(prog, tag=""):
    import glob
    return [r for f in sorted(glob.glob(os.path.join(group_dir(tag), f"{prog}_*.json"))) for r in json.load(open(f))]


def fitted(prog, key, tag="", **variant):
    """Rebuild the law stage's result for one AMPL group without refitting: per AMPL, (capture, model, usable end), on
    the K0 grid, band-limited, with the stored clock and gains. variant: coefs= / F= as in Model.parts (default: the
    law). The residual y - model past the first pass is what Q3b works on (ml_structure §8.4)."""
    rows = {r["amp"]: r for r in json.load(open(group_path(prog, key, tag)))}
    M = Model(prog)
    g = cc.groups(prog)[key]
    caps, _, amps = cc.captures(prog, key, max(len(cc._IR(f)) for f in g.values()) + 200)
    n = len(caps[amps[-1]])
    C, R = M.parts(n, key, **variant)
    out = {}
    for amp in amps:
        r = rows[amp]
        Cw, Rw = warp(C, *r["clock"]), warp(R, *r["clock"])
        mdl = np.stack([r["gains"][ch][0] * Cw[:, ch] + r["gains"][ch][1] * Rw[:, ch] for ch in (0, 1)], 1)
        out[amp] = (caps[amp], mdl, usable_end(caps[amp], g[amp], prog))
    return out


def rt_law_bands(prog, key):
    lf, rt, hf = (float(x) for x in key.split("_"))
    return rt_of(line_coefs(cc.LINES[prog], rt, lf, hf), cc.LINES[prog], OCT)


def _split(model_db, repeat_db):
    """Capture-to-capture variance v and the model's own error e (both re the capture, dB), assuming each capture is
    a common response plus an independent random part: repeat = 2 v, model = v + e."""
    m, r = 10 ** (np.asarray(model_db) / 10), 10 ** (np.asarray(repeat_db) / 10)
    v = r / 2
    return 10 * np.log10(v), 10 * np.log10(np.maximum(m - v, 1e-6))


def report_stage(progs, tag=""):
    """Breakdowns of law.json: nulls by control against the repeat captures; band RT ratios by LF and HF step."""
    for p in progs:
        rows = load_rows(p, tag)
        if not rows:
            continue
        lab = lambda r: dict(zip(("lf", "rt", "hf"), (float(x) for x in r["key"].split("_"))))
        print(f"\n## {p}: {len({r['key'] for r in rows})} groups, {len(rows)} captures")
        print("- nulls, median L / R; [repeat]; v / e = capture-to-capture variance / the model's own error (R), dB")
        for w in windows(p)[0]:
            sel = [r for r in rows if w in r["law"] and r["law"][w][1] is not None and r["repeat"][w][1] is not None]
            if not sel:
                continue
            m = np.array([r["law"][w] for r in sel], float)
            q = np.array([r["repeat"][w] for r in sel], float)
            v, e = _split(np.median(m[:, 1]), np.median(q[:, 1]))
            u = np.array([r["uncoupled"][w][1] for r in sel])
            nc = np.array([r["nocaps"][w][1] for r in sel])
            print(f"  {w:11s} law {np.nanmedian(m[:, 0]):6.1f} / {np.median(m[:, 1]):6.1f} [{np.nanmedian(q[:, 0]):6.1f} / {np.median(q[:, 1]):6.1f}]"
                  f"  v {v:6.1f} e {e:6.1f}; no limits R {np.median(nc):6.1f}; uncoupled R {np.median(u):6.1f}; n {len(sel)}")
        for ctl in ("rt", "lf", "hf"):
            vals = sorted({lab(r)[ctl] for r in rows})
            print(f"- by {ctl.upper()} (R; law [repeat]; recirc1 / recirc2; no-limits recirc1):")
            for x in vals:
                sel = [r for r in rows if lab(r)[ctl] == x]
                cells = []
                for w in ("recirc1", "recirc2"):
                    a = [r["law"][w][1] for r in sel if w in r["law"] and r["law"][w][1] is not None and r["repeat"][w][1] is not None]
                    b = [r["repeat"][w][1] for r in sel if w in r["law"] and r["law"][w][1] is not None and r["repeat"][w][1] is not None]
                    cells.append(f"{np.median(a):6.1f} [{np.median(b):6.1f}]" if a else "     -         ")
                nc = [r["nocaps"]["recirc1"][1] for r in sel if r["nocaps"].get("recirc1", [None, None])[1] is not None]
                print(f"  {x:5}: " + " / ".join(cells) + (f"; {np.median(nc):6.1f}" if nc else "") + f"  ({len(sel)} captures)")
        print("- RT ratio measured / law per octave (median over captures and channels), by LF and HF step:")
        print("  " + " " * 12 + " ".join(f"{f:>6d}" for f in OCT))
        for ctl in ("lf", "hf"):
            for x in sorted({lab(r)[ctl] for r in rows}):
                sel = [r for r in rows if lab(r)[ctl] == x]
                cells = []
                for f in OCT:
                    rr = []
                    for r in sel:
                        rtl = rt_law_bands(p, r["key"])[OCT.index(f)]
                        for ch in (0, 1):
                            b = r["bands"][str(f)][ch]
                            if b is not None and 60 / rtl - b[0] > 0:
                                rr.append(60 / (60 / rtl - b[0]) / rtl)
                    cells.append(f"{np.median(rr):6.3f}" if rr else "     -")
                print(f"  {ctl.upper()} {x:<8} " + " ".join(cells))
        print("- band level re the model at the onset (dB, median), by HF step:")
        for x in sorted({lab(r)["hf"] for r in rows}):
            sel = [r for r in rows if lab(r)["hf"] == x]
            print(f"  HF {x:<6} " + " ".join(f"{np.median([r['bands'][str(f)][ch][1] for r in sel for ch in (0, 1) if r['bands'][str(f)][ch] is not None]):+6.2f}" for f in OCT))


def report_caps(progs, tag=""):
    """The coefficient-limit scans: recirculation nulls per limit value, over the groups where that limit can act."""
    for p in progs:
        rows = load_rows(p, tag)
        for name, scan, act in (("HF pole limit", [f"p{x}" for x in P_SCAN], "hf"), ("per-pass DC gain floor", [f"dc{x}" for x in DC_SCAN], "dc"),
                                ("beta", [f"beta{x}" for x in BETA_SCAN], "dc")):
            sel = [r for r in rows if r["active"][act]]
            if not sel:
                continue
            print(f"\n## {p}: {name}, {len({r['key'] for r in sel})} groups where it can act; recirc1 / recirc2 / recirc3 null, "
                  f"median over captures and channels (dB)")
            for vn in scan:
                cells = []
                for w in ("recirc1", "recirc2", "recirc3"):
                    x = [r[vn][w][ch] for r in sel for ch in (0, 1) if w in r[vn] and r[vn][w][ch] is not None]
                    cells.append(f"{np.median(x):6.2f}" if x else "     -")
                print(f"  {vn:9s} " + " / ".join(cells))


# ---------------------------------------------------------------- fit stage
T_MIN = 0.15                    # a band's decay counts only if it was followed this long (s) before its floor or the fade


def _fit_job(a):
    """One group: the line law's three numbers v = (rt_hi, dc ratio, HF factor), no limits, fitted to the measured
    decay rate per octave of every capture (target = law rate - measured slope; the band-centre approximation cancels
    against the law model the slopes were measured on). Returns v, its residual, and the per-line coefficients."""
    from scipy.optimize import least_squares
    p, key, rows = a
    lf, rt, hf = (float(x) for x in key.split("_"))
    P = cc.LINES[p]
    law_rate = 60 / rt_law_bands(p, key)
    f, tgt = [], []
    for r in rows:
        for i, fc in enumerate(OCT):
            for ch in (0, 1):
                b = r["bands"][str(fc)][ch]
                if b is not None and b[2] >= T_MIN:
                    f.append(fc)
                    tgt.append(law_rate[i] - b[0])
    f, tgt = np.array(f), np.array(tgt)
    if len(set(f)) < 4:
        return p, key, None
    res = lambda v: (-60 / rt_of(line_coefs(P, rt, lf, hf, caps=False, v=v), P, f) - (-tgt)) / tgt
    L = lf_law(lf)
    x0 = np.array([L["hi"] * rt, L["dc"] / L["hi"], min(hf, 0.999)])
    s = least_squares(res, x0, bounds=([0.05 * rt, 0.02, 0.02], [5 * rt, 5.0, 1.0]), loss="soft_l1", f_scale=0.02)
    co = line_coefs(P, rt, lf, hf, caps=False, v=s.x)
    law = line_coefs(P, rt, lf, hf, caps=True)
    return p, key, dict(v=s.x.tolist(), rms=float(np.sqrt(np.mean(s.fun ** 2))), n=int(len(f)), bands=sorted(set(f.tolist())),
                        dc_gain=[g * r0 for g, r0, _ in co], hf_pole=[pp for _, _, pp in co],
                        law_dc_gain=[g * r0 for g, r0, _ in law], law_hf_pole=[pp for _, _, pp in law],
                        rt_fit=rt_of(co, P, OCT).tolist(), rt_law=rt_law_bands(p, key).tolist())


def fit_stage(progs, workers, tag=""):
    jobs = []
    for p in progs:
        by = {}
        for r in load_rows(p, tag):
            by.setdefault(r["key"], []).append(r)
        jobs += [(p, k, rr) for k, rr in by.items()]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_fit_job, jobs))
    out = {p: {k: v for q, k, v in res if q == p and v is not None} for p in progs}
    json.dump(out, open(os.path.join(OUT, f"fit{tag}.json"), "w"), indent=0)
    for p in progs:
        o = out[p]
        lab = {k: dict(zip(("lf", "rt", "hf"), (float(x) for x in k.split("_")))) for k in o}
        print(f"\n## {p}: line law fitted per setting to the band decay ({len(o)} groups; median rms rate error "
              f"{np.median([x['rms'] for x in o.values()]) * 100:.2f} %)")
        print("- per LF step, settings with HF at its top step: Nyquist-end RT / label (S law: L['hi']), DC RT / Nyquist-end RT "
              "(S law: dc / hi); median [IQR] over RT; RT where no limit is active in the law")
        hf_top = max(l["hf"] for l in lab.values())
        for lf in sorted({l["lf"] for l in lab.values()}):
            ks = [k for k in o if lab[k]["lf"] == lf and lab[k]["hf"] == hf_top]
            free = [k for k in ks if min(o[k]["law_dc_gain"]) > DC_MIN * 1.0001]
            hi = np.array([o[k]["v"][0] / lab[k]["rt"] for k in free])
            dc = np.array([o[k]["v"][1] for k in free])
            L = lf_law(lf)
            if len(hi):
                print(f"  LF {lf}: hi {np.median(hi):.3f} [{np.percentile(hi, 25):.3f}–{np.percentile(hi, 75):.3f}] (S {L['hi']:.3f}); "
                      f"dc/hi {np.median(dc):.3f} [{np.percentile(dc, 25):.3f}–{np.percentile(dc, 75):.3f}] (S {L['dc'] / L['hi']:.3f}); "
                      f"{len(free)} of {len(ks)} settings uncapped")
        print("- per-pass DC gain (min over lines) fitted vs law-with-limits, settings where the law's DC limit is active:")
        for k in sorted(o, key=lambda k: (lab[k]["lf"], lab[k]["rt"], lab[k]["hf"])):
            if min(o[k]["law_dc_gain"]) <= DC_MIN * 1.0001:
                print(f"  {k}: fitted {min(o[k]['dc_gain']):.3f} (law {min(o[k]['law_dc_gain']):.3f}); RT fit / law per octave "
                      + " ".join(f"{a / b:.2f}" for a, b in zip(o[k]["rt_fit"], o[k]["rt_law"])))
        print("- HF pole (max over lines) fitted vs law-with-limits, per HF step (median over settings) and where the law's pole is capped:")
        for h in sorted({l["hf"] for l in lab.values()}):
            ks = [k for k in o if lab[k]["hf"] == h]
            capped = [k for k in ks if max(o[k]["law_hf_pole"]) >= P_MAX - 1e-9]
            print(f"  HF {h}: fitted {np.median([max(o[k]['hf_pole']) for k in ks]):.3f} (law {np.median([max(o[k]['law_hf_pole']) for k in ks]):.3f}); "
                  f"capped in the law: {len(capped)} settings, fitted max pole there {np.median([max(o[k]['hf_pole']) for k in capped]) if capped else float('nan'):.3f}")
        print("- RT fit / law per octave (median over all groups), by HF step:")
        for h in sorted({l["hf"] for l in lab.values()}):
            ks = [k for k in o if lab[k]["hf"] == h]
            print(f"  HF {h:<5}" + " ".join(f"{np.median([o[k]['rt_fit'][i] / o[k]['rt_law'][i] for k in ks]):6.3f}" for i in range(len(OCT))))


# ---------------------------------------------------------------- main
def main():
    a = argparse.ArgumentParser()
    a.add_argument("--stage", default="law", help="law | report | fit | all")
    a.add_argument("--prog", default="ML")
    a.add_argument("--workers", type=int, default=4)
    a.add_argument("--only", type=int, default=None, help="smoke test: this many groups per program")
    a.add_argument("--chunk", default=None, help="i/n: run only every n-th group from the i-th")
    o = a.parse_args()
    os.makedirs(OUT, exist_ok=True)
    tag = "_smoke" if o.only else ""
    if o.stage in ("law", "all"):
        law_stage(o.prog, o.workers, o.only, tuple(int(x) for x in o.chunk.split("/")) if o.chunk else None)
    if o.stage in ("report", "all"):
        report_stage(o.prog, tag)
        report_caps(o.prog, tag)
    if o.stage in ("fit", "all"):
        fit_stage(o.prog, o.workers, tag)


if __name__ == "__main__":
    main()
