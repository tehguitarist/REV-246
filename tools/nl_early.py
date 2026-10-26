"""Program N (Non-Lin): the early section, fitted tap by tap on the Tomkin 96 kHz library (docs/NonLin_program.md §9).

    .venv/bin/python tools/nl_early.py --stage fit [--late 0.0] [--workers 6]

Structure (found in NonLin_program.md §9, device samples, the front-tap frame):
  front   three plain taps per channel, no diffuser: L 41 / 410 / 463, R 41 / 278 / 463 (S's cluster's place);
  early   the input through AP(124), then the main taps (TAPS, R = L - 221, per-channel signs), then AP(180) in L and
          AP(188) in R, all at g = -1/sqrt(2) (S's diffuser gain). EARLY sets one gain per stage (STAGES); SPAN sets how
          many stages run;
  jitter  the AP(124) branch (early and late) sits o samples from the front taps, o per capture (0 in most, else about
          -64, -84, -87): a capture-level time variance, fitted as an alignment.

fit: per file (LATE 0 by default, so nothing but the early section is there), one sub-sample shift for the front taps,
one for the AP(124) branch (o, searched around offsets of a coarse scan), and free least-squares gains for the three
front taps and every main tap of the active stages, per channel, on [0, last tap + 4000) in 100 Hz-10 kHz. Writes
analysis/nl/early_fit.json: per file the shifts, the gains and the residual (in-sample) per channel.
"""
import argparse, json, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.signal import butter, fftconvolve, lfilter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nl_structure as nl  # noqa: E402

FS, K0 = nl.FSD, nl.K0
G_AP = -1 / np.sqrt(2)
BAND = butter(4, [100, 10000], btype="band", fs=FS, output="sos")
FRONT = {0: (41, 410, 463), 1: (41, 278, 463)}
AP_IN, AP_OUT = 124, {0: 180, 1: 188}
R_OFF = -221
R_FIRST_DIRECT = 465            # front-frame position of R's first tap's direct term
TAPS = [685, 1230, 1767, 2523, 2817, 3747, 3835, 4821, 4902, 5775, 5988, 6697, 7005, 7587, 7953, 8445, 8832, 9271,
        9642, 10065, 10383, 10827, 11055, 11557, 11658, 12192, 12255]
SIGNS = {0: [1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1],
         1: [1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1]}
# stage j = the taps sharing one gain (the first stage has one tap)
STAGES = [[0], [1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14], [15, 16], [17, 18], [19, 20], [21, 22],
          [23, 24], [25, 26]]
STAGE_S = 0.0345                # law time per stage, s (stage gains at EARLY -3: 0.69 dB per stage)
N_STAGES = {0.0: 0, 0.1: 2, 0.2: 4, 0.3: 7, 0.4: 10, 0.5: 14}     # stages that run at each library SPAN
OUT = os.path.join(nl.OUT, "early_fit.json")
KFILE = os.path.join(nl.OUT, "out_kernel.json")
KH = int(os.environ.get("NL_KH", 31))   # output kernel half-length: 63 taps (the 12 kHz LOW PASS edge needs more than M / L.s 11)
KTRAIN = [(0.4, 10.0, 0.0), (0.4, -3.0, 0.0), (0.4, 1.0, 0.0), (0.5, 10.0, 0.0), (0.5, -3.0, 0.0), (0.5, 1.0, 0.0)]


def _load_kernel():
    """The per-channel output kernel (stage kernel); NL_KERNEL=0 leaves it out."""
    if os.environ.get("NL_KERNEL") == "0" or not os.path.exists(KFILE):
        return None
    k = json.load(open(KFILE))
    return {0: np.array(k["L"]), 1: np.array(k["R"])}


KER = _load_kernel()


def ap(x, M, g=G_AP):
    return lfilter(np.r_[-g, np.zeros(M - 1), 1.0], np.r_[1.0, np.zeros(M - 1), -g], x)


def fshift(x, s):
    N = 1 << int(np.ceil(np.log2(len(x) + 64)))
    return np.fft.irfft(np.fft.rfft(x, N) * np.exp(-2j * np.pi * np.fft.rfftfreq(N) * s), N)[:len(x)]


def diffuser(ch, n):
    d = np.zeros(n)
    d[0] = 1.0
    return ap(ap(d, AP_IN), AP_OUT[ch])


def taps(ch, n_stage):
    idx = [i for j in range(n_stage) for i in STAGES[j]]
    return [(TAPS[i] + (R_OFF if ch else 0), SIGNS[ch][i], i) for i in idx]


def columns(ch, n, n_stage, o_front, o_early, kappa=0.0):
    """Design matrix (n, 3 + taps): front spikes at o_front + f (1 + kappa), diffused main taps at o_early +
    t (1 + kappa) (kappa: the capture's clock ratio - 1), every position band-limited and fractional."""
    N = 1 << int(np.ceil(np.log2(2 * n + 2 * TAPS[-1])))   # room for taps shifted before the file start
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    cols = [np.fft.irfft(np.exp(ph * (K0 + o_front + f * (1 + kappa))), N)[:n] for f in FRONT[ch]]
    d = diffuser(ch, N // 2)
    Df = np.fft.rfft(d, N)
    for t, s, i in taps(ch, n_stage):
        if ch == 1 and i == 0:
            # R's first tap: its direct (feed-forward) term sits in the front frame at R_FIRST_DIRECT, the delayed
            # terms in the branch frame (NonLin_program.md 9.4)
            c = np.fft.irfft((Df - d[0]) * np.exp(ph * (K0 + o_early + t * (1 + kappa))), N)[:n]
            c += d[0] * np.fft.irfft(np.exp(ph * (K0 + o_front + R_FIRST_DIRECT * (1 + kappa))), N)[:n]
            cols.append(s * c)
            continue
        cols.append(s * np.fft.irfft(Df * np.exp(ph * (K0 + o_early + t * (1 + kappa))), N)[:n])
    A = np.array(cols).T
    if KER is not None:
        A = fftconvolve(A, KER[ch][:, None], axes=0)[KH:KH + n]
    return sosfiltfilt(BAND, A, axis=0)


def solve(y, A):
    g, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ g
    return g, float(np.sum(r * r) / np.sum(y * y))


FRONT_BIAS = -0.0055   # the peak picker's bias on the 410 tap (its neighbours pull it): the front taps fit -0.005...-0.006
                       # samples earlier in every capture and both channels (NonLin_program 11.1); 0-50 ms -42 -> -48...-56 dB


def front_shift(h):
    """Sub-sample position of the L front tap (410) re K0 + 410: the file's trim offset."""
    from core import subsample_peak
    i = K0 + 380 + int(np.argmax(np.abs(h[K0 + 380:K0 + 460, 0])))
    return subsample_peak(h[:, 0], i)[0] - (K0 + 410) + FRONT_BIAS


def stage_db(early, n_stage):
    """Nominal per-stage gain (dB) of the label's law: 60 / EARLY dB/s over STAGE_S per stage, falling from 0 dB at
    the first stage (EARLY > 0), rising to 0 dB at the last (EARLY < 0); EARLY 10 is flat."""
    step = 0.0 if early >= 10 else -60 * STAGE_S / early
    return [step * (j - (n_stage - 1 if early < 0 else 0)) for j in range(n_stage)]


TICK = 0.004                    # the law's time quantum, s: front taps at FRONT_TICKS, stage j at 4 + 9 j
FRONT_TICKS = {0: (0, 1, 3), 1: (0, 2, 3)}   # L 41 / 410 / 463, R 41 / 278 / 463
FRONT_AMP = {0: (1, 2, 1), 1: (1, 2, -1)}        # front taps in main-tap units (signs included)


def law_db(span, early, ch=0):
    """The EARLY law (NonLin_program.md 9.5): gain in dB of the three front taps and of each stage, re a main tap at
    EARLY 10 (flat). Law time tau: front taps at FRONT_TICKS (L 0 / 4 / 12 ms, R 0 / 8 / 12 ms); EARLY > 0: stage j at 16 + 36 j ms, gain -60 tau / EARLY;
    EARLY < 0: the last stage at tau_end = SPAN + 12 ms, the others 36 ms apart before it, gain
    +60 (tau - tau_end) / |EARLY| (the ramp rises over the whole SPAN to the flat level)."""
    ns = N_STAGES[span]
    tf = np.array(FRONT_TICKS[ch]) * TICK
    if early >= 10:
        return np.zeros(3), np.zeros(ns)
    if early > 0:
        ts = (4 + 9 * np.arange(ns)) * TICK
        return -60 * tf / early, -60 * ts / early
    t_end = span + 3 * TICK
    ts = t_end - 9 * TICK * (ns - 1 - np.arange(ns))
    return 60 * (tf - t_end) / -early, 60 * (ts - t_end) / -early


def law_column(ch, n, span, early, o_front, o_early, kappa=0.0):
    """The whole early section of one channel at the law's gains, one column (n,)."""
    ns = N_STAGES[span]
    fdb, sdb = law_db(span, early, ch)
    A = columns(ch, n, ns, o_front, o_early, kappa)
    w = [FRONT_AMP[ch][i] * 10 ** (fdb[i] / 20) for i in range(3)]
    for t, s, i in taps(ch, ns):
        w.append(10 ** (sdb[[j for j, st in enumerate(STAGES) if i in st][0]] / 20))
    return A @ np.array(w)


def coarse_offsets(h, key, fr, k=8):
    """The k best branch offsets of an FFT match of L (both diffuser stages undone, front taps blanked) with the tap
    train at the label's stage gains."""
    ns = N_STAGES[key[0]]
    z = h[:, 0][::-1].copy()
    for M in (AP_IN, AP_OUT[0]):
        z = ap(z, M)
    z = z[::-1]
    z[:K0 + int(fr) + 480] = 0
    db = stage_db(key[1], ns)
    T = np.zeros(TAPS[-1] + 1)
    for t, s, i in taps(0, ns):
        T[t] = s * 10 ** (db[[j for j, st in enumerate(STAGES) if i in st][0]] / 20)
    N = 1 << int(np.ceil(np.log2(len(z) + len(T) + 16)))
    c = np.abs(np.fft.irfft(np.fft.rfft(z, N) * np.conj(np.fft.rfft(T, N)), N))
    lags = np.r_[np.arange(-TAPS[-1], 0), np.arange(0, 400)]
    v = c[lags % N]
    out = []
    for i in np.argsort(v)[::-1]:
        o = int(lags[i]) - K0 - int(round(fr))
        if all(abs(o - q) > 3 for q in out):
            out.append(o)
        if len(out) == k:
            break
    return out


def fit_file(key):
    span = key[0]
    ns = N_STAGES[span]
    h = sosfiltfilt(BAND, nl.ir(key), axis=0)
    fr = front_shift(h)
    last = taps(0, ns)[-1][0] if ns else 600
    n = min(len(h), K0 + int(fr) + last + 4000)
    y = h[:n]

    def cost(o, kappa=0.0, chs=(0, 1)):
        return sum(solve(y[:, ch], columns(ch, n, ns, fr, fr + o, kappa))[1] for ch in chs)
    out = {"front_shift": fr}
    kappa = 0.0
    if ns:
        # coarse: candidates from a template match (thousands of samples apart for steep negative EARLY) and the
        # common small offsets, scored on L; then offset and clock refined on both channels
        grid = sorted({o + d for o in coarse_offsets(h, key, fr) + [0, -64, -84, -87, -163] for d in (-1, 0, 1)})
        cs = [cost(o, chs=(0,)) for o in grid]
        o0 = grid[int(np.argmin(cs))]
        o0 = minimize_scalar(lambda v: cost(v), bounds=(o0 - 1, o0 + 1), method="bounded", options={"xatol": 0.01}).x
        r = minimize(lambda v: np.log(cost(v[0], v[1] * 1e-6)), [o0, 0.0], method="Nelder-Mead",
                     options={"xatol": 1e-3, "fatol": 1e-5, "initial_simplex": [[o0, 0], [o0 + 0.1, 0], [o0, 20]]})
        o, kappa = float(r.x[0]), float(r.x[1] * 1e-6)
    else:
        o = 0.0
    out["kappa_ppm"] = kappa * 1e6
    out["o"] = float(o)
    for ch in (0, 1):
        A = columns(ch, n, ns, fr, fr + o, kappa)
        g, r = solve(y[:, ch], A)
        out[f"front{ch}"] = g[:3].tolist()
        out[f"taps{ch}"] = g[3:].tolist()
        out[f"res{ch}_db"] = 10 * np.log10(r)
    return key, out


def kernel_stage(iters=6):
    """One 11-tap output kernel per channel, shared by every setting, fitted on KTRAIN only (alternating least squares:
    the tap gains given the kernel, the kernel given the gains), with each file's shifts and clock from early_fit.json
    (fitted without a kernel). Every other setting is then held out for it. Writes out_kernel.json."""
    global KER
    KER = None
    fits = json.load(open(OUT))
    data = []
    for key in KTRAIN:
        v, ns = fits[str(key)], N_STAGES[key[0]]
        h = sosfiltfilt(BAND, nl.ir(key), axis=0)
        n = min(len(h), K0 + int(v["front_shift"]) + taps(0, ns)[-1][0] + 4000)
        data.append((h[:n], [columns(ch, n, ns, v["front_shift"], v["front_shift"] + v["o"], v["kappa_ppm"] * 1e-6)
                             for ch in (0, 1)]))
    ker = {}
    for ch in (0, 1):
        k = np.zeros(2 * KH + 1)
        k[KH] = 1.0
        for it in range(iters):
            X, Y, tot, res = [], [], 0.0, 0.0
            for y, As in data:
                A = fftconvolve(As[ch], k[:, None], axes=0)[KH:KH + len(y)]
                g, r = solve(y[:, ch], A)
                tot += np.sum(y[:, ch] ** 2)
                res += r * np.sum(y[:, ch] ** 2)
                m = As[ch] @ g
                X.append(np.array([np.roll(m, j - KH) for j in range(2 * KH + 1)]).T)
                Y.append(y[:, ch])
            print(f"ch {ch} iter {it}: train residual {10 * np.log10(res / tot):.2f} dB", flush=True)
            k = np.linalg.lstsq(np.vstack(X), np.concatenate(Y), rcond=None)[0]
            k /= k.sum()
        ker["LR"[ch]] = k.tolist()
    json.dump(ker, open(KFILE, "w"), indent=1)
    print(ker)


LAWOUT = os.path.join(nl.OUT, "early_law.json")


def law_file(key):
    """The law model held out on one capture: per capture only the alignment (front shift, branch offset o, clock)
    and one gain per channel; nulls (100 Hz-10 kHz) over the front taps [0, 600) and the early section
    [600, last tap + 3000), plus the free per-tap fit's residual for comparison."""
    span, early, _ = key
    ns = N_STAGES[span]
    free = json.load(open(OUT))[str(key)]
    h = sosfiltfilt(BAND, nl.ir(key), axis=0)
    fr = free["front_shift"]
    last = taps(0, ns)[-1][0] if ns else 600
    n = min(len(h), K0 + int(fr) + last + 3000)
    y = h[:n]

    def resid(o, kappa):
        out = []
        for ch in (0, 1):
            m = law_column(ch, n, span, early, fr, fr + o, kappa)
            a = (y[:, ch] @ m) / (m @ m)
            out.append((y[:, ch] - a * m, a))
        return out

    def cost(v):
        return np.log(sum(np.sum(r ** 2) for r, _ in resid(v[0], v[1] * 1e-6)))
    o, kp = free["o"], free["kappa_ppm"]
    if ns:
        r = minimize(cost, [o, kp], method="Nelder-Mead",
                     options={"xatol": 1e-3, "fatol": 1e-6, "initial_simplex": [[o, kp], [o + 0.05, kp], [o, kp + 2]]})
        o, kp = float(r.x[0]), float(r.x[1])
    out = {"o": o, "kappa_ppm": kp}
    i0, i1 = K0 + int(fr), K0 + int(fr) + 600
    for ch, (res, a) in enumerate(resid(o, kp * 1e-6)):
        e = y[:, ch]
        nul = lambda s0, s1: float(10 * np.log10(np.sum(res[s0:s1] ** 2) / np.sum(e[s0:s1] ** 2)))
        out[f"gain{ch}"] = float(a)
        out[f"front{ch}_db"] = nul(i0 + FRONT[0][0] - 20, i1)
        out[f"early{ch}_db"] = nul(i1, n) if ns else None
        out[f"all{ch}_db"] = nul(i0 + FRONT[0][0] - 20, n)   # from just before the first front tap: the trim cut its pre-ringing
        out[f"free{ch}_db"] = free[f"res{ch}_db"]
    return key, out


def law_stage(workers):
    keys = [eval(k) for k in json.load(open(OUT))]
    res = {}
    with ProcessPoolExecutor(workers) as ex:
        for k, out in ex.map(law_file, keys):
            res[str(k)] = out
            print(k, f"law L {out['all0_db']:.1f} R {out['all1_db']:.1f} (free {out['free0_db']:.1f} / "
                     f"{out['free1_db']:.1f}) o {out['o']:+.2f}", flush=True)
    json.dump(res, open(LAWOUT, "w"), indent=1)


def fit_stage(late, workers, only):
    keys = sorted(k for k in nl.table() if k[2] == late)
    if only:
        keys = keys[:: max(1, len(keys) // only)][:only]
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    todo = [k for k in keys if str(k) not in res]
    with ProcessPoolExecutor(workers) as ex:
        for k, out in ex.map(fit_file, todo):
            res[str(k)] = out
            print(k, f"o {out['o']:+.2f} clock {out['kappa_ppm']:+.1f} ppm res L {out['res0_db']:.1f} R {out['res1_db']:.1f}", flush=True)
    os.makedirs(nl.OUT, exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--stage", default="fit")
    ap_.add_argument("--late", type=float, default=0.0)
    ap_.add_argument("--workers", type=int, default=6)
    ap_.add_argument("--only", type=int, default=0)
    a = ap_.parse_args()
    if a.stage == "fit":
        fit_stage(a.late, a.workers, a.only)
    elif a.stage == "kernel":
        kernel_stage()
    elif a.stage == "law":
        law_stage(a.workers)


if __name__ == "__main__":
    main()
