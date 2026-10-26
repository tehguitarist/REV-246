"""Program E (the 246's EMT 250 program) as a structure, fitted to the Altiverb REV-250 presets and scored held out by RT step
(docs/EMT250_program.md §11).

    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/e250/e_model.py --stage prep
    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/e250/e_model.py --stage static --tfit 1.0 [--law 246] [--held 2.5]

Static model, per channel, at the device rate 26,373.6 Hz (R's tap 1 at K0):
  y = K_ch * [ s_ch direct + AP46 AP146 AP392 sum_k c_k G_k z^-d_k / (1 + G_k z^-P_k) ],
  direct = taps 0 / 429 / 718 (L: 314 later, inverted), loops P = 2215 / 2507 / 2776 with negative feedback, first reads
  d = 2223 / 805 / 1802 (L's reverb 315 later), G_k one pass of the loop filter (T60 = RT m(f), or the 246's law at panel
  LF / HF), K_ch the chain kernel (tap 1's own response, from the training presets). 16 numbers (11 with --law 246).
The direct path is fitted first on its own windows, then the reverb. The reverb is time variant (§11.1: each preset carries a
fresh random part, ~72 % of its energy), so the reverb is scored as the fraction of a held-out capture's energy explained (rho^2)
against the leave-one-out mean of the other presets; the direct path as a null.
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import butter, sosfiltfilt, firwin, fftconvolve

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FS = 19.2e6 / 728
R48 = 48000
RTS = [0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5]
K0 = 64                                   # device index of R's first tap
OCT = np.array([125, 250, 500, 1000, 2000, 4000, 8000.0])
MB = np.array([1.40, 1.26, 1.11, 1.00, 0.92, 0.82, 0.67])          # UAD250_test_plan §2.9 (R; L the same to 0.01)
BAND = butter(4, [100, 10000], btype="band", fs=FS, output="sos")
CACHE = os.path.join(ROOT, "analysis", "altiverb", "cache")
OUT = os.path.join(ROOT, "analysis", "e250")


# ---------------------------------------------------------------- data
def to_device(x48, ppm, t0, n):
    """48 kHz (n48, 2) -> device grid, sample j at 48 kHz time t0 + j R48 / FS (1 + ppm 1e-6); low-passed at 12.4 kHz first."""
    lp = firwin(1023, 12400, fs=R48, window=("kaiser", 12.0))
    x = fftconvolve(x48, lp[:, None], axes=0)[511:511 + len(x48)]
    t = t0 + np.arange(n) * R48 / FS * (1 + ppm * 1e-6)
    i = np.floor(t).astype(int); fr = t - i; half = 64; k = np.arange(-half + 1, half + 1)
    xp = np.concatenate([np.zeros((half, 2)), x, np.zeros((half + 2, 2))])
    out = np.zeros((n, 2))
    for s0 in range(0, n, 4096):
        dd = fr[s0:s0 + 4096, None] - k[None, :]
        w = np.sinc(dd) * np.i0(10 * np.sqrt(np.clip(1 - (dd / (half + 0.5)) ** 2, 0, None))) / np.i0(10)
        idx = np.clip(i[s0:s0 + 4096, None] + k[None, :] + half, 0, len(xp) - 1)
        out[s0:s0 + 4096] = np.einsum("nk,nkc->nc", w, xp[idx])
    return out


def load48(rt):
    h = np.load(os.path.join(CACHE, f"REV-250 {rt}s.npy")).T          # (n, 2): ch 0 = L (late onset), ch 1 = R
    a = np.abs(h[:, 1]); i = int(np.argmax(a > 0.5 * a.max()))
    i = i - 3 + int(np.argmax(a[i - 3:i + 4]))
    # sub-sample peak of R's first tap (parabolic on the band-limited peak)
    y0, y1, y2 = h[i - 1:i + 2, 1]; frac = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    return h, i + frac


def stage_prep(ppm=42.8, tlen=6.0):
    os.makedirs(OUT, exist_ok=True)
    for rt in RTS:
        h, tp = load48(rt)
        n = int(tlen * FS)
        y = to_device(h, ppm, tp - K0 * R48 / FS, n)
        np.save(os.path.join(OUT, f"dev_{rt}.npy"), y)
        print(f"RT {rt}: R tap at {tp:.2f} (48k), device peak R {y[K0, 1]:+.1f} L-onset {np.argmax(np.abs(y[:, 0]) > 0.5 * np.abs(y[:, 0]).max()) - K0}")


_DEV = {}


def dev(rt):
    """Device-rate IR normalised by R's first tap: the presets carry a per-preset gain (tap 367 -> 264 from RT 0.4 to 4.5
    while the early reverb / tap ratio stays at -33.5 +- 0.4 dB), i.e. a capture / preset level, not part of the program."""
    if rt not in _DEV:
        y = np.load(os.path.join(OUT, f"dev_{rt}.npy")); _DEV[rt] = y / y[K0, 1]
    return _DEV[rt]


# ---------------------------------------------------------------- loop filter
def loss_db(f, rt, n, law):
    """Per-pass loss in dB of a loop of length n at RT (the label) at frequency f: T60(f) = RT m(f)."""
    m = np.interp(np.log(np.clip(f, OCT[0], OCT[-1])), np.log(OCT), law)
    return 60.0 * n / (FS * rt * m)


def minphase(ldb, nt=96, M=2048):
    """Minimum-phase FIR with magnitude 10^(-ldb/20), ldb on rfftfreq(M)."""
    c = np.fft.irfft(-ldb / 20 * np.log(10), M)
    fold = np.zeros(M); fold[0] = c[0]; fold[1:M // 2] = 2 * c[1:M // 2]; fold[M // 2] = c[M // 2]
    return np.fft.irfft(np.exp(np.fft.rfft(fold)), M)[:nt]


# ---------------------------------------------------------------- static structure (EMT250_program §11)
P = (2215, 2507, 2776)                    # loop lengths, negative feedback (the old 4430 / 5014 / 5552 are 2 P)
DK = (2223, 805, 1802)                    # each loop's first read after the direct tap (A, B, C; all-pass chain and chain kernel divided out)
AP_D = (46, 146, 392)                     # all-pass chain, outside the loops, shared by static and random parts (g ~ 0.70 each)
TAPS = (0, 429, 718)                      # direct path
LDEL = 314                                # L's direct path: R's 314 samples later, inverted
LREV = 315                                # L's reverb: R's static reverb 315 samples later, same sign
KL = 410                                  # chain kernel length (tap 1's response; tap 2 starts at 429 - 16)
NFS = 1 << 18                             # 9.9 s: no circular wrap of RT 4.5 tails into the direct window
W = np.fft.rfftfreq(NFS)
FK = np.fft.rfftfreq(2048, 1 / FS)


KLONG = 540                               # §16: the kernel past tap 2 (tap 2's copy removed recursively) up to the first moving read


def kernels(rts, a2=-0.709):
    """Per channel, tap 1's response (chain), from the mean of the given presets; normalised to 1 at the tap. Up to KL
    samples it is the data itself; from tap 2 (429) on, tap 2's own copy of the kernel is subtracted recursively,
    k[n] = d[n] - a2 k[n - 429], which extends it to KLONG (the first moving reads start ~540; §16)."""
    m = np.mean([dev(r) for r in rts], 0); out = []
    for ch, o in ((0, LDEL), (1, 0)):
        d = m[K0 + o - 16:K0 + o + KLONG, ch].copy(); d /= d[16]; k = d.copy()
        for i in range(TAPS[1], len(k)):        # tap 2's copy starts 16 samples before it (index 16 = the tap)
            k[i] = d[i] - a2 * k[i - TAPS[1]]
        out.append(k)
    if TAIL and os.path.exists(TAIL_FILE):              # §16: the chain's slow tail, a cubic B-spline from KLONG to 4000
        tl = json.load(open(TAIL_FILE)); lags = np.arange(tl["lo"], tl["hi"])
        for ch, key in ((0, "L"), (1, "R")):
            tail = np.sum([b * _bspl(tl, j)(lags) for j, b in enumerate(tl[key])], 0)
            out[ch] = np.r_[out[ch], np.zeros(16 + tl["hi"] - len(out[ch]))]; out[ch][16 + tl["lo"]:16 + tl["hi"]] += tail
    return out


TAIL = os.environ.get("E_TAIL", "1") == "1"
TAIL_FILE = os.path.join(OUT, "kernel_tail.json")


def _bspl(tl, j):
    from scipy.interpolate import BSpline
    t = np.r_[[tl["lo"]] * 3, np.linspace(tl["lo"], tl["hi"], tl["nb"] - 2), [tl["hi"]] * 3]
    return BSpline(t, np.eye(tl["nb"])[j], 3)


def fit_tail(x, rts, lo=540, hi=4000, nb=16, secs=1.0):
    """The kernel's slow tail (§16): per channel, nb cubic B-spline coefficients on lags lo..hi, fitted by least squares
    below 250 Hz on E's expected response (render_mean) over the first `secs`, with the kernel without its tail."""
    global TAIL
    from scipy.signal import butter, sosfiltfilt
    keep = TAIL; TAIL = False
    try:
        kern = kernels(rts); n = K0 + int(secs * FS); lp = butter(4, 250, fs=FS, output="sos"); lags = np.arange(lo, hi)
        tl = dict(lo=lo, hi=hi, nb=nb); Bm = np.array([_bspl(tl, j)(lags) for j in range(nb)])
        for ch, key in ((0, "L"), (1, "R")):
            A, b = [], []
            for rt in rts:
                d = np.zeros(KLONG + 16); d[16] = 1.0
                u = render_mean(x, rt, [d, d.copy()], n)[:, ch]; y = dev(rt)[:n, ch]; m = render_mean(x, rt, kern, n)[:, ch]
                A.append(np.array([sosfiltfilt(lp, np.convolve(u, np.r_[np.zeros(lo), Bm[j]])[:n])[K0 - 20:] for j in range(nb)]).T)
                b.append(sosfiltfilt(lp, y - m)[K0 - 20:])
            tl[key] = np.linalg.lstsq(np.vstack(A), np.concatenate(b), rcond=None)[0].tolist()
    finally:
        TAIL = keep
    return tl


LAWMODE = "mb"                            # "mb": T60(f) = RT m(f), 7 knots; "246": the 246's loop law (ml_lines.line_coefs) at panel LF / HF


def unpack(x):
    return dict(a=x[0:2], gL=x[2], c=x[3:6], g=x[6:9], m=x[9:16] if LAWMODE == "mb" else x[9:11])


def loop_G(rt, Pk, law):
    """One pass of loop Pk's filter on the rfft grid W."""
    if LAWMODE == "246":
        import ml_lines as mll
        return mll.line_resp(mll.line_coefs([Pk], rt, law[0], law[1]), np.exp(-2j * np.pi * W))[:, 0]
    return np.fft.rfft(minphase(loss_db(FK, rt, Pk, law)), NFS)


X0 = np.r_[-0.709, 0.59, 1.1146, 0.324, 0.414, -0.419, 0.7, 0.7, 0.7, MB]


def render(x, rt, kern, n):
    """Static model at RT, (n, 2) device samples, R's tap 1 at K0."""
    q = unpack(x); z = lambda d: np.exp(-2j * np.pi * W * d)
    direct = 1 + q["a"][0] * z(TAPS[1]) + q["a"][1] * z(TAPS[2])
    ap = np.ones_like(W, dtype=complex)
    for D, g in zip(AP_D, q["g"]):
        ap *= (-g + z(D)) / (1 - g * z(D))
    rev = 0
    for Pk, dk, ck in zip(P, DK, q["c"]):
        G = loop_G(rt, Pk, q["m"])
        rev = rev + ck * G * z(dk) / (1 + G * z(Pk))
    rev = rev * ap
    out = np.zeros((n, 2))
    for ch, (sd, gc, o, orv) in enumerate(((-1.0, q["gL"], LDEL, LREV), (1.0, 1.0, 0, 0))):
        Hk = np.fft.rfft(kern[ch], NFS) * z(o - 16 + K0)
        out[:, ch] = gc * np.fft.irfft(Hk * (sd * direct + rev * z(orv - o)), NFS)[:n]
    return out


def render_mean(x, rt, kern, n, reads=None, box=70):
    """E's expected response (§16): the static network plus each fixed moving read averaged over its glide (a centred boxcar of
    `box` samples = the triangle's full excursion), with its own R / L amplitude (L 314 samples later). The deterministic part
    of every capture; what a single capture can be nulled against."""
    if reads is None:
        reads = json.load(open(os.path.join(OUT, "fixed_reads.json")))["reads"]
    q = unpack(x); z = lambda d: np.exp(-2j * np.pi * W * d)
    ap = np.ones_like(W, dtype=complex)
    for D, g in zip(AP_D, q["g"]):
        ap *= (-g + z(D)) / (1 - g * z(D))
    bx = np.fft.rfft(np.ones(box) / box, NFS) * z(-(box - 1) / 2)
    out = render(x, rt, kern, n)
    G = {k: loop_G(rt, P[k], q["m"]) for k in range(3)}
    for ch, o, key in ((0, LDEL, "L"), (1, 0, "R")):
        Hk = np.fft.rfft(kern[ch], NFS) * z(o - 16 + K0)
        acc = 0
        for r in reads:
            k = "ABC".index(r["loop"]); acc = acc + r[key] * G[k] / (1 + G[k] * z(P[k])) * z(r["centre"])
        out[:, ch] += np.fft.irfft(Hk * ap * bx * acc, NFS)[:n]
    return out


def resid(x, rts, kern, n, w):
    return np.concatenate([((render(x, rt, kern, n) - dev(rt)[:n]) * w[:, None]).ravel() for rt in rts])


def fitstatic(rts, tfit, x0=None):
    x0 = X0 if x0 is None else x0
    from scipy.optimize import least_squares
    kern = kernels(rts); n = K0 + int(tfit * FS)
    t = np.maximum(np.arange(n) - K0, 0) / FS; w = 10 ** (np.minimum(t / tfit, 1) * 20 / 20)   # mild late weight (+20 dB at tfit)
    # 1. the direct path (a, gL) on its own window (to 19 ms after each channel's tap 1: random reads start ~20 ms)
    # (and tap 3's own 702-760 window, where it is ~30 dB above the random reads)
    nd = K0 + LDEL + 760; wd = np.zeros((nd, 2))
    for ch, o in ((0, LDEL), (1, 0)):
        wd[K0 + o - 20:K0 + o + 500, ch] = 1; wd[K0 + o + 702:K0 + o + 760, ch] = 1
    def rdir(v):
        xx = x0.copy(); xx[0:3] = v
        return np.concatenate([((render(xx, rt, kern, nd) - dev(rt)[:nd]) * wd).ravel() for rt in rts])
    x0 = x0.copy(); x0[0:3] = least_squares(rdir, x0[0:3], diff_step=1e-5).x
    # 2. the reverb (c, g, m) with the direct path fixed (LF held at E_LF when set: §16, LF from the decay statistics)
    free = np.arange(3, len(x0))
    if LAWMODE == "246" and os.environ.get("E_LF"):
        x0[9] = float(os.environ["E_LF"]); free = free[free != 9]
    def rrev(v):
        xx = x0.copy(); xx[free] = v
        return resid(xx, rts, kern, n, w)
    sol = least_squares(rrev, x0[free], x_scale=0.1, diff_step=1e-4, max_nfev=60)
    x = x0.copy(); x[free] = sol.x
    return x, kern, sol


def wnull(y, m, i0, i1, band=None):
    band = BAND if band is None else band
    yb, mb = sosfiltfilt(band, y), sosfiltfilt(band, m)
    return float(10 * np.log10(np.sum((yb[i0:i1] - mb[i0:i1]) ** 2) / np.sum(yb[i0:i1] ** 2)))


def expl(y, m, i0, i1):
    """Fraction of y's band energy in [i0, i1) the predictor m explains (optimal gain): rho^2."""
    yb, mb = sosfiltfilt(BAND, y)[i0:i1], sosfiltfilt(BAND, m)[i0:i1]
    return float(np.dot(yb, mb) ** 2 / (np.dot(yb, yb) * np.dot(mb, mb)))


def stage_static(args):
    rts = [r for r in RTS if r >= args.rtmin]
    res = {}
    held = rts if args.held == "all" else [float(v) for v in args.held.split(",")]
    for rh in held:
        tr = [r for r in rts if r != rh]
        x, kern, sol = fitstatic(tr, args.tfit)
        n = K0 + int(2.0 * FS); y = dev(rh)[:n]; m = render_mean(x, rh, kern, n) if os.path.exists(os.path.join(OUT, "fixed_reads.json")) else render(x, rh, kern, n)
        loo = np.mean([dev(r)[:n] for r in tr], 0)
        r = dict(x=x.tolist(), nfev=int(sol.nfev))
        r["direct_null"] = [wnull(y[:, 0], m[:, 0], K0 + LDEL - 20, K0 + LDEL + 500), wnull(y[:, 1], m[:, 1], K0 - 20, K0 + 500)]   # to 19 ms: random reads start ~20 ms
        fb = butter(4, [20, 12400], btype="band", fs=FS, output="sos")
        r["direct_null_fullband"] = [wnull(y[:, 0], m[:, 0], K0 + LDEL - 20, K0 + LDEL + 500, fb), wnull(y[:, 1], m[:, 1], K0 - 20, K0 + 500, fb)]
        for a, b in ((0.03, 0.165), (0.17, 0.5), (0.5, 1.5)):
            i0, i1 = K0 + int(a * FS), K0 + int(b * FS)
            r[f"rho2_{int(a*1e3)}-{int(b*1e3)}"] = dict(model=[expl(y[:, c], m[:, c], i0, i1) for c in (0, 1)],
                                                         loo_mean=[expl(y[:, c], loo[:, c], i0, i1) for c in (0, 1)])
        res[str(rh)] = r
        q = unpack(x)
        print(f"held {rh}: direct null L/R {r['direct_null'][0]:.1f} / {r['direct_null'][1]:.1f} dB (20 Hz-12.4 kHz {r['direct_null_fullband'][0]:.1f} / {r['direct_null_fullband'][1]:.1f}) | rho2 model vs LOO-mean (L,R) "
              + " | ".join(f"{k[5:]} {v['model'][0]:.2f},{v['model'][1]:.2f} vs {v['loo_mean'][0]:.2f},{v['loo_mean'][1]:.2f}" for k, v in r.items() if k.startswith("rho2"))
              + f" | a {q['a'].round(3)} gL {q['gL']:.4f} c {q['c'].round(3)} g {q['g'].round(3)} m {q['m'].round(2)}", flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, f"static_{args.tag}.json"), "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--law", default="mb")
    ap.add_argument("--tfit", type=float, default=2.0)
    ap.add_argument("--held", default="all")
    ap.add_argument("--rtmin", type=float, default=1.0)
    ap.add_argument("--ppm", type=float, default=42.8)    # the S presets' clock; the direct taps confirm it (§13.7: +55.1 is worse by 4-5 dB)
    ap.add_argument("--tag", default="base")
    a = ap.parse_args()
    if a.law == "246":
        sys.path.insert(0, os.path.join(ROOT, "tools")); LAWMODE = "246"; X0 = np.r_[X0[:9], 1.1, 0.7]
    if a.stage == "prep":
        stage_prep(a.ppm)
    elif a.stage == "static":
        stage_static(a)
    elif a.stage == "tail":                      # all-data kernel tail on the 246-law static fit (held-out evidence: §16)
        sys.path.insert(0, os.path.join(ROOT, "tools")); LAWMODE = "246"; x = np.array(json.load(open(os.path.join(OUT, "static_all_246.json")))["x"])
        json.dump(fit_tail(x, [r for r in RTS if r >= 1.0]), open(TAIL_FILE, "w"), indent=1); print("wrote", TAIL_FILE)
