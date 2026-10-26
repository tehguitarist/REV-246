"""Task 1: is the Altiverb "dark chain" (REV-250 and LONG RT presets) the 246's front-panel LP?
(docs/EMT250_program.md section 5, EMT246_unit_facts LP rows, NEXT_STEPS section 7.)

    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/e250/lp_chain.py --stage ratio
    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/e250/lp_chain.py --stage fit

The first tap of each program's direct / cluster path is one single-sample impulse through the analogue output chain, so
(dark first-tap spectrum) / (bright first-tap spectrum), per channel, is the extra filtering of the dark presets.
Within a family (S / MED / LONG / REV-250) the first tap is the same in every preset (agreement -60...-70 dB, the capture
noise), so a family mean is used; the 48 kHz deconvolved Altiverb cache (not committable) is read, derived numbers go to
analysis/e250/lp_chain.json.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CACHE = os.path.join(ROOT, "analysis", "altiverb", "cache")
OUT = os.path.join(ROOT, "analysis", "e250")
FS = 19.2e6 / 728
FA = 48000.0
PRE = 48                                   # 1 ms before the tap
NMAX = 150000
FAM = {
    "S": ["NON-LIN SHORT RT %ss" % x for x in ["0.2", "0.6", "1.0", "1.6", "2.0", "2.2"]],
    "M": ["MED RT %ss" % x for x in ["0.5", "1.0", "1.4", "1.8", "2.0", "2.6", "3.0", "3.5", "4.0", "5.0"]],
    "L": ["LONG RT %ss" % x for x in ["2.0", "3.0", "4.0", "5.0", "6.0", "8.0", "10.0"]],
    "E": ["REV-250 %ss" % x for x in ["0.4", "1.0", "1.5", "2.0", "2.5", "3.0", "3.6", "4.0", "4.5"]],
}
# next tap after the first one (48 kHz samples after the first tap), from the cluster / E taps (checked by --stage taps)
NEXT = {("S", 0): 98, ("S", 1): 260, ("M", 0): 197, ("M", 1): 518, ("L", 0): 393, ("L", 1): 1037, ("E", 0): 781, ("E", 1): 781}
NFFT = 1 << 15
FREQS = np.geomspace(100, 20000, 400)


def load_tap(name, ch, anchor=None, post=480):
    h = np.load(os.path.join(CACHE, name + ".npy"), mmap_mode="r")[ch, :NMAX]
    if anchor is None:
        a = np.abs(h); i = int(np.argmax(a > a.max() * 0.1)); anchor = i - 3 + int(np.argmax(a[i - 3:i + 4]))
    return np.array(h[anchor - PRE: anchor + post]), anchor


def family_tap(fam, ch, post=480):
    """Mean of the family's first taps (each scaled to the first member by least squares over 90 samples), same anchor."""
    _, a0 = load_tap(FAM[fam][0], ch, None, post)
    anchors = []
    for n in FAM[fam]:
        _, a = load_tap(n, ch, None, post); anchors.append(a)
    anc = int(np.bincount(anchors).argmax())           # the mode (some files pick the neighbouring sample of a tie)
    xs = [load_tap(n, ch, anc, post)[0] for n in FAM[fam]]
    ref = xs[0]; acc = np.zeros_like(ref); gains = []
    for x in xs:
        g = np.dot(x[:PRE + 90], ref[:PRE + 90]) / np.dot(ref[:PRE + 90], ref[:PRE + 90])
        acc += x / g; gains.append(g)
    return acc / len(xs), anc, xs, np.array(gains)


def window(x, w_ms, taper_ms=0.25):
    """Rectangular window from -1 ms to +w_ms around the tap (index PRE), raised-cosine edges of taper_ms."""
    n = PRE + int(round(w_ms * FA / 1000)) + 1
    w = np.ones(n)
    nt = max(2, int(round(taper_ms * FA / 1000)))
    ramp = 0.5 * (1 - np.cos(np.pi * (np.arange(nt) + 0.5) / nt))
    w[:nt] = ramp; w[-nt:] = ramp[::-1]
    y = np.zeros(NFFT); y[:n] = x[:n] * w
    return y


def spec(x, w_ms, f=FREQS):
    X = np.fft.rfft(window(x, w_ms)); fr = np.fft.rfftfreq(NFFT, 1 / FA)
    k = np.arange(NFFT // 2 + 1)
    return np.interp(f, fr, X.real) + 1j * np.interp(f, fr, X.imag)


def tap_time(x, w_ms=3.0):
    """Sub-sample tap time (48 kHz samples from the window start) from the phase slope over 400-3000 Hz."""
    X = np.fft.rfft(window(x, w_ms)); fr = np.fft.rfftfreq(NFFT, 1 / FA); m = (fr > 400) & (fr < 3000)
    return -np.polyfit(2 * np.pi * fr[m], np.unwrap(np.angle(X[m])), 1)[0] * FA


def stage_taps():
    for fam in FAM:
        for ch in (0, 1):
            x, anc, xs, g = family_tap(fam, ch)
            d = [20 * np.log10(np.linalg.norm(a[:PRE + 90] / gg - x[:PRE + 90]) / np.linalg.norm(x[:PRE + 90]) + 1e-30)
                 for a, gg in zip(xs, g)]
            a = np.abs(x); 
            print(fam, "LR"[ch], "anchor", anc, "members vs mean (dB) max", round(max(d), 1), "min", round(min(d), 1),
                  "t =", round(anc - PRE + tap_time(x), 3), "(= device sample %.3f)" % ((anc - PRE + tap_time(x)) / (FA / FS * (1 + 42.8e-6))),
                  "gain spread dB", round(20 * np.log10(np.max(np.abs(g)) / np.min(np.abs(g))), 2))


def stage_phase():
    for ch in (0, 1):
        for dark, bright, w in (("E", "M", 3.0), ("L", "M", 3.0), ("E", "L", 3.0)):
            R, lev, ph = ratio(dark, bright, ch, w)
            if dark == "E" and ch == 1:
                R = -R
            phs = np.unwrap(np.angle(R[(FREQS > 200)]))
            fr = FREQS[FREQS > 200]
            gd = -np.gradient(phs, 2 * np.pi * fr) * 1e6             # microseconds
            print(f"{'LR'[ch]} {dark}/{bright}")
            for f in PR:
                i = int(np.argmin(np.abs(fr - f)))
                print(f"   {f:6d} Hz  mag {20*np.log10(abs(R[FREQS>200][i])):6.2f} dB  phase {np.degrees(phs[i]):8.2f} deg  GD {gd[i]:7.2f} us ({gd[i]*1e-6*FS:5.2f} dev samples)")


PR = [250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000]


def fam_spec(fam, ch, w_ms):
    x, anc, _, _ = family_tap(fam, ch)
    return spec(x, w_ms), anc


KDEV = FA / FS * (1 + 42.8e-6)             # 48 kHz samples per device sample, capture clock included
# nominal first-tap time in integer device samples (the DELAY: 10 / 20 / 40 ms plus the constant chain; L = R's, E's L 314 later)
NDEV = {("S", 0): 268, ("S", 1): 268, ("M", 0): 536, ("M", 1): 536, ("L", 0): 1066, ("L", 1): 1066, ("E", 0): 850, ("E", 1): 536}


def abs_phase(S, anc, fam, ch, f=FREQS):
    """Refer the window spectrum to the tap's nominal time: absolute time of the window start, minus NDEV device samples."""
    return S * np.exp(-2j * np.pi * f * ((anc - PRE) - NDEV[(fam, ch)] * KDEV) / FA)


def ratio(dark, bright, ch, w_ms):
    """Complex ratio dark / bright on the absolute time axis, normalised to unit magnitude at 200-1000 Hz; returns
    (R, level_dB, phase at 200-1000 Hz in degrees: 0 = same polarity, 180 = inverted, other = a time offset)."""
    Sd, ad = fam_spec(dark, ch, w_ms); Sb, ab = fam_spec(bright, ch, w_ms)
    R = abs_phase(Sd, ad, dark, ch) / abs_phase(Sb, ab, bright, ch)
    lf = (FREQS > 200) & (FREQS < 1000)
    g = np.mean(np.abs(R[lf]))
    ph = np.degrees(np.angle(np.mean(R[lf] / np.abs(R[lf]))))
    return R / g, 20 * np.log10(g), ph


def stage_ratio():
    for ch in (0, 1):
        for dark in ("E", "L"):
            for bright in ("S", "M"):
                for w in (1.0, 1.7, 3.0, 4.0, 8.0):
                    if w > min(NEXT[(dark, ch)], NEXT[(bright, ch)]) / FA * 1000 - 0.15:
                        continue
                    R, lev, ph = ratio(dark, bright, ch, w)
                    row = []
                    for f in PR:
                        i = int(np.argmin(np.abs(FREQS - f)))
                        row.append(f"{20*np.log10(abs(R[i])):6.2f}")
                    print(f"{'LR'[ch]} {dark}/{bright} W={w:3.1f}ms LF phase {ph:+7.1f} deg | " + " ".join(row))
    print("freqs", PR)


# ------------------------------------------------------------------ fits
from scipy import signal
from scipy.optimize import least_squares

FMAX_FIT = 12000.0
FIT = FREQS[(FREQS >= 100) & (FREQS <= FMAX_FIT)]
IFIT = np.nonzero((FREQS >= 100) & (FREQS <= FMAX_FIT))[0]


def get_ratio(dark, bright, ch, w):
    R, lev, ph = ratio(dark, bright, ch, w)
    if abs(ph) > 90:                           # E's R is inverted relative to the others: polarity is not part of the filter
        R = -R
    return R[IFIT]


def noise_floor(fam, ch, w):
    """Ratio of two disjoint halves of one family's members (should be 1): the estimate noise, dB and degrees, per band."""
    _, anc, xs, g = family_tap(fam, ch)
    xa = np.mean([x / gg for x, gg in list(zip(xs, g))[0::2]], axis=0); xb = np.mean([x / gg for x, gg in list(zip(xs, g))[1::2]], axis=0)
    r = spec(xa, w)[IFIT] / spec(xb, w)[IFIT]
    return 20 * np.log10(np.abs(r)), np.degrees(np.angle(r))


def model_resp(kind, p, f):
    """Complex response of the candidate filter; p = parameter vector. Returns H(f) (delay included where the form has it)."""
    if kind.startswith("fir"):
        k = int(kind[3:]); a = np.asarray(p[:k - 1]); h = np.concatenate([[1 - a.sum()], a])
        return np.sum(h[:, None] * np.exp(-2j * np.pi * f[None, :] * np.arange(k)[:, None] / FS), axis=0)
    typ, n = kind.split("_")[0], int(kind.split("_")[1]); fc = p[0]; tau = p[-1]
    if typ in ("butter", "bessel", "cheby1", "cheby2", "ellip"):
        if typ == "butter": z, pl, k = signal.buttap(n)
        elif typ == "bessel": z, pl, k = signal.besselap(n, norm="mag")
        elif typ == "cheby1": z, pl, k = signal.cheb1ap(n, p[1])
        elif typ == "cheby2": z, pl, k = signal.cheb2ap(n, p[1])
        else: z, pl, k = signal.ellipap(n, p[1], p[2])
        w, h = signal.freqs_zpk(z, pl, k, f / fc)
    elif typ == "dbutter":                               # bilinear Butterworth at the device rate
        b, a_ = signal.butter(n, fc, fs=FS); w, h = signal.freqz(b, a_, worN=f, fs=FS)
    return h * np.exp(-2j * np.pi * f * tau)


# (kind, initial p, lower, upper); fc in Hz, tau in seconds (free; the 'notau' variants fix it at zero)
def candidates():
    c = []
    for n in range(1, 7):
        c.append((f"butter_{n}", [9e3 * (1 + 0.3 * n), 0], [1e3, -60e-6], [8e4, 60e-6]))
        c.append((f"bessel_{n}", [9e3, 0], [1e3, -60e-6], [2e5, 60e-6]))
        c.append((f"dbutter_{n}", [9e3 * (1 + 0.3 * n), 0], [1e3, -60e-6], [1.31e4, 60e-6]))
        c.append((f"cheby1_{n}", [1.2e4, 0.5, 0], [1e3, 0.01, -60e-6], [8e4, 3.0, 60e-6]))
        c.append((f"cheby2_{n}", [1.2e4, 30, 0], [1e3, 8, -60e-6], [8e4, 90, 60e-6]))
        if n >= 2: c.append((f"ellip_{n}", [1.2e4, 0.5, 30, 0], [1e3, 0.01, 8, -60e-6], [8e4, 3.0, 90, 60e-6]))
    for k in (2, 3, 4):
        c.append((f"fir{k}", [0.1] * (k - 1), [-1.0] * (k - 1), [1.0] * (k - 1)))
    return c


def resid(kind, p, f, R):
    e = np.log(model_resp(kind, p, f) / R)
    e = e.real + 1j * np.angle(np.exp(1j * e.imag))
    return np.concatenate([e.real, e.imag])


def fit_one(kind, p0, lo, hi, f, R, nstart=6, fix_tau=False):
    rng = np.random.default_rng(1)
    best = None
    lo = np.array(lo, float); hi = np.array(hi, float)
    for i in range(nstart):
        x0 = np.array(p0, float) if i == 0 else lo + (hi - lo) * rng.random(len(lo))
        if kind.startswith(("butter", "bessel", "cheby", "ellip", "dbutter")) and i:
            x0[0] = np.exp(rng.uniform(np.log(4e3), np.log(5e4 if not kind.startswith("dbutter") else 1.3e4)))
        if fix_tau and not kind.startswith("fir"):
            x0[-1] = 0.0
        l2, h2 = lo.copy(), hi.copy()
        if fix_tau and not kind.startswith("fir"):
            l2[-1], h2[-1] = -1e-12, 1e-12
        try:
            r = least_squares(lambda q: resid(kind, q, f, R), np.clip(x0, l2, h2), bounds=(l2, h2), x_scale=np.maximum(np.abs(x0), 1e-3))
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best


def score(kind, p, f, R):
    H = model_resp(kind, p, f); e = np.log(H / R)
    db = 8.685889638 * e.real; ph = np.degrees(np.angle(np.exp(1j * e.imag)))
    return float(np.sqrt(np.mean(db ** 2))), float(np.max(np.abs(db))), float(np.sqrt(np.mean(ph ** 2))), float(np.max(np.abs(ph)))


def stage_abs():
    """Absolute first-tap magnitude response of each family, dB re 1 kHz (window 1.7 ms, channels L / R)."""
    fs_ = [500, 1000, 2000, 4000, 6000, 8000, 9000, 10000, 11000, 11500, 12000, 12100, 12200, 12400, 12800, 13500, 14500, 16000, 18000]
    print("freq (Hz):", fs_)
    for fam in ("S", "M", "L", "E"):
        for ch in (0, 1):
            x, anc, _, _ = family_tap(fam, ch)
            fr = np.fft.rfftfreq(NFFT, 1 / FA); X = np.abs(np.fft.rfft(window(x, 1.7)))
            ref = X[np.argmin(np.abs(fr - 1000))]
            print(f"{fam} {'LR'[ch]}", " ".join(f"{20*np.log10(X[np.argmin(np.abs(fr-f))]/ref):6.2f}" for f in fs_))


def stage_fit(args):
    sets = {}
    for ch in (0, 1):
        for dark in ("E", "L"):
            sets[(dark, "M", ch)] = get_ratio(dark, "M", ch, 3.0)
            sets[(dark, "S", ch)] = get_ratio(dark, "S", ch, 1.7)
        sets[("M", "S", ch)] = get_ratio("M", "S", ch, 1.7)
        sets[("E", "L", ch)] = get_ratio("E", "L", ch, 3.0)
    nf = {}
    for fam in ("E", "M", "L"):
        db, ph = noise_floor(fam, 1, 3.0)
        nf[fam] = (float(np.sqrt(np.mean(db ** 2))), float(np.sqrt(np.mean(ph ** 2))))
    print("ratio noise floor from split halves (rms dB, rms deg) E/M/L:", {k: tuple(round(x, 3) for x in v) for k, v in nf.items()})
    train = ("E", "M", 1)                               # fit: REV-250 R against MED R (W = 3 ms)
    held_M = [k for k in sets if k[1] == "M" and k != train]
    held_S = [k for k in sets if k[1] == "S" and k[0] in ("E", "L")]
    Rtr = sets[train]
    print("fit on", train, "100 Hz-12 kHz; held out: the other 3 dark/MED ratios (L ch, LONG) and the 4 dark/S ratios (different bright reference)")
    print(f"{'form':14s} {'np':>2s} | fit dB/deg  | held vs M (mean dB, max dB, mean deg) | held vs S dB/deg | params")
    res = {}
    for fix_tau in (False, True):
        for kind, p0, lo, hi in candidates():
            if fix_tau and kind.startswith("fir"):
                continue
            r = fit_one(kind, p0, lo, hi, FIT, Rtr, fix_tau=fix_tau)
            if r is None: continue
            p = r.x; tr = score(kind, p, FIT, Rtr)
            hm = [score(kind, p, FIT, sets[k]) for k in held_M]; hs = [score(kind, p, FIT, sets[k]) for k in held_S]
            name = kind + ("/tau0" if fix_tau else "")
            res[name] = {"params": list(map(float, p)), "train_rms_db_deg": [tr[0], tr[2]],
                         "held_vs_M": [float(np.mean([h[0] for h in hm])), float(np.max([h[1] for h in hm])), float(np.mean([h[2] for h in hm]))],
                         "held_vs_S": [float(np.mean([h[0] for h in hs])), float(np.mean([h[2] for h in hs]))]}
            q = res[name]
            print(f"{name:14s} {len(p) - (1 if fix_tau and not kind.startswith('fir') else 0):2d} | {tr[0]:6.3f} {tr[2]:5.2f} | {q['held_vs_M'][0]:7.3f} {q['held_vs_M'][1]:7.3f} {q['held_vs_M'][2]:6.2f} | {q['held_vs_S'][0]:6.3f} {q['held_vs_S'][1]:5.2f} | "
                  + " ".join(f"{x:.4g}" for x in p))
    # the one-parameter two-tap form refitted on every pair: a = fractional device-sample delay with linear interpolation
    print("two-tap [1-a, a] refit per pair (a; resulting rms dB, deg); LF delay = a / FS:")
    pa = {}
    for k, R in sets.items():
        r = fit_one("fir2", [0.1], [-1.0], [1.0], FIT, R)
        sc = score("fir2", r.x, FIT, R)
        pa["%s/%s %s" % (k[0], k[1], "LR"[k[2]])] = [float(r.x[0]), sc[0], sc[2]]
        print(f"   {k[0]}/{k[1]} {'LR'[k[2]]}: a = {r.x[0]:+.4f} ({r.x[0] / FS * 1e6:+.2f} us)  rms {sc[0]:.3f} dB {sc[2]:.2f} deg")
    res["fir2_per_pair"] = pa; res["noise_floor_rms_db_deg"] = nf
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, "lp_chain_fits.json"), "w"), indent=1)


def stage_shape(args):
    """What sets the scale of the droop? Two-tap form with a free tap spacing T' (a, T'), and an extra zero-order hold."""
    R = get_ratio("E", "M", 1, 3.0); Rl = get_ratio("E", "M", 0, 3.0); Rg = get_ratio("L", "M", 1, 3.0)
    def h2(q, f):
        return (1 - q[0]) + q[0] * np.exp(-2j * np.pi * f * q[1] * 1e-6)
    def hold(q, f):                                   # ratio of a hold of length T + h to a hold of length T (h in device periods)
        T = 1 / FS; hh = q[0] * T
        return np.sinc(f * (T + hh)) / np.sinc(f * T) * np.exp(-1j * np.pi * f * hh)
    for name, fn, p0, lo, hi in (("two-tap, free spacing (a, T' us)", h2, [0.2, 38.0], [-1, 5.0], [1, 200.0]),
                                  ("extra zero-order hold h/T", hold, [0.4], [0.0], [3.0])):
        best = None
        for s0 in ([0.2, 38.0], [0.5, 25.0], [0.3, 60.0], [0.1, 120.0]) if len(p0) == 2 else (p0, [1.0]):
            r = least_squares(lambda q: np.concatenate([(lambda e: [e.real, np.angle(np.exp(1j * e.imag))])(np.log(fn(q, FIT) / R))[0],
                                                       (lambda e: [e.real, np.angle(np.exp(1j * e.imag))])(np.log(fn(q, FIT) / R))[1]]), s0, bounds=(lo, hi))
            if best is None or r.cost < best.cost: best = r
        e = np.log(fn(best.x, FIT) / R)
        print(f"{name}: params {np.round(best.x, 4)}  rms {8.686 * np.sqrt(np.mean(e.real ** 2)):.3f} dB {np.degrees(np.sqrt(np.mean(np.angle(np.exp(1j * e.imag)) ** 2))):.2f} deg"
              f"  (device period = {1e6 / FS:.2f} us; 44.1 kHz = 22.68, 48 kHz = 20.83, 96 kHz = 10.42)")
        # scan T' to show how sharply it is fixed
        if len(p0) == 2:
            for Tp in (20.83, 22.68, 30.0, 34.0, 37.92, 42.0, 50.0, 76.0):
                rr = least_squares(lambda a_: np.concatenate([np.log(np.abs(h2([a_[0], Tp], FIT) / R)), np.angle(h2([a_[0], Tp], FIT) / R)]), [0.2], bounds=([-1], [1]))
                e = np.log(h2([rr.x[0], Tp], FIT) / R)
                print(f"      T' = {Tp:6.2f} us: a = {rr.x[0]:+.3f}  rms {8.686 * np.sqrt(np.mean(e.real ** 2)):.3f} dB {np.degrees(np.sqrt(np.mean(np.angle(np.exp(1j * e.imag)) ** 2))):.2f} deg")


def stage_kernel():
    """Compare with the Tomkin 96 kHz M / L output kernel (ml_structure 7.5; M [.0235 .9518 .0276], L [.1134 .7834 .1217], symmetric,
    S = delta): fit symmetric [p/2, 1-p, p/2] (+ a free delay), the one-sided [1-a, a], and a general taps -1..+2 form to
    M/S, E/S, LONG/S per channel (the bright reference S is the Tomkin S's twin)."""
    T = 1 / FS
    def sym(q, f, free_tau):
        h = q[0] / 2 * np.exp(2j * np.pi * f * T) + (1 - q[0]) + q[0] / 2 * np.exp(-2j * np.pi * f * T)
        return h * (np.exp(-2j * np.pi * f * q[1] * 1e-6) if free_tau else 1)
    def gen(q, f):                                   # taps at -1, 0, +1, +2 (sum 1)
        t = np.array([q[0], 1 - q[0] - q[1] - q[2], q[1], q[2]]); k = np.array([-1, 0, 1, 2])
        return np.sum(t[:, None] * np.exp(-2j * np.pi * f[None, :] * k[:, None] * T), axis=0)
    def fit(fn, p0, R):
        r = least_squares(lambda q: np.concatenate([(lambda e: np.concatenate([e.real, np.angle(np.exp(1j * e.imag))]))(np.log(fn(q, FIT) / R))]), p0)
        e = np.log(fn(r.x, FIT) / R)
        return r.x, 8.686 * np.sqrt(np.mean(e.real ** 2)), np.degrees(np.sqrt(np.mean(np.angle(np.exp(1j * e.imag)) ** 2)))
    out = {}
    print("pair   ch | one-sided a (dB,deg) | symmetric p, zero delay (dB,deg) | symmetric p + free delay us (dB,deg) | general taps -1..+2")
    for (d, b) in (("M", "S"), ("E", "S"), ("L", "S")):
        for ch in (0, 1):
            R = get_ratio(d, b, ch, 1.7)
            x1, m1, p1 = fit(lambda q, f: (1 - q[0]) + q[0] * np.exp(-2j * np.pi * f * T), [0.1], R)
            x2, m2, p2 = fit(lambda q, f: sym(q, f, False), [0.1, 0], R)
            x3, m3, p3 = fit(lambda q, f: sym(q, f, True), [0.1, 0.0], R)
            x4, m4, p4 = fit(gen, [0.05, 0.2, 0.0], R)
            out[f"{d}/{b} {'LR'[ch]}"] = dict(onesided_a=float(x1[0]), sym_p=float(x2[0]), sym_free_p=float(x3[0]), sym_free_tau_us=float(x3[1]), general=list(map(float, x4)))
            print(f"{d}/{b} {'LR'[ch]} | a={x1[0]:.3f} ({m1:.3f},{p1:.2f}) | p={x2[0]:.3f} ({m2:.3f},{p2:.2f}) | p={x3[0]:.3f} tau={x3[1]:+.1f} ({m3:.3f},{p3:.2f}) | "
                  f"taps -1..2 = {x4[0]:+.3f} {1 - x4.sum():+.3f} {x4[1]:+.3f} {x4[2]:+.3f} ({m4:.3f},{p4:.2f})")
    # Tomkin kernels
    for nm, k in (("Tomkin M", [-0.0093, 0.0235, 0.9518, 0.0276, 0.0063]), ("Tomkin L", [-0.0392, 0.1134, 0.7834, 0.1217, 0.0173])):
        H = np.sum([c * np.exp(-2j * np.pi * PR_ * (i - 2) * T) for i, c in enumerate(k)], axis=0)
        print(nm, "mag dB at", list(PR_), np.round(20 * np.log10(np.abs(H)), 2))
    for (d, ch) in (("M", 1), ("L", 1), ("E", 1)):
        R = ratio(d, "S", ch, 1.7)[0]
        print(f"Altiverb {d}/S R mag dB:", [round(20 * np.log10(abs(R[np.argmin(np.abs(FREQS - f))])), 2) for f in PR_])
    json.dump(out, open(os.path.join(OUT, "lp_chain_kernel.json"), "w"), indent=1)


PR_ = np.array([2000, 4000, 6000, 8000, 10000, 12000.0])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", default="taps"); args = ap.parse_args()
    {"taps": lambda: stage_taps(), "ratio": lambda: stage_ratio(), "phase": lambda: stage_phase(), "fit": lambda: stage_fit(args), "abs": lambda: stage_abs(), "shape": lambda: stage_shape(args), "kernel": lambda: stage_kernel()}[args.stage]()
