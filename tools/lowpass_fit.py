"""The 246's LOW PASS lever's filter shape, inferred from the Altiverb presets read below 12 kHz (presets.json).

    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/lowpass_fit.py

Per preset (S / M / L / E / N with a first tap): the first tap (cluster tap 0; E's direct tap 0; N's front tap 41),
both channels, windowed (sqrt-Hann, +-HALF of tools/altiverb_presets.py), divided by its own family's first tap in the
same channel (the family-relative transfer: the Altiverb's capture chain and the program's output kernel cancel).
Check of that cancellation: each family's first tap over the S family's must be the program's output kernel
(EMT250_program.md 14.1: L / E [0.77, 0.23], M 0.047).

Measured per 1/6 octave 1-13 kHz (magnitude dB, phase deg after the fitted delay). Fitted 1-12 kHz on the complex
log ratio (1 dB = 0.115 Np against 1 rad; equal weight per log frequency), each form with a free gain and a free
delay (nuisance, not counted), at the device rate 26,373.6 Hz:
  fir2       [1 - a, a]                                   1
  repl       fir2(a_tot) / fir2(a_family): the program's output kernel replaced, not cascaded      1
  pole1      (1 - p) / (1 - p z^-1)                       1
  bw2        the plugin's DeviceLowPass (2nd-order Butterworth, bilinear, prewarped) at fc       1
  biquad     general biquad, DC gain 1                    4 (+ the gain = 5 coefficients)
  fir2pole   fir2(a) * pole1(p)                           2
  bwN        N-th order Butterworth at fc (PRESET23 only) 2
Also per preset: the noise floor at the file's end over the family's (uninformative: the files end in digital silence).

Stage `extra` (analysis/lowpass/extra.json, extra.txt): the free tap spacing of the 2-tap replacement (device period or
not); PRESET23 in detail (linear-phase test, notches, Chebyshev / elliptic IIR fits, symmetric FIRs at the device rate
and 48 kHz); the candidate law (PRESET23's device-rate 13-tap shape scaled in frequency to each panel step, H) against
the plugin's Butterworth per step (-3 dB point, level at 5 / 8 / 10 / 12 kHz).

    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/lowpass_fit.py [--stage fit|extra|all]
Output: analysis/lowpass/fit.json, report.txt (stage fit); extra.json, extra.txt (stage extra).
"""
import json, os, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import altiverb_presets as ap  # noqa: E402

FS, FSD = ap.FS, ap.FSD
OUT = os.path.join(ROOT, "analysis", "lowpass")
NFFT = 1 << 15
F = np.fft.rfftfreq(NFFT, 1 / FS)
FIT = (F >= 1000) & (F <= 12000)
W = 1 / F[FIT]
W = W / W.sum()
BANDS = [1000 * 2 ** (i / 6) for i in range(23)]          # 1.0 ... 12.7 kHz
FAM_A = {"S": 0.0, "M": 0.047, "L": 0.23, "E": 0.23, "N": 0.0}


def zinv(f):
    return np.exp(-2j * np.pi * f / FSD)


def h_fir2(f, a):
    return (1 - a) + a * zinv(f)


def h_pole1(f, p):
    return (1 - p) / (1 - p * zinv(f))


def h_bw(f, fc, order=2):
    """Bilinear N-th order Butterworth at the device rate, prewarped (order 2 = DeviceLowPass)."""
    k = np.tan(np.pi * fc / FSD)
    s = 1j * np.tan(np.pi * f / FSD) / k          # bilinear frequency map
    h = np.ones_like(s)
    for m in range(order // 2):
        th = np.pi * (2 * m + 1) / (2 * order)
        h = h / (s * s + 2 * np.sin(th) * s + 1)
    if order % 2:
        h = h / (s + 1)
    return h * (1 + 0 * f)


def h_biquad(f, q):
    b1, b2, a1, a2 = q
    z = zinv(f)
    g = (1 + a1 + a2) / (1 + b1 + b2)
    return g * (1 + b1 * z + b2 * z * z) / (1 + a1 * z + a2 * z * z)


# ---------------------------------------------------------------- measurement
def first_tap(name, prog, onset_ms):
    """(index in 48 kHz samples per channel, sign) of the first tap."""
    ir = ap.load(name)
    k = int(round(onset_ms * 1e-3 * FS + (41 / FSD * FS if prog == "N" else 0)))
    ks = []
    for c in (0, 1):
        kc = k + (int(round(314 / FSD * FS)) if (prog == "E" and c == 0) else 0)
        kc = kc + int(np.argmax(np.abs(np.array(ir[c][kc - 4: kc + 5])))) - 4
        ks.append(kc)
    return ir, ks


def spectrum(ir, c, k, half):
    n = int(half * FS)
    seg = np.array(ir[c][k - n: k + n + 1]) * np.hanning(2 * n + 1) ** 0.5
    x = np.zeros(NFFT)
    x[: n + 1] = seg[n:]
    x[-n:] = seg[:n]                                  # time origin at the window centre
    return np.fft.rfft(x)


def noise_spec(ir, c):
    """Mean power spectrum of the file's last 2 s (Welch, 4096)."""
    x = np.array(ir[c][-2 * FS:])
    segs = [x[i: i + 4096] * np.hanning(4096) for i in range(0, len(x) - 4096, 2048)]
    return np.mean([np.abs(np.fft.rfft(s)) ** 2 for s in segs], axis=0), np.fft.rfftfreq(4096, 1 / FS)


def band_vals(h, tau):
    """Per 1/6 octave: magnitude dB (power mean) and phase deg (after removing the delay tau, seconds)."""
    ph = np.unwrap(np.angle(h * np.exp(2j * np.pi * F * tau)))
    ph -= 2 * np.pi * np.round(np.interp(1000, F, ph) / (2 * np.pi))
    mag, pha = [], []
    for fc in BANDS:
        m = (F >= fc * 2 ** (-1 / 12)) & (F < fc * 2 ** (1 / 12))
        mag.append(float(10 * np.log10(np.mean(np.abs(h[m]) ** 2))))
        pha.append(float(np.degrees(np.mean(ph[m]))))
    return mag, pha


# ---------------------------------------------------------------- fitting
def fit_form(h, form, fam_a, x0=None, order=None):
    f = F[FIT]
    lh = np.log(h[FIT])
    lh = np.real(lh) + 1j * np.unwrap(np.imag(lh))
    sw = np.sqrt(W)

    def model(q):
        if form == "fir2":
            return h_fir2(f, q[0])
        if form == "repl":
            return h_fir2(f, q[0]) / h_fir2(f, fam_a)
        if form == "pole1":
            return h_pole1(f, q[0])
        if form == "bw2":
            return h_bw(f, q[0] * 1e3, 2)
        if form == "bwN":
            return h_bw(f, q[0] * 1e3, order)
        if form == "biquad":
            return h_biquad(f, q)
        if form == "fir2pole":
            return h_fir2(f, q[0]) * h_pole1(f, q[1])
        if form == "open":
            return np.ones_like(f, dtype=complex)

    def resid(qq):
        q, lg, tau = qq[:-2], qq[-2], qq[-1]
        hm = model(q)
        d = lh - (np.log(hm) + lg - 2j * np.pi * f * tau * 1e-6)
        dr, di = np.real(d), np.imag(d)
        di = di - 2 * np.pi * np.round(np.sum(W * di) / (2 * np.pi))
        return np.concatenate([sw * dr, sw * di])

    starts = {"fir2": [[0.1], [0.3]], "repl": [[0.3], [0.45]], "pole1": [[0.1], [0.4]],
              "bw2": [[10.0], [7.0], [11.5]], "bwN": [[9.0], [10.0]],
              "biquad": [[2.0, 1.0, -0.5, 0.2], [1.0, 0.0, 0.0, 0.0], [1.6, 1.0, -0.3, 0.3]],
              "fir2pole": [[0.1, 0.1], [0.3, 0.2]], "open": [[]]}[form]
    lb = {"fir2": [-0.5], "repl": [0.0], "pole1": [-0.9], "bw2": [1.0], "bwN": [1.0], "biquad": [-4, -4, -1.99, -0.99],
          "fir2pole": [-0.5, -0.9], "open": []}[form]
    ub = {"fir2": [0.5], "repl": [0.5], "pole1": [0.95], "bw2": [13.1], "bwN": [13.1], "biquad": [4, 4, 1.99, 0.99],
          "fir2pole": [0.5, 0.95], "open": []}[form]
    # the delay: start from the phase slope at 1-4 kHz
    m = (f < 4000)
    sl = np.polyfit(f[m], np.unwrap(np.imag(lh))[m], 1)[0]
    tau0 = -sl / (2 * np.pi) * 1e6
    best = None
    for s in starts:
        try:
            r = least_squares(resid, list(s) + [0.0, tau0], bounds=(lb + [-np.inf, -1e4], ub + [np.inf, 1e4]))
        except ValueError:
            continue
        if best is None or r.cost < best.cost:
            best = r
    q, lg, tau = best.x[:-2], best.x[-2], best.x[-1]
    hm = model(q)
    d = lh - (np.log(hm) + lg - 2j * np.pi * f * tau * 1e-6)
    dr = np.real(d) * 20 / np.log(10)
    di = np.degrees(np.imag(d))
    di -= 360 * np.round(np.sum(W * di) / 360)
    return {"params": [float(v) for v in q], "delay_us": float(tau), "gain_db": float(lg * 20 / np.log(10)),
            "rms_db": float(np.sqrt(np.sum(W * dr ** 2))), "rms_deg": float(np.sqrt(np.sum(W * di ** 2))),
            "max_db": float(np.abs(dr).max())}


def model_db(form, q, f, fam_a=0.0, order=2):
    if form == "fir2":
        h = h_fir2(f, q[0])
    elif form == "repl":
        h = h_fir2(f, q[0]) / h_fir2(f, fam_a)
    elif form == "pole1":
        h = h_pole1(f, q[0])
    elif form == "bw2":
        h = h_bw(f, q[0] * 1e3, 2)
    elif form == "bwN":
        h = h_bw(f, q[0] * 1e3, order)
    elif form == "biquad":
        h = h_biquad(f, q)
    else:
        h = h_fir2(f, q[0]) * h_pole1(f, q[1])
    return 20 * np.log10(np.abs(h))


def stage_fit():
    os.makedirs(OUT, exist_ok=True)
    pj = json.load(open(os.path.join(ROOT, "analysis", "altiverb", "presets.json")))["presets"]
    names = [n for n, v in pj.items() if v.get("program") in ("S", "M", "L", "E", "N")]
    spec, info = {}, {}
    for n in names:
        v = pj[n]
        prog = v["program"]
        if prog == "L" and n == "PRESET90 VERY LG":
            continue                       # no cluster
        ir, ks = first_tap(n, prog, v["onset_ms"])
        spec[n] = [spectrum(ir, c, ks[c], ap.HALF[prog]) for c in (0, 1)]
        ns, fn = noise_spec(ir, 0)
        info[n] = {"prog": prog, "lp_set": v["settings"]["lowpass"]["value"], "hf_set": v["settings"].get("hf", {}).get("value"),
                   "noise": ns, "fn": fn, "peak": [float(np.abs(np.array(ir[c][ks[c] - 2: ks[c] + 3])).max()) for c in (0, 1)]}
        del ir
    res = {"method": __doc__, "bands_hz": BANDS, "kernel_check": {}, "presets": {}}
    lines = []
    P = lines.append

    # ---- 1. cancellation check: each family's first tap over S's (needs the same window: S's HALF is short, use L's)
    P("== Kernel check: family first tap / S family first tap (window +-HALF of the numerator's program) ==")
    for prog in "MLEN":
        fam = ap.FAMILY_REF[prog]
        out = {}
        for c in (0, 1):
            ir_f, ks_f = first_tap(fam, prog, pj[fam]["onset_ms"])
            ir_s, ks_s = first_tap(ap.FAMILY_REF["S"], "S", pj[ap.FAMILY_REF["S"]]["onset_ms"])
            hw = ap.HALF["S"]
            h = spectrum(ir_f, c, ks_f[c], hw) / spectrum(ir_s, c, ks_s[c], hw)
            r = {k: fit_form(h, k, 0.0) for k in ("fir2", "pole1", "bw2")}
            out["LR"[c]] = r
            mag, _ = band_vals(h, r["fir2"]["delay_us"] * 1e-6)
            P(f"  {prog} {'LR'[c]}: fir2 a = {r['fir2']['params'][0]:.3f} ({r['fir2']['rms_db']:.3f} dB / {r['fir2']['rms_deg']:.2f} deg)  "
              f"pole1 p = {r['pole1']['params'][0]:.3f} ({r['pole1']['rms_db']:.3f} dB / {r['pole1']['rms_deg']:.2f} deg)  "
              f"bw2 {r['bw2']['params'][0]:.2f} kHz ({r['bw2']['rms_db']:.3f} dB / {r['bw2']['rms_deg']:.2f} deg)  "
              f"| 8/10/12 kHz {np.interp(8000, BANDS, mag) - mag[0]:+.2f} {np.interp(10000, BANDS, mag) - mag[0]:+.2f} {np.interp(12000, BANDS, mag) - mag[0]:+.2f} dB")
        res["kernel_check"][prog] = out

    # ---- 2. every preset over its family
    P("")
    P("== Presets / family: shape params; rms dB / rms deg (1-12 kHz, gain + delay free); per channel L | R ==")
    forms = ["fir2", "repl", "pole1", "bw2", "biquad", "fir2pole"]
    for n in names:
        if n not in spec:
            continue
        prog = info[n]["prog"]
        fam = ap.FAMILY_REF[prog]
        if n == fam:
            continue
        r = {"prog": prog, "lp_set": info[n]["lp_set"], "hf_set": info[n]["hf_set"], "ch": {}}
        for c in (0, 1):
            h = spec[n][c] / spec[fam][c]
            fr = {k: fit_form(h, k, FAM_A[prog]) for k in forms + ["open"] if not (k == "repl" and FAM_A[prog] == 0)}
            if n.startswith("PRESET23"):
                for o in (2, 4, 6, 8, 10):
                    fr[f"bw{o}"] = fit_form(h, "bwN", 0.0, order=o)
            mag, pha = band_vals(h, fr["fir2"]["delay_us"] * 1e-6)
            mag = [m - mag[0] for m in mag]
            fr["bands_db"], fr["bands_deg"] = mag, pha
            r["ch"]["LR"[c]] = fr
        # the noise floor at the file's end over the family's (L channel), 4-11 kHz relative to 1-3 kHz
        nf = info[n]["noise"] / info[fam]["noise"]
        fn = info[n]["fn"]
        lo = np.mean(nf[(fn > 1000) & (fn < 3000)])
        r["noise_rel_db"] = {k: float(10 * np.log10(np.mean(nf[(fn > k * 1e3 * 0.95) & (fn < k * 1e3 * 1.05)]) / lo)) for k in (5, 8, 9, 10, 11, 12)}
        r["noise_abs_db"] = {k: float(10 * np.log10(np.mean(info[n]["noise"][(fn > k * 1e3 * 0.95) & (fn < k * 1e3 * 1.05)]) / info[n]["peak"][0] ** 2)) for k in (2, 8, 10, 11)}
        res["presets"][n] = r
        fl = r["ch"]["L"]
        P(f"{n:22s} {prog} LPset {r['lp_set']:>2} HF {r['hf_set']}  8/10/11/12k {np.interp(8000, BANDS, fl['bands_db']):+6.2f} "
          f"{np.interp(10000, BANDS, fl['bands_db']):+6.2f} {np.interp(11000, BANDS, fl['bands_db']):+6.2f} {np.interp(12000, BANDS, fl['bands_db']):+6.2f}")
        for k in forms + ["open"]:
            if k not in fl:
                continue
            a, b = r["ch"]["L"][k], r["ch"]["R"][k]
            ps = ",".join(f"{x:.3f}" for x in a["params"]) + " | " + ",".join(f"{x:.3f}" for x in b["params"])
            P(f"    {k:9s} {ps:34s} {a['rms_db']:6.3f} / {a['rms_deg']:5.2f}  | {b['rms_db']:6.3f} / {b['rms_deg']:5.2f}")
        if n.startswith("PRESET23"):
            for o in (2, 4, 6, 8, 10):
                a, b = r["ch"]["L"][f"bw{o}"], r["ch"]["R"][f"bw{o}"]
                P(f"    bw{o:<7d} {a['params'][0]:.3f} | {b['params'][0]:.3f} kHz   {a['rms_db']:6.3f} / {a['rms_deg']:6.2f}  | {b['rms_db']:6.3f} / {b['rms_deg']:6.2f}")
        P(f"    noise floor / family's (re 1-3 kHz) 5/8/9/10/11/12k: " + " ".join(f"{r['noise_rel_db'][k]:+.1f}" for k in (5, 8, 9, 10, 11, 12))
          + "   abs re tap peak 2/8/10/11k: " + " ".join(f"{r['noise_abs_db'][k]:.0f}" for k in (2, 8, 10, 11)))
    json.dump(res, open(os.path.join(OUT, "fit.json"), "w"), indent=1)
    open(os.path.join(OUT, "report.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------------------------------------------------------------- stage extra
LP_STEPS = [2, 3, 4, 5, 6, 8, 10, 12]          # plugin/src/params/Parameters.cpp lowpassStepsKHz


def tap_h(n, fam, prog, c=0, half=None):
    pj = json.load(open(os.path.join(ROOT, "analysis", "altiverb", "presets.json")))["presets"]
    hw = half or ap.HALF[prog]
    a = first_tap(n, prog, pj[n]["onset_ms"])
    b = first_tap(fam, prog, pj[fam]["onset_ms"])
    return spectrum(a[0], c, a[1][c], hw) / spectrum(b[0], c, b[1][c], hw)


def sym_fir(f, t, K, rate):
    M = (K - 1) // 2
    X = np.stack([np.ones_like(f)] + [2 * np.cos(2 * np.pi * f * k / rate) for k in range(1, M + 1)], 1)
    return np.linalg.lstsq(X, t, rcond=None)[0]


def sym_amp(c, f, rate):
    return c[0] + 2 * sum(c[k] * np.cos(2 * np.pi * f * k / rate) for k in range(1, len(c)))


def m3db(fun, lo=500.0, hi=13180.0):
    g = np.linspace(lo, hi, 20000)
    v = fun(g)
    i = np.argmax(v < -3.0103)
    return float(g[i]) if v[i] < -3.0103 else None


def stage_extra():
    from scipy.signal import ellip, cheby2, freqz
    os.makedirs(OUT, exist_ok=True)
    out, L = {}, []
    P = L.append
    # 1. free spacing of the 2-tap replacement
    P("== Free tap spacing of the replacement 2-tap (us; device period 37.92, 48 kHz 20.83) ==")
    f = F[FIT]
    sw = np.sqrt(W)
    out["spacing"] = {}
    for n, fam, prog in (("PRESET9 CONCERT", "LONG RT 6.0s", "L"), ("PRESET4 STRINGS", "LONG RT 6.0s", "L"),
                         ("PRESET2 BRASS", "REV-250 2.0s", "E"), ("PRESET86 MED WARM", "REV-250 2.0s", "E"),
                         ("PRESET82 SM REF", "MED RT 3.0s", "M"), ("PRESET81 SM LIVE", "MED RT 3.0s", "M")):
        h = tap_h(n, fam, prog)
        lh = np.log(h[FIT])
        lh = np.real(lh) + 1j * np.unwrap(np.imag(lh))
        fa = FAM_A[prog]

        def r(q):
            a, T, lg, tau = q
            hm = ((1 - a) + a * np.exp(-2j * np.pi * f * T * 1e-6)) / h_fir2(f, fa)
            d = lh - (np.log(hm) + lg - 2j * np.pi * f * tau * 1e-6)
            di = np.imag(d)
            di = di - 2 * np.pi * np.round(np.sum(W * di) / (2 * np.pi))
            return np.concatenate([sw * np.real(d), sw * di])
        best = min((least_squares(r, [a0, T0, 0, 0]) for a0 in (0.05, 0.3, 0.45) for T0 in (30, 37.9, 45)), key=lambda z: z.cost)
        out["spacing"][n] = {"a": float(best.x[0]), "spacing_us": float(best.x[1])}
        P(f"  {n:20s} a {best.x[0]:.4f}  spacing {best.x[1]:.2f} us")
    # 2. PRESET23
    P("")
    P("== PRESET23 / NON-LIN, tap 41, L ==")
    h = tap_h("PRESET23 SM SECT1", "NON-LIN", "N")
    i1 = int(1000 / F[1])
    mag = 20 * np.log10(np.abs(h) / np.abs(h[i1]))
    m = (F > 1000) & (F < 5000)
    ph = np.unwrap(np.angle(h))
    cf = np.polyfit(F[m], ph[m], 1)
    dev = np.degrees(ph - np.polyval(cf, F))
    rows = {}
    for fx in (2000, 4000, 6000, 7000, 8000, 8500, 9000, 9500, 10000, 10500, 10800, 11000, 11500, 12000):
        i = int(round(fx / F[1]))
        rows[fx] = (float(mag[i]), float(dev[i]))
    P("  f Hz: magnitude dB / phase from a pure delay (fitted 1-5 kHz), deg")
    P("  " + "  ".join(f"{k}: {v[0]:+.1f}/{v[1]:+.1f}" for k, v in rows.items()))
    mm = (F > 10500) & (F < 13000)
    d = mag[mm]
    ff = F[mm]
    notch = [(float(ff[i]), float(d[i])) for i in range(1, len(d) - 1) if d[i] < d[i - 1] and d[i] < d[i + 1] and d[i] < -35]
    P("  notches: " + ", ".join(f"{a:.0f} Hz {b:.1f} dB" for a, b in notch))
    out["p23"] = {"points": rows, "notches": notch}
    lh = np.log(h[FIT])
    lh = np.real(lh) + 1j * np.unwrap(np.imag(lh))
    tau0 = -cf[0] / (2 * np.pi) * 1e6

    def iir_fit(kind, N, starts, lb, ub):
        def H(q):
            if kind == "ellip":
                bb, aa = ellip(N, q[0], q[1], q[2] * 1e3, fs=FSD)
            else:
                bb, aa = cheby2(N, q[0], q[1] * 1e3, fs=FSD)
            return freqz(bb, aa, worN=f, fs=FSD)[1]

        def r(qq):
            d = lh - (np.log(H(qq[:-2]) + 1e-30) + qq[-2] - 2j * np.pi * f * qq[-1] * 1e-6)
            di = np.imag(d)
            di = di - 2 * np.pi * np.round(np.sum(W * di) / (2 * np.pi))
            return np.concatenate([sw * np.real(d), sw * di])
        best = None
        for s0 in starts:
            try:
                rr = least_squares(r, list(s0) + [0, tau0], bounds=(lb + [-np.inf, -1e4], ub + [np.inf, 1e4]))
            except ValueError:
                continue
            if best is None or rr.cost < best.cost:
                best = rr
        q = best.x
        d = lh - (np.log(H(q[:-2]) + 1e-30) + q[-2] - 2j * np.pi * f * q[-1] * 1e-6)
        dr = np.real(d) * 20 / np.log(10)
        di = np.degrees(np.imag(d))
        di -= 360 * np.round(np.sum(W * di) / 360)
        return [float(v) for v in q[:-2]], float(np.sqrt(np.sum(W * dr ** 2))), float(np.sqrt(np.sum(W * di ** 2)))
    out["p23"]["iir"] = {}
    for N in (2, 3, 4, 6):
        q, e1, e2 = iir_fit("ellip", N, [[0.5, 30, 8.0], [0.2, 40, 8.5], [1, 25, 7.5]], [0.01, 10, 3], [3, 80, 12.5])
        out["p23"]["iir"][f"ellip{N}"] = (q, e1, e2)
        P(f"  elliptic order {N} (3 numbers): rp {q[0]:.2f} rs {q[1]:.1f} fp {q[2]:.2f} kHz  {e1:.2f} dB / {e2:.1f} deg rms")
    for N in (3, 4, 6):
        q, e1, e2 = iir_fit("cheby2", N, [[30, 10.8], [40, 11]], [10, 8], [80, 13])
        out["p23"]["iir"][f"cheby2_{N}"] = (q, e1, e2)
        P(f"  Chebyshev II order {N} (2 numbers): rs {q[0]:.1f} fs {q[1]:.2f} kHz  {e1:.2f} dB / {e2:.1f} deg rms")
    # symmetric (linear-phase) FIRs, zero-phase amplitude signed past the first notch
    A = np.where(F < notch[0][0], np.abs(h) / np.abs(h[i1]), -np.abs(h) / np.abs(h[i1]))
    mf = (F >= 200) & (F <= 12300)
    pb = F[mf] < 10500
    out["p23"]["fir"] = {}
    P("  symmetric FIR (zero phase), K taps = (K+1)/2 numbers; rms / max dB 0.2-10.5 kHz, linear rms to 12.3 kHz:")
    for rate, nm in ((FSD, "device"), (48000.0, "48k")):
        row = []
        for K in (5, 7, 9, 11, 13, 15, 17, 21):
            c = sym_fir(F[mf], A[mf], K, rate)
            y = sym_amp(c, F[mf], rate)
            e = 20 * np.log10(np.maximum(np.abs(y[pb]), 1e-9) / A[mf][pb])
            out["p23"]["fir"][f"{nm}{K}"] = {"rms_db": float(np.sqrt(np.mean(e ** 2))), "max_db": float(np.abs(e).max()),
                                            "lin_rms": float(np.sqrt(np.mean((y - A[mf]) ** 2))), "coef": c.tolist()}
            row.append(f"{K}: {np.sqrt(np.mean(e ** 2)):.3f}/{np.abs(e).max():.2f}/{np.sqrt(np.mean((y - A[mf]) ** 2)):.4f}")
        P(f"   {nm:6s} " + "  ".join(row))
    c13 = np.array(out["p23"]["fir"]["device13"]["coef"])
    P(f"  device 13-tap (centre..edge): {', '.join(f'{v:.4f}' for v in c13)}; 48 kHz 17-tap coefficients reach "
      f"{np.abs(out['p23']['fir']['48k17']['coef']).max():.0f} (ill-conditioned)")
    # 3. the candidate law: PRESET23's device 13-tap shape scaled to each step (H: one observation, taken as step 8)
    P("")
    P("== Candidate law (H): PRESET23 shape at step 8, frequency-scaled to each step (held at its Nyquist stop-band level, -27 dB, past the scaled Nyquist), against the plugin Butterworth ==")
    P("  step kHz | law: -3 dB kHz, level at 5 / 8 / 10 / 12 kHz | Butterworth: -3 dB kHz, level at 5 / 8 / 10 / 12 kHz")
    out["law"] = {}

    def law_db(fx, step):
        return 20 * np.log10(np.maximum(np.abs(sym_amp(c13, np.minimum(np.asarray(fx) * 8.0 / step, FSD / 2), FSD)) / c13.dot([1] + [2] * 6), 1e-6))
    for st in LP_STEPS[:-1]:
        lw = [float(law_db(np.array([x * 1e3]), st)[0]) for x in (5, 8, 10, 12)]
        bw = [float(20 * np.log10(np.abs(h_bw(np.array([x * 1e3]), st * 1e3, 2))[0])) for x in (5, 8, 10, 12)]
        l3 = m3db(lambda g: law_db(g, st))
        b3 = m3db(lambda g: 20 * np.log10(np.abs(h_bw(g, st * 1e3, 2))))
        out["law"][st] = {"law_m3db_khz": l3 and l3 / 1e3, "law_db": lw, "bw2_m3db_khz": b3 and b3 / 1e3, "bw2_db": bw}
        P(f"  {st:>4} | {l3 / 1e3 if l3 else float('nan'):5.2f}  " + " ".join(f"{v:+7.1f}" for v in lw)
          + f" | {b3 / 1e3 if b3 else float('nan'):5.2f}  " + " ".join(f"{v:+7.1f}" for v in bw))
    P("    12 | open (exact bypass) in both")
    # the 2-tap totals for reference (what the dark presets are)
    P("")
    P("  the dark presets' total 2-tap kernels [1 - a, a] for reference: -3 dB kHz, level at 5 / 8 / 10 / 12 kHz")
    for a in (0.047, 0.061, 0.091, 0.23, 0.333, 0.411, 0.498):
        lv = [float(20 * np.log10(np.abs(h_fir2(np.array([x * 1e3]), a))[0])) for x in (5, 8, 10, 12)]
        a3 = m3db(lambda g: 20 * np.log10(np.abs(h_fir2(g, a))))
        P(f"   a {a:.3f}: {a3 / 1e3 if a3 else float('nan'):5.2f}  " + " ".join(f"{v:+6.1f}" for v in lv))
    json.dump(out, open(os.path.join(OUT, "extra.json"), "w"), indent=1)
    open(os.path.join(OUT, "extra.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    import argparse
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--stage", default="all", choices=["fit", "extra", "all"])
    st_ = ap_.parse_args().stage
    if st_ in ("fit", "all"):
        stage_fit()
    if st_ in ("extra", "all"):
        stage_extra()
