"""Fork C: the complex residual above 11 kHz measured on the FRONT TAPS only (no reverb), per channel, pooled over captures.
Per capture and channel: the cached device-rate IR from file start to just before the first branch (main) tap carries only the
front taps. Ideal = sum_i a_i exp(-i w T_i) * R_ch(f) * G_dg(f) (nl_structure.processing_response x resampler_gain), a_i from the
front law, T_i free to sub-sample and one real gain g per channel, fitted on 2-10 kHz where the residual is ~0; the complex ratio
capture / ideal is then read in 10-12.8 kHz on the same knots as resid_eq.py (hat functions, weights |ideal|^2).
usage: c_front.py N_CAPTURES [outfile]   (no cache, no shared file is touched)"""
import os, sys, json, random
os.environ["NL_WARP2"] = "1"
os.environ.pop("NL_BK", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.optimize import least_squares
import defects

KNOTS = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
NF = 4096                       # FFT size of the segment spectra
NB = 1 << 14                    # size of the model's time grid
FITBAND = (2000.0, 10000.0)
SPAN_BEFORE_BRANCH = 45         # samples kept clear before the first main tap (its band-limited pre-ringing)
B0 = {0: 685, 1: 464}           # front-frame position of the first branch tap, L / R (R's 463 front tap sits on it: excluded)
TAPS_USED = {0: (0, 1, 2), 1: (0, 1)}       # indices into ne.FRONT[ch]; variant 'no41' drops index 0


def hats(f):
    W = []
    for kk in range(len(KNOTS)):
        lo = KNOTS[kk - 1] if kk else KNOTS[0] - 1; hi = KNOTS[kk + 1] if kk < len(KNOTS) - 1 else KNOTS[-1] + 1
        W.append(np.clip(np.minimum((f - lo) / (KNOTS[kk] - lo), (hi - f) / (hi - KNOTS[kk])), 0, 1))
    return np.array(W)


def usable(key, b):
    """front taps visible: not trimmed, and every used front tap within 10 dB of the strongest, both channels"""
    if b.get("trimmed"): return False
    for ch in (0, 1):
        fdb, _ = ne.law_db(key[0], key[1], ch)
        amp = np.abs([ne.FRONT_AMP[ch][i] * 10 ** (fdb[i] / 20) for i in range(3)])
        used = [amp[i] for i in TAPS_USED[ch]]
        if min(used) < 0.25 * max(amp.max(), 1e-9) or max(used) < 0.5: return False
    return True


_PROC = {}
def proc(ch):
    if "f" not in _PROC:
        _PROC["f"] = np.fft.rfftfreq(NB, 1 / FS)
        _PROC["R"] = nl.processing_response(_PROC["f"]) * nl.resampler_gain(_PROC["f"])[:, None]
    return _PROC["f"], _PROC["R"][:, ch]


def fit_channel(key, ch, b, h, used):
    """returns None or dict(num, den, taps, g, nulls...) for one channel"""
    fr = b["front_shift"]
    fdb, _ = ne.law_db(key[0], key[1], ch)
    a = np.array([ne.FRONT_AMP[ch][i] * 10 ** (fdb[i] / 20) for i in used])
    s1 = int(M.K0 + fr + B0[ch] - SPAN_BEFORE_BRANCH)
    s0 = 0 if used[0] == 0 else int(M.K0 + fr + ne.FRONT[ch][used[0]] - 100)      # without tap 41 the window starts after it
    y = h[s0:s1, ch]
    f, P = proc(ch)
    fnf = np.fft.rfftfreq(NF, 1 / FS)
    Y = np.fft.rfft(y, NF)
    T0 = np.array([M.K0 + fr + ne.FRONT[ch][i] for i in used])

    def model(d):
        S = sum(ai * np.exp(-2j * np.pi * f / FS * (T0[j] + d[j])) for j, ai in enumerate(a)) * P
        m = np.fft.irfft(S, NB)
        return np.fft.rfft(m[s0:s1], NF)
    sel = (fnf >= FITBAND[0]) & (fnf < FITBAND[1])

    def resid(d):
        Mm = model(d)[sel]
        g = np.real(np.vdot(Mm, Y[sel])) / np.real(np.vdot(Mm, Mm))
        r = Y[sel] - g * Mm
        return np.r_[r.real, r.imag]
    sol = least_squares(resid, np.zeros(len(a)), x_scale=0.1, bounds=(-1.5, 1.5))
    d = sol.x
    Mm = model(d)
    g = np.real(np.vdot(Mm[sel], Y[sel])) / np.real(np.vdot(Mm[sel], Mm[sel]))
    Mg = g * Mm
    W = hats(fnf)
    num = (W * (Y * np.conj(Mg))[None, :]).sum(1); den = (W * (np.abs(Mg) ** 2)[None, :]).sum(1)
    inb = (fnf >= 3000) & (fnf < 9000)
    fitnull = 10 * np.log10(np.sum(np.abs(Y[inb] - Mg[inb]) ** 2) / np.sum(np.abs(Y[inb]) ** 2))
    return dict(num=num, den=den, d=d.tolist(), g=float(g), fitnull=float(fitnull))


def run(keys, variant):
    out = []
    for key in keys:
        b = common.bench(key); h = nl.ir(nl.table()[key])
        row = dict(key=key)
        for ch in (0, 1):
            used = TAPS_USED[ch] if variant == "all" else tuple(i for i in TAPS_USED[ch] if i != 0)
            r = fit_channel(key, ch, b, h, used)
            row[ch] = r
        out.append(row)
    return out


def pooled(rows, ch):
    num = sum(r[ch]["num"] for r in rows); den = sum(r[ch]["den"] for r in rows)
    return num / den


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
    random.seed(23)
    ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
    ks = [k for k in ks if usable(k, common.bench(k))]
    print(len(ks), "captures with usable front taps")
    sample = random.sample(ks, min(n, len(ks)))
    res = {}
    for variant in ("all", "no41"):
        rows = run(sample, variant)
        res[variant] = rows
        print("== variant", variant, "n", len(rows))
        for ch in (0, 1):
            fn = np.median([r[ch]["fitnull"] for r in rows])
            c = pooled(rows, ch)
            # bootstrap over captures for the spread
            rng = np.random.default_rng(0); bs = []
            for _ in range(100):
                idx = rng.integers(0, len(rows), len(rows))
                nn = sum(rows[i][ch]["num"] for i in idx); dd = sum(rows[i][ch]["den"] for i in idx); bs.append(nn / dd)
            bs = np.array(bs)
            print(f" ch{ch}  median 3-9 kHz fit null {fn:.1f} dB; tap offsets (median, samples) {np.round(np.median([r[ch]['d'] for r in rows], 0), 3)}")
            print("   knots kHz      ", np.round(KNOTS / 1000, 2))
            print("   level dB       ", np.round(20 * np.log10(np.abs(c)), 3))
            print("   phase mrad     ", np.round(np.angle(c) * 1e3, 1))
            print("   phase mrad s.e.", np.round(np.std(np.angle(bs), 0) * 1e3, 1), " level s.e. dB", np.round(np.std(20 * np.log10(np.abs(bs)), 0), 3))
    import pickle; pickle.dump(res, open(OUT + "/" + (sys.argv[2] if len(sys.argv) > 2 else "c_front_smoke.pkl"), "wb"))
