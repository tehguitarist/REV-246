"""Program N (Non-Lin): the late section, from the Tomkin 96 kHz library (docs/NonLin_program.md §10).

    .venv/bin/python tools/nl_late.py --stage excite

The late section is fed from the main-tap line (NonLin_program.md 9.3; feed tap per SPAN in FEED) through AP(124) and
depends on LATE only: RT = LATE, flat in frequency. Every (SPAN, EARLY) is therefore another capture of the same late
network at each LATE.

excite: the SPAN 0 captures (EARLY 10 and -3: the late at 0 dB, nothing but the front taps before it), all ten LATE
values, both channels. Each capture is unwarped by the session clock (CLOCK_PPM), its front taps blanked, aligned and
scaled on the late's first pass (which no loop reaches), exponentially windowed, and then, per frequency bin, solved
by least squares for

    Y_c = X [F + sum_i E_i / (1 - g_i(LATE_c) z^-P_i)],   g_i = 10^(-3 P_i / (fs LATE_c)),

X = AP(124), P_i = LOOPS: one LATE-independent part F and one excitation E_i per loop, shared by every capture.
Writes analysis/nl/late_excite.npz (time-domain F and E_i per channel, re the feed) and prints the solve's residual
per capture and the excitations' largest taps.
"""
import argparse, functools, os, sys
import numpy as np
from scipy.signal import firwin, lfilter, resample_poly, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nl_early as ne  # noqa: E402
import nl_structure as nl  # noqa: E402

FS, K0 = nl.FSD, nl.K0
BAND = ne.BAND
LOOPS = [1108, 1254, 1386]
FEED = {0.0: 684, 0.1: 1766, 0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}   # one sample before a main tap (NonLin_program 10)
R_LATE = -220                                    # R's late section re L's
CLOCK_PPM = 11.6                                 # the session's capture clock (median of the early fits, +-0.5)
LATES = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
OUT = os.path.join(nl.OUT, "late_excite.npz")
AP_LATE = 314                                    # the late path's second all-pass (after AP124), g = -1/sqrt(2)
# four equal read taps per loop (re the late feed), the same positions in both channels, per-channel signs
READS = {1108: (301, 712, 923, 1107), 1254: (398, 689, 1022, 1253), 1386: (582, 897, 1179, 1385)}
READ_SIGNS = {0: {1108: (-1, 1, 1, 1), 1254: (1, -1, 1, 1), 1386: (-1, -1, 1, -1)},
              1: {1108: (1, -1, 1, 1), 1254: (1, -1, -1, -1), 1386: (1, -1, 1, 1)}}


def loop_gain(P, late):
    """Jot gain of a flat loop of P samples for RT = LATE (s)."""
    return 10 ** (-3 * P / (FS * late))


def late_columns(ch, n, span, late, o_branch, kappa=0.0, r_off=R_LATE):
    """(n, 12) band-limited columns of the late section, one per (loop, read tap), unsigned and at unit gain: the input
    through AP(124) and AP(AP_LATE), delayed to the feed (FEED[span] + o_branch, R r_off earlier), through the loop
    (P, Jot gain), read at the tap. Positions scale with the capture clock (1 + kappa)."""
    N = 1 << int(np.ceil(np.log2(n + 4 * FEED[0.5])))
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    u = np.zeros(N)
    u[0] = 1.0
    u = ne.ap(ne.ap(u, ne.AP_IN), AP_LATE)
    cols = []
    for P, qs in READS.items():
        v = np.zeros(N)
        v[0] = 1.0
        g = loop_gain(P, late)
        v = lfilter([1.0], np.r_[1.0, np.zeros(P - 1), -g], v)              # the comb's impulse response
        V = np.fft.rfft(np.fft.irfft(np.fft.rfft(u) * np.fft.rfft(v), N), N)
        for q in qs:
            pos = K0 + o_branch + (FEED[span] + (r_off if ch else 0) + q) * (1 + kappa)
            cols.append(np.fft.irfft(V * np.exp(ph * pos), N)[:n])
    return sosfiltfilt(BAND, np.array(cols).T, axis=0)


def late_column(ch, n, span, late, o_branch, kappa=0.0, r_off=R_LATE):
    """The late section of one channel at unit read-tap gain, the read taps' signs applied (n,)."""
    sg = [s for P in READS for s in READ_SIGNS[ch][P]]
    return late_columns(ch, n, span, late, o_branch, kappa, r_off) @ np.array(sg, float)


# The clock warp with a long interpolation kernel (NonLin_program 11.6; default, NL_WARP2=0 restores the old one). scipy's default Kaiser(5) kernel
# (10 input samples per side) has a 3.9 kHz transition band around Nyquist, which leaves the image of 11-13 kHz in the
# 32x signal: a fractional-position resample then gives G(f) + G(fs - f) exp(i 2 pi frac), a phase error growing with
# file time and an amplitude loss at 11-12.8 kHz. The long kernel (150 input samples per side, Kaiser 9) is flat to
# 12.8 kHz and rejects from 13.6 kHz, where the captures' analysis band ends. The old warp is kept for the legacy benches.
WARP2 = os.environ.get("NL_WARP2", "1") == "1"
WARP_HALF, WARP_BETA = 150, 9.0


@functools.lru_cache(maxsize=1)
def _warp_kernel(up):
    return firwin(2 * WARP_HALF * up + 1, 1.0 / up, window=("kaiser", WARP_BETA))


def upsample(x, up=32):
    """x (n, c) on the 32x grid; the part of unwarp that does not depend on the clock ratio (cache it per capture)."""
    return resample_poly(x, up, 1, axis=0, window=_warp_kernel(up)) if WARP2 else resample_poly(x, up, 1, axis=0)


def unwarp_from(u, n, ppm, up=32):
    """The first n samples of the clock-corrected capture from its 32x grid u (see unwarp)."""
    t = np.arange(n) * (1 + ppm * 1e-6) * up
    return np.stack([np.interp(t, np.arange(len(u)), u[:, c]) for c in range(u.shape[1])], 1)


def unwarp(x, ppm, up=32):
    """x resampled from the capture clock to the device clock: sample n of the result is x at n (1 + ppm 1e-6)."""
    return unwarp_from(upsample(x, up), len(x), ppm, up)


def late_capture(key, n):
    """Band-limited, unwarped device-rate IR with the late's first L arrival near K0 + FEED + 301 (front frame),
    front taps blanked, n samples."""
    h = sosfiltfilt(BAND, unwarp(nl.ir(key), CLOCK_PPM), axis=0)
    fr = ne.front_shift(h)
    s = int(round(fr))
    y = np.zeros((n, 2))
    m = min(n, len(h) - s)
    y[:m] = h[s:s + m]
    y[:K0 + FEED[key[0]] + 250] = 0.0
    return y


def excite_stage(n=40000, r_db_s=60.0):
    keys = [(0.0, e, l) for e in (10.0, -3.0) for l in LATES]
    ys = {k: late_capture(k, n) for k in keys}
    ref = ys[(0.0, 10.0, 2.0)]
    i0, i1 = K0 + 685 + 250, K0 + 685 + 1050                 # first pass, before any recirculation (301 + 1108)
    # align (sub-sample) and scale every capture onto the reference on the first pass, per channel
    for k in keys:
        for ch in (0, 1):
            a = ref[i0:i1, ch]
            best = None
            for s in np.arange(-3, 3.01, 0.05):
                b = ne.fshift(ys[k][:, ch], -s)[i0:i1]
                g = (a @ b) / (b @ b)
                e = np.sum((a - g * b) ** 2)
                if best is None or e < best[0]:
                    best = (e, s, g)
            ys[k][:, ch] = best[2] * ne.fshift(ys[k][:, ch], -best[1])
            print(k, ch, f"shift {best[1]:+.2f} gain {20 * np.log10(best[2]):+.2f} dB "
                         f"first-pass null {10 * np.log10(best[0] / np.sum(a * a)):.1f}", flush=True)
    r = 10 ** (-r_db_s / 20 / FS)                              # exponential window: -60 dB per second
    w = r ** np.arange(n)
    N = 1 << int(np.ceil(np.log2(2 * n)))
    om = np.fft.rfftfreq(N)
    x = np.zeros(n)
    x[K0] = 1.0
    X = np.fft.rfft(ne.ap(x, ne.AP_IN) * w, N)
    out = {}
    for ch in (0, 1):
        feed = 685 + (R_LATE if ch else 0)
        Y = np.array([np.fft.rfft(ys[k][:, ch] * w, N) for k in keys])           # (C, B)
        z = np.exp(-2j * np.pi * om)
        cols = [np.ones_like(z)[None, :].repeat(len(keys), 0)]
        for P in LOOPS:
            g = np.array([10 ** (-3 * P / (FS * k[2])) for k in keys])[:, None]
            cols.append(1 / (1 - g * r ** P * z[None, :] ** P))
        A = np.stack(cols, 2) * X[None, :, None]                                  # (C, B, 1 + loops)
        sol = np.zeros((len(om), A.shape[2]), complex)
        res = np.zeros(len(keys))
        for b in range(len(om)):
            s_, *_ = np.linalg.lstsq(A[:, b, :], Y[:, b], rcond=None)
            sol[b] = s_
            res += np.abs(Y[:, b] - A[:, b, :] @ s_) ** 2
        for c, k in enumerate(keys):
            print(f"ch {ch} {k}: solve residual {10 * np.log10(res[c] / np.sum(np.abs(Y[c]) ** 2)):.1f} dB")
        parts = np.fft.irfft(sol, N, axis=0)[:n] / w[:, None]                     # time domain, unweighted
        parts = np.roll(parts, -(K0 + feed), axis=0)[:6000]                       # re the feed
        out[f"ch{ch}"] = parts
        for j, name in enumerate(["F"] + [f"E{P}" for P in LOOPS]):
            p = parts[:, j]
            m = np.abs(p).max()
            idx = np.argsort(np.abs(p))[::-1]
            taps = []
            for i in idx:
                if all(abs(i - q) > 3 for q in taps):
                    taps.append(int(i))
                if len(taps) == 14:
                    break
            print(f"  ch{ch} {name}: energy {10 * np.log10(np.sum(p ** 2)):.1f} dB, taps "
                  + " ".join(f"{i}:{20 * np.log10(abs(p[i]) / m):.0f}{'+' if p[i] > 0 else '-'}" for i in sorted(taps)))
    np.savez(OUT, **out)


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--stage", default="excite")
    a = ap_.parse_args()
    if a.stage == "excite":
        excite_stage()


if __name__ == "__main__":
    main()
