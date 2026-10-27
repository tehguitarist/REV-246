"""Shared helpers for the 11-12.3 kHz drift investigation (read-only use of tools/nl_*; writes only here)."""
import json, os, sys
import numpy as np
ROOT = "/Users/leighpierce/Git/246"
sys.path.insert(0, os.path.join(ROOT, "tools"))
import nl_model as M  # noqa: E402  (imports nl_early / nl_late / nl_structure; no side effects on bench dirs)
nl, ne, nlt = M.nl, M.ne, M.nlt
FS = M.FS
OUT = os.path.join(ROOT, "analysis", "nl_hfdrift")
BENCH = os.path.join(ROOT, "analysis", "nl", "bench")


def bench(key):
    return json.load(open(os.path.join(BENCH, "_".join(str(v) for v in key) + ".json")))


def pair(key, full=True, proc=True):
    """(y, m, info): the capture (clock-unwarped, full band) and the model at the capture's own bench fit, both (n, 2),
    gain applied to the model, plus the front-frame start index."""
    b = bench(key)
    raw = nl.ir(nl.table()[key])
    h = M.sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]
    end = M.usable_end(h)
    n = min(end, M.K0 + int(fr) + int(1.0 * FS)) if False else end
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = nlt.unwarp(rawp, b["ppm"])[:n]
    m = M.render(key, n, b["o"], fr, band=False, proc=proc) * np.array(b["gain"])
    return y, m, dict(fr=fr, n=n, b=b, i0=M.window_start(h, fr)[0])


def local_delay(y, m, f_lo, f_hi, win=2048, hop=1024, start=0, stop=None, thr_db=-60):
    """Per STFT frame: delay (device samples, y relative to m: positive = y later) from the cross spectrum in the band,
    phase-delay (through zero) and group-delay (free intercept) estimates, weights = |m|^2 energy. Returns arrays
    (t_centre_s, tau_phase, tau_group, coh) for one channel pair (1-D inputs)."""
    w = np.hanning(win)
    f = np.fft.rfftfreq(win)
    sel = (f * FS >= f_lo) & (f * FS < f_hi)
    out = []
    stop = len(y) if stop is None else stop
    for a in range(start, stop - win, hop):
        Y = np.fft.rfft(y[a:a + win] * w)[sel]
        Mm = np.fft.rfft(m[a:a + win] * w)[sel]
        x = Y * np.conj(Mm)
        e = np.sum(np.abs(Mm) ** 2)
        if e <= 0 or np.sum(np.abs(Mm) ** 2) < 1e-12:
            continue
        ph = np.angle(x)
        wt = np.abs(x)
        ff = f[sel]
        # phase delay: y(t) ~ m(t - tau)  =>  Y = M e^{-i 2 pi f tau}  => arg(Y conj M) = -2 pi f tau
        tau_p = -np.sum(wt * ph * ff) / (2 * np.pi * np.sum(wt * ff ** 2))
        A = np.vstack([ff, np.ones_like(ff)]).T * np.sqrt(wt)[:, None]
        sol = np.linalg.lstsq(A, ph * np.sqrt(wt), rcond=None)[0]
        tau_g = -sol[0] / (2 * np.pi)
        coh = np.abs(np.sum(x)) / np.sum(wt)
        out.append(((a + win / 2) / FS, tau_p, tau_g, coh, 10 * np.log10(e + 1e-30)))
    return np.array(out)


def unwarp2(x, ppm, up=32, half=150, beta=9.0):
    """unwarp with a long interpolation kernel: transition ~0.03 cycles/sample (12.8 -> 13.6 kHz) instead of the default
    Kaiser(5) with 10 taps per side (about 3.9 kHz wide, which leaves the image of 11-13 kHz unsuppressed)."""
    from scipy.signal import firwin, resample_poly
    h = firwin(2 * half * up + 1, 1.0 / up, window=("kaiser", beta))
    u = resample_poly(x, up, 1, axis=0, window=h)
    t = np.arange(len(x)) * (1 + ppm * 1e-6) * up
    return np.stack([np.interp(t, np.arange(len(u)), u[:, c]) for c in range(x.shape[1])], 1)


def unwarp_kernel_response(up=32, f_hz=None, half=10, beta=5.0):
    """Gain of the interpolation filter at original-rate frequencies f_hz (default scipy kernel)."""
    from scipy.signal import firwin
    h = firwin(2 * half * up + 1, 1.0 / up, window=("kaiser", beta)) * up
    w = 2 * np.pi * np.asarray(f_hz) / FS / up
    k = np.arange(len(h)) - half * up
    return np.array([abs(np.sum(h * np.exp(-1j * wi * k))) / up for wi in w])
