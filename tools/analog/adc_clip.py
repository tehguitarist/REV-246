"""The A/D's view of the input and its clipping, for the plugin's 16-bit mode (October 2026).

The unit's A/D (service manual pp. 37-39) is a discrete 16-bit SAR behind the input filter 7 246 131 and an
integrate/hold S&H. A SAR cannot wrap: past full scale every comparison gives the same answer, so the code sticks at
the end of the range, a hard clip of the held sample at the device rate. What the A/D sees is the input filter's
output, with its +5.1 dB HF pre-emphasis, so the clip point is not where a plain level reading puts it.

This script:
  1. extracts the exact poles and zeros of the transcribed input filter (netlists.py) from the MNA matrices
     (generalised eigenvalues of G + sC), checked against mna.response;
  2. designs the reinjection kernel g0 = IFT{T(w) / H_band(jw)}: the clip error found in the A/D's domain, brought
     back to the input domain so that a later input filter (the engine's analogue chain) turns it into the A/D's own
     error. T is a taper that zeroes it where the filter stops passing (and doubles as the device-rate reconstruction
     low-pass);
  3. runs the whole scheme as the plugin does (host-rate filter, sinc sampling at the device instants, clip, sparse
     reconstruction, delayed dry) and checks the input filter applied to its output against the ideal: the clipped
     A/D samples, delayed.

Writes analysis/analog/adc_clip.json. Usage: .venv/bin/python tools/analog/adc_clip.py
"""
import json
import os
import sys

import numpy as np
import scipy.linalg as sl
from scipy import signal

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from mna import response  # noqa: E402
from netlists import INPUT_FILTER  # noqa: E402

FS_DEV = 19.2e6 / 728
TWO_PI = 2 * np.pi
# Kernel and interpolator choices (mirrored in plugin/src/dsp/AdcClipper.h).
TAPER_LO, TAPER_HI = 11.5e3, 12.6e3   # T(w): 1 below, raised cosine, 0 above (device Nyquist 13.19 kHz)
PRE, POST = 1.2e-3, 0.3e-3            # g0 support around its centre, seconds (it is mostly anticausal)
INTERP_HALF, INTERP_FC, INTERP_BETA = 0.2e-3, 22e3, 7.0  # sampling the host-rate A/D signal at device instants
COUPLING_HZ = 3.35                    # C106 10 uF / R117 4.75 k into the A/D summing point (unit facts 2)


def pole_zero(net):
    elems, amps = net["elements"], net["opamps"]
    nodes = sorted(({n for e in elems for n in e[1:3]} | {n for a in amps for n in a[:3]}) - {"gnd"})
    idx = {n: i for i, n in enumerate(nodes)}
    nn, na = len(nodes), len(amps)
    N = nn + na + 1
    G, C = np.zeros((N, N)), np.zeros((N, N))
    for kind, a, b, val, _ in elems:
        M, y = (G, 1 / val) if kind == "R" else (C, val)
        for p_, q_ in ((a, b), (b, a)):
            if p_ != "gnd":
                M[idx[p_], idx[p_]] += y
                if q_ != "gnd":
                    M[idx[p_], idx[q_]] -= y
    for j, (pl, mi, ou, _) in enumerate(amps):
        r = nn + j
        G[idx[ou], r] += 1
        if pl != "gnd":
            G[r, idx[pl]] += 1
        if mi != "gnd":
            G[r, idx[mi]] -= 1
    G[idx["in"], nn + na] += 1
    G[nn + na, idx["in"]] = 1
    w = sl.eigvals(G, -C)
    poles = w[np.isfinite(w) & (np.abs(w) < 1e9)]
    Gb, Cb = np.zeros((N + 1, N + 1)), np.zeros((N + 1, N + 1))
    Gb[:N, :N], Cb[:N, :N] = G, C
    Gb[N, idx[net["out"]]] = 1
    Gb[N - 1, N] = 1
    w = sl.eigvals(Gb, -Cb)
    zeros = w[np.isfinite(w) & (np.abs(w) < 1e9)]
    clean = lambda r: np.where(np.abs(r.real) < 1e-6 * np.abs(r), 1j * r.imag, r)
    return clean(poles), clean(zeros)


def upper(r):
    """One of each conjugate pair (imag >= 0), sorted by frequency."""
    r = r[r.imag >= -1e-9]
    return r[np.argsort(np.abs(r))]


def h_analog(s, poles, zeros):
    return np.prod([s - z for z in zeros], 0) / np.prod([s - p for p in poles], 0)


def matched_sos(poles, zeros, fs, f_norm=1e3, f_fit=12.6e3):
    """Matched-z (z = e^{sT}: every pole and zero at its own frequency; the bilinear transform warps the steep 12 kHz
    skirt by up to 1 kHz at 44.1 kHz) plus one zero at z = -alpha, alpha in [0, 1] chosen for the smallest peak
    magnitude error to f_fit. Normalised to unity at f_norm (the plugin's method). Returns (sos, alpha, max_err_db, delay_s)."""
    T = 1 / fs
    f = np.linspace(20, f_fit, 64)
    ha = 20 * np.log10(np.abs(h_analog(TWO_PI * 1j * f, poles, zeros)))
    ha -= 20 * np.log10(np.abs(h_analog(TWO_PI * 1j * np.array([f_norm]), poles, zeros)))[0]
    best = None
    for alpha in np.linspace(0, 1, 201):
        sos = signal.zpk2sos(np.r_[np.exp(zeros * T), -alpha], np.exp(poles * T), 1.0)
        _, h1 = signal.sosfreqz(sos, [f_norm], fs=fs)
        sos[0, :3] /= np.abs(h1[0])
        _, hd = signal.sosfreqz(sos, f, fs=fs)
        e = float(np.max(np.abs(20 * np.log10(np.abs(hd)) - ha)))
        if best is None or e < best[2]:
            best = (sos, float(alpha), e)
    # Matched-z lags the analogue filter by a near-constant delay (phase linear in f to 0.1 deg): measure it at
    # 1 kHz. The A/D samples then belong to analogue time t - delay, and the reinjection is placed accordingly.
    _, h1 = signal.sosfreqz(best[0], [f_norm], fs=fs)
    ph = np.angle(h1[0] / h_analog(TWO_PI * 1j * np.array([f_norm]), poles, zeros)[0])
    return best + (float(-ph / (TWO_PI * f_norm)),)


def taper(f):
    f = np.abs(f)
    t = np.clip((f - TAPER_LO) / (TAPER_HI - TAPER_LO), 0, 1)
    return np.where(f <= TAPER_LO, 1.0, 0.5 * (1 + np.cos(np.pi * t)))


def kernel_g0(poles_b, zeros_b, norm, tau):
    """g0(tau) = (1/pi) int_0^inf Re{T(w) / H_band(jw) e^{jw tau}} dw, H_band normalised by `norm`."""
    f = np.arange(0, TAPER_HI + 5, 5.0)
    G = taper(f) * norm / h_analog(TWO_PI * 1j * f, poles_b, zeros_b)
    wts = np.full(len(f), 5.0)
    wts[0] = wts[-1] = 2.5                                      # trapezoid in Hz
    ph = np.exp(1j * TWO_PI * np.outer(tau, f))
    return 2 * np.real(ph @ (G * wts))                          # 2 int_0 Re{} df = (1/pi) int_0 Re{} dw


def kaiser_sinc(tau):
    x = np.clip(tau / INTERP_HALF, -1, 1)
    w = np.i0(INTERP_BETA * np.sqrt(1 - x * x)) / np.i0(INTERP_BETA)
    return 2 * INTERP_FC * np.sinc(2 * INTERP_FC * tau) * w * (np.abs(tau) < INTERP_HALF)


def run(x, fs, full_scale, sos, g0_tau, g0, t0=0.0, lag=0.0):
    """The plugin's scheme. Returns (out, delay_samples, ad_samples, clip_errors, device_times)."""
    v = signal.sosfilt(sos, x)
    a = np.exp(-TWO_PI * COUPLING_HZ / fs)                     # A/D coupling high-pass (one-pole, matched DC)
    v = signal.lfilter([(1 + a) / 2, -(1 + a) / 2], [1, -a], v)
    T = 1 / fs
    D = int(np.ceil((PRE + INTERP_HALF) / T)) + 1               # dry delay, host samples
    tk = t0 + np.arange(int((len(x) / fs - 2 * INTERP_HALF - t0) * FS_DEV)) / FS_DEV
    vk = np.empty(len(tk))
    half = int(np.ceil(INTERP_HALF / T))
    for k, t in enumerate(tk):
        j0 = int(np.floor(t / T))
        j = np.arange(j0 - half, j0 + half + 2)
        j = j[(j >= 0) & (j < len(x))]
        vk[k] = T * np.sum(v[j] * kaiser_sinc(t - j * T))
    hi = full_scale * 32767 / 32768
    ek = np.clip(vk, -full_scale, hi) - vk
    out = np.r_[np.zeros(D), x][: len(x)]
    for k in np.nonzero(ek)[0]:
        tc = tk[k] + D * T - lag
        m = np.arange(int(np.ceil((tc - PRE) / T)), int(np.floor((tc + POST) / T)) + 1)
        m = m[(m >= 0) & (m < len(x))]
        out[m] += ek[k] / FS_DEV * np.interp(m * T - tc, g0_tau, g0)
    return out, D, vk, ek, tk


def main():
    poles, zeros = pole_zero(INPUT_FILTER)
    f = np.array([20, 100, 1e3, 2e3, 5e3, 10e3, 10.5e3, 12e3, 12.2e3, 13e3, 14e3])
    s = TWO_PI * 1j * f
    h_mna = response(INPUT_FILTER, f)
    k = h_mna[2] / h_analog(s[2:3], poles, zeros)[0]
    assert np.max(np.abs(k * h_analog(s, poles, zeros) / h_mna - 1)) < 1e-6
    # The band part: drop the filter's own sub-Hz coupling (pole near 0.6 Hz, zero at DC), which is not inverted.
    hp_p = poles[np.abs(poles) < TWO_PI * 5]
    hp_z = zeros[np.abs(zeros) < TWO_PI * 1e-3]
    pb = poles[np.abs(poles) >= TWO_PI * 5]
    zb = zeros[np.abs(zeros) >= TWO_PI * 1e-3]
    norm = np.abs(h_analog(TWO_PI * 1j * np.array([1e3]), pb, zb))[0]
    db = lambda h: 20 * np.log10(np.abs(h))
    resp = {f"{fi:.0f}": round(float(d), 3) for fi, d in zip(f, db(h_analog(s, poles, zeros) / h_analog(s[2:3], poles, zeros)))}
    print("input filter re 1 kHz:", resp)

    # g0 support: how much of the ideal (untruncated) kernel lies outside [-PRE, POST]?
    tau_wide = np.arange(-4e-3, 4e-3, 2e-6)
    gw = kernel_g0(pb, zb, norm, tau_wide)
    inside = (tau_wide >= -PRE) & (tau_wide <= POST)
    leak_db = 10 * np.log10(np.sum(gw[~inside] ** 2) / np.sum(gw ** 2))
    print(f"g0 energy outside [-{PRE*1e3} ms, {POST*1e3} ms]: {leak_db:.1f} dB")
    tau = np.arange(-PRE, POST + 1e-9, 0.5e-6)
    g0 = kernel_g0(pb, zb, norm, tau)
    # Truncation: a short raised-cosine fade over the last 50 us of each end.
    fade = np.clip(np.minimum(tau + PRE, POST - tau) / 50e-6, 0, 1)
    g0 = g0 * 0.5 * (1 - np.cos(np.pi * fade))

    # Interpolator flatness (sampling the host-rate A/D signal at the device instants).
    out = {"input_filter": {
        "poles_rad_s": [[float(p.real), float(p.imag)] for p in upper(poles)],
        "zeros_rad_s": [[float(z.real), float(z.imag)] for z in upper(zeros)],
        "gain_re_1khz_db": resp,
        "group_delay_us": {}}, "kernel": {"pre_ms": PRE * 1e3, "post_ms": POST * 1e3,
                                         "energy_outside_db": round(float(leak_db), 1),
                                         "taper_hz": [TAPER_LO, TAPER_HI]}, "rates": {}}
    fg = np.linspace(10, 12e3, 4000)
    hg = h_analog(TWO_PI * 1j * fg, pb, zb)
    gd = -np.gradient(np.unwrap(np.angle(hg)), TWO_PI * fg)
    for fq in (100, 1e3, 5e3, 10e3, 11.5e3, 12e3):
        out["input_filter"]["group_delay_us"][f"{fq:.0f}"] = round(float(np.interp(fq, fg, gd) * 1e6), 1)

    rng = np.random.default_rng(1)
    for fs in (44100.0, 48000.0, 96000.0):
        sos, alpha, dig_err, lag = matched_sos(pb, zb, fs)
        # Interpolator (continuous kernel): response at 14 kHz, rejection of the first image of 14 kHz at fs - 14 kHz.
        tt = np.linspace(-INTERP_HALF, INTERP_HALF, 20001)
        ir = lambda fq: abs(np.trapezoid(kaiser_sinc(tt) * np.exp(-TWO_PI * 1j * fq * tt), tt))
        interp = {"14k_db": round(20 * np.log10(ir(14e3)), 3), "image_db": round(20 * np.log10(ir(fs - 14e3) + 1e-12), 1)}
        res = {"matched_z_alpha": alpha, "matched_z_lag_us": round(lag * 1e6, 2), "matched_z_max_err_db_to_12.6k": round(dig_err, 3), "interp": interp, "tests": {}}
        n = int(0.25 * fs)
        tests = {
            "sine_1k_+3dB": (1e3, 10 ** (3 / 20)),
            "sine_9k_0dB": (9e3, 1.0),          # 9 kHz is +4.7 dB through the filter: clips although under FS
            "noise_+6dB": (None, 10 ** (6 / 20)),
        }
        for name, (fq, amp) in tests.items():
            t = np.arange(n) / fs
            if fq is None:
                b, a = signal.butter(4, 8e3, fs=fs)
                x = signal.lfilter(b, a, rng.standard_normal(n))
                x *= amp * 0.3 / np.max(np.abs(x[n // 4:]))
            else:
                x = amp * 0.3 * np.sin(TWO_PI * fq * t)
            fsc = 0.3                                                    # full scale = the sine peak that just fits
            y, D, vk, ek, tk = run(x, fs, fsc, sos, tau, g0, lag=lag)
            # The check: the input filter applied to the output, sampled at the device instants, against the clipped
            # A/D samples delayed by D. Skip the start (filter settling).
            _, _, vy, _, _ = run(y, fs, 1e9, sos, tau, g0, t0=D / fs)
            ref = np.clip(vk, -fsc, fsc * 32767 / 32768)[: len(vy)]
            m = tk[: len(vy)] > 0.05
            ek, vk = ek[: len(vy)], vk[: len(vy)]
            resid = vy[m] - ref[m]
            err_db = 10 * np.log10(np.sum(resid ** 2) / np.sum(ek[np.nonzero(m)[0]] ** 2) + 1e-30)
            # What an exact implementation still leaves: the error's parts the taper removes (above 11.5 kHz).
            F = np.fft.rfftfreq(len(ek), 1 / FS_DEV)
            ek_t = np.fft.irfft(np.fft.rfft(ek) * taper(F), len(ek))
            ideal_db = 10 * np.log10(np.sum((ek - ek_t)[m] ** 2) / np.sum(ek[m] ** 2) + 1e-30)
            peak_db = 20 * np.log10(np.max(np.abs(vy[m])) / fsc)
            clipped = float(np.mean(ek[m] != 0))
            res["tests"][name] = {"clipped_fraction": round(clipped, 4),
                                  "residual_re_clip_error_db": round(float(err_db), 1),
                                  "taper_floor_db": round(float(ideal_db), 1),
                                  "ad_peak_after_db_re_fs": round(float(peak_db), 2),
                                  "ad_peak_before_db_re_fs": round(float(20 * np.log10(np.max(np.abs(vk[m])) / fsc)), 2)}
        res["dry_delay_samples"] = D
        out["rates"][f"{fs:.0f}"] = res
        print(fs, json.dumps(res))
    np.savetxt(os.path.join(ROOT, "analysis", "analog", "adc_clip_g0.txt"), np.c_[tau[::20], g0[::20]], fmt="%.6e",
               header="tau_s g0 (every 10 us; the plugin computes its own table)")
    with open(os.path.join(ROOT, "analysis", "analog", "adc_clip.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
