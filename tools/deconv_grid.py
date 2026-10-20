"""Recover the EMT 246's internal (device-rate) impulse response from a 96 kHz capture.

Assumes: the unit runs at fs_dev (≈19.2 MHz/728 = 26373.6 Hz), the cluster's first tap is a
single-sample impulse, and the analogue output chain is identical on both channels apart from
a fractional-sample delay. Writes <out>.npy with shape (n, 2) at fs_dev, first cluster tap lands exactly on an integer index.
"""
import argparse, json
import numpy as np, soundfile as sf
from scipy.interpolate import CubicSpline

def locate_first_tap(s, fs, thresh_db=-40):
    thr = np.abs(s).max() * 10 ** (thresh_db / 20)
    return int(np.argmax(np.abs(s) > thr))

def recover(x, fs, fs_dev, kernel_ch=1, kernel_ms=5.2, pre_ms=0.1, f_hi=12300.0, f_roll=11000.0, eps_rel=1e-3):
    on = locate_first_tap(x[:, kernel_ch], fs)
    i0 = on - int(pre_ms / 1000 * fs)
    L = int(kernel_ms / 1000 * fs)
    g = x[i0:i0 + L, kernel_ch].copy()
    nt = int(0.0015 * fs); g[-nt:] *= np.hanning(2 * nt)[nt:]
    N = 1 << int(np.ceil(np.log2(len(x) + L)) + 1)
    f = np.fft.rfftfreq(N, 1 / fs)
    G = np.fft.rfft(g, N)
    eps = eps_rel * np.max(np.abs(G)) ** 2
    win = np.where(f < f_roll, 1.0, np.where(f < f_hi, 0.5 * (1 + np.cos(np.pi * (f - f_roll) / (f_hi - f_roll))), 0.0))
    # inter-channel fractional delay relative to kernel channel (estimated from first tap phase slope)
    out, offsets = [], []
    for c in range(x.shape[1]):
        Y = np.fft.rfft(x[:, c], N)
        tau = 0.0
        if c != kernel_ch:
            Lw = int(0.0019 * fs)
            A = np.fft.rfft(x[i0:i0 + Lw, c], 1 << 16); B = np.fft.rfft(x[i0:i0 + Lw, kernel_ch], 1 << 16)
            ff = np.fft.rfftfreq(1 << 16, 1 / fs); m = (ff > 300) & (ff < 9000)
            tau = -np.polyfit(ff[m], np.unwrap(np.angle(A[m] / B[m])), 1)[0] / (2 * np.pi)
            Y = Y * np.exp(2j * np.pi * f * tau)  # advance by tau
        offsets.append(tau)
        D = Y * np.conj(G) / (np.abs(G) ** 2 + eps) * win
        d = np.fft.irfft(D, N)[: len(x)]
        cs = CubicSpline(np.arange(len(d)), d)
        # deconvolution puts each tap exactly at the kernel start (i0), so sample on a grid through i0
        k0 = int(np.floor(i0 / fs * fs_dev))
        n = int((len(x) - 1 - i0) / fs * fs_dev) + k0
        out.append(cs(i0 + (np.arange(n) - k0) / fs_dev * fs))
    h = np.stack(out, 1)
    h /= np.abs(h).max()
    return h, i0 / fs, offsets

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("wav"); ap.add_argument("out")
    ap.add_argument("--fs-dev", type=float, default=19.2e6 / 728)
    a = ap.parse_args()
    x, fs = sf.read(a.wav, always_2d=True)
    h, t0, offs = recover(x, fs, a.fs_dev)
    np.save(a.out, h)
    print(json.dumps({"fs_dev": a.fs_dev, "kernel_start_s": t0, "channel_offsets_us": [o * 1e6 for o in offs], "n": len(h)}))
