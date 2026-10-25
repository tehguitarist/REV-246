"""Per-loop delay tracking from long white-noise renders (C_loops; docs/UAD250_test_plan.md §2.5).

    .venv/bin/python tools/uad250/analyze_loops.py [case ...]

The output of a feedback comb correlates with itself at its loop lag. The global autocorrelation of the (band-limited)
output locates the loops; then for each loop the lag of the correlation peak is tracked in short sliding windows
(0.25 s, hop 0.05 s, sub-sample parabolic peak), giving that loop's delay as a function of absolute time. The delay series
is analysed for lines (CLEAN), depth and, between the two seeds at Time 4.5, for reproducibility (a delay that is a function of
absolute time must agree for any input).
"""
import os, sys
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt, find_peaks
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "C_loops")
OUT = os.path.join(ROOT, "analysis", "uad250", "loops"); os.makedirs(OUT, exist_ok=True)
r = 48000; FDEV = 19.2e6 / 728; W = 24000; HOP = 2400; RANGE = 60; BAND = tuple(float(v) for v in os.environ.get('LOOP_BAND', '100,400').split(','))

def prep(path, ch):
    y, _ = sf.read(path); y = y[:, ch].astype(np.float64)[5 * r:]
    return sosfiltfilt(butter(4, list(BAND), "band", fs=r, output="sos"), y)

def global_acf(y, lo=0.10, hi=0.45):
    N = 1 << int(np.ceil(np.log2(2 * len(y)))); S = np.abs(np.fft.rfft(y, N)) ** 2
    ac = np.fft.irfft(S, N)[:int(hi * r)]; return ac / ac[0]

def track(y, L, rng=RANGE):
    """rho(lag) in sliding windows for lags L-rng..L+rng; returns lags (samples, sub-sample) and peak rho per window."""
    n = len(y); e = np.concatenate([[0], np.cumsum(y * y)]); starts = np.arange(rng + 1 + int(L) + 1, n - W, HOP)
    # windows start at s; compare y[s:s+W] with y[s-lag:s-lag+W]
    rho = np.zeros((2 * rng + 1, len(starts)))
    for k, m in enumerate(range(-rng, rng + 1)):
        lag = int(L) + m; p = y[lag:] * y[:n - lag]  # p[i] = y[i+lag]*y[i]
        cs = np.concatenate([[0], np.cumsum(p)])
        a = starts - lag  # window of the early copy y[a:a+W] pairs with y[a+lag:...]
        c = cs[a + W] - cs[a]; e1 = e[a + W] - e[a]; e2 = e[a + lag + W] - e[a + lag]
        rho[k] = c / np.sqrt(e1 * e2 + 1e-30)
    # peak position per window (parabolic). The shifts (about +-20 samples) exceed the width of a high-band peak, so the linear
    # matched-filter estimate fails there; the 100-400 Hz peak is broad and its position is unambiguous within +-RANGE.
    i = rho.argmax(0); ii = np.clip(i, 1, 2 * rng - 1); cols = np.arange(rho.shape[1])
    a_, b_, c_ = rho[ii - 1, cols], rho[ii, cols], rho[ii + 1, cols]
    d = 0.5 * (a_ - c_) / (a_ - 2 * b_ + c_ + 1e-30)
    return int(L) - rng + ii + d, rho.max(0), (starts + W / 2) / r + 5.0

def clean(t, d, n=6, fmax=3.0):
    d = d - d.mean(); res = d.copy(); out = []; fs = np.arange(0.02, fmax, 0.001)
    for _ in range(n):
        best = None
        for f in fs:
            A = np.c_[np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)]; c = np.linalg.lstsq(A, res, rcond=None)[0]; e = ((res - A @ c) ** 2).sum()
            if best is None or e < best[0]: best = (e, f, c)
        _, f, c = best; res = res - np.c_[np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)] @ c
        out.append((f, np.hypot(*c), np.degrees(np.arctan2(c[1], c[0]))))
    return out, 1 - (res ** 2).sum() / (d ** 2).sum()

def main():
    names = sys.argv[1:] or ["noise240s_T4.5_seed11", "noise240s_T4.5_seed12", "noise240s_T2.5_seed11", "noise240s_T1.5_seed11"]
    series = {}
    for n in names:
        for ch, cn in ((1, "R"), (0, "L")):
            y = prep(os.path.join(CAP, n + ".wav"), ch); ac = global_acf(y)
            lo = int(0.12 * r); pk, _ = find_peaks(ac[lo:int(0.45 * r)], height=0.05, distance=200)
            pk = pk + lo; pk = pk[np.argsort(-ac[pk])][:8]
            print(f"\n== {n} {cn}: global ACF peaks 120-450 ms (ms, height): " + ", ".join(f"{p / r * 1e3:.2f} ({ac[p]:.2f})" for p in sorted(pk)))
            loops = [p for p in sorted(pk) if min(abs(p / r * 1e3 - q) for q in (167.42, 190.0, 210.25)) < 1.5]
            for L in loops:
                lag, rh, t = track(y, L)
                dly = (lag - np.median(lag)) / r * 1e3  # ms
                o, v = clean(t, dly)
                series[(n, cn, round(L / r * 1e3, 1))] = (t, dly)
                print(f"  loop {L / r * 1e3:7.2f} ms ({L / r * FDEV:6.1f} dev smp): window rho {np.median(rh):.2f}; delay rms {dly.std():.3f} ms = {dly.std() * FDEV / 1e3:4.1f} dev smp, p-p {np.ptp(dly):.2f} ms; "
                      f"lines ({v * 100:.0f}% expl.) Hz/ms(sinc-corrected)/deg: " + "; ".join(f"{f:.3f}/{A / abs(np.sinc(f * W / r)):.3f}/{p:.0f}" for f, A, p in o[:5]))
    # reproducibility: same loop, two seeds, same channel
    for cn in "RL":
        a = [k for k in series if k[0] == "noise240s_T4.5_seed11" and k[1] == cn]
        for k in a:
            k2 = ("noise240s_T4.5_seed12", cn, k[2])
            if k2 in series:
                t1, d1 = series[k]; t2, d2 = series[k2]; m = min(len(d1), len(d2))
                print(f"reproducibility {cn} loop {k[2]} ms, seed 11 vs 12: corr {np.corrcoef(d1[:m], d2[:m])[0, 1]:+.3f}, rms diff {np.std(d1[:m] - d2[:m]):.3f} ms (series rms {d1.std():.3f})")
    # cross-loop relations (T4.5 seed 11)
    for cn in "RL":
        ks = sorted(k for k in series if k[0] == "noise240s_T4.5_seed11" and k[1] == cn)
        for i, a in enumerate(ks):
            for b in ks[i + 1:]:
                m = min(len(series[a][1]), len(series[b][1])); print(f"cross-loop {cn} {a[2]} vs {b[2]} ms: corr {np.corrcoef(series[a][1][:m], series[b][1][:m])[0, 1]:+.2f}")
    np.savez(os.path.join(OUT, "series.npz"), **{"|".join(map(str, k)): np.c_[v[0], v[1]] for k, v in series.items()})

if __name__ == "__main__":
    main()
