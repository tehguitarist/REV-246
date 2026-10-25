"""Per-loop modulation line spectrum by direct model fit (C_loops; docs/UAD250_test_plan.md §2.5).

    LOOP_BAND=100,400 .venv/bin/python tools/uad250/analyze_loops2.py [case ...]

Peak picking per window is dominated by wrong-lag outliers, so instead: rho[m, t] (correlation of the output with itself at
lag L_j + m in 0.5 s windows) is computed once per loop; for a trial line (f, A, phase) the loop's delay is
L_j + A sin(2 pi f t + phase) and the objective is the summed rho at that lag over all windows. The line spectrum S_j(f) is
the best gain over (A, phase) relative to no modulation, in units of mean rho; its peaks name the loop's modulation lines, and
the fitted (A, phase) must agree between two seeds at the same Time (the delay depends on absolute time, not the input).
"""
import os, sys
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "C_loops")
r = 48000; FDEV = 19.2e6 / 728; W = 24000; HOP = 2400; RNG = 60
BAND = tuple(float(v) for v in os.environ.get("LOOP_BAND", "100,400").split(","))
LOOPS = {"167": 167.42, "190": 190.00, "210": 210.25}

def prep(path, ch):
    y, _ = sf.read(path); y = y[:, ch].astype(np.float64)[5 * r:]
    return sosfiltfilt(butter(4, list(BAND), "band", fs=r, output="sos"), y)

def rho_table(y, L):
    n = len(y); e = np.concatenate([[0], np.cumsum(y * y)]); starts = np.arange(RNG + 1 + int(L) + 1, n - W, HOP)
    rho = np.zeros((2 * RNG + 1, len(starts)))
    for k, m in enumerate(range(-RNG, RNG + 1)):
        lag = int(L) + m; cs = np.concatenate([[0], np.cumsum(y[lag:] * y[:n - lag])]); a = starts - lag
        rho[k] = (cs[a + W] - cs[a]) / np.sqrt((e[a + W] - e[a]) * (e[a + lag + W] - e[a + lag]) + 1e-30)
    return rho, (starts + W / 2) / r + 5.0

def objective(rho, off):
    """sum over windows of rho at fractional lag index (RNG + off[t])."""
    x = np.clip(RNG + off, 0, 2 * RNG - 1e-9); i = x.astype(int); fr = x - i; c = np.arange(rho.shape[1])
    return (rho[i, c] * (1 - fr) + rho[i + 1, c] * fr).sum()

def spectrum(rho, t, fs, amps, nph=12):
    j0 = objective(rho, np.zeros(len(t))); out = []
    for f in fs:
        best = (-1e9, 0, 0)
        for A in amps:
            for ph in np.linspace(0, 2 * np.pi, nph, endpoint=False):
                j = objective(rho, A * np.sin(2 * np.pi * f * t + ph))
                if j > best[0]: best = (j, A, ph)
        out.append(((best[0] - j0) / rho.shape[1], best[1], best[2]))
    return np.array(out)

def main():
    names = sys.argv[1:] or ["noise240s_T4.5_seed11", "noise240s_T4.5_seed12"]
    fs = np.arange(0.30, 1.30, 0.01); amps = np.array([4, 8, 12, 16, 20, 26]) * 1.0  # samples at 48 kHz (0.083 .. 0.54 ms)
    res = {}
    for n in names:
        for ch, cn in ((1, "R"), (0, "L")):
            y = prep(os.path.join(CAP, n + ".wav"), ch)
            for lk, Lms in LOOPS.items():
                rho, t = rho_table(y, Lms * r / 1000)
                S = spectrum(rho, t, fs, amps); res[(n, cn, lk)] = (S, rho.max(0).mean())
                top = np.argsort(-S[:, 0])[:5]
                print(f"{n} {cn} loop {lk}: mean peak rho {rho.max(0).mean():.2f}; strongest lines Hz (gain, A ms, phase deg): " +
                      "; ".join(f"{fs[i]:.2f} ({S[i, 0]:.4f}, {S[i, 1] / r * 1e3:.2f}, {np.degrees(S[i, 2]):.0f})" for i in sorted(top, key=lambda i: fs[i])), flush=True)
    np.savez(os.path.join(ROOT, "analysis", "uad250", "loops", f"spectra_{int(BAND[0])}_{int(BAND[1])}.npz"), fs=fs,
             **{"|".join(k): v[0] for k, v in res.items()})
    print("\nreproducibility of the line spectra (corr over f of S_j(f), seed 11 vs 12):")
    for cn in "RL":
        for lk in LOOPS:
            a, b = ("noise240s_T4.5_seed11", cn, lk), ("noise240s_T4.5_seed12", cn, lk)
            if a in res and b in res: print(f"  {cn} loop {lk}: corr {np.corrcoef(res[a][0][:, 0], res[b][0][:, 0])[0, 1]:+.2f}")

if __name__ == "__main__":
    main()


def refine(names=("noise240s_T4.5_seed11", "noise240s_T4.5_seed12")):
    """Significance of each loop's lines against the off-line floor, then a fine (f, A, phase) fit of the two main lines."""
    z = np.load(os.path.join(ROOT, "analysis", "uad250", "loops", f"spectra_{int(BAND[0])}_{int(BAND[1])}.npz")); fs = z["fs"]
    print(f"\nband {BAND}: peak significance (gain / median over f, and gain at 0.74 and 0.55 Hz re median), seeds 11 | 12")
    for cn in "RL":
        for lk in LOOPS:
            row = []
            for n in names:
                S = z["|".join((n, cn, lk))][:, 0]; med = np.median(S); i74 = np.argmin(abs(fs - 0.74)); i55 = np.argmin(abs(fs - 0.555))
                row.append(f"top {fs[np.argmax(S)]:.2f} Hz ({S.max() / med:4.1f}x med), 0.74 Hz {S[i74 - 1:i74 + 2].max() / med:4.1f}x, 0.55 Hz {S[i55 - 2:i55 + 3].max() / med:4.1f}x")
            print(f"  {cn} loop {lk}: " + " | ".join(row))
    print("\nfine fit of the two main lines (joint amplitude/phase, each f searched +-0.012 Hz):")
    f1s = np.arange(0.540, 0.572, 0.004); f2s = np.arange(0.728, 0.760, 0.004); A_ = np.arange(1, 25, 1.5); P_ = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    for n in names:
        for ch, cn in ((1, "R"), (0, "L")):
            y = prep(os.path.join(CAP, n + ".wav"), ch)
            for lk, Lms in LOOPS.items():
                rho, t = rho_table(y, Lms * r / 1000); j0 = objective(rho, np.zeros(len(t))); out = []
                for fgrid in (f1s, f2s):
                    best = (-1e9, 0, 0, 0)
                    for f in fgrid:
                        for A in A_:
                            for ph in P_:
                                j = objective(rho, A * np.sin(2 * np.pi * f * t + ph))
                                if j > best[0]: best = (j, f, A, ph)
                    out.append(best)
                print(f"  {n[-6:]} {cn} loop {lk}: " + " | ".join(f"{b[1]:.3f} Hz A {b[2] / r * 1e3:.2f} ms = {b[2] / r * FDEV:.1f} dev smp ph {np.degrees(b[3]):3.0f} (gain {(b[0] - j0) / rho.shape[1]:.4f})" for b in out), flush=True)

if os.environ.get("REFINE"):
    refine()
