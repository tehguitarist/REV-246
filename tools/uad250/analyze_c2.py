"""C2: modulation spectrum of the 250 algorithm from the long-tone renders (docs/UAD250_test_plan.md §3 C2).

    .venv/bin/python tools/uad250/analyze_c2.py [case ...]

Complex-demodulates the steady-state output at the tone frequency: amplitude (AM, gain modulation) and phase (PM/FM,
delay modulation) versus time, their spectra (modulation rate), depth (rms), and for FM how depth scales with f
(delay modulation: phase deviation proportional to f). Uses 5-55 s of the 60 s files.
"""
import os, sys
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt, hilbert, find_peaks
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "C_mod")
r = 48000

def demod(x, f):
    sos = butter(4, [f - 40, f + 40], "band", fs=r, output="sos")
    z = hilbert(sosfiltfilt(sos, x)); ph = np.unwrap(np.angle(z * np.exp(-2j * np.pi * f * np.arange(len(z)) / r)))
    return np.abs(z), ph

def spec_peaks(v, fs, n=6, fmax=50):
    v = v - v.mean(); V = np.abs(np.fft.rfft(v * np.hanning(len(v)))) / len(v) * 4; fr = np.fft.rfftfreq(len(v), 1 / fs)
    m = (fr > 0.02) & (fr < fmax); pk, _ = find_peaks(V * m, distance=5)
    top = pk[np.argsort(-V[pk])[:n]]
    return [(round(float(fr[i]), 3), float(V[i])) for i in top], V, fr

names = sys.argv[1:] or ["tone250_60s_T2.5_noise0", "tone1000_60s_T2.5_noise0", "tone4000_60s_T2.5_noise0",
                         "tone1000_60s_T0.4_noise1", "tone1000_60s_T4.5_noise1", "tone1000_60s_T2.5_0VU_noise1"]
for n in names:
    p = os.path.join(CAP, n + ".wav")
    if not os.path.exists(p): continue
    f = float(n.split("tone")[1].split("_")[0]); y, _ = sf.read(p)
    print(f"\n== {n}  ({f:.0f} Hz)")
    for c, ch in enumerate("LR"):
        a, ph = demod(y[5 * r:55 * r, c].astype(float), f)
        a, ph = a[r:-r], ph[r:-r]
        d = 10  # decimate to ~4.8 kHz
        a, ph = a[::d * 10], ph[::d * 10]; fs = r / (d * 10)
        am = a / a.mean()
        print(f"  {ch}: AM rms {np.std(am) * 100:6.2f} % ({20 * np.log10(am.max() / am.min()):5.2f} dB p-p) | PM rms {np.std(ph):7.4f} rad "
              f"(= {np.std(ph) / (2 * np.pi * f) * 1e6:7.2f} us of delay) p-p {ph.max() - ph.min():.3f} rad")
        pa, _, _ = spec_peaks(am, fs); pp, _, _ = spec_peaks(ph, fs)
        print("     AM peaks Hz (rel amp): " + ", ".join(f"{fr:.3f} ({v:.3f})" for fr, v in pa))
        print("     PM peaks Hz (rad):     " + ", ".join(f"{fr:.3f} ({v:.4f})" for fr, v in pp))
