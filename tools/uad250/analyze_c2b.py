"""C2b: carrier vs sideband power of the long-tone renders (modulation index without phase unwrapping).

    .venv/bin/python tools/uad250/analyze_c2b.py

A steady tone through a static (LTI) reverb returns one steady tone. The complex baseband around the tone is analysed:
carrier = power at |nu| < 0.05 Hz (what a static reverb would give), sidebands = the rest. Reports the carrier fraction,
the modulation bandwidth (nu below which 50 / 90 / 99 % of the sideband power lies) and the strongest sideband lines.
For delay modulation of depth d (seconds) the FM index is 2*pi*f*d, so the sideband fraction grows with f until it saturates.
"""
import os
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt, hilbert, decimate
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "C_mod")
r = 48000
for n in ["tone250_60s_T2.5_noise0", "tone1000_60s_T2.5_noise0", "tone4000_60s_T2.5_noise0",
          "tone1000_60s_T0.4_noise1", "tone1000_60s_T4.5_noise1"]:
    f = float(n.split("tone")[1].split("_")[0]); y, _ = sf.read(os.path.join(CAP, n + ".wav"))
    print(f"\n== {n}")
    for c, ch in enumerate("LR"):
        x = y[5 * r:55 * r, c].astype(float)
        z = hilbert(sosfiltfilt(butter(4, [f - 100, f + 100], "band", fs=r, output="sos"), x))
        z = z * np.exp(-2j * np.pi * f * np.arange(len(z)) / r); z = z[r:-r]
        z = decimate(z.real, 96, zero_phase=True) + 1j * decimate(z.imag, 96, zero_phase=True); fs = r / 96  # 500 Hz
        Z = np.abs(np.fft.fft(z * np.hanning(len(z)))) ** 2; nu = np.fft.fftfreq(len(z), 1 / fs)
        tot = Z.sum()
        # the carrier of a static reverb may sit a few mHz off (clock ratio): take the strongest bin within +-0.1 Hz
        m0 = np.abs(nu) < 0.1; k0 = np.argmax(np.where(m0, Z, 0)); w = np.abs(nu - nu[k0]) < 0.05
        car = Z[w].sum() / tot; sb = np.where(w, 0, Z); sbt = sb.sum()
        order = np.argsort(np.abs(nu)); cum = np.cumsum(sb[order]) / sbt
        bw = [float(np.abs(nu[order])[np.searchsorted(cum, q)]) for q in (0.5, 0.9, 0.99)]
        top = np.argsort(-sb)[:5]
        print(f"  {ch}: carrier {10 * np.log10(car):6.1f} dB of total ({car * 100:5.1f} %), carrier offset {nu[k0] * 1e3:6.1f} mHz; sideband power 50/90/99 % within |nu| < "
              f"{bw[0]:.2f}/{bw[1]:.2f}/{bw[2]:.2f} Hz; strongest lines (Hz, dB re carrier): " +
              ", ".join(f"{nu[i]:+.2f} ({10 * np.log10(sb[i] / Z[w].sum() + 1e-30):.0f})" for i in top))
