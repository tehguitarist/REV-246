"""B1/B2/B3: level behaviour of the UAD EMT 250 (docs/UAD250_test_plan.md §3 B). Reference level: 0 VU = -18 dBFS RMS.

    .venv/bin/python tools/uad250/analyze_b1.py

B1: each impulse-ladder IR is divided by its input gain and compared with the -30 dB case (inside the linear range, above
the floor); the residual energy in time windows (re the reference) shows non-linearity, the tail floor shows quantisation.
B2: harmonics 2-5 of the sine ladder (re the fundamental), by level. B3: idle/floor noise and the Power-off/Mix-0 paths.
"""
import glob, os
import numpy as np, soundfile as sf
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "B_level")
LAT, AT = 523, 100
def rd(n): return sf.read(os.path.join(CAP, n + ".wav"))
def e(x): return float((x ** 2).sum())
def dbr(a, b): return 10 * np.log10((a + 1e-40) / (b + 1e-40))

print("B1 impulse ladder: residual vs the -30 dB IR scaled to equal input (dB re reference energy in the window; L+R)")
W = [(0, .03), (.03, .3), (.3, 1), (1, 2), (2, 3)]
for T in ("0.4", "4.5"):
    ref, r = rd(f"imp_-30dB_T{T}"); ref = ref * 10 ** (30 / 20)
    print(f"\nTime {T} s   windows " + "  ".join(f"{a}-{b}s" for a, b in W) + "   peak gain re -30 dB (dB)   floor dBFS(last 1s)")
    for lv in [12, 6, 0, -6, -12, -18, -24, -30, -36, -42, -48, -54, -60, -66, -72, -78, -84, -90]:
        y, _ = rd(f"imp_{lv:+d}dB_T{T}"); g = 10 ** (-lv / 20); y = y * g
        res = []
        for a, b in W:
            s = slice(int((a) * r) + AT + LAT, int(b * r) + AT + LAT)
            res.append(dbr(e(y[s] - ref[s]), e(ref[s])))
        pk = dbr(abs(y).max() ** 2, abs(ref).max() ** 2)
        fl = 10 * np.log10(np.mean(y[-r:] ** 2) / g ** 2 + 1e-40)
        print(f"{lv:+4d} dBFS  " + "  ".join(f"{v:7.1f}" for v in res) + f"   {pk:6.2f}   {fl:7.1f}")
print("\nB2 sine ladder (Time 0.4, Noise off): output level, harmonics 2-5 re fundamental (dB), 3 s..7 s steady state")
for f in (250, 1000, 4000):
    print(f"  {f} Hz")
    for nm, lv in [("-60dB", -60), ("-40dB", -40), ("-20dB", -20), ("0VU", -15), ("-6dB", -6), ("+6VU", -9), ("+12VU", -3), ("0dB", 0), ("+6dB", 6)]:
        n = f"sine{f}_{lv:+d}dB" if nm[0] in "-+" and "VU" not in nm else f"cal_sine{f}_{nm}"
        if not os.path.exists(os.path.join(CAP, n + ".wav")): continue
        y, r = rd(n); x = y[3 * r:7 * r, 0]; N = len(x); X = np.abs(np.fft.rfft(x * np.hanning(N))) ** 2
        fb = lambda h: X[int(h * f * N / r) - 3:int(h * f * N / r) + 4].sum()
        print(f"   in {lv:+3d} dBFS pk: out {10 * np.log10((x ** 2).mean() + 1e-40):7.1f} dBFS rms; H2..H5 " + " ".join(f"{dbr(fb(h), fb(1)):6.1f}" for h in (2, 3, 4, 5)))
print("\nB3 floor / bypass paths")
for n in sorted(glob.glob(os.path.join(CAP, "silence_*.wav")) + glob.glob(os.path.join(CAP, "whitenoise_*.wav"))):
    y, r = sf.read(n); b = os.path.basename(n)[:-4]
    x = y[r:, 0]; S = np.abs(np.fft.rfft(x * np.hanning(len(x)))) ** 2
    fr = np.fft.rfftfreq(len(x), 1 / r)
    bands = [(20, 200), (200, 2e3), (2e3, 8e3), (8e3, 12e3), (12e3, 20e3)]
    print(f"{b:42s} rms {10 * np.log10((x ** 2).mean() + 1e-40):7.1f} dBFS  band share(dB re total) " +
          " ".join(f"{dbr(S[(fr >= a) & (fr < c)].sum(), S.sum()):6.1f}" for a, c in bands))
