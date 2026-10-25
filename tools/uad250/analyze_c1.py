"""C1: modulation snapshots from the impulse trains (docs/UAD250_test_plan.md §3 C).

    .venv/bin/python tools/uad250/analyze_c1.py [T]

Cuts each train into per-impulse IR snapshots (Noise on) and compares them with each other and with the Noise-off IR:
windowed, octave-band normalised correlation versus time since onset (the capture-to-capture decorrelation curve of
ml_structure.md §4-5), plus exact-repeat detection (does the random state repeat within the train?).
"""
import os, sys
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "C_mod")
SP, LEAD, N_ON = 12.0, 0.1, 24
BANDS = [125, 250, 500, 1000, 2000, 4000, 8000]
TIMES = [0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]

def snaps(name, n, extra=0):
    y, r = sf.read(os.path.join(CAP, name + ".wav")); out = []
    for k in range(n):
        a = int(round((LEAD + k * SP) * r)) + k * extra
        out.append(y[a:a + int((SP - 0.5) * r)])
    return out, r

def rho(a, b, r, t, w=0.1):
    s = slice(int(t * r), int((t + w) * r)); x, z = a[s], b[s]
    return float((x * z).sum() / np.sqrt((x ** 2).sum() * (z ** 2).sum() + 1e-40))

def table(title, A, B_list, r):
    print(f"\n{title}  (rho in 100 ms windows starting at t after onset; L channel; mean over pairs)")
    print("band  " + " ".join(f"{t:5.2f}" for t in TIMES))
    for f in BANDS:
        sos = butter(4, [f / 1.414, min(f * 1.414, 0.45 * r)], "band", fs=r, output="sos")
        fa = sosfiltfilt(sos, A[:, 0])
        row = []
        for t in TIMES:
            row.append(np.mean([rho(fa, sosfiltfilt(sos, b[:, 0]), r, t) for b in B_list]))
        print(f"{f:5d} " + " ".join(f"{v:5.2f}" for v in row))

for T in (sys.argv[1:] or ["1.0", "2.5", "4.5"]):
    n1, n0 = f"train12s_T{T}_noise1", f"train12s_T{T}_noise0"
    if not os.path.exists(os.path.join(CAP, n1 + ".wav")): continue
    S1, r = snaps(n1, N_ON); S0, _ = snaps(n0, 4)
    print(f"\n===== Time {T} s =====")
    d0 = [10 * np.log10(((S0[0] - s) ** 2).sum() / (S0[0] ** 2).sum() + 1e-30) for s in S0[1:]]
    print("Noise-off snapshots vs snapshot 0 (time invariance; dB):", np.round(d0, 1))
    pairs = [(i, j) for i in range(N_ON) for j in range(i + 1, N_ON)]
    nul = np.array([10 * np.log10(((S1[i] - S1[j]) ** 2).sum() / (S1[i] ** 2).sum() + 1e-30) for i, j in pairs])
    print(f"Noise-on snapshot pairs: best null {nul.min():.1f} dB (pair {pairs[int(nul.argmin())]}), median {np.median(nul):.1f} dB -> "
          + ("a pair REPEATS exactly" if nul.min() < -120 else "no exact repeats within the train"))
    table("Noise-on snapshot k vs snapshot 0", S1[0], S1[1:], r)
    table("Noise-on snapshots vs the Noise-off (static) IR", S0[0], S1, r)
