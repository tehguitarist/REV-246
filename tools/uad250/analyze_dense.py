"""Delay-modulation time series from the dense onset sampling (C_dense, Time 0.4 s, 0.25 s steps over 240 s).

    .venv/bin/python tools/uad250/analyze_dense.py

For each onset snapshot (R channel by default) and each lag window of the response, the shift d_i(window) against an
iteratively refined template is found by cross-correlation of the 60-300 Hz band (coherence there survives shifts of
about 1 ms), with a parabolic sub-sample peak. Reports the spread of d, its autocorrelation, spectrum (period of the
modulation) and the correlation between lag windows (common-mode versus independent per loop).
"""
import os, sys
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "C_dense")
r = 48000; FDEV = 19.2e6 / 728; OFF = 11  # ms from onset to first arrival, R channel (latency 523 smp)

def load(prefix, n_per, spacing, offs, ch):
    H, T = [], []
    for k, o in enumerate(offs):
        y, _ = sf.read(os.path.join(CAP, prefix.format(o)))
        for j in range(n_per):
            t = 0.1 + o / 1000 + j * spacing; a = int(round(t * r)); H.append(y[a:a + int(0.6 * r), ch]); T.append(t)
    i = np.argsort(T); return np.array(H)[i], np.array(T)[i]

def shifts(X, a_ms, b_ms, lo=60, hi=300, maxlag_ms=3.0, iters=3):
    sos = butter(4, [lo, hi], "band", fs=r, output="sos"); F = sosfiltfilt(sos, X, axis=1)
    pad = int(maxlag_ms * 1e-3 * r) + 2; s0 = int((a_ms + OFF) * r / 1000); s1 = int((b_ms + OFF) * r / 1000)
    ref = F[0, s0 - pad:s1 + pad].copy(); w = np.hanning(s1 - s0); d = np.zeros(len(F))
    for _ in range(iters):
        sh = []
        for i in range(len(F)):
            seg = F[i, s0:s1] * w; best = None
            c = np.array([(seg * ref[pad + k:pad + k + (s1 - s0)]).sum() for k in range(-pad + 1, pad - 1)])
            m = int(np.argmax(c)); k = m - (pad - 1)
            if 0 < m < len(c) - 1:
                den = c[m - 1] - 2 * c[m] + c[m + 1]; k += 0.5 * (c[m - 1] - c[m + 1]) / den if den != 0 else 0
            sh.append(-k)  # shift of snapshot relative to ref (samples at 48k)
        d = np.array(sh) / r * 1e3  # ms
        # refine template: mean of snapshots aligned by the estimated shifts (integer rounding)
        al = np.mean([np.roll(F[i], int(round(d[i] * r / 1e3)))[s0 - pad:s1 + pad] for i in range(len(F))], axis=0)
        ref = al
    return d

def report(name, X, T, step):
    print(f"\n=== {name}: {len(X)} snapshots, step {step} s")
    wins = [(40, 80), (80, 160), (160, 260)]; D = {}
    for a, b in wins:
        d = shifts(X, a, b); d = d - np.median(d); D[(a, b)] = d
        ac = np.correlate(d - d.mean(), d - d.mean(), "full")[len(d) - 1:] / ((d - d.mean()) ** 2).sum()
        sp = np.abs(np.fft.rfft((d - d.mean()) * np.hanning(len(d)))); fr = np.fft.rfftfreq(len(d), step)
        top = np.argsort(-sp[1:])[:4] + 1
        print(f" lag {a:3d}-{b:3d} ms: shift rms {d.std():6.3f} ms = {d.std() * FDEV / 1e3:5.1f} dev smp (p-p {d.max() - d.min():.2f} ms); "
              f"autocorr at 0.25/0.5/1/2/5 s: " + " ".join(f"{ac[int(t / step)]:5.2f}" for t in (0.25, 0.5, 1, 2, 5) if int(t / step) < len(ac)) +
              f" | spectrum peaks Hz: " + ", ".join(f"{fr[i]:.3f}" for i in top))
    ks = list(D)
    print(" correlation of shift series between lag windows: " + ", ".join(f"{a}-{b}ms~{c}-{e}ms {np.corrcoef(D[(a, b)], D[(c, e)])[0, 1]:+.2f}" for i, (a, b) in enumerate(ks) for (c, e) in ks[i + 1:]))
    return D

for ch, nm in ((1, "R"), (0, "L")):
    X, T = load("dense_T0.4_off{:03d}ms.wav", 240, 1.0, [0, 250, 500, 750], ch)
    D = report(f"T0.4 dense 0.25 s, {nm}", X, T, 0.25)
    if nm == "R": np.save(os.path.join(ROOT, "analysis", "uad250", "dense_shift_R.npy"), np.array([D[k] for k in D]))


# ---- harmonic fit of the shift series: frequency, amplitude (ms and device samples), phase, harmonics, residual
from scipy.optimize import least_squares
def fit(t, d, f0):
    def model(p, nh=3):
        f = p[0]; out = np.full_like(t, p[1])
        for h in range(1, nh + 1): out += p[2 * h] * np.sin(2 * np.pi * h * f * t) + p[2 * h + 1] * np.cos(2 * np.pi * h * f * t)
        return out
    fs = np.linspace(f0 - 0.03, f0 + 0.03, 241)
    best = min(fs, key=lambda f: ((d - np.polyval([0], t) - model(np.r_[f, 0, 1, 0, 0, 0, 0, 0], 3) * 0).sum() * 0 + np.linalg.lstsq(np.c_[np.ones_like(t), np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)], d, rcond=None)[1].sum()
                                  if True else 0))
    p0 = np.r_[best, 0, np.zeros(6)]
    A = np.c_[np.ones_like(t), np.sin(2 * np.pi * best * t), np.cos(2 * np.pi * best * t)]; c = np.linalg.lstsq(A, d, rcond=None)[0]
    p0[1:4] = c
    sol = least_squares(lambda p: model(p) - d, p0).x
    return sol, model(sol)

print("\n=== sinusoid fit of the shift series (3 harmonics), fundamental near 0.55 Hz")
for ch, nm in ((1, "R"), (0, "L")):
    X, T = load("dense_T0.4_off{:03d}ms.wav", 240, 1.0, [0, 250, 500, 750], ch)
    for a, b in [(40, 80), (80, 160), (160, 260)]:
        d = shifts(X, a, b); d = d - np.median(d); sol, m = fit(T, d, 0.55)
        amp = [np.hypot(sol[2 * h], sol[2 * h + 1]) for h in (1, 2, 3)]
        ph = np.degrees(np.arctan2(sol[3], sol[2]))
        res = 1 - ((d - m) ** 2).sum() / ((d - d.mean()) ** 2).sum()
        print(f" {nm} {a:3d}-{b:3d} ms: f = {sol[0]:.4f} Hz (period {1 / sol[0]:.3f} s); amplitude H1/H2/H3 = {amp[0]:.3f}/{amp[1]:.3f}/{amp[2]:.3f} ms "
              f"(H1 = {amp[0] * FDEV / 1e3:.1f} dev smp peak); phase {ph:6.1f} deg; variance explained {res * 100:.0f} %")


# ---- CLEAN-style line extraction (up to the 2 Hz Nyquist of the 0.25 s sampling)
def clean(t, d, n=6, fmax=2.0):
    d = d - d.mean(); out = []; res = d.copy(); fs = np.arange(0.01, fmax, 0.0005)
    for _ in range(n):
        best = None
        for f in fs:
            A = np.c_[np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)]; c = np.linalg.lstsq(A, res, rcond=None)[0]; e = ((res - A @ c) ** 2).sum()
            if best is None or e < best[0]: best = (e, f, c)
        e, f, c = best; res = res - np.c_[np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)] @ c
        out.append((f, np.hypot(*c), np.degrees(np.arctan2(c[1], c[0]))))
    return out, 1 - (res ** 2).sum() / (d ** 2).sum()

print("\n=== six strongest lines of the shift series (Hz, amplitude ms, phase deg)")
for ch, nm in ((1, "R"), (0, "L")):
    X, T = load("dense_T0.4_off{:03d}ms.wav", 240, 1.0, [0, 250, 500, 750], ch)
    for a, b in [(40, 80), (80, 160), (160, 260)]:
        d = shifts(X, a, b); d -= np.median(d); o, v = clean(T, d)
        print(f" {nm} {a:3d}-{b:3d} ms: {v * 100:3.0f}% explained: " + "; ".join(f"{f:.3f}, {A:.3f}, {p:.0f}" for f, A, p in o))


# ---- coherence versus onset separation, per octave band: a sinusoidal delay of amplitude A gives rho(dt=T/2) ~ J0(4 pi f A)
from scipy.signal import butter as _b, sosfiltfilt as _sf
from scipy.special import j0
from scipy.optimize import minimize_scalar
print("\n=== pair coherence rho(f) versus onset separation (R, lag window 80-260 ms; Time 0.4 s)")
X, T = load("dense_T0.4_off{:03d}ms.wav", 240, 1.0, [0, 250, 500, 750], 1)
s0, s1 = int((80 + OFF) * r / 1000), int((260 + OFF) * r / 1000); w = np.hanning(s1 - s0)
bands = [125, 250, 500, 1000, 2000]; seps = [0.25, 0.5, 0.75, 0.9, 1.0, 1.25, 1.5, 1.8, 2.0, 2.5, 3.6, 5.4, 9.0]
curves = {}
for f in bands:
    F = _sf(_b(4, [f / 1.414, f * 1.414], "band", fs=r, output="sos"), X[:, :int(0.4 * r)], axis=1)[:, s0:s1] * w
    row = []
    for dt in seps:
        k = int(round(dt / 0.25)); num = sum((F[i] * F[i + k]).sum() for i in range(len(F) - k)); den = sum((F[i] ** 2).sum() for i in range(len(F) - k)); row.append(num / den)
    curves[f] = np.array(row)
print("sep s   " + " ".join(f"{d:5.2f}" for d in seps))
for f in bands: print(f"{f:5d}Hz " + " ".join(f"{v:5.2f}" for v in curves[f]))
k = seps.index(0.9)
print("anti-phase (0.9 s) fit of A with rho = a + (1-a) J0(4 pi f A) using the in-phase (1.8 s) value as 'a'-free check:")
for f in bands:
    rho_anti, rho_in = curves[f][k], curves[f][seps.index(1.8)]
    res = minimize_scalar(lambda A: (j0(4 * np.pi * f * A) - (rho_anti)) ** 2, bounds=(1e-6, 4e-4 if f >= 250 else 2e-3), method="bounded")
    print(f" {f:5d} Hz: rho(0.9 s) {rho_anti:5.2f}, rho(1.8 s) {rho_in:5.2f} -> A (J0 model) {res.x * 1e3:.3f} ms = {res.x * FDEV:.1f} device samples")
