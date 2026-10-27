"""Same front section => same device output in t < 25 ms (causality). Compare captures of different SPAN (or LATE) at identical EARLY, in the front frame
(t = 0 at the front-frame zero; fractional shifts removed). dY = data difference, dM = model difference; both after a common per-capture gain (bench)."""
import sys
from lb_common import *
B = bench_all()
bA, bB = bands(20, 100), bands(100, 1000)
def frame(y, fr, n=3000, pre=50):
    N = 8192; f = np.fft.rfftfreq(N)
    Y = np.fft.rfft(np.vstack([y[:n], np.zeros((N - n, 2))]), axis=0)
    Y = Y * np.exp(2j * np.pi * f * fr)[:, None]          # y(t + fr)
    return np.fft.irfft(Y, N, axis=0)[:n]
def getr(key):
    r = B[key]; y, m, _ = load(key, r, n=3000); return r, frame(y, r["front_shift"]), frame(m, r["front_shift"])
def rep(ref, keys, label):
    rr, yr, mr = getr(ref)
    pk = np.abs(mr[K0 + 0: K0 + 800, 0]).max()
    ib = K0 + int(0.025 * FS); ia = 10
    q = lambda v: 20 * np.log10(np.sqrt(np.mean(v[ia:ib] ** 2)) / pk + 1e-12)
    print(f"\n{label}: ref {ref} fs {rr['front_shift']:.2f};  rms dB re 410-peak over front-frame t in [-40, 25 ms]:   full | 20-100 | 100-1k    (dY data / dM model), ch L then R")
    for k in keys:
        if k not in B: continue
        r, y, m = getr(k)
        line = f"  {str(k):20s} fs {r['front_shift']:6.2f} o {r['o']:7.2f} |"
        for c in (0, 1):
            dy, dm = y[:, c] - yr[:, c], m[:, c] - mr[:, c]
            line += " ch%d: " % c + " ".join(f"{q(f(dy)):6.1f}/{q(f(dm)):6.1f}" for f in (lambda v: v, lambda v: sosfiltfilt(bA, v), lambda v: sosfiltfilt(bB, v))) + " |"
        print(line)
for e in (10.0, 3.0):
    rep((0.0, e, 0.0), [(s, e, 0.0) for s in (0.1, 0.2, 0.3, 0.4, 0.5)], f"EARLY {e}, LATE 0: SPAN s vs SPAN 0")
rep((0.3, 10.0, 0.0), [(0.3, 10.0, l) for l in (0.2, 0.4, 0.8, 1.2, 2.0)], "SPAN 0.3 EARLY 10: LATE vs LATE 0")
