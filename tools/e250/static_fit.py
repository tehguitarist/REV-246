"""Static fit of the E program (Altiverb REV-250 presets): decay law and early taps, scored held-out (docs/UAD250_test_plan.md §2.9).

    .venv/bin/python tools/e250/static_fit.py

1. Decay law. Per octave band b the decay rate r_b = 60 / T60_b (dB/s) is modelled as r_b = p_b + q_b r_1k, r_1k = 60 / RT_label, which is what a fixed
   per-pass loop filter (RT-independent loss) plus a gain set by RT gives. Fitted on 7 of the 9 RT steps, predicted on the 2 left out (all 36
   pairs); 2 parameters per band (7 bands, 14 in all). Reported: T60 error in % per band, held out, against the in-sample error.
2. Early taps. The cluster is modelled as three copies of one kernel (the first tap's response, 1.25 ms) at fitted fractional delays and gains,
   per channel; reported: tap times (ms and device samples at 26,373.6 Hz), gains, and the residual over each tap's window. Parameters: 2 delays + 2 gains
   per channel (the kernel is the data).
"""
import itertools, os, sys
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import resample_poly
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "tools", "harness"))
from harness.core import band_rt  # noqa: E402
r = 48000; FDEV = 19.2e6 / 728; OCT = [125, 250, 500, 1000, 2000, 4000, 8000]
RTS = [0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5]
def first(x, th=0.02): return int(np.argmax(np.abs(x) > th * np.abs(x).max()))
def load(RT): return np.load(os.path.join(ROOT, "analysis", "altiverb", "cache", f"REV-250 {RT}s.npy"))

H = {RT: load(RT) for RT in RTS}
lead = {RT: 0 if first(H[RT][0]) < first(H[RT][1]) else 1 for RT in RTS}   # channel index of R (leads)

def t60_table(pick):
    return np.array([[band_rt(pick(RT)[:int(min(3 * RT + 2, 12) * r)], r, f / 1.414, f * 1.414)["T"] for f in OCT] for RT in RTS])

MODELS = {  # name: (features of RT, target transform, inverse)
    "rate-linear   r = p + q r1k":      (lambda RT: np.c_[np.ones_like(RT), 60 / RT], lambda T: 60 / T, lambda y: 60 / y),
    "log-power     lnT = a + c lnRT":   (lambda RT: np.c_[np.ones_like(RT), np.log(RT)], np.log, np.exp),
    "log-quad      lnT = a + c l + d l^2": (lambda RT: np.c_[np.ones_like(RT), np.log(RT), np.log(RT) ** 2], np.log, np.exp),
    "ratio-linear  T/RT = a + c RT":    (lambda RT: np.c_[np.ones_like(RT), RT], lambda T: T, None),
}

def score(T, feat, tf, inv):
    RT = np.array(RTS); F = feat(RT); out_in, out_held = [], []
    for b in range(T.shape[1]):
        y = tf(T[:, b] / (RT if inv is None else 1)); c = np.linalg.lstsq(F, y, rcond=None)[0]
        pred = (F @ c) * (RT if inv is None else 1) if inv is None else inv(F @ c)
        out_in.append(np.sqrt(np.mean((pred / T[:, b] - 1) ** 2)) * 100); errs = []
        for a, bb in itertools.combinations(range(len(RTS)), 2):
            tr = [i for i in range(len(RTS)) if i not in (a, bb)]; c2 = np.linalg.lstsq(F[tr], y[tr], rcond=None)[0]
            for k in (a, bb):
                pk = (F[k] @ c2) * (RT[k] if inv is None else 1) if inv is None else inv(F[k] @ c2); errs.append(pk / T[k, b] - 1)
        out_held.append(np.sqrt(np.mean(np.square(errs))) * 100)
    return out_in, out_held

print("1. DECAY LAW: octave T60 versus the RT label, nine presets; in-sample and held-out (leave two RT steps out, all 36 pairs) T60 error %")
print("   band                                    " + " ".join(f"{f:6d}" for f in OCT))
TT = {}
for cn, pick in (("R", lambda RT: H[RT][lead[RT]]), ("L", lambda RT: H[RT][1 - lead[RT]])):
    T = t60_table(pick); TT[cn] = T
    print(f"  {cn}: T60/label   RT 0.4: " + " ".join(f"{T[0, b] / 0.4:6.2f}" for b in range(7)) + "\n               RT 1.0: " + " ".join(f"{T[1, b] / 1.0:6.2f}" for b in range(7)) + "\n               RT 4.5: " + " ".join(f"{T[8, b] / 4.5:6.2f}" for b in range(7)))
    for nm, (feat, tf, inv) in MODELS.items():
        i, h = score(T, feat, tf, inv); npar = feat(np.array(RTS)).shape[1] * 7
        print(f"   {nm:36s} held-out " + " ".join(f"{v:6.2f}" for v in h) + f"   ({npar} params; in-sample mean {np.mean(i):.2f} %)")

print("\n2. EARLY TAPS: three copies of the first tap's kernel (1.25 ms), fractional delays and gains fitted per channel and preset")
UP = 8
def taps(h):
    x = resample_poly(h, UP, 1); ia = first(x); n = 63 * UP; k0 = x[ia - 3 * UP: ia - 3 * UP + n]; res = []
    for guess_ms in (16.4, 27.3):
        o = int(guess_ms * 1e-3 * r * UP); a0 = ia + o - 3 * UP; W = 10 * UP; w = x[a0 - W: a0 + n + W]
        def f(p): return p[0] * np.interp(np.arange(len(w)) - W - p[1] * UP, np.arange(n), k0, left=0, right=0) - w
        best = min((least_squares(f, [0.7, d0]) for d0 in np.arange(-8, 8.1, 0.5)), key=lambda s: s.cost)
        g, d = best.x; res.append(((o + d * UP) / UP / r * 1e3, g, 20 * np.log10(np.sqrt(2 * best.cost / (w ** 2).sum()) + 1e-12)))
    return ia / UP / r * 1e3, res
print("   RT   ch  first tap ms |  tap 2: ms (device smp) gain resid dB |  tap 3: ms (device smp) gain resid dB")
tapdat = []
for RT in RTS:
    for cn, ch in (("R", lead[RT]), ("L", 1 - lead[RT])):
        t1, res = taps(H[RT][ch]); tapdat.append((RT, cn, res))
        print(f"   {RT:3.1f}  {cn}   {t1:7.3f}     | " + " | ".join(f"{t:6.3f} ({t * 1e-3 * FDEV:6.1f}) {g:+.3f} {e:6.1f}" for t, g, e in res))
g2 = np.array([d[2][0][1] for d in tapdat]); g3 = np.array([d[2][1][1] for d in tapdat]); t2 = np.array([d[2][0][0] for d in tapdat]); t3 = np.array([d[2][1][0] for d in tapdat])
print(f"\n   over 9 presets x 2 channels: tap 2 gain {g2.mean():+.3f} +- {g2.std():.3f}, at {t2.mean():.3f} +- {t2.std():.3f} ms ({t2.mean() * 1e-3 * FDEV:.1f} device samples); "
      f"tap 3 gain {g3.mean():+.3f} +- {g3.std():.3f}, at {t3.mean():.3f} +- {t3.std():.3f} ms ({t3.mean() * 1e-3 * FDEV:.1f})")

print("\n3. DECAY LAW, constant per-band multiplier T60_b = RT * m_b (the card's RT x LF / RT x HF form), RT >= 1.0 (8 presets), held-out leave-two-out")
print("   band                       " + " ".join(f"{f:6d}" for f in OCT))
RT_ = np.array(RTS); sel = RT_ >= 1.0
for cn in "RL":
    T = TT[cn]; ratio = T[sel] / RT_[sel, None]; idx = np.where(sel)[0]; errs = []
    for a, bb in itertools.combinations(range(len(idx)), 2):
        tr = [i for i in range(len(idx)) if i not in (a, bb)]; m = ratio[tr].mean(0)
        for k in (a, bb): errs.append(m / ratio[k] - 1)
    m = ratio.mean(0)
    print(f"  {cn}: m_b (all 8)             " + " ".join(f"{v:6.3f}" for v in m))
    print("     held-out T60 error %      " + " ".join(f"{v:6.2f}" for v in np.sqrt(np.mean(np.square(errs), 0)) * 100) + "   (7 parameters)")
    print("     RT 0.4 measured / (RT m_b)" + " ".join(f"{v:6.2f}" for v in T[0] / (0.4 * m)) + "   (the short preset departs from the law, mostly in the top bands)")
