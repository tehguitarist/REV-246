"""Parametric LF corrections of the model (families with 1-3 numbers), fitted on one session, tested on the other + leave-EARLY-out:
A: low-edge window re-shaped: T = W(f; a, b) / W(f; 3, 12), W = sin^2(pi/2 clip((f-a)/(b-a))) ; B: first-order HP T = jf/(jf+fc) ; C: A x B."""
import pickle, numpy as np, itertools
exec(open("lb_tfcorr.py").read().split("def interp")[0])
from scipy.optimize import minimize
fs0 = [k for k, d in S.items() if abs(d["fs"]) < .5]; fs23 = [k for k in S if k not in fs0]
def W(f, a, b): return np.sin(np.pi / 2 * np.clip((f - a) / max(b - a, 1e-3), 0, 1)) ** 2
W0 = W(fr, 3, 12)
def Tf(th, fam):
    f = np.maximum(fr, 1e-3)
    if fam == "A": return (W(fr, th[0], th[1]) + 1e-9) / (W0 + 1e-9)
    if fam == "B": return 1j * f / (1j * f + th[0])
    if fam == "C": return (W(fr, th[0], th[1]) + 1e-9) / (W0 + 1e-9) * 1j * f / (1j * f + th[2])
def cost(th, fam, keys, lo):
    mk = (fr >= lo) & (fr < 200); t = Tf(th, fam)[mk]
    n = d = 0.0
    for k in keys:
        for c in (0, 1):
            n += np.sum(abs(sp[k][0][mk, c] - t * sp[k][1][mk, c]) ** 2); d += np.sum(abs(sp[k][0][mk, c]) ** 2)
    return n / d
def fit(fam, keys, lo):
    starts = {"A": [[3, 12], [2, 10], [4, 14], [1, 8], [5, 16]], "B": [[0.3], [1.0], [3.0], [0.05]], "C": [[3, 12, 0.3], [3, 12, 1.0], [2, 10, 0.1]]}[fam]
    best = None
    for s in starts:
        r = minimize(cost, s, args=(fam, keys, lo), method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-10, "maxiter": 600})
        if best is None or r.fun < best.fun: best = r
    return best.x
def nulls(keys, fam, th, lo_hi):
    out = []
    for lo, hi in lo_hi:
        mk = (fr >= lo) & (fr < hi); t = Tf(th, fam)[mk] if th is not None else 1
        out.append(np.median([10 * np.log10(np.sum(abs(sp[k][0][mk, c] - t * sp[k][1][mk, c]) ** 2) / np.sum(abs(sp[k][0][mk, c]) ** 2)) for k in keys for c in (0, 1)]))
    return out
bands_ = [(20, 100), (100, 1000)]
for lo_fit in (5, 16):
    print(f"--- cost band {lo_fit}-200 Hz")
    for fam in "ABC":
        for lab, tr, te in (("fs23->fs0", fs23, fs0), ("fs0->fs23", fs0, fs23), ("LOO E1.0", [k for k in fs23 if k[1] != 1.0], [k for k in fs23 if k[1] == 1.0]), ("LOO E2.0", [k for k in fs23 if k[1] != 2.0], [k for k in fs23 if k[1] == 2.0])):
            th = fit(fam, tr, lo_fit); b = nulls(te, fam, None, bands_); a = nulls(te, fam, th, bands_)
            print(f"  {fam} {lab:10s} th={np.array2string(th, precision=3):22s} 20-100 {b[0]:6.1f} -> {a[0]:6.1f} | 100-1k {b[1]:6.1f} -> {a[1]:6.1f}")
