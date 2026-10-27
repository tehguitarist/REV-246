"""(a) cont.: (i) session x EARLY-class table; (ii) few-parameter physical forms for T(f<18) fitted on some groups, tested on others (LS-metric energy reduction, held out)."""
from solve import *
from scipy.optimize import minimize
import collections
D = pickle.load(open(os.path.join(HERE, "g2_all.pkl"), "rb")); v = "w60"; lam = 1e-2
Ts = LJ.response(KNOTS, 2)
KI = [i for i in range(len(KNOTS)) if KNOTS[i] <= 18]; cols = np.array([2 * i + j for i in KI for j in (0, 1)]); fk = KNOTS[KI]
keys = sorted(D); ecl = lambda k: "E<=1" if k[1] <= 1.0 else ("E1.5-3" if k[1] <= 3.0 else "E10")
def est(ks, ch):
    G, b, yy, n = pool(D, ks, v, ch); Gs = G[np.ix_(cols, cols)]
    c = np.linalg.solve(Gs + lam * np.mean(np.diag(Gs)) * np.eye(len(cols)), b[cols]); return (1 + c[0::2] + 1j * c[1::2]) * Ts[KI, ch]
print("(i) T_total (dB / deg) at 10, 12, 14 Hz by session x EARLY class (L | R)")
for s in (0, 16, 20, 23, 36):
    for e in ("E<=1", "E1.5-3", "E10"):
        ks = [k for k in keys if D[k]["sess"] == s and ecl(k) == e]
        if len(ks) < 8: continue
        out = []
        for ch in (0, 1):
            T = est(ks, ch); out.append(" ".join(f"{20*np.log10(abs(T[list(fk).index(f)])):+5.1f}" for f in (6, 10, 12, 14)))
        print(f"fs{s:3d} {e:7s} n{len(ks):3d}  L(6,10,12,14Hz): {out[0]}   R: {out[1]}")
# (ii) forms
edge = lambda f, a, b: np.sin(np.pi / 2 * np.clip((f - a) / (b - a), 0, 1)) ** 2
def F_edge(p): a, b = p; return edge(fk, a, b) / np.maximum(edge(fk, 3.0, 12.0), 1e-9) + 0j
def F_corner(p): fc = p[0]; H = lambda fcc: (1j * fk / fcc) / (1 + 1j * fk / fcc); return H(fc) / H(3.35)
def F_both(p): return F_edge(p[:2]) * F_corner(p[2:])
def F_const(p): A, tau = p; return A * np.exp(-2j * np.pi * fk * tau / 26373.6 * 26373.6 * 0 - 2j * np.pi * fk * tau * 1e-3)   # tau in ms
FORMS = {"edge(a,b)": (F_edge, [3.0, 12.0], 2), "HP corner fc": (F_corner, [3.35], 1), "edge+corner": (F_both, [3.0, 12.0, 3.35], 3), "gain+delay": (F_const, [1.0, 0.0], 2)}
def red(ks, ch, T):
    G, b, yy, n = pool(D, ks, v, ch); c = np.zeros(NC); d = T - 1                  # T applied on top of the shipped base
    c[cols[0::2]] = d.real; c[cols[1::2]] = d.imag; return 2 * c @ b - c @ G @ c, yy
def nonpar(ks_tr, ks_te, ch):
    G, b, yy, n = pool(D, ks_tr, v, ch); Gs = G[np.ix_(cols, cols)]; c = np.zeros(NC); c[cols] = np.linalg.solve(Gs + lam * np.mean(np.diag(Gs)) * np.eye(len(cols)), b[cols])
    G2, b2, yy2, _ = pool(D, ks_te, v, ch); return (2 * c @ b2 - c @ G2 @ c) / yy2
tests = {"train fs 0,16,36 -> test fs 20,23": ([0, 16, 36], [20, 23]), "train fs 20,23 -> test fs 0,16,36": ([20, 23], [0, 16, 36]), "train fs 23 -> test fs 20": ([23], [20]),
         "train fs 20 -> test fs 23": ([20], [23]), "train fs 0 -> test fs 23": ([0], [23]), "train fs 16 -> test fs 36": ([16], [36])}
print("\n(ii) held-out LS-metric reduction of the whole-spectrum residual energy (percent of residual energy removed; positive = better), mean of L,R")
print(f"{'':36s}" + "".join(f"{n:>14s}" for n in list(FORMS) + ["nonpar 18/ch"]))
for lab, (tr, te) in tests.items():
    ktr = [k for k in keys if D[k]["sess"] in tr]; kte = [k for k in keys if D[k]["sess"] in te]; row = []
    for nm, (F, p0, npar) in FORMS.items():
        r = []
        for ch in (0, 1):
            obj = lambda p: -red(ktr, ch, F(p) * 1.0)[0] / pool(D, ktr, v, ch)[2]
            p = minimize(obj, p0, method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-9, maxiter=400)).x
            rr, yy = red(kte, ch, F(p)); r.append(100 * rr / yy)
        row.append(np.mean(r))
    row.append(np.mean([100 * nonpar(ktr, kte, ch) for ch in (0, 1)]))
    print(f"{lab:36s}" + "".join(f"{x:14.3f}" for x in row))
# fitted parameter values on everything (in-sample, for the record)
print("\nfitted on all captures (L, R):")
for nm, (F, p0, npar) in FORMS.items():
    out = []
    for ch in (0, 1):
        obj = lambda p: -red(keys, ch, F(p))[0] / pool(D, keys, v, ch)[2]
        p = minimize(obj, p0, method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-9, maxiter=600)).x; out.append((np.round(p, 2), round(100 * red(keys, ch, F(p))[0] / pool(D, keys, v, ch)[2], 3)))
    print(f"  {nm:14s} L {out[0]}  R {out[1]}   (nonparametric all: L {100*nonpar(keys, keys, 0):.3f}  R {100*nonpar(keys, keys, 1):.3f})")
