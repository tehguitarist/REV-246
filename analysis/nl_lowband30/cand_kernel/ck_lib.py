"""Held-out evaluation framework for compact T(f) candidates (mirrors analysis/nl_lowband30/lb_joint.py; read-only on its inputs).
Nothing here touches the shipped tables."""
import glob, os, pickle
import numpy as np
from ck_common import *
NL = ROOT + '/analysis/nl_lowband30'
TBIN, TNB = 2.0, 6150
fc = (np.arange(TNB) + 0.5) * TBIN
FOLDS = {"A": (0.3, 0.8, 1.4, 2.2), "B": (0.6, 1.0, 1.8, 2.5)}
BD = [(20, 40), (40, 80), (80, 160), (160, 270)]
S, N = {}, {}
for p in sorted(glob.glob(ROOT + "/analysis/fullband/s/*.npz")):
    k = os.path.basename(p)[:-4]; lf, rt, hf = (float(x) for x in k.split("_")); zz = np.load(p)
    S[k] = dict(lf=lf, rt=rt, sym=zz["sym"].astype(complex), smm=zz["smm"].astype(float), syy=zz["syy"].astype(float))
N = {k: dict(span=v[0], sym=v[1], smm=v[2], syy=v[3]) for k, v in pickle.load(open(NL + "/joint_N_sums.pkl", "rb")).items()}
band = (fc >= 20) & (fc < 1600)
def tot(ds): return (sum(d["sym"] for d in ds), sum(d["smm"] for d in ds))
def norm_dt(smm): return np.array([smm[c][band].sum() for c in (0, 1)])[:, None]
def pooled(parts):
    """parts: [(sym, smm, w)] -> (num, den) per channel per bin, each program normalised by its own 20-1600 Hz |M|^2."""
    num = 0; den = 0
    for sym, smm, w in parts:
        d = norm_dt(smm); num = num + w * sym / d; den = den + w * smm / d
    return num, den
def smooth(x, SM=2):
    ker = np.ones(2 * SM + 1); return np.stack([np.convolve(x[c], ker, "same") for c in (0, 1)])
def T_table(parts, SM=2):
    num, den = pooled(parts); return smooth(num, SM) / np.maximum(smooth(den, SM), 1e-30)
def taper(T):
    """the shipped recipe: 1 -> T over 12-16 Hz, T -> 1 over 100-130 Hz."""
    w = np.clip((fc - 12) / 4, 0, 1) * np.clip((130 - fc) / 30, 0, 1)
    return 1 + w[None, :] * (T - 1)
S_keys = lambda fold=None: [k for k, v in S.items() if v["lf"] >= 1.0 and (fold is None or v["rt"] in FOLDS[fold])]
def err(d, T, m, c):
    t = T[c][m]; return np.sum(d["syy"][c][m] + abs(t) ** 2 * d["smm"][c][m] - 2 * np.real(t * np.conj(d["sym"][c][m])))
one = np.ones((2, TNB), complex)
def folds(prog):
    """yield (training parts-builder tag, test dicts, training parts) for held-out evaluation, exactly as lb_joint.py 'joint'."""
    if prog == "S":
        sN = tot(list(N.values()))
        for other, mine in (("B", "A"), ("A", "B")):
            sS = tot([S[q] for q in S_keys(other)])
            yield f"S fold {mine}", [S[k] for k, v in S.items() if v["rt"] in FOLDS[mine]], [(*sS, 1.0), (*sN, 1.0)]
    else:
        Sall = tot([S[q] for q in S_keys()])
        for sp in (0.2, 0.3, 0.4, 0.5):
            sN = tot([d for d in N.values() if d["span"] != sp])
            yield f"N span {sp}", [d for d in N.values() if d["span"] == sp], [(*Sall, 1.0), (*sN, 1.0)]
def evaluate(builder, progs=("N", "S"), verbose=False):
    """builder(parts) -> T (2, TNB) complex. Returns {prog: {(c, bi): median dB}} (median over all test captures/settings)."""
    out = {}
    for prog in progs:
        r = {}
        for tag, tests, parts in folds(prog):
            T = builder(parts)
            for d in tests:
                for c in (0, 1):
                    for bi, (lo, hi) in enumerate(BD):
                        m = (fc >= lo) & (fc < hi); den = d["syy"][c][m].sum()
                        r.setdefault((c, bi), []).append(10 * np.log10(max(err(d, T, m, c), 1e-30) / den))
        out[prog] = {k: float(np.median(v)) for k, v in r.items()}
    return out
def fmt(res, name):
    s = f"{name:34s}"
    for prog in res:
        s += f" | {prog}: " + "  ".join(f"{res[prog][(0,bi)]:6.1f}/{res[prog][(1,bi)]:6.1f}" for bi in range(len(BD)))
    return s
# ------------------------------------------------------------------ compact basis
def hfac(f):
    H = np.exp(2j * np.pi * f * LEAD / FS)
    for h in HP: H = H * (1j * f / h) / (1 + 1j * f / h)
    return H
FMAX = 400.0
fb = fc[fc <= FMAX]
Gc = np.stack([G_of(fb, c) for c in (0, 1)])           # (2, nb) chain3 response on the bin centres
Hc = hfac(fb)
def bspline_basis_ft(t0_ms, t1_ms, d_ms, dec=4, nmax_ms=160.0):
    """FT (on fb) of cubic B-splines with uniform knot spacing d_ms spanning lags [t0, t1] ms (FIR origin; the FIR ends at 38.9 ms),
    in units of kernel tap amplitude (columns scaled by 1e-5 so ridge is O(1)). Returns (nb, K) complex and the knot info."""
    from scipy.interpolate import BSpline
    n = np.arange(0, int(nmax_ms * 96), dec); lag = n / 96.0
    k0 = np.arange(t0_ms - 3 * d_ms, t1_ms + 3 * d_ms + 1e-9, d_ms)
    K = len(k0) - 4
    cols = []
    for i in range(K):
        b = BSpline(k0, np.eye(len(k0) - 4)[i], 3, extrapolate=False)(lag); b = np.nan_to_num(b)
        if np.abs(b).max() > 0: cols.append(b)
    B = np.stack(cols, 1)
    E = np.exp(-2j * np.pi * np.outer(fb, n) / FS)
    return (E @ B) * dec * 1e-5, B, lag
def design(ft):
    """phi[c] (nb, K): T-1 = phi @ a  -> (FT * H / G)."""
    return [ft * (Hc / Gc[c])[:, None] for c in (0, 1)]
def fit_compact(parts, phi, lo=14.0, hi=130.0, ridge=1e-3, SMW=0):
    """weighted complex LS (real coefficients) of the pooled per-bin T-1 on bins lo..hi; returns a (2, K)."""
    num, den = pooled(parts)
    a = []
    for c in (0, 1):
        m = (fb >= lo) & (fb <= hi)
        t = (num[c][:len(fb)] / np.maximum(den[c][:len(fb)], 1e-30) - 1)[m]; w = np.sqrt(den[c][:len(fb)][m] / den[c][:len(fb)][m].sum())
        P = phi[c][m]
        A = np.vstack([(w[:, None] * P).real, (w[:, None] * P).imag]); y = np.r_[w * t.real, w * t.imag]
        K = A.shape[1]
        a.append(np.linalg.solve(A.T @ A + ridge * np.trace(A.T @ A) / K * np.eye(K), A.T @ y))
    return np.array(a)
def T_from_a(a, phi, cut=None):
    T = np.ones((2, TNB), complex)
    for c in (0, 1):
        T[c, :len(fb)] = 1 + phi[c] @ a[c]
    if cut is not None:
        w = np.clip((cut[1] - fc) / (cut[1] - cut[0]), 0, 1)[None, :]; T = 1 + w * (T - 1)
    return T
