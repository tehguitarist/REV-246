"""Few-parameter time-domain forms for the kernel error dg(t): per channel and per form, nonlinear shape parameters by grid
search on the TRAINING sums (variable projection; amplitude linear), held out as in lb_joint. Parameter counts are per channel x2."""
import sys, itertools; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
dec = 4; nms = 600.0
n = np.arange(0, int(nms * 96), dec); tms = n / 96.0
E = np.exp(-2j * np.pi * np.outer(fb, n) / FS)
def ft(b): return (E @ b) * dec * 1e-5
def f_gauss(c, s): return np.exp(-0.5 * ((tms - c) / s) ** 2)
def f_exp(ts, tau): return np.where(tms >= ts, np.exp(-(tms - ts) / tau), 0.0)
def f_gamma(k, tau): return (tms / tau) ** k * np.exp(-tms / tau)
def f_hp(): return None
forms = {
 'gauss(c,s)': ([(c, s) for c in (10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80) for s in (3, 5, 8, 12, 18, 25, 35)], lambda p: f_gauss(*p)),
 'exp-tail(ts,tau)': ([(ts, tau) for ts in (0, 5, 10, 20, 30, 38.9) for tau in (5, 10, 20, 30, 50, 80, 120, 200)], lambda p: f_exp(*p)),
 'gamma(k,tau)': ([(k, tau) for k in (1, 2, 3) for tau in (5, 8, 12, 18, 25, 35, 50)], lambda p: f_gamma(*p)),
}
def fit1(parts, cols_list, lo=18.0, hi=130.0):
    """cols_list[c]: (nb, k) design for channel c. Returns per-channel (coef, weighted residual)."""
    num, den = pooled(parts); out = []
    for c in (0, 1):
        m = (fb >= lo) & (fb <= hi)
        t = (num[c][:len(fb)] / np.maximum(den[c][:len(fb)], 1e-30) - 1)[m]; w = np.sqrt(den[c][:len(fb)][m] / den[c][:len(fb)][m].sum())
        P = cols_list[c][m]
        A = np.vstack([(w[:, None] * P).real, (w[:, None] * P).imag]); y = np.r_[w * t.real, w * t.imag]
        co, *_ = np.linalg.lstsq(A, y, rcond=None); out.append((co, float(np.sum((A @ co - y) ** 2))))
    return out
def make_builder(name, ncomp=1):
    grid, fn = forms[name]
    cache = {}
    def col(p, c):
        if (p, c) not in cache: cache[(p, c)] = ft(fn(p)) * Hc / Gc[c]
        return cache[(p, c)]
    def bld(parts):
        T = np.ones((2, TNB), complex)
        for c in (0, 1):
            best = None
            if ncomp == 1:
                for p in grid:
                    r = fit1(parts, [None, None] if False else [col(p, 0)[:, None], col(p, 1)[:, None]])[c]
                    if best is None or r[1] < best[1]: best = (p, r[1], r[0])
                ps = [best[0]]
            else:   # greedy second component
                for p in grid:
                    r = fit1(parts, [col(p, 0)[:, None], col(p, 1)[:, None]])[c]
                    if best is None or r[1] < best[1]: best = (p, r[1], r[0])
                p1 = best[0]; best2 = None
                for p in grid:
                    if p == p1: continue
                    M = [np.stack([col(p1, 0), col(p, 0)], 1), np.stack([col(p1, 1), col(p, 1)], 1)]
                    r = fit1(parts, M)[c]
                    if best2 is None or r[1] < best2[1]: best2 = (p, r[1], r[0])
                best = best2; ps = [p1, best2[0]]
            Pm = np.stack([col(p, c) for p in ps], 1)
            T[c, :len(fb)] = 1 + Pm @ best[2]
            bld.last = getattr(bld, 'last', {}); bld.last[c] = (ps, best[2])
        w = np.clip((130 - fc) / 30, 0, 1)[None, :]; return 1 + w * (T - 1)
    return bld
if __name__ == '__main__':
    smoke = len(sys.argv) > 1 and sys.argv[1] == 'smoke'
    for name in forms:
        for nc in (1, 2):
            if smoke and nc == 2: continue
            b = make_builder(name, nc)
            r = evaluate(b, ("N",) if smoke else ("N", "S"))
            print(fmt(r, f'{name} x{nc} ({3*nc if "exp" in name or "gauss" in name else 2*nc}+/ch)'), '  last fit:', {c: (b.last[c][0], np.round(b.last[c][1], 2).tolist()) for c in (0, 1)}, flush=True)
