"""Smallest spline tail within ~1 dB of the 2 Hz table (20-40, 40-80 Hz; N, S; L, R held out)."""
import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
ref = evaluate(lambda p: taper(T_table(p)))
none = evaluate(lambda p: one)
def short(r):
    d = [(r[p][(c, bi)] - ref[p][(c, bi)]) for p in ("N", "S") for bi in (0, 1) for c in (0, 1)]
    g = [(none[p][(c, bi)] - r[p][(c, bi)]) / (none[p][(c, bi)] - ref[p][(c, bi)]) for p in ("N", "S") for bi in (0, 1) for c in (0, 1)]
    return np.mean(d), np.max(d), np.mean(g)
print('table:', fmt(ref, 'table'))
rows = []
for (t0, t1) in ((0, 60), (0, 80), (0, 100), (5, 100), (10, 80)):
    for d in (5.0, 6.0, 7.5, 9.0):
        ft, B, lag = bspline_basis_ft(t0, t1, d, nmax_ms=max(160.0, t1 + 20)); phi = design(ft); K = phi[0].shape[1]
        for ridge in (1e-6, 1e-4):
            def bld(parts):
                a = fit_compact(parts, phi, lo=18.0, ridge=ridge); return T_from_a(a, phi, cut=(100, 130))
            r = evaluate(bld); m, mx, gfrac = short(r)
            print(f'{t0:>3}-{t1:<4} d{d:<4} r{ridge:g} K{K:>2}x2={2*K:>2}  mean shortfall {m:5.2f} dB  worst {mx:5.2f} dB  gain fraction of table {gfrac:5.2f}', flush=True)
