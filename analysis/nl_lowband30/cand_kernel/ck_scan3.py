"""Where in lag does the fixed error live? Support-restricted spline tails (lag from the FIR origin; FIR ends at 38.9 ms)."""
import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
print(fmt(evaluate(lambda p: one), 'none'))
print(fmt(evaluate(lambda p: taper(T_table(p))), 'table shipped recipe (232 nums)'))
for (t0, t1) in ((0, 10), (10, 38.9), (38.9, 80), (0, 38.9), (10, 80), (38.9, 150), (60, 150), (0, 100)):
    ft, B, lag = bspline_basis_ft(t0, t1, 5.0, nmax_ms=max(160.0, t1 + 20)); phi = design(ft); K = phi[0].shape[1]
    def bld(parts):
        a = fit_compact(parts, phi, lo=18.0, ridge=1e-6); return T_from_a(a, phi, cut=(100, 130))
    print(fmt(evaluate(bld), f'bs {t0}-{t1} d5 lo18 K{K}x2'), flush=True)
