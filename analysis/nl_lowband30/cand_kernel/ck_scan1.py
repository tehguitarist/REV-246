import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
print(fmt(evaluate(lambda p: one), 'none'))
print(fmt(evaluate(lambda p: taper(T_table(p))), 'table shipped recipe (232 nums)'))
print(fmt(evaluate(lambda p: T_table(p)), 'table raw 2Hz (lb_joint joint)'))
for (t0, t1) in ((0, 60), (5, 70), (10, 60), (10, 100), (0, 100)):
    for d in (4.0, 5.0, 7.5, 10.0):
        ft, B, lag = bspline_basis_ft(t0, t1, d); phi = design(ft); K = phi[0].shape[1]
        for cut in ((100, 130), None):
            def bld(parts):
                a = fit_compact(parts, phi); return T_from_a(a, phi, cut=cut)
            print(fmt(evaluate(bld), f'bs {t0}-{t1}ms d{d} K{K}x2 cut{cut}'), flush=True)
