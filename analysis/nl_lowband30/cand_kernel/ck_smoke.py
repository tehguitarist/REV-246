import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
print(len(S), len(N))
# smoke: baseline none and shipped-table recipe, one prog N only
print(fmt(evaluate(lambda p: one, ("N",)), 'none'))
print(fmt(evaluate(lambda p: taper(T_table(p)), ("N",)), 'table shipped recipe'))
ft, B, lag = bspline_basis_ft(5, 70, 5.0); phi = design(ft); print('K', ft.shape)
def bld(parts):
    a = fit_compact(parts, phi); return T_from_a(a, phi, cut=(100, 130))
print(fmt(evaluate(bld, ("N",)), 'bspline 5-70ms d5 (cut100-130)'))
