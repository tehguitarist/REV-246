import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
ft, B, lag = bspline_basis_ft(0, 100, 6.0); phi = design(ft)
for nm, bld in (('compact 38', lambda p: T_from_a(fit_compact(p, phi, lo=18.0, ridge=1e-4), phi, cut=(100, 130))), ('table', lambda p: taper(T_table(p)))):
    cnt = {b: [0, 0] for b in ((20, 40), (40, 80), (80, 160))}
    for prog in ("N", "S"):
        for tag, tests, parts in folds(prog):
            T = bld(parts)
            for d in tests:
                for c in (0, 1):
                    for b in cnt:
                        m = (fc >= b[0]) & (fc < b[1]); cnt[b][1] += 1; cnt[b][0] += (10 * np.log10(err(d, T, m, c) / err(d, one, m, c)) > 0.5)
    print(nm, {b: f'{v[0]}/{v[1]}' for b, v in cnt.items()})
