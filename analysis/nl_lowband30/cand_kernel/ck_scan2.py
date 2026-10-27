import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
def hpcols(extra):
    if not extra: return [np.zeros((len(fb), 0), complex)] * 2
    x = lambda fc_: 1.0 / (1 + 1j * fb / fc_)
    cols = [x(3.35), x(0.5)]
    return [np.stack(cols, 1)] * 2
for (t0, t1) in ((5, 70), (0, 100), (0, 150)):
    d = 5.0
    ft, B, lag = bspline_basis_ft(t0, t1, d, nmax_ms=max(160.0, t1 + 20)); phi0 = design(ft)
    for extra in (False, True):
        hc = hpcols(extra)
        phi = [np.hstack([phi0[c], hc[c]]) for c in (0, 1)]; K = phi[0].shape[1]
        for lo in (14.0, 18.0):
            for ridge in (1e-3, 1e-6):
                def bld(parts):
                    a = fit_compact(parts, phi, lo=lo, ridge=ridge); return T_from_a(a, phi, cut=(100, 130))
                print(fmt(evaluate(bld), f'bs {t0}-{t1} d5 hp{int(extra)} lo{lo:.0f} r{ridge:g} K{K}x2'), flush=True)
