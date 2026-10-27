"""Fitted time-domain error kernel dg(t) (all data, in-sample, for shape only) against the chain3 tail g(t)."""
import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
ft, B, lag = bspline_basis_ft(0, 100, 5.0); phi = design(ft)
parts = [(*tot([S[q] for q in S_keys()]), 1.0), (*tot(list(N.values())), 1.0)]
a = fit_compact(parts, phi, lo=18.0, ridge=1e-6)
np.save('ck_a_0_100_d5_lo18.npy', a)
dg = np.stack([B @ a[c] * 1e-5 for c in (0, 1)], 1)       # at dec 4 -> per-tap amplitude
print('lag(ms)   dg_L(e-6)  g_L(e-6; 5ms mean)  | dg_R  g_R   [per-tap amplitude]')
n = np.arange(len(lag)) * 4
for l in range(0, 100, 5):
    sel = (lag >= l) & (lag < l + 5)
    gsel = (np.arange(len(g96)) / 96.0 >= l) & (np.arange(len(g96)) / 96.0 < l + 5)
    row = [dg[sel, c].mean() * 1e6 for c in (0, 1)] 
    gm = [g96[gsel, c].mean() * 1e6 if gsel.any() else float('nan') for c in (0, 1)]
    print('%4d-%-4d  %8.2f  %8.2f   | %8.2f %8.2f' % (l, l + 5, row[0], gm[0], row[1], gm[1]))
# integrated dg (area) vs the DC sum of g
print('area of dg (sum over taps, unit: tap amplitude):', [float((dg[:, c].sum() * 4)) for c in (0, 1)], ' FIR sum:', g96.sum(0))
print('rms dg (0-100ms) / rms g(10-38.9ms):', [float(np.sqrt(np.mean(dg[(lag>=10)&(lag<38.9), c]**2)) / np.sqrt(np.mean(g96[960:3731, c]**2))) for c in (0, 1)])
