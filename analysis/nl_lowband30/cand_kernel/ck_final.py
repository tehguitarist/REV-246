"""Final head-to-head (held out) of the best compact tail (0-100 ms, d6 ms, fit 18-130 Hz, ridge 1e-4; 19 coefficients per channel)
against the shipped 2 Hz table recipe and no correction; also writes tail_candidate.json (coefficients fitted on ALL N+S sums)."""
import sys, json; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
res = {}
for name, cfg in (('0-100 d6 (38 numbers)', (0, 100, 6.0)), ('0-80 d6 (32 numbers)', (0, 80, 6.0))):
    t0, t1, d = cfg
    ft, B, lag = bspline_basis_ft(t0, t1, d, nmax_ms=max(160.0, t1 + 20)); phi = design(ft)
    bld = lambda parts, phi=phi: T_from_a(fit_compact(parts, phi, lo=18.0, ridge=1e-4), phi, cut=(100, 130))
    res[name] = evaluate(bld)
res['none (0)'] = evaluate(lambda p: one)
res['table shipped recipe (232)'] = evaluate(lambda p: taper(T_table(p)))
for k, v in res.items(): print(fmt(v, k))
print('bands: 20-40 | 40-80 | 80-160 | 160-270 Hz, L/R, median dB re |Y|^2')
# per-capture counts: compact vs none, worse by > 0.5 dB in 20-80 Hz
ft, B, lag = bspline_basis_ft(0, 100, 6.0); phi = design(ft)
worse = tot_ = 0
for prog in ("N", "S"):
    for tag, tests, parts in folds(prog):
        T = T_from_a(fit_compact(parts, phi, lo=18.0, ridge=1e-4), phi, cut=(100, 130))
        for d in tests:
            for c in (0, 1):
                for (lo, hi) in ((20, 40), (40, 80), (80, 160)):
                    m = (fc >= lo) & (fc < hi); e1 = err(d, T, m, c); e0 = err(d, one, m, c); tot_ += 1; worse += (10 * np.log10(e1 / e0) > 0.5)
print('channel-bands worse than none by >0.5 dB:', worse, 'of', tot_)
parts = [(*tot([S[q] for q in S_keys()]), 1.0), (*tot(list(N.values())), 1.0)]
a = fit_compact(parts, phi, lo=18.0, ridge=1e-4)
json.dump(dict(note='chain3 kernel-tail error dg(t) = sum_k a[c][k] 1e-5 B_k(lag_ms) per FIR tap (96 kHz, FIR origin = index 0, FIR ends 38.9 ms); '
               'cubic B-splines, uniform knots every 6 ms over lags 0-100 ms (knots from -18 to 118 ms, scipy BSpline, nonzero basis kept); '
               'T(f) = 1 + FT[dg](f) * H(f) / G(f), H = lead phase * analytic HP corners; fitted 18-130 Hz on joint N+S sums (all data, in-sample); '
               'apply with T -> 1 linearly over 100-130 Hz.', a=a.tolist()), open('tail_candidate.json', 'w'))
