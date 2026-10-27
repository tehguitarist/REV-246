"""N's own 130-1000 Hz structure (fineT_all_sm6): is it the same kernel-tail correction? (a) extrapolate the 18-130 Hz tail fit
upward; (b) refit a lag-domain tail directly on 130-1000 Hz and compare lag profiles. Unweighted per 0.4 Hz bin (N only)."""
import sys; sys.path.insert(0, '/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_kernel')
from ck_lib import *
from scipy.interpolate import BSpline
z = np.load(ROOT + '/analysis/nl/lowband_N_fine.npz'); fr = z['fr']; TF = [z['T0'], z['T1']]
dec = 4
def basis(f, t0, t1, d, nmax=200.0):
    n = np.arange(0, int(nmax * 96), dec); lag = n / 96.0
    k0 = np.arange(t0 - 3 * d, t1 + 3 * d + 1e-9, d)
    cols = []
    for i in range(len(k0) - 4):
        b = np.nan_to_num(BSpline(k0, np.eye(len(k0) - 4)[i], 3, extrapolate=False)(lag))
        if np.abs(b).max() > 0: cols.append(b)
    B = np.stack(cols, 1); E = np.exp(-2j * np.pi * np.outer(f, n) / FS)
    return (E @ B) * dec * 1e-5, B, lag
def solve(P, t, ridge):
    A = np.vstack([P.real, P.imag]); y = np.r_[t.real, t.imag]; K = A.shape[1]
    return np.linalg.solve(A.T @ A + ridge * np.trace(A.T @ A) / K * np.eye(K), A.T @ y)
mlo = (fr >= 18) & (fr <= 130); mhi = (fr > 130) & (fr <= 1000)
flo, fhi = fr[mlo], fr[mhi]
bands = ((130, 200), (200, 300), (300, 500), (500, 1000))
def frac(t, tp, f):
    out = []
    for lo, hi in bands:
        m = (f > lo) & (f <= hi); out.append(1 - np.sum(abs(t[m] - tp[m]) ** 2) / np.sum(abs(t[m]) ** 2))
    return out
# sanity: fineT vs joint-N column at 18-130
j = J(); print('rms |T_fine - T_jointN| / rms|T_jointN - 1| on 18-130 Hz:', [float(np.sqrt(np.mean(abs(np.interp(j['f'][(j['f']>=18)&(j['f']<=130)], fr, TF[c].real) + 1j*np.interp(j['f'][(j['f']>=18)&(j['f']<=130)], fr, TF[c].imag) - j[('N', c)][(j['f']>=18)&(j['f']<=130)]) ** 2)) / np.sqrt(np.mean(abs(j[('N', c)][(j['f']>=18)&(j['f']<=130)] - 1) ** 2))) for c in (0, 1)])
print('rms |T-1| by band (fine, ch0/ch1):', [(b, round(float(np.sqrt(np.mean(abs(TF[0][(fr>b[0])&(fr<=b[1])]-1)**2))),5), round(float(np.sqrt(np.mean(abs(TF[1][(fr>b[0])&(fr<=b[1])]-1)**2))),5)) for b in ((18,40),(40,80),(80,130))+bands])
for (t0, t1, d) in ((0, 100, 6.0), (0, 100, 3.0)):
    flo_ft, B, lag = basis(flo, t0, t1, d)
    fhi_ft, _, _ = basis(fhi, t0, t1, d)
    print(f'\n--- spline 0-100 ms d{d} (K={B.shape[1]}): fit on N fineT 18-130 Hz, extrapolate to 130-1000 Hz')
    for c in (0, 1):
        Plo = flo_ft * (hfac(flo) / G_of(flo, c))[:, None]; Phi = fhi_ft * (hfac(fhi) / G_of(fhi, c))[:, None]
        a = solve(Plo, TF[c][mlo] - 1, 1e-6)
        fin = 1 - np.sum(abs(TF[c][mlo] - 1 - Plo @ a) ** 2) / np.sum(abs(TF[c][mlo] - 1) ** 2)
        ex = frac(TF[c][mhi] - 1, Phi @ a, fhi)
        print(f' ch{c}: in-band 18-130 explained {fin:.3f}; extrapolated explained fraction (130-200,200-300,300-500,500-1000 Hz): {np.round(ex, 2)}; predicted rms/actual rms 130-200: {np.sqrt(np.mean(abs(Phi@a)[(fhi<=200)]**2))/np.sqrt(np.mean(abs(TF[c][mhi]-1)[(fhi<=200)]**2)):.2f}')
        np.save(f'ck_b_fineN_ch{c}_d{d:g}.npy', a)
# refit on 130-1000 directly (finer knots), compare lag profiles with the 18-130 Hz fit
print('\n--- refit directly on 130-1000 Hz (0-100 ms, d1.5, K=%d); interleaved 10 Hz blocks held out' % basis(fhi[:3], 0, 100, 1.5)[0].shape[1])
fhi_ft, B15, lag15 = basis(fhi, 0, 100, 1.5)
blk = (np.floor(fhi / 10).astype(int) % 2 == 0)
for c in (0, 1):
    Phi = fhi_ft * (hfac(fhi) / G_of(fhi, c))[:, None]; t = TF[c][mhi] - 1
    for ridge in (1e-3, 1e-1):
        a = solve(Phi[blk], t[blk], ridge)
        tr = 1 - np.sum(abs(t[blk] - Phi[blk] @ a) ** 2) / np.sum(abs(t[blk]) ** 2); te = 1 - np.sum(abs(t[~blk] - Phi[~blk] @ a) ** 2) / np.sum(abs(t[~blk]) ** 2)
        print(f' ch{c} ridge {ridge:g}: train explained {tr:.2f}, held-out blocks {te:.2f}')
    a = solve(Phi, t, 1e-1); dgh = B15 @ a * 1e-5
    # lag profile (5 ms rms) and correlation with the low-band-fit profile
    flo_ft6, B6, lag6 = basis(flo, 0, 100, 6.0); Plo = flo_ft6 * (hfac(flo) / G_of(flo, c))[:, None]; a6 = solve(Plo, TF[c][mlo] - 1, 1e-6); dgl = B6 @ a6 * 1e-5
    lg = lag6; 
    prof = lambda dg, lagx: [float(np.sqrt(np.mean(dg[(lagx >= l) & (lagx < l + 10)] ** 2))) * 1e6 for l in range(0, 100, 10)]
    print(f' ch{c} rms dg per 10 ms lag bin (e-6): low-band fit {np.round(prof(dgl, lag6), 1)}; high-band fit {np.round(prof(dgh, lag15), 1)}')
