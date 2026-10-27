"""Is the 100-1k (and 20-100) residual an additive waveform common to captures of different EARLY (same SPAN, session, LATE), or tied to the model's content?
Correlation and LS amplitude ratio of e between EARLY values; also vs e of the same EARLY at another LATE (control)."""
import itertools
from lb_common import *
B = bench_all()
for bn, (lo, hi) in (("100-1k", (100, 1000)), ("20-100", (20, 100)), ("1-5k", (1000, 5000))):
    sos = bands(lo, hi)
    E = {}; Mdl = {}
    for k in sorted(B):
        if k[0] == 0.5 and k[2] in (0.0, 0.2, 0.4) and 0.7 <= k[1] <= 2.0 and abs(B[k]["front_shift"] - 22.94) < .3 and not B[k]["trimmed"] and abs(B[k]["o"] - 1.0) < 0.05:
            r = B[k]; z = K0 + r["front_shift"]; a, b = int(z + 0.05 * FS), int(z + 0.45 * FS)
            y, m, _ = load(k, r, n=b + 64, band=False)
            for c in (0, 1):
                E[(k, c)] = (sosfiltfilt(sos, y[:, c]) - sosfiltfilt(sos, m[:, c]))[a:b]; Mdl[(k, c)] = sosfiltfilt(sos, m[:, c])[a:b]
    print(f"--- {bn}: {len(E)//2} captures")
    cross = []; same = []
    for (k1, c1), (k2, c2) in itertools.combinations(sorted(E), 2):
        if c1 != c2: continue
        rho = np.corrcoef(E[(k1, c1)], E[(k2, c2)])[0, 1]; amp = (E[(k1, c1)] @ E[(k2, c2)]) / (E[(k2, c2)] @ E[(k2, c2)])
        mrat = np.sqrt(np.sum(Mdl[(k1, c1)] ** 2) / np.sum(Mdl[(k2, c2)] ** 2))
        (same if k1[1] == k2[1] else cross).append((rho, amp, mrat, np.corrcoef(Mdl[(k1, c1)], Mdl[(k2, c2)])[0, 1]))
    for nm, L in (("same EARLY, other LATE", same), ("different EARLY", cross)):
        L = np.array(L)
        print(f"  {nm:24s} n{len(L):3d}: corr(e1,e2) median {np.median(L[:,0]):5.2f}   LS amplitude e1/e2 median {np.median(L[:,1]):5.2f}   model rms ratio {np.median(L[:,2]):5.2f}   corr of models {np.median(L[:,3]):5.2f}")
