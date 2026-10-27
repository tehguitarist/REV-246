"""Fine time structure of the 100-1k residual: relative error energy vs local model energy in 1-ms bins (slope of log e^2 vs log m^2), and the error energy
split into 'near a tap' (+-8 samples of any main-tap position, front-frame) vs elsewhere, SPAN 0.5, EARLY>=0.7, LATE<=0.4, window 0.05-0.45 s."""
from lb_common import *
B = bench_all(); sos = bands(100, 1000)
keys = sorted(k for k in B if k[0] == 0.5 and k[1] >= 0.7 and k[2] <= 0.4 and not B[k]["trimmed"])
nb = int(0.001 * FS)
rs = []; near_frac = []; tapfrac = []
for k in keys:
    r = B[k]; z = K0 + r["front_shift"]; a, b = int(z + 0.05 * FS), int(z + 0.45 * FS)
    y, m, _ = load(k, r, n=b + 64, band=False)
    for c in (0, 1):
        yb, mb = sosfiltfilt(sos, y[:, c]), sosfiltfilt(sos, m[:, c]); e = (yb - mb)[a:b]; mm = mb[a:b]
        E_ = e[:len(e) // nb * nb].reshape(-1, nb); M_ = mm[:len(mm) // nb * nb].reshape(-1, nb)
        ee, mmm = np.sum(E_ ** 2, 1), np.sum(M_ ** 2, 1)
        ok = mmm > 1e-12
        sl = np.polyfit(np.log(mmm[ok]), np.log(ee[ok] + 1e-30), 1)[0]
        rs.append(sl)
        # near-tap mask: main taps (device samples from the front-frame zero) + o
        taps = np.array([t + (ne.R_OFF if c else 0) for t in ne.TAPS]) + r["o"] + K0 + r["front_shift"]
        mask = np.zeros(len(e), bool)
        for t in taps:
            lo_, hi_ = int(t - 8 - a), int(t + 8 - a)
            if hi_ > 0 and lo_ < len(e): mask[max(lo_, 0):hi_] = True
        near_frac.append((np.sum(e[mask] ** 2) / np.sum(e ** 2), mask.mean(), np.sum(mm[mask] ** 2) / np.sum(mm ** 2)))
print("log-log slope of 1-ms error energy vs model energy: median %.2f (1 = error proportional to local level; 0 = uniform)" % np.median(rs))
nf = np.array(near_frac)
print("error energy within +-8 samples of tap positions: %.3f of total (those samples are %.3f of the window, and hold %.3f of the model energy)" % tuple(np.median(nf, 0)))
