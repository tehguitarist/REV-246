"""Task 1: inverse transform of the joint T-1 (and N's, S's) to the FIR time domain: dg(n) = IFT[(T-1) G / H], band 14-130 Hz.
Prints the analytic-signal energy envelope by lag (ms, from the FIR origin; FIR ends at 38.9 ms) and the delay slope."""
from ck_common import *
j = J(); f = j['f']
m = (f >= 14) & (f <= 130)
ff = f[m]
lags_ms = np.arange(0, 200, 2.5)
for prog in ('joint', 'N', 'S'):
    print('\n==', prog)
    for ch in (0, 1):
        G = G_of(ff, ch)
        H = np.exp(2j*np.pi*ff*LEAD/FS)
        for fc in HP: H = H*(1j*ff/fc)/(1+1j*ff/fc)
        D = (j[(prog, ch)][m]-1)*G/H          # FT of the FIR-domain error kernel
        # analytic signal envelope: sum_f D e^{+i2pi f t} * df  (one-sided), Hann taper over band for low sidelobes
        w = np.hanning(len(ff)+2)[1:-1]
        t = lags_ms*1e-3
        a = (D*w)[None, :] * np.exp(2j*np.pi*np.outer(t, ff)) 
        env = np.abs(a.sum(1))*2*2.0/FS*FS  # scale arbitrary
        env /= env.max()
        # group delay slope of D*(phase) : fit phase vs f after removing nothing
        ph = np.unwrap(np.angle(D)); sl = np.polyfit(ff, ph, 1)[0]/(2*np.pi)
        print(' ch%d  phase-slope delay %.1f ms ; envelope by lag (ms:rel) ' % (ch, -sl*1e3), ' '.join('%g:%.2f' % (l, e) for l, e in zip(lags_ms[::2], env[::2]) if e > 0.1))
        # energy-weighted centroid of |env|^2
        e2 = env**2; print('     centroid %.1f ms, 25/50/75%% of energy at %s' % ((e2*lags_ms).sum()/e2.sum(), np.interp([.25,.5,.75], np.cumsum(e2)/e2.sum(), lags_ms)))
