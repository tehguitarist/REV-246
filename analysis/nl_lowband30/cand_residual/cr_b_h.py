"""(b) The LTI part of the error as an impulse response: h = IFFT(T - 1) from the fine-bin (4.8 Hz) cross-spectral T of the captures of one (EARLY, session) cell,
band-limited 100-1000 Hz (and 20-100 Hz separately); peaks of |h| in lags +-3000 samples."""
from cr_common import *
B = bench_all()
EAR = (1.0, 1.5, 2.0, 3.0, 10.0)
NF = 16384
fr = np.fft.rfftfreq(NF, 1 / FS)
PRE = butter(6, [12, 3000], btype='band', fs=FS, output="sos")
def db(x): return 20 * np.log10(x)
cells = {}
for k in sorted(B):
    if k[0] == 0.5 and k[1] in EAR and not B[k]["trimmed"]:
        cells.setdefault((k[1], round(B[k]["front_shift"] / 10)), []).append(k)
def Tcell(ks, grp=3):
    num = np.zeros(NF // 2 + 1, complex); den = np.zeros(NF // 2 + 1)
    for k in ks:
        y, m = cached(k, B[k]); a, b = window(k, B[k]); w = np.hanning(b - a)
        out = []
        for c in (0, 1):
            Y = np.fft.rfft(sosfiltfilt(PRE, y[:, c])[a:b] * w, NF); M = np.fft.rfft(sosfiltfilt(PRE, m[:, c])[a:b] * w, NF)
            sc = 1 / np.sum(np.abs(Y) ** 2)
            num += sc * Y * np.conj(M); den += sc * np.abs(M) ** 2
    # smooth over grp bins
    ker = np.ones(grp)
    n_ = np.convolve(num, ker, "same"); d_ = np.convolve(den, ker, "same")
    return n_ / np.maximum(d_, 1e-30)
def hresp(T, lo, hi):
    msk = ((fr >= lo) & (fr < hi)).astype(float)
    # raised-cosine edges 10 % to avoid ringing from the mask
    from scipy.ndimage import uniform_filter1d
    msk = uniform_filter1d(msk, max(3, int(0.1 * (hi - lo) / (fr[1]))))
    h = np.fft.irfft((T - 1) * msk, NF)
    return np.fft.fftshift(h), np.arange(-NF // 2, NF // 2)
res = {}
for cell, ks in sorted(cells.items()):
    T = Tcell(ks); res[cell] = T
    print(f"EARLY {cell[0]} session {cell[1]} (n={len(ks)}):")
    for lo, hi in ((20, 100), (100, 1000)):
        h, lag = hresp(T, lo, hi)
        sel = np.abs(lag) <= 4000
        hh = h[sel]; ll = lag[sel]
        tot = np.sum(h ** 2)
        order = np.argsort(-np.abs(hh))
        picked = []
        for i in order:
            if all(abs(ll[i] - p) > 12 for p, _ in picked): picked.append((ll[i], hh[i]))
            if len(picked) == 6: break
        print(f"   {lo}-{hi} Hz: |T-1| rms {db(np.sqrt(np.mean(np.abs((T-1)[(fr>=lo)&(fr<hi)])**2))):6.1f} dB; peaks (lag samples, dB re 1): " + "  ".join(f"{l:+d}:{db(abs(v)):.1f}" for l, v in picked) + f"   peak share of energy {np.sum([v*v for _, v in picked])/tot:.2f}")
np.savez(os.path.join(HERE, "b_h_T.npz"), f=fr, **{f"T_{c[0]}_{c[1]}": v for c, v in res.items()})
