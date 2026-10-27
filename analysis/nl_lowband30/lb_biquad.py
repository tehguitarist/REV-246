"""Compact LF correction: model x T(f), T = (1 + c1 jw - c2 w^2)/(1 + d1 jw - d2 w^2) (4 numbers, both channels), fitted on part of the captures by the spectral
error in 8-200 Hz (Hann window z+0.05..0.45 s), tested held out by session / EARLY on the 20-100 Hz and 100 Hz-1 kHz band error (spectral, and the bench's time-domain
butter(4) bands)."""
import pickle, numpy as np, json
from scipy.optimize import minimize
from scipy.signal import butter, sosfiltfilt
exec(open("lb_tfcorr.py").read().split("def interp")[0])
fs0 = [k for k, d in S.items() if abs(d["fs"]) < .5]; fs23 = [k for k in S if k not in fs0]
w = 2 * np.pi * fr
def T(th, f=None):
    ww = 2 * np.pi * (fr if f is None else f)
    c1, c2, d1, d2 = th
    return (1 + 1j * ww * c1 - ww ** 2 * c2) / (1 + 1j * ww * d1 - ww ** 2 * d2)
band_fit = (fr >= 16) & (fr < 200)
def cost(th, keys):
    t = T(th)[band_fit]
    num = den = 0.0
    for k in keys:
        for c in (0, 1):
            Y, Mo = sp[k][0][band_fit, c], sp[k][1][band_fit, c]
            num += np.sum(abs(Y - t * Mo) ** 2); den += np.sum(abs(Y) ** 2)
    return num / den
def fit(keys):
    best = None
    for th0 in ([0.003, 1e-4, 0.002, 1e-4], [0.02, 0.0, 0.02, 0.0], [0.0, 2e-4, 0.01, 5e-4], [0.01, 1e-3, 0.004, 1e-4], [-0.01, 1e-3, 0.004, 1e-3]):
        r = minimize(cost, th0, args=(keys,), method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-12, "maxiter": 4000, "maxfev": 4000})
        if best is None or r.fun < best.fun: best = r
    return best.x, best.fun
def specnull(keys, th, lo, hi):
    mk = (fr >= lo) & (fr < hi); t = T(th)[mk] if th is not None else 1
    return np.median([10 * np.log10(np.sum(abs(sp[k][0][mk, c] - t * sp[k][1][mk, c]) ** 2) / np.sum(abs(sp[k][0][mk, c]) ** 2)) for k in keys for c in (0, 1)])
def timenull(keys, th, lo, hi, order=4):
    sos = butter(order, [lo, hi], btype="band", fs=FS, output="sos"); v = []
    for k in keys:
        d = S[k]; a, b = int(d["z"] + 0.05 * FS), int(d["z"] + 0.45 * FS)
        for c in (0, 1):
            m = d["m"][:, c]
            if th is not None:
                N = 1 << 15; Mf = np.fft.rfft(m, N); f = np.fft.rfftfreq(N, 1 / FS); m = np.fft.irfft(Mf * T(th, f), N)[:len(m)]
            yb, mb = sosfiltfilt(sos, d["y"][:, c]), sosfiltfilt(sos, m)
            v.append(10 * np.log10(np.sum((yb[a:b] - mb[a:b]) ** 2) / np.sum(yb[a:b] ** 2)))
    return np.median(v)
def rep(train, test, label):
    th, f = fit(train)
    print(f"{label:30s} th={np.array2string(th, precision=5)} train cost {10*np.log10(f):6.1f} dB | test spectral 20-100: {specnull(test,None,20,100):6.1f} -> {specnull(test,th,20,100):6.1f}  100-1k: {specnull(test,None,100,1000):6.1f} -> {specnull(test,th,100,1000):6.1f} | time-domain butter 20-100: {timenull(test,None,20,100):6.1f} -> {timenull(test,th,20,100):6.1f}  100-1k: {timenull(test,None,100,1000):6.1f} -> {timenull(test,th,100,1000):6.1f}")
    return th
th1 = rep(fs23, fs0, "train fs23 -> test fs0")
th2 = rep(fs0, fs23, "train fs0 -> test fs23")
for e in (0.7, 1.0, 1.5, 2.0):
    rep([k for k in fs23 if k[1] != e], [k for k in fs23 if k[1] == e], f"fs23 leave-out EARLY {e}")
thall = rep(list(S), list(S), "in-sample (reference)")
json.dump({"train_fs23": list(th1), "train_fs0": list(th2), "all": list(thall)}, open("biquad.json", "w"))
