import pickle, numpy as np
FS = 26373.626373626373
D = pickle.load(open("tf_all2.pkl", "rb")); F = D["F"]; rows = [r for r in D["rows"] if abs(r["fs"] - 26.2) > 1]
N = 1 << 17; f = np.fft.rfftfreq(N, 1 / FS)
np.set_printoptions(linewidth=200, precision=4, suppress=True)
for c in (0, 1):
    t = sum(r["num"][:, c] for r in rows) / sum(r["den"][:, c] for r in rows)
    m = (F >= 17) & (F <= 150); ff = np.r_[17 / 1.4, F[m], 150 * 1.4]
    lm = np.r_[0, np.log(abs(t[m])), 0]; ph = np.r_[0, np.unwrap(np.angle(t[m])), 0]
    x = np.log(np.clip(f, 1e-3, None))
    T = np.exp(np.interp(x, np.log(ff), lm, left=0, right=0) + 1j * np.interp(x, np.log(ff), ph, left=0, right=0))
    h = np.fft.irfft(T - 1, N)
    tt = (np.arange(N) + N // 2) % N - N // 2
    print(f"ch{c}: impulse response of (T-1): energy fraction within |t|<=ms:", {ms: round(float(np.sum(h[abs(tt) <= ms * 1e-3 * FS] ** 2) / np.sum(h ** 2)), 3) for ms in (5, 10, 20, 40, 80, 160)})
    # print envelope in 10 ms bins, signed peak
    for a in range(-160, 161, 20):
        sel = (tt >= a * 1e-3 * FS) & (tt < (a + 20) * 1e-3 * FS)
        print(f"   t {a:5d}..{a+20:<4d} ms  rms {np.sqrt(np.mean(h[sel]**2))*1e3:7.3f}e-3  peak {h[sel][np.argmax(abs(h[sel]))]*1e3:8.3f}e-3")
