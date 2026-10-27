"""Does a fixed LTI correction T(f) of the model (data/model ratio pooled from OTHER captures, log-spaced +-12 % bands, smooth interpolation) remove the
20-100 Hz error?  Held out by session (fs0: EARLY 3, 10 | fs23: EARLY 0.7-2) and by EARLY (leave-one-EARLY-out).  Null = spectral, Hann window z+0.05..0.45 s."""
import pickle, collections, numpy as np
FS = 26373.626373626373
S = pickle.load(open("low_sel.pkl", "rb"))
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
def spec(d, t0=0.05, t1=0.45):
    a, b = int(d["z"] + t0 * FS), int(d["z"] + t1 * FS); w = np.hanning(b - a)[:, None]
    return np.fft.rfft(d["y"][a:b] * w, NF, axis=0), np.fft.rfft(d["m"][a:b] * w, NF, axis=0)
sp = {k: spec(d) for k, d in S.items()}
F = np.geomspace(6, 400, 40)
def est(keys, c):
    out = []
    for f in F:
        mk = (fr >= f / 1.12) & (fr < f * 1.12)
        out.append(sum(np.sum(sp[k][0][mk, c] * sp[k][1][mk, c].conj()) for k in keys) / sum(np.sum(abs(sp[k][1][mk, c]) ** 2) for k in keys))
    return np.array(out)
def interp(Tf, f):                      # interpolate log-magnitude and unwrapped phase on log f; 1 outside [F0, F1] edges (hold)
    lm = np.interp(np.log(f), np.log(F), np.log(abs(Tf))); ph = np.interp(np.log(f), np.log(F), np.unwrap(np.angle(Tf)))
    return np.exp(lm + 1j * ph)
def null(k, c, Tf, lo, hi):
    Y, Mo = sp[k][0][:, c], sp[k][1][:, c]
    mk = (fr >= lo) & (fr < hi)
    Mc = Mo * interp(Tf, fr[mk].clip(F[0], F[-1])).__class__(np.ones(1)) if False else None
    T = np.ones(len(fr), complex); T[mk] = interp(Tf, fr[mk]) if Tf is not None else 1
    return 10 * np.log10(np.sum(abs(Y[mk] - T[mk] * Mo[mk]) ** 2) / np.sum(abs(Y[mk]) ** 2))
def report(train, test, label):
    row = {}
    for lo, hi in ((20, 100), (100, 1000)):
        a = []; b = []
        for c in (0, 1):
            Tf = est(train, c)
            for k in test:
                a.append(null(k, c, None, lo, hi)); b.append(null(k, c, Tf, lo, hi))
        row[(lo, hi)] = (np.median(a), np.median(b))
    print(f"{label:34s} n_train {len(train):2d} n_test {len(test):2d} | 20-100: {row[(20,100)][0]:6.1f} -> {row[(20,100)][1]:6.1f} | 100-1k: {row[(100,1000)][0]:6.1f} -> {row[(100,1000)][1]:6.1f}")
fs0 = [k for k, d in S.items() if abs(d["fs"]) < .5]; fs23 = [k for k in S if k not in fs0]
report(fs23, fs0, "train fs23 -> test fs0")
report(fs0, fs23, "train fs0  -> test fs23")
for e in sorted({k[1] for k in fs23}):
    report([k for k in fs23 if k[1] != e], [k for k in fs23 if k[1] == e], f"fs23 leave-out EARLY {e}")
report(list(S), list(S), "in-sample (all, 40 bands x 2 numbers)")
