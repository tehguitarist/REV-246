"""Leave-one-SPAN-out evaluation of the empirical LF correction T(f) (per channel, log grid) applied to the model's complex spectrum.
Evaluation captures: untrimmed, SPAN >= 0.2, LATE 0 or 0.2 (late small in the window), EARLY >= 0.3; window front-zero + 0.05 s ... min(z+0.6 s, window before the late
feed - 10 ms, usable end); null in bands from Y,M spectra (Hann) and, in time domain, from butter(4) band-passes of y and of T-corrected m.
env: NL_BK=1 NL_BK_ALPHA=0.65 LB_BENCH=bench_bk_a0.65"""
import pickle, json, collections
from lb_common import *
D = pickle.load(open("tf_all2.pkl", "rb")); Fg = D["F"]; rows = [r for r in D["rows"] if abs(r["fs"] - 26.2) > 1]
B = bench_all()
import sys
LO = float(sys.argv[1]) if len(sys.argv) > 1 else 17.0; HI = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
def table(excl_span, lo=17.0, hi=150.0):
    rs = [r for r in rows if r["key"][0] != excl_span]
    T = {}
    for c in (0, 1):
        t = sum(r["num"][:, c] for r in rs) / sum(r["den"][:, c] for r in rs)
        T[c] = t
    return T
def apply(T, fq, c, lo=None, hi=None):
    lo = LO if lo is None else lo; hi = HI if hi is None else hi
    """complex T on frequency array fq (Hz): log-interp on [lo,hi] of the grid, blended to 1 at lo/1.4 and hi*1.4"""
    m = (Fg >= lo) & (Fg <= hi); f = np.r_[lo / 1.4, Fg[m], hi * 1.4]
    lm = np.r_[0, np.log(abs(T[c][m])), 0]; ph = np.r_[0, np.unwrap(np.angle(T[c][m])), 0]
    x = np.log(np.clip(fq, 1e-3, None))
    return np.exp(np.interp(x, np.log(f), lm, left=0, right=0) + 1j * np.interp(x, np.log(f), ph, left=0, right=0))
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
BD = [(20, 32), (32, 50), (50, 100), (100, 300), (300, 1000)]
sos20 = bands(20, 100); sos100 = bands(100, 1000)
res = collections.defaultdict(list)
keys = [k for k in sorted(B) if k[0] >= 0.2 and k[2] in (0.0, 0.2) and k[1] >= 0.3 and not B[k]["trimmed"] and abs(B[k]["front_shift"] - 26.2) > 1]
for k in keys:
    r = B[k]; z = K0 + r["front_shift"]
    feed = {0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}[k[0]] if k[2] > 0 else 10 ** 9
    raw = nl.ir(nl.table()[k]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h)
    b = int(min(z + 0.6 * FS, end, K0 + r["front_shift"] + feed - 0.01 * FS + (r["o"] if abs(r["o"]) < 5 else 0))); a = int(z + 0.05 * FS)
    if b - a < int(0.12 * FS): continue
    y, m, _ = load(k, r, n=b + 64, band=False)
    T = table(k[0])
    w = np.hanning(b - a)[:, None]
    Y = np.fft.rfft(y[a:b] * w, NF, axis=0); Mo = np.fft.rfft(m[a:b] * w, NF, axis=0)
    # time domain: correct the whole model record
    N = 1 << int(np.ceil(np.log2(len(m)))); ff = np.fft.rfftfreq(N, 1 / FS)
    mc = np.stack([np.fft.irfft(np.fft.rfft(m[:, c], N) * apply(T, ff, c), N)[:len(m)] for c in (0, 1)], 1)
    for c in (0, 1):
        Tc = apply(T, fr, c)
        for lo, hi in BD:
            mk = (fr >= lo) & (fr < hi)
            den = np.sum(abs(Y[mk, c]) ** 2)
            res[("spec", lo, hi, "raw")].append(10 * np.log10(np.sum(abs(Y[mk, c] - Mo[mk, c]) ** 2) / den))
            res[("spec", lo, hi, "cor")].append(10 * np.log10(np.sum(abs(Y[mk, c] - Tc[mk] * Mo[mk, c]) ** 2) / den))
        for nm, sos in (("20-100", sos20), ("100-1k", sos100)):
            yb = sosfiltfilt(sos, y[:, c])[a:b]
            for tag, mm in (("raw", m), ("cor", mc)):
                mb = sosfiltfilt(sos, mm[:, c])[a:b]
                res[("time", nm, tag)].append(10 * np.log10(np.sum((yb - mb) ** 2) / np.sum(yb ** 2)))
            cst = (y[8:40, c] - mc[8:40, c]).mean(); ys = y[:, c] - cst * (0.5 * (1 - np.tanh((np.arange(len(y)) - 45.5) / 1.5)))
            ybs = sosfiltfilt(sos, ys)[a:b]; mb = sosfiltfilt(sos, mc[:, c])[a:b]
            res[("time", nm, "cor+step")].append(10 * np.log10(np.sum((ybs - mb) ** 2) / np.sum(ybs ** 2)))
    res[("span",)].append(k[0])
spans = np.array(res[("span",)]); n = len(spans)
print(n, "captures; leave-one-SPAN-out; median over captures and channels (raw -> corrected)")
def med(key, sel=None):
    v = np.array(res[key]).reshape(n, 2)
    if sel is not None: v = v[sel]
    return np.median(v)
for lo, hi in BD:
    print(f"  spectral {lo:4d}-{hi:<4d}: " + "  ".join(f"SPAN {s}: {med(('spec',lo,hi,'raw'), spans==s):6.1f} -> {med(('spec',lo,hi,'cor'), spans==s):6.1f}" for s in (0.2, 0.3, 0.4, 0.5)) + f" | all {med(('spec',lo,hi,'raw')):6.1f} -> {med(('spec',lo,hi,'cor')):6.1f}")
for nm in ("20-100", "100-1k"):
    print(f"  cor+step {nm}: " + "  ".join(f"SPAN {s}: {med(('time',nm,'cor+step'), spans==s):6.1f}" for s in (0.2, 0.3, 0.4, 0.5)) + f" | all {med(('time',nm,'cor+step')):6.1f}")
    print(f"  time butter {nm:7s}: " + "  ".join(f"SPAN {s}: {med(('time',nm,'raw'), spans==s):6.1f} -> {med(('time',nm,'cor'), spans==s):6.1f}" for s in (0.2, 0.3, 0.4, 0.5)) + f" | all {med(('time',nm,'raw')):6.1f} -> {med(('time',nm,'cor')):6.1f}")
v = np.array(res[("time", "20-100", "cor")]) - np.array(res[("time", "20-100", "raw")])
print("  20-100 time-domain: worse by >0.5 dB in %d of %d channel-captures, better by >3 dB in %d" % ((v > 0.5).sum(), len(v), (v < -3).sum()))
