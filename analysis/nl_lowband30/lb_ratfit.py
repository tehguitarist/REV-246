"""Rational s-domain fit (real coefficients, order n, T(inf)=1, Sanathanan-Koerner) of the empirical LF correction T(f) per channel, and its leave-one-SPAN-out
evaluation on the same captures / metrics as lb_tf_eval.py.  Parameters: 2n per channel.  usage: lb_ratfit.py n [fit_lo fit_hi]
env: NL_BK=1 NL_BK_ALPHA=0.65 LB_BENCH=bench_bk_a0.65"""
import sys, pickle, collections
from lb_common import *
D = pickle.load(open("tf_all2.pkl", "rb")); Fg = D["F"]; rows = [r for r in D["rows"] if abs(r["fs"] - 26.2) > 1]
W0 = 2 * np.pi * 30.0
def fit_rat(f, T, n, iters=40):
    s = 1j * 2 * np.pi * f / W0; a = np.zeros(n)
    for it in range(iters):
        Dp = 1 + sum(a[k - 1] * s ** k for k in range(1, n + 1))
        cols = [s ** k / Dp for k in range(n)] + [-T * s ** k / Dp for k in range(1, n)] + [(s ** n - T * s ** n) / Dp]
        A = np.stack(cols, 1); rhs = T / Dp
        x = np.linalg.lstsq(np.vstack([A.real, A.imag]), np.r_[rhs.real, rhs.imag], rcond=None)[0]
        a = x[n:]
        b = x[:n]
    return np.r_[b, a[-1]], np.r_[1.0, a]
def H(num, den, fq):
    s = 1j * 2 * np.pi * np.asarray(fq) / W0
    return sum(num[k] * s ** k for k in range(len(num))) / sum(den[k] * s ** k for k in range(len(den)))
def make(excl_span, n, lo, hi):
    out = {}
    rs = [r for r in rows if r["key"][0] != excl_span]
    for c in (0, 1):
        t = sum(r["num"][:, c] for r in rs) / sum(r["den"][:, c] for r in rs)
        m = (Fg >= lo) & (Fg <= hi)
        out[c] = fit_rat(Fg[m], t[m], n)
    return out
if __name__ == "__main__":
    n = int(sys.argv[1]); lo = float(sys.argv[2]) if len(sys.argv) > 2 else 14; hi = float(sys.argv[3]) if len(sys.argv) > 3 else 150
    B = bench_all()
    NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
    sos20 = bands(20, 100); sos100 = bands(100, 1000)
    res = collections.defaultdict(list); spans = []
    full = make(None, n, lo, hi)
    for c in (0, 1):
        t = sum(r["num"][:, c] for r in rows) / sum(r["den"][:, c] for r in rows); m = (Fg >= lo) & (Fg <= hi)
        e = H(*full[c], Fg[m]) - t[m]; pol = np.roots(full[c][1][::-1])
        print(f"n={n} ch{c}: fit error max {abs(e).max():.4f} ({20*np.log10(abs(e).max()+1e-12):.1f} dB), poles (s/W0 units) max real {pol.real.max():.2f}")
    keys = [k for k in sorted(B) if k[0] >= 0.2 and k[2] in (0.0, 0.2) and k[1] >= 0.3 and not B[k]["trimmed"] and abs(B[k]["front_shift"] - 26.2) > 1]
    for k in keys:
        r = B[k]; z = K0 + r["front_shift"]
        feed = {0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}[k[0]] if k[2] > 0 else 10 ** 9
        raw = nl.ir(nl.table()[k]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h)
        b = int(min(z + 0.6 * FS, end, K0 + r["front_shift"] + feed - 0.01 * FS + (r["o"] if abs(r["o"]) < 5 else 0))); a = int(z + 0.05 * FS)
        if b - a < int(0.12 * FS): continue
        y, m, _ = load(k, r, n=b + 64, band=False)
        P = make(k[0], n, lo, hi)
        N = 1 << int(np.ceil(np.log2(len(m)))); ff = np.fft.rfftfreq(N, 1 / FS)
        mc = np.stack([np.fft.irfft(np.fft.rfft(m[:, c], N) * H(*P[c], ff), N)[:len(m)] for c in (0, 1)], 1)
        for c in (0, 1):
            for nm, sos in (("20-100", sos20), ("100-1k", sos100)):
                yb = sosfiltfilt(sos, y[:, c])[a:b]
                for tag, mm in (("raw", m), ("cor", mc)):
                    mb = sosfiltfilt(sos, mm[:, c])[a:b]
                    res[(nm, tag)].append(10 * np.log10(np.sum((yb - mb) ** 2) / np.sum(yb ** 2)))
        spans.append(k[0])
    sp = np.array(spans)
    for nm in ("20-100", "100-1k"):
        v = {t: np.array(res[(nm, t)]).reshape(-1, 2) for t in ("raw", "cor")}
        print(f"  n={n} time butter {nm:7s}: " + "  ".join(f"SPAN {s}: {np.median(v['raw'][sp==s]):6.1f} -> {np.median(v['cor'][sp==s]):6.1f}" for s in (0.2, 0.3, 0.4, 0.5)) + f" | all {np.median(v['raw']):6.1f} -> {np.median(v['cor']):6.1f}")
