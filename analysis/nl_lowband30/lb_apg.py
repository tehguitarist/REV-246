"""Do the 100 Hz-1 kHz (and 20-100 Hz) errors come from the all-pass coefficients?  Per capture/channel regress the windowed spectral residual on the
derivatives of the branch with respect to g(AP124) and g(AP180/188) (relative derivative of the all-pass response times the branch spectrum).
Reports in-sample null, the fitted dg per capture, and a leave-one-EARLY-out with pooled dg."""
import sys, collections
from lb_common import *
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
BD = [(20, 100), (100, 1000), (1000, 5000)]
keys = sorted(k for k in B if k[0] == 0.5 and k[1] >= 0.7 and k[2] <= 1.0 and not B[k]["trimmed"] and abs(B[k]["front_shift"] - 22.94) < .3)
print(len(keys), "captures; EARLY", sorted({k[1] for k in keys}))
g0 = ne.G_AP
def rel(M_, f, g=g0):
    z = np.exp(-2j * np.pi * f / FS * M_)           # z^-M
    return (z ** 2 - 1) / ((1 - g * z) * (-g + z))
D = {}
for k in keys:
    r = B[k]; z0 = K0 + r["front_shift"]; a, b = int(z0 + 0.05 * FS), int(z0 + 0.45 * FS); n = b + 64
    y, m, _ = load(k, r, n=n, band=False)
    N = 1 << int(np.ceil(np.log2(n + 2 * ne.TAPS[-1] + 4096))); ff = np.fft.rfftfreq(N, 1 / FS)
    win = np.hanning(b - a)
    for c in (0, 1):
        cols = stage_columns(k, c, N, r, with_deriv=False)
        Br = sum(v for nme, v in cols.items() if nme != "front")
        out = {}
        for nm, M_ in (("g124", ne.AP_IN), ("gout", ne.AP_OUT[c])):
            col = np.fft.irfft(Br * rel(M_, ff), N)[a:b] * r["gain"][c]
            out[nm] = np.fft.rfft(col * win, NF)
        D[(k, c)] = (np.fft.rfft(y[a:b, c] * win, NF), np.fft.rfft(m[a:b, c] * win, NF), out)
print("built")
def solve(keys_, c, mk, names):
    A = []; bb = []
    for k in keys_:
        Y, Mo, cl = D[(k, c)]
        A.append(np.stack([cl[nm][mk] for nm in names], 1)); bb.append(Y[mk] - Mo[mk])
    A = np.vstack(A); bb = np.concatenate(bb)
    x = np.linalg.lstsq(np.vstack([A.real, A.imag]), np.r_[bb.real, bb.imag], rcond=None)[0]
    return x
def null(k, c, mk, names, x):
    Y, Mo, cl = D[(k, c)]
    corr = sum(xi * cl[nm][mk] for xi, nm in zip(x, names)) if x is not None else 0
    return 10 * np.log10(np.sum(abs(Y[mk] - Mo[mk] - corr) ** 2) / np.sum(abs(Y[mk]) ** 2))
for names in (["g124"], ["gout"], ["g124", "gout"]):
    for lo, hi in BD:
        mk = (fr >= lo) & (fr < hi)
        row = []
        for c in (0, 1):
            raw = []; ins = []; lo_ = []; xs = []
            x_all = solve(keys, c, mk, names); xs.append(x_all)
            for k in keys:
                raw.append(null(k, c, mk, names, None)); ins.append(null(k, c, mk, names, x_all))
                x_lo = solve([q for q in keys if q[1] != k[1]], c, mk, names)
                lo_.append(null(k, c, mk, names, x_lo))
            row.append(f"ch{c}: {np.median(raw):6.1f} -> in-sample {np.median(ins):6.1f}, leave-EARLY-out {np.median(lo_):6.1f}, dg = " + np.array2string(x_all, precision=5))
        print(f"{'+'.join(names):10s} {lo:5d}-{hi:<5d} " + " | ".join(row))
