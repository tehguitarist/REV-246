"""Pooled fine-resolution T(f) (per FFT bin, then smoothed over +-w bins) from all SPAN >= 0.2 captures with LATE <= 0.2 and EARLY >= 0.3 (window front-zero + 0.05 s ...
min(0.45 s, late feed - 10 ms, usable end)); leave-one-SPAN-out evaluation in the time domain (4th-order band nulls) and a final pooled table saved to fineT_all.npz.
usage: lb_fineT3.py [smooth_bins=3]"""
import sys, pickle, collections
from lb_common import *
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
SM = int(sys.argv[1]) if len(sys.argv) > 1 else 3
FEED = {0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}
keys = []
for k in sorted(B):
    r = B[k]
    if k[0] < 0.2 or k[1] < 0.3 or k[2] > 0.2 or r["trimmed"] or abs(r["front_shift"] - 26.2) < 1 or abs(r["front_shift"] - 29.5) < 1: continue
    raw = nl.ir(nl.table()[k]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h)
    z = K0 + r["front_shift"]
    b = min(z + 0.45 * FS, end)
    if k[2] > 0: b = min(b, z + FEED[k[0]] + r["o"] - 0.01 * FS)
    a = z + 0.05 * FS
    if b - a >= 0.1 * FS: keys.append((k, int(a), int(b)))
print(len(keys), "captures:", collections.Counter(k[0][0] for k in keys))
S = {}; REC = {}
for k, a, b in keys:
    r = B[k]; y, m, _ = load(k, r, n=b + 64, band=False); w = np.hanning(b - a)[:, None]
    S[k] = (np.fft.rfft(y[a:b] * w, NF, axis=0), np.fft.rfft(m[a:b] * w, NF, axis=0)); REC[k] = (y, m, a, b)
def estimate(excl_span):
    ks = [k for k, _, _ in keys if k[0] != excl_span]
    T = []
    for c in (0, 1):
        num = sum(S[k][0][:, c] * S[k][1][:, c].conj() for k in ks); den = sum(abs(S[k][1][:, c]) ** 2 for k in ks)
        if SM > 0:
            ker = np.hanning(2 * SM + 3)[1:-1]; ker /= ker.sum()
            num = np.convolve(num, ker, "same"); den = np.convolve(den, ker, "same")
        t = num / (den + 1e-30)
        wgt = np.clip(np.minimum((fr - 14) / 4, (3500 - fr) / 700), 0, 1)
        T.append(1 + (t - 1) * wgt)
    return T
BDs = [(20, 100), (100, 300), (300, 1000), (1000, 5000), (5000, 10000)]
SOS = {b: bands(*b) for b in BDs}
res = collections.defaultdict(lambda: [[], []])
for span in (0.2, 0.3, 0.4, 0.5):
    T = estimate(span)
    for k, a, b in keys:
        if k[0] != span: continue
        y, m, a, b = REC[k]
        mc = np.stack([np.fft.irfft(np.fft.rfft(m[:, c], NF) * T[c], NF)[:len(m)] for c in (0, 1)], 1)
        for bd in BDs:
            for c in (0, 1):
                yb = sosfiltfilt(SOS[bd], y[:, c])[a:b]
                for j, mm in enumerate((m, mc)):
                    mb = sosfiltfilt(SOS[bd], mm[:, c])[a:b]
                    res[(span, bd)][j].append(10 * np.log10(np.sum((yb - mb) ** 2) / np.sum(yb ** 2)))
print(f"leave-one-SPAN-out, smoothing +-{SM} bins ({SM*FS/NF:.1f} Hz): median null raw -> corrected")
for span in (0.2, 0.3, 0.4, 0.5):
    print(f"  SPAN {span} n{len(res[(span, BDs[0])][0])//2:3d} | " + " | ".join(f"{b[0]}-{b[1]}: {np.median(res[(span,b)][0]):6.1f}->{np.median(res[(span,b)][1]):6.1f}" for b in BDs))
allr = {b: [sum((res[(s, b)][j] for s in (0.2, 0.3, 0.4, 0.5)), []) for j in (0, 1)] for b in BDs}
print("  all           | " + " | ".join(f"{b[0]}-{b[1]}: {np.median(allr[b][0]):6.1f}->{np.median(allr[b][1]):6.1f}" for b in BDs))
v = np.array(allr[(100, 300)][1]) - np.array(allr[(100, 300)][0]); print("  100-300: worse by >0.5 dB in %d of %d, better by >3 dB in %d" % ((v > 0.5).sum(), len(v), (v < -3).sum()))
v = np.array(allr[(20, 100)][1]) - np.array(allr[(20, 100)][0]); print("  20-100: worse by >0.5 dB in %d of %d, better by >3 dB in %d" % ((v > 0.5).sum(), len(v), (v < -3).sum()))
Tall = estimate(None)
np.savez(f"fineT_all_sm{SM}.npz", fr=fr, T0=Tall[0], T1=Tall[1], n=len(keys))
