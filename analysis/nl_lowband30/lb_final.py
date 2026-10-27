"""Final pooled estimate of the low-band correction on a log grid (band = +-3 % around each grid frequency), leave-one-SPAN-out evaluation of the
SAME procedure (grid, interpolation) in the time domain, and the saved table lowband_T.json.   usage: lb_final.py"""
import sys, json, collections
from lb_common import *
import lowband_fix as lf
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
FEED = {0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}
G = lf.grid(); print(len(G), "grid points", G[0], G[-1])
keys = []
for k in sorted(B):
    r = B[k]
    if k[0] < 0.2 or k[1] < 0.3 or k[2] > 0.2 or r["trimmed"] or abs(r["front_shift"] - 26.2) < 1 or abs(r["front_shift"] - 29.5) < 1: continue
    raw = nl.ir(nl.table()[k]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h); z = K0 + r["front_shift"]
    b = min(z + 0.45 * FS, end)
    if k[2] > 0: b = min(b, z + FEED[k[0]] + r["o"] - 0.01 * FS)
    a = z + 0.05 * FS
    if b - a >= 0.1 * FS: keys.append((k, int(a), int(b)))
S = {}; REC = {}
for k, a, b in keys:
    r = B[k]; y, m, _ = load(k, r, n=b + 64, band=False); w = np.hanning(b - a)[:, None]
    S[k] = (np.fft.rfft(y[a:b] * w, NF, axis=0), np.fft.rfft(m[a:b] * w, NF, axis=0)); REC[k] = (y, m, a, b)
masks = [(fr >= f / 1.03) & (fr < f * 1.03) for f in G]
def estimate(excl):
    ks = [k for k, _, _ in keys if k[0] != excl]; out = []
    for c in (0, 1):
        t = []
        for mk in masks:
            num = sum(np.sum(S[k][0][mk, c] * S[k][1][mk, c].conj()) for k in ks); den = sum(np.sum(abs(S[k][1][mk, c]) ** 2) for k in ks)
            t.append(num / den)
        out.append(np.array(t))
    return out
def to_full(T):
    lf._T = {"f_hz": list(G), "db0": list(20 * np.log10(abs(T[0]))), "deg0": list(np.degrees(np.angle(T[0]))), "db1": list(20 * np.log10(abs(T[1]))), "deg1": list(np.degrees(np.angle(T[1])))}
    return [lf.interp(fr, c) for c in (0, 1)]
BDs = [(20, 100), (100, 300), (300, 1000), (1000, 5000), (5000, 10000)]
SOS = {b: bands(*b) for b in BDs}
res = collections.defaultdict(lambda: [[], [], []])
for span in (0.2, 0.3, 0.4, 0.5):
    T = to_full(estimate(span))
    for k, _, _ in keys:
        if k[0] != span: continue
        y, m, a, b = REC[k]
        mc = np.stack([np.fft.irfft(np.fft.rfft(m[:, c], NF) * T[c], NF)[:len(m)] for c in (0, 1)], 1)
        for bd in BDs:
            for c in (0, 1):
                yb = sosfiltfilt(SOS[bd], y[:, c])[a:b]
                for j, mm in enumerate((m, mc)):
                    mb = sosfiltfilt(SOS[bd], mm[:, c])[a:b]
                    res[(span, bd)][j].append(10 * np.log10(np.sum((yb - mb) ** 2) / np.sum(yb ** 2)))
print("leave-one-SPAN-out (log grid, 2 x %d complex numbers): median null raw -> corrected" % len(G))
for span in (0.2, 0.3, 0.4, 0.5):
    print(f"  SPAN {span} n{len(res[(span, BDs[0])][0])//2:3d} | " + " | ".join(f"{b[0]}-{b[1]}: {np.median(res[(span,b)][0]):6.1f}->{np.median(res[(span,b)][1]):6.1f}" for b in BDs))
allr = {b: [sum((res[(s, b)][j] for s in (0.2, 0.3, 0.4, 0.5)), []) for j in (0, 1)] for b in BDs}
print("  all           | " + " | ".join(f"{b[0]}-{b[1]}: {np.median(allr[b][0]):6.1f}->{np.median(allr[b][1]):6.1f}" for b in BDs))
for b in BDs[:3]:
    v = np.array(allr[b][1]) - np.array(allr[b][0]); print("  %d-%d: worse by >0.5 dB in %d of %d, better by >3 dB in %d" % (b[0], b[1], (v > 0.5).sum(), len(v), (v < -3).sum()))
T = estimate(None)
json.dump({"f_hz": list(map(float, G)), "db0": list(map(float, 20 * np.log10(abs(T[0])))), "deg0": list(map(float, np.degrees(np.angle(T[0])))),
           "db1": list(map(float, 20 * np.log10(abs(T[1])))), "deg1": list(map(float, np.degrees(np.angle(T[1])))),
           "n_captures": len(keys), "note": "data/model complex transfer ratio, SPAN>=0.2 LATE<=0.2 EARLY>=0.3, window front-zero+0.05 s .. min(0.45 s, late feed); bench_bk_a0.65 alignment"}, open("lowband_T.json", "w"), indent=1)
print("saved lowband_T.json")
