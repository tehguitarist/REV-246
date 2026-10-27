"""Fine-resolution pooled transfer T(f) = <Y M*>/<|M|^2> per FFT bin (Hann window 0.4 s -> ~3-bin smoothing), 20-2000 Hz, estimated from the other EARLY values
and applied to the held-out EARLY (SPAN 0.5, session fs 23, LATE <= 0.4).  Shows whether the 100 Hz-1 kHz error is a fine-structured LTI error (echo / comb)."""
from lb_common import *
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
keys = sorted(k for k in B if k[0] == 0.5 and 0.7 <= k[1] <= 2.0 and k[2] <= 0.4 and abs(B[k]["front_shift"] - 22.94) < .3 and not B[k]["trimmed"])
print(len(keys), "captures", [(k[1], k[2], round(B[k]["o"], 2)) for k in keys][:30])
S = {}
for k in keys:
    r = B[k]; z = K0 + r["front_shift"]; a, b = int(z + 0.05 * FS), int(z + 0.45 * FS)
    y, m, _ = load(k, r, n=b + 64, band=False)
    w = np.hanning(b - a)[:, None]
    S[k] = (np.fft.rfft(y[a:b] * w, NF, axis=0), np.fft.rfft(m[a:b] * w, NF, axis=0))
BD = [(20, 100), (100, 300), (300, 1000), (1000, 2000)]
for held in (0.7, 1.0, 1.5, 2.0):
    tr = [k for k in keys if k[1] != held]; te = [k for k in keys if k[1] == held]
    out = []
    for c in (0, 1):
        T = sum(S[k][0][:, c] * S[k][1][:, c].conj() for k in tr) / (sum(abs(S[k][1][:, c]) ** 2 for k in tr) + 1e-30)
        for lo, hi in BD:
            mk = (fr >= lo) & (fr < hi)
            for tag, TT in (("raw", 1.0), ("fineT", T[mk])):
                v = [10 * np.log10(np.sum(abs(S[k][0][mk, c] - TT * S[k][1][mk, c]) ** 2) / np.sum(abs(S[k][0][mk, c]) ** 2)) for k in te]
                out.append((c, lo, hi, tag, np.median(v)))
    print(f"held out EARLY {held}: " + " | ".join(f"{lo}-{hi}: " + "/".join(f"{np.mean([o[4] for o in out if o[1]==lo and o[3]==t]):6.1f}" for t in ("raw", "fineT")) for lo, hi in BD))
