"""Cross-session / cross-SPAN test of the fine-resolution transfer T(f): estimated on SPAN 0.5, session fs 23, EARLY 0.7-2, LATE <= 0.4, applied (as an LTI filter on the model
record) to other SPANs and sessions.  Metric: time-domain 4th-order band nulls over front-zero + 0.05 s ... min(0.45 s, late feed - 10 ms).  T = 1 outside 15 Hz-3 kHz."""
import pickle
from lb_common import *
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
def spec_of(k, t0=0.05, t1=0.45):
    r = B[k]; z = K0 + r["front_shift"]; a, b = int(z + t0 * FS), int(z + t1 * FS)
    y, m, _ = load(k, r, n=b + 64, band=False); w = np.hanning(b - a)[:, None]
    return np.fft.rfft(y[a:b] * w, NF, axis=0), np.fft.rfft(m[a:b] * w, NF, axis=0)
tr = sorted(k for k in B if k[0] == 0.5 and 0.7 <= k[1] <= 2.0 and k[2] <= 0.4 and abs(B[k]["front_shift"] - 22.94) < .3 and not B[k]["trimmed"])
S = {k: spec_of(k) for k in tr}
T = [sum(S[k][0][:, c] * S[k][1][:, c].conj() for k in tr) / (sum(abs(S[k][1][:, c]) ** 2 for k in tr) + 1e-30) for c in (0, 1)]
sel = (fr >= 15) & (fr <= 3000); wgt = np.clip(np.minimum((fr - 10) / 5, (4000 - fr) / 1000), 0, 1)
T = [1 + (t - 1) * wgt for t in T]
pickle.dump(dict(fr=fr, T=T), open("fineT.pkl", "wb"))
BDs = [(20, 100), (100, 300), (300, 1000), (1000, 5000)]
SOS = {b: bands(*b) for b in BDs}
FEED = {0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}
groups = {
    "train set itself (SPAN .5 fs23, in-sample)": [k for k in tr],
    "SPAN .5 fs23 EARLY .3/.5 (other EARLY)": [k for k in B if k[0] == 0.5 and k[1] in (0.3, 0.5) and k[2] <= 0.4 and abs(B[k]["front_shift"] - 22.94) < .3 and not B[k]["trimmed"]],
    "SPAN .5 fs0 EARLY 3/10 (other session)": [k for k in B if k[0] == 0.5 and k[1] in (3.0, 10.0) and k[2] <= 0.4 and abs(B[k]["front_shift"]) < .3 and not B[k]["trimmed"]],
    "SPAN .4 fs20 EARLY>=.7": [k for k in B if k[0] == 0.4 and k[1] >= 0.7 and k[2] <= 0.2 and not B[k]["trimmed"]],
    "SPAN .3 fs0 EARLY>=.7": [k for k in B if k[0] == 0.3 and k[1] >= 0.7 and k[2] <= 0.2 and not B[k]["trimmed"]],
    "SPAN .2 fs23 EARLY>=.7": [k for k in B if k[0] == 0.2 and k[1] >= 0.7 and k[2] <= 0.2 and not B[k]["trimmed"]],
}
for gname, ks in groups.items():
    res = {b: [[], []] for b in BDs}
    for k in ks:
        r = B[k]; z = K0 + r["front_shift"]
        end = 0.45 * FS if k[2] == 0 or k[0] == 0.5 else (FEED[k[0]] + r["o"] - 0.01 * FS)
        if k[2] > 0 and k[0] < 0.5: end = min(end, FEED[k[0]] + r["o"] - 0.01 * FS)
        a, b = int(z + 0.05 * FS), int(z + min(0.45 * FS, end))
        if b - a < int(0.08 * FS): continue
        y, m, _ = load(k, r, n=b + 64, band=False)
        mc = np.stack([np.fft.irfft(np.fft.rfft(m[:, c], NF) * T[c], NF)[:len(m)] for c in (0, 1)], 1)
        for bd in BDs:
            for c in (0, 1):
                yb = sosfiltfilt(SOS[bd], y[:, c])[a:b]
                for j, mm in enumerate((m, mc)):
                    mb = sosfiltfilt(SOS[bd], mm[:, c])[a:b]
                    res[bd][j].append(10 * np.log10(np.sum((yb - mb) ** 2) / np.sum(yb ** 2)))
    print(f"{gname:46s} n{len(ks):3d} | " + " | ".join(f"{b[0]}-{b[1]}: {np.median(res[b][0]):6.1f} -> {np.median(res[b][1]):6.1f}" for b in BDs))
