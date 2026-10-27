"""Is the first-25 ms residual e the same in captures that share the front section? Correlation of e (front frame, band) between captures of one EARLY,
within and across sessions (front-shift clusters). Captures with a branch offset o not near 0/1 excluded (first main tap inside the window)."""
import sys, itertools
from lb_common import *

B = bench_all()
band = {"20-300": bands(20, 300), "100-1k": bands(100, 1000)}
def collect(early, lates=(0.0, 0.2, 0.4, 0.6)):
    out = []
    for k, r in sorted(B.items()):
        if k[1] != early or k[2] not in lates or r["trimmed"] or abs(r["o"] - 0.5) > 1.2: continue
        y, m, _ = load(k, r, n=3000)
        yf, mf = frame(y, r["front_shift"]), frame(m, r["front_shift"])
        pk = np.abs(mf[K0:K0 + 800, 0]).max()
        out.append((k, r["front_shift"], (yf - mf) / pk))
    return out
ia, ib = K0 - 40, K0 + int(0.025 * FS)
for early in (10.0, 3.0):
    caps = collect(early)
    print(f"\nEARLY {early}: {len(caps)} captures (LATE<=0.6, o in 0..1.7):")
    for bn, sos in band.items():
        for c in (0, 1):
            E = np.array([sosfiltfilt(sos, e[:, c])[ia:ib] for _, _, e in caps])
            names = [f"{k[0]}/{k[2]}@{fs:.0f}" for k, fs, _ in caps]
            amp = 20 * np.log10(np.sqrt((E ** 2).mean(1)))
            C = np.corrcoef(E)
            same = [C[i, j] for i, j in itertools.combinations(range(len(E)), 2) if abs(caps[i][1] - caps[j][1]) < 1.0]
            diff = [C[i, j] for i, j in itertools.combinations(range(len(E)), 2) if abs(caps[i][1] - caps[j][1]) >= 1.0]
            print(f"  {bn} ch{c}: e rms dB re peak  min {amp.min():.1f} med {np.median(amp):.1f} max {amp.max():.1f};  corr(e_i,e_j) same fs-cluster: mean {np.mean(same):.2f} (n{len(same)}) | different cluster: mean {np.mean(diff):.2f} (n{len(diff)})")
            # mean e per cluster and similarity of cluster means
            cl = {}
            for (k, fs, _), e in zip(caps, E): cl.setdefault(round(fs), []).append(e)
            mm = {f: np.mean(v, 0) for f, v in cl.items()}
            fl = sorted(mm)
            print("     cluster-mean e rms (dB) :", {f: round(20 * np.log10(np.sqrt((mm[f] ** 2).mean())), 1) for f in fl}, " n:", {f: len(cl[f]) for f in fl})
            print("     corr of cluster means   :", np.round(np.corrcoef([mm[f] for f in fl]), 2).tolist() if len(fl) > 1 else "")
