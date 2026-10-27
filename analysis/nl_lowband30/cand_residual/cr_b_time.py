"""(b) time structure: local error / local signal in 25 ms bins (SPAN 0.5, window from the front-frame zero), 20-100 / 100-400 / 400-1000 Hz, per EARLY; plus
the deterministic part (error averaged over the LATE repeats of a setting, noise-reduced)."""
from cr_common import *
B = bench_all()
EAR = (1.0, 2.0, 3.0, 10.0)
BD = [(20, 100), (100, 400), (400, 1000), (1000, 5000)]
SOS = {b: band_sos(*b) for b in BD}
def db(x): return 10 * np.log10(x)
W = int(0.025 * FS)
for e_ in EAR:
    ks = [k for k in sorted(B) if k[0] == 0.5 and k[1] == e_ and not B[k]["trimmed"]]
    ks = [k for k in ks if abs(B[k]["o"] - B[ks[0]]["o"]) < 0.02 and abs(B[k]["front_shift"] - B[ks[0]]["front_shift"]) < 0.3]
    print(f"\nSPAN 0.5 EARLY {e_}: {len(ks)} captures (same session / o)")
    z = K0 + int(B[ks[0]]["front_shift"])
    Y = {}; E = {}
    for k in ks:
        y, m = cached(k, B[k])
        for bd in BD:
            Y[(k, bd)] = sosfiltfilt(SOS[bd], y, axis=0); E[(k, bd)] = sosfiltfilt(SOS[bd], y - m, axis=0)
    print("   t(ms): " + " ".join(f"{1000*(t0)/FS:6.0f}" for t0 in range(int(0.0 * FS), int(0.45 * FS), W)[1:]))
    for bd in BD:
        eavg = np.mean([E[(k, bd)] for k in ks], axis=0)                     # deterministic part (noise averaged down)
        line_t, line_d, line_n = [], [], []
        for t0 in range(int(0.0 * FS), int(0.45 * FS), W)[1:]:
            a, b = z + t0, z + t0 + W
            sig = np.mean([fpow(Y[(k, bd)][a:b]) for k in ks])
            tot = np.mean([fpow(E[(k, bd)][a:b]) for k in ks])
            det = fpow(eavg[a:b]) * (1 - 1 / len(ks)) - 0       # |mean e|^2 includes noise/N
            line_t.append(db(tot / sig)); line_d.append(db(max(fpow(eavg[a:b]) / sig, 1e-9)))
        print(f"  {bd[0]:>4}-{bd[1]:<4} total   " + " ".join(f"{v:6.1f}" for v in line_t))
        print(f"            mean-e  " + " ".join(f"{v:6.1f}" for v in line_d))
