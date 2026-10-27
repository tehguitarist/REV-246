"""(a) non-parametric sub-18 Hz T on top of the shipped base: per session / SPAN / file-length group, bootstrap sd, equality z-scores. v = window variant."""
from solve import *
import collections
D = pickle.load(open(os.path.join(HERE, "g2_all.pkl"), "rb")); v = sys.argv[1] if len(sys.argv) > 1 else "w60"; lam = 1e-2
# total T = shipped table x (1+delta): add the shipped table at the knots
Ts = LJ.response(KNOTS, 2)           # (nk, 2)
rng = np.random.default_rng(0)
KI = [i for i in range(len(KNOTS)) if KNOTS[i] <= 18]; cols = np.array([2 * i + j for i in KI for j in (0, 1)])
def est(ks, ch):
    G, b, yy, n = pool(D, ks, v, ch); Gs = G[np.ix_(cols, cols)]
    c = np.linalg.solve(Gs + lam * np.mean(np.diag(Gs)) * np.eye(len(cols)), b[cols]); return 1 + c[0::2] + 1j * c[1::2]
def boot(ks, ch, nb=40):
    return np.array([est([ks[i] for i in rng.integers(0, len(ks), len(ks))], ch) for _ in range(nb)])
def show(lab, ks):
    for ch in (0, 1):
        T = est(ks, ch) * Ts[KI, ch]; bs = boot(ks, ch) * Ts[KI, ch]; sd = np.sqrt(np.mean(abs(bs - bs.mean(0)) ** 2, 0))
        print(f"{lab:30s} ch{ch} n{len(ks):3d} " + " ".join(f"{KNOTS[i]:.0f}:{20*np.log10(abs(T[j])):+5.1f}/{np.degrees(np.angle(T[j])):+4.0f}(±{sd[j]:.2f})" for j, i in enumerate(KI) if KNOTS[i] in (6, 8, 10, 12, 14, 16, 18)))
    return
keys = sorted(D); print(collections.Counter((D[k]["sess"]) for k in keys), "variant", v)
print("T_total = shipped x (1+delta); (±) = bootstrap sd of |T| linear")
for s in sorted(set(D[k]["sess"] for k in keys)):
    ks = [k for k in keys if D[k]["sess"] == s]
    if len(ks) >= 20: show(f"fs{s}", ks)
print("--- by SPAN within session")
for s in (0, 20, 23, 16, 36):
    for sp in sorted(set(k[0] for k in keys if D[k]["sess"] == s)):
        ks = [k for k in keys if D[k]["sess"] == s and k[0] == sp]
        if len(ks) >= 20: show(f"fs{s} SPAN{sp}", ks)
print("--- by file length (LATE) within session")
for s in (0, 20, 23, 16, 36):
    for lab, f in (("LATE<=0.2", lambda k: k[2] <= 0.2), ("LATE 0.4-0.8", lambda k: 0.4 <= k[2] <= 0.8), ("LATE>=1.0", lambda k: k[2] >= 1.0)):
        ks = [k for k in keys if D[k]["sess"] == s and f(k)]
        if len(ks) >= 15: show(f"fs{s} {lab}", ks)
print("--- by EARLY")
for lab, f in (("EARLY 0.3-1", lambda k: k[1] <= 1.0), ("EARLY 1.5-3", lambda k: 1.0 < k[1] <= 3.0), ("EARLY 10", lambda k: k[1] >= 10)):
    show(lab, [k for k in keys if f(k)])
