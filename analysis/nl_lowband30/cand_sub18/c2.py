"""(c) continued: sub-18 Hz table fitted on the 20-100 Hz band objective (sum of per-capture relative band errors), held out by group. Exact scoring from band normal equations."""
from solve import *
D = pickle.load(open(os.path.join(HERE, "g2_all.pkl"), "rb"))
lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-2; kmax = float(sys.argv[2]) if len(sys.argv) > 2 else 18
KN = [i for i in range(len(KNOTS)) if KNOTS[i] <= kmax]; cols = np.array([2 * i + j for i in KN for j in (0, 1)])
keys = sorted(D); db = lambda x, y: 10 * np.log10(max(x, 1e-30) / y)
def bd(k, ch): return D[k]["g"]["band"][("20-100", ch)]
def null(k, ch, c, cor):
    b = bd(k, ch); cc = np.zeros(NC)
    if c is not None: cc[cols] = c
    bb = b["bb_cor" if cor else "bb_raw"]; return db(b["ee_cor" if cor else "ee_raw"] - 2 * cc @ bb + cc @ b["Gb"] @ cc, b["yy"])
def fit(ks, ch, cor):
    G = 0; r = 0
    for k in ks:
        b = bd(k, ch); w = 1 / b["yy"]; G = G + w * b["Gb"][np.ix_(cols, cols)]; r = r + w * b["bb_cor" if cor else "bb_raw"][cols]
    return np.linalg.solve(G + lam * np.mean(np.diag(G)) * np.eye(len(cols)), r)
groups = {"leave-one-SPAN-out": lambda k: k[0], "leave-one-(session,SPAN)-out": lambda k: (D[k]["sess"], k[0]), "leave-one-session-out": lambda k: D[k]["sess"],
          "leave-one-EARLY-out": lambda k: k[1]}
print("lam", lam, "knots <=", kmax, "params", len(cols), "per channel")
for cor in (False, True):
    base = np.array([[null(k, ch, None, cor) for ch in (0, 1)] for k in keys])
    print(("step-subtracted" if cor else "raw"), "shipped: median L/R %.1f / %.1f" % tuple(np.median(base, 0)))
    for lab, grp in groups.items():
        cache = {}
        def f(k, ch):
            g = (grp(k), ch)
            if g not in cache: cache[g] = fit([q for q in keys if grp(q) != grp(k)], ch, cor)
            return null(k, ch, cache[g], cor)
        v_ = np.array([[f(k, ch) for ch in (0, 1)] for k in keys]); ch_ = v_ - base
        print(f"   {lab:30s} median L/R {np.median(v_[:,0]):6.1f} / {np.median(v_[:,1]):6.1f}   change median {np.median(ch_):+.2f}  worse>0.5: {np.sum(ch_>0.5)}/{ch_.size}  better>0.5: {np.sum(ch_<-0.5)}")
    # in-sample pooled
    c = [fit(keys, ch, cor) for ch in (0, 1)]
    v_ = np.array([[null(k, ch, c[ch], cor) for ch in (0, 1)] for k in keys]); print(f"   {'in-sample pooled':30s} median L/R {np.median(v_[:,0]):6.1f} / {np.median(v_[:,1]):6.1f}")
    # by subgroup of the shipped null: which captures reach -50
    print("   shipped: fraction of channel-windows <= -50 dB: %.2f" % np.mean(base <= -50))
