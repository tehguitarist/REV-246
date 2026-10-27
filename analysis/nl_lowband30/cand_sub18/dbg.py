from solve import *
D = pickle.load(open(os.path.join(HERE, "g2_all.pkl"), "rb"))
bad = []
for k, d in D.items():
    for ch in (0, 1):
        g = d["g"]["ls"].get("w60", {}).get(ch)
        if g is None: bad.append((k, "missing")); break
        if not (np.isfinite(g[0]).all() and np.isfinite(g[1]).all()) or g[4] < 1e-12: bad.append((k, "nonfinite/sc", g[4])); break
        b = d["g"]["band"][("20-100", ch)]
        if not np.isfinite(b["bb_raw"]).all(): bad.append((k, "bandnan")); break
print(len(D), len(bad)); print(bad[:10])
KN = [i for i in range(len(KNOTS)) if KNOTS[i] <= 18]; cols = np.array([2*i+j for i in KN for j in (0,1)])
import collections
print(collections.Counter(k[0] for k in D))
for sp in sorted(set(k[0] for k in D)):
    ks = [q for q in D if q[0] != sp]; G, b, yy, n = pool(D, ks, "w60", 0)
    Gs = G[np.ix_(cols, cols)]; print(sp, len(ks), np.isfinite(Gs).all(), np.mean(np.diag(Gs)), np.linalg.cond(Gs + 1e-2*np.mean(np.diag(Gs))*np.eye(len(cols))))
