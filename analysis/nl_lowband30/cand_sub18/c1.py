"""(c) 20-100 Hz null of the shipped model (lb2 alignment), raw and step-subtracted, and with a sub-18 Hz table on top (exact, from band normal equations).
Held out: table for capture k estimated WITHOUT that capture's SPAN (leave-one-SPAN-out, pooled over sessions) or without that capture's session-and-SPAN group."""
from solve import *
D = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb")); v = sys.argv[2] if len(sys.argv) > 2 else "w60"
lam = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-2
KN = [i for i in range(len(KNOTS)) if KNOTS[i] <= 18]           # sub-18 Hz knots (2..18)
cols = np.array([2 * i + j for i in KN for j in (0, 1)])
def db(x, y): return 10 * np.log10(max(x, 1e-30) / y)
keys = sorted(D)
def nullset(k, c, ch, cor):
    b = D[k]["g"]["band"][("20-100", ch)]; ee = b["ee_cor" if cor else "ee_raw"]; bb = b["bb_cor" if cor else "bb_raw"]
    cc = np.zeros(NC); 
    if c is not None: cc[cols] = c
    # band-null scaling: the LS columns were scaled by sc per capture; the band columns are NOT scaled (unit gain), so c must be rescaled
    return db(ee - 2 * cc @ bb + cc @ b["Gb"] @ cc, b["yy"])
def fit(ks, ch):
    G, b, yy, n = pool(D, ks, v, ch)
    # coefficient c refers to unscaled columns after dividing by sc; the LS used res/sc and cols/sc, so c is dimensionless and shared: OK
    Gs = G[np.ix_(cols, cols)]; return np.linalg.solve(Gs + lam * np.mean(np.diag(Gs)) * np.eye(len(cols)), b[cols])
def summ(name, f):
    v_ = np.array([[f(k, ch) for ch in (0, 1)] for k in keys]); print(f"{name:55s} median L/R {np.median(v_[:,0]):6.1f} / {np.median(v_[:,1]):6.1f}   mean {v_.mean():6.1f}  n={len(keys)}")
    return v_
print("captures", len(keys), "variant", v, "lam", lam)
raw = summ("shipped lb2, raw", lambda k, ch: nullset(k, None, ch, False))
cor = summ("shipped lb2, pre-roll step subtracted (2 numbers/capture)", lambda k, ch: nullset(k, None, ch, True))
print("step effect: median change %.2f dB, better>0.5 dB in %d of %d" % (np.median(cor - raw), np.sum((cor - raw) < -0.5), cor.size))
# oracle per capture (in-sample upper bound; 2*len(KN) params per capture and channel)
def per_cap(k, ch, cor_):
    b = D[k]["g"]["band"][("20-100", ch)]; bb = b["bb_cor" if cor_ else "bb_raw"][cols]; G = b["Gb"][np.ix_(cols, cols)]
    c = np.linalg.solve(G + 1e-6 * np.mean(np.diag(G)) * np.eye(len(cols)), bb); return nullset(k, c, ch, cor_)
summ("ORACLE per-capture sub-18 table (in-sample, %d params), raw" % len(cols), lambda k, ch: per_cap(k, ch, False))
summ("ORACLE per-capture sub-18 table (in-sample), step-subtracted", lambda k, ch: per_cap(k, ch, True))
# held out tables
for cor_ in (False, True):
    for grp, lab in ((lambda k: D[k]["r"] and k[0], "leave-one-SPAN-out"), (lambda k: (D[k]["sess"], k[0]), "leave-one-(session,SPAN)-out"), (lambda k: k, "leave-one-capture-out (pooled)"), (lambda k: D[k]["sess"], "leave-one-session-out")):
        cache = {}
        def f(k, ch):
            g = grp(k); key = (g, ch)
            if key not in cache: cache[key] = fit([q for q in keys if grp(q) != g], ch)
            return nullset(k, cache[key], ch, cor_)
        if lab.startswith("leave-one-capture"):
            # cheap: pooled all minus this capture
            Gt, bt = [], []
            for ch in (0, 1):
                G, b, yy, n = pool(D, keys, v, ch); Gt.append(G); bt.append(b)
            def f(k, ch):
                g = D[k]["g"]["ls"][v][ch]; G = Gt[ch] - g[0]; b = bt[ch] - g[1]
                Gs = G[np.ix_(cols, cols)]; c = np.linalg.solve(Gs + lam * np.mean(np.diag(Gs)) * np.eye(len(cols)), b[cols]); return nullset(k, c, ch, cor_)
        summ(f"held out {lab}, {'step-sub' if cor_ else 'raw'} (+{len(cols)} params shared)", f)
