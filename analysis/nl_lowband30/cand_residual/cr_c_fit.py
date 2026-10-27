"""(c) Regression of the residual e = y - m on derivative columns of candidate perturbations, from the per-capture sufficient statistics (stats2/).
Per channel parameters (reported count = per channel x 2). Held out by EARLY (SPAN 0.5, leave one EARLY out) and by SPAN (0.5 <-> 0.2/0.3/0.4).
Metric: median over test captures x channels of 10 log10(|e - X beta|^2 / |y|^2) in the band (4th-order band-pass, window 0.05-0.45 s or to the last tap + 20 ms)."""
from cr_common import *
import glob
B = bench_all()
files = sorted(glob.glob(os.path.join(HERE, "stats2", "*.npz")))
D = {}
for f in files:
    k = tuple(float(v) for v in os.path.basename(f)[:-4].split("_"))
    z = np.load(f); D[k] = {n: z[n] for n in z.files}
names = list(D[next(iter(D))]["names"]); P = len(names); ix = {n: i for i, n in enumerate(names)}
NUIS = [ix["nuis_gain"], ix["nuis_shift"], ix["nuis_ppm"]]
def db(x): return 10 * np.log10(x)

def lin(*groups):
    """transform matrix G (P x q) from lists of (name list); each group is a list of column names summed into one parameter."""
    G = np.zeros((P, len(groups)))
    for j, g in enumerate(groups):
        for n in g: G[ix[n], j] += 1.0
    return G
def stage_of(i): return [j for j, st in enumerate(ne.STAGES) if i in st][0]
tapn = lambda kind: [f"tap{i}_{kind}" for i in range(27)]
stageG = lambda kind, pre="tap": lin(*[[f"{pre}{i}_{kind}" if pre == "tap" else f"bk_tap{i}" for i in st] for st in ne.STAGES])
CAND = {}
CAND["branch_gain (1)"] = lin(["branch_gain"])
CAND["branch_dt (1)"] = lin(["branch_dt"])
CAND["branch gain+dt (2)"] = lin(["branch_gain"], ["branch_dt"])
CAND["branch drift: dt, gain ramps (2)"] = lin(["branch_dt_ramp"], ["branch_gain_ramp"])
CAND["front gain+dt+rdir_dt (3)"] = lin(["front_gain"], ["front_dt"], ["rdir_dt"])
CAND["front taps gain+dt (6)"] = lin(*[[f"fr{j}_gain"] for j in range(3)], *[[f"fr{j}_dt"] for j in range(3)])
CAND["per-tap gain (27)"] = lin(*[[n] for n in tapn("gain")])
CAND["per-tap dt (27)"] = lin(*[[n] for n in tapn("dt")])
CAND["per-tap gain+dt (54)"] = lin(*[[n] for n in tapn("gain") + tapn("dt")])
CAND["per-stage gain (14)"] = stageG("gain")
CAND["per-stage dt (14)"] = stageG("dt")
CAND["BK alpha, all taps (1)"] = lin(["bk_alpha_all"])
CAND["BK alpha per stage (14)"] = stageG("bk", pre="bk")
CAND["read shape +-14 free (29)"] = lin(*[[f"shape{j:+d}"] for j in range(-14, 15)])
CAND["AP124, APout g (2)"] = lin(["ap124"], ["apout"])
allc = [n for n in names if not n.startswith("nuis")]
CAND["ALL columns (%d)" % len(allc)] = lin(*[[n] for n in allc])

def reduce(d, ch, bd, G, nuis):
    A = d[f"XtX_{ch}_{bd}"]; b = d[f"Xte_{ch}_{bd}"]; c = d[f"ete_{ch}_{bd}"]
    if nuis:
        Nn = NUIS; Pp = [i for i in range(P) if i not in Nn]
        Ann = A[np.ix_(Nn, Nn)]; Ainv = np.linalg.pinv(Ann, rcond=1e-12)
        Apn = A[np.ix_(Pp, Nn)]
        A2 = np.zeros_like(A); b2 = np.zeros_like(b)
        A2[np.ix_(Pp, Pp)] = A[np.ix_(Pp, Pp)] - Apn @ Ainv @ Apn.T
        b2[Pp] = b[Pp] - Apn @ Ainv @ b[Nn]
        c2 = c - b[Nn] @ Ainv @ b[Nn]
        A, b, c = A2, b2, c2
    return G.T @ A @ G, G.T @ b, c
def fold_scores(train, test, G, bd, nuis, lam=1e-7):
    out = []
    for ch in (0, 1):
        q = G.shape[1]
        if q:
            Aw = np.zeros((q, q)); bw = np.zeros(q)
            for k in train:
                A, b, c = reduce(D[k], ch, bd, G, nuis); w = 1.0 / D[k][f"yty_{ch}_{bd}"]
                Aw += w * A; bw += w * b
            beta = np.linalg.solve(Aw + lam * np.trace(Aw) / q * np.eye(q), bw)
        for k in test:
            A, b, c = reduce(D[k], ch, bd, G, nuis)
            e = c - (2 * beta @ b - beta @ A @ beta if q else 0.0)
            out.append(db(max(e, 1e-30) / D[k][f"yty_{ch}_{bd}"]))
    return np.array(out)
E5 = [k for k in D if k[0] == 0.5]; EL = sorted({k[1] for k in E5})
def run(G, bd, nuis):
    """held-out scores: dict fold -> median dB, plus pooled over (leave-EARLY-out) and (SPAN) test sets"""
    res = {}
    sc = []
    for e_ in EL:
        tr = [k for k in E5 if k[1] != e_]; te = [k for k in E5 if k[1] == e_]
        s = fold_scores(tr, te, G, bd, nuis); res[f"E{e_:g}"] = np.median(s); sc.append(s)
    res["EARLY-pooled"] = np.median(np.concatenate(sc))
    lo = [k for k in D if k[0] < 0.5]
    res["SPAN .5->.2-.4"] = np.median(fold_scores(E5, lo, G, bd, nuis))
    res["SPAN .2-.4->.5"] = np.median(fold_scores(lo, E5, G, bd, nuis))
    return res
if __name__ == "__main__":
    bands = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [100, 20, 1000]
    Gnone = np.zeros((P, 0))
    for nuis in (False, True):
        for bd in bands:
            print(f"\n=== band {bd}-{ {20:100,100:1000,1000:5000}[bd] } Hz, nuisance (gain, shift, ppm per capture) {'FREE' if nuis else 'fixed at the bench alignment'}; {len(D)} captures ===")
            base = run(Gnone, bd, nuis)
            def line(nm, r):
                return f"{nm:36s} " + " ".join(f"{r[k]:6.1f}" for k in [f"E{e:g}" for e in EL] + ["EARLY-pooled", "SPAN .5->.2-.4", "SPAN .2-.4->.5"])
            print(f"{'':36s} " + " ".join(f"{k:>6s}" for k in [f"E{e:g}" for e in EL] + ["Epool", "S>lo", "S>.5"]))
            print(line("baseline (no parameters)", base))
            for nm, G in CAND.items():
                r = run(G, bd, nuis)
                print(line(nm, r) + f"   d(EARLY-pool) {r['EARLY-pooled']-base['EARLY-pooled']:+5.2f}  d(S>lo) {r['SPAN .5->.2-.4']-base['SPAN .5->.2-.4']:+5.2f}  d(S>.5) {r['SPAN .2-.4->.5']-base['SPAN .2-.4->.5']:+5.2f}", flush=True)
