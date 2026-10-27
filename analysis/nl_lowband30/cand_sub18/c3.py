"""Is a per-session (or per-(session,SPAN), per-LATE-class) sub-18 table better than the pooled one? leave-one-capture-out inside the group, band objective."""
import os, sys
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "c2.py")).read().split("groups =")[0])
import itertools
print("lam", lam, "params", len(cols))
for cor in (True,):
    base = np.array([[null(k, ch, None, cor) for ch in (0, 1)] for k in keys])
    print("step-sub shipped: median L/R %.1f / %.1f" % tuple(np.median(base, 0)))
    defs = {"pooled (all others)": lambda k: 0, "own session": lambda k: D[k]["sess"], "own (session,SPAN)": lambda k: (D[k]["sess"], k[0]),
            "own (session,SPAN,LATE>0)": lambda k: (D[k]["sess"], k[0], k[2] > 0), "own (session,SPAN,LATE)": lambda k: (D[k]["sess"], k[0], k[2]),
            "own (session,SPAN,EARLY)": lambda k: (D[k]["sess"], k[0], k[1])}
    for lab, grp in defs.items():
        out = np.full((len(keys), 2), np.nan); cnt = []
        # sum G, r once per group, subtract own capture
        for ch in (0, 1):
            tot = {}
            for k in keys:
                b = bd(k, ch); w = 1 / b["yy"]; g = grp(k)
                G0, r0, n0 = tot.get(g, (0, 0, 0)); tot[g] = (G0 + w * b["Gb"][np.ix_(cols, cols)], r0 + w * b["bb_cor"][cols], n0 + 1)
            for i, k in enumerate(keys):
                G0, r0, n0 = tot[grp(k)]
                if n0 < 4: out[i, ch] = base[i, ch]; continue
                b = bd(k, ch); w = 1 / b["yy"]; G = G0 - w * b["Gb"][np.ix_(cols, cols)]; r = r0 - w * b["bb_cor"][cols]
                c = np.linalg.solve(G + lam * np.mean(np.diag(G)) * np.eye(len(cols)), r); out[i, ch] = null(k, ch, c, cor)
        d = out - base
        print(f"   {lab:30s} median L/R {np.median(out[:,0]):6.1f} / {np.median(out[:,1]):6.1f}  change {np.median(d):+.2f}  worse>0.5 {np.sum(d>0.5)}/{d.size}  better>0.5 {np.sum(d<-0.5)}   groups {len(tot)}")
