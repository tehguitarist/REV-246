import json, numpy as np, collections
d = json.load(open("oracle.json"))
g = collections.defaultdict(lambda: collections.defaultdict(list))
for k, v in d.items():
    s, e, l = (float(x) for x in k.split("_"))
    ses = (s, round(v["fs"]) if abs(v["fs"]) > 0.5 else 0)
    g[ses][l].append((v["c0"] * 1e4, v["c1"] * 1e4, e))
print("pre-roll step c (1e-4 re 410 peak), L / R: median over EARLY>=0.3 captures by LATE; spread = std over EARLY")
Ls = [0.0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]
print("SPAN fs   " + "".join(f"  LATE {l:3.1f}     " for l in Ls))
for ses in sorted(g):
    row = f"{ses[0]:3.1f} {ses[1]:3d}  "
    for l in Ls:
        a = np.array(g[ses].get(l, [(np.nan,) * 3]))
        row += f" {np.nanmedian(a[:,0]):5.1f}/{np.nanmedian(a[:,1]):5.1f}  "
    print(row)
# correlations between L and R c, and a session-level model c = a_s + b_s*LATE
print("\nper session: corr(c_L, c_R); R slope vs LATE (1e-4 per unit LATE); std of c_R over captures at fixed LATE (mean) ")
for ses in sorted(g):
    allc = [(l, c0, c1) for l, L in g[ses].items() for c0, c1, e in L]
    A = np.array(allc)
    sl = np.polyfit(A[:, 0], A[:, 2], 1)[0]
    within = np.mean([np.std([c1 for _, c1, _ in L]) for l, L in g[ses].items() if len(L) > 3])
    print(f"  SPAN {ses[0]} fs {ses[1]:3d}: corr(L,R) {np.corrcoef(A[:,1],A[:,2])[0,1]:5.2f}  R slope {sl:5.2f}  within-LATE std R {within:4.1f}  (n {len(A)})")
