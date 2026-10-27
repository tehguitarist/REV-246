"""Joint N + S low-band transfer.  T from N alone, S alone, or both (equal program weight), each estimated WITHOUT the test captures
(N: leave-one-SPAN-out; S: other RT fold, LF >= 1.0 as sml_fullband.fold_T), scored on both programs from the stored cross sums
(|Y - T M|^2 = |Y|^2 + |T|^2 |M|^2 - 2 Re(T conj(Y M*)), per 2 Hz bin, summed in bands), median over captures / settings.   usage: lb_joint.py [smooth=2] [wN=1]"""
import sys, glob, os, pickle
import numpy as np
SM = int(sys.argv[1]) if len(sys.argv) > 1 else 2; WN = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
TBIN, TNB = 2.0, 6150; fc = (np.arange(TNB) + 0.5) * TBIN
FOLDS = {"A": (0.3, 0.8, 1.4, 2.2), "B": (0.6, 1.0, 1.8, 2.5)}
S = {}
for p in sorted(glob.glob("/Users/leighpierce/Git/246/analysis/fullband/s/*.npz")):
    k = os.path.basename(p)[:-4]; lf, rt, hf = (float(x) for x in k.split("_")); z = np.load(p)
    S[k] = dict(lf=lf, rt=rt, sym=z["sym"].astype(complex), smm=z["smm"].astype(float), syy=z["syy"].astype(float))
N = {k: dict(span=v[0], sym=v[1], smm=v[2], syy=v[3]) for k, v in pickle.load(open("joint_N_sums.pkl", "rb")).items()}
print(len(S), "S settings,", len(N), "N captures; smoothing +-%d bins, N weight %.2f" % (SM, WN))
def tot(ds): return (sum(d["sym"] for d in ds), sum(d["smm"] for d in ds))
band = (fc >= 20) & (fc < 1600)
def norm_dt(smm): return np.array([smm[c][band].sum() for c in (0, 1)])[:, None]
def sm(x):
    ker = np.ones(2 * SM + 1)
    return np.stack([np.convolve(x[c], ker, "same") for c in (0, 1)])
def T_from(parts):
    """parts: list of (sym, smm, weight); each program normalised by its total smm in 20-1600 Hz before pooling."""
    num = 0; den = 0
    for sym, smm, w in parts:
        d = norm_dt(smm); num = num + w * sym / d; den = den + w * smm / d
    return sm(num) / np.maximum(sm(den), 1e-30)
S_keys = lambda fold=None: [k for k, v in S.items() if v["lf"] >= 1.0 and (fold is None or v["rt"] in FOLDS[fold])]
BD = [(20, 40), (40, 80), (80, 160), (160, 270), (270, 400), (400, 800), (800, 1600), (1600, 3200)]
def err(d, T, m, c):
    t = T[c][m]; return np.sum(d["syy"][c][m] + abs(t) ** 2 * d["smm"][c][m] - 2 * np.real(t * np.conj(d["sym"][c][m])))
one = np.ones((2, TNB), complex)
res = {"S": {}, "N": {}}
# ---- S held out by RT fold
for name in ("none", "N", "S", "joint"):
    for k, v in S.items():
        other = "B" if v["rt"] in FOLDS["A"] else "A"
        sS = tot([S[q] for q in S_keys(other)]); sN = tot(list(N.values()))
        T = {"none": one, "N": T_from([(*sN, 1.0)]), "S": T_from([(*sS, 1.0)]), "joint": T_from([(*sS, 1.0), (*sN, WN)])}[name]
        for c in (0, 1):
            for bi, (lo, hi) in enumerate(BD):
                m = (fc >= lo) & (fc < hi); den = v["syy"][c][m].sum()
                res["S"].setdefault((name, c, bi), []).append(10 * np.log10(max(err(v, T, m, c), 1e-30) / den))
# ---- N held out by SPAN
Sall = tot([S[q] for q in S_keys()])
for name in ("none", "N", "S", "joint"):
    for sp in (0.2, 0.3, 0.4, 0.5):
        sN = tot([d for d in N.values() if d["span"] != sp])
        T = {"none": one, "N": T_from([(*sN, 1.0)]), "S": T_from([(*Sall, 1.0)]), "joint": T_from([(*Sall, 1.0), (*sN, WN)])}[name]
        for k, d in N.items():
            if d["span"] != sp: continue
            for c in (0, 1):
                for bi, (lo, hi) in enumerate(BD):
                    m = (fc >= lo) & (fc < hi); den = d["syy"][c][m].sum()
                    res["N"].setdefault((name, c, bi), []).append(10 * np.log10(max(err(d, T, m, c), 1e-30) / den))
for prog in ("N", "S"):
    print(f"\n=== scored on {prog}  ({'leave-one-SPAN-out' if prog=='N' else 'other RT fold'}); median band error, dB re |Y|^2, L / R   [T estimated from: none | N | S | joint N+S]")
    for bi, (lo, hi) in enumerate(BD):
        print(f"  {lo:5d}-{hi:<5d} " + "  ".join(f"{nm:>5s} {np.median(res[prog][(nm,0,bi)]):6.1f}/{np.median(res[prog][(nm,1,bi)]):6.1f}" for nm in ("none", "N", "S", "joint")))
pickle.dump(res, open("joint_res.pkl", "wb"))
