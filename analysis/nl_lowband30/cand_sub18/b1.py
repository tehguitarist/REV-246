"""(b) estimation-design sensitivity of the shared low-band table (as lb_joint.py: joint N+S, equal program weight unless varied).  Held out: N leave-one-SPAN-out, S other RT fold.
Reports median band error energy (dB re |Y|^2, L/R) for 20-40, 40-80, 80-160 Hz, scored on (T0) the original 51 N captures at w45f and (Tall) every N capture with EARLY>=0.3 / enough signal at w45f.
Compared with the baseline estimation design on the SAME test captures."""
import pickle, glob, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
TBIN, TNB = 2.0, 6150; fc = (np.arange(TNB) + 0.5) * TBIN
FOLDS = {"A": (0.3, 0.8, 1.4, 2.2), "B": (0.6, 1.0, 1.8, 2.5)}
S = {}
for p in sorted(glob.glob("/Users/leighpierce/Git/246/analysis/fullband/s/*.npz")):
    k = os.path.basename(p)[:-4]; lf, rt, hf = (float(x) for x in k.split("_")); z = np.load(p)
    S[k] = dict(lf=lf, rt=rt, sym=z["sym"].astype(complex), smm=z["smm"].astype(float), syy=z["syy"].astype(float))
NN = pickle.load(open(os.path.join(HERE, "bsums.pkl"), "rb")); NN = {k: v for k, v in NN.items() if v}
N0 = pickle.load(open("/Users/leighpierce/Git/246/analysis/nl_lowband30/joint_N_sums.pkl", "rb"))
BD = [(20, 40), (40, 80), (80, 160)]
band = (fc >= 20) & (fc < 1600)
def nd(smm): return np.array([smm[c][band].sum() for c in (0, 1)])[:, None]
def sm(x, SM):
    if SM == 0: return x
    ker = np.ones(2 * SM + 1); return np.stack([np.convolve(x[c], ker, "same") for c in (0, 1)])
def T_from(parts, SM):
    num = 0; den = 0
    for sym, smm, w in parts:
        d = nd(smm); num = num + w * sym / d; den = den + w * smm / d
    return sm(num, SM) / np.maximum(sm(den, SM), 1e-30)
def tot(ds): return (sum(d[0] for d in ds), sum(d[1] for d in ds))
def err(sym, smm, syy, T, m, c):
    t = T[c][m]; return np.sum(syy[c][m] + abs(t) ** 2 * smm[c][m] - 2 * np.real(t * np.conj(sym[c][m])))
def enough(v):
    f = v["floor"].get("0.05-0.15") or [None, None]
    return all(x is not None and x <= -30 for x in f)
def sums(v, win):
    w = v["win"].get(win); return None if w is None else (w[0], w[1], w[2])
# ---- training designs: name -> (selector(key, v) -> window name or None)
DES = {
 "base (SPAN>=.2, LATE<=.2, EARLY>=.3, w45)": lambda k, v: "w45f" if k[0] >= 0.2 and k[2] <= 0.2 and k[1] >= 0.3 else None,
 "+LATE>0.2 (window before late feed)":        lambda k, v: "w45f" if k[0] >= 0.2 and k[1] >= 0.3 else None,
 "+EARLY<0 ramps with enough signal":          lambda k, v: "w45f" if k[0] >= 0.2 and k[2] <= 0.2 and (k[1] >= 0.3 or enough(v)) else None,
 "+SPAN 0.1 (LATE<=.2, pre-feed)":             lambda k, v: "w45f" if k[0] >= 0.1 and k[2] <= 0.2 and k[1] >= 0.3 else None,
 "+SPAN 0.1 incl. LATE>.2 (window NOT limited)": lambda k, v: ("w45n" if k[0] < 0.2 and k[2] > 0.2 else "w45f") if k[0] >= 0.1 and k[1] >= 0.3 else None,
 "window to 0.6 s":                            lambda k, v: "w60f" if k[0] >= 0.2 and k[2] <= 0.2 and k[1] >= 0.3 else None,
 "all of the above (w60f)":                    lambda k, v: "w60f" if k[0] >= 0.1 and (k[1] >= 0.3 or enough(v)) else None,
 "LATE 0 only":                                lambda k, v: "w45f" if k[0] >= 0.2 and k[2] == 0 and k[1] >= 0.3 else None,
 "EARLY>=1.5 only":                            lambda k, v: "w45f" if k[0] >= 0.2 and k[2] <= 0.2 and k[1] >= 1.5 else None,
 "EARLY 0.3-1 only":                           lambda k, v: "w45f" if k[0] >= 0.2 and k[2] <= 0.2 and 0.3 <= k[1] <= 1.0 else None,
}
def train_set(des):
    out = {}
    for k, v in NN.items():
        w = DES[des](k, v)
        if w and sums(v, w) is not None: out[k] = sums(v, w)
    return out
T0 = list(N0)                                                             # original 51
Tall = [k for k, v in NN.items() if k[0] >= 0.2 and k[1] >= 0.3 and "w45f" in v["win"]]
S_keys = lambda fold=None: [k for k, v in S.items() if v["lf"] >= 1.0 and (fold is None or v["rt"] in FOLDS[fold])]
Sall = tot([(S[q]["sym"], S[q]["smm"]) for q in S_keys()])
def evaluate(des, SM=2, WN=1.0, WS=1.0):
    TR = train_set(des); res = {"N0": {}, "Nall": {}, "S": {}}
    for sp in (0.2, 0.3, 0.4, 0.5):
        sN = tot([(a, b) for k, (a, b, c) in TR.items() if k[0] != sp])
        T = T_from([(*Sall, WS), (*sN, WN)], SM)
        for tag, tk in (("N0", T0), ("Nall", Tall)):
            for k in tk:
                if k[0] != sp: continue
                s = NN[k]["win"]["w45f"]; sym, smm, syy = s[0], s[1], s[2]
                for c in (0, 1):
                    for bi, (lo, hi) in enumerate(BD):
                        m = (fc >= lo) & (fc < hi); res[tag].setdefault((c, bi), []).append(10 * np.log10(max(err(sym, smm, syy, T, m, c), 1e-30) / syy[c][m].sum()))
    sN = tot([(a, b) for k, (a, b, c) in TR.items()])
    for k, v in S.items():
        other = "B" if v["rt"] in FOLDS["A"] else "A"; sS = tot([(S[q]["sym"], S[q]["smm"]) for q in S_keys(other)])
        T = T_from([(*sS, WS), (*sN, WN)], SM)
        for c in (0, 1):
            for bi, (lo, hi) in enumerate(BD):
                m = (fc >= lo) & (fc < hi); res["S"].setdefault((c, bi), []).append(10 * np.log10(max(err(v["sym"], v["smm"], v["syy"], T, m, c), 1e-30) / v["syy"][c][m].sum()))
    return res, len(TR)
def row(res):
    return {t: [float(np.median(res[t][(c, bi)])) for bi in range(len(BD)) for c in (0, 1)] for t in res}
def pr(name, r, n, ref):
    s = f"{name:46s} n{n:3d} "
    for t in ("N0", "Nall", "S"):
        s += f" {t}:" + " ".join(f"{x:6.1f}" for x in r[t][:4]) + ("  d%+.2f" % np.max(np.abs(np.array(r[t][:4]) - np.array(ref[t][:4]))) if ref else "")
    print(s, flush=True)
if __name__ == "__main__":
    print("columns per test set: 20-40 L R, 40-80 L R (dB re |Y|^2, median); d = max |change| vs baseline on that test set; N0 = original 51 captures (n=%d), Nall n=%d" % (len(T0), len(Tall)))
    base_res, nb = evaluate(list(DES)[0]); ref = row(base_res); pr(list(DES)[0], ref, nb, None)
    # no table
    for des in list(DES)[1:]:
        r, n = evaluate(des); pr(des, row(r), n, ref)
    print("--- smoothing (+-bins) on the base design")
    for SM in (0, 1, 4, 8):
        r, n = evaluate(list(DES)[0], SM=SM); pr(f"smoothing +-{SM} bins ({2*SM+1} x 2 Hz)", row(r), n, ref)
    print("--- program weights N:S on the base design")
    for lab, wn, ws in (("N:S 1:2", 1, 2), ("N:S 2:1", 2, 1), ("N only", 1, 1e-9), ("S only", 1e-9, 1)):
        r, n = evaluate(list(DES)[0], WN=wn, WS=ws); pr(lab, row(r), n, ref)
    print("--- weights on the all-of-the-above design")
    for lab, wn, ws in (("N:S 1:1", 1, 1), ("N:S 1:2", 1, 2), ("N:S 2:1", 2, 1)):
        r, n = evaluate("all of the above (w60f)", WN=wn, WS=ws); pr(lab, row(r), n, ref)
