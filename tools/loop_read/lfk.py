"""An LF-dependent HF pole term (model_s_results §25.20): the law's per-LF factor k (T = 10^(-3 P / (fs RT) k (1/hf - 1)))
refitted on the per-loop poles (hfloop_read.json), leaving each RT step out, with the §25.18 ramp (that fold's c, p0 0.62)
and ceiling 0.8475 on top; held-out pole error (precision-weighted, well-determined poles) against the §25.18 baseline
(law + ramp + ceiling + HF table) and the same structure with the law's k. Writes analysis/lowband/chain3_lfk.json
(k per LF per fold).

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/loop_read/lfk.py"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.optimize import least_squares
import dense_mid as dm
q1, ls = dm.q1, dm.ls
law = dm.q2.final_injection(compact=ls.lb.lf_.COMPACT)[1]
law["lf"] = {float(k): x for k, x in law["lf"].items()}
R = json.load(open(os.path.join(ROOT, "analysis/lowband/chain3_hfloop_read.json")))["rows"]
CAP, P0 = 0.8475, 0.62
folds = {float(f[5:-5]): json.load(open(os.path.join(ROOT, "analysis/loop_cv/cap0.8475_folds", f))) for f in os.listdir(os.path.join(ROOT, "analysis/loop_cv/cap0.8475_folds"))}
P = np.array([q1.P[r["loop"]] for r in R]); rt = np.array([r["rt"] for r in R]); hf = np.array([r["hf"] for r in R]); lf = np.array([r["lf"] for r in R])
p = np.array([r["p"] for r in R]); pb = np.array([r["p_base"] for r in R]); sig = np.abs(np.array([r["pA"] - r["pB"] for r in R])) / 2 + 1e-5
good = (sig < 3e-4) & (np.array([r["p_law"] for r in R]) < 0.86)
w = 1 / sig ** 2
LFS = sorted(set(lf))
def pred(kv, sel, c):
    k = np.array([kv[LFS.index(l)] for l in lf[sel]])
    T = 10 ** (-3 * P[sel] / q1.FS / rt[sel] * k * (1 / hf[sel] - 1))
    pl = q1.onepole_for(T, law["f_ref"])
    return np.minimum(pl + c * np.maximum(0, pl - P0), CAP)
k0 = np.array([law["lf"][l]["k"] for l in LFS])
out = {"A: §25.18 baseline (law + ramp + ceiling + HF table)": pb - p}
eB, eC, KF = np.zeros_like(p), np.zeros_like(p), {}
for r_ in sorted(set(rt)):
    c = folds[r_]["pramp"][0]
    tr, te = good & (rt != r_), rt == r_
    eB[te] = pred(k0, te, c) - p[te]
    f = least_squares(lambda kv: (pred(kv, tr, c) - p[tr]) * np.sqrt(w[tr] / w[tr].mean()), k0, bounds=(0.5, 1.5))
    eC[te] = pred(f.x, te, c) - p[te]; KF[str(r_)] = f.x.tolist()
out["B: law k + ramp + ceiling, no table"] = eB
out["C: k per LF refitted (4 numbers) + ramp + ceiling, no table"] = eC
rms = lambda e, s: 1e3 * np.sqrt(np.average(e[s] ** 2, weights=w[s]))
print(f"held-out pole error (rms x 1e-3, precision-weighted, {good.sum()} well-determined poles); law k {np.round(k0, 4).tolist()}")
print("candidate                                                   | all  | " + " | ".join(f"HF {h:g}" for h in sorted(set(hf))) + " | " + " | ".join(f"LF {l:g}" for l in LFS))
for nm, e in out.items():
    print(f"{nm:60s}| {rms(e, good):4.2f} | " + " | ".join(f"{rms(e, good & (hf == h)):5.2f}" for h in sorted(set(hf))) + " | " + " | ".join(f"{rms(e, good & (lf == l)):5.2f}" for l in LFS))
print("k per LF refitted, per held-out RT fold: " + "; ".join(f"RT {r_}: {np.round(v, 4).tolist()}" for r_, v in KF.items()))
json.dump(dict(lfs=LFS, k_law=k0.tolist(), folds=KF), open(os.path.join(ROOT, "analysis/lowband/chain3_lfk.json"), "w"), indent=1)
