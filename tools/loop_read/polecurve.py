"""Fit p_true = F(p_law) (polelaw.py's rows) as one curve in the coefficient domain, leave-one-RT-out, weighted by the
cost curvature; settings at the grid edge excluded. Candidates: (a) the baseline's ramp + cap (fold constants, no table);
(b) a smooth curve: piecewise-linear Delta(p_law) on fixed knots + cap; (c) the law with its own reference frequency
refitted (onepole_for at f') + cap. Reports the held-out pole error (rms, x 1e-3, weighted) per RT step and by p_law bin.
model_s_results §25.17.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/polecurve.py"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.optimize import least_squares
import dense_mid as dm
q1, ls = dm.q1, dm.ls
R = [r for r in json.load(open(os.path.join(ROOT, "analysis/lowband/chain3_polelaw.json"))) if not r["at_edge"] and r["curv"] > 0]
x = np.array([r["p_law"] for r in R]); y = np.array([r["p_true"] for r in R]); wt = np.sqrt(np.array([r["curv"] for r in R]))
wt = wt / wt.mean(); rt = np.array([r["rt"] for r in R]); T = np.array([r["T"] for r in R])
KN = np.array([0.0, 0.3, 0.45, 0.55, 0.62, 0.7, 0.78, 0.9])
CAP0 = 0.84
def curve(th, xx):
    d = np.interp(xx, KN, th[:len(KN)])
    return np.minimum(xx + d, th[len(KN)])
def fref(th, TT):
    return np.minimum(q1.onepole_for(TT ** th[1], th[0]), th[2])
def ramp(xx, rtv):
    out = []
    for xi, r_ in zip(xx, rtv):
        c, p0 = ls.cv_terms(r_)["pramp"]
        out.append(min(xi + c * max(0.0, xi - p0), CAP0))
    return np.array(out)
res = {"law + cap": np.minimum(x, CAP0) - y, "ramp + cap (no table)": ramp(x, rt) - y}
pc, pf = np.zeros_like(y), np.zeros_like(y)
TH = {}
for r_ in sorted(set(rt)):
    tr, te = rt != r_, rt == r_
    f1 = least_squares(lambda th: (curve(th, x[tr]) - y[tr]) * wt[tr], np.r_[np.zeros(len(KN)), CAP0])
    pc[te] = curve(f1.x, x[te]); TH[r_] = f1.x.tolist()
    f2 = least_squares(lambda th: (fref(th, T[tr]) - y[tr]) * wt[tr], [6400.0, 1.0, CAP0], bounds=([2000, 0.5, 0.7], [12000, 2.0, 0.95]))
    pf[te] = fref(f2.x, T[te])
res["smooth curve (11 numbers)"] = pc - y
res["law at refitted f_ref + exponent + cap (3)"] = pf - y
full = least_squares(lambda th: (curve(th, x) - y) * wt, np.r_[np.zeros(len(KN)), CAP0])
f2a = least_squares(lambda th: (fref(th, T) - y) * wt, [6400.0, 1.0, CAP0], bounds=([2000, 0.5, 0.7], [12000, 2.0, 0.95]))
rms = lambda e, s=slice(None): 1e3 * np.sqrt(np.sum((e[s] * wt[s]) ** 2) / np.sum(wt[s] ** 2))
print(f"{len(R)} (setting, loop) rows; held-out pole error, rms x 1e-3 (weighted), leave-one-RT-out")
print("candidate                                   | all  | " + " | ".join(f"RT {r_:g}" for r_ in sorted(set(rt))))
for nm, e in res.items():
    print(f"{nm:44s}| {rms(e):4.2f} | " + " | ".join(f"{rms(e, rt == r_):5.2f}" for r_ in sorted(set(rt))))
bins = [(0, 0.3), (0.3, 0.45), (0.45, 0.55), (0.55, 0.62), (0.62, 0.7), (0.7, 0.8), (0.8, 1.0)]
print("\nby p_law bin (mean signed error x 1e-3 / n): " + "; ".join(f"{a}-{b}" for a, b in bins))
for nm, e in res.items():
    print(f"{nm:44s} " + "  ".join(f"{1e3 * np.average(e[(x >= a) & (x < b)], weights=wt[(x >= a) & (x < b)]):+5.2f}/{np.sum((x >= a) & (x < b))}" if np.any((x >= a) & (x < b)) else "--" for a, b in bins))
print("\nall-RT fit, smooth curve: Delta at knots " + " ".join(f"{k:.2f}:{1e3 * v:+.2f}" for k, v in zip(KN, full.x[:len(KN)])) + f" (x1e-3); cap {full.x[len(KN)]:.4f}")
print(f"all-RT fit, refitted law: f_ref {f2a.x[0]:.0f} Hz, exponent {f2a.x[1]:.4f}, cap {f2a.x[2]:.4f}")
json.dump(dict(knots=KN.tolist(), folds=TH, all=full.x.tolist(), fref_all=f2a.x.tolist()), open(os.path.join(ROOT, "analysis/lowband/chain3_polecurve.json"), "w"), indent=1)
