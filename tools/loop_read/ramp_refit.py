"""The ramp refitted around another HF cap (model_s_results §25.18): per RT fold, c of p = min(p_law + c max(0, p_law -
p0), cap), p0 = 0.62, by weighted least squares on the other RT steps' own best poles (polelaw.py's rows, the cost
curvature as weight, grid-edge settings out); the fold files of loop_cv's lf04dc set with only pramp replaced ->
analysis/loop_cv/cap<cap>_folds/. Also the same fit at cap 0.84 as a check against §22.8's c.

    .venv/bin/python tools/loop_read/ramp_refit.py 0.8475"""
import json, os, sys
import numpy as np
from scipy.optimize import minimize_scalar
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
cap, P0 = float(sys.argv[1]), 0.62
R = [r for r in json.load(open(os.path.join(ROOT, "analysis/lowband/chain3_polelaw.json"))) if not r["at_edge"] and r["curv"] > 0]
x = np.array([r["p_law"] for r in R]); y = np.array([r["p_true"] for r in R]); w = np.array([r["curv"] for r in R]); rt = np.array([r["rt"] for r in R])
def fit(sel, cp):
    f = lambda c: np.sum(w[sel] * (np.minimum(x[sel] + c * np.maximum(0, x[sel] - P0), cp) - y[sel]) ** 2)
    return minimize_scalar(f, bounds=(0, 0.5), method="bounded").x
src, out = os.path.join(ROOT, "analysis/loop_cv/lf04dc"), os.path.join(ROOT, f"analysis/loop_cv/cap{sys.argv[1]}_folds")
os.makedirs(out, exist_ok=True)
for r_ in sorted(set(rt)):
    c_new, c_old = fit(rt != r_, cap), fit(rt != r_, 0.84)
    f = json.load(open(os.path.join(src, f"fold_{r_:g}.json" if os.path.exists(os.path.join(src, f"fold_{r_:g}.json")) else f"fold_{r_}.json")))
    old = f["pramp"]
    f["pramp"] = [float(c_new), P0]
    f["cap"] = cap
    json.dump(f, open(os.path.join(out, f"fold_{f['hold']}.json"), "w"), indent=1)
    print(f"fold RT {r_:g}: c at cap {cap} = {c_new:.4f} (same fit at cap 0.84: {c_old:.4f}; §22.8's fold value {old[0]:.4f})")
