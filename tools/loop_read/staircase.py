"""Is the HF pole a staircase of the requested loss (a firmware table indexed by the per-pass target)? Pooled S (polelaw.py)
and M / L (ml_poles.py) best poles minus S's smooth curve (polecurve.py, all RT); the remainder scanned for a periodic
component against two candidate indices (p_law; the per-pass target loss at f_ref in dB), periods 0.002-0.05 (p_law) /
0.05-2 dB; significance of the best period by permutation (shuffled remainders, 500x). M / L shortest RT and settings
without a drift zero left out. model_s_results §25.18.

    .venv/bin/python tools/loop_read/staircase.py"""
import json, os
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pc = json.load(open(os.path.join(ROOT, "analysis/lowband/chain3_polecurve.json")))
kn, th = np.array(pc["knots"]), np.array(pc["all"])
curve = lambda x: np.minimum(x + np.interp(x, kn, th[:len(kn)]), th[len(kn)])
X, Y, S = [], [], []
for r in json.load(open(os.path.join(ROOT, "analysis/lowband/chain3_polelaw.json"))):
    if not r["at_edge"] and r["curv"] > 0 and r["p_law"] < 0.84:
        X.append(r["p_law"]); Y.append(r["p_true"]); S.append("S")
SHORT = {"M": 0.6, "L": 2.2}
for prog in ("M", "L"):
    for r in json.load(open(os.path.join(ROOT, f"analysis/lowband/chain3_ml_poles_{prog}.json"))):
        if float(r["key"].split("_")[1]) == SHORT[prog] or not np.isfinite(r["dp_best"]):
            continue
        for pl, pb in zip(r["p_law"], r["p_best"]):
            if pl < 0.84:
                X.append(pl); Y.append(pb); S.append(prog)
X, Y, S = np.array(X), np.array(Y), np.array(S)
res = Y - curve(X)
c = np.cos(2 * np.pi * 6400 / (19.2e6 / 728))
loss = -20 * np.log10(np.abs((1 - X) / (1 - X * np.exp(-1j * 2 * np.pi * 6400 / (19.2e6 / 728)))))   # the pole's loss at f_ref
print(f"{len(X)} points (S {np.sum(S == 'S')}, M {np.sum(S == 'M')}, L {np.sum(S == 'L')}); remainder rms {1e3 * res.std():.2f} x 1e-3")
rng = np.random.default_rng(0)
for nm, u, periods in (("p_law", X, np.geomspace(0.002, 0.05, 400)), ("target loss dB", loss, np.geomspace(0.05, 2.0, 400))):
    def power(v):
        return np.array([abs(np.mean(v * np.exp(2j * np.pi * u / P))) for P in periods])
    pw = power(res - res.mean())
    j = int(np.argmax(pw))
    null = np.array([power(rng.permutation(res - res.mean())).max() for _ in range(500)])
    pval = np.mean(null >= pw[j])
    print(f"index {nm}: best period {periods[j]:.4g}, amplitude {2e3 * pw[j]:.2f} x 1e-3; permutation p = {pval:.3f} (null 95 % {2e3 * np.quantile(null, 0.95):.2f})")
    for prog in ("S", "M", "L"):
        sel = S == prog
        if sel.sum() > 5:
            ph = np.angle(np.mean((res[sel] - res[sel].mean()) * np.exp(2j * np.pi * u[sel] / periods[j])))
            print(f"   {prog}: amplitude at that period {2e3 * abs(np.mean((res[sel] - res[sel].mean()) * np.exp(2j * np.pi * u[sel] / periods[j]))):.2f} x 1e-3, phase {np.degrees(ph):+.0f} deg")
