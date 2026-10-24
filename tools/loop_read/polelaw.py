"""Toward an HF pole law in the coefficient domain (model_s_results §25.17). Per S setting with HF < 1 and loop: the
law's uncapped pole p_law (onepole_for of the per-pass target at f_ref), the §25.6 baseline's pole p_base (ramp, fold
hfdp, cap) and the setting's own best common shift dp_own from lowband_shelf's hfcell grid (parabolic minimum, with the
cost curvature as a weight) -> p_true = p_base + dp_own. Writes analysis/lowband/chain3_polelaw.json.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/polelaw.py"""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import dense_mid as dm
q2, q1, ls = dm.q2, dm.q1, dm.ls
law = q2.final_injection(compact=ls.lb.lf_.COMPACT)[1]
law["lf"] = {float(k): x for k, x in law["lf"].items()}
rows = []
for p in sorted(glob.glob(os.path.join(ROOT, "analysis/lowband/chain3_shelf_hfcell/*.json"))):
    d = json.load(open(p))
    lf, rt, hf = d["lf"], d["rt"], d["hf"]
    dp = np.array([g["dp"] for g in d["grid"]]); c = 10 * np.log10(np.array([g["cost"] for g in d["grid"]]) / d["base"])
    j = int(np.argmin(c))
    if 0 < j < len(c) - 1:
        a, b, cc = np.polyfit(dp[j - 1:j + 2], c[j - 1:j + 2], 2)
        own, curv = (-b / (2 * a), a) if a > 0 else (dp[j], 0.0)
    else:
        own, curv = dp[j], 0.0
    t = ls.cv_terms(rt)
    L = law["lf"][lf]
    for i, Pi in enumerate(q1.P):
        kk = -3 * Pi / q1.FS
        T = float(10 ** (kk / rt * L["k"] * (1 / hf - 1)))
        p_law = float(q1.onepole_for(np.array(T), law["f_ref"]))
        pr = t["pramp"]
        p_base = min(p_law + pr[0] * max(0.0, p_law - pr[1]), q2.P_MAX)
        hd = t["hfdp"]
        if hd and str(hf) in hd and hd.get("p_above", 0.0) <= p_base < hd["p_below"] and lf not in hd.get("exclude_lf", []):
            p_base += hd[str(hf)]
        rows.append(dict(key=d["key"], lf=lf, rt=rt, hf=hf, loop=i, P=int(Pi), T=T, p_law=p_law, p_base=p_base, dp_own=float(own),
                         p_true=p_base + float(own), curv=float(curv), at_edge=bool(j in (0, len(c) - 1)), gain_db=float(c[j])))
json.dump(rows, open(os.path.join(ROOT, "analysis/lowband/chain3_polelaw.json"), "w"), indent=0)
print(f"{len(rows)} (setting, loop) rows from {len(rows) // 4} settings; at grid edge: {sum(r['at_edge'] for r in rows) // 4} settings")
for hf in sorted({r["hf"] for r in rows}):
    for rt in sorted({r["rt"] for r in rows}):
        u = [r for r in rows if r["hf"] == hf and r["rt"] == rt and r["loop"] == 0]
        if u:
            print(f"HF {hf} RT {rt}: p_law(loop 0) " + " ".join(f"{r['p_law']:.4f}" for r in sorted(u, key=lambda r: r['lf'])) +
                  " | p_true " + " ".join(f"{r['p_true']:.4f}" for r in sorted(u, key=lambda r: r['lf'])) +
                  " | dp_own " + " ".join(f"{r['dp_own']:+.4f}{'*' if r['at_edge'] else ''}" for r in sorted(u, key=lambda r: r['lf'])))
