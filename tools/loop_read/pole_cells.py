"""The HF pole in the pole-limited cells as one curve of the law's pole with a drop-back (model_s_results 32.4): a variant
fold set for loop_cv --retake --tag pcell, the §32.3 baseline's folds and cells with each pole-limited setting's four poles
(cells.json hf_poles) replacing its HF cell.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds \\
        SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_lowrt/cells.json \\
        SHELF_TAKEFIT=$PWD/analysis/lowband/chain3_shelf_takefit_cap0.8475_lowrt \\
        .venv/bin/python tools/loop_read/pole_cells.py    # seconds -> analysis/loop_cv/hfcells_chain3_cap0.8475_pcell/,
                                                          #            analysis/lowband/chain3_pole_cells.{log,json}
    for rt in 0.3 0.6 0.8; do EMT_HFRT_MIN=0.12 EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py \\
        --limit-gb 9 -- .venv/bin/python tools/loop_cv.py --retake --tag pcell --stage run --rt $rt --workers 8; done
    <same> .venv/bin/python tools/loop_read/pole_cells.py --all     # -> ..._pall/, chain3_pole_cells_all.*: every HF < 1 setting,
                                                                    #    all eight RT folds to rerun (--tag pall)

The unit's per-loop poles (baseline pole + takefit's pole shift, mean of the takes; settings whose takes disagree by more
than 2e-3 left out) against each loop's law pole (uncapped, the HF floor applied) lie on one curve for every LF step,
loop and RT step up to a law pole of about 0.84 (HF 0.2 / RT 0.8 and HF 0.4 / RT 0.3 have the same law poles and the same
unit poles to 3e-4); above about 0.856 the unit gives one value, about 0.8465, below the curve's top (the drop-back).
Pole-limited cells (HF 0.2 at RT <= 0.8, HF 0.4 at RT 0.3; §32.3's set): each loop's pole is the curve at its law pole,
a piecewise-linear curve on fixed knots, flat beyond the last, fitted on the other RT steps' points (held out by RT step),
or, at a law pole >= FOLD, the baseline's pole there (the 0.8475 ceiling + the (HF, RT) cell's shift, chosen held out by
LF step on the bench: the drop-back value). The takefit drop-back poles are not used for it: LF 0.4's (0.841-0.844, the
shelf knee's cell) pull a mean of the other LF steps 2e-3 under the unit at LF 2.0, which costs 12 dB there.
--all: the same curve (knots extended down to 0) for every HF < 1 setting at every RT step, in place of the ramp, the
ceiling, the HF cells and their slopes (the drop-back loops keep the baseline's pole), each RT step's curve fitted on the
other RT steps (model_s_results 32.5). --all --by-lf: each LF step's curve fitted on the other LF steps instead (32.6):
leaving an RT step out removes every point in a stretch of law poles only that RT step samples (0.64-0.75: HF 0.2 /
RT 1.0), while the claim to test is that the curve is the same for every LF step. --all --by-lf --lfrt-dp: also each
(LF, RT) cell's common offset from the curve (mean over the loops of unit - curve), taken from the cell's other HF steps
(held out by HF step), added to every curve loop (32.7).
"""
import glob, json, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import lowband_shelf as ls  # noqa: E402
import q1_generalise as q1  # noqa: E402
import q2_injection as q2  # noqa: E402

ALL = "--all" in sys.argv
BY_LF = "--by-lf" in sys.argv                             # with --all: hold out by LF step, not RT step (32.6)
LFRT_DP = "--lfrt-dp" in sys.argv                         # with --all --by-lf: + the (LF, RT) cell's offset (32.7)
GROUP = "lf" if BY_LF else "rt"
SRC = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_lowrt")
OUT = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_" + ("pallfd" if LFRT_DP else "pallf" if BY_LF else "pall" if ALL else "pcell"))
LOG = os.path.join(ROOT, "analysis", "lowband", "chain3_pole_cells" + ("_allfd" if LFRT_DP else "_allf" if BY_LF else "_all" if ALL else ""))
KNOTS = np.array(([0.0, 0.15, 0.25, 0.35, 0.45, 0.55] if ALL else []) + [0.62, 0.68, 0.74, 0.78, 0.81, 0.837])   # flat beyond the last
FOLD = 0.856                                              # between the highest law pole on the curve (0.8515) and the lowest dropped (0.8602)
LAM = 1e-6                                                # second-difference penalty on the knot values
POLE = (lambda lab: lab["hf"] < 1.0) if ALL else (lambda lab: (lab["hf"] == 0.2 and lab["rt"] <= 0.8) or (lab["hf"] == 0.4 and lab["rt"] == 0.3))


def poles():
    """Per HF < 1 setting with a takefit: each loop's law pole, baseline pole and unit pole."""
    v, law0, _ = q2.final_injection(compact=ls.lb.lf_.COMPACT)
    w = np.linspace(0, np.pi, 513)
    rows = []
    for p in sorted(glob.glob(os.path.join(ls.TAKEFIT_DIR, "*.json"))):
        t = json.load(open(p))
        if t["hf"] >= 1.0:                                      # no takefit (one clean take): law and baseline poles only
            continue
        lab = dict(lf=t["lf"], rt=t["rt"], hf=t["hf"])
        tt, lfdc, kw = ls.base_kw(t["key"], ls.HF_CELLS_JSON)
        law = json.loads(json.dumps(law0))
        law["lf"] = {float(k): x for k, x in law["lf"].items()}
        for lf_k, x in lfdc.items():
            law["lf"][float(lf_k)]["dc"] = float(x)
        args = dict(lfrt=tt["lfrt"], gamma=tt["gamma"], hfdp=tt["hfdp"], pramp=tt["pramp"], **kw)
        G, G1 = q2.final_Gs(lab, law, w, **args), q2.final_Gs(dict(lab, hf=1.0), law, w, **args)
        pb = np.array([(1 - h) / (1 + h) for h in [float(np.real(g[-1] / g1[-1])) for g, g1 in zip(G, G1)]])
        L = law["lf"][lab["lf"]]
        hf_e = max(lab["hf"], q2.HFRT_MIN / lab["rt"]) if q2.HFRT_MIN else lab["hf"]
        pl = [float(q1.onepole_for(np.array(10 ** (-3 * Pi / q1.FS / lab["rt"] * L["k"] * (1 / hf_e - 1))), law["f_ref"])) for Pi in q1.P]
        u = {k: pb + np.array(x[8:12]) for k, x in t["fit"].items() if len(x) > 8}
        rows.append(dict(key=t["key"], **lab, p_law=pl, pb=pb.tolist(), unit={k: x.tolist() for k, x in u.items()}))
    return rows


def hat(x, K=KNOTS):
    """Piecewise-linear basis on the knots, flat outside them."""
    x = np.clip(np.atleast_1d(np.asarray(x, float)), K[0], K[-1])
    B = np.zeros((len(x), len(K)))
    for j in range(len(K) - 1):
        m = (x >= K[j]) & (x <= K[j + 1])
        f = (x[m] - K[j]) / (K[j + 1] - K[j])
        B[m, j], B[m, j + 1] = 1 - f, f
    return B / B.sum(1, keepdims=True)


def fit_curve(pts):
    B, y = hat([p["pl"] for p in pts]), np.array([p["u"] for p in pts])
    D = np.diff(np.eye(len(KNOTS)), 2, axis=0)
    return np.linalg.solve(B.T @ B + LAM * D.T @ D, B.T @ y)


def main():
    rows = poles()
    pts = []
    for r in rows:
        u = r["unit"]
        if not u:
            continue
        if len(u) == 2 and np.abs(np.array(u["A"]) - np.array(u["B"])).max() > 2e-3:   # takes disagree (e.g. 2.0_0.8_0.2)
            continue
        uu = np.mean(list(u.values()), 0)
        pts += [dict(pl=r["p_law"][i], u=float(uu[i]), pb=r["pb"][i], rt=r["rt"], lf=r["lf"], hf=r["hf"], key=r["key"], loop=i) for i in range(4)]
    os.makedirs(OUT, exist_ok=True)
    for f in glob.glob(os.path.join(SRC, "fold_*.json")):
        rt = float(os.path.basename(f).split("_")[1].replace(".json", ""))
        if not f.endswith("_rows.json") or (rt > 0.8 and not ALL):   # RT >= 1.0 folds are unchanged: their rows too
            open(os.path.join(OUT, os.path.basename(f)), "w").write(open(f).read())
    cells = json.load(open(os.path.join(SRC, "cells.json")))
    info = q2.q1_tables()[0]
    curves, out, lines = {}, [], []
    for v in sorted({g[GROUP] for g in info.values()}) if ALL else (0.3, 0.6, 0.8):
        c = fit_curve([p for p in pts if p[GROUP] != v and KNOTS[0] - 0.02 <= p["pl"] < FOLD])
        curves[str(v)] = c.tolist()
        lines.append(f"{GROUP.upper()} {v} held out: curve at knots {KNOTS.tolist()} = " + " ".join(f"{x:.4f}" for x in c))
    by_key = {r["key"]: r for r in rows}
    for key, g in sorted(info.items()):
        lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
        if not POLE(lab) or key not in by_key:
            continue
        pl = np.array(by_key[key]["p_law"])
        c = np.array(curves[str(lab[GROUP])])
        dp = 0.0
        if LFRT_DP:                                         # the (LF, RT) cell's offset from its other HF steps (32.7)
            o = [np.mean([p["u"] - float((hat(p["pl"]) @ np.array(curves[str(p["lf"])]))[0]) for p in g4])
                 for g4 in ([p for p in pts if p["key"] == k2] for k2 in {p["key"] for p in pts if p["lf"] == lab["lf"]
                                                                       and p["rt"] == lab["rt"] and p["hf"] != lab["hf"]})
                 if len(g4) == 4 and all(p["pl"] < FOLD for p in g4)]
            dp = float(np.mean(o)) if o else 0.0
        hp = [float((hat(x) @ c)[0]) + dp if x < FOLD else pb for x, pb in zip(pl, by_key[key]["pb"])]   # dropped: the HF cell
        cells.setdefault(key, {})["hf_poles"] = [round(x, 6) for x in hp]
        unit = by_key[key]["unit"]
        uu = np.mean(list(unit.values()), 0) if unit else None
        out.append(dict(key=key, **lab, p_law=pl.tolist(), pb=by_key[key]["pb"], hf_poles=hp, unit=None if uu is None else uu.tolist()))
        lines.append(f"  {key}: law " + " ".join(f"{x:.4f}" for x in pl) + " | base " + " ".join(f"{x:.4f}" for x in by_key[key]["pb"])
                     + " | model " + " ".join(f"{x:.4f}" for x in hp) + ("" if uu is None else " | unit " + " ".join(f"{x:.4f}" for x in uu)))
    json.dump(cells, open(os.path.join(OUT, "cells.json"), "w"), indent=1)
    c_all = fit_curve([p for p in pts if KNOTS[0] - 0.02 <= p["pl"] < FOLD]).tolist()     # all data: the plugin's / M and L's map
    lines.append(f"all data: curve at knots {KNOTS.tolist()} = " + " ".join(f"{x:.4f}" for x in c_all))
    json.dump(dict(knots=KNOTS.tolist(), fold=FOLD, curves=curves, all=c_all, rows=out), open(LOG + ".json", "w"), indent=1)
    lines.insert(0, f"{len(out)} pole-limited settings with hf_poles -> {OUT}; drop-back above a law pole of {FOLD}")
    open(LOG + ".log", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
