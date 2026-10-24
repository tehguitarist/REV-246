"""NEXT_STEPS §5 item 5: the transferable loop terms of docs/model_s_results.md §22.3-22.8 together, every constant
chosen leaving one RT step out (8 folds), each setting scored with the constants of the fold that held out its RT step.
Results: docs/model_s_results.md §22.9.

    EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --stage constants      # -> analysis/loop_cv/fold_<rt>.json, constants.log
    for rt in 0.3 0.6 0.8 1.0 1.4 1.8 2.2 2.5; do EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --stage run --rt $rt; done
    EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --stage report         # -> report.log

--hfcells (any stage, implies --cells): also the HF pole shift per (HF, RT) cell (model_s_results §25.7: lowband_shelf.py's
hfcell grid, chosen per cell on its other LF steps) in cells.json; outputs in analysis/loop_cv/hfcells/, reported against cells/.
--cells (any stage, implies --lf04dc): also the DC floor's constants per (LF, RT) cell under the floor (model_s_results
§25.6: lowband_shelf.py's floor grid with gamma on the Nyquist end, chosen per cell on its other HF steps, so held out by
setting) -> cells.json; outputs in analysis/loop_cv/cells/, reported against the lf04dc folds.
--lf04dc (any stage): also LF 0.4's shelf DC ratio per fold (model_s_results §25: the median of the per-setting best
of lowband_shelf.py's dc stage, lfrt kept, chain3, over the fold's training RT steps); outputs in analysis/loop_cv/lf04dc/,
reported against the §22.9 folds above.
--retake (run / report; implies --hfcells): the §25.7 baseline re-taken on this kernel (model_s_results §25.11): the
hfcells folds' constants and cells.json unchanged, the compact numbers refitted on this kernel (q2_injection diffuser
stage, diffuser_tied_fix-0.707107_ap230new_<kernel>); outputs in analysis/loop_cv/hfcells_<kernel>[_<xadd>]/, reported
against hfcells/. --xadd TAG: also dense_mid.py's correction x_<fold>_<TAG>.npy, the fold whose training RT steps
exclude the setting's (so held out by RT step, 2 folds); --refs: the repeat-capture references (default nulls_free on
this kernel, else chain2's).

The terms (G_capped / q2_injection final options) and how each fold fits them, from the other RT steps only:
  gamma   a constant extra broadband loss where the DC floor acts (LF 0.1): the floor grid (loop_fir floor stage), DC_MIN
          0.18 / beta 0.019, the gamma with the lowest mean cost re gamma = 0 over the LF 0.1 settings;
  lfrt    LF 2.0: rate x (1 + a + b RT) and LF 0.4: rate x (1 + a max(0, 1.8 - RT)) (gated, §22.7), least squares on the
          broadband corrections of the loop_fir fits (median over HF 0.6 / 1.0, 63 Hz-2 kHz, the four loops);
  hfdp    HF 0.2: the pole + dp where 0.5 <= p < 0.62 and LF != 0.1 (§22.5, §22.7), the median of the loop_fir corrections
          read as a pole change (loops with that law pole);
  pramp   p = min(p_law + c max(0, p_law - 0.62), 0.84) (§22.8): per cell above the knee its best constant offset from the
          hfcap grid (cap 0.84), c = the median of dp / (median law pole - 0.62) over the cells.
"""
import argparse, glob, json, os, subprocess, sys

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import loop_fir as lf_  # noqa: E402

q2 = lf_.q2
q1 = q2.q1
FS = q1.FS
RTS = (0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5)
KNEE, LF04_END, BAND = 0.62, 1.8, (0.5, 0.62)
OUT = os.path.join(ROOT, "analysis", "loop_cv")
LFIR = lf_.OUT                                    # analysis/loop_fir/<kernel>
BASE_OUT = OUT
SHELF_DC = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_dcK")   # lowband_shelf.py --stage dc --keep-lfrt
SHELF_FLOOR = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_floorN")   # lowband_shelf.py --stage floor --gamma-nyq
SHELF_HF = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_hfcell")      # lowband_shelf.py --stage hfcell
CELLS = False
HFCELLS = False
RETAKE, XADD, REFS, PCURVE = False, None, None, False
DENSE = os.path.join(ROOT, "analysis", "dense_mid")


def hfcells_table():
    """Per setting with HF < 1: the HF pole shift of its (HF, RT) cell chosen on the cell's other LF steps (lowest mean
    100 Hz-10 kHz cost change, lowband_shelf.report_hfcell's rule)."""
    D = [json.load(open(p)) for p in glob.glob(os.path.join(SHELF_HF, "*.json"))]
    rel = {d["key"]: np.array([10 * np.log10(g["cost"] / d["base"]) for g in d["grid"]]) for d in D}
    out = {}
    for d in D:
        oth = [rel[e["key"]] for e in D if e["hf"] == d["hf"] and e["rt"] == d["rt"] and e["key"] != d["key"]]
        out[d["key"]] = float(d["grid"][int(np.argmin(np.mean(oth, 0)))]["dp"])
    return out


def cells_table():
    """Per setting under the floor (LF 0.1, LF 0.4 / RT 0.3): the floor constants of its cell chosen on the cell's other
    HF steps (lowband_shelf.cell_choice)."""
    import lowband_shelf as ls
    D = [json.load(open(p)) for p in glob.glob(os.path.join(SHELF_FLOOR, "*.json"))]
    out = {}
    for d in D:
        (dc, dc_min, beta, gamma), _ = ls.cell_choice(D, d["key"])
        out[d["key"]] = dict(dc=float(dc), dc_min=float(dc_min), beta=float(beta), gamma=float(gamma), gamma_dc=False)
    return out


def lf04_dc(tr):
    """LF 0.4's DC ratio from the training RT steps tr: the median of the per-setting best (25-1000 Hz cost)."""
    best = []
    for path in glob.glob(os.path.join(SHELF_DC, "*.json")):
        d = json.load(open(path))
        if d["rt"] in tr:
            best.append(d["grid"][int(np.argmin([g["cost"] for g in d["grid"]]))]["dc"])
    return float(np.median(best))


def law():
    lw = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    lw["lf"] = {float(k): x for k, x in lw["lf"].items()}
    return lw


def tables():
    """Per setting from the loop_fir fits: the broadband correction (fraction of the law's broadband rate) and, per loop,
    the correction read as a pole change (dp) with the law's pole."""
    lw = law()
    F = np.array([63, 125, 250, 500, 1000, 2000]); w = 2 * np.pi * F / FS
    f2 = np.linspace(100, 9000, 200); w2 = 2 * np.pi * f2 / FS; z1 = np.exp(-1j * w2)
    op = lambda p: (1 - p) / (1 - p * z1)
    bb, dp = {}, []
    for path in glob.glob(os.path.join(LFIR, "*.json")):
        s = json.load(open(path))
        lab = dict(lf=s["lf"], rt=s["rt"], hf=s["hf"])
        if s["hf"] in (0.6, 1.0):
            v = [np.median(20 * np.log10(np.abs(np.exp(-1j * np.outer(w, s["lags"])) @ np.array(c)))) * FS / Pi for c, Pi in zip(s["fir"], q1.P)]
            bb[(s["lf"], s["hf"], s["rt"])] = np.median(v) / (60 / (lw["lf"][s["lf"]]["hi"] * s["rt"]))
        if s["hf"] < 1.0 and s["lf"] != 0.1:
            E = np.exp(-1j * np.outer(w2, s["lags"]))
            for i, c in enumerate(s["fir"]):
                p = lf_.law_coefs(lab, lw)[i][3]
                C = E @ np.array(c)

                def r(u):
                    d = np.log((1 + u[1]) * op(p + u[0]) / op(p)) - np.log(C)
                    return np.r_[d.real, d.imag]
                u = least_squares(r, [0.0, 0.0], bounds=([-0.3, -0.2], [min(0.3, 0.999 - p), 0.2])).x
                dp.append(dict(**lab, p=p, dp=float(u[0])))
    return bb, dp


def constants_stage():
    os.makedirs(OUT, exist_ok=True)
    bb, dp = tables()
    floor = json.load(open(lf_.OUT + "_floor.json"))
    hfcap = json.load(open(lf_.OUT + "_hfcap.json"))
    lw = law()
    log = []
    for hold in RTS + ("all",):                    # "all": every RT step, the constants of the consolidated model
        tr = [r for r in RTS if r != hold]
        # gamma
        gs = sorted({g["gamma"] for g in floor[0]["grid"]})
        cost = lambda s, ga: next(g["cost"] for g in s["grid"] if (g["dc_min"], g["beta"], g["gamma"]) == (0.18, 0.019, ga))
        tot = {ga: np.mean([cost(s, ga) / cost(s, 0.0) for s in floor if s["lf"] == 0.1 and s["rt"] in tr]) for ga in gs}
        gamma = min(tot, key=tot.get)
        # LF forms
        y = lambda lf, rts: np.array([np.median([bb[(lf, hf, rt)] for hf in (0.6, 1.0)]) for rt in rts])
        r_ = np.array(tr)
        a20, b20 = np.linalg.lstsq(np.stack([np.ones_like(r_), r_], 1), y(2.0, tr), rcond=None)[0]
        x4 = np.maximum(0, LF04_END - r_); a04 = float(np.dot(x4, y(0.4, tr)) / np.dot(x4, x4))
        # HF 0.2 shift in the band
        band = [d["dp"] for d in dp if d["hf"] == 0.2 and d["rt"] in tr and BAND[0] <= d["p"] < BAND[1]]
        dpb = float(np.median(band)) if band else 0.0
        # ramp: cells above the knee (by their median law pole), best constant offset at cap 0.84
        cells = {}
        for s in hfcap:
            if s["rt"] not in tr:
                continue
            pm = np.median([float(q1.onepole_for(np.array(10 ** (-3 * Pi / FS / s["rt"] * lw["lf"][s["lf"]]["k"] * (1 / s["hf"] - 1))), lw["f_ref"]))
                            for Pi in q1.P])                       # the law's uncapped pole
            if KNEE + 0.03 < pm < 0.84:
                base = next(g["cost"] for g in s["grid"] if g["pmax"] is None)
                cells.setdefault((s["hf"], s["rt"]), dict(p=pm, rel={}))
                for g in s["grid"]:
                    if g["pmax"] == 0.84:
                        cells[(s["hf"], s["rt"])]["rel"].setdefault(g["dp"], []).append(g["cost"] / base)
        cs = []
        for (hf, rt), c in cells.items():
            best = min(c["rel"], key=lambda d: np.mean(c["rel"][d]))
            cs.append(best / (c["p"] - KNEE))
        cramp = float(np.median(cs)) if cs else 0.0
        fold = dict(hold=hold, train_rt=tr, gamma=gamma, **({"lfdc": {"0.4": lf04_dc(tr)}} if OUT != BASE_OUT else {}),
                    lfrt=dict(form="a + b RT", coefs={"0.4": [LF04_END * a04, -a04], "2.0": [float(a20), float(b20)]}, clip=["0.4"]),
                    hfdp={"0.2": dpb, "p_above": BAND[0], "p_below": BAND[1], "exclude_lf": [0.1]},
                    pramp=[cramp, KNEE], ramp_cells={f"{hf}_{rt}": c["p"] for (hf, rt), c in cells.items()})
        json.dump(fold, open(os.path.join(OUT, f"fold_{hold}.json"), "w"), indent=1)
        line = (f"{'all RT steps' if hold == 'all' else f'hold RT {hold}'}: gamma {gamma:+.2f}; LF 0.4 a {a04:+.5f}; LF 2.0 ({a20:+.5f}, {b20:+.5f}); HF 0.2 dp {dpb:+.4f} "
                f"(n {len(band)}); ramp c {cramp:.3f} from {sorted(fold['ramp_cells'])}" + (f"; LF 0.4 dc {fold['lfdc']['0.4']:.5f}" if "lfdc" in fold else ""))
        print(line); log.append(line)
    if CELLS:
        t = cells_table()
        if HFCELLS:
            for k, dp in hfcells_table().items():
                t.setdefault(k, {})["hf_dp"] = dp
        json.dump(t, open(os.path.join(OUT, "cells.json"), "w"), indent=1)
        for k in sorted(t):
            line = f"cell {k}: " + ", ".join(f"{n} {x:g}" for n, x in t[k].items() if n != "gamma_dc")
            print(line); log.append(line)
    open(os.path.join(OUT, "constants.log"), "w").write("\n".join(log) + "\n")


def run_stage(rt, workers):
    fold = json.load(open(os.path.join(OUT, f"fold_{rt}.json")))
    for k in ("lfrt", "hfdp"):
        json.dump(fold[k], open(os.path.join(OUT, f"fold_{rt}_{k}.json"), "w"), indent=1)
    info = json.load(open(os.path.join(q2.Q1, q1.DATA)))
    keys = ",".join(k for k, g in info.items() if g["rt"] == rt)
    compact = f"diffuser_tied_fix-0.707107_ap230new_{q1.ms.KTAG}" if RETAKE else lf_.COMPACT
    cmd = [sys.executable, os.path.join(HERE, "q2_injection.py"), "--stage", "final", "--compact", compact,
           "--ap230", os.path.join(q2.OUT, lf_.AP230), "--gamma", str(fold["gamma"]),
           "--lfrt", os.path.join(OUT, f"fold_{rt}_lfrt.json"), "--hfdp", os.path.join(OUT, f"fold_{rt}_hfdp.json"),
           "--pramp", f"{fold['pramp'][0]},{fold['pramp'][1]}", "--only", keys, "--workers", str(workers)]
    if "lfdc" in fold:
        cmd += ["--lfdc", ",".join(f"{k}:{x}" for k, x in fold["lfdc"].items())]
    if CELLS:
        cmd += ["--floorcells", os.path.join(OUT, "cells.json")]
    if PCURVE:                                      # §25.17: the pole curve fitted without this RT step replaces pramp, hfdp, HF cells
        pc = json.load(open(os.path.join(ROOT, "analysis", "lowband", "chain3_polecurve.json")))
        th = pc["folds"][str(rt)] if str(rt) in pc["folds"] else pc["folds"][repr(rt)]
        json.dump(dict(knots=pc["knots"], values=th[:len(pc["knots"])], cap=th[len(pc["knots"])]), open(os.path.join(OUT, f"fold_{rt}_pcurve.json"), "w"), indent=1)
        i = cmd.index("--hfdp"); del cmd[i:i + 2]
        i = cmd.index("--pramp"); del cmd[i:i + 2]
        cmd += ["--pcurve", os.path.join(OUT, f"fold_{rt}_pcurve.json")]
    if XADD:                                        # dense_mid's 2 folds: A trains on RT 0.3 / 0.8 / 1.4 / 2.2
        fd = "B" if rt in (0.3, 0.8, 1.4, 2.2) else "A"
        cmd += ["--xadd", os.path.join(DENSE, q1.ms.KTAG, f"x_{fd}_{XADD}.npy")]
    with open(os.path.join(OUT, f"fold_{rt}.log"), "w") as f:
        subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
    src = sorted(glob.glob(os.path.join(q2.OUT, q1.ms.tagged("final_s") + "_refit_lfrt_gamma*" + ("" if PCURVE else "_hfdp_pramp*") + ("_lfdc" if "lfdc" in fold else "") + ("_cells" if CELLS else "") + ("_xadd" if XADD else "") + ("_pcurve" if PCURVE else "") + "_smoke.json")),
                 key=os.path.getmtime)[-1]
    os.replace(src, os.path.join(OUT, f"fold_{rt}_rows.json"))


def report_stage():
    if RETAKE:                                      # §25.7 (chain2) -> re-taken on this kernel; with XADD: the re-take -> + the correction
        ref = os.path.join(BASE_OUT, "hfcells" + (f"_{q1.ms.KTAG}" if XADD or PCURVE else ""))
        a = {(r["key"], r["amp"]): r for p in glob.glob(os.path.join(ref, "fold_*_rows.json")) for r in json.load(open(p))["rows"]}
    elif OUT == BASE_OUT:                           # §21 -> the CV loop terms
        a = {(r["key"], r["amp"]): r for r in json.load(open(os.path.join(q2.OUT, q1.ms.tagged("final_s") + "_refit.json")))["rows"]}
    else:                                           # the §22.9 folds -> the same plus LF 0.4's dc (--cells: the lf04dc folds -> plus the floor cells)
        ref = os.path.join(BASE_OUT, "cells" if HFCELLS else "lf04dc") if CELLS else BASE_OUT
        a = {(r["key"], r["amp"]): r for p in glob.glob(os.path.join(ref, "fold_*_rows.json")) for r in json.load(open(p))["rows"]}
    c = {(r["key"], r["amp"]): r for p in glob.glob(os.path.join(OUT, "fold_*_rows.json")) for r in json.load(open(p))["rows"]}
    rp = REFS or os.path.join(q2.OUT, q1.ms.tagged("nulls_free") + ".json")
    if not os.path.exists(rp):
        rp = os.path.join(q2.OUT, "nulls_free_chain2.json")
    rep = {(r["key"], r["amp"]): r for r in json.load(open(rp))}
    ws = ["0-50", "50-150", "150-300", "300-600", "600-1000", "1000-1500"]

    def med(ks, src, w, ch):
        x = [src[k]["null_comp"][w][ch] for k in ks if src[k]["null_comp"].get(w)]
        return np.median(x) if x else np.nan

    def ex(ks, src, w):
        e = [src[k]["null_comp"][w][1] - rep[k]["repeat_comp"][w][1] for k in ks if src[k]["null_comp"].get(w) and rep.get(k) and rep[k]["repeat_comp"].get(w)]
        return np.median(e) if e else np.nan
    ks = list(c)
    out = [f"## All loop terms, leave-one-RT-out ({len({k[0] for k in ks})} settings, {len(ks)} captures; every setting scored with constants "
           f"fitted without its RT step); R null [excess over repeat]: " + ((f"re-taken on {q1.ms.KTAG} -> pole curve (no ramp, hfdp, HF cells)" if PCURVE else f"re-taken on {q1.ms.KTAG} -> + {XADD}" if XADD else f"§25.7 -> re-taken on {q1.ms.KTAG}") if RETAKE else "§21 -> CV" if OUT == BASE_OUT else "§25.6 -> + HF cells" if HFCELLS else "§25.4 -> + floor cells" if CELLS else "§22.9 -> + LF 0.4 dc") + "; kernel " + str(q1.ms.KTAG) + "; refs " + os.path.basename(rp)]
    out.append("| Settings | " + " | ".join(ws[1:5]) + " |\n|---|" + "---|" * 4)
    groups = [("All 128", lambda r: True), ("Flat", lambda r: r["lf"] == 1.0 and r["hf"] == 1.0), ("LF 0.1", lambda r: r["lf"] == 0.1),
              ("LF 0.4", lambda r: r["lf"] == 0.4), ("LF 2.0", lambda r: r["lf"] == 2.0), ("HF 0.2", lambda r: r["hf"] == 0.2),
              ("HF 0.4", lambda r: r["hf"] == 0.4), ("RT <= 1.0", lambda r: r["rt"] <= 1.0), ("RT >= 1.4", lambda r: r["rt"] >= 1.4)] + \
             [(f"RT {rt}", lambda r, rt=rt: r["rt"] == rt) for rt in RTS]
    for g, sel in groups:
        kk = [k for k in ks if sel(c[k])]
        out.append(f"| {g} | " + " | ".join(f"{med(kk, a, w, 1):.1f} [{ex(kk, a, w):+.1f}] -> {med(kk, c, w, 1):.1f} [{ex(kk, c, w):+.1f}]" for w in ws[1:5]) + " |")
    out.append("- All 128, L / R: " + "; ".join(f"{w}: {med(ks, a, w, 0):.1f} / {med(ks, a, w, 1):.1f} -> {med(ks, c, w, 0):.1f} / {med(ks, c, w, 1):.1f}" for w in ws))
    sets = sorted({k[0] for k in ks})
    d = {s: ex([k for k in ks if k[0] == s], c, "300-600") - ex([k for k in ks if k[0] == s], a, "300-600") for s in sets}
    d = {s: v for s, v in d.items() if np.isfinite(v)}                  # settings with no 300-600 ms window (RT 0.3) drop out
    out.append(f"- per setting, change of the 300-600 ms excess: median {np.median(list(d.values())):+.1f} dB; better (< -0.5) {sum(v < -0.5 for v in d.values())}, "
               f"worse (> +0.5) {sum(v > 0.5 for v in d.values())} of {len(d)}; worst: " +
               ", ".join(f"{s} {v:+.1f}" for s, v in sorted(d.items(), key=lambda kv: -kv[1])[:6]))
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["constants", "run", "report"], required=True)
    ap.add_argument("--rt", type=float, default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--lf04dc", action="store_true")
    ap.add_argument("--cells", action="store_true")
    ap.add_argument("--hfcells", action="store_true", help="also the HF pole shift per (HF, RT) cell (implies --cells)")
    ap.add_argument("--retake", action="store_true", help="the hfcells constants on this kernel's refitted compact numbers (implies --hfcells)")
    ap.add_argument("--xadd", default=None, help="with --retake: dense_mid x tag (after the fold letter) added to the injection")
    ap.add_argument("--pcurve", action="store_true", help="with --retake: the §25.17 pole curve (fitted without the fold's RT step) in place of pramp, hfdp and the HF cells")
    ap.add_argument("--refs", default=None, help="report: the repeat-capture reference file (default nulls_free on this kernel)")
    a = ap.parse_args()
    global OUT, CELLS, HFCELLS, RETAKE, XADD, REFS, PCURVE
    RETAKE, XADD, REFS, PCURVE = a.retake, a.xadd, a.refs, a.pcurve
    if a.lf04dc or a.cells or a.hfcells or a.retake:
        CELLS, HFCELLS = a.cells or a.hfcells or a.retake, a.hfcells or a.retake
        OUT = os.path.join(BASE_OUT, "hfcells" if a.hfcells or a.retake else "cells" if a.cells else "lf04dc")
        if a.retake:
            src, OUT = OUT, OUT + f"_{q1.ms.KTAG}" + (f"_{a.xadd}" if a.xadd else "") + ("_pcurve" if a.pcurve else "")
            os.makedirs(OUT, exist_ok=True)
            for f in [f"fold_{rt}.json" for rt in RTS] + ["cells.json"]:      # the constants, unchanged
                if not os.path.exists(os.path.join(OUT, f)):
                    open(os.path.join(OUT, f), "w").write(open(os.path.join(src, f)).read())
        os.makedirs(OUT, exist_ok=True)
    if a.stage == "constants":
        constants_stage()
    elif a.stage == "run":
        run_stage(a.rt, a.workers)
    else:
        report_stage()


if __name__ == "__main__":
    main()
