"""The DC floor's compact form (model_s_results 34): one plateau per RT step plus one soft shoulder and one Nyquist-end
rule shared by every setting and loop, in place of the per-setting floor cells (dc, dc_min, beta, gamma) of §25.6.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds \\
        SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_pallf/cells.json \\
        SHELF_TAKEFIT=$PWD/analysis/lowband/chain3_shelf_takefit_cap0.8475_pallf \\
        .venv/bin/python tools/loop_read/floor_form.py --stage fit     # seconds -> analysis/lowband/chain3_floor_form.{log,json},
                                                                      #   <OUT>/dcfloor.json (EMT_DCFLOOR's table)
    <same> .venv/bin/python tools/loop_read/floor_form.py --stage folds   # -> <OUT>: _pallf's folds, cells.json without the
                                                                          #   floor cells (shelf cells still _pallf's)
    then takefit on <OUT> (EMT_DCFLOOR=<OUT>/dcfloor.json SHELF_CELLS=<OUT>/cells.json SHELF_TAG=_pfloor), and
    tools/loop_read/lf_cells.py --delta for both halves (the shelf cells refitted on the new floor).

Data: each loop's DC end and Nyquist end as the unit has them (the _pallf model's, times takefit's per-loop broadband
gain and DC section, mean of the two takes), against the law's request x (per-pass DC gain in dB, law DC ratio, no floor).
The form (q2.G_capped floor_form), on a cut shelf (LF 0.1, 0.4):
    DC end    D = max(x, F - c) + c exp(-max(x - F + c, 0) / tau)     F: the plateau per RT step; c, tau: the shoulder (dB)
    Nyquist   where x < F0 (F0: the median of the per-RT plateaus, no extra number): -beta (F0 - x) + gamma . B dB,
              B the law's broadband loss per pass (dB), gamma < 0
Hold-out: c, tau, beta, gamma (and F0) fitted without the scored RT step (leave-one-RT-out); F per RT step fitted on that
RT step's other HF steps (leave-one-HF-out, as §25.6 chose the cells). RT steps with no loop requesting under 0.2 per pass
(RT 2.2, 2.5) take the nearest such RT step's plateau. beta / gamma are fitted on LF 0.1's floored loops at RT <= 1.8
(LF 0.4 / RT 0.3, the only other floored cell, sits at RT <= 0.8, where the shelf cells keep their common gain).
Also: whether the per-RT plateaus follow a law in RT (constant, linear in RT, in 1 / RT), held out by RT step.

--stage shoulder (model_s_results 35): the shoulder alone where no floor cell acts (LF 0.4 at RT >= 0.6), the floor cells
kept everywhere else. c, tau and the per-RT plateaus are fitted on the floor-cell settings only (LF 0.1, LF 0.4 / RT 0.3:
held out by setting, nothing from LF 0.4 / RT >= 0.6) and written to <OUT2>/dcfloor.json for the 28 LF 0.4 / RT >= 0.6
settings (beta = gamma = 0: no loop there requests under the plateau, so the Nyquist end is the baseline's), and as "floor_form" into <OUT2>/cells.json (read without EMT_DCFLOOR); <OUT2> gets
_pallf's folds, cells.json (floor cells kept) and RT 0.3's rows (that fold is unchanged). Then takefit on the 28 keys with
SHELF_PREP_CELLS, lf_cells.py --delta (both halves) on <OUT2>, and loop_cv --retake --tag pshoulder for RT 0.6-2.5.
"""
import argparse, collections, glob, json, os, shutil, sys

import numpy as np
from scipy.optimize import least_squares

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import lowband_shelf as ls  # noqa: E402

q2, q1 = ls.q2, ls.q1
SRC = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_pallf")
OUT = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_pfloor")
OUT2 = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_pshoulder")      # --stage shoulder (35)
LOG = os.path.join(ROOT, "analysis", "lowband", "chain3_floor_form")
TF = os.environ.get("SHELF_TAKEFIT") or os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_takefit_cap0.8475_pallf")
RTS, HFS = [0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5], [0.2, 0.4, 0.6, 1.0]
XFL = 20 * np.log10(0.2)          # an RT step "has floored loops" if any loop requests less than 0.2 per pass
FLOOR_KEYS = ("dc", "dc_min", "beta", "gamma", "gamma_dc")
db = lambda x: 20 * np.log10(x)
rms = lambda a: float(np.sqrt(np.mean(np.square(a)))) if len(a) else float("nan")


def data():
    """Per (setting, loop) at LF 0.1 / 0.4: the law's request, the unit's DC and Nyquist ends, the current model's."""
    law0 = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    info = q2.q1_tables()[0]
    cells = json.load(open(ls.HF_CELLS_JSON))
    w = np.array([0.0, np.pi])
    rows = []
    for key in sorted(info):
        lab = dict(lf=info[key]["lf"], rt=info[key]["rt"], hf=info[key]["hf"])
        if lab["lf"] not in (0.1, 0.4):
            continue
        tf = json.load(open(os.path.join(TF, key + ".json")))
        fits = [np.array(tf["fit"][k]) for k in ("A", "B") if k in tf["fit"]]
        if not fits:
            continue
        t, lfdc, kw = ls.base_kw(key, ls.HF_CELLS_JSON)
        law, lawL = (json.loads(json.dumps(law0)) for _ in range(2))
        for x in (law, lawL):
            x["lf"] = {float(k): v for k, v in x["lf"].items()}
        for k_, x in lfdc.items():
            law["lf"][float(k_)]["dc"] = float(x)
        for k_, x in ls.cv_fold(lab["rt"])["lfdc"].items():          # the law as the form sees it (no cell dc)
            lawL["lf"][float(k_)]["dc"] = float(x)
        l1, args = dict(lab, hf=1.0), dict(lfrt=t["lfrt"], gamma=t["gamma"], hfdp=None, pramp=t["pramp"])
        G = q2.final_Gs(l1, law, w, **args, **kw)                                       # the current model, HF 1.0 form
        Gf = q2.final_Gs(l1, law, w, **args, **{k: v for k, v in kw.items() if k != "loopcell"})   # its floor alone
        Gl = q2.final_Gs(l1, lawL, w, **args, dc_min=0.0, beta=0.0)                    # the law (lfrt included)
        L = lawL["lf"][lab["lf"]]
        for i, Pi in enumerate(q1.P):
            kk = -3 * Pi / q1.FS
            g = 10 ** (kk / (L["hi"] * lab["rt"]))
            req = g * 10 ** (kk * (1 / (L["dc"] * lab["rt"]) - 1 / (L["hi"] * lab["rt"])))
            m_dc, m_ny = abs(G[i][0]), abs(G[i][1])
            rows.append(dict(key=key, **lab, loop=i, x=db(req), B=-db(g), lfr=db(abs(Gl[i][0]) / req), law_ny=db(abs(Gl[i][1])),
                             y=float(np.mean([db(m_dc) + u[i] + u[4 + i] for u in fits])),
                             yn=float(np.mean([db(m_ny) + u[i] for u in fits])) - db(abs(Gl[i][1])),
                             cur=db(m_dc), cur_f=db(abs(Gf[i][0])), cur_fn=db(abs(Gf[i][1])) - db(abs(Gl[i][1])),
                             cell="dc_min" in cells.get(key, {})))
    return rows


def dc_pred(d, F, c, tau):
    return d["lfr"] + max(d["x"], F - c) + c * np.exp(-max(d["x"] - F + c, 0.0) / tau)


def ny_pred(d, F0, beta, gamma):
    u = F0 - d["x"]
    return (-beta * u + gamma * d["B"]) if u > 0 else 0.0


def nearest(Fs):
    """Per-RT plateaus completed: RT steps without floored loops take the nearest RT step's value."""
    return {rt: Fs[rt] if rt in Fs else Fs[min(Fs, key=lambda r: abs(r - rt))] for rt in RTS}


def fit_shape(ds):
    """Per-RT plateaus (RT steps with floored loops) and the shoulder (c, tau), jointly on ds."""
    rts = sorted({d["rt"] for d in ds if d["x"] < XFL})
    n = len(rts)

    def res(v):
        F = nearest(dict(zip(rts, v[:n])))
        return [dc_pred(d, F[d["rt"]], v[n], v[n + 1]) - d["y"] for d in ds]
    f = least_squares(res, [-14.68] * n + [0.3, 2.0], bounds=([-17.0] * n + [0.0, 0.2], [-12.0] * n + [1.0, 10.0]))
    return nearest(dict(zip(rts, f.x[:n]))), float(f.x[n]), float(f.x[n + 1])


def fit_plateau(ds, c, tau):
    """One RT step's plateau with the shoulder fixed (None if the RT step has no floored loop)."""
    if not any(d["x"] < XFL for d in ds):
        return None
    return float(least_squares(lambda v: [dc_pred(d, v[0], c, tau) - d["y"] for d in ds], [-14.68], bounds=(-17.0, -12.0)).x[0])


def fit_ny(ds, F0):
    s = [d for d in ds if d["lf"] == 0.1 and d["rt"] <= 1.8 and d["x"] < F0]
    return tuple(float(x) for x in least_squares(lambda v: [ny_pred(d, F0, v[0], v[1]) - d["yn"] for d in s], [0.02, -0.03]).x)


def fit_stage():
    D = data()
    log = [f"# The DC floor's compact form (model_s_results 34); unit = the current model x takefit ({os.path.basename(TF)}), "
           f"cells {os.path.relpath(ls.HF_CELLS_JSON, ROOT)}", ""]
    Fa, ca, ta = fit_shape(D)
    F0a = float(np.median([Fa[rt] for rt in sorted({d['rt'] for d in D if d['x'] < XFL})]))
    ba, ga = fit_ny(D, F0a)
    log.append("All data: plateau F per RT step " + " ".join(f"{rt}:{10 ** (Fa[rt] / 20):.4f}" for rt in RTS)
               + f"; shoulder c {ca:.3f} dB, tau {ta:.2f} dB; F0 {10 ** (F0a / 20):.4f}; beta {ba:.4f}, gamma {ga:+.4f}")
    out, folds = [], {}
    for rt in RTS:                                    # c, tau, beta, gamma, F0 without this RT step
        tr = [d for d in D if d["rt"] != rt]
        Ft, c, tau = fit_shape(tr)
        F0 = float(np.median([Ft[r] for r in sorted({d['rt'] for d in tr if d['x'] < XFL})]))
        beta, gamma = fit_ny(tr, F0)
        folds[rt] = dict(c=c, tau=tau, F0=F0, beta=beta, gamma=gamma, F={})
        for hf in HFS:                                # the plateau from this RT step's other HF steps
            Fh = fit_plateau([d for d in D if d["rt"] == rt and d["hf"] != hf], c, tau)
            F = Fh if Fh is not None else Ft[rt]
            folds[rt]["F"][hf] = F
            out += [dict(d, p=dc_pred(d, F, c, tau), pn=ny_pred(d, F0, beta, gamma)) for d in D if d["rt"] == rt and d["hf"] == hf]
        log.append(f"fold RT {rt}: c {c:.3f} tau {tau:.2f} F0 {10 ** (F0 / 20):.4f} beta {beta:.4f} gamma {gamma:+.4f}; "
                   "F by HF " + " ".join(f"{hf}:{10 ** (folds[rt]['F'][hf] / 20):.4f}" for hf in HFS))
    # -- the table EMT_DCFLOOR reads: every LF 0.1 / 0.4 setting with its held-out values
    info = q2.q1_tables()[0]
    tab = {}
    for k, g in info.items():
        if g["lf"] in (0.1, 0.4):
            f = folds[g["rt"]]
            tab[k] = dict(F=round(10 ** (f["F"][g["hf"]] / 20), 6), c=round(f["c"], 5), tau=round(f["tau"], 5),
                          F0=round(10 ** (f["F0"] / 20), 6), beta=round(f["beta"], 6), gamma=round(f["gamma"], 6))
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(form="model_s_results 34", all=dict(F={str(rt): round(10 ** (Fa[rt] / 20), 6) for rt in RTS}, c=ca, tau=ta,
                                                        F0=10 ** (F0a / 20), beta=ba, gamma=ga), settings=tab),
              open(os.path.join(OUT, "dcfloor.json"), "w"), indent=1)
    # -- offline scores
    log += ["", "Held out (F by HF step, the rest by RT step), dB per pass rms over loops: DC end against the unit; "
            "Nyquist end's loop mean per setting (what the shelf cells' zero-mean gain cannot take at RT >= 1.0)",
            f"{'cell':14s} {'DC: form':>9s} {'cells':>7s} {'model':>7s} | {'Nyq mean: form':>15s} {'cells':>7s} | unit Nyq, form, cells (dB)"]
    groups = [("LF 0.1 + 0.4", lambda d: True), ("LF 0.1", lambda d: d["lf"] == 0.1), ("LF 0.4", lambda d: d["lf"] == 0.4),
              ("RT >= 1.0", lambda d: d["rt"] >= 1.0)] + \
             [(f"LF {lf} RT {rt}", (lambda lf, rt: lambda d: d["lf"] == lf and d["rt"] == rt)(lf, rt)) for lf in (0.1, 0.4) for rt in RTS]
    summ = {}
    for name, sel in groups:
        s = [d for d in out if sel(d)]
        bk = collections.defaultdict(list)
        for d in s:
            bk[d["key"]].append(d)
        nm = lambda f: rms([np.mean([f(d) for d in v]) for v in bk.values()])
        r = dict(dc_form=rms([d["p"] - d["y"] for d in s]), dc_cells=rms([d["cur_f"] - d["y"] for d in s]),
                 dc_model=rms([d["cur"] - d["y"] for d in s]), ny_form=nm(lambda d: d["pn"] - d["yn"]), ny_cells=nm(lambda d: d["cur_fn"] - d["yn"]))
        summ[name] = r
        log.append(f"{name:14s} {r['dc_form']:9.4f} {r['dc_cells']:7.4f} {r['dc_model']:7.4f} | {r['ny_form']:15.4f} {r['ny_cells']:7.4f} | "
                   f"{np.mean([d['yn'] for d in s]):+.3f} {np.mean([d['pn'] for d in s]):+.3f} {np.mean([d['cur_fn'] for d in s]):+.3f}")
    log.append("('cells' = the per-setting floor cells alone, 'model' = _pallf with its held-out shelf cells: the target the refitted shelf cells reach)")
    # -- per-RT plateaus as a law in RT, held out by RT step
    rts = [rt for rt in RTS if any(d["x"] < XFL and d["rt"] == rt for d in D)]
    Fr = {rt: fit_plateau([d for d in D if d["rt"] == rt], ca, ta) for rt in rts}
    log += ["", "Per-RT plateaus (all HF steps, shoulder fixed at the all-data fit): " + " ".join(f"{rt}:{10 ** (Fr[rt] / 20):.4f}" for rt in rts),
            "As a law in RT, held out by RT step (predicted - fitted plateau, dB):"]
    laws = {"constant": lambda r: np.ones((len(r), 1)), "a + b RT": lambda r: np.c_[np.ones(len(r)), r], "a + b / RT": lambda r: np.c_[np.ones(len(r)), 1 / np.array(r)]}
    lawres = {}
    for nm_, X in laws.items():
        e = []
        for rt in rts:
            tr = [r for r in rts if r != rt]
            co = np.linalg.lstsq(X(tr), np.array([Fr[r] for r in tr]), rcond=None)[0]
            e.append(float((X([rt]) @ co)[0] - Fr[rt]))
        lawres[nm_] = e
        log.append(f"  {nm_:10s} " + " ".join(f"{rt}:{x:+.3f}" for rt, x in zip(rts, e)) + f"  rms {rms(e):.3f}")
    log.append(f"  (the per-RT plateau held out by HF step instead: rms {rms([folds[rt]['F'][hf] - Fr[rt] for rt in rts for hf in HFS]):.4f} dB)")
    open(LOG + ".log", "w").write("\n".join(log) + "\n")
    json.dump(dict(all=dict(F={str(rt): Fa[rt] for rt in RTS}, c=ca, tau=ta, F0=F0a, beta=ba, gamma=ga),
                   folds={str(rt): dict(v, F={str(h): x for h, x in v["F"].items()}) for rt, v in folds.items()},
                   plateaus={str(rt): Fr[rt] for rt in rts}, laws=lawres, scores=summ), open(LOG + ".json", "w"), indent=1)
    print("\n".join(log))


def folds_stage():
    """<OUT>: _pallf's fold constants, cells.json without the floor cells (shelf cells and hf_poles unchanged)."""
    os.makedirs(OUT, exist_ok=True)
    for f in glob.glob(os.path.join(SRC, "fold_*.json")):
        if not f.endswith(("_rows.json", "_lfrt.json", "_hfdp.json")):
            shutil.copy(f, OUT)
    cells = json.load(open(os.path.join(SRC, "cells.json")))
    n = sum(1 for v in cells.values() if "dc_min" in v)
    cells = {k: {n_: x for n_, x in v.items() if n_ not in FLOOR_KEYS} for k, v in cells.items()}
    json.dump({k: v for k, v in cells.items() if v}, open(os.path.join(OUT, "cells.json"), "w"), indent=1)
    print(f"{n} floor cells removed -> {OUT}/cells.json")


def shoulder_stage():
    """--stage shoulder (model_s_results 35): c, tau and the plateaus fitted on the floor-cell settings, applied to LF 0.4 /
    RT >= 0.6 only; offline scores against the unit's DC ends there and the variant fold set <OUT2>."""
    D = data()
    fit = [d for d in D if d["cell"]]                                  # LF 0.1 and LF 0.4 / RT 0.3: the floor-cell settings
    tst = [d for d in D if d["lf"] == 0.4 and d["rt"] >= 0.6]
    F, c, tau = fit_shape(fit)
    pl = lambda rt: 10 ** (F[rt] / 20)
    log = [f"# The shoulder alone at LF 0.4 / RT >= 0.6 (model_s_results 35); unit = the current model x takefit ({os.path.basename(TF)})",
           f"Fitted on the {len({d['key'] for d in fit})} floor-cell settings ({len(fit)} loops), held out by setting: c {c:.4f} dB, "
           f"tau {tau:.4f} dB; plateau per RT step " + " ".join(f"{rt}:{pl(rt):.4f}" for rt in RTS),
           f"min request at LF 0.4 / RT >= 0.6: {min(d['x'] for d in tst):.2f} dB per pass (plateau {F[0.6]:.2f}): no loop under it", "",
           "DC end against the unit, dB per pass (LF 0.4, per RT step): rms and mean of form - unit; the baseline's hard floor without "
           "its shelf cells; _pallf with its shelf cells", f"{'RT':>4s} {'form rms':>9s} {'mean':>8s} {'floor':>7s} {'_pallf':>7s} | unit / form lift over the law, loops 1-4 (HF mean)"]
    sc = {}
    for rt in RTS[1:]:
        s_ = [d for d in tst if d["rt"] == rt]
        e = [dc_pred(d, F[rt], c, tau) - d["y"] for d in s_]
        lift = lambda f: " ".join(f"{np.mean([f(d) for d in s_ if d['loop'] == i]):+.3f}" for i in range(4))
        sc[str(rt)] = dict(form=rms(e), mean=float(np.mean(e)), floor=rms([d["cur_f"] - d["y"] for d in s_]), pallf=rms([d["cur"] - d["y"] for d in s_]))
        log.append(f"{rt:4.1f} {sc[str(rt)]['form']:9.4f} {sc[str(rt)]['mean']:+8.4f} {sc[str(rt)]['floor']:7.4f} {sc[str(rt)]['pallf']:7.4f} | "
                   + lift(lambda d: d["y"] - d["lfr"] - d["x"]) + " / " + lift(lambda d: dc_pred(d, F[rt], c, tau) - d["lfr"] - d["x"]))
    os.makedirs(OUT2, exist_ok=True)
    for f in glob.glob(os.path.join(SRC, "fold_*.json")):
        if not f.endswith(("_rows.json", "_lfrt.json", "_hfdp.json")) or os.path.basename(f) == "fold_0.3_rows.json":
            shutil.copy(f, OUT2)
    info = q2.q1_tables()[0]
    tab = {k: dict(F=round(pl(g["rt"]), 6), c=round(c, 5), tau=round(tau, 5), F0=round(pl(g["rt"]), 6), beta=0.0, gamma=0.0)
           for k, g in info.items() if g["lf"] == 0.4 and g["rt"] >= 0.6}
    cells = json.load(open(os.path.join(SRC, "cells.json")))        # the form also in cells.json ("floor_form", read by every
    for k, x in tab.items():                                         # reader where no floor cell acts): no EMT_DCFLOOR needed
        cells.setdefault(k, {})["floor_form"] = x
    json.dump(cells, open(os.path.join(OUT2, "cells.json"), "w"), indent=1)
    json.dump(dict(form="model_s_results 35: the shoulder alone at LF 0.4 / RT >= 0.6", all=dict(F={str(rt): pl(rt) for rt in RTS}, c=c, tau=tau),
                   settings=tab), open(os.path.join(OUT2, "dcfloor.json"), "w"), indent=1)
    log.append(f"\n{len(tab)} settings -> {os.path.relpath(OUT2, ROOT)}/dcfloor.json; _pallf's folds and cells.json (floor cells kept) copied")
    open(LOG.replace("floor_form", "floor_shoulder") + ".log", "w").write("\n".join(log) + "\n")
    json.dump(dict(c=c, tau=tau, F={str(rt): F[rt] for rt in RTS}, scores=sc), open(LOG.replace("floor_form", "floor_shoulder") + ".json", "w"), indent=1)
    print("\n".join(log))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=["fit", "folds", "shoulder"], required=True)
    a = ap.parse_args()
    {"fit": fit_stage, "folds": folds_stage, "shoulder": shoulder_stage}[a.stage]()
