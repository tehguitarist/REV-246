"""S's AP230 step 0 at b = 0 with the loops re-adapted (model_s_results.md section 39; the probe was section 38.2).

The 92 undamped S settings run AP230 at b = 0.00069 (the chain2 read, section 16.9); b = 0 alone helps 50-300 ms but costs
HF 0.2 / long RT late, because the loop tables (pole curve, shelf cells) were fitted around it. takefit (lowband_shelf.py) on
every undamped setting at both b gives, per setting, the 12 per-loop numbers (gain dB per pass, DC section dB, HF pole) the
loops want on top of the current model; their difference fit(b = 0) - fit(b = 0.00069) is the loops' share of the b change
(each setting's own loop error cancels in it). Here that change is made compact and held out, and written into a variant
fold set: lf_gain / lf_dc / hf_poles of each undamped setting moved by its predicted change (settings without lf_gain /
lf_dc get them), with AP230's step 0 at b = 0 (analysis/q2/ap230_mltable_b0.json, EMT_AP230).

    .venv/bin/python tools/s_step0.py --stage shift   -> analysis/s_step0/shift.{json,log}: the change per setting, and the
        compact forms' prediction error: per HF step (mean of the other settings of the step, leave one setting out),
        per (HF, RT) cell (mean of the cell's other LF steps, leave one LF step out), and per HF step linear in log RT
    .venv/bin/python tools/s_step0.py --stage cells --form <name>   -> analysis/loop_cv/hfcells_chain3_cap0.8475_apml0/
        (_apml's folds, its cells.json with the held-out change added; forms: hf_step, hf_rt_cell, hf_logrt, hf_rt_cell_nogain, zero)
Result (model_s_results 39): tried, not adopted; the variant fold set is not kept (its cells.json: analysis/s_step0/cells_<form>.json)
    then: for rt in 0.3 0.6 0.8 1.0 1.4 1.8 2.2 2.5; do EMT_AP230=ap230_mltable_b0.json <bench env> loop_cv.py --retake --tag apml0
        --stage run --rt $rt; done; foldcmp.py <_apml> <_apml0> --per-setting
"""
import argparse, collections, glob, json, os, shutil, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TF = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_takefit_cap0.8475_")
OUT = os.path.join(ROOT, "analysis", "s_step0")
BASE = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_apml")
VAR = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_apml0")
PARTS = (("gain", slice(0, 4)), ("dc", slice(4, 8)), ("pole", slice(8, 12)))


def deltas():
    """{key: dict(lf, rt, hf, d (12,), with the pole entries nan at HF 1.0)} from the two takefit runs (mean of the takes)."""
    out = {}
    for f in sorted(glob.glob(TF + "rb069/*.json")):
        a = json.load(open(f))
        g = TF + "rb0/" + os.path.basename(f)
        if not os.path.exists(g) or not a["fit"]:
            continue
        b = json.load(open(g))
        fit = lambda t: np.r_[np.mean([t["fit"][x] for x in t["fit"]], 0), np.full(12 - len(next(iter(t["fit"].values()))), np.nan)]
        out[a["key"]] = dict(lf=a["lf"], rt=a["rt"], hf=a["hf"], d=fit(b) - fit(a))
    return out


def forms(D):
    """Each compact form's held-out prediction per setting: {form: {key: (12,)}}."""
    keys = sorted(D)
    pred = collections.defaultdict(dict)
    for k in keys:
        d = D[k]
        same_hf = [j for j in keys if j != k and D[j]["hf"] == d["hf"]]
        pred["hf_step"][k] = np.nanmean([D[j]["d"] for j in same_hf], 0) if same_hf else np.zeros(12)
        cell = [j for j in keys if j != k and D[j]["hf"] == d["hf"] and D[j]["rt"] == d["rt"] and D[j]["lf"] != d["lf"]]
        pred["hf_rt_cell"][k] = np.nanmean([D[j]["d"] for j in cell], 0) if cell else pred["hf_step"][k]
        if len(same_hf) >= 3:                       # per HF step, each number linear in log RT, fitted without the setting
            x = np.log([D[j]["rt"] for j in same_hf])
            Y = np.array([D[j]["d"] for j in same_hf])
            p = np.full(12, np.nan)
            for i in range(12):
                ok = np.isfinite(Y[:, i])
                if ok.sum() >= 3:
                    c = np.polyfit(x[ok], Y[ok, i], 1)
                    p[i] = np.polyval(c, np.log(d["rt"]))
            pred["hf_logrt"][k] = p
        else:
            pred["hf_logrt"][k] = pred["hf_step"][k]
        pred["zero"][k] = np.zeros(12)
        pred["hf_rt_cell_nogain"][k] = np.r_[np.zeros(4), pred["hf_rt_cell"][k][4:]]      # the common gain left out (it overshoots late)
    return pred


def shift_stage():
    os.makedirs(OUT, exist_ok=True)
    D = deltas()
    P = forms(D)
    lines = [f"{len(D)} undamped settings with both takefit runs; change fit(b = 0) - fit(b = 0.00069), mean of the two takes"]
    for hf in sorted({d["hf"] for d in D.values()}):
        x = np.array([d["d"] for d in D.values() if d["hf"] == hf])
        lines.append(f"- HF {hf} ({len(x)}): gain mdB {np.round(np.nanmean(x[:, :4], 0) * 1e3, 2)} (sd {np.round(np.nanstd(x[:, :4], 0) * 1e3, 2)}); "
                     f"DC mdB {np.round(np.nanmean(x[:, 4:8], 0) * 1e3, 2)} (sd {np.round(np.nanstd(x[:, 4:8], 0) * 1e3, 2)}); "
                     f"pole 1e-4 {np.round(np.nanmean(x[:, 8:], 0) * 1e4, 2) if np.isfinite(x[:, 8:]).any() else '--'} (sd {np.round(np.nanstd(x[:, 8:], 0) * 1e4, 2) if np.isfinite(x[:, 8:]).any() else '--'})")
        for rt in sorted({d["rt"] for d in D.values() if d["hf"] == hf}):
            y = np.array([d["d"] for d in D.values() if d["hf"] == hf and d["rt"] == rt])
            lines.append(f"    RT {rt} ({len(y)}): gain {np.round(np.nanmean(y[:, :4]) * 1e3, 2)} mdB, DC {np.round(np.nanmean(y[:, 4:8]) * 1e3, 2)} mdB"
                         + (f", pole {np.round(np.nanmean(y[:, 8:]) * 1e4, 2)}e-4" if np.isfinite(y[:, 8:]).any() else ""))
    lines.append("\nheld-out prediction error, rms over settings and loops (zero = do nothing), by part:")
    for name, Pf in P.items():
        s = []
        for part, sl in PARTS:
            e = np.array([Pf[k][sl] - D[k]["d"][sl] for k in D])
            a = np.array([D[k]["d"][sl] for k in D])
            ok = np.isfinite(e) & np.isfinite(a)
            scale = 1e4 if part == "pole" else 1e3
            s.append(f"{part} {np.sqrt(np.mean(e[ok] ** 2)) * scale:.2f} (of {np.sqrt(np.mean(a[ok] ** 2)) * scale:.2f})")
        lines.append(f"- {name}: " + "; ".join(s) + "  [gain / DC mdB, pole 1e-4]")
    json.dump(dict(delta={k: dict(D[k], d=D[k]["d"].tolist()) for k in D},
                   pred={n: {k: v.tolist() for k, v in Pf.items()} for n, Pf in P.items()}),
              open(os.path.join(OUT, "shift.json"), "w"), indent=1)
    open(os.path.join(OUT, "shift.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def cells_stage(form):
    S = json.load(open(os.path.join(OUT, "shift.json")))
    pred = S["pred"][form]
    os.makedirs(VAR, exist_ok=True)
    for f in glob.glob(os.path.join(BASE, "fold_*.json")) + [os.path.join(BASE, "dcfloor.json")]:
        if not f.endswith("_rows.json"):
            shutil.copy2(f, VAR)
    for f in glob.glob(os.path.join(BASE, "fold_*_rows.json")):  # rows: replaced when each fold is rerun
        shutil.copy2(f, VAR)
    cells = json.load(open(os.path.join(BASE, "cells.json")))
    n = collections.Counter()
    for k, p in pred.items():
        p = np.nan_to_num(np.array(p))
        c = cells.setdefault(k, {})
        c["lf_gain"] = (np.array(c.get("lf_gain", np.zeros(4))) + p[:4]).tolist()
        c["lf_dc"] = (np.array(c.get("lf_dc", np.zeros(4))) + p[4:8]).tolist()
        n["loopcell"] += 1
        if "hf_poles" in c:
            c["hf_poles"] = (np.array(c["hf_poles"]) + p[8:12]).tolist()
            n["poles"] += 1
    json.dump(cells, open(os.path.join(VAR, "cells.json"), "w"), indent=1)
    json.dump(dict(form=form, ap230="ap230_mltable_b0.json", base=os.path.basename(BASE)), open(os.path.join(VAR, "step0.json"), "w"), indent=1)
    print(f"{VAR}: {form}, {n['loopcell']} settings' lf_gain / lf_dc moved, {n['poles']} with hf_poles")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["shift", "cells"], required=True)
    ap.add_argument("--form", default="hf_step")
    a = ap.parse_args()
    shift_stage() if a.stage == "shift" else cells_stage(a.form)


if __name__ == "__main__":
    main()
