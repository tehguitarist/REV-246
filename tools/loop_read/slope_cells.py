"""The HF cells' slope across the loops as a baseline candidate (model_s_results 32): a variant fold set for
loop_cv --retake --tag slope, the §29 baseline's folds and cells with each HF < 1 setting's hf_slope added.

    .venv/bin/python tools/loop_read/slope_cells.py          # -> analysis/loop_cv/hfcells_chain3_cap0.8475_slope/
    for rt in 0.3 0.6 0.8 1.0 1.4 1.8 2.2 2.5; do EMT_HFRT_MIN=0.12 EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py \
        --limit-gb 9 -- .venv/bin/python tools/loop_cv.py --retake --tag slope --stage run --rt $rt --workers 8; done

hf_slope = the mean, over the same (HF, RT) cell's other LF steps, of the slope of (per-loop pole - baseline pole) against
(P_i - Pbar) / Pbar from lowband_shelf --stage hfloop (all captures): held out by LF step. No slope (0) in cells where a
baseline pole is within 0.005 of the 0.8475 ceiling (the ceiling clips every loop alike) and at HF 0.2 / RT 0.3 (its hfloop
poles were read against the pre-§29 base). G_capped adds hf_slope (P_i - Pbar) / Pbar to each loop's pole with the cell's hf_dp.
"""
import glob, json, os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_hfrt")
OUT = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_slope")
HL = [json.load(open(p)) for p in glob.glob(os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_hfloop_cap0.8475", "*.json"))]
P = np.array([1321, 1389, 1461, 1587.])
X = (P - P.mean()) / P.mean()
CEIL = 0.8475


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in glob.glob(os.path.join(SRC, "fold_*.json")):
        if not f.endswith(("_rows.json",)):
            open(os.path.join(OUT, os.path.basename(f)), "w").write(open(f).read())
    cells = json.load(open(os.path.join(SRC, "cells.json")))
    noslope = {(h["hf"], h["rt"]) for h in HL if max(h["p_base"]) > CEIL - 0.005 or h["hf"] * h["rt"] < 0.12 - 1e-9}
    slope = {h["key"]: float(np.polyfit(X, np.array(h["all"]["p"]) - np.array(h["p_base"]), 1)[0]) for h in HL}
    n = 0
    for h in HL:
        if (h["hf"], h["rt"]) in noslope:
            continue
        o = [slope[g["key"]] for g in HL if g["hf"] == h["hf"] and g["rt"] == h["rt"] and g["lf"] != h["lf"]]
        if o and h["key"] in cells and "hf_dp" in cells[h["key"]]:
            cells[h["key"]]["hf_slope"] = round(float(np.mean(o)), 5)
            n += 1
    json.dump(cells, open(os.path.join(OUT, "cells.json"), "w"), indent=1)
    print(f"{n} settings with a slope; no slope in {sorted(noslope)}; -> {OUT}")
    for hf in (0.2, 0.4, 0.6):
        print(f"  HF {hf}: " + "  ".join(f"RT {rt}: " + ",".join(f"{1e3 * cells[k].get('hf_slope', 0):+.1f}" for k in sorted(cells) if k.endswith(f"_{rt}_{hf}"))
                                         for rt in (0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5)))


if __name__ == "__main__":
    main()
