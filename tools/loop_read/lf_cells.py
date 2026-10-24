"""The (LF, RT) cells' per-loop shelf terms as a baseline candidate (model_s_results 32.2): a variant fold set for
loop_cv --retake --tag lfcell, the §32.1 baseline's folds and cells with each RT >= 1.0 setting's lf_gain / lf_dc added.

    .venv/bin/python tools/loop_read/lf_cells.py [--zero-mean-gain]   # -> analysis/loop_cv/hfcells_chain3_cap0.8475_lfcell[zm]/
    for rt in 1.0 1.4 1.8 2.2 2.5; do EMT_HFRT_MIN=0.12 EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py \
        --limit-gb 9 -- .venv/bin/python tools/loop_cv.py --retake --tag lfcell --stage run --rt $rt --workers 8; done
    (RT 0.3-0.8 unchanged: copy their _rows.json from the slope fold set)

lf_gain / lf_dc = each loop's broadband gain (dB per pass) and DC section (dB, the shelf's pole), the mean over the same
(LF, RT) cell's other HF steps of lowband_shelf --stage takefit on §32.1 (both takes): held out by HF step. Only at RT >= 1.0,
where §32 found the corrections separable into LF and HF parts (at RT <= 0.8 they are not, and the cell means hurt, §32.2).
"""
import glob, json, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_slope")
ZM = "--zero-mean-gain" in sys.argv          # lf_gain minus its loop mean: takefit keeps each capture's gain fixed, so a common
                                             # broadband gain can be capture-gain error (NEXT_STEPS 0, rules); -> _lfcellzm
OUT = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_lfcell" + ("zm" if ZM else ""))
TF = [json.load(open(p)) for p in glob.glob(os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_takefit_cap0.8475_s321", "*.json"))]
RT_MIN = 1.0


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in glob.glob(os.path.join(SRC, "fold_*.json")):
        rt = float(os.path.basename(f).split("_")[1].replace(".json", ""))
        if not f.endswith("_rows.json") or rt < RT_MIN:          # RT < 1.0 folds are unchanged: their rows too
            open(os.path.join(OUT, os.path.basename(f)), "w").write(open(f).read())
    cells = json.load(open(os.path.join(SRC, "cells.json")))
    n = 0
    for t in TF:
        if t["rt"] < RT_MIN:
            continue
        o = [np.mean([np.array(c["fit"][k][:8]) for k in ("A", "B") if k in c["fit"]], 0)
             for c in TF if c["lf"] == t["lf"] and c["rt"] == t["rt"] and c["hf"] != t["hf"] and c["fit"]]
        if o:
            u = np.mean(o, 0)
            if ZM:
                u[:4] -= u[:4].mean()
            cells.setdefault(t["key"], {}).update(lf_gain=[round(float(x), 5) for x in u[:4]], lf_dc=[round(float(x), 5) for x in u[4:8]])
            n += 1
    json.dump(cells, open(os.path.join(OUT, "cells.json"), "w"), indent=1)
    print(f"{n} settings with per-loop shelf terms (RT >= {RT_MIN}) -> {OUT}")
    for lf in (0.1, 0.4, 1.0, 2.0):
        for rt in (1.0, 1.4, 1.8, 2.2, 2.5):
            k = next(k for k in sorted(cells) if k.startswith(f"{lf}_{rt}_") and "lf_gain" in cells[k])
            print(f"  LF {lf} RT {rt} ({k}): gain " + " ".join(f"{x:+.4f}" for x in cells[k]["lf_gain"]) + "  dc " + " ".join(f"{x:+.4f}" for x in cells[k]["lf_dc"]))


if __name__ == "__main__":
    main()
