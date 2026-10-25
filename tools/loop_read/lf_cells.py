"""The (LF, RT) cells' per-loop shelf terms as a baseline candidate (model_s_results 32.2): a variant fold set for
loop_cv --retake --tag lfcell, the §32.1 baseline's folds and cells with each RT >= 1.0 setting's lf_gain / lf_dc added.

    .venv/bin/python tools/loop_read/lf_cells.py [--zero-mean-gain]   # -> analysis/loop_cv/hfcells_chain3_cap0.8475_lfcell[zm]/
    .venv/bin/python tools/loop_read/lf_cells.py --low [--zero-mean-gain]   # -> ..._lowrt[zm]/: the same at RT <= 0.8 on _lfcellzm
    for rt in 1.0 1.4 1.8 2.2 2.5; do EMT_HFRT_MIN=0.12 EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py \
        --limit-gb 9 -- .venv/bin/python tools/loop_cv.py --retake --tag lfcell --stage run --rt $rt --workers 8; done
    (RT 0.3-0.8 unchanged: copy their _rows.json from the slope fold set)

lf_gain / lf_dc = each loop's broadband gain (dB per pass) and DC section (dB, the shelf's pole), the mean over the same
(LF, RT) cell's other HF steps of lowband_shelf --stage takefit on §32.1 (both takes): held out by HF step. Only at RT >= 1.0,
where §32 found the corrections separable into LF and HF parts (at RT <= 0.8 they are not, and the cell means hurt, §32.2).

--low (model_s_results 32.3): the same terms at RT <= 0.8 on top of the §32.2 baseline, taken only from the cell's HF steps
outside the pole-limited regime (HF 0.2, and HF 0.4 at RT 0.3, where the ceiling / HF floor act and the takefit gain trades
with the pole): with those left out, the cell's other HF steps predict the shelf side there as well as at RT >= 1.0.

--delta NEW OLD --dir D (model_s_results 34): refit the shelf cells of fold set D in place after a model change below them
(the DC floor's compact form): each setting with a NEW takefit gets lf_gain / lf_dc + (mean of NEW - mean of OLD) over the
same sources as above (the cell's other HF steps; --low without the pole-limited ones; --zero-mean-gain: the gain's loop
mean removed from the change). NEW and OLD must be takefit runs that differ only in the loop responses (SHELF_PREP_CELLS:
the captures fitted on the same model), so the change is the cells' share of the model change and nothing else. Run it
twice, as the two halves above:
    .venv/bin/python tools/loop_read/lf_cells.py --zero-mean-gain --delta <NEW> <OLD> --dir <D>
    .venv/bin/python tools/loop_read/lf_cells.py --low --delta <NEW> <OLD> --dir <D>
"""
import glob, json, os, sys

import numpy as np

LOW = "--low" in sys.argv
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_" + ("lfcellzm" if LOW else "slope"))
ZM = "--zero-mean-gain" in sys.argv          # lf_gain minus its loop mean: takefit keeps each capture's gain fixed, so a common
                                             # broadband gain can be capture-gain error (NEXT_STEPS 0, rules); -> _lfcellzm
OUT = os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_" + ("lowrt" if LOW else "lfcell") + ("zm" if ZM else ""))
TF = [json.load(open(p)) for p in glob.glob(os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_takefit_cap0.8475_s321", "*.json"))]
RT_MIN = 1.0
POLE = lambda c: c["hf"] == 0.2 or (c["rt"] == 0.3 and c["hf"] == 0.4)   # --low: pole-limited HF steps, not used as sources
in_set = (lambda rt: rt < RT_MIN) if LOW else (lambda rt: rt >= RT_MIN)


def delta(new, old, d):
    """--delta: the cells of fold set d moved by the change between two takefit runs (model_s_results 34)."""
    load = lambda p: {t["key"]: t for t in (json.load(open(f)) for f in glob.glob(os.path.join(p, "*.json")))}
    N, O = load(new), load(old)
    cells = json.load(open(os.path.join(d, "cells.json")))
    fit = lambda c: np.mean([np.array(c["fit"][k][:8]) for k in ("A", "B") if k in c["fit"]], 0)
    n, log = 0, []
    for key, t in sorted(N.items()):
        if not in_set(t["rt"]) or "lf_gain" not in cells.get(key, {}):
            continue
        src = [k for k, c in N.items() if c["lf"] == t["lf"] and c["rt"] == t["rt"] and c["hf"] != t["hf"] and c["fit"]
               and k in O and O[k]["fit"] and not (LOW and POLE(c))]
        if not src:
            continue
        u = np.mean([fit(N[k]) - fit(O[k]) for k in src], 0)
        if ZM:
            u[:4] -= u[:4].mean()
        c = cells[key]
        c["lf_gain"] = [round(float(x + y), 5) for x, y in zip(c["lf_gain"], u[:4])]
        c["lf_dc"] = [round(float(x + y), 5) for x, y in zip(c["lf_dc"], u[4:8])]
        log.append(f"  {key}: gain " + " ".join(f"{x:+.4f}" for x in u[:4]) + "  dc " + " ".join(f"{x:+.4f}" for x in u[4:8]))
        n += 1
    json.dump(cells, open(os.path.join(d, "cells.json"), "w"), indent=1)
    print(f"{n} settings' shelf cells moved (RT {'<' if LOW else '>='} {RT_MIN}{', zero-mean gain change' if ZM else ''}) -> {d}/cells.json")
    print("\n".join(log))


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in glob.glob(os.path.join(SRC, "fold_*.json")):
        rt = float(os.path.basename(f).split("_")[1].replace(".json", ""))
        if not f.endswith("_rows.json") or not in_set(rt):      # folds outside the set are unchanged: their rows too
            open(os.path.join(OUT, os.path.basename(f)), "w").write(open(f).read())
    cells = json.load(open(os.path.join(SRC, "cells.json")))
    n = 0
    for t in TF:
        if not in_set(t["rt"]):
            continue
        o = [np.mean([np.array(c["fit"][k][:8]) for k in ("A", "B") if k in c["fit"]], 0)
             for c in TF if c["lf"] == t["lf"] and c["rt"] == t["rt"] and c["hf"] != t["hf"] and c["fit"] and not (LOW and POLE(c))]
        if o:
            u = np.mean(o, 0)
            if ZM:
                u[:4] -= u[:4].mean()
            cells.setdefault(t["key"], {}).update(lf_gain=[round(float(x), 5) for x in u[:4]], lf_dc=[round(float(x), 5) for x in u[4:8]])
            n += 1
    json.dump(cells, open(os.path.join(OUT, "cells.json"), "w"), indent=1)
    print(f"{n} settings with per-loop shelf terms (RT {'<' if LOW else '>='} {RT_MIN}) -> {OUT}")
    for lf in (0.1, 0.4, 1.0, 2.0):
        for rt in ((0.3, 0.6, 0.8) if LOW else (1.0, 1.4, 1.8, 2.2, 2.5)):
            k = next(k for k in sorted(cells) if k.startswith(f"{lf}_{rt}_") and "lf_gain" in cells[k])
            print(f"  LF {lf} RT {rt} ({k}): gain " + " ".join(f"{x:+.4f}" for x in cells[k]["lf_gain"]) + "  dc " + " ".join(f"{x:+.4f}" for x in cells[k]["lf_dc"]))


if __name__ == "__main__":
    if "--delta" in sys.argv:
        i = sys.argv.index("--delta")
        delta(sys.argv[i + 1], sys.argv[i + 2], sys.argv[sys.argv.index("--dir") + 1])
    else:
        main()
