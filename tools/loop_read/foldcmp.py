"""Two loop_cv fold sets compared on clean captures (tools/defects.py), R and L medians and the excess over the chain3
repeat references (analysis/q2/nulls_free_chain3.json, §30), by group, and per setting (model_s_results 32).

    .venv/bin/python tools/loop_read/foldcmp.py [base_dir] [candidate_dir] [--per-setting]   # default: the §29 baseline -> the slope variant
--per-setting: also every setting's change in every window from 50 to 1000 ms, R and L, listing those beyond +-0.5 dB (34).
"""
import glob, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from defects import DEFECTS  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(HERE))
CV = os.path.join(ROOT, "analysis", "loop_cv")
A = sys.argv[1] if len(sys.argv) > 1 else os.path.join(CV, "hfcells_chain3_cap0.8475_hfrt")
B = sys.argv[2] if len(sys.argv) > 2 else os.path.join(CV, "hfcells_chain3_cap0.8475_slope")
WS = ["0-50", "50-150", "150-300", "300-600", "600-1000"]


def rows(d):
    return {(r["key"], r["amp"]): r for p in glob.glob(os.path.join(d, "fold_*_rows.json")) for r in json.load(open(p))["rows"]
            if (r["key"], r["amp"]) not in DEFECTS}


a, b = rows(A), rows(B)
rep = {(r["key"], r["amp"]): r for r in json.load(open(os.path.join(ROOT, "analysis", "q2", "nulls_free_chain3.json")))}
ks = sorted(set(a) & set(b))
rts = sorted({a[k]["rt"] for k in ks})
print(f"{os.path.basename(A)} -> {os.path.basename(B)}: {len({k[0] for k in ks})} settings, {len(ks)} clean captures, RT folds {rts}")
groups = [("All", lambda r: True), ("Flat", lambda r: r["lf"] == 1.0 and r["hf"] == 1.0), ("LF 0.1", lambda r: r["lf"] == 0.1),
          ("LF 0.4", lambda r: r["lf"] == 0.4),
          ("HF 0.2", lambda r: r["hf"] == 0.2), ("HF 0.4", lambda r: r["hf"] == 0.4), ("HF 0.6", lambda r: r["hf"] == 0.6),
          ("RT <= 1.0", lambda r: r["rt"] <= 1.0), ("RT >= 1.4", lambda r: r["rt"] >= 1.4)]
for ch, nm in ((1, "R"), (0, "L")):
    print(f"\n{nm}: median null [excess over repeat]")
    print("| Settings | " + " | ".join(WS) + " |\n|---|" + "---|" * len(WS))
    for g, f in groups:
        kk = [k for k in ks if f(a[k])]
        cells = []
        for w in WS:
            va = [a[k]["null_comp"][w][ch] for k in kk if a[k]["null_comp"].get(w) and b[k]["null_comp"].get(w)]
            vb = [b[k]["null_comp"][w][ch] for k in kk if a[k]["null_comp"].get(w) and b[k]["null_comp"].get(w)]
            ea = [a[k]["null_comp"][w][ch] - rep[k]["repeat_comp"][w][ch] for k in kk if a[k]["null_comp"].get(w) and k in rep and rep[k]["repeat_comp"].get(w)]
            eb = [b[k]["null_comp"][w][ch] - rep[k]["repeat_comp"][w][ch] for k in kk if b[k]["null_comp"].get(w) and k in rep and rep[k]["repeat_comp"].get(w)]
            cells.append(f"{np.median(va):.1f} [{np.median(ea):+.1f}] → {np.median(vb):.1f} [{np.median(eb):+.1f}]" if va else "")
        print(f"| {g} | " + " | ".join(cells) + " |")
for w in ("150-300", "300-600"):
    per = {}
    for k in ks:
        if a[k]["null_comp"].get(w) and b[k]["null_comp"].get(w):
            per.setdefault(k[0], []).append(b[k]["null_comp"][w][1] - a[k]["null_comp"][w][1])
    dl = sorted((np.mean(v), s) for s, v in per.items())
    print(f"\nR {w}, per setting: better > 0.5 dB {sum(x < -0.5 for x, _ in dl)}, worse > 0.5 {sum(x > 0.5 for x, _ in dl)} of {len(dl)}; "
          f"best {', '.join(f'{s} {x:+.1f}' for x, s in dl[:5])}; worst {', '.join(f'{s} {x:+.1f}' for x, s in dl[-5:])}")
if "--per-setting" in sys.argv:                     # every window 50-1000 ms, R and L: each setting's change (mean over its
    print("\nPer setting, mean change over the setting's clean captures (dB), every setting beyond +-0.5 dB:")   # clean captures)
    for ch, nm in ((1, "R"), (0, "L")):
        for w in WS[1:]:
            per = {}
            for k in ks:
                if a[k]["null_comp"].get(w) and b[k]["null_comp"].get(w):
                    per.setdefault(k[0], []).append(b[k]["null_comp"][w][ch] - a[k]["null_comp"][w][ch])
            dl = sorted((float(np.mean(v)), s_) for s_, v in per.items())
            out = [f"{s_} {x:+.2f}" for x, s_ in dl if abs(x) > 0.5]
            print(f"  {nm} {w}: max worse {dl[-1][1]} {dl[-1][0]:+.2f}, max better {dl[0][1]} {dl[0][0]:+.2f}; beyond +-0.5: " + (", ".join(out) or "none"))
