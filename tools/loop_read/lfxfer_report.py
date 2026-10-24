"""Report for lowband_shelf.py --stage lfxfer (model_s_results 32.2): each loop's gain and DC section per (LF, RT) cell from
the cell's other HF steps, as per-loop means and as common + slope across the loops, held out by HF step and by take,
against the §32.1 baseline and takefit's per-setting refit on it (same captures), R and L medians; per-setting tails; the
cells where the form helps or hurts.

    .venv/bin/python tools/loop_read/lfxfer_report.py      # analysis/lowband/chain3_lfxfer_report.log
"""
import glob, json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LX = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_lfxfer_cap0.8475_s321")
TF = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_takefit_cap0.8475_s321")
WIN = ["50-150", "150-300", "300-600"]
D = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(LX, "*.json")))]
T = {json.load(open(p))["key"]: json.load(open(p)) for p in glob.glob(os.path.join(TF, "*.json"))}
out = []


def vals(sel, src, w, ch):
    v = []
    for d in sel:
        h = {"base": d["base"], "loop": d["held"].get("loop", {}), "cs": d["held"].get("cs", {}), "takefit": T.get(d["key"], {}).get("held", {})}[src]
        amps = set(d["held"].get("cs", {})) & set(T.get(d["key"], {}).get("held", {}))
        v += [h[a][w][ch] for a in amps if w in h.get(a, {})]
    return v


groups = [("all 128", lambda d: True), ("LF 0.1", lambda d: d["lf"] == 0.1), ("LF 0.4", lambda d: d["lf"] == 0.4), ("LF 1.0", lambda d: d["lf"] == 1.0),
          ("LF 2.0", lambda d: d["lf"] == 2.0), ("RT <= 0.8", lambda d: d["rt"] <= 0.8), ("RT >= 1.0", lambda d: d["rt"] >= 1.0),
          ("LF 0.1 RT >= 1.4", lambda d: d["lf"] == 0.1 and d["rt"] >= 1.4), ("HF 0.2", lambda d: d["hf"] == 0.2)]
for ch, name in ((1, "R"), (0, "L")):
    out.append(f"\n{name}: §32.1 base -> per-loop means (8 / cell) | common + slope (4 / cell) | takefit 12 / setting; median dB [n]")
    for g, f in groups:
        sel = [d for d in D if f(d)]
        cells = [f"{w}: {np.median(vals(sel, 'base', w, ch)):.1f} -> {np.median(vals(sel, 'loop', w, ch)):.1f} | {np.median(vals(sel, 'cs', w, ch)):.1f} | "
                 f"{np.median(vals(sel, 'takefit', w, ch)):.1f} [{len(vals(sel, 'base', w, ch))}]" for w in WIN]
        out.append(f"  {g:17s} " + "   ".join(cells))
for form in ("cs", "loop"):
    for w in ("150-300", "300-600"):
        dl = []
        for d in D:
            h = d["held"].get(form, {})
            amps = [a for a in h if w in h[a] and w in d["base"][a]]
            if amps:
                dl.append((np.mean([h[a][w][1] - d["base"][a][w][1] for a in amps]), d["key"]))
        dl.sort()
        out.append(f"\nR {w} ({form}): better > 0.5 dB {sum(x < -0.5 for x, _ in dl)}, worse > 0.5 {sum(x > 0.5 for x, _ in dl)} of {len(dl)}; "
                   f"best {', '.join(f'{k} {x:+.1f}' for x, k in dl[:4])}; worst {', '.join(f'{k} {x:+.1f}' for x, k in dl[-6:])}")
out.append("\nper (LF, RT) cell, R 300-600 mean change (cs) and the cell's numbers (take A, dB per pass: gain common / slope, DC common / slope):")
x = (np.array([1321, 1389, 1461, 1587.]) - 1439.5) / 1439.5
for lf in (0.1, 0.4, 1.0, 2.0):
    row = []
    for rt in (0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5):
        sel = [d for d in D if d["lf"] == lf and d["rt"] == rt]
        ch = [np.mean([d["held"]["cs"][a]["300-600"][1] - d["base"][a]["300-600"][1] for a in d["held"].get("cs", {}) if "300-600" in d["held"]["cs"][a]])
              for d in sel if d["held"].get("cs")]
        ch = [c for c in ch if np.isfinite(c)]
        u = np.array(sel[0]["u"]["A"]["cs"]) if sel and "A" in sel[0]["u"] else None
        nums = "" if u is None else f" g {u[:4].mean():+.3f}/{np.polyfit(x, u[:4], 1)[0]:+.3f} dc {u[4:].mean():+.3f}/{np.polyfit(x, u[4:], 1)[0]:+.3f}"
        row.append(f"RT {rt}: {np.mean(ch) if ch else float('nan'):+.1f}{nums}")
    out.append(f"  LF {lf}: " + "; ".join(row))
open(os.path.join(ROOT, "analysis", "lowband", "chain3_lfxfer_report.log"), "w").write("\n".join(out).lstrip() + "\n")
print("\n".join(out))
