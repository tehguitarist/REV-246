"""Report for lowband_shelf.py --stage dcxfer and addxfer (model_s_results 32): the per-(LF, RT, loop) DC table held out
by HF step and by take, and the additive LF + HF table per (RT, loop) held out by setting and by take, against the §29
baseline and takefit's per-setting 12-number refit (same captures), R and L medians.

    .venv/bin/python tools/loop_read/dcxfer_report.py [dcxfer_dir] [takefit_dir] [addxfer_dir]
"""
import glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_dcxfer_cap0.8475_s29")
TF = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_takefit_cap0.8475_s29")
WIN = ["50-150", "150-300", "300-600"]
D = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(DX, "*.json")))]
AX = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_addxfer_cap0.8475_s29")
X = {json.load(open(p))["key"]: json.load(open(p)) for p in glob.glob(os.path.join(AX, "*.json"))}
T = {json.load(open(p))["key"]: json.load(open(p)) for p in glob.glob(os.path.join(TF, "*.json"))}


def med(sel, src, w, ch):
    v = []
    for d in sel:
        tab = {"base": d["base"], "loop": d["held"].get("loop", {}), "cell": d["held"].get("cell", {}),
               "takefit": T.get(d["key"], {}).get("held", {}), "add": X.get(d["key"], {}).get("held", {})}[src]
        amps = set(d["held"].get("loop", {})) & set(T.get(d["key"], {}).get("held", {})) & set(X.get(d["key"], {}).get("held", {}))   # the same captures for every column
        v += [tab[a][w][ch] for a in amps if w in tab.get(a, {})]
    return np.median(v) if v else np.nan, len(v)


groups = [("all 128", lambda d: True), ("LF 0.1", lambda d: d["lf"] == 0.1), ("LF 0.4", lambda d: d["lf"] == 0.4),
          ("LF 1.0", lambda d: d["lf"] == 1.0), ("LF 2.0", lambda d: d["lf"] == 2.0), ("RT <= 1.0", lambda d: d["rt"] <= 1.0), ("RT 1.0", lambda d: d["rt"] == 1.0),
          ("RT >= 1.4", lambda d: d["rt"] >= 1.4), ("HF 0.2", lambda d: d["hf"] == 0.2), ("LF 0.1 RT 1.4-1.8", lambda d: d["lf"] == 0.1 and d["rt"] in (1.4, 1.8)),
          ("LF 0.4 RT 0.6", lambda d: d["lf"] == 0.4 and d["rt"] == 0.6)]
for ch, name in ((1, "R"), (0, "L")):
    print(f"\n{name}: base (§29) -> DC table per loop | per cell (loop mean) | additive LF + HF | takefit 12/setting, median dB over held-out captures [n]")
    for g, f in groups:
        sel = [d for d in D if f(d)]
        cells = []
        for w in WIN:
            b, n = med(sel, "base", w, ch)
            cells.append(f"{w}: {b:.1f} -> {med(sel, 'loop', w, ch)[0]:.1f} | {med(sel, 'cell', w, ch)[0]:.1f} | {med(sel, 'add', w, ch)[0]:.1f} | {med(sel, 'takefit', w, ch)[0]:.1f} [{n}]")
        print(f"  {g:18s} " + "   ".join(cells))
# per setting gains and losses at 150-300 / 300-600 (R)
for form, src in (("DC table", lambda d: d["held"].get("loop", {})), ("additive", lambda d: X.get(d["key"], {}).get("held", {}))):
    for w in ("150-300", "300-600"):
        dl = []
        for d in D:
            h = src(d)
            amps = [a for a in h if w in h[a] and w in d["base"][a]]
            if amps:
                dl.append((np.mean([h[a][w][1] - d["base"][a][w][1] for a in amps]), d["key"]))
        dl.sort()
        print(f"\nR {w}: per-setting change ({form} - base), dB: better > 0.5: {sum(x < -0.5 for x, _ in dl)}, worse > 0.5: {sum(x > 0.5 for x, _ in dl)}; "
              f"best {', '.join(f'{k} {x:+.1f}' for x, k in dl[:4])}; worst {', '.join(f'{k} {x:+.1f}' for x, k in dl[-4:])}")

# the combination: additive LF + HF at RT >= 1.0, the DC table below (the split the per-frequency check gives, §32)
print("\ncombined (additive at RT >= 1.0, DC table at RT <= 0.8): R / L median, held out by setting and take; per setting R")
for g, f in groups[:8]:
    sel = [d for d in D if f(d)]
    for ch, name in ((1, "R"), (0, "L")):
        cells = []
        for w in WIN:
            vb, vc = [], []
            for d in sel:
                h = X.get(d["key"], {}).get("held", {}) if d["rt"] >= 1.0 else d["held"].get("loop", {})
                amps = [a for a in set(h) & set(T.get(d["key"], {}).get("held", {})) & set(X.get(d["key"], {}).get("held", {})) if w in h[a]]
                vb += [d["base"][a][w][ch] for a in amps]
                vc += [h[a][w][ch] for a in amps]
            cells.append(f"{w}: {np.median(vb):.1f} -> {np.median(vc):.1f}")
        print(f"  {g:12s} {name} " + "   ".join(cells))
for w in ("150-300", "300-600"):
    dl = []
    for d in D:
        h = X.get(d["key"], {}).get("held", {}) if d["rt"] >= 1.0 else d["held"].get("loop", {})
        amps = [a for a in h if w in h[a] and w in d["base"][a]]
        if amps:
            dl.append((np.mean([h[a][w][1] - d["base"][a][w][1] for a in amps]), d["key"]))
    dl.sort()
    print(f"R {w}: better > 0.5: {sum(x < -0.5 for x, _ in dl)}, worse > 0.5: {sum(x > 0.5 for x, _ in dl)} of {len(dl)}; "
          f"worst {', '.join(f'{k} {x:+.1f}' for x, k in dl[-5:])}")

# the HF-pole slope across the loops per (HF, RT) cell (slopexfer), held out by LF step and by take
SX = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_slopexfer_cap0.8475_s29")
SL = {json.load(open(p))["key"]: json.load(open(p)) for p in glob.glob(os.path.join(SX, "*.json"))}
# No slope where it has no meaning: cells with a baseline pole within 0.005 of the 0.8475 ceiling (the ceiling clips every
# loop alike), and HF 0.2 / RT 0.3, whose hfloop poles were read against the pre-§29 base (§29 models it as HF 0.4).
HL = [json.load(open(p)) for p in glob.glob(os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_hfloop_cap0.8475", "*.json"))]
NOSLOPE = {(h["hf"], h["rt"]) for h in HL if max(h["p_base"]) > 0.8475 - 0.005 or h["hf"] * h["rt"] < 0.12 - 1e-9}
for d in SL.values():
    if (d["hf"], d["rt"]) in NOSLOPE:
        d["held"] = {a: d["base"][a] for a in d["held"]}
if SL:
    print(f"\npole slope table (no slope in the cells at the ceiling and HF 0.2 / RT 0.3: {sorted(NOSLOPE)}; per (HF, RT) cell from its other LF steps), HF < 1: base -> slope; R / L median")
    for g, f in [("HF < 1", lambda d: True), ("HF 0.2", lambda d: d["hf"] == 0.2), ("HF 0.4", lambda d: d["hf"] == 0.4),
                 ("HF 0.6", lambda d: d["hf"] == 0.6), ("RT <= 1.0", lambda d: d["rt"] <= 1.0), ("RT >= 1.4", lambda d: d["rt"] >= 1.4)]:
        sel = [d for d in SL.values() if f(d)]
        for ch, name in ((1, "R"), (0, "L")):
            cells = []
            for w in WIN:
                vb = [d["base"][a][w][ch] for d in sel for a in d["held"] if w in d["held"][a]]
                vs = [d["held"][a][w][ch] for d in sel for a in d["held"] if w in d["held"][a]]
                cells.append(f"{w}: {np.median(vb):.1f} -> {np.median(vs):.1f}")
            print(f"  {g:10s} {name} " + "   ".join(cells) + f"  [{len(vb)}]")
    for w in ("150-300", "300-600"):
        dl = sorted((np.mean([d["held"][a][w][1] - d["base"][a][w][1] for a in d["held"] if w in d["held"][a]]), d["key"])
                    for d in SL.values() if any(w in d["held"][a] for a in d["held"]))
        print(f"R {w}: better > 0.5: {sum(x < -0.5 for x, _ in dl)}, worse > 0.5: {sum(x > 0.5 for x, _ in dl)} of {len(dl)}; "
              f"best {', '.join(f'{k} {x:+.1f}' for x, k in dl[:4])}; worst {', '.join(f'{k} {x:+.1f}' for x, k in dl[-4:])}")
    print("slopes s (take A) per cell (x 1e-3):")
    for hf in (0.2, 0.4, 0.6):
        print(f"  HF {hf}: " + "  ".join(f"RT {rt}: {1e3 * np.mean([d['s']['A'] for d in SL.values() if d['hf'] == hf and d['rt'] == rt and 'A' in d['s']]):+.1f}"
                                         for rt in sorted({d["rt"] for d in SL.values() if d["hf"] == hf})))
