"""S part of the full-band gate with the shared sub-100 Hz transfer (gate_lowband) against the reference (gate_ref) and the control (gate_lbcontrol = switch off):
median null (dB) per band and window over the 16 gate settings' captures, L / R."""
import json, os, sys
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
ROOT = "/Users/leighpierce/Git/246/analysis/fullband"
GATE_S = [f"{lf}_{rt}_{hf}" for i, rt in enumerate((0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5)) for lf, hf in []]
import pulse_check as pc
EDGES = [20, 40, 80, 160, 270, 400, 800, 1600, 3200, 6400, 10000, 11000, 12300]; BN = [f"{a}-{b}" for a, b in zip(EDGES[:-1], EDGES[1:])]
def db(n, d): return 10 * np.log10(n / d) if n > 0 and d > 0 else None
def val(root, b, w, ch):
    xs = [db(c["bands"][b][w]["e"][ch], c["bands"][b][w]["y"][ch]) for k in pc.GATE_S for c in json.load(open(os.path.join(root, "s", k + ".json")))["caps"].values() if w in c["bands"].get(b, {})]
    xs = [x for x in xs if x is not None]
    return np.median(xs) if xs else np.nan
roots = {n: os.path.join(ROOT, d) for n, d in (("ref", "gate_ref"), ("control", "gate_lbcontrol"), ("lowband", "gate_lowband")) if os.path.isdir(os.path.join(ROOT, d, "s"))}
print("roots:", list(roots))
for w in ("50-150", "150-300", "300-600"):
    print(f"\nS {w} ms: median null per band (dB), L / R:   ref | control | lowband   (lowband - ref)")
    for b in BN[:8]:
        row = [f"{val(r, b, w, 0):6.1f}/{val(r, b, w, 1):6.1f}" for r in roots.values()]
        d = [val(roots["lowband"], b, w, c) - val(roots["ref"], b, w, c) for c in (0, 1)]
        print(f"  {b:>9s}  " + " | ".join(row) + f"   ({d[0]:+.2f} / {d[1]:+.2f})")
