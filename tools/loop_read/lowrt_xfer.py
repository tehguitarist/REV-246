"""Do the per-loop shelf terms carry across HF steps at RT <= 0.8 once the pole-limited HF steps are left out? (model_s_results
32.3). Offline, on lowband_shelf --stage takefit on §32.1 (_s321) and takefit_read's per-pass changes: take B's per-pass change
per loop (all 12 numbers, at 125 Hz - 8 kHz) predicted by the shelf side (gain + DC section, 8 numbers) of take A's fits at the
same (LF, RT) cell's other HF steps, as explained variance, per RT step and for pole-limited / other targets. Source rules:
none, the mean of all other HF steps (§32.2's rule), the mean without the pole-limited steps (HF 0.2, and HF 0.4 at RT 0.3),
the median; each also with the gain's loop mean removed (/zm).

    .venv/bin/python tools/loop_read/lowrt_xfer.py      # seconds -> stdout
"""
import glob, json, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import model_s as m  # noqa: E402

FS = 19.2e6 / 728
F = np.array([125, 250, 500, 1000, 2000, 4000, 8000])
w = 2 * np.pi * F / FS
A_POLE = 0.9167                                     # the shelf's pole (law["a"])
LB = os.path.join(ROOT, "analysis", "lowband")
TF = {(t["lf"], t["rt"], t["hf"]): t for t in (json.load(open(p)) for p in glob.glob(os.path.join(LB, "chain3_shelf_takefit_cap0.8475_s321", "*.json")))}
RD = {}
for r in json.load(open(os.path.join(LB, "chain3_takefit_read_s321.json"))):
    RD.setdefault((r["lf"], r["rt"], r["hf"]), {})[r["loop"]] = r


def resp(u8):
    """Per-pass dB at F of each loop's gain u8[i] and DC section u8[4 + i] (loops x freqs)."""
    return np.array([u8[i] + 20 * np.log10(np.abs(m.shelf(w, 10 ** (u8[4 + i] / 20), 1.0, A_POLE))) for i in range(4)])


def zm(u):
    u = np.array(u, float)
    u[:4] -= u[:4].mean()
    return u


pole = lambda lf, rt, hf: hf == 0.2 or (rt == 0.3 and hf == 0.4)
VAR = {"none": lambda src: None,
       "all": lambda src: np.mean([u for _, u in src], 0),
       "nopole": lambda src: np.mean([u for k, u in src if not pole(*k)], 0) if any(not pole(*k) for k, _ in src) else None,
       "median": lambda src: np.median([u for _, u in src], 0)}

for rts in ((0.3,), (0.6,), (0.8,), (1.0, 1.4, 1.8, 2.2, 2.5)):
    for tgt in ("other", "pole"):
        out = {}
        for v, f in VAR.items():
            for z in (False, True):
                err, tot = [], []
                for (lf, rt, hf), t in TF.items():
                    if rt not in rts or (lf, rt, hf) not in RD or pole(lf, rt, hf) != (tgt == "pole"):
                        continue
                    for tr, te in (("A", "B"), ("B", "A")):
                        src = [((c["lf"], c["rt"], c["hf"]), np.array(c["fit"][tr][:8])) for k, c in TF.items()
                               if k[0] == lf and k[1] == rt and k[2] != hf and tr in c["fit"]]
                        if not src:
                            continue
                        u = f(src)
                        p = np.zeros((4, 7)) if u is None else resp(zm(u) if z else u)
                        y = np.array([RD[(lf, rt, hf)][lp][te] for lp in range(4)])
                        err.append((y - p) ** 2)
                        tot.append(y ** 2)
                if err:
                    out[v + ("/zm" if z else "")] = 1 - np.sum(err, 0).sum(0) / np.sum(tot, 0).sum(0)
        if out:
            print(f"RT {rts}, {tgt} targets: explained variance at {F.tolist()} Hz")
            for k, v in out.items():
                if k != "none/zm":
                    print(f"   {k:10}", " ".join(f"{x:+.2f}" for x in v))
