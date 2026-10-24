"""Do the read HF coefficients sit on a fixed-point grid? Rayleigh test of c = a (1 - p), p and 1 - p (poles.py's
read values) against steps 2^-6 .. 2^-11; the baseline's own (smooth-law) values as a false-positive check.
model_s_results §25.14.

    .venv/bin/python tools/loop_read/grid.py"""
import json, os
import numpy as np
r = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                                "analysis/loop_fir/chain3_current_kc8_read/poles.json")))
for nm, sel in (("one-pole fits <= -8 dB", [x for x in r if x["fit_db"] <= -8]), ("all", r)):
    c = np.array([x["a"] * (1 - x["p"]) for x in sel]); p = np.array([x["p"] for x in sel]); cb = np.array([x["ab"] * (1 - x["pb"]) for x in sel])
    print(f"{nm}: n={len(sel)}; c = a(1-p) {c.min():.4f}..{c.max():.4f}; c/(1-p)-1 median {np.median(c / (1 - p) - 1):+.4f}")
    for lab, v in (("c read", c), ("p read", p), ("c baseline", cb)):
        line = []
        for n in range(6, 12):
            R = abs(np.mean(np.exp(2j * np.pi * v * 2 ** n)))
            line.append(f"2^-{n}: R {R:.2f} p~{np.exp(-len(v) * R ** 2):.2g}")
        print(f"   {lab:10s} " + "; ".join(line))
