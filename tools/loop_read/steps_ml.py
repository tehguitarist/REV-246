"""Lag tracks (warp.py) of the M / L captures against the law model (ml_lines.fitted: stored clocks and gains), per
capture: median lag in time bands and the fitted slope (ppm), at the strong-HF short-RT settings and high-HF controls.
model_s_results §25.16.

    .venv/bin/python tools/loop_read/steps_ml.py M|L key,key,..."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ml_lines as ml
import warp as W
prog, keys = sys.argv[1], sys.argv[2].split(",")
BANDS = [(60, 100), (120, 200), (200, 300), (300, 450), (450, 600)]
for k in keys:
    out = []
    for amp, (y, m, end) in ml.fitted(prog, k).items():
        i0, i1 = ml.K0 + int(0.06 * ml.FS), min(end[1], ml.K0 + int(0.6 * ml.FS), len(y) - W.MAXLAG)
        T, L, E = W.lag_track(y[:, 1], m[:, 1], ml.FS, i0, i1)
        t = T / ml.FS * 1000 + 60
        med = [np.median(L[(t >= a) & (t < b)]) if np.any((t >= a) & (t < b)) else np.nan for a, b in BANDS]
        out.append(f"A{amp}: " + "/".join(f"{x:+.2f}" for x in med) + f" ({W.slope(T, L, E)[0]:+.0f} ppm)")
    print(f"{prog} {k}: median lag (smp) at " + "/".join(f"{a}-{b}" for a, b in BANDS) + " ms | " + " | ".join(out), flush=True)
