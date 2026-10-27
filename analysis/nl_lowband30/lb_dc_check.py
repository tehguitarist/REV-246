"""Spot check of the pipeline agent's claim: the pre-roll step c is proportional to the file's DC in its first 30 samples (c in 1e-4 of the 410-tap peak from oracle.json)."""
import json, glob, os
import numpy as np, soundfile as sf
from lb_common import *
O = json.load(open("oracle.json")); B = bench_all()
rows = []
for ks, v in list(O.items())[::6]:
    k = tuple(float(x) for x in ks.split("_"))
    x, fs = sf.read(nl.table()[k], always_2d=True)
    dc = x[:30].mean(0); pk = np.abs(x[:, :]).max(0)
    r = B[k]; z = int(K0 + r["front_shift"])
    raw = nl.ir(nl.table()[k]); pky = np.abs(raw[z:z + 800]).max(0)
    rows.append((dc[0] / pky[0], dc[1] / pky[1], v["c0"], v["c1"]))
A = np.array(rows)
for c in (0, 1):
    a, b = A[:, c] * 1.0, A[:, 2 + c]
    k = (a @ b) / (a @ a)
    print(f"ch{c}: n {len(A)} corr(DC/peak, c) {np.corrcoef(a, b)[0,1]:.3f}  slope c/(DC/peak) {k:.2f}  residual rms / c rms {np.sqrt(np.mean((b - k*a)**2)) / np.sqrt(np.mean(b**2)):.2f}")
