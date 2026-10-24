"""Lag steps in the Tomkin S captures: per capture (dense_mid prep: the §25.7 baseline, clocks per capture), the lag
track of warp.py; the step = median lag 120-300 ms minus median lag 60-100 ms (samples), flagged at |step| >= 0.5,
and which captures of a setting share their late reverb (lag tracks equal to 0.05 samples). model_s_results §25.16.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/steps.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import dense_mid as dm
import warp as W
info = dm.q2.q1_tables()[0]
keys = sorted(p[:-4] for p in os.listdir(os.path.join(dm.OUT, "prep")))
flag = []
for k in keys:
    z = np.load(os.path.join(dm.OUT, "prep", k + ".npz"))
    row = []
    for j in sorted(int(n[3:]) for n in z.files if n.startswith("res")):
        y, m = z[f"y{j}"][:, 1].astype(float), z[f"m{j}"][:, 1].astype(float)
        i0, i1 = dm.PRE + int(0.06 * dm.FS), min(dm.PRE + int(min(float(z[f"te{j}"]), 0.3) * dm.FS), len(y) - W.MAXLAG)
        T, L, E = W.lag_track(y, m, dm.FS, i0, i1)
        t = T / dm.FS * 1000 + 60
        early, late = L[(t >= 60) & (t <= 100)], L[(t >= 120) & (t <= 300)]
        st = float(np.median(late) - np.median(early)) if len(early) and len(late) else np.nan
        row.append((int(z[f"amp{j}"]), st))
    bad = [(a, s) for a, s in row if abs(s) >= 0.5]
    print(f"{k}: " + " ".join(f"A{a}:{s:+.2f}" for a, s in row) + ("   <-- step" if bad else ""), flush=True)
    flag += [(k, a, s) for a, s in bad]
print(f"\n{len(flag)} captures with |step| >= 0.5 samples in {len({f[0] for f in flag})} of {len(keys)} settings; by AMPL: " +
      ", ".join(f"{a}: {sum(1 for f in flag if f[1] == a)}" for a in (0, 30, 60, 90)))
