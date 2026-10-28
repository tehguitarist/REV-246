"""Repeat-capture level for N per band, incl. 20-100 Hz: captures at SPAN 0.5 and EARLY <= 0.3 differ only in LATE,
whose level is -145...-97 dB re a main tap there (late_db), so they are repeats of one early section."""
import os, sys, json, itertools
from pathlib import Path
import numpy as np
from scipy.signal import butter, sosfiltfilt
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
import nl_model as nm, nl_early as ne, nl_structure as nl
FS, K0 = nm.FS, nm.K0
tab = nl.table()
bands = [(20, 100), (100, 1000), (1000, 5000), (5000, 10000)]
span = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
rows = {b: [] for b in bands}; mrows = {b: [] for b in bands}
for early in (0.1, 0.2, 0.3):
    keys = sorted(k for k in tab if k[0] == span and k[1] == early and k not in nm.defects.NL_DEFECTS)
    print(span, early, "captures:", len(keys), flush=True)
    sig = {}
    for k in keys:
        raw = nl.ir(tab[k]); h = sosfiltfilt(ne.BAND, raw, axis=0); fr = ne.front_shift(h)
        sig[k] = (nl.ir(tab[k]), fr)
    for ka, kb in itertools.combinations(keys, 2):
        (ya, fa), (yb, fb) = sig[ka], sig[kb]
        n = min(len(ya), len(yb)) - 200
        a = K0 + int(fa) + int(0.05 * FS); b = min(K0 + int(fa) + int(0.45 * FS), n - 100)
        sh = int(round(fb - fa)); yb2 = np.roll(yb, -sh, axis=0)
        for lo, hi in bands:
            sos = butter(4, [lo, hi], btype="band", fs=FS, output="sos")
            A = sosfiltfilt(sos, ya[:n], axis=0)[a:b]
            Bfull = sosfiltfilt(sos, yb2[:n], axis=0)
            t0 = np.arange(a - 60, b + 60); out = []
            for c in (0, 1):
                cs = CubicSpline(np.arange(len(Bfull)), Bfull[:, c])
                tt = np.arange(a, b)

                def cost(p):
                    Bs = cs(tt + p[0] + p[1] * 1e-6 * (tt - a)); g = (A[:, c] @ Bs) / (Bs @ Bs)
                    return np.sum((A[:, c] - g * Bs) ** 2) / np.sum(A[:, c] ** 2)
                r = minimize(cost, [0.0, 0.0], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-10, "initial_simplex": [[0, 0], [0.05, 0], [0, 2.0]]})
                out.append(10 * np.log10(r.fun) - 3.0)
            rows[(lo, hi)].append(out)
        for lo, hi in bands:      # the model's own null on both captures, from the bench (0.05-0.6 s window)
            for kk in (ka, kb):
                f = os.path.join(nm.nl.OUT, "bench_bk_a0.65", "_".join(str(v) for v in kk) + ".json")
                if os.path.exists(f):
                    bd = json.load(open(f)).get("bands", {}).get(f"{lo}-{hi}")
                    if bd: mrows[(lo, hi)].append(bd)
print(f"SPAN {span}, EARLY 0.1/0.2/0.3, pairs differing only in LATE (window 0.05-0.45 s, gain + fractional shift fitted per band)")
print("band          repeat level (pairs, -3 dB)   model (bench alpha 0.65, 0.05-0.6 s)   model - repeat")
for bd in bands:
    r = np.array(rows[bd]); m = np.array(mrows[bd])
    print(f"{bd[0]:>5}-{bd[1]:<5} Hz   {np.median(r):6.1f}  (n {len(r)})            {np.median(m):6.1f}  (n {len(m)})            {np.median(m) - np.median(r):+5.1f}")
