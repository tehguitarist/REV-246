"""The common residual waveform (mean over undamped settings and captures, prep: y - model, 100 Hz-10 kHz): its envelope
in time, its spectrum re the model's in the same window, and its cross-correlation with the mean model output (an echo /
extra tap / missing stage shows as a peak at a fixed lag).
model_s_results §25.12.

    EMT_KERNEL=chain3 .venv/bin/python tools/resid_common/echo.py"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt
import dense_mid as dm
FS, PRE = dm.FS, dm.PRE
sos = butter(4, (100, 10000), btype="band", fs=FS, output="sos")
ap = json.load(open(os.path.join(dm.q2.OUT, dm.ls.lb.lf_.AP230)))
bstep = lambda k: ap["steps"][ap["settings"][k]["step"]]
info = dm.q2.q1_tables()[0]
n_end = PRE + int(0.16 * FS)
res, mod = {0: [], 1: []}, {0: [], 1: []}
for p in sorted(os.listdir(os.path.join(dm.OUT, "prep"))):
    k = p[:-4]
    if bstep(k) >= 0.01 or info[k]["rt"] < 1.0:
        continue
    z = np.load(os.path.join(dm.OUT, "prep", p))
    for j in sorted(int(n[3:]) for n in z.files if n.startswith("res")):
        y, m = z[f"y{j}"].astype(float), z[f"m{j}"].astype(float)
        if len(y) < n_end + 2000 or int(z[f"amp{j}"]) == 0:
            continue
        for c in (0, 1):
            res[c].append(sosfiltfilt(sos, (y - m)[:n_end + 2000, c])[:n_end]); mod[c].append(sosfiltfilt(sos, m[:n_end + 2000, c])[:n_end])
for c in (0, 1):
    r, m = np.mean(res[c], axis=0), np.mean(mod[c], axis=0)
    print(f"== {'LR'[c]}: {len(res[c])} captures (undamped, RT >= 1.0, AMPL > 0)")
    t = lambda ms: PRE + int(ms / 1000 * FS)
    env = [10 * np.log10(np.sum(r[t(a):t(a + 10)] ** 2) / np.sum(m[t(a):t(a + 10)] ** 2)) for a in range(0, 150, 10)]
    print("  residual / model per 10 ms from 0 ms (dB): " + " ".join(f"{e:.0f}" for e in env))
    i0, i1 = t(30), t(100)
    R_, M_ = np.fft.rfft(r[i0:i1] * np.hanning(i1 - i0)), np.fft.rfft(m[i0:i1] * np.hanning(i1 - i0))
    f = np.fft.rfftfreq(i1 - i0, 1 / FS)
    print("  30-100 ms, residual / model by octave (dB): " + " ".join(f"{fc}:{10 * np.log10(np.sum(np.abs(R_[(f > fc / 1.41) & (f < fc * 1.41)]) ** 2) / np.sum(np.abs(M_[(f > fc / 1.41) & (f < fc * 1.41)]) ** 2)):.0f}" for fc in (125, 250, 500, 1000, 2000, 4000, 8000)))
    # cross-correlation of the residual (30-100 ms) with the model (whole 0-160 ms), lag = residual time - model time
    seg = r[i0:i1]
    L = 3000
    xc = np.array([seg @ m[i0 - s:i1 - s] / np.sqrt(seg @ seg * (m[i0 - s:i1 - s] @ m[i0 - s:i1 - s])) if 0 <= i0 - s and i1 - s <= len(m) else 0 for s in range(-L, L + 1)])
    top = np.argsort(-np.abs(xc))[:8]
    print("  residual vs model shifted by s samples (residual = a m(t - s)), top |corr|: " + ", ".join(f"{s - L:+d} ({xc[s]:+.2f})" for s in top))
    # the same against the derivative of the model (a delay / phase error)
    md = np.gradient(m)
    print(f"  at s = 0: corr with model {xc[L]:+.3f}, with its derivative {seg @ md[i0:i1] / np.sqrt(seg @ seg * (md[i0:i1] @ md[i0:i1])):+.3f}")
