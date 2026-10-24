"""S's pole curve applied to M / L with nothing fitted on them (model_s_results §25.17): per setting the capture-vs-model
lag drift (300-600 ms minus 60-100 ms, mean over captures, R) with the law's capped poles and with the curve's poles
(polecurve.py's all-RT fit), and the residual energy re the capture at 150-600 ms (177 Hz - 2.8 kHz, stored clocks).

    .venv/bin/python tools/loop_read/ml_curve.py M|L key,key,..."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt
import ml_lines as ml
import warp as W
prog, keys = sys.argv[1], sys.argv[2].split(",")
pc = json.load(open(os.path.join(ROOT, "analysis/lowband/chain3_polecurve.json")))
kn, th = np.array(pc["knots"]), np.array(pc["all"])
curve = lambda x: min(x + float(np.interp(x, kn, th[:len(kn)])), float(th[len(kn)]))
M = ml.Model(prog)
sos = butter(4, W.BAND, btype="band", fs=ml.FS, output="sos")
def score(key, coefs):
    dr, rs = [], []
    for amp, (y, m, end) in ml.fitted(prog, key, coefs=coefs).items():
        i0, i1 = ml.K0 + int(0.06 * ml.FS), min(end[1], ml.K0 + int(0.6 * ml.FS), len(y) - W.MAXLAG)
        T, L, E = W.lag_track(y[:, 1], m[:, 1], ml.FS, i0, i1)
        t = T / ml.FS * 1000 + 60
        e, l = L[(t >= 60) & (t < 100)], L[(t >= 300) & (t < 600)]
        if len(e) and len(l):
            dr.append(np.median(l) - np.median(e))
        a, b = ml.K0 + int(0.15 * ml.FS), i1
        yb, mb = sosfiltfilt(sos, y[:, 1])[a:b], sosfiltfilt(sos, m[:, 1])[a:b]
        rs.append(10 * np.log10(np.sum((yb - mb) ** 2) / np.sum(yb ** 2)))
    return (float(np.mean(dr)) if dr else np.nan), float(np.median(rs))
for key in keys:
    lf, rt, hf = (float(x) for x in key.split("_"))
    used = ml.line_coefs(M.P, rt, lf, hf)
    raw = ml.line_coefs(M.P, rt, lf, hf, caps=False)
    cur = [(g, r0, curve(pr[2]) if hf < 1.0 else 0.0) for (g, r0, _), pr in zip(used, raw)]
    (d0, r0_), (d1, r1) = score(key, used), score(key, cur)
    print(f"{prog} {key}: drift {d0:+.2f} -> {d1:+.2f} smp | residual 150-600 ms {r0_:.1f} -> {r1:.1f} dB | poles " +
          " ".join(f"{u[2]:.4f}->{c[2]:.4f}" for u, c in zip(used, cur)), flush=True)
