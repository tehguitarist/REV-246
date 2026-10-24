"""M / L's own best HF pole at strong HF, against S's pole curve (model_s_results §25.17). Per setting: the law's poles
per line (uncapped p_law, capped as used), a common pole shift dp on a grid applied to every line, and for each dp the
capture-vs-model lag drift (warp.py's track, R, median lag 300-600 ms minus 60-100 ms, mean over captures); the dp
that zeroes the drift (linear interpolation) -> p_best = p_used + dp, compared with S's curve at the same p_law.

    .venv/bin/python tools/loop_read/ml_poles.py M|L key,key,..."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
os.environ.setdefault("ML_POLE_CURVE", "0")         # these compare against the capped law (§25.16-25.17)
import ml_lines as ml
import warp as W
prog, keys = sys.argv[1], sys.argv[2].split(",")
pc = json.load(open(os.path.join(ROOT, "analysis/lowband/chain3_polecurve.json")))
kn, th = np.array(pc["knots"]), np.array(pc["all"])
curve = lambda x: min(x + float(np.interp(x, kn, th[:len(kn)])), float(th[len(kn)]))
DP = np.round(np.arange(-0.060, 0.0201, 0.005), 4)
M = ml.Model(prog)


def drift(key, coefs):
    out = []
    for amp, (y, m, end) in ml.fitted(prog, key, coefs=coefs).items():
        i0, i1 = ml.K0 + int(0.06 * ml.FS), min(end[1], ml.K0 + int(0.6 * ml.FS), len(y) - W.MAXLAG)
        T, L, E = W.lag_track(y[:, 1], m[:, 1], ml.FS, i0, i1)
        t = T / ml.FS * 1000 + 60
        e, l = L[(t >= 60) & (t < 100)], L[(t >= 300) & (t < 600)]
        if len(e) and len(l):
            out.append(np.median(l) - np.median(e))
    return float(np.mean(out)) if out else np.nan


rows = []
for key in keys:
    lf, rt, hf = (float(x) for x in key.split("_"))
    used = ml.line_coefs(M.P, rt, lf, hf)
    raw = ml.line_coefs(M.P, rt, lf, hf, caps=False)
    p_law = [c[2] for c in raw]; p_used = [c[2] for c in used]
    d = [drift(key, [(g, r0, p + dp) for g, r0, p in used]) for dp in DP]
    dd = np.array(d)
    ok = np.isfinite(dd)
    s = np.nan
    if ok.sum() >= 2:
        j = np.where(np.diff(np.sign(dd[ok])) != 0)[0]
        if len(j):
            x0, x1, y0, y1 = DP[ok][j[0]], DP[ok][j[0] + 1], dd[ok][j[0]], dd[ok][j[0] + 1]
            s = float(x0 - y0 * (x1 - x0) / (y1 - y0))
    pb = [p + s for p in p_used] if np.isfinite(s) else [np.nan] * len(p_used)
    pcur = [curve(x) for x in p_law]
    rows.append(dict(key=key, p_law=p_law, p_used=p_used, dp_best=s, p_best=pb, p_curve=pcur, drift=d))
    print(f"{prog} {key}: lines P {list(M.P)}; p_law " + " ".join(f"{x:.4f}" for x in p_law) + " | used " + " ".join(f"{x:.4f}" for x in p_used) +
          f" | drift by dp {' '.join(f'{x:+.2f}' for x in d)} -> dp {s:+.4f} | best " + " ".join(f"{x:.4f}" for x in pb) +
          " | S curve " + " ".join(f"{x:.4f}" for x in pcur), flush=True)
json.dump(rows, open(os.path.join(ROOT, f"analysis/lowband/chain3_ml_poles_{prog}.json"), "w"), indent=0)
