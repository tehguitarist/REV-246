"""Model null per band in the repeat window (z+0.05 .. z+0.45 s), SPAN given, bench alignment of bench_bk_a0.65 (env: NL_BK=1 NL_BK_ALPHA=0.65 LB_BENCH=bench_bk_a0.65).
fixed = bench o / clock / gain; fit = gain + fractional shift + clock re-fitted in the band (the repeat method's freedom).  usage: lb_model_band.py SPAN [late_stride]"""
import sys, json
from lb_common import *
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
B = bench_all()
span = float(sys.argv[1]); stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1
BANDS = [(20, 100), (100, 1000), (1000, 5000), (5000, 10000), (10000, 12000)]
SOS = {b: bands(*b) for b in BANDS}
t0s, t1s = 0.05, 0.45
keys = sorted(k for k in B if k[0] == span and not B[k]["trimmed"] and k[1] > 0)
keys = [k for i, k in enumerate(keys) if (i % stride == 0)]
out = {}
for k in keys:
    r = B[k]; n = K0 + int(abs(r["front_shift"])) + int(0.5 * FS)
    y, m, _ = load(k, r, n=n, band=False)           # y: unwarped by the bench clock, m: bench alignment and gain
    raw = nl.ir(nl.table()[k])[:n + 64]
    z = K0 + r["front_shift"]; a = int(z + t0s * FS); b = int(z + t1s * FS)
    res = {}
    for bd in BANDS:
        yb, mb = sosfiltfilt(SOS[bd], y, axis=0), sosfiltfilt(SOS[bd], m, axis=0)
        mu = M.render(k, n, r["o"], r["front_shift"], band=False)          # unit-gain model for the free fit
        mub = sosfiltfilt(SOS[bd], mu, axis=0)
        rawb = sosfiltfilt(SOS[bd], raw, axis=0)
        row = []
        for c in (0, 1):
            fx = 10 * np.log10(np.sum((yb[a:b, c] - mb[a:b, c]) ** 2) / np.sum(yb[a:b, c] ** 2))
            cs = CubicSpline(np.arange(len(mub)), mub[:, c]); tt = np.arange(a, b); Ac = rawb[a:b, c]
            def cost(p):
                Ms = cs(tt + p[0] + p[1] * 1e-6 * (tt - a)); g = (Ac @ Ms) / (Ms @ Ms)
                return np.sum((Ac - g * Ms) ** 2) / np.sum(Ac ** 2)
            f = minimize(cost, [0.0, r["ppm"] * 0 + 0.0], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-10, "initial_simplex": [[0, 0], [0.05, 0], [0, 2.0]]}).fun
            row.append((float(fx), float(10 * np.log10(f))))
        res[f"{bd[0]}-{bd[1]}"] = row
    out["_".join(map(str, k))] = res
    print(k, {b: [round(v[0], 1) for v in res[b]] for b in ("20-100", "100-1000")}, flush=True)
json.dump(out, open(f"modelband_{span}.json", "w"))
