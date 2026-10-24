"""The per-loop HF poles (lowband_shelf.py --stage hfloop, §25.18 base) read out (model_s_results §25.19): precision from
the two takes (AMPL 0 + 90 against 30 + 60), the fit's gain, the poles against the law (p_law) with a fresh smooth curve
(piecewise linear, leave-one-RT-out), a periodic-remainder scan (staircase) and a grid test on the absolute poles
(Rayleigh, 2^-6 .. 2^-12), both weighted by precision and on the well-determined poles only; permutation p-values.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/loop_read/hfloop_read.py"""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.optimize import least_squares
import dense_mid as dm
q1, ls = dm.q1, dm.ls
law = dm.q2.final_injection(compact=ls.lb.lf_.COMPACT)[1]
law["lf"] = {float(k): x for k, x in law["lf"].items()}
rows = []
for p in sorted(glob.glob(os.path.join(ROOT, "analysis/lowband/chain3_shelf_hfloop_cap0.8475/*.json"))):
    d = json.load(open(p))
    L = law["lf"][d["lf"]]
    for i, Pi in enumerate(q1.P):
        T = float(10 ** (-3 * Pi / q1.FS / d["rt"] * L["k"] * (1 / d["hf"] - 1)))
        rows.append(dict(key=d["key"], lf=d["lf"], rt=d["rt"], hf=d["hf"], loop=i, p_law=float(q1.onepole_for(np.array(T), law["f_ref"])),
                         p_base=d["p_base"][i], p=d["all"]["p"][i], pA=d["take_0_90"]["p"][i], pB=d["take_30_60"]["p"][i],
                         gain=10 * np.log10(d["all"]["cost"] / d["all"]["base"])))
x = np.array([r["p_law"] for r in rows]); p = np.array([r["p"] for r in rows]); dAB = np.array([r["pA"] - r["pB"] for r in rows])
rt = np.array([r["rt"] for r in rows]); pb = np.array([r["p_base"] for r in rows])
sig = np.abs(dAB) / 2 + 1e-5                                         # each take's estimate ~ |pA - pB| / sqrt 2; the mean ~ / 2
print(f"{len(rows)} loop poles from {len(rows) // 4} settings")
print("precision (|take A - take B|, median) by p_law: " + "; ".join(f"{a}-{b}: {1e4 * np.median(np.abs(dAB[(x >= a) & (x < b)])):.1f}e-4 (n {np.sum((x >= a) & (x < b))})"
      for a, b in ((0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.75), (0.75, 1.0))))
g = np.array([r["gain"] for r in rows[::4]])
print(f"fit gain over the §25.18 baseline (100 Hz-10 kHz from 50 ms, all captures, in sample): median {np.median(g):+.2f} dB, 10 % {np.percentile(g, 10):+.2f}")
print(f"p - p_base: median {1e3 * np.median(p - pb):+.2f}e-3, rms {1e3 * np.sqrt(np.mean((p - pb) ** 2)):.2f}e-3")
cap = q1.ms and float(os.environ.get("EMT_PMAX", 0.84))
hi = x > 0.86
print(f"poles where the law asks past the ceiling (p_law > 0.86, n {hi.sum()}): median {np.median(p[hi]):.4f}, IQR {np.percentile(p[hi], 25):.4f}-{np.percentile(p[hi], 75):.4f} (ceiling {cap})")
KN = np.array([0.0, 0.15, 0.3, 0.45, 0.55, 0.62, 0.7, 0.78, 0.9, 1.0])
def curve(th, xx):
    return np.minimum(xx + np.interp(xx, KN, th[:-1]), th[-1])
good = (sig < 3e-4) & (x < 0.86)
w = 1 / sig ** 2
pred = np.zeros_like(p)
for r_ in sorted(set(rt)):
    tr, te = good & (rt != r_), rt == r_
    f = least_squares(lambda th: (curve(th, x[tr]) - p[tr]) * np.sqrt(w[tr]) / np.sqrt(w[tr]).mean(), np.r_[np.zeros(len(KN)), cap])
    pred[te] = curve(f.x, x[te])
res = p - pred
print(f"well-determined poles (|A-B|/2 < 3e-4, p_law < 0.86): n {good.sum()}; held-out curve remainder rms {1e3 * np.sqrt(np.average(res[good] ** 2, weights=w[good])):.2f}e-3 "
      f"(take precision rms {1e3 * np.sqrt(np.mean(sig[good] ** 2)):.2f}e-3)")
rng = np.random.default_rng(1)
def scan(v, u, periods, wt):
    return np.array([abs(np.sum(wt * v * np.exp(2j * np.pi * u / P)) / np.sum(wt)) for P in periods])
loss = -20 * np.log10(np.abs((1 - x) / (1 - x * np.exp(-2j * np.pi * 6400 / q1.FS))))
for nm, u, periods in (("p_law", x, np.geomspace(0.001, 0.05, 500)), ("target loss dB", loss, np.geomspace(0.02, 2.0, 500))):
    v, uu, wt = res[good] - np.average(res[good], weights=w[good]), u[good], w[good] / w[good].mean()
    pw = scan(v, uu, periods, wt); j = int(np.argmax(pw))
    null = np.array([scan(rng.permutation(v), uu, periods, wt).max() for _ in range(300)])
    print(f"staircase vs {nm}: best period {periods[j]:.4g}, amplitude {2e3 * pw[j]:.3f}e-3, permutation p = {np.mean(null >= pw[j]):.3f} (null 95 % {2e3 * np.quantile(null, 0.95):.3f}e-3)")
print("grid test on the absolute poles (well-determined), Rayleigh R and permutation-free p ~ exp(-n_eff R^2):")
pg, wg = p[good], w[good] / w[good].mean()
neff = wg.sum() ** 2 / np.sum(wg ** 2)
for n in range(6, 13):
    R = abs(np.sum(wg * np.exp(2j * np.pi * pg * 2 ** n)) / wg.sum())
    print(f"   2^-{n} ({2.0 ** -n:.5f}): R {R:.3f}, p ~ {np.exp(-neff * R ** 2):.3g}  (n_eff {neff:.0f}; resolvable if take precision << step: {1e3 * np.median(sig[good]):.2f}e-3)")
json.dump(dict(rows=rows), open(os.path.join(ROOT, "analysis/lowband/chain3_hfloop_read.json"), "w"))
print("\nunweighted grid test (well-determined poles, n %d):" % good.sum())
for n in range(6, 13):
    R = abs(np.mean(np.exp(2j * np.pi * p[good] * 2 ** n)))
    print(f"   2^-{n}: R {R:.3f}, p ~ {np.exp(-good.sum() * R ** 2):.3g}")
keys = np.array([r["key"] for r in rows])
com, own = [], []
for k in sorted(set(keys[good])):
    s = good & (keys == k)
    if s.sum() == 4:
        m = res[s].mean(); com.append(m); own += list(res[s] - m)
com, own = np.array(com), np.array(own)
print(f"\ncurve remainder split (settings with 4 well-determined loops, n {len(com)}): common to the setting's loops rms {1e3 * com.std():.2f}e-3; "
      f"loop-specific rms {1e3 * own.std() * np.sqrt(4 / 3):.2f}e-3 (take precision ~0.1e-3)")
lo = [r for r in rows]
byc = {}
for k_, cm in zip(sorted(set(keys[good])), com):
    lf, rt_, hf = k_.split("_"); byc.setdefault((hf, rt_), []).append(cm)
spread = [np.std(v) for v in byc.values() if len(v) >= 3]
print(f"common part by (HF, RT) cell: between-cell rms {1e3 * np.std([np.mean(v) for v in byc.values()]):.2f}e-3, within-cell (across LF) rms {1e3 * np.sqrt(np.mean(np.square(spread))):.2f}e-3")
