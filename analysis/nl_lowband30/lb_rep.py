"""Repeat-capture level per band: captures of one SPAN / EARLY that differ only in LATE (the late starts after the window), same front-shift session.
Per pair, band, channel: gain + fractional shift + clock fitted (as n_low_rep.py), null re capture A; reported raw and with -3 dB (noise of one capture).
usage: lb_rep.py SPAN EARLY[,EARLY...] [pairs_per_early] [t0 t1]"""
import sys, itertools, json
from lb_common import *
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
B = bench_all()
span = float(sys.argv[1]); earlies = [float(x) for x in sys.argv[2].split(",")]
maxp = int(sys.argv[3]) if len(sys.argv) > 3 else 20
t0s, t1s = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (0.05, 0.45)
BANDS = [(20, 100), (100, 1000), (1000, 5000), (5000, 10000), (10000, 12000)]
SOS = {b: bands(*b) for b in BANDS}
rng = np.random.default_rng(0)
out = {}
def fit(A, Bf, a, b, c):
    cs = CubicSpline(np.arange(len(Bf)), Bf[:, c]); tt = np.arange(a, b); Ac = A[a:b, c]
    def cost(p):
        Bs = cs(tt + p[0] + p[1] * 1e-6 * (tt - a)); g = (Ac @ Bs) / (Bs @ Bs)
        return np.sum((Ac - g * Bs) ** 2) / np.sum(Ac ** 2)
    r = minimize(cost, [0.0, 0.0], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-10, "initial_simplex": [[0, 0], [0.05, 0], [0, 2.0]]})
    return r.fun, r.x
for early in earlies:
    keys = sorted(k for k in B if k[0] == span and k[1] == early and not B[k]["trimmed"])
    if len(keys) < 2: print("skip", early); continue
    sig = {k: nl.ir(nl.table()[k]) for k in keys}
    pairs = [(a, b) for a, b in itertools.combinations(keys, 2) if abs(B[a]["front_shift"] - B[b]["front_shift"]) < 0.3 and abs(B[a]["o"] - B[b]["o"]) < 0.02]
    if len(pairs) > maxp: pairs = [pairs[i] for i in sorted(rng.choice(len(pairs), maxp, replace=False))]
    res = {bd: [[], []] for bd in BANDS}
    filt = {}
    for ka, kb in pairs:
        fa, fb = B[ka]["front_shift"], B[kb]["front_shift"]
        n = min(len(sig[ka]), len(sig[kb])) - 200
        a = K0 + int(fa) + int(t0s * FS); b = min(K0 + int(fa) + int(t1s * FS), n - 100)
        sh = int(round(fb - fa))
        for bd in BANDS:
            for kk, shf in ((ka, 0), (kb, sh)):
                if (kk, bd, shf) not in filt:
                    filt[(kk, bd, shf)] = sosfiltfilt(SOS[bd], np.roll(sig[kk], -shf, axis=0)[:n + 100], axis=0) if False else sosfiltfilt(SOS[bd], np.roll(sig[kk], -shf, axis=0)[:n], axis=0)
            Af, Bf = filt[(ka, bd, 0)], filt[(kb, bd, sh)]
            for c in (0, 1):
                f, _ = fit(Af, Bf, a, b, c)
                res[bd][c].append(10 * np.log10(f))
    out[early] = {f"{bd[0]}-{bd[1]}": [float(np.median(res[bd][0])), float(np.median(res[bd][1])), len(res[bd][0])] for bd in BANDS}
    out[str(early) + "_pairs"] = {f"{bd[0]}-{bd[1]}": res[bd] for bd in BANDS}
    print(f"SPAN {span} EARLY {early}: {len(pairs)} pairs of {len(keys)} captures; median pair null (L, R; raw, subtract 3 dB for one capture):", flush=True)
    for bd in BANDS:
        v = out[early][f"{bd[0]}-{bd[1]}"]
        print(f"   {bd[0]:>5}-{bd[1]:<5}: {v[0]:6.1f} {v[1]:6.1f}  (-3 dB: {v[0]-3:6.1f} {v[1]-3:6.1f})", flush=True)
json.dump(out, open(f"rep_{span}_{t0s}_{t1s}.json", "w"))
