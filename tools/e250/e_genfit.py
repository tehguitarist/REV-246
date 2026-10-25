"""Fit the moving-read generator of e_sim to the random part's statistics, held out by RT step (docs/EMT250_program.md §12).

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_genfit.py [--law tri|walk] [--M 6] [--tsw 10]   (or --sinegrid)

Fitted on the statistics of RT 1.0 / 2.0 / 3.0 / 4.0, scored on 1.5 / 2.5 / 3.6 / 4.5 (and the reverse split). Free: b, lo, hi,
v (up), vdn (down), vsd, share (7 numbers); M, the motion law and the sweep length are fixed per run (the sweep length is the Altiverb capture's,
unknown, so it is scanned). Common random numbers (fixed seeds) keep the cost smooth for Nelder-Mead.
"""
import argparse, json, os, sys
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import e_model as E  # noqa: E402
import e_sim as S  # noqa: E402

NAMES = ("b", "lo", "hi", "v", "vdn", "vsd", "share")
X0 = (0.8, 400.0, 3600.0, 60.0, 80.0, 0.3, 0.5)
SPLITS = (((1.0, 2.0, 3.0, 4.0), (1.5, 2.5, 3.6, 4.5)), ((1.5, 2.5, 3.6, 4.5), (1.0, 2.0, 3.0, 4.0)))


def data_ref(rts):
    S.static_params(); kern = E.kernels([r for r in E.RTS if r >= 1.0]); n = E.K0 + int(2.5 * E.FS)
    return S.stats([E.dev(r)[:n] for r in rts], [E.render(S.X_STATIC, r, kern, n) for r in rts])


def prm_of(x, fixed):
    p = dict(S.DEF) | fixed | dict(zip(NAMES, x))
    p["vsd"] = float(np.clip(p["vsd"], 0, 0.95)); p["share"] = float(np.clip(p["share"], 0, 1))
    return p


def sinegrid():
    """§12.2: the sinusoidal law (the UAD plugin's) on a rate x depth grid, the other parameters at the triangle fit's values;
    cost on split A (fit half) and split B (other half): the choice made on A is held out on B."""
    (tr, te), _ = SPLITS
    rA, rB = data_ref(list(tr)), data_ref(list(te)); base = dict(S.DEF, law="sine", b=0.81, lo=405, hi=3600, share=0.5, M=6, tsw=10.0, seed=7)
    rows = []
    for f in (0.3, 0.55, 0.8, 1.2):
        for A in (8, 15, 25, 40):
            p = base | dict(f=f, A=A)
            cA = S.cost(S.sim_stats(p, list(tr)), rA); cB = S.cost(S.sim_stats(p | dict(seed=11), list(te)), rB)
            rows.append(dict(f=f, A=A, peak_speed=2 * np.pi * f * A, cost_A=cA, cost_B=cB))
            print(f"sine f {f} A {A} (peak {2 * np.pi * f * A:.0f} samples/s): cost A {cA:.1f}, B {cB:.1f}", flush=True)
    json.dump(rows, open(os.path.join(E.OUT, "genfit_sinegrid.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sinegrid", action="store_true")
    ap.add_argument("--law", default="tri"); ap.add_argument("--M", type=int, default=6); ap.add_argument("--tsw", type=float, default=10.0)
    ap.add_argument("--iters", type=int, default=120); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    if a.sinegrid:
        return sinegrid()
    fixed = dict(law=a.law, M=a.M, tsw=a.tsw, seed=7)
    out = dict(fixed=fixed, splits=[])
    for tr, te in SPLITS:
        ref_tr, ref_te = data_ref(list(tr)), data_ref(list(te))
        f = lambda x: S.cost(S.sim_stats(prm_of(x, fixed), list(tr)), ref_tr) if x[1] < x[2] - 200 and x[0] > 0 and x[3] > 0 and x[4] > 0 else 1e9
        res = minimize(f, X0, method="Nelder-Mead", options=dict(maxfev=a.iters, xatol=1e-2, fatol=0.5, initial_simplex=None))
        p = prm_of(res.x, fixed)
        st_te = S.sim_stats(p | dict(seed=11), list(te)); c_te = S.cost(st_te, ref_te)
        # reference: the cost two disjoint halves of the data have against each other
        c_dd = S.cost(ref_tr, ref_te)
        r = dict(train=tr, test=te, params={k: p[k] for k in NAMES}, cost_train=float(res.fun), cost_test=c_te, cost_data_halves=c_dd,
                 sim_test=st_te, data_test=ref_te)
        out["splits"].append(r)
        print(f"{a.law} M {a.M} tsw {a.tsw}: train {tr} cost {res.fun:.1f} -> test {te} cost {c_te:.1f} (data half vs half {c_dd:.1f}) | "
              + " ".join(f"{k} {p[k]:.3g}" for k in NAMES), flush=True)
        S.show(st_te, "   sim  test"); S.show(ref_te, "   data test")
    tag = a.tag or f"{a.law}_M{a.M}_t{a.tsw:g}"
    json.dump(out, open(os.path.join(E.OUT, f"genfit_{tag}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
