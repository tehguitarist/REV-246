"""Q9 / Q10: do S's coefficients sit on a fixed-point grid, and is fixed-point arithmetic visible in the tails?
(model_s_results 36, NEXT_STEPS §3 Q9 / Q10).

    <the S environment> .venv/bin/python tools/loop_read/fixedpoint.py coefs     # seconds -> analysis/lowband/chain3_fixedpoint_coefs.json
    .venv/bin/python tools/loop_read/fixedpoint.py lattice                      # -> analysis/lowband/chain3_fixedpoint_lattice.{log,json}
    .venv/bin/python tools/loop_read/fixedpoint.py consts                       # -> chain3_fixedpoint_consts.log (scalar constants, AP230 gap)
    .venv/bin/python tools/loop_read/fixedpoint.py snap <n> [<offset in steps>]  # -> analysis/loop_cv/hfcells_chain3_cap0.8475_snap<n>[o<off>]/
                                                                                 #    (the baseline's fold set, every hf_poles entry on the 2^-n grid)

<the S environment> = EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds
                      SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_pshoulder/cells.json

coefs: per (setting, loop, take) the unit's loop coefficients as takefit reads them on the baseline (model_s_results 31,
32.9; `_pallf`'s takefit, `_pshoulder`'s for the 28 LF 0.4 / RT >= 0.6 settings; takes holding a GLITCH capture dropped,
DEFECTS were never in): the HF pole p (baseline pole + takefit's shift), the per-pass DC gain D and the Nyquist gain N1 of the
HF 1.0 form (broadband gain x shelf), and the first-order shelf's numerator with the broadband gain folded in, c0 / c1
(D (1 - a) = c0 + c1, N1 (1 + a) = c0 - c1, a the law's shelf pole). Also the baseline's own values (smooth law and
curve), the false-positive check.

lattice: for each coefficient form and grid 2^-n (n = 6..16; 2^-12 is also the TMS32010 MPYK 13-bit immediate's grid),
a precision-weighted Rayleigh statistic Z = |sum w e^{2 pi i v / q}|^2 / sum w^2, w = exp(-2 pi^2 sigma^2 / q^2) the
coherence a value of precision sigma could show on a true grid (sigma from the take-to-take spread, per RT step and form);
analytic p = exp(-Z), and an empirical p against 2000 random spacings (log-uniform within x1.6 of q, powers of two
excluded), which keeps whatever regularity the panel steps put into the values. Z_on = 1 + sum w^2 is what a true grid
would give: where it is below ~5 the grid can't be tested with this precision. Held out: by take (the phase found on
take A tested on take B and vice versa: a reproducible phase, not measurement noise) and by setting (the phase found on
LF 0.1 / 1.0 tested on LF 0.4 / 2.0 and vice versa: a lattice, not a coincidence of the true values). Bonferroni over
all forms x grids.
"""
import glob, json, os, shutil, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools"))
LB = os.path.join(ROOT, "analysis", "lowband")
CV = os.path.join(ROOT, "analysis", "loop_cv")
BASE = os.path.join(CV, "hfcells_chain3_cap0.8475_pshoulder")
TF = os.path.join(LB, "chain3_shelf_takefit_cap0.8475_pallf")
TF_SH = os.path.join(LB, "chain3_shelf_takefit_cap0.8475_pshoulder")
COEFS = os.path.join(LB, "chain3_fixedpoint_coefs.json")
NS = list(range(6, 17))
FORMS = ["p", "D", "N1", "c0", "c1", "ap", "n0", "n1"]   # ap = a p, n0 / n1 = c0 (1 - p), c1 (1 - p): a direct-form loop filter's coefficients
rng = np.random.default_rng(246)


class Tee:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.f.write(s + "\n")
        self.f.flush()


def coefs():
    import lowband_shelf as ls
    from defects import GLITCH
    q2, q1 = ls.q2, ls.q1
    assert os.path.samefile(ls.HF_CELLS_JSON, os.path.join(BASE, "cells.json")), "set SHELF_CELLS to the _pshoulder cells.json"
    law0 = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    a = float(law0["a"])
    info = q2.q1_tables()[0]
    w = np.array([0.0, np.pi])
    rows = []
    for key in sorted(info):
        lab = dict(lf=info[key]["lf"], rt=info[key]["rt"], hf=info[key]["hf"])
        src = TF_SH if os.path.exists(os.path.join(TF_SH, key + ".json")) else TF
        tf = json.load(open(os.path.join(src, key + ".json")))
        bad = {"A" if amp in (0, 90) else "B" for k, amp in GLITCH if k == key}
        fits = {k: np.array(x) for k, x in tf["fit"].items() if k not in bad}
        t, lfdc, kw = ls.base_kw(key, ls.HF_CELLS_JSON)
        law = json.loads(json.dumps(law0))
        law["lf"] = {float(k): x for k, x in law["lf"].items()}
        for k_, x in lfdc.items():
            law["lf"][float(k_)]["dc"] = float(x)
        args = dict(lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], **kw)
        G1 = q2.final_Gs(dict(lab, hf=1.0), law, w, **args)
        hf1 = lab["hf"] >= 1.0
        G = G1 if hf1 else q2.final_Gs(lab, law, w, **args)
        pb = [0.0] * 4 if hf1 else [float((1 - h) / (1 + h)) for h in [np.real(g[-1] / g1[-1]) for g, g1 in zip(G, G1)]]
        for i in range(4):
            d0, n0 = float(np.real(G1[i][0])), float(np.real(G1[i][1]))
            r = dict(key=key, **lab, loop=i, a=a, base=dict(D=d0, N1=n0, p=pb[i] if not hf1 else None), take={})
            for k, u in fits.items():
                D = d0 * 10 ** ((u[i] + u[4 + i]) / 20)
                N1 = n0 * 10 ** (u[i] / 20)
                r["take"][k] = dict(D=D, N1=N1, p=None if hf1 else pb[i] + float(u[8 + i]))
            rows.append(r)
    for r in rows:
        for x in [r["base"]] + list(r["take"].values()):
            x["c0"] = (x["D"] * (1 - a) + x["N1"] * (1 + a)) / 2
            x["c1"] = (x["D"] * (1 - a) - x["N1"] * (1 + a)) / 2
    json.dump(dict(a=a, rows=rows), open(COEFS, "w"), indent=0)
    print(f"{len(rows)} loop rows, {sum(len(r['take']) == 2 for r in rows)} with two takes; shelf pole a = {a:.6f} -> {COEFS}")


def zstat(v, sig, q):
    w = np.exp(-2 * np.pi ** 2 * (sig / q) ** 2)
    s = np.sum(w * np.exp(2j * np.pi * v / q))
    return float(abs(s) ** 2 / np.sum(w ** 2)), float(np.angle(s)), float(1 + np.sum(w ** 2)), w


def aligned(v, sig, q, phi):
    """Held-out phase test: sum w cos(theta - phi) / sqrt(sum w^2 / 2), N(0, 1) under no lattice."""
    w = np.exp(-2 * np.pi ** 2 * (sig / q) ** 2)
    th = 2 * np.pi * v / q
    return float(np.sum(w * np.cos(th - phi)) / np.sqrt(np.sum(w ** 2) / 2))


def rand_p(v, sig, q, z, K=2000):
    qs = q * np.exp(rng.uniform(-np.log(1.6), np.log(1.6), 3 * K))
    l2 = np.log2(qs)
    qs = qs[np.abs(l2 - np.round(l2)) > 0.1][:K]
    zz = np.array([zstat(v, sig, x)[0] for x in qs])
    return float((1 + np.sum(zz >= z)) / (1 + len(zz)))


def is_drop(p):
    return p is not None and min(abs(p - 0.8465), abs(p - 0.847)) < 1e-6


def group_key(r, form):
    """One true coefficient per key: the shelf side (D, N1, c0, c1) is set per (LF, RT, loop), so a cell's HF steps are
    repeat measurements of it; loops on the DC floor (one plateau per RT step, tested as constants) and drop-back poles (one
    value) are left out; the HF pole is per (LF, RT, effective HF, loop): RT 0.3's HF 0.2 and HF 0.4 are one filter (29)."""
    if form in ("p", "ap", "n0", "n1"):
        if r["base"]["p"] is None or is_drop(r["base"]["p"]):
            return None
        return (r["lf"], r["rt"], round(max(r["hf"], 0.12 / r["rt"]), 3), r["loop"])
    if r["base"]["D"] < 0.1905:
        return None
    return (r["lf"], r["rt"], r["loop"])


def table(rows, form, base=False):
    """Two-take values of one form, one per true coefficient (group_key): the group's take A / B means, the mean, and the
    precision of the mean: the take-to-take spread (robust, per RT step) and the spread across the group's settings."""
    g = {}
    for r in rows:
        k = group_key(r, form)
        if k is not None and len(r["take"]) == 2:
            g.setdefault(k, []).append(r)
    keys = sorted(g)
    A = np.array([np.mean([r["take"]["A"][form] for r in g[k]]) for k in keys])
    B = np.array([np.mean([r["take"]["B"][form] for r in g[k]]) for k in keys])
    v = (A + B) / 2
    spread = np.array([np.std([(r["take"]["A"][form] + r["take"]["B"][form]) / 2 for r in g[k]], ddof=1) / np.sqrt(len(g[k]))
                       if len(g[k]) > 1 else 0.0 for k in keys])
    rt = np.array([k[1] for k in keys])
    sig = np.zeros(len(keys))
    for x in set(rt):
        m = rt == x
        sig[m] = 0.741 * np.median(np.abs(A[m] - B[m]))        # sigma of the mean of two takes (Gaussian, robust)
    sig = np.sqrt(sig ** 2 + spread ** 2)
    sel = [dict(lf=k[0], rt=k[1], base=dict(form=np.mean([r["base"][form] for r in g[k]]))) for k in keys]
    for s in sel:
        s["base"][form] = s["base"].pop("form")
    return sel, v, A, B, sig


def lattice():
    from scipy.stats import norm
    d = json.load(open(COEFS))
    rows = d["rows"]
    for r in rows:                                            # the direct-form coefficients of shelf x one-pole
        for x in [r["base"]] + list(r["take"].values()):
            hf = x["p"] is not None
            x["ap"] = d["a"] * x["p"] if hf else None
            x["n0"], x["n1"] = (x["c0"] * (1 - x["p"]), x["c1"] * (1 - x["p"])) if hf else (None, None)
    log = Tee(os.path.join(LB, "chain3_fixedpoint_lattice.log"))
    out = {}
    log(f"Lattice tests of S's per-loop coefficients (model_s_results 36); {len(rows)} loop rows, shelf pole a = {d['a']:.6f}")
    log("Z: precision-weighted Rayleigh (any offset); p_an = exp(-Z); p_rand: random spacings within x1.6; Z_on: expected Z on a true grid;")
    log("phi: the phase (0 = grid points at k 2^-n); held out: by take (phase of A tested on B, B on A) and by setting (LF 0.1/1.0 <-> 0.4/2.0), z-scores")
    ntest = len(FORMS) * len(NS)
    for form in FORMS:
        sel, v, A, B, sig = table(rows, form)
        lf = np.array([r["lf"] for r in sel])
        base = np.array([r["base"][form] for r in sel])
        log(f"\n{form}: n = {len(v)} values ({v.min():.4f} .. {v.max():.4f}); precision of the mean of two takes, median {np.median(sig):.2e} "
            f"(per RT step {', '.join(f'{x}: {np.median(sig[np.array([r['rt'] for r in sel]) == x]):.1e}' for x in sorted({r['rt'] for r in sel}))})")
        log(f"   {'grid':>6s} {'q':>9s} {'q/sig':>6s} {'Z_on':>6s} {'Z':>6s} {'p_an':>7s} {'p_rand':>7s} {'phi':>6s} | {'take z':>12s} | {'setting z':>12s} | baseline Z")
        out[form] = []
        sA = sig * np.sqrt(2)                                     # single take
        for n in NS:
            q = 2.0 ** -n
            z, phi, zon, _ = zstat(v, sig, q)
            if zon < 2:                                           # sum w^2 < 1: a true grid could not show; not a test
                log(f"   2^-{n:<3d} {q:9.2e} {q / np.median(sig):6.1f} {zon:6.1f}  (not testable)")
                continue
            pr = rand_p(v, sig, q, z)
            phA, phB = zstat(A, sA, q)[1], zstat(B, sA, q)[1]
            tz = (aligned(B, sA, q, phA), aligned(A, sA, q, phB))
            m1 = np.isin(lf, [0.1, 1.0])
            sz = (aligned(v[~m1], sig[~m1], q, zstat(v[m1], sig[m1], q)[1]), aligned(v[m1], sig[m1], q, zstat(v[~m1], sig[~m1], q)[1]))
            zb = zstat(base, np.full(len(base), 1e-9), q)[0]
            log(f"   2^-{n:<3d} {q:9.2e} {q / np.median(sig):6.1f} {zon:6.1f} {z:6.2f} {np.exp(-z):7.3f} {pr:7.3f} {phi:+6.2f} | {tz[0]:+5.2f} {tz[1]:+5.2f} | {sz[0]:+5.2f} {sz[1]:+5.2f} | {zb:6.2f}")
            out[form].append(dict(n=n, q=q, z=z, p_an=float(np.exp(-z)), p_rand=pr, phi=phi, z_on=zon, take_z=tz, setting_z=sz,
                                  setting_p=[float(norm.sf(x)) for x in sz], base_z=zb))
    ntest = sum(len(out[f]) for f in out)
    best = min((x["p_an"], f, x["n"]) for f in out for x in out[f])
    log(f"\nSmallest analytic p over {ntest} testable form x grid tests: {best[0]:.3g} ({best[1]}, 2^-{best[2]}); Bonferroni {min(1, best[0] * ntest):.3g}")
    log("Upper limit on an on-grid share f (|mean coherence| / expected, 95 %): " + "; ".join(
        f"{f}: " + ", ".join(f"2^-{x['n']} {min(1, (np.sqrt(x['z']) + 1.73) / np.sqrt(x['z_on'] - 1)):.2f}" for x in out[f] if x["z_on"] >= 5) for f in out))
    testable = {f: [x["n"] for x in out[f] if x["z_on"] >= 5] for f in out}
    log("Grids a true lattice would show at Z_on >= 5 (testable): " + "; ".join(f"{f}: 2^-{min(v)}..2^-{max(v)}" if v else f"{f}: none" for f, v in testable.items()))
    json.dump(out, open(os.path.join(LB, "chain3_fixedpoint_lattice.json"), "w"), indent=1)


def snap(n, off=0.0):
    """A variant fold set: the baseline's, every cells.json hf_poles entry snapped to the 2^-n grid (shifted by off steps)."""
    q = 2.0 ** -n
    tag = f"snap{n}" + (f"o{off:g}" if off else "")
    out = os.path.join(CV, f"hfcells_chain3_cap0.8475_{tag}")
    os.makedirs(out, exist_ok=True)
    for f in glob.glob(os.path.join(BASE, "*.json")):
        if not f.endswith("_rows.json"):
            shutil.copy(f, out)
    cells = json.load(open(os.path.join(BASE, "cells.json")))
    dev = []
    for k, c in cells.items():
        if "hf_poles" in c:
            new = [float((np.round(p / q - off) + off) * q) for p in c["hf_poles"]]
            dev += [x - p for x, p in zip(new, c["hf_poles"])]
            c["hf_poles"] = new
    json.dump(cells, open(os.path.join(out, "cells.json"), "w"), indent=1)
    dev = np.array(dev)
    print(f"{tag}: {len(dev)} poles moved by rms {np.sqrt(np.mean(dev ** 2)):.2e}, max {np.abs(dev).max():.2e} -> {out} (loop_cv --retake --tag {tag})")


if __name__ == "__main__":
    st = sys.argv[1]
    if st == "coefs":
        coefs()
    elif st == "lattice":
        lattice()
    elif st == "snap":
        snap(int(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 0.0)
