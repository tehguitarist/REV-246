"""Q9 / Q10: do S's coefficients sit on a fixed-point grid, and is fixed-point arithmetic visible in the tails?
(model_s_results 36, NEXT_STEPS §3 Q9 / Q10).

    <the S environment> .venv/bin/python tools/loop_read/fixedpoint.py coefs [<takefit dir> <out json>]   # seconds -> chain3_fixedpoint_coefs.json
    .venv/bin/python tools/loop_read/fixedpoint.py lattice      # -> analysis/lowband/chain3_fixedpoint_lattice.{log,json}
    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/fixedpoint.py consts   # -> chain3_fixedpoint_consts.{log,json} (scalar constants, AP230 gap and steps)
    .venv/bin/python tools/loop_read/fixedpoint.py panel        # -> chain3_fixedpoint_panel.{log,json} (a panel-step table or a computation?)
    .venv/bin/python tools/loop_read/fixedpoint.py snap <n> [<offset in steps>]  # -> analysis/loop_cv/hfcells_chain3_cap0.8475_snap<n>[o<off>]/
                                                                 #    (the baseline's fold set, every hf_poles entry on the 2^-n grid)
    .venv/bin/python tools/loop_read/fixedpoint.py q10 [keys]   # -> chain3_fixedpoint_q10.{log,json} (float vs 16-bit loops, c/s_loops_fx.c,
                                                                 #    compiled with cc into $TMPDIR)

<the S environment> = EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds
                      SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_pshoulder/cells.json

coefs: per (setting, loop, take) the unit's loop coefficients as takefit reads them on the baseline (model_s_results 31,
32.9; `_pallf`'s takefit, `_pshoulder`'s for the 28 LF 0.4 / RT >= 0.6 settings; takes holding a GLITCH capture dropped,
DEFECTS were never in): the HF pole p (baseline pole + takefit's shift), the per-pass DC gain D and the Nyquist gain N1 of the
HF 1.0 form (broadband gain x shelf), and the first-order shelf's numerator with the broadband gain folded in, c0 / c1
(D (1 - a) = c0 + c1, N1 (1 + a) = c0 - c1, a the law's shelf pole). Also the baseline's own values (smooth law and
curve), the false-positive check, and the law's D / N1 without the shelf cells (`nocell`, for panel). With a takefit dir and an
output (and SHELF_CELLS = that baseline's cells.json) it reads another baseline's takefit: the read-bias check (_lowrt).

lattice: for each coefficient form and grid 2^-n (n = 6..16; 2^-12 is also the TMS32010 MPYK 13-bit immediate's grid),
a precision-weighted Rayleigh statistic Z = |sum w e^{2 pi i v / q}|^2 / sum w^2, w = exp(-2 pi^2 sigma^2 / q^2) the
coherence a value of precision sigma could show on a true grid. sigma = the take-to-take spread (per (RT, effective HF) cell
for the pole forms, per RT step for the shelf) + the read bias both takes share (the sd of this read minus the _lowrt
baseline's read of the same takes, per RT step; a lower bound). p_an = exp(-Z); p_rand against 2000 random spacings
(log-uniform within x1.6 of q, powers of two excluded), the reported p. Z_on = 1 + sum w^2 is what a true grid would give:
where it is below ~5 the grid can't be tested. Held out by setting (the phase found on LF 0.1 / 1.0 tested on LF 0.4 / 2.0 and
vice versa); the take A -> B phase z is ~ sqrt(2 Z) whenever q >> the take noise, so it only shows reproducibility. Bonferroni
over all forms x grids.
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
FORMS = ["p", "D", "N1", "c0", "c1", "ap", "n0", "n1", "gp", "g1p"]   # ap = a p, n0 / n1 = c0 (1 - p), c1 (1 - p): a direct-form loop filter's coefficients; gp / g1p = N1 p, N1 (1 - p)
rng = np.random.default_rng(246)


class Tee:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.f.write(s + "\n")
        self.f.flush()


def coefs(tf_dir=None, out=COEFS):
    """tf_dir / out: another baseline's takefit (with SHELF_CELLS = that baseline's cells.json), e.g. the bias check's _lowrt."""
    import lowband_shelf as ls
    from defects import GLITCH
    q2, q1 = ls.q2, ls.q1
    if tf_dir is None:
        assert os.path.samefile(ls.HF_CELLS_JSON, os.path.join(BASE, "cells.json")), "set SHELF_CELLS to the _pshoulder cells.json"
    law0 = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    a = float(law0["a"])
    info = q2.q1_tables()[0]
    w = np.array([0.0, np.pi])
    rows = []
    for key in sorted(info):
        lab = dict(lf=info[key]["lf"], rt=info[key]["rt"], hf=info[key]["hf"])
        src = tf_dir or (TF_SH if os.path.exists(os.path.join(TF_SH, key + ".json")) else TF)
        if not os.path.exists(os.path.join(src, key + ".json")):
            continue
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
        GL = q2.final_Gs(dict(lab, hf=1.0), law, w, **{k_: x for k_, x in args.items() if k_ != "loopcell"})   # without the shelf cells
        for i in range(4):
            d0, n0 = float(np.real(G1[i][0])), float(np.real(G1[i][1]))
            r = dict(key=key, **lab, loop=i, a=a, base=dict(D=d0, N1=n0, p=pb[i] if not hf1 else None), take={},
                     nocell=dict(D=float(np.real(GL[i][0])), N1=float(np.real(GL[i][1]))))
            for k, u in fits.items():
                D = d0 * 10 ** ((u[i] + u[4 + i]) / 20)
                N1 = n0 * 10 ** (u[i] / 20)
                r["take"][k] = dict(D=D, N1=N1, p=None if hf1 else pb[i] + float(u[8 + i]))
            rows.append(r)
    for r in rows:
        for x in [r["base"]] + list(r["take"].values()):
            x["c0"] = (x["D"] * (1 - a) + x["N1"] * (1 + a)) / 2
            x["c1"] = (x["D"] * (1 - a) - x["N1"] * (1 + a)) / 2
    json.dump(dict(a=a, rows=rows), open(out, "w"), indent=0)
    print(f"{len(rows)} loop rows, {sum(len(r['take']) == 2 for r in rows)} with two takes; shelf pole a = {a:.6f} -> {out}")


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


POLE_FORMS = ("p", "ap", "n0", "n1", "gp", "g1p")


def pole_limited(r):
    """HF 0.2 at RT <= 0.8 and HF 0.4 at RT 0.3: their shelf reads carry the pole's error (32.3), not shelf repeats."""
    return r["rt"] <= 0.8 and (r["hf"] <= 0.2 or (r["rt"] <= 0.3 and r["hf"] <= 0.4))


def group_key(r, form):
    """One true coefficient per key: the shelf side (D, N1, c0, c1) is set per (LF, RT, loop), so a cell's HF steps are
    repeat measurements of it (the pole-limited HF steps left out, as the shelf cells do, 32.3); loops on the DC floor (one
    plateau per RT step, tested as constants) and drop-back poles (one value) are left out; the HF pole is per (LF, RT,
    effective HF, loop): RT 0.3's HF 0.2 and HF 0.4 are one filter (29)."""
    if form in POLE_FORMS:
        if r["base"]["p"] is None or is_drop(r["base"]["p"]):
            return None
        return (r["lf"], r["rt"], round(max(r["hf"], 0.12 / r["rt"]), 3), r["loop"])
    if r["base"]["D"] < 0.1905 or pole_limited(r):
        return None
    return (r["lf"], r["rt"], r["loop"])


def add_forms(d):
    """The products a TMS32010 multiplier would see: the direct-form loop filter (shelf x one-pole) a p, c0 (1 - p), c1 (1 - p)
    (a + p is p's lattice shifted by a constant: the same test), and the one-pole's two multipliers with the loop's HF gain
    N1 folded in, N1 p and N1 (1 - p) (1 - p is p's lattice)."""
    for r in d["rows"]:
        for x in [r["base"]] + list(r["take"].values()):
            hf = x["p"] is not None
            x["ap"] = d["a"] * x["p"] if hf else None
            x["n0"], x["n1"] = (x["c0"] * (1 - x["p"]), x["c1"] * (1 - x["p"])) if hf else (None, None)
            x["gp"], x["g1p"] = (x["N1"] * x["p"], x["N1"] * (1 - x["p"])) if hf else (None, None)


def table(rows, form, bias=None):
    """Two-take values of one form, one per true coefficient (group_key): the group's take A / B means, the mean, and the
    precision of the mean: the take-to-take spread (robust; per (RT, effective HF) cell for the pole forms, whose precision
    varies tenfold with the pole, per RT step for the shelf), the spread across the group's settings, and (bias: a dict
    of per-RT-step sd) the model-dependence of the read, which both takes share and |A - B| cannot show."""
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
    cell = [(k[1], k[2]) if form in POLE_FORMS else k[1] for k in keys]
    sig = np.zeros(len(keys))
    for x in set(cell):
        m = np.array([c == x for c in cell])
        sig[m] = 0.741 * np.median(np.abs(A[m] - B[m]))        # sigma of the mean of two takes (Gaussian, robust)
    sig = np.sqrt(sig ** 2 + spread ** 2)
    sb = np.array([bias.get(k[1], 0.0) for k in keys]) if bias else np.zeros(len(keys))
    sel = [dict(lf=k[0], rt=k[1], base={form: np.mean([r["base"][form] for r in g[k]])}) for k in keys]
    return sel, v, A, B, sig, sb


def read_bias(d, alt_path):
    """Per form and RT step: the robust sd of (this read - another baseline's read of the same captures), take means. Both
    reads fit the same takes, so their take noise is shared and cancels: what is left is how much the read moves when the
    rest of the model changes, a lower bound on the bias common to both takes."""
    if not os.path.exists(alt_path):
        return {}
    alt = json.load(open(alt_path))
    add_forms(alt)
    o = {(r["key"], r["loop"]): r for r in alt["rows"]}
    out = {}
    for form in FORMS:
        dd = {}
        for r in d["rows"]:
            r2 = o.get((r["key"], r["loop"]))
            if group_key(r, form) is None or len(r["take"]) < 2 or r2 is None or len(r2["take"]) < 2:
                continue
            x = (r["take"]["A"][form] + r["take"]["B"][form] - r2["take"]["A"][form] - r2["take"]["B"][form]) / 2
            dd.setdefault(r["rt"], []).append(x)
        out[form] = {rt: float(1.4826 * np.median(np.abs(np.array(x) - np.median(x)))) for rt, x in dd.items()}
    return out


def lattice():
    from scipy.stats import norm
    d = json.load(open(COEFS))
    rows = d["rows"]
    add_forms(d)
    bias = read_bias(d, COEFS.replace(".json", "_lowrt.json"))
    log = Tee(os.path.join(LB, "chain3_fixedpoint_lattice.log"))
    out = {}
    log(f"Lattice tests of S's per-loop coefficients (model_s_results 36); {len(rows)} loop rows, shelf pole a = {d['a']:.6f}")
    log("Z: precision-weighted Rayleigh (any offset), weights from take noise + read bias; p_an = exp(-Z); p_rand: random spacings within x1.6")
    log("(powers of two excluded), the look-elsewhere-safe p; Z_on: the Z a true grid would give with take noise only / with the read bias too")
    log("(bias: the sd of this read minus the _lowrt baseline's read of the same takes, per RT step; a lower bound); b/q: that bias over the step.")
    log("phi: the phase (0 = grid points at k 2^-n). Held out by setting: the phase found on LF 0.1 / 1.0 tested on LF 0.4 / 2.0 and vice versa, z.")
    log("A / B phase z (take A's phase on take B): ~ sqrt(2 Z) whenever q >> take noise, so it shows only that the values are reproducible, not a lattice.")
    for form in FORMS:
        sel, v, A, B, sig, sb = table(rows, form, bias.get(form))
        st = np.sqrt(sig ** 2 + sb ** 2)
        lf = np.array([r["lf"] for r in sel])
        rts = np.array([r["rt"] for r in sel])
        base = np.array([r["base"][form] for r in sel])
        log(f"\n{form}: n = {len(v)} values ({v.min():.4f} .. {v.max():.4f}); take precision of the mean, median {np.median(sig):.2e}; read bias median {np.median(sb):.2e} "
            f"(per RT step take / bias: {', '.join(f'{x}: {np.median(sig[rts == x]):.1e} / {np.median(sb[rts == x]):.1e}' for x in sorted(set(rts)))})")
        log(f"   {'grid':>6s} {'q':>9s} {'b/q':>5s} {'Z_on t':>6s} {'Z_on':>6s} {'Z':>6s} {'p_an':>7s} {'p_rand':>7s} {'phi':>6s} | {'setting z':>12s} | {'A/B z':>11s} | baseline Z")
        out[form] = []
        sA = sig * np.sqrt(2)                                     # single take
        for n in NS:
            q = 2.0 ** -n
            zon_t = zstat(v, sig, q)[2]
            z, phi, zon, _ = zstat(v, st, q)
            if zon < 2:                                           # sum w^2 < 1: a true grid could not show; not a test
                log(f"   2^-{n:<3d} {q:9.2e} {np.median(sb) / q:5.2f} {zon_t:6.1f} {zon:6.1f}  (not testable)")
                continue
            pr = rand_p(v, st, q, z)
            phA, phB = zstat(A, sA, q)[1], zstat(B, sA, q)[1]
            tz = (aligned(B, sA, q, phA), aligned(A, sA, q, phB))
            m1 = np.isin(lf, [0.1, 1.0])
            sz = (aligned(v[~m1], st[~m1], q, zstat(v[m1], st[m1], q)[1]), aligned(v[m1], st[m1], q, zstat(v[~m1], st[~m1], q)[1]))
            zb = zstat(base, np.full(len(base), 1e-9), q)[0]
            log(f"   2^-{n:<3d} {q:9.2e} {np.median(sb) / q:5.2f} {zon_t:6.1f} {zon:6.1f} {z:6.2f} {np.exp(-z):7.3f} {pr:7.3f} {phi:+6.2f} | {sz[0]:+5.2f} {sz[1]:+5.2f} | {tz[0]:+5.2f} {tz[1]:+5.2f} | {zb:6.2f}")
            out[form].append(dict(n=n, q=q, z=z, p_an=float(np.exp(-z)), p_rand=pr, phi=phi, z_on=zon, z_on_take=zon_t, bias_q=float(np.median(sb) / q),
                                  take_z=tz, setting_z=sz, setting_p=[float(norm.sf(x)) for x in sz], base_z=zb))
    ntest = sum(len(out[f]) for f in out)
    best = min((x["p_rand"], f, x["n"], x["p_an"]) for f in out for x in out[f])
    log(f"\nSmallest p_rand over {ntest} testable form x grid tests (with the bias): {best[0]:.3g} ({best[1]}, 2^-{best[2]}, p_an {best[3]:.3g}); "
        f"Bonferroni {min(1, best[0] * ntest):.3g} (conservative: the pole forms are functions of one another)")
    log("Upper limit on an on-grid share f (|mean coherence| / expected, 95 %), Z_on >= 5: " + "; ".join(
        f"{f}: " + ", ".join(f"2^-{x['n']} {min(1, (np.sqrt(x['z']) + 1.73) / np.sqrt(x['z_on'] - 1)):.2f}" for x in out[f] if x["z_on"] >= 5) for f in out))
    testable = {f: [x["n"] for x in out[f] if x["z_on"] >= 5] for f in out}
    log("Grids a true lattice would show at Z_on >= 5 (testable, bias included): " + "; ".join(f"{f}: 2^-{min(v)}..2^-{max(v)}" if v else f"{f}: none" for f, v in testable.items()))
    testable_t = {f: [x["n"] for x in out[f] if x["z_on_take"] >= 5] for f in out}
    log("... with take noise only (the first read's claim): " + "; ".join(f"{f}: 2^-{min(v)}..2^-{max(v)}" if v else f"{f}: none" for f, v in testable_t.items()))
    json.dump(dict(bias=bias, tests=out), open(os.path.join(LB, "chain3_fixedpoint_lattice.json"), "w"), indent=1)


CSD2 = np.unique(np.abs([s1 * 2.0 ** -a + s2 * 2.0 ** -b for a in range(0, 13) for b in range(0, 13) for s1 in (1,) for s2 in (-1, 0, 1)]))


def near(c, pts):
    j = np.argmin(np.abs(pts - c))
    return pts[j], c - pts[j]


def chance(c, d, pts, sig, K=20000):
    """How often a constant drawn uniformly around c (+-max(0.05 |c|, 20 sig)) lands at least as close to the set pts as c does."""
    h = max(0.05 * abs(c), 20 * sig)
    u = c + rng.uniform(-h, h, K)
    m = (pts > c - 2 * h) & (pts < c + 2 * h)
    dd = np.min(np.abs(u[:, None] - pts[None, m if m.sum() >= 2 else slice(None)]), axis=1)
    return float((1 + np.sum(dd <= abs(d))) / (1 + K))


def rationals(c, qmax=12):
    lo, hi = np.floor(abs(c)) - 1, np.ceil(abs(c)) + 1
    return np.unique([np.sign(c) * k / q for q in range(1, qmax + 1) for k in range(int(lo * q), int(hi * q) + 1)])


def consts():
    """Scalar constants of S (and the AP230 damping steps of S, M, L) against 2^-n grids, simple rationals (denominator <= 12) and
    two-term shift-and-add forms (+-2^-a +- 2^-b, a, b <= 12), with each constant's precision and the look-elsewhere chance."""
    import q1_generalise as q1
    import q2_injection as q2
    log = Tee(os.path.join(LB, "chain3_fixedpoint_consts.log"))
    C = []                                                    # (name, value, sigma, sigma's source, group)
    # AP230's feedforward / feedback gains: dense_mid's fits (25.12-25.13), 3 variants x 2 folds; A = (-g_ff + z^-230) / (1 - g_fb z^-230)
    g0 = -0.70710678
    th = [np.load(f)[10:12] for f in sorted(glob.glob(os.path.join(ROOT, "analysis", "dense_mid", "chain3", "theta_*_combo-ck2_ap230*.npy")))]
    th = np.array(th)
    ff, fb = -(g0 + th[:, 0]), -(g0 + th[:, 1])
    for nm, v in (("AP230 |g_ff|", ff), ("AP230 |g_fb|", fb), ("AP230 (|g_ff| + |g_fb|) / 2", (ff + fb) / 2), ("AP230 gap |g_fb| - |g_ff|", fb - ff)):
        C.append((nm, float(np.mean(v)), float(max(np.std(v, ddof=1), 2e-6)), f"sd of {len(v)} fits (3 bases x 2 folds; same data: a lower bound)", "ap"))
    log(f"AP230 fits (dense_mid, 25.13): d g_ff {np.round(th[:, 0], 6).tolist()}, d g_fb {np.round(th[:, 1], 6).tolist()}")
    log(f"   |g_ff| {np.round(ff, 6).tolist()}; |g_fb| {np.round(fb, 6).tolist()}; gap {np.round(fb - ff, 6).tolist()} (25.12's 2-term fits: 4.82e-4 / 4.88e-4)")
    log(f"   reference: 1/sqrt2 = 0.7071068; Q15 23170 / 32768 = {23170 / 32768:.7f}, 23171 = {23171 / 32768:.7f}; Q12 (MPYK) 2896 / 4096 = {2896 / 4096:.7f}, "
        f"2897 = {2897 / 4096:.7f}; 2^-11: 1448 / 2048 = {1448 / 2048:.7f}, 1449 = {1449 / 2048:.7f}")
    # the AP230 damping steps, three programs
    ap = json.load(open(os.path.join(ROOT, "analysis", "q2", "ap230_g0.7071_chain2.json")))
    st = {}
    for v in ap["settings"].values():
        st.setdefault(v["step"], []).append(v["b"])
    steps = {"S": {k: (float(np.median(x)), float(np.std(x, ddof=1) / np.sqrt(len(x))), len(x)) for k, x in st.items()}}
    fp = json.load(open(os.path.join(ROOT, "analysis", "ml_carryover", "firstpass.json")))
    bs = json.load(open(os.path.join(ROOT, "analysis", "ml_lines", "b_steps.json")))
    for prog in "ML":
        st = {}
        for r in fp[prog]["rows"]:
            lf_, rt_, hf_ = r["key"].split("_")
            st.setdefault(bs[prog]["cells"][f"{rt_}_{hf_}"], {})[r["key"]] = r["b"]
        steps[prog] = {k: (float(np.mean(list(x.values()))), float(np.std(list(x.values()), ddof=1) / np.sqrt(len(x))), len(x)) for k, x in st.items()}
    log("\nAP230 damping steps b (S: median of the per-setting fits, chain2, 16.9; M / L: mean, firstpass.json, ml_structure 19), +- SEM, n:")
    for prog in "SML":
        log(f"   {prog}: " + "; ".join(f"step {k}: {steps[prog][k][0]:.5f} +- {steps[prog][k][1]:.1e} ({steps[prog][k][2]})" for k in sorted(steps[prog])))
    for k in (1, 2, 3):
        v = np.array([steps[p][k][0] for p in "SML"])
        sem = np.array([steps[p][k][1] for p in "SML"])
        chi2 = float(np.sum((v - np.average(v, weights=sem ** -2)) ** 2 / sem ** 2))
        log(f"   step {k}: S / M / L spread sd {np.std(v, ddof=1):.1e} against SEMs {sem.max():.0e} or less (chi2 {chi2:.0f} on 2 dof): the programs disagree beyond "
            f"their statistics; offset by each program's step-0 value (S {steps['S'][0][0]:.4f}, M {steps['M'][0][0]:.4f}, L {steps['L'][0][0]:.4f}): "
            f"{', '.join(f'{x:.5f}' for x in v - np.array([steps[p][0][0] for p in 'SML']))} (sd {np.std(v - np.array([steps[p][0][0] for p in 'SML']), ddof=1):.1e})")
        C.append((f"AP step {k}, mean of S / M / L", float(np.mean(v)), float(np.std(v, ddof=1)), "sd across the three programs (per-program bias)", "step"))
        for p in "SML":
            C.append((f"AP step {k}, {p}", steps[p][k][0], float(np.hypot(steps[p][k][1], np.std(v, ddof=1))), "SEM (+) the S/M/L spread", "step_p"))
    sp = {p: (steps[p][2][0] - steps[p][1][0], steps[p][3][0] - steps[p][2][0]) for p in "SML"}
    d12 = np.array([sp[p][0] for p in "SML"]); d23 = np.array([sp[p][1] for p in "SML"])
    log(f"   spacings step 2 - 1 / 3 - 2: " + "; ".join(f"{p} {sp[p][0]:.5f} / {sp[p][1]:.5f}" for p in "SML") +
        f"; mean {d12.mean():.5f} / {d23.mean():.5f}, difference {np.mean(d23 - d12):+.5f} +- {np.std(d23 - d12, ddof=1) / np.sqrt(3):.5f} (S alone: "
        f"{sp['S'][1] - sp['S'][0]:+.5f} +- {np.sqrt(2 * steps['S'][2][1] ** 2 + steps['S'][1][1] ** 2 + steps['S'][3][1] ** 2):.5f}): not equally spaced, so no k / 6 or k d (+ c) table")
    # the law constants (q1.fit_laws, 13 numbers): bootstrap over the law files
    fits = json.load(open(os.path.join(q2.Q1, "loops.json")))
    law = q1.fit_laws(fits)
    names = [(f"law {lf} {x}", (lambda L, lf=lf, x=x: L["lf"][lf][x])) for lf in law["lf"] for x in ("hi", "dc", "k") if not (x == "dc" and lf == 1.0)]
    names.append(("shelf pole a", lambda L: L["a"]))
    bt = []
    for b in range(300):
        L = q1.fit_laws([fits[i] for i in rng.integers(0, len(fits), len(fits))])
        bt.append([f(L) for _, f in names])
    bt = np.array(bt)
    for j, (nm, f) in enumerate(names):
        C.append((nm, float(f(law)), float(np.std(bt[:, j])), "bootstrap over the law files (300)", "law"))
    log(f"\nLaw constants: {law['n_files']} files; f_ref = {law['f_ref']:.0f} Hz is a scan on a 100 Hz grid (not tested); LF 1.0's dc is fixed at 1")
    # the DC floor, the pole ceiling and drop-back, the injection ratio, the M / L output kernels
    sh = json.load(open(os.path.join(LB, "chain3_floor_shoulder.json")))
    for rt in ("0.6", "1.4", "1.8"):
        v = 10 ** (sh["F"][rt] / 20)
        C.append((f"DC floor plateau RT {rt}{' (0.6-1.0)' if rt == '0.6' else ''}", float(v), float(v * (10 ** (0.003 / 20) - 1)), "0.003 dB (34: by HF step)", "floor"))
    log(f"DC floor plateaus (35.1): RT 0.6 / 0.8 / 1.0 " + " / ".join(f"{10 ** (sh['F'][r] / 20):.5f}" for r in ("0.6", "0.8", "1.0")) + " (one value to 3e-5)")
    for nm, v, s_, src in (("pole ceiling", 0.8475, 1e-3, "25.17-25.19: S upper quartile 0.8472, M / L 0.8443-0.8517"),
                           ("drop-back pole LF 0.1 / 0.4 / 2.0", 0.8465, 5e-4, "32.4: the ceiling + an old cell shift"),
                           ("drop-back pole LF 1.0", 0.847, 5e-4, "32.4"), ("drop-back threshold (law pole)", 0.856, 2e-3, "32.4: between law-pole samples"),
                           ("cluster gain / tap magnitude", 26.315, 1e-3, "25.11: chain2 26.316 -> chain3 26.315"),
                           ("tap magnitude / cluster gain", 1 / 26.315, 1e-3 / 26.315 ** 2, "as above")):
        C.append((nm, v, s_, src, "misc"))
    kf = json.load(open(os.path.join(ROOT, "analysis", "ml_lines", "kernel_fix.json")))
    for prog in "ML":
        k, o = np.array(kf[prog]["kernel"]), np.array(kf[prog]["old"])
        C.append((f"{prog} output kernel centre tap", float(k[5]), float(abs(k[5] - o[5])), "18: refit change (library chain inside)", "kernel"))
        C.append((f"{prog} output kernel sum", float(k.sum()), float(abs(k.sum() - o.sum())), "18: refit change", "kernel"))
    # the tests
    log(f"\nPer constant: value +- sigma; on 2^-n grids testable where q >= 4 sigma: the finest testable n*, the offset at n* in LSB, and the")
    log("smallest chance p = 2 |d| / q over the testable grids, Sidak-corrected over their number (p_c); the nearest rational k / q (q <= 12) and two-term")
    log("shift-and-add form (+-2^-a +- 2^-b), each with its chance p (a constant drawn around the value lands as close).")
    log(f"{'constant':38s} {'value':>11s} {'sigma':>8s} | {'n*':>3s} {'d/q@n*':>7s} {'best n':>6s} {'p':>6s} {'p_c':>6s} | {'rational':>10s} {'d':>9s} {'p':>5s} | {'2-term':>10s} {'d':>9s} {'p':>5s}")
    res = []
    for nm, v, s_, src, grp in C:
        tn = [n for n in range(6, 17) if 2.0 ** -n >= 4 * s_]
        gp = []
        for n in tn:
            q = 2.0 ** -n
            dq = (v - np.round(v / q) * q) / q
            gp.append((min(1.0, 2 * abs(dq)), n, dq))
        if gp:
            pbest, nb, _ = min(gp)
            pc = 1 - (1 - pbest) ** len(gp)
            dqs = gp[-1][2]
        else:
            pbest = pc = nb = dqs = np.nan
        R = rationals(v)
        rv, rd = near(v, R)
        rp = chance(v, rd, R, s_)
        cv, cd = near(abs(v), CSD2)
        cp = chance(abs(v), cd, CSD2, s_)
        rq = next(q for q in range(1, 13) if abs(round(rv * q) - rv * q) < 1e-9)
        log(f"{nm:38s} {v:11.6f} {s_:8.1e} | {tn[-1] if tn else '-':>3} {dqs:+7.2f} {nb if gp else '-':>6} {pbest:6.3f} {pc:6.3f} | "
            f"{f'{int(round(rv * rq))}/{rq}':>10s} {rd:+9.1e} {rp:5.2f} | {cv:10.6f} {cd:+9.1e} {cp:5.2f}" + ("  (|d| > 2 sigma)" if abs(rd) > 2 * s_ and abs(cd) > 2 * s_ else ""))
        res.append(dict(name=nm, value=v, sigma=s_, sigma_src=src, group=grp, testable_n=tn, grid_p=pbest, grid_p_corr=pc, grid_best_n=nb,
                        dq_finest=dqs, rational=[float(rv), int(rq)], rational_d=float(rd), rational_p=rp, csd2=float(cv), csd2_d=float(cd), csd2_p=cp))
    from scipy.stats import chi2
    ind = [r for r in res if r["group"] in ("ap", "step", "law", "floor", "misc", "kernel") and r["name"] not in
           ("AP230 (|g_ff| + |g_fb|) / 2", "tap magnitude / cluster gain") and np.isfinite(r["grid_p_corr"])]
    for key in ("grid_p_corr", "rational_p", "csd2_p"):
        ps = np.clip([r[key] for r in ind], 1e-6, 1)
        X = -2 * np.sum(np.log(ps))
        log(f"Fisher over {len(ps)} independent constants, {key}: chi2 {X:.1f} on {2 * len(ps)} dof, p {chi2.sf(X, 2 * len(ps)):.3f}; "
            f"smallest {min(ps):.3f} ({ind[int(np.argmin(ps))]['name']}), Bonferroni {min(1, min(ps) * len(ps)):.2f}")
    # the gap as a number of LSBs of some 2^-n (1, 2 or 3 LSB), on a log-uniform null
    gap = next(r for r in res if r["name"].startswith("AP230 gap"))
    l2 = np.log2(gap["value"])
    d1 = abs(l2 - np.round(l2))
    dl = min(d1, np.min(np.abs(l2 - (np.arange(-20, 0) + np.log2(3)))))
    lo, hi = np.log2(4.82e-4), np.log2(gap["value"] + 2 * gap["sigma"])
    log(f"\nAP230 gap = {gap['value']:.4e} +- {gap['sigma']:.1e} (6 fits; 25.12's 4.82e-4 / 4.88e-4 widen it to 4.8-5.2e-4) = {gap['value'] / 2 ** -11:.3f} LSB of 2^-11 "
        f"(log2 {l2:.3f}; the range of fits spans log2 {lo:.3f} .. {hi:.3f}). Chance that a gap drawn log-uniformly lands as close to one LSB of some 2^-n: "
        f"{min(1, 2 * d1):.2f}; to 1, 2 or 3 LSB: {min(1, 4 * dl):.2f}")
    json.dump(dict(constants=res, steps=steps), open(os.path.join(LB, "chain3_fixedpoint_consts.json"), "w"), indent=1)


def panel():
    """Table indexed by panel steps, or computed? Per loop, the unit's coefficient minus the baseline's smooth part (the HF pole
    minus the held-out pole curve; the shelf's DC / Nyquist gains in dB per pass minus the law without shelf cells), predicted
    for each setting from OTHER settings' take A and scored on its own take B (and B -> A): a table per (HF, RT) step (from the
    other LF steps), per (LF, RT) (other HF steps), per (LF, HF) (other RT steps), or a smooth function of the panel values
    (least squares in the law pole / RT and log LF, from all other settings). The ceiling is the setting's own other take.
    For the shelf also the loop structure: one dB-per-pass offset per setting (a coefficient table entry shared by the loops),
    one offset proportional to the loop length (one RT-like quantity the firmware converts per loop), or per loop."""
    from defects import DEFECTS
    bad = {k for k, _ in DEFECTS}
    d = json.load(open(COEFS))
    P = np.array([1321, 1389, 1461, 1587], float)
    log = Tee(os.path.join(LB, "chain3_fixedpoint_panel.log"))
    log("Panel table or computed (model_s_results 36): held-out R^2 = 1 - sum (r_test - pred)^2 / sum r_test^2, pred from other settings' other take;")
    log("'repeat' = the setting's own other take (the reproducible share); r = unit - baseline's smooth part")
    out = {}

    def score(rows, r_of, cells, smooth_x):
        """rows: list of dict(key, lf, rt, hf, loop, A, B); cells: name -> cell-key function; smooth_x: regressors."""
        res = {}
        for name, fk in list(cells.items()) + [("smooth", None), ("repeat", None)]:
            se, ss = 0.0, 0.0
            for tr, te in (("A", "B"), ("B", "A")):
                for j, r in enumerate(rows):
                    if name == "repeat":
                        pred = r[tr]
                    elif name == "smooth":
                        m = np.array([x["key"] != r["key"] for x in rows])
                        X = np.array([smooth_x(x) for x in rows])
                        y = np.array([x[tr] for x in rows])
                        c = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
                        pred = float(X[j] @ c)
                    else:
                        o = [x[tr] for x in rows if fk(x) == fk(r) and x["key"] != r["key"]]
                        if not o:
                            continue
                        pred = float(np.mean(o))
                    se += (r[te] - pred) ** 2
                    ss += r[te] ** 2
            res[name] = 1 - se / ss
        return res

    # the HF pole against its curve
    rows = []
    for r in d["rows"]:
        if r["base"]["p"] is None or is_drop(r["base"]["p"]) or len(r["take"]) < 2 or r["key"] in bad:
            continue
        A, B = r["take"]["A"]["p"] - r["base"]["p"], r["take"]["B"]["p"] - r["base"]["p"]
        if abs(A - B) > 3e-4:
            continue
        rows.append(dict(key=r["key"], lf=r["lf"], rt=r["rt"], hf=round(max(r["hf"], 0.12 / r["rt"]), 3), loop=r["loop"], A=A, B=B, pb=r["base"]["p"]))
    cells = {"(HF, RT) table": lambda x: (x["hf"], x["rt"], x["loop"]), "(LF, RT) table": lambda x: (x["lf"], x["rt"], x["loop"]),
             "(LF, HF) table": lambda x: (x["lf"], x["hf"], x["loop"])}
    sm = lambda x: [1, x["pb"], x["pb"] ** 2, np.log(x["lf"]), np.log(x["lf"]) * x["pb"], x["loop"]]
    res = score(rows, None, cells, sm)
    rms = np.sqrt(np.mean([(x["A"] + x["B"]) ** 2 / 4 for x in rows]))
    log(f"\nHF pole - curve: {len(rows)} loops (|A - B| < 3e-4, no drop-back, no DEFECTS setting), rms {rms:.2e}")
    log("   " + "; ".join(f"{k} {v:+.2f}" for k, v in res.items()))
    out["pole"] = res
    # the shelf in dB per pass, against the law without cells
    for form in ("D", "N1"):
        rows = []
        for r in d["rows"]:
            if r["base"]["D"] < 0.1905 or pole_limited(r) or len(r["take"]) < 2 or r["key"] in bad:
                continue
            l0 = 20 * np.log10(r["nocell"][form])
            rows.append(dict(key=r["key"], lf=r["lf"], rt=r["rt"], hf=r["hf"], loop=r["loop"], A=20 * np.log10(r["take"]["A"][form]) - l0,
                             B=20 * np.log10(r["take"]["B"][form]) - l0))
        cells = {"(LF, RT) table": lambda x: (x["lf"], x["rt"], x["loop"]), "(HF, RT) table": lambda x: (x["hf"], x["rt"], x["loop"]),
                 "(LF, HF) table": lambda x: (x["lf"], x["hf"], x["loop"])}
        sm = lambda x: [1, 1 / x["rt"], 1 / x["rt"] ** 2, np.log(x["lf"]), np.log(x["lf"]) / x["rt"], x["loop"]]
        res = score(rows, None, cells, sm)
        rms = np.sqrt(np.mean([(x["A"] + x["B"]) ** 2 / 4 for x in rows]))
        log(f"\nShelf {form} (dB per pass) - law without cells: {len(rows)} loops (no floor, no pole-limited HF step), rms {rms:.4f} dB")
        log("   " + "; ".join(f"{k} {v:+.2f}" for k, v in res.items()))
        # loop structure inside a setting: shares of the reproducible part (A x B) carried by a common dB/pass, a P-proportional term, the rest
        keys = sorted({x["key"] for x in rows})
        num = {"common dB/pass": 0.0, "proportional to P": 0.0}
        den = 0.0
        for k in keys:
            g = [x for x in rows if x["key"] == k]
            if len(g) < 4:
                continue
            a_, b_ = np.array([x["A"] for x in g]), np.array([x["B"] for x in g])
            pp = P[[x["loop"] for x in g]] / P.mean()
            for nm, X in (("common dB/pass", np.ones(4)), ("proportional to P", pp)):
                ca = X @ a_ / (X @ X)
                num[nm] += float((ca * X) @ b_)          # held out by take: fitted on A, the part of B it carries
            den += float(a_ @ b_)
        log(f"   loop structure (A fitted, B scored; share of the reproducible sum A.B): " + "; ".join(f"{k} {v / den:.2f}" for k, v in num.items()))
        out[form] = dict(res, **{k: v / den for k, v in num.items()})
    json.dump(out, open(os.path.join(LB, "chain3_fixedpoint_panel.json"), "w"), indent=1)


FS_DEV = 19.2e6 / 728
FXBIN = os.path.join(os.environ.get("TMPDIR", "/tmp"), "s_loops_fx")


def fx_build():
    import subprocess
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "c", "s_loops_fx.c")
    if not os.path.exists(FXBIN) or os.path.getmtime(FXBIN) < os.path.getmtime(src):
        subprocess.run(["cc", "-O2", "-o", FXBIN, src, "-lm"], check=True)


def fx_run(par, mode, x, tmp):
    """par: 4 tuples (P, c0, c1, a, p, b, sign); mode 0 float, 1 truncation, 2 rounding; x: input (full scale 1)."""
    import subprocess
    pf, xf, of = (os.path.join(tmp, f) for f in ("par.txt", "x.f64", "y.f64"))
    with open(pf, "w") as f:
        f.write(f"{len(x)}\n" + "".join(" ".join(f"{v:.12g}" for v in r) + "\n" for r in par))
    np.asarray(x, np.float64).tofile(xf)
    r = subprocess.run([FXBIN, pf, str(mode), xf, of], check=True, capture_output=True, text=True)
    nclip = int(r.stdout.strip() or 0)
    if nclip:
        print(f"   [{nclip} samples saturated, mode {mode}]")
    return np.fromfile(of, np.float64)


def q10(keys=("1.0_2.5_1.0", "1.0_1.0_1.0", "2.0_0.8_0.2")):
    """Q10: S's four loops (the baseline's per-loop shelf, gain and HF pole; one read tap per loop) in float and in 16-bit fixed point
    with a 32-bit accumulator, storage truncated or rounded (c/s_loops_fx.c). Internal scaling (an assumption): the loops' input gain
    set so the loops' impulse response has unit energy (white noise passes at unity power gain). (1) Impulse response at -6 dBFS
    peak: the error fixed - float per window, re the 50-300 ms level, and the fixed tail itself (floor, dead band, limit cycle, DC);
    (2) idle noise: ADC-like white noise at -94 dBFS rms in (the DELAY program's 88 dB below nominal, nominal -6 dBFS), the output
    noise re nominal against the REVERB spec (73 dB); (3) what a sweep measurement leaves: an exponential sweep (10 Hz - 13 kHz,
    time constant 3.4 s, about 24 s, as ml_structure 9.5's simulation; Tomkin's sweep is unknown) at -6 dBFS peak through float and
    fixed models with the same input noise, deconvolved: the fixed-point and input-noise floors re the deconvolved 50-300 ms level."""
    import tempfile
    from scipy.signal import butter, sosfilt
    HP = butter(2, 20, btype="high", fs=FS_DEV, output="sos")         # the analogue chain's high-pass: DC never reaches a capture
    hp = lambda v: sosfilt(HP, v)
    fx_build()
    d = json.load(open(COEFS))
    P = [1321, 1389, 1461, 1587]
    sg = [1, -1, 1, -1]
    log = Tee(os.path.join(LB, "chain3_fixedpoint_q10.log"))
    log("Q10 (model_s_results 36): S's four loops, float vs 16-bit fixed point (32-bit accumulator, Q15 coefficients), truncation vs rounding")
    rng2 = np.random.default_rng(10)
    db = lambda v: 10 * np.log10(np.mean(np.asarray(v) ** 2) + 1e-300)
    out = {}
    W = [(0.05, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0), (5.0, 6.0)]
    with tempfile.TemporaryDirectory() as tmp:
        for key in keys:
            rows = [r for r in d["rows"] if r["key"] == key]
            base = [r["base"] for r in sorted(rows, key=lambda r: r["loop"])]
            par0 = [(P[i], base[i]["c0"], base[i]["c1"], d["a"], base[i]["p"] or 0.0, 1.0, sg[i]) for i in range(4)]
            N = int(6.0 * FS_DEV)
            imp = np.zeros(N); imp[0] = 1.0
            h = fx_run(par0, 0, imp, tmp)
            gmax = np.array([max(b_["D"], b_["N1"]) for b_ in base])          # the loop filter's largest gain (shelf x one-pole: DC or Nyquist)
            b_A = 1 / np.sum(1 / (1 - gmax))                                   # A: no loop write and no output can overflow for any input <= FS
            b_B = 1 / np.sqrt(np.sum(h ** 2))                                  # B: the loops' IR has unit energy (white noise at unity power gain)
            log(f"\n{key}: loops (c0, c1, p) " + "; ".join(f"({t[1]:.4f}, {t[2]:.4f}, {t[4]:.4f})" for t in par0) + f"; input gain A (no overflow) {b_A:.4f}, "
                f"B (unit-energy IR) {b_B:.4f}: B is {20 * np.log10(b_B / b_A):+.1f} dB hotter")
            o = dict(b_A=b_A, b_B=b_B)
            for scl, bin_ in (("A", b_A), ("B", b_B)):
                par = [t[:5] + (bin_,) + t[6:] for t in par0]
                log(f"  scaling {scl}:")
                o[scl] = oo = {}
                # (1) impulse at -6 dBFS
                x = 0.5 * imp
                yr = {m: fx_run(par, m, x, tmp) for m in (0, 1, 2)}
                y = {m: hp(v) for m, v in yr.items()}
                ref = db(y[0][int(0.05 * FS_DEV):int(0.3 * FS_DEV)])
                seg = lambda v, a_, b_: v[int(a_ * FS_DEV):int(b_ * FS_DEV)]
                log(f"   (1) impulse -6 dBFS peak, outputs high-passed at 20 Hz; 50-300 ms level {ref:.1f} dBFS; per window, dB re that: float | trunc: error, fixed tail, "
                    f"tail DC before the high-pass | round: error, fixed tail")
                oo["impulse"] = {}
                for a_, b_ in W:
                    f0, t1, r2 = seg(y[0], a_, b_), seg(y[1], a_, b_), seg(y[2], a_, b_)
                    row = dict(float=db(f0) - ref, trunc_err=db(t1 - f0) - ref, trunc=db(t1) - ref, trunc_dc=20 * np.log10(abs(np.mean(seg(yr[1], a_, b_))) + 1e-300) - ref,
                               round_err=db(r2 - f0) - ref, round=db(r2) - ref)
                    oo["impulse"][f"{a_}-{b_}"] = row
                    log(f"      {a_:.2f}-{b_:.2f} s: {row['float']:7.1f} | {row['trunc_err']:7.1f} {row['trunc']:7.1f} {row['trunc_dc']:7.1f} | {row['round_err']:7.1f} {row['round']:7.1f}")
                t_last, r_last = seg(yr[1], 5.0, 6.0), seg(yr[2], 5.0, 6.0)
                for nm, v in (("trunc", t_last), ("round", r_last)):
                    lv = np.round(v * 32768).astype(int)
                    sp = np.abs(np.fft.rfft(v - v.mean())) ** 2
                    fpk = np.fft.rfftfreq(len(v), 1 / FS_DEV)[np.argmax(sp)] if np.any(v != v.mean()) else 0.0
                    log(f"      5-6 s, {nm}: output values {sorted(set(lv.tolist()))[:6]}{'...' if len(set(lv.tolist())) > 6 else ''} (LSB), mean {v.mean() * 32768:+.2f} LSB, "
                        f"rms about the mean {np.std(v) * 32768:.2f} LSB" + (f", spectral peak {fpk:.0f} Hz (a limit cycle)" if np.std(v) > 0 else ""))
                    oo[f"tail_{nm}"] = dict(mean_lsb=float(v.mean() * 32768), rms_lsb=float(np.std(v) * 32768), peak_hz=float(fpk))
                # (2) idle noise
                M = int(10 * FS_DEV)
                xn = rng2.standard_normal(M) * 10 ** (-94 / 20)
                ynr = {m: fx_run(par, m, xn, tmp) for m in (0, 1, 2)}
                yn = {m: hp(v) for m, v in ynr.items()}
                nom = -6.0
                sn = {m: nom - db(yn[m][M // 2:]) for m in (0, 1, 2)}
                add = {m: nom - db((yn[m] - yn[0])[M // 2:]) for m in (1, 2)}
                fs = lambda x: "silent (dead band)" if x > 200 else f"{x:.1f}"   # db()'s 1e-300 floor: an all-zero output
                log(f"   (2) idle noise in, -94 dBFS rms (88 dB below nominal): output S/N re nominal: float {fs(sn[0])} dB, trunc {fs(sn[1])}, round {fs(sn[2])} "
                    f"(the fixed-point part alone: trunc {fs(add[1])}, round {fs(add[2])}; trunc output DC {np.mean(ynr[1][M // 2:]) * 32768:+.1f} LSB, removed by the 20 Hz high-pass); spec REVERB 73, DELAY 88")
                oo["idle"] = dict(float=sn[0], trunc=sn[1], round=sn[2], trunc_added=add[1], round_added=add[2], trunc_dc_lsb=float(np.mean(ynr[1][M // 2:]) * 32768))
                # (3) through a sweep, deconvolved
                Ls, f1, f2 = 3.4, 10.0, 13000.0
                T = Ls * np.log(f2 / f1)
                t = np.arange(int(T * FS_DEV)) / FS_DEV
                sw = 0.5 * np.sin(2 * np.pi * f1 * Ls * (np.exp(t / Ls) - 1))
                sw[:200] *= np.hanning(400)[:200]; sw[-200:] *= np.hanning(400)[200:]
                K = int(7 * FS_DEV)
                xs = np.r_[sw, np.zeros(K)]
                noise = rng2.standard_normal(len(xs)) * 10 ** (-94 / 20)
                Nf = 1 << int(np.ceil(np.log2(len(xs) + K)))
                S = np.fft.rfft(sw, Nf)
                f = np.fft.rfftfreq(Nf, 1 / FS_DEV)
                band = np.clip((f - 15) / 10, 0, 1) * np.clip((12500 - f) / 1500, 0, 1)
                band = 0.5 - 0.5 * np.cos(np.pi * band)                # raised-cosine band edges (15-25 Hz, 11-12.5 kHz)
                inv = band / np.where(np.abs(S) > 0, S, 1)
                dec = lambda v: np.fft.irfft(np.fft.rfft(v, Nf) * inv, Nf)[:K]
                if scl == "B":                                   # B clips on the sweep's resonant sines: the sweep is scored on A only
                    continue
                yf = fx_run(par, 0, xs, tmp)
                hf = dec(yf)
                ref3 = db(seg(hf, 0.05, 0.3))
                log(f"   (3) the sweep's output peak {20 * np.log10(np.abs(yf).max()):.1f} dBFS (the headroom a hotter sweep or a hotter internal scaling could use)")
                oo["sweep"] = {}
                log(f"   (3) sweep {T:.1f} s at -6 dBFS peak, deconvolved 25 Hz-11 kHz (tapered), re its 50-300 ms level: input noise alone | trunc: total floor, fixed-point part | round: total, part")
                hn = dec(fx_run(par, 0, xs + noise, tmp))
                ht = dec(fx_run(par, 1, xs + noise, tmp))
                hr = dec(fx_run(par, 2, xs + noise, tmp))
                for a_, b_ in W:
                    row = dict(noise=db(seg(hn - hf, a_, b_)) - ref3, trunc=db(seg(ht - hf, a_, b_)) - ref3, trunc_part=db(seg(ht - hn, a_, b_)) - ref3,
                               round=db(seg(hr - hf, a_, b_)) - ref3, round_part=db(seg(hr - hn, a_, b_)) - ref3, ir=db(seg(hf, a_, b_)) - ref3)
                    oo["sweep"][f"{a_}-{b_}"] = row
                    log(f"      {a_:.2f}-{b_:.2f} s: IR {row['ir']:7.1f} | {row['noise']:7.1f} | {row['trunc']:7.1f} {row['trunc_part']:7.1f} | {row['round']:7.1f} {row['round_part']:7.1f}")
            out[key] = o
    json.dump(out, open(os.path.join(LB, "chain3_fixedpoint_q10.json"), "w"), indent=1)


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
        coefs(*([os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3])] if len(sys.argv) > 3 else []))
    elif st == "lattice":
        lattice()
    elif st == "consts":
        consts()
    elif st == "panel":
        panel()
    elif st == "q10":
        q10(*([sys.argv[2:]] if len(sys.argv) > 2 else []))
    elif st == "snap":
        snap(int(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 0.0)
