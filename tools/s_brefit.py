"""S's AP230 damping steps b refitted on the current S model (model_s_results.md section 38; the S analogue of ml_structure 22).

S's AP230 table (analysis/q2/ap230_g0.7071_chain2.json, steps 0.00069 / 0.17109 / 0.33160 / 0.50192) was read on the chain2
IRs with the injection and loops of that time (section 16.9 / 21). Since then the IRs went to chain3 and the loops to the
section 35 baseline; M / L's steps, refitted the same way on their current model, moved by 0.6-3.4e-3 and fell on one table
(ml_structure 22). Here each S setting is rendered as pulse_check's s stage does (the fold that left its RT step out, the
cells, per capture clock and gain refitted for every b) on a grid of AP230 b, and every window's null (q2.W_COMP, L and R)
is kept per capture, so any candidate is scored without rendering again.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 .venv/bin/python tools/memguard.py --limit-gb 9 -- \\
        .venv/bin/python tools/s_brefit.py --stage fit [--only N] [--workers 8]
        per setting -> analysis/s_brefit/<key>.json (resumable): damped settings (step 1-3) on b0 - 0.006 ... b0 + 0.002
        (0.0005), undamped (step 0) on 0 / 0.00035 / 0.00069; pipeline check: the stored b against analysis/pulse/s
    .venv/bin/python tools/s_brefit.py --stage cand   -> analysis/s_brefit/cand.{json,log}: per-setting b fitted on one
        take's 0-50 ms (AMPL 0 + 90, or 30 + 60) and scored on the other take; step tables (per-setting optima averaged
        over the step, all data and leave-one-LF-column-out); step 0 at b = 0; the M / L table (ml_structure 22).
        Every window 0-1000 ms, R and L, clean captures (defects.DEFECTS), against the stored b, with the strict rule;
        writes analysis/s_brefit/b_<name>.json {key: b} for each candidate
"""
import argparse, collections, glob, json, os, sys

for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
OUT = os.path.join(ROOT, "analysis", "s_brefit")
WINS = ["0-50", "50-150", "150-300", "300-600", "600-1000"]
TAKES = ((0, 90), (30, 60))
ML_TABLE = [0.0, 0.169859, 0.329683, 0.499595]      # ml_structure 22


def grid_of(step, b0):
    if step == 0:
        return [0.0, 0.00035, b0]
    return [b0 + 0.0005 * i for i in range(-12, 5)]              # i = 0 is b0 itself (the stored value, exactly)


def _job(a):
    key, path = a
    import pulse_check as pc
    q1, q2 = pc.q1, pc.q2
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    g = info[key]
    fold = json.load(open(os.path.join(pc.CV, f"fold_{g['rt']}.json")))
    lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
    caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(g)]
    step = ap["settings"][key]["step"]
    b0 = ap["steps"][step]
    grid = grid_of(step, b0)
    N = {}
    for b in grid:
        G = pc.setting_G(lab, b, law, fold, cells.get(key))
        row = {}
        for c in caps:
            y, mdl, gg, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
            nl = q1.nulls(y, mdl * gg, te, q2.W_COMP, narrow=True)
            row[str(c["amp"])] = {w: nl[w] for w in WINS}
        N[repr(b)] = row
    ref = os.path.join(ROOT, "analysis", "pulse", "s", key + ".json")
    dev = None
    if os.path.exists(ref):                          # the pipeline check: the stored b against the pulse check
        R = {str(r["amp"]): r["null"] for r in json.load(open(ref))["rows"]}
        dev = max(abs(R[a][w][ch] - N[repr(b0)][a][w][ch]) for a in R for w in WINS if R[a].get(w) for ch in (0, 1))
    json.dump(dict(key=key, **lab, step=step, b0=b0, grid=grid, nulls=N, pulse_check_dev=dev), open(path, "w"))
    print(f"{key}: step {step}, {len(grid)} b, pulse dev {dev}", flush=True)
    return key


def fit_stage(workers, only=0):
    from concurrent.futures import ProcessPoolExecutor
    import q2_injection as q2
    os.makedirs(OUT, exist_ok=True)
    info = q2.q1_tables()[0]
    ap = json.load(open(os.path.join(q2.OUT, "ap230_g0.7071_chain2.json")))
    keys = sorted(info, key=lambda k: (-ap["settings"][k]["step"], k))      # damped first (the long jobs)
    if only:
        keys = [k for k in keys if ap["settings"][k]["step"] > 0][:only] + [k for k in keys if ap["settings"][k]["step"] == 0][:1]
    jobs = [(k, os.path.join(OUT, f"{k}.json")) for k in keys if not os.path.exists(os.path.join(OUT, f"{k}.json"))]
    print(f"{len(jobs)} settings to scan", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_job, jobs):
            pass


def load():
    from defects import DEFECTS
    D = {}
    for p in sorted(glob.glob(os.path.join(OUT, "[0-9]*_*_*.json"))):
        d = json.load(open(p))
        d["clean"] = [a for a in d["nulls"][repr(d["b0"])] if (d["key"], int(a)) not in DEFECTS]
        D[d["key"]] = d
    return D


def at(d, b, a, w, ch):
    """Null (dB) of capture a, window w, channel ch at b: cubic in linear energy over the grid (exact at grid points)."""
    from scipy.interpolate import CubicSpline
    g = np.array(d["grid"])
    e = np.array([10 ** (d["nulls"][repr(x)][a][w][ch] / 10) if d["nulls"][repr(x)][a][w] else np.nan for x in d["grid"]])
    if not np.all(np.isfinite(e)):
        return None
    if not g[0] - 1e-9 <= b <= g[-1] + 1e-9:
        raise ValueError(f"{d['key']}: b {b} outside {g[0]}..{g[-1]}")
    if len(g) < 4:
        return float(10 * np.log10(np.interp(b, g, e)))
    return float(10 * np.log10(CubicSpline(g, e)(b)))


def best_b(d, amps, wins=("0-50",)):
    """b minimising the summed L + R energy of the captures amps in wins (fine scan of the spline)."""
    g = np.array(d["grid"])
    bs = np.linspace(g[0], g[-1], 1 + int(round((g[-1] - g[0]) / 1e-5)))
    tot = np.zeros(len(bs))
    for a in amps:
        for w in wins:
            for ch in (0, 1):
                v = [at(d, b, a, w, ch) for b in bs] if at(d, g[0], a, w, ch) is not None else None
                if v is not None:
                    tot += 10 ** (np.array(v) / 10)
    return float(bs[int(np.argmin(tot))])


def cand_stage():
    D = load()
    lines, res, cand = [], {}, collections.defaultdict(dict)
    dev = [d["pulse_check_dev"] for d in D.values() if d["pulse_check_dev"] is not None]
    lines.append(f"{len(D)} settings; pipeline check (stored b vs analysis/pulse/s, every window L / R) max |dev| {max(dev):.4f} dB")
    damped = {k: d for k, d in D.items() if d["step"] > 0 and d["clean"]}      # 1.0_0.3_0.4 has no clean capture
    per_take, both = {}, {}
    for k, d in damped.items():                     # per setting: fitted on one take's 0-50 ms, both takes, and per take
        tk = [[str(a) for a in t if str(a) in d["clean"]] for t in TAKES]
        per_take[k] = [best_b(d, t) if t else None for t in tk]
        both[k] = best_b(d, d["clean"])
    steps = sorted({d["step"] for d in damped.values()})
    st = {s: float(np.mean([both[k] for k, d in damped.items() if d["step"] == s])) for s in steps}
    lfs = sorted({d["lf"] for d in damped.values()})
    fold = {lf: {s: float(np.mean([both[k] for k, d in damped.items() if d["step"] == s and d["lf"] != lf] or [st[s]])) for s in steps} for lf in lfs}
    lines.append("damped settings, b fitted on 0-50 ms (L + R, clean captures of both takes): stored -> refit per step: " + "; ".join(
        f"step {s}: {D[next(k for k in damped if damped[k]['step'] == s)]['b0']:.5f} -> {st[s]:.5f} (sd {np.std([both[k] for k in damped if damped[k]['step'] == s]):.5f}, "
        f"n {sum(1 for k in damped if damped[k]['step'] == s)})" for s in steps))
    lines.append("  take A - take B optimum: median " + f"{np.median([p[0] - p[1] for p in per_take.values() if None not in p]):+.5f}, "
                 f"MAD {np.median(np.abs([p[0] - p[1] for p in per_take.values() if None not in p])):.5f}")
    for s in steps:
        x = sorted(((d["lf"], d["rt"], d["hf"], round(both[k], 5)) for k, d in damped.items() if d["step"] == s))
        lines.append(f"  step {s}: " + " ".join(f"{lf}/{rt}/{hf}:{b}" for lf, rt, hf, b in x))
    edge = [k for k, d in damped.items() if min(abs(both[k] - d["grid"][0]), abs(both[k] - d["grid"][-1])) < 2e-4]
    lines.append(f"  optima at the grid's edge (widen the grid): {edge}")
    for k, d in D.items():                          # candidates: {key: (b per capture)}; held out by take for per-setting
        cand["stored"][k] = {a: d["b0"] for a in d["clean"]}
        cand["step0_zero"][k] = {a: (0.0 if d["step"] == 0 else d["b0"]) for a in d["clean"]}
        if d["step"] > 0 and k in damped:
            ho = {}
            for i, t in enumerate(TAKES):           # capture in take i scored with the other take's fit
                for a in t:
                    if str(a) in d["clean"]:
                        o = per_take[k][1 - i]
                        ho[str(a)] = o if o is not None else both[k]
            cand["setting_heldout_take"][k] = ho
            cand["steps_all"][k] = {a: st[d["step"]] for a in d["clean"]}
            cand["steps_cv"][k] = {a: fold[d["lf"]][d["step"]] for a in d["clean"]}
            cand["ml_table"][k] = {a: ML_TABLE[d["step"]] for a in d["clean"]}
        for nm in ("setting_heldout_take", "steps_all", "steps_cv", "ml_table"):
            if d["step"] == 0:
                cand[nm][k] = {a: 0.0 for a in d["clean"]}
                cand[nm + "_keep0"][k] = {a: d["b0"] for a in d["clean"]}       # the damped steps alone, step 0 as stored
            elif k in cand[nm]:
                cand[nm + "_keep0"][k] = cand[nm][k]
    for nm, C in cand.items():
        if nm == "stored":
            continue
        json.dump({k: float(np.mean(list(v.values()))) if v else None for k, v in C.items()}, open(os.path.join(OUT, f"b_{nm}.json"), "w"), indent=1)
        lines.append(f"\n## candidate {nm}: clean captures, stored b -> candidate; medians, dB (R / L)")
        rec = collections.defaultdict(list)
        for k, d in D.items():
            for a in d["clean"]:
                for w in WINS:
                    x0 = [at(d, cand["stored"][k][a], a, w, ch) for ch in (0, 1)]
                    x1 = [at(d, C[k][a], a, w, ch) for ch in (0, 1)]
                    if None in x0 or None in x1:
                        continue
                    rec[w].append(dict(key=k, amp=a, step=d["step"], lf=d["lf"], rt=d["rt"], hf=d["hf"], b=x0, a=x1))
        bad = collections.defaultdict(dict)
        for w in WINS:
            r = rec[w]
            for sel, f in (("all", lambda x: True), ("damped", lambda x: x["step"] > 0), ("undamped", lambda x: x["step"] == 0)):
                x = [y for y in r if f(y)]
                if not x:
                    continue
                b_, a_ = np.array([y["b"] for y in x]), np.array([y["a"] for y in x])
                lines.append(f"- {w} ms {sel} ({len(x)}): R {np.median(b_[:, 1]):.2f} -> {np.median(a_[:, 1]):.2f}, L {np.median(b_[:, 0]):.2f} -> "
                             f"{np.median(a_[:, 0]):.2f}; worst capture R {np.max(a_[:, 1] - b_[:, 1]):+.2f} L {np.max(a_[:, 0] - b_[:, 0]):+.2f}; "
                             f"better > 0.3 dB R {int(np.sum(a_[:, 1] - b_[:, 1] < -0.3))}, worse > 0.3 dB R {int(np.sum(a_[:, 1] - b_[:, 1] > 0.3))}")
            g = collections.defaultdict(list)
            for y in r:
                g[y["key"]].append(max(y["a"][1] - y["b"][1], y["a"][0] - y["b"][0]))
            for k, v in g.items():
                if np.mean(v) > 0.5:
                    bad[k][w] = round(float(np.mean(v)), 2)
            for gname in ("lf", "rt", "hf"):
                gg = collections.defaultdict(lambda: ([], []))
                for y in r:
                    gg[y[gname]][0].append(y["b"][1])
                    gg[y[gname]][1].append(y["a"][1])
                med_worse = [v for v, (b_, a_) in gg.items() if np.median(a_) > np.median(b_) + 0.1]
                if med_worse:
                    lines.append(f"  - {w} ms: {gname} group medians worse by > 0.1 dB (R): {sorted(med_worse)}")
        lines.append(f"  - settings worse than 0.5 dB (mean over clean captures, max of L / R): {len(bad)} {dict(sorted(bad.items()))}")
        res[nm] = {w: rec[w] for w in WINS}
    json.dump(dict(steps_all=st, fold=fold, both=both, per_take=per_take), open(os.path.join(OUT, "cand.json"), "w"), indent=1)
    open(os.path.join(OUT, "cand.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["fit", "cand"], required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", type=int, default=0)
    a = ap.parse_args()
    if a.stage == "fit":
        fit_stage(a.workers, a.only)
    else:
        cand_stage()


if __name__ == "__main__":
    main()
