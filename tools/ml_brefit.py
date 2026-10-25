"""M / L's damped all-pass steps b refitted on the current model (ml_structure.md section 22; NEXT_STEPS audit item 6).

The 96 damped settings' b (M 56, L 40; firstpass.json) were fitted on 22 October by ml_carryover's firstpass stage: per
setting on the high-AMPL capture (M AMPL 99, L AMPL 80), the summed L + R first-pass residual energy, one shift and a
cluster and a tap gain per capture, with the kernel and IRs of that time (bounded Brent, xatol 0.002). Since then the IRs
went to chain3 and section 18 set the undamped step to b = 0 and refitted the output kernel on the undamped settings only.
Here b is refitted the same way on the current model: ml_lines.Model (kernel_fix kernel, held fixed: it never saw a
damped setting), each capture with the law stage's clock and gains (as pulse_check / ml_lines.fitted), one gain per window
and channel (as pulse_check._ml_nulls), b minimising the high-AMPL capture's L + R first-pass residual energy.

    .venv/bin/python tools/memguard.py --limit-gb 9 -- .venv/bin/python tools/ml_brefit.py --stage fit [--only N] [--workers 8]
        per damped setting -> analysis/ml_lines/brefit/<prog>_<key>.json (resumable): the stored b, the refit (hi capture,
        exact Brent), per-capture refits, the first-pass residual energy per capture and channel on a b grid around both,
        cluster-and-gap nulls at both b
    .venv/bin/python tools/ml_brefit.py --stage write      -> analysis/ml_lines/b_refit.json {prog: {key: b}} (the refit;
        ml_lines reads it with ML_BREFIT=1) and b_refit.log (how b moved, per cell, per LF column, step means)
    .venv/bin/python tools/ml_brefit.py --stage cand       -> analysis/ml_lines/b_refit_cand.{json,log}: candidates scored on
        the grids (first pass R and L per capture; refit, step tables all-data / leave-one-LF-column-out, S's table, one shared
        table) against the stored b, with the strict rule; writes each candidate's b file analysis/ml_lines/b_cand_<name>.json
"""
import argparse, collections, glob, json, os, sys

os.environ["ML_BREFIT"] = "0"                      # the stored b are the reference (b0)
os.environ["ML_BSTEP"] = "0"
os.environ.setdefault("EMT_KERNEL", "chain3")       # the IRs the law files and analysis/pulse/ml were made with
for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
OUT = os.path.join(ROOT, "analysis", "ml_lines")
DIR = os.path.join(OUT, "brefit")
PAD = 3000                                          # render past the first pass's end (band filter, warp)
STEP, HALF = 0.00025, 0.004                         # b grid: spacing, and its reach beyond [min, max] of stored / refit
S_STEPS = [0.0, 0.17109, 0.33160, 0.50192]          # S's AP230 table (model_s_results 16.9, ml_structure 19-20)


def step_of(b):
    return 0 if b < 0.05 else (1 if b < 0.25 else 2 if b < 0.42 else 3)


def damped(prog):
    import ml_lines as mll
    M = mll.Model(prog)
    return sorted(k for k, b in M.b.items() if b >= 0.05)


def _job(a):
    prog, key, path = a
    import ml_carryover as cc
    import ml_lines as mll
    from scipy.optimize import minimize_scalar
    K0 = mll.K0
    M = mll.Model(prog)
    b0 = M.b[key]
    rows = {r["amp"]: r for r in json.load(open(mll.group_path(prog, key)))}
    g = cc.groups(prog)[key]
    caps, _, amps = cc.captures(prog, key, max(len(cc._IR(f)) for f in g.values()) + 200)
    W = mll.windows(prog)
    n = K0 + max(W[ch]["first_pass"][1] for ch in (0, 1)) + PAD
    ends = {a: mll.usable_end(caps[a], g[a], prog) for a in amps}
    ys = {a: np.array(caps[a][:n]) for a in amps}
    del caps
    hi = amps[-1]
    C = M.parts(n, key)[0]
    Cw = {a: mll.warp(C, *rows[a]["clock"]) for a in amps}
    cache = {}

    def model(b):
        M.b[key] = b
        R = M.parts(n, key)[1]
        out = {}
        for a in amps:
            Rw = mll.warp(R, *rows[a]["clock"])
            gs = rows[a]["gains"]
            out[a] = np.stack([gs[ch][0] * Cw[a][:, ch] + gs[ch][1] * Rw[:, ch] for ch in (0, 1)], 1)
        return out

    def rel(a, mdl, w, ch):
        a0, a1 = K0 + W[ch][w][0], min(K0 + W[ch][w][1], ends[a][ch])
        if a1 - a0 <= 1000:                         # as pulse_check._ml_nulls: too little before the fade (M AMPL 10)
            return float("nan")
        x1, x2 = ys[a][a0:a1, ch], mdl[a0:a1, ch]
        return float(np.sum((x1 - np.dot(x1, x2) / np.dot(x2, x2) * x2) ** 2) / np.sum(x1 ** 2))

    def E(b):                                       # first-pass residual energy (relative) per capture and channel
        b = float(b)
        if b not in cache:
            m = model(b)
            cache[b] = {a: [rel(a, m[a], "first_pass", ch) for ch in (0, 1)] for a in amps}
        return cache[b]
    unscored = not np.all(np.isfinite(E(b0)[hi]))  # no usable first pass in the fit capture (M 0.1_0.6_0.7: none in any)
    bf = b0 if unscored else float(minimize_scalar(lambda b: sum(E(b)[hi]), bounds=(b0 - 0.01, b0 + 0.01), method="bounded",
                                                   options={"xatol": 1e-7}).x)
    lo, up = min(b0, bf) - HALF, max(b0, bf) + HALF
    grid = np.round(np.arange(np.floor(lo / STEP) * STEP, up + STEP / 2, STEP), 8)
    Eg = np.array([[E(b)[a] for a in amps] for b in grid])          # (nb, ncap, 2)
    per_cap = {}
    for i, a in enumerate(amps):                    # each capture's own optimum, L + R and R alone, from the grid (quadratic)
        if not np.all(np.isfinite(Eg[:, i])):
            continue
        for name, s in (("LR", Eg[:, i].sum(1)), ("R", Eg[:, i, 1])):
            j = int(np.argmin(s))
            j = min(max(j, 1), len(grid) - 2)
            c = np.polyfit(grid[j - 1:j + 2] - grid[j], s[j - 1:j + 2], 2)
            per_cap.setdefault(str(a), {})[name] = float(grid[j] - c[1] / (2 * c[0])) if c[0] > 0 else float(grid[j])
    cg = {}
    for lab, b in (("b0", b0), ("bf", bf)):
        m = model(b)
        cg[lab] = {str(a): [10 * np.log10(rel(a, m[a], "cluster_gap", ch)) for ch in (0, 1)] for a in amps}
    fp = {lab: {str(a): [10 * np.log10(x) for x in E(b)[a]] for a in amps} for lab, b in (("b0", b0), ("bf", bf))}
    pc = os.path.join(ROOT, "analysis", "pulse", "ml", f"{prog}_{key}.json")
    check = None
    if os.path.exists(pc):                          # the pulse check's R first pass at the stored b: the pipeline check
        P = {r["amp"]: r["null"]["first_pass"] for r in json.load(open(pc))["rows"]}
        check = max((abs(P[a] - fp["b0"][str(a)][1]) for a in amps if P.get(a) is not None), default=None)
    json.dump(dict(prog=prog, key=key, kernel=os.environ["EMT_KERNEL"], amps=amps, hi=hi, b0=b0, bf=bf, unscored=unscored, per_cap=per_cap, grid=grid.tolist(),
                   E=Eg.tolist(), cluster_gap=cg, first_pass=fp, pulse_check_dev=check, renders=len(cache)), open(path, "w"))
    print(f"{prog} {key}: b {b0:.5f} -> {bf:.5f} ({bf - b0:+.5f}); R first pass hi {fp['b0'][str(hi)][1]:.2f} -> "
          f"{fp['bf'][str(hi)][1]:.2f}; pulse dev {check}", flush=True)
    return key


def fit_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    os.makedirs(DIR, exist_ok=True)
    jobs = []
    for prog in ("M", "L"):
        ks = damped(prog)
        if only:
            ks = ks[::max(1, len(ks) // only)][:only]
        jobs += [(prog, k, os.path.join(DIR, f"{prog}_{k}.json")) for k in ks]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    print(f"{len(jobs)} settings to fit", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=2) as ex:
        for _ in ex.map(_job, jobs):
            pass


def load():
    return {(d["prog"], d["key"]): d for d in (json.load(open(p)) for p in sorted(glob.glob(os.path.join(DIR, "*.json"))))}


def split(key):
    lf, rt, hf = key.split("_")
    return lf, rt, hf


def write_stage():
    D = load()
    out, lines = {"M": {}, "L": {}}, []
    for prog in ("M", "L"):
        rr = [d for (p, _), d in D.items() if p == prog]
        for d in rr:
            out[prog][d["key"]] = d["bf"]
        lines.append(f"\n{prog}: unscored (no usable first pass, stored b kept, left out below): {[d['key'] for d in rr if d.get('unscored')]}")
        rr = [d for d in rr if not d.get("unscored")]
        db = np.array([d["bf"] - d["b0"] for d in rr])
        dev = [d["pulse_check_dev"] for d in rr if d["pulse_check_dev"] is not None]
        lines.append(f"\n## {prog}: {len(rr)} damped settings; refit - stored: mean {db.mean():+.5f}, sd {db.std():.5f}, "
                     f"range {db.min():+.5f} .. {db.max():+.5f}; pipeline check (R first pass at stored b vs pulse_check) max |dev| {max(dev):.3f} dB")
        lines.append("| step | n | stored mean +- SEM | refit mean +- SEM | refit sd | refit - stored |")
        lines.append("|---|---|---|---|---|---|")
        for s in (1, 2, 3):
            a = np.array([d["b0"] for d in rr if step_of(d["b0"]) == s])
            b = np.array([d["bf"] for d in rr if step_of(d["b0"]) == s])
            lines.append(f"| {s} | {len(a)} | {a.mean():.5f} +- {a.std(ddof=1) / np.sqrt(len(a)):.5f} | {b.mean():.5f} +- "
                         f"{b.std(ddof=1) / np.sqrt(len(b)):.5f} | {b.std(ddof=1):.5f} | {b.mean() - a.mean():+.5f} |")
        lines.append("- per (RT, HF) cell: step, mean(refit - stored), sd over LF of stored / refit:")
        cells = collections.defaultdict(list)
        for d in rr:
            lf, rt, hf = split(d["key"])
            cells[(float(rt), float(hf))].append(d)
        for c in sorted(cells):
            x = cells[c]
            lines.append(f"  - RT {c[0]} HF {c[1]}: step {step_of(x[0]['b0'])}, n {len(x)}, move {np.mean([d['bf'] - d['b0'] for d in x]):+.5f}, "
                         f"sd stored {np.std([d['b0'] for d in x]):.5f} / refit {np.std([d['bf'] for d in x]):.5f}")
        for lab in ("b0", "bf"):
            sm = {s: np.mean([d[lab] for d in rr if step_of(d["b0"]) == s]) for s in (1, 2, 3)}
            off = collections.defaultdict(list)
            for d in rr:
                lf, rt, hf = split(d["key"])
                o = d[lab] - sm[step_of(d["b0"])]
                off[lf].append(o)
                if prog == "L" and lf == "2.0":
                    off["2.0, RT >= 2.6" if float(rt) >= 2.6 else "2.0, RT 2.2"].append(o)
            lines.append(f"- offsets from the step means, {'stored' if lab == 'b0' else 'refit'}: " + "; ".join(
                f"LF {k}: {np.mean(v):+.5f} (sd {np.std(v):.5f}, n {len(v)})" for k, v in sorted(off.items())))
        pc = {s: [] for s in ("LR", "R")}
        for d in rr:
            for a, v in d["per_cap"].items():
                if int(a) != d["hi"]:
                    for s in pc:
                        pc[s].append(v[s] - d["bf"])
        lines.append(f"- other captures' own optimum minus the refit (hi capture): L + R median {np.median(pc['LR']):+.5f}, "
                     f"MAD {np.median(np.abs(np.array(pc['LR']) - np.median(pc['LR']))):.5f}; R alone median {np.median(pc['R']):+.5f}")
    json.dump(out, open(os.path.join(OUT, "b_refit.json"), "w"), indent=1)
    open(os.path.join(OUT, "b_refit.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def interp_db(d, b):
    """First-pass null (dB) per capture and channel at b, cubic in the grid (exact at grid points)."""
    from scipy.interpolate import CubicSpline
    g, E = np.array(d["grid"]), np.array(d["E"])
    if not g[0] - 1e-9 <= b <= g[-1] + 1e-9:
        raise ValueError(f"{d['key']}: b {b} outside the grid {g[0]}..{g[-1]}")
    out = np.full(E.shape[1:], np.nan)
    ok = np.all(np.isfinite(E), axis=0)                         # unusable captures stay nan
    out[ok] = 10 * np.log10(CubicSpline(g, E[:, ok], axis=0)(b))
    return out                                                  # (ncap, 2)


def cand_stage():
    D = load()
    cands = {}
    for prog in ("M", "L"):
        rr = [d for (p, _), d in D.items() if p == prog]
        sc = [d for d in rr if not d.get("unscored")]  # step means from the settings with a usable first pass
        st = {s: np.mean([d["bf"] for d in sc if step_of(d["b0"]) == s]) for s in (1, 2, 3)}
        fold = {lf: {s: np.mean([d["bf"] for d in sc if step_of(d["b0"]) == s and split(d["key"])[0] != lf]) for s in (1, 2, 3)}
                for lf in {split(d["key"])[0] for d in rr}}
        cands.setdefault("refit", {})[prog] = {d["key"]: d["bf"] for d in rr}
        cands.setdefault("steps_all", {})[prog] = {d["key"]: st[step_of(d["b0"])] for d in rr}
        cands.setdefault("steps_cv", {})[prog] = {d["key"]: fold[split(d["key"])[0]][step_of(d["b0"])] for d in rr}
        cands.setdefault("steps_S", {})[prog] = {d["key"]: S_STEPS[step_of(d["b0"])] for d in rr}
        cands.setdefault("_steps", {})[prog] = dict(all=st, fold=fold)
    shared = {s: float(np.mean([S_STEPS[s]] + [cands["_steps"][p]["all"][s] for p in ("M", "L")])) for s in (1, 2, 3)}
    ml = {s: float(np.mean([d["bf"] for d in D.values() if step_of(d["b0"]) == s and not d.get("unscored")])) for s in (1, 2, 3)}
    for prog in ("M", "L"):
        rr = [d for (p, _), d in D.items() if p == prog]
        other = cands["_steps"]["L" if prog == "M" else "M"]["all"]
        cands.setdefault("steps_shared", {})[prog] = {d["key"]: shared[step_of(d["b0"])] for d in rr}
        cands.setdefault("steps_ml", {})[prog] = {d["key"]: ml[step_of(d["b0"])] for d in rr}
        cands.setdefault("steps_cross", {})[prog] = {d["key"]: other[step_of(d["b0"])] for d in rr}
    lines = [f"shared table (mean of S and the refit M / L step means): {[round(shared[s], 5) for s in (1, 2, 3)]}",
             f"M / L table (all 96 refits, by step): {[round(ml[s], 5) for s in (1, 2, 3)]}; steps_cross: each program on the other's step means"]
    for prog in ("M", "L"):
        lines.append(f"{prog} refit step means {[round(float(cands['_steps'][prog]['all'][s]), 5) for s in (1, 2, 3)]}; "
                     f"fold max |dev| {max(abs(cands['_steps'][prog]['fold'][lf][s] - cands['_steps'][prog]['all'][s]) for lf in cands['_steps'][prog]['fold'] for s in (1, 2, 3)):.5f}")
    res = {}
    for name in ("refit", "steps_all", "steps_cv", "steps_S", "steps_shared", "steps_ml", "steps_cross"):
        json.dump({p: {k: float(v) for k, v in cands[name][p].items()} for p in ("M", "L")},
                  open(os.path.join(OUT, f"b_cand_{name}.json"), "w"), indent=1)
        lines.append(f"\n## candidate {name}: first pass, damped captures, before (stored b) -> after; dB")
        for prog in ("M", "L"):
            rec = []
            for (p, key), d in D.items():
                if p != prog:
                    continue
                a0 = interp_db(d, d["b0"])
                a1 = interp_db(d, cands[name][prog][key])
                for i, amp in enumerate(d["amps"]):
                    if not np.all(np.isfinite(a0[i])):      # unusable capture (pulse_check scores none)
                        continue
                    rec.append(dict(key=key, amp=amp, hi=amp == d["hi"], before=a0[i].tolist(), after=a1[i].tolist()))
            res.setdefault(name, {})[prog] = rec
            for ch, chn in ((1, "R"), (0, "L")):
                for sel, f in (("all", lambda r: True), ("fit capture", lambda r: r["hi"]), ("held-out captures", lambda r: not r["hi"])):
                    x = [r for r in rec if f(r)]
                    b_, a_ = np.array([r["before"][ch] for r in x]), np.array([r["after"][ch] for r in x])
                    lines.append(f"- {prog} {chn} {sel} ({len(x)}): median {np.median(b_):.2f} -> {np.median(a_):.2f}, mean change "
                                 f"{np.mean(a_ - b_):+.3f}, worst capture {np.max(a_ - b_):+.2f}, better > 0.3 dB {int(np.sum(a_ - b_ < -0.3))}, worse > 0.3 dB {int(np.sum(a_ - b_ > 0.3))}")
                g = collections.defaultdict(list)
                for r in rec:
                    g[r["key"]].append(r["after"][ch] - r["before"][ch])
                gm = {k: float(np.mean(v)) for k, v in g.items()}
                bad = {k: round(v, 2) for k, v in sorted(gm.items()) if v > 0.5}
                lines.append(f"  - {prog} {chn} settings (mean over captures) worse > 0.5 dB: {len(bad)} {bad}; worst {max(gm.values()):+.2f} "
                             f"({max(gm, key=gm.get)}); best {min(gm.values()):+.2f}")
                for gname, gi in (("LF", 0), ("RT", 1), ("HF", 2)):
                    gg = collections.defaultdict(lambda: ([], []))
                    for r in rec:
                        v = split(r["key"])[gi]
                        gg[v][0].append(r["before"][ch])
                        gg[v][1].append(r["after"][ch])
                    lines.append(f"  - {prog} {chn} by {gname}: " + "; ".join(
                        f"{v}: {np.median(b_):.2f} -> {np.median(a_):.2f}" for v, (b_, a_) in sorted(gg.items(), key=lambda t: float(t[0]))))
    json.dump(res, open(os.path.join(OUT, "b_refit_cand.json"), "w"))
    open(os.path.join(OUT, "b_refit_cand.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["fit", "write", "cand"], required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", type=int, default=0)
    a = ap.parse_args()
    if a.stage == "fit":
        fit_stage(a.workers, a.only)
    elif a.stage == "write":
        write_stage()
    else:
        cand_stage()


if __name__ == "__main__":
    main()
