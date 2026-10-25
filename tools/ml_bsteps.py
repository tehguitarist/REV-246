"""M / L's last all-pass damping b as three step values per program (ml_structure.md section 19, NEXT_STEPS audit item 6).

The 96 damped settings (M 56, L 40) carried each its own fitted b (firstpass.json), in three clusters (0.1707-0.1711 /
0.3303-0.3308 / 0.5000-0.5006 in M, 0.172-0.173 / 0.332 / 0.502 in L). Like S's AP230 (model_s_results 14.3) the step is
set by the RT and HF controls, not LF: every (RT, HF) cell has the same step in all four LF columns (checked here). This
writes the cell -> step table and, per program, the step means: from all damped settings ("steps"), and leave-one-LF-step-out
("fold": the means of the other three LF columns, used for the settings of the held-out LF) for the held-out score.

    .venv/bin/python tools/ml_bsteps.py        -> analysis/ml_lines/b_steps.json   (reads analysis/ml_carryover/firstpass.json)
    .venv/bin/python tools/ml_bsteps.py score BASE_PULSE_DIR CAND_PULSE_DIR [OUT.json]
        before -> after of two `pulse_check --stage ml` runs (their <dir>/ml/*.json) on the damped settings: first pass and
        recirculation 1 per capture, per step and per group, with the groups that regress by > 0.3 dB listed

ml_lines.Model reads it with ML_BSTEP=1 (all-data means) or ML_BSTEP=cv (leave-LF-out means); unset: the stored per-setting b.
"""
import collections, json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FP = os.path.join(ROOT, "analysis", "ml_carryover", "firstpass.json")
OUT = os.path.join(ROOT, "analysis", "ml_lines", "b_steps.json")


def step_of(b):
    return 0 if b < 0.05 else (1 if b < 0.25 else 2 if b < 0.42 else 3)


def main():
    fp = json.load(open(FP))
    res = {}
    for prog in ("M", "L"):
        b = {}
        for r in fp[prog]["rows"]:
            b[r["key"]] = r["b"]
        cells, S = collections.defaultdict(set), collections.defaultdict(list)
        for k, v in b.items():
            lf, rt, hf = k.split("_")
            cells[f"{rt}_{hf}"].add(step_of(v))
            S[step_of(v)].append((lf, v))
        assert all(len(s) == 1 for s in cells.values()), "a cell's step depends on LF"
        steps = [0.0] + [float(np.mean([v for _, v in S[s]])) for s in (1, 2, 3)]
        lfs = sorted({lf for s in (1, 2, 3) for lf, _ in S[s]}, key=float)
        fold = {lf: [0.0] + [float(np.mean([v for l, v in S[s] if l != lf])) for s in (1, 2, 3)] for lf in lfs}
        dmax = max(abs(f[s] - steps[s]) for f in fold.values() for s in (1, 2, 3))
        res[prog] = dict(steps=steps, fold=fold, cells={c: int(next(iter(s))) for c, s in sorted(cells.items())},
                         n={str(s): len(S[s]) for s in (1, 2, 3)}, max_fold_dev=dmax,
                         max_dev_stored=max(abs(v - steps[s]) for s in (1, 2, 3) for _, v in S[s]))
        print(prog, "steps", [round(x, 5) for x in steps], "n", res[prog]["n"], f"max |fold - all| {dmax:.5f}",
              f"max |stored b - step| {res[prog]['max_dev_stored']:.5f}")
    json.dump(res, open(OUT, "w"), indent=1)


def score(base, cand, out=None):
    import glob
    cells = json.load(open(OUT))
    rec = []
    for prog in ("M", "L"):
        for p in sorted(glob.glob(os.path.join(cand, "ml", f"{prog}_*.json"))):
            A, B = json.load(open(os.path.join(base, "ml", os.path.basename(p)))), json.load(open(p))
            lf, rt, hf = B["key"].split("_")
            ra = {r["amp"]: r["null"] for r in A["rows"]}
            for r in B["rows"]:
                a, b = ra[r["amp"]], r["null"]
                rec.append(dict(prog=prog, key=B["key"], amp=r["amp"], step=cells[prog]["cells"][f"{rt}_{hf}"],
                                first_pass=[a["first_pass"], b["first_pass"]], recirc1=[a["recirc1"], b["recirc1"]]))
    for prog in ("M", "L"):
        for st in (0, 1, 2, 3):
            x = [r for r in rec if r["prog"] == prog and r["step"] == st and r["first_pass"][0] is not None]
            if not x:
                continue
            a, b = np.array([r["first_pass"] for r in x]).T
            print(f"{prog} step {st}: {len(x)} captures, first pass median {np.median(a):.2f} -> {np.median(b):.2f}, mean change {np.mean(b - a):+.3f}, "
                  f"worse by > 0.3 dB: {int(np.sum(b - a > 0.3))}, better: {int(np.sum(b - a < -0.3))}, max worse {np.max(b - a):+.2f}")
        g = collections.defaultdict(list)
        for r in rec:
            if r["prog"] == prog and r["step"] > 0 and r["first_pass"][0] is not None:
                g[r["key"]].append(r["first_pass"][1] - r["first_pass"][0])
        print(f"{prog} groups regressing > 0.3 dB (mean over captures):", {k: round(float(np.mean(v)), 2) for k, v in sorted(g.items()) if np.mean(v) > 0.3})
    if out:
        json.dump(rec, open(out, "w"))


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "score":
        score(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    else:
        main()
