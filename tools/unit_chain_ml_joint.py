"""§37.x: the M / L first pass on the Altiverb MED / LONG presets (unit_chain.ml_stage's protocol, copied; only the
variants differ) through CLA's, the Altiverb's and the joint low-band cascades (10 and 13 numbers).
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 3 -- \
        .venv/bin/python tools/unit_chain_ml_joint.py        # -> analysis/unit_chain/ml_joint.{json,log}
Run from unit_chain.py's environment (same as `--stage ml`)."""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import unit_chain as uc  # noqa: E402

pc, FS, ROOT, OUT = uc.pc, uc.FS, uc.ROOT, uc.OUT


def main():
    import altiverb_ml as am
    import ml_lines as mll
    S = json.load(open(os.path.join(OUT, "s.json")))["irs"]
    bq = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))
    J = json.load(open(os.path.join(OUT, "joint.json")))
    variants = {"none": None, "CLA10": bq["CLA"]["x"], "Alt10": bq["Altiverb RT 1.6"]["x"],
                "CLA13": S["CLA"]["chain@150-300"]["chain"]["bq"], "Alt13": S["Altiverb RT 1.6"]["chain@150-300"]["chain"]["bq"],
                **{k: v["x"] for k, v in J.items()}}
    if os.environ.get("EMT_ML_PLUS1"):                 # §37.2: CLA's 10 numbers + one peak (lb_resp's 13-number form)
        P = json.load(open(os.path.join(OUT, "plus1.json")))
        variants = {"none": None, "CLA10": bq["CLA"]["x"], **{k: bq["CLA"]["x"] + v["x"] for k, v in P.items()}}
    K0, res = am.K0, []
    for prog, (pre, lf, hf) in am.SET.items():
        names = sorted((os.path.basename(p)[:-4] for p in glob.glob(os.path.join(am.CACHE, pre + " *.npy"))), key=lambda t: float(t.split(" ")[-1][:-1]))
        for name in names:
            rt = float(name.split(" ")[-1][:-1])
            n = K0 + int(min(1.5, 1.2 * rt + 0.5) * FS)
            y = am.device_ir(name, n)
            i1 = min(am.usable(y), n - 1)
            key = f"{lf}_{rt}_{hf}"
            M = mll.Model(prog)
            row = dict(prog=prog, preset=name, rt=rt)
            fb = None
            for b in am.B_STEPS:
                M.b[key] = b
                C, Rr = M.parts(n, key)
                fpn = pc._ml_nulls(prog, y, am.fit(prog, y, C, Rr, i1), i1)["first_pass"]
                if fb is None or fpn < fb[1]:
                    fb = (b, fpn, C, Rr)
            b, _, C0, R0 = fb
            best = None
            for d in np.arange(-6.0, 6.01, 1.0):
                C, Rr = (mll.warp(C0, 0.0, d), mll.warp(R0, 0.0, d)) if d else (C0, R0)
                nl = pc._ml_nulls(prog, y, am.fit(prog, y, C, Rr, i1), i1)
                if best is None or (nl["recirc1"] or 0) < (best[1]["recirc1"] or 0):
                    best = (d, nl, C, Rr)
            d, _, C, Rr = best
            row.update(b=b)
            W = mll.windows(prog)[1]
            for vn, x in variants.items():
                mdl = am.fit(prog, y, uc.lowband(C, x), uc.lowband(Rr, x), i1)
                row[vn] = dict(null=dict(cluster_gap=uc.cluster_gap(prog, y, mdl), **pc._ml_nulls(prog, y, mdl, i1)), fp_oct=uc.fp_oct(y, mdl, K0, W))
            res.append(row)
            print(prog, name, " | ".join(f"{vn} {row[vn]['null']['first_pass']:.1f}" for vn in variants), flush=True)
    json.dump(res, open(os.path.join(OUT, ("ml_plus1.json" if os.environ.get("EMT_ML_PLUS1") else "ml_joint.json")), "w"))
    lines = []
    for prog in ("M", "L"):
        rr = [r for r in res if r["prog"] == prog and r["b"] == 0.0]
        lines.append(f"{prog} median over {len(rr)} b=0 presets, R: first pass [cluster+gap] (recirc1)")
        for vn in variants:
            lines.append(f"  {vn:10s} fp {np.median([r[vn]['null']['first_pass'] for r in rr]):.1f} cl {np.median([r[vn]['null']['cluster_gap'] for r in rr]):.1f} r1 {np.median([r[vn]['null']['recirc1'] for r in rr]):.1f}"
                         + "   oct " + " ".join(f"{f}:{np.median([r[vn]['fp_oct'][str(f)] for r in rr]):.1f}" for f in uc.OCT))
    open(os.path.join(OUT, ("ml_plus1.log" if os.environ.get("EMT_ML_PLUS1") else "ml_joint.log")), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
