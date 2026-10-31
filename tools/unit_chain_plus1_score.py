"""§37.2: fullband scoring of the "CLA 10 + 1" candidates, (a) refitted per variant as unit_transfer.fullband_stage and (b) FROZEN:
each capture's clock, shift, gain and 14-number FIR as fitted for the shipped CLA 10, applied unchanged to the candidate's model
(so the refit's own movement is removed). Same IRs, model and bands as fullband_stage (copied loop). 11-12.3 kHz is skipped.
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/unit_chain_plus1_score.py
 -> analysis/unit_chain/plus1_score.json ; print with --report"""
import glob, json, os, sys
import numpy as np
from scipy.signal import sosfiltfilt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "harness"))
import unit_transfer as ut  # noqa: E402
import unit_chain_plus1 as p1  # noqa: E402

ROOT, OUT = p1.ROOT, os.path.join(p1.ROOT, "analysis", "unit_chain")
P = json.load(open(os.path.join(OUT, "plus1.json")))
VARS = ["plain", "bq_CLA"] + list(P)
pc, q1, q2, FS, PRE = ut.pc, ut.q1, ut.q2, ut.FS, ut.PRE


def H_of(vn):
    if vn == "plain":
        return None
    f = (lambda f: p1.ut.bq_resp(p1.X10, f)) if vn == "bq_CLA" else (lambda f, x=P[vn]["x"]: p1.resp(x, f))
    return lambda ff: np.repeat(f(ff)[:, None], 2, 1)


def run():
    import altiverb_predict as av
    import sml_fullband as fb
    q1.EXT = 1700
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    fold = json.load(open(os.path.join(ROOT, "analysis", "loop_cv", "fold_all.json")))
    fold["pramp"] = [float(np.median([json.load(open(p))["pramp"][0] for p in
                                      glob.glob(os.path.join(ROOT, "analysis", "loop_cv", f"cap{os.environ['EMT_PMAX']}_folds", "fold_*.json"))])), 0.62]
    irs = [("CLA", 1.0, av.cla_h())] + [(f"Altiverb RT {k}", rt, av.altiverb_h(k)[0]) for k, rt in av.PRESETS.items()]
    bands = [b for b in fb.BNAME if b != fb.BNAME[-1]]
    res = {"bands": bands, "wins": [f"{a}-{b}" for a, b in fb.WINS], "irs": {}}
    for name, rt, h in irs:
        G = pc.setting_G(dict(lf=1.0, rt=rt, hf=0.9), ap["steps"][0], law, fold, None)
        hs = h * (np.sign(h[PRE].sum()) or 1.0)
        yb = sosfiltfilt(q1.LFBAND, hs, axis=0)
        n = min(len(yb) - 600, G["N"])
        m0 = q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=80)])]).render(v, G, n, 80)
        res["irs"][name] = {"rt": rt}
        frozen = None
        for vn in VARS:                                    # bq_CLA second in VARS... frozen needs it first
            pass
        order = ["bq_CLA"] + [v_ for v_ in VARS if v_ != "bq_CLA"]
        for vn in order:
            H = H_of(vn)
            mdl = m0 if H is None else ut.through(m0, H)
            eg, shg, gain = av.fit_clock(yb, mdl, starts=[98.4] if name == "CLA" else [42.8])
            y = q1.apply_clock(yb, eg, shg)[:n] / gain
            mdl2, fir = av.channel_fir(y, mdl)
            eg2, shg2, gain2 = av.fit_clock(yb, mdl2, starts=[eg])
            y2 = q1.apply_clock(yb, eg2, shg2)[:n] / gain2
            te = min(1.5, q2.comp_end(y2, np.ones_like(y2)))
            if vn == "bq_CLA":
                frozen = (eg2, shg2, gain2, fir, y2, te)
            _, _, _, fr_fir, fr_y, fr_te = frozen
            K = 3
            mf = np.zeros_like(mdl)                         # the candidate through the shipped CLA 10's FIR, clock, shift, gain
            for c in (0, 1):
                mf[:, c] = np.stack([np.roll(mdl[:, c], k) for k in range(-K, K + 1)], 1) @ fr_fir[c]
            res["irs"][name][vn] = dict(refit=ut.band_nulls(y2, mdl2, te, bands), frozen=ut.band_nulls(fr_y, mf, fr_te, bands))
            print(name, vn, "50-150 R 270-400 refit/frozen", round(res["irs"][name][vn]["refit"]["270-400"]["50-150"][1], 1),
                  round(res["irs"][name][vn]["frozen"]["270-400"]["50-150"][1], 1), flush=True)
    json.dump(res, open(os.path.join(OUT, "plus1_score.json"), "w"))


def report():
    R = json.load(open(os.path.join(OUT, "plus1_score.json")))
    W = ["50-150", "150-300", "300-600", "600-1000"]
    for mode in ("frozen", "refit"):
        for sn, us in (("CLA", ["CLA"]), ("Altiverb RT 1.6", ["Altiverb RT 1.6"]), ("Altiverb other 4 (held out of P_ho), median", [u for u in R["irs"] if u.startswith("Altiverb") and u != "Altiverb RT 1.6"])):
            print(f"## {mode}: {sn}: null R minus plain (dB, mean of windows 50-1000 ms; negative = better) | candidate minus shipped CLA 10 (+ = worse)")
            cs = [v for v in VARS if v not in ("plain", "bq_CLA")]
            print(f"{'band':12s}{'CLA10':>8s}" + "".join(f"{c:>8s}" for c in cs) + "  |" + "".join(f"{c:>8s}" for c in cs))
            for b in R["bands"]:
                def d(vn):
                    return np.median([np.mean([R["irs"][u][vn][mode][b][w][1] - R["irs"][u]["plain"][mode][b][w][1] for w in W if R["irs"][u][vn][mode][b][w]]) for u in us])
                a = d("bq_CLA")
                print(f"{b:12s}{a:+8.1f}" + "".join(f"{d(c):+8.1f}" for c in cs) + "  |" + "".join(f"{d(c) - a:+8.2f}" for c in cs))
            print()


if __name__ == "__main__":
    report() if "--report" in sys.argv else run()
