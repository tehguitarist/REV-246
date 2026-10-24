"""Two read-only diagnostics on one S setting, through the pulse check's render (analysis/pulse/s is not touched):

    bscan   the AP230 damping b scanned (the table's four steps and off-grid values), R nulls of the AMPL 0 and 90 captures
    resid   the clean AMPL 0 / 90 residuals (capture - model): energy re the signal per window, their correlation (a shared
            error is the model's; an unshared one is a capture's or its fit's), and the 0-60 ms per 5 ms bin

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 4 -- \\
        .venv/bin/python tools/loop_read/setting_diag.py resid 1.0_0.3_0.4 [more keys]
    ... setting_diag.py bscan 1.0_0.3_0.4 0.17,0.33,0.42,0.6      # also runs the stock b; NEXT_STEPS §0 item 5"""
import json, os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pulse_check as pc

q2, FS, PRE = pc.q2, pc.FS, pc.PRE


def setup(key):
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    g = info[key]
    fold = json.load(open(os.path.join(pc.CV, f"fold_{g['rt']}.json")))
    lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
    b0 = ap["steps"][ap["settings"][key]["step"]]
    files = [(a, nm) for a, nm in q2.files_of(g)]
    return dict(info=info, fade=fade, eps0=eps0, ends=ends, v=v, law=law, cells=cells, fold=fold, lab=lab, b0=b0, files=files, steps=ap["steps"])


def bscan(key, extra):
    s = setup(key)
    caps = [dict(amp=a, file=nm, eps=s["eps0"][nm], fade2=s["fade"][nm]) for a, nm in s["files"]]
    print(f"{key}: stock b {s['b0']:.4f}, steps {np.round(s['steps'], 4).tolist()}")
    for b in [s["b0"]] + extra:
        path = os.path.join(tempfile.mkdtemp(), "x.json")
        pc._s_job((key, s["lab"], b, s["law"], s["fold"], s["cells"].get(key), s["v"], caps, s["ends"], path))
        for r in json.load(open(path))["rows"]:
            if r["amp"] in (0, 90):
                print(f"  b {b:.4f} amp {r['amp']:2d}: " + " ".join(f"{w}:{x[1]:.1f}" for w, x in r["null"].items() if x))


def resid(key):
    s = setup(key)
    G = pc.setting_G(s["lab"], s["b0"], s["law"], s["fold"], s["cells"].get(key))
    R = {}
    for amp, nm in s["files"]:
        if amp in (0, 90):
            y, mdl, g, te = q2.fit_capture(G, s["v"], amp, nm, s["fade"][nm], s["eps0"][nm], s["ends"][nm], 1700)
            R[amp] = (pc.prep(y, g)["bb"], pc.prep(mdl * g)["bb"])
    print("##", key)
    sl = lambda w0, w1: slice(PRE + int(w0 * FS), PRE + int(w1 * FS))
    for w0, w1 in ((0, .05), (.05, .15), (.15, .30)):
        r = {a: R[a][0][sl(w0, w1)] - R[a][1][sl(w0, w1)] for a in R}
        e = {a: 10 * np.log10(np.sum(r[a][:, 1] ** 2) / np.sum(R[a][0][sl(w0, w1), 1] ** 2)) for a in R}
        rho = [np.corrcoef(r[0][:, c], r[90][:, c])[0, 1] for c in (0, 1)]
        print(f"  {int(w0 * 1000):3d}-{int(w1 * 1000):3d} ms: residual/signal R amp0 {e[0]:6.1f} amp90 {e[90]:6.1f} dB | correlation L {rho[0]:+.2f} R {rho[1]:+.2f}")
    print("  per 5 ms, 0-60 ms (R, residual / signal dB: amp0 | amp90; amp 0 has no cluster, so its first 15 ms are empty):")
    for k in range(12):
        w = sl(k * .005, (k + 1) * .005)
        e = {a: 10 * np.log10(np.sum((R[a][0][w, 1] - R[a][1][w, 1]) ** 2) / max(np.sum(R[a][0][w, 1] ** 2), 1e-30)) for a in R}
        print(f"   {k * 5:2d}-{k * 5 + 5:2d} ms: {e[0]:6.1f} | {e[90]:6.1f}")


if __name__ == "__main__":
    if sys.argv[1] == "bscan":
        bscan(sys.argv[2], [float(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [])
    else:
        for k in sys.argv[2:]:
            resid(k)
