"""A second capture defect at 99.9 ms file time: a waveform departure without a lag step (model_s_results 35.2).

    EMT_HFRT_MIN=0.12 EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 9 -- \\
        .venv/bin/python tools/loop_read/glitch_scan.py scan       # 9 min, 3.8 GB -> analysis/lowband/chain3_glitch_scan.{json,log}
    <same> .venv/bin/python tools/loop_read/glitch_scan.py pairs     # 1 min -> analysis/lowband/chain3_glitch_pairs.log

scan: every S capture (512) against the baseline model (pulse_check's render, PULSE_CV), residual / signal, R, 100 Hz-10 kHz,
at 88-96 and 98-106 ms after tap 0 (tap 0 is 2.67 ms into the file, so 97 ms = 99.9 ms file time, §28's slip time); a
capture whose residual jumps by more than 6 dB there is listed, with whether tools/defects.py already has it. §28's lag
detector (stepband.py) finds only lag steps; this finds any departure from the model at that time.
pairs: at RT 0.3, the HF 0.2 and HF 0.4 captures of each LF step (the same filter in the unit, model_s_results 29, and the
same model) against each other and against the model (50-150 / 150-300 ms), and the departure's onset per 1 ms for the
suspect captures against a sibling capture of the same filter."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import sosfiltfilt
import pulse_check as pc
from defects import DEFECTS, BAD, SLIPPED

q2, q1 = pc.q2, pc.q2.q1
FS, PRE = pc.FS, pc.PRE
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LOG = os.path.join(ROOT, "analysis", "lowband", "chain3_glitch_")


def ctx():
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    return dict(info=info, fade=fade, eps0=eps0, ap=ap, ends=ends, v=v, law=law, cells=cells)


def G_of(c, key):
    g = c["info"][key]
    lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
    fold = json.load(open(os.path.join(pc.CV, f"fold_{g['rt']}.json")))
    return pc.setting_G(lab, c["ap"]["steps"][c["ap"]["settings"][key]["step"]], c["law"], fold, c["cells"].get(key))


def cap(c, G, key, a):
    nm = dict(q2.files_of(c["info"][key]))[a]
    y, mdl, g, te = q2.fit_capture(G, c["v"], a, nm, c["fade"][nm], c["eps0"][nm], c["ends"][nm], 1700)
    return y, mdl * g, te, nm


def rel(x, y, t0, t1):
    i0, i1 = PRE + int(t0 * FS), PRE + int(t1 * FS)
    return float(10 * np.log10(np.sum((x[i0:i1] - y[i0:i1]) ** 2) / np.sum(x[i0:i1] ** 2)))


def scan_job(key):
    c = ctx()
    G = G_of(c, key)
    out = []
    for a, nm in q2.files_of(c["info"][key]):
        y, m, te, _ = cap(c, G, key, a)
        yb, mb = sosfiltfilt(q1.ms.BAND, y[:, 1]), sosfiltfilt(q1.ms.BAND, m[:, 1])
        out.append(dict(key=key, amp=a, file=nm, before=rel(yb, mb, 0.088, 0.096), after=rel(yb, mb, 0.098, 0.106), defect=(key, a) in DEFECTS))
    return out


def report(res):
    jump = sorted([r for r in res if r["after"] - r["before"] > 6], key=lambda r: r["before"] - r["after"])
    d = np.array([r["after"] - r["before"] for r in res if r["after"] - r["before"] <= 6])
    lines = [f"# model {os.path.relpath(pc.CV, ROOT)}; {len(res)} captures; residual / signal (R, 100 Hz-10 kHz) at 88-96 -> 98-106 ms after tap 0",
             f"jump > 6 dB ({len(jump)}):"]
    lines += [f"  {r['key']:12s} AMPL {r['amp']:2d} {r['file']:28s} {r['before']:6.1f} -> {r['after']:6.1f}   {'in DEFECTS' if r['defect'] else 'NOT LISTED'}" for r in jump]
    lines.append(f"the other {len(d)} captures: jump median {np.median(d):+.2f} dB, 99th percentile {np.percentile(d, 99):+.2f}, max {d.max():+.2f}")
    return "\n".join(lines)


def pairs():
    c = ctx()
    W = [(50, 150), (150, 300)]
    nl = lambda y, m, te: " / ".join(f"{x[1]:.1f}" if x else "-" for x in q1.nulls(y, m, te, W, narrow=True).values())
    L = ["# RT 0.3: HF 0.2 and HF 0.4 captures against the shared model and against each other; R, 100 Hz-10 kHz, 50-150 / 150-300 ms"]
    for lf in (0.1, 0.4, 1.0, 2.0):
        k2, k4 = f"{lf}_0.3_0.2", f"{lf}_0.3_0.4"
        G = G_of(c, k4)
        Y = {(k, a): cap(c, G, k, a) for k in (k2, k4) for a, _ in q2.files_of(c["info"][k]) if (k, a) not in BAD | SLIPPED}   # the GLITCH files stay in: they are what pairs examines
        L.append(f"## LF {lf}")
        L += [f"   {k} AMPL {a:2d} vs model: {nl(y, m, te)}" for (k, a), (y, m, te, _) in Y.items()]
        L += [f"   {k2} vs {k4}, AMPL {a:2d}: {nl(Y[(k2, a)][0], Y[(k4, a)][0], min(Y[(k2, a)][2], Y[(k4, a)][2]))}"
              for a in (0, 30, 60, 90) if (k2, a) in Y and (k4, a) in Y]
        L += [f"   in-take repeat {k} AMPL 0 vs 90: {nl(Y[(k, 0)][0], Y[(k, 90)][0], min(Y[(k, 0)][2], Y[(k, 90)][2]))}"
              for k in (k2, k4) if (k, 0) in Y and (k, 90) in Y]
    L.append("\n# onset: suspect capture against a sibling of the same filter, residual / signal per 1 ms after tap 0 (dB), R")
    for mk, ref, sus in [("1.0_0.3_0.2", ("1.0_0.3_0.2", 0), ("1.0_0.3_0.4", 0)), ("1.0_0.3_0.2", ("1.0_0.3_0.2", 90), ("1.0_0.3_0.4", 90)),
                         ("2.0_0.3_0.4", ("2.0_0.3_0.4", 60), ("2.0_0.3_0.4", 30)), ("2.0_0.8_0.2", ("2.0_0.8_0.2", 60), ("2.0_0.8_0.2", 30)),
                         ("2.0_0.3_0.4", ("2.0_0.3_0.4", 60), ("2.0_0.3_0.4", 90))]:
        G = G_of(c, mk)
        a = sosfiltfilt(q1.ms.BAND, cap(c, G, *ref)[0][:, 1])
        b = sosfiltfilt(q1.ms.BAND, cap(c, G, *sus)[0][:, 1])
        L.append(f"  {sus[0]} AMPL {sus[1]} vs AMPL {ref[1]} of {ref[0]}: " + " ".join(f"{k}:{rel(a, b, k / 1e3, (k + 1) / 1e3):.0f}" for k in range(88, 106)))
    return "\n".join(L)


if __name__ == "__main__":
    if sys.argv[1] == "scan":
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(8) as ex:
            res = [r for rr in ex.map(scan_job, sorted(q2.q1_tables()[0])) for r in rr]
        json.dump(res, open(LOG + "scan.json", "w"), indent=0)
        txt = report(res)
        open(LOG + "scan.log", "w").write(txt + "\n")
    elif sys.argv[1] == "report":
        txt = report(json.load(open(LOG + "scan.json")))
        open(LOG + "scan.log", "w").write(txt + "\n")
    else:
        txt = pairs()
        open(LOG + "pairs.log", "w").write(txt + "\n")
    print(txt)
