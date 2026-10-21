"""NEXT_STEPS Q2: the per-loop injection/read-out structure of program S (Tomkin 96 kHz library).
Results: docs/model_s_results.md §15.

    .venv/bin/python tools/q2_injection.py [--stage windows|free|joint|report] [--workers 10] [--only keys] [--out analysis/q2]

Builds on tools/q1_generalise.py (its analysis/q1 cache, laws, AP230 table and clock fits).

Stages:
  windows  corrected window ends (fade2.json): q1's fade.json lets files shorter than about 1.5 s run into Tomkin's
           fade (its guard is min(0.70 s, 45 % of the file); the fade starts 0.65 s before the last non-zero sample
           in every file, short ones included). Here the end is min(q1's end, last non-zero sample - FADE_GUARD).
  free     per setting, a free injection of its own (no sharing across settings), loops from the laws:
           does each setting fit to the repeat-capture level when its injection is free? Three tests per setting:
             in-sample   fit on all its captures, up to each fade;
             capture     fit on AMPL 30 + 90, null the AMPL 0 and 60 captures (other captures, same reverb);
             time        fit on 0-T_FIT only, predict T_FIT to the fade (the loops and diffuser carry it).
           If these reach the repeat-capture null, the loops and diffuser are right and Q2 is only how the
           injection depends on the controls; if not, the recursion itself is wrong where it fails.
  ramp     Tomkin's fade-out, for compensation only (never part of the reverb model): every file ends with a linear
           amplitude ramp of FADE_L = 0.700 s (127 file/partner pairs agree within 0.1 dB; per-file fits 0.699 s,
           IQR 0.693-0.703). ramp.json: each raw file's end mapped onto q1's device grid per channel. With L fixed, the
           ramp ends FADE_END0 = -0.65 ms from the raw end in every file (IQR -0.90..-0.36 ms, from the gains against an
           unfaded AMPL partner), so the compensation has no free parameter and does not depend on the model.
  nulls    re-evaluate saved models (free per setting, joint per split) with the fade compensated: clock and gain fitted
           before the fade as before, then nulls of the data against model x ramp (the model itself stays unfaded),
           windows running on to where the compensated level reaches FLOOR_DB below 50-300 ms or the ramp falls below
           RAMP_MIN. The repeat-capture null is compensated the same way (y_a g_b against y_b g_a).
  joint    q1's held-out bench (run_split: shared injection x0 + G.x1, laws, AP230 table) re-run with the
           corrected windows, so Q2 candidates are compared against a fade-free baseline. --only limits it to a
           few settings for a smoke test (training = the --only keys not in the split's test set).
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import q1_generalise as q1  # noqa: E402

FS, PRE = q1.FS, q1.PRE
Q1 = os.path.join(ROOT, "analysis", "q1")
T_FIT = 0.25                    # time test: fit window 0-T_FIT s after tap 0
WINDOWS = [(0, 50), (50, 150), (150, 300), (300, 600), (600, 1000), (1000, 1500)]
OCT = [125, 250, 500, 1000, 2000, 4000, 8000]
OCT_WIN = [(50, 150), (150, 300), (300, 600)]
OCT_SOS = {f: butter(3, [f / np.sqrt(2), min(f * np.sqrt(2), 11000)], btype="band", fs=q1.FS, output="sos") for f in OCT}


def octave_nulls(y, mdl, t_end):
    """Null per octave band and window (R channel), each band relative to its own energy; None past t_end."""
    out = {}
    for f, sos in OCT_SOS.items():
        yb, mb = sosfiltfilt(sos, y[:, 1]), sosfiltfilt(sos, mdl[:, 1])
        row = []
        for a, b in OCT_WIN:
            i0, i1 = PRE + int(a / 1000 * FS), min(PRE + int(b / 1000 * FS), PRE + int(t_end * FS), len(y))
            row.append(None if i1 - i0 < int(0.05 * FS) else
                       float(10 * np.log10(np.sum((yb[i0:i1] - mb[i0:i1]) ** 2) / np.sum(yb[i0:i1] ** 2))))
        out[f] = row
    return out


FADE_GUARD = 0.70               # Tomkin's fades start 0.65 s before the last non-zero sample (q1 FADE_GUARD)


def windows_stage():
    """fade2.json: per file min(q1 window end, last non-zero sample - FADE_GUARD); 0 where nothing is left."""
    fade = json.load(open(os.path.join(Q1, "fade.json")))
    out = {}
    for nm, te in fade.items():
        h = q1.load(Q1, nm)
        last = (np.nonzero(np.any(h != 0, 1))[0][-1] - PRE) / FS
        out[nm] = float(max(min(te, last - FADE_GUARD), 0.0))
    json.dump(out, open(os.path.join(OUT, "fade2.json"), "w"), indent=0)
    d = np.array([fade[k] - out[k] for k in out])
    print(f"fade2: {np.sum(d > 0.005)} of {len(d)} windows shortened, by median {np.median(d[d > 0.005]):.3f} s (max {d.max():.3f})")
    return out


FADE_L = 0.700                  # Tomkin's fade: linear amplitude ramp over the last 0.700 s
FADE_END0 = -0.65e-3            # ramp end re the raw file end (all files, L fixed)
RAMP_MIN = 0.1                  # compensated windows stop where the ramp is below this (-20 dB)
FLOOR_DB = 60.0                 # ... or where the compensated level is this far below the 50-300 ms level (24-bit floor)


def raw_end_dev(nm):
    """q1 device index (L, R) of the raw 96 kHz file's end, through ml_structure's device grid and q1's shift_LR."""
    import soundfile as sf
    import ml_structure as ms
    from core import PAD
    info = json.load(open(os.path.join(Q1, "data.json")))
    x, fs = sf.read(next(p for p in ms.files("S") if os.path.basename(p) == nm), always_2d=True)
    a = np.abs(x).max(1)
    main = int(np.argmax(a > a.max() * 0.1))
    lo = max(main - int(0.002 * fs), 0)
    onset = lo + int(np.argmax(a[lo:main + 1] > a.max() * 1e-2))       # as ml_structure._make
    anchor = float(onset + PAD)
    k0 = int(np.floor(anchor / fs * ms.FSD))
    sh = next(f["shift_LR"] for g in info.values() for k, f in g["files"].items() if not k.startswith("_") and f["file"] == nm)
    return [float(k0 + (len(x) + PAD - anc) * ms.FSD / fs - max(k0 - ms.K0, 0) + sh[c])
            for c, anc in ((0, anchor + ms.LR_US * 1e-6 * fs), (1, anchor))]


def ramp_stage(workers):
    from concurrent.futures import ProcessPoolExecutor
    names = sorted(json.load(open(os.path.join(OUT, "fade2.json"))))
    with ProcessPoolExecutor(workers) as ex:
        ends = list(ex.map(raw_end_dev, names, chunksize=8))
    json.dump({nm: e for nm, e in zip(names, ends)}, open(os.path.join(OUT, "ramp.json"), "w"), indent=0)
    print(f"ramp: {len(names)} files mapped")


def ramp_gain(n, end_dev):
    """(n, 2) fade gain: 1 before the ramp, linear to 0 at end_dev + FADE_END0 (per channel, device samples)."""
    k = np.arange(n)[:, None]
    e = np.array(end_dev)[None] + FADE_END0 * FS
    return np.clip((e - k) / (FADE_L * FS), 0.0, 1.0)


def comp_end(y, g):
    """End (s after tap 0) of the compensated window: ramp >= RAMP_MIN and level within FLOOR_DB of 50-300 ms."""
    nf = int(0.025 * FS)
    lv = np.array([np.mean((y[j:j + nf] / np.maximum(g[j:j + nf], 1e-3)) ** 2) for j in range(PRE, len(y) - nf, nf)])
    lv = 10 * np.log10(lv / lv[2:12].mean() + 1e-30)
    k = next((k for k in range(2, len(lv)) if lv[k] < -FLOOR_DB), len(lv))
    j = int(np.argmax(np.min(g, 1) < RAMP_MIN)) if np.any(np.min(g, 1) < RAMP_MIN) else len(g)
    return min(k * 0.025, (j - PRE) / FS)


def q1_tables():
    info = json.load(open(os.path.join(Q1, "data.json")))
    fade = json.load(open(os.path.join(OUT, "fade2.json")))
    fits = json.load(open(os.path.join(Q1, "loops.json")))
    law = json.load(open(os.path.join(Q1, "laws.json")))
    law["lf"] = {float(k): v for k, v in law["lf"].items()}
    ap = json.load(open(os.path.join(Q1, "ap230.json")))
    eps0, _ = q1.file_eps(fits)
    return info, fade, law, ap, eps0


def files_of(g):
    return [(int(k), f["file"]) for k, f in sorted(g["files"].items()) if not k.startswith("_")]


def _band(name):
    return sosfiltfilt(q1.LFBAND, q1.load(Q1, name), axis=0)


def fit_setting(args):
    """Free injection for one setting (variant x: excitations P_i + EXT, first-pass FIR, cluster gains)."""
    key, g, fade, law, b230, eps0, iters = args
    lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
    G = q1.group_spectra(lab, law, q1.grid_n(g["rt"]), b230=b230)
    fl = files_of(g)
    eps = {nm: eps0[nm] for _, nm in fl}
    shift = {nm: 0.0 for _, nm in fl}

    def build(use, t_end):
        files = []
        for amp, nm in fl:
            if amp in use:
                y, w = q1.prep(Q1, nm, g["rt"], eps[nm], min(t_end, fade[nm]), shift=shift[nm])
                files.append(dict(y=y, w=w, amp=amp, name=nm))
        return q1.Joint([dict(G, files=files)])

    def evaluate(J, v, test, keep=None):
        out = {}
        for amp, nm in fl:
            if amp not in test:
                continue
            yb = _band(nm)
            n = min(PRE + int(fade[nm] * FS), len(yb))
            mdl = J.render(v, G, n, amp)
            _, (e, sh, gain) = q1._clock_job((yb, mdl, fade[nm], eps[nm], amp == 0))
            y = q1.apply_clock(yb, e, sh)[:n] / gain
            if keep is not None:
                keep[amp] = y
            out[amp] = dict(eps=e, gain_db=float(20 * np.log10(gain)),
                            null=q1.nulls(y, mdl, fade[nm], WINDOWS, narrow=True), octave=octave_nulls(y, mdl, fade[nm]))
        return out

    def ceiling(ys):
        """Repeat-capture null per window: AMPL partners against each other, each re-gridded and scaled by
        its own clock and gain fit against the in-sample model (as q1 run_split)."""
        out = {}
        amps = {a: nm for a, nm in fl}
        for a_, b_ in ((30, 60), (60, 90), (90, 30), (0, 30)):
            if a_ in ys and b_ in ys:
                n = min(len(ys[a_]), len(ys[b_]))
                te = min(fade[amps[a_]], fade[amps[b_]])
                out[f"{a_}-{b_}"] = dict(q1.nulls(ys[a_][:n], ys[b_][:n], te, WINDOWS[1:], narrow=True),
                                         octave=octave_nulls(ys[a_][:n], ys[b_][:n], te))
        return out

    # clock and AMPL-0 shift: one refit round against this setting's own in-sample model
    J = build({0, 30, 60, 90}, 9.0)
    v = J.solve(iters)
    for amp, nm in fl:
        yb = _band(nm)
        n = min(PRE + int(fade[nm] * FS), len(yb))
        (e, sh), _ = q1._clock_job((yb, J.render(v, G, n, amp), fade[nm], eps[nm], amp == 0))
        eps[nm], shift[nm] = e, sh
    res = dict(key=key, **lab, b230=b230, fade={nm: fade[nm] for _, nm in fl})
    J = build({0, 30, 60, 90}, 9.0)
    v = J.solve(iters, x0=v)
    ys = {}
    res["in_sample"] = evaluate(J, v, {0, 30, 60, 90}, keep=ys)
    res["repeat"] = ceiling(ys)
    np.save(os.path.join(OUT, "v", key + ".npy"), v.astype(np.float32))
    J = build({30, 90}, 9.0)
    res["capture"] = evaluate(J, J.solve(iters, x0=v), {0, 60})
    J = build({0, 30, 60, 90}, T_FIT)
    res["time"] = evaluate(J, J.solve(iters, x0=v), {0, 30, 60, 90})
    res["eps"] = eps
    res["shift"] = shift
    print(f"{key}: " + "  ".join(f"{t} " + " ".join(
        f"{w}:{'--' if x is None else '%.1f/%.1f' % tuple(x)}" for w, x in res[t][min(res[t])]["null"].items() if x is not None)
        for t in ("in_sample", "capture", "time")), flush=True)
    return res


def free_stage(workers, iters, only=None):
    from concurrent.futures import ProcessPoolExecutor
    info, fade, law, ap, eps0 = q1_tables()
    os.makedirs(os.path.join(OUT, "v"), exist_ok=True)
    keys = [k for k in info if only is None or k in only]
    jobs = [(k, info[k], fade, law, ap["steps"][ap["settings"][k]["step"]], eps0, iters) for k in keys]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(fit_setting, jobs))
    path = os.path.join(OUT, "free.json")
    old = {r["key"]: r for r in json.load(open(path))} if os.path.exists(path) else {}
    old.update({r["key"]: r for r in res})
    json.dump(list(old.values()), open(path, "w"), indent=0)
    return res


def _rows(test, w, c, sel=lambda r: True):
    out = []
    for r in json.load(open(os.path.join(OUT, "free.json"))):
        if not sel(r):
            continue
        for amp, x in r[test].items():
            x = x if test == "repeat" else x["null"]
            if x.get(w) is not None:
                out.append(x[w][c])
    return out


def _excess(test, w, c, sel=lambda r: True):
    """Per setting: median null of `test` minus the median repeat-capture null (same window); one value per setting."""
    out = []
    for r in json.load(open(os.path.join(OUT, "free.json"))):
        if not sel(r) or "repeat" not in r:
            continue
        a = [x["null"][w][c] for x in r[test].values() if x["null"].get(w) is not None]
        b = [x[w][c] for x in r["repeat"].values() if x.get(w) is not None]
        if a and b:
            out.append(np.median(a) - np.median(b))
    return out


def report_free():
    res = json.load(open(os.path.join(OUT, "free.json")))
    ws = [f"{a}-{b}" for a, b in WINDOWS]
    print(f"\n## Free injection per setting ({len(res)} settings; 100 Hz-10 kHz, one clock + gain per capture; median L / R)")
    print("| window (ms) | in-sample | other captures (AMPL 0, 60) | time held out (fit 0-250 ms) | repeat-capture | excess over repeat: capture / time (R) |")
    print("|---|---|---|---|---|---|")
    for w in ws:
        cells = []
        for t in ("in_sample", "capture", "time", "repeat"):
            a, b = _rows(t, w, 0), _rows(t, w, 1)
            cells.append(f"{np.median(a):.1f} / {np.median(b):.1f} (n {len(a)})" if a else "--")
        ex = [_excess(t, w, 1) for t in ("capture", "time")]
        cells.append(" / ".join(f"{np.median(e):+.1f}" if e else "--" for e in ex))
        print(f"| {w} | " + " | ".join(cells) + " |")
    for ctrl in ("rt", "hf", "lf"):
        vals = sorted({r[ctrl] for r in res})
        for t in ("capture", "time"):
            line = []
            for v_ in vals:
                sel = lambda r, v_=v_: r[ctrl] == v_
                a = [np.median(_excess(t, w, 1, sel)) if _excess(t, w, 1, sel) else np.nan for w in ("50-150", "150-300", "300-600")]
                line.append(f"{ctrl} {v_}: " + "/".join(f"{x:+.1f}" for x in a))
            print(f"- excess over repeat (R), {t}, by {ctrl} (50-150/150-300/300-600 ms): " + "; ".join(line))
    for wi, w in enumerate(OCT_WIN):
        print(f"\nExcess over repeat per octave (R, other captures, {w[0]}-{w[1]} ms; median over settings)")
        print("| settings | n | " + " | ".join(str(f) for f in OCT) + " |\n|---|---|" + "---|" * len(OCT))
        groups = [("all", lambda r: True)] + [(f"{c} {v}", lambda r, c=c, v=v: r[c] == v)
                                               for c in ("rt", "hf", "lf") for v in sorted({r[c] for r in res})]
        for name, sel in groups:
            ex = {f: [] for f in OCT}
            n = 0
            for r in res:
                if not sel(r) or "repeat" not in r:
                    continue
                n += 1
                for f in OCT:
                    a = [x["octave"][str(f)][wi] for x in r["capture"].values() if x["octave"][str(f)][wi] is not None]
                    b = [x["octave"][str(f)][wi] for x in r["repeat"].values() if x["octave"][str(f)][wi] is not None]
                    if a and b:
                        ex[f].append(np.median(a) - np.median(b))
            print(f"| {name} | {n} | " + " | ".join(f"{np.median(ex[f]):+.1f}" if ex[f] else "--" for f in OCT) + " |")


def joint_stage(splits, workers, iters, rounds, only=None):
    """q1 run_split with the corrected windows; outputs split_<name>.json and v_<name>.npy in OUT (cache -> analysis/q1)."""
    info, fade, _, ap, _ = q1_tables()
    fits = json.load(open(os.path.join(Q1, "loops.json")))
    link = os.path.join(OUT, "cache")
    if not os.path.exists(link):
        os.symlink(os.path.join(Q1, "cache"), link)
    for sp in splits:
        keys = [k for k in info if only is None or k in only]
        test = [k for k in keys if q1.SPLITS[sp](info[k])]
        train = [k for k in keys if k not in test]
        q1.run_split(OUT, info, fits, fade, f"{sp}_x0+Gx1" + ("_smoke" if only else ""), train, test, "x0+Gx1",
                     rounds=rounds, iters=iters, workers=workers, ap230=ap)


W_COMP = [(0, 50), (50, 150), (150, 300), (300, 600), (600, 1000), (1000, 1500), (1500, 2500)]


def _nulls_setting(job):
    """One setting: render the saved model per capture (unfaded), clock + gain before the fade, ramp end per capture,
    nulls cut at the fade (fade2) and compensated (model x ramp, to comp_end); compensated repeat-capture nulls."""
    key, lab, b230, law, through_g, vpath, caps, ends = job
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    G = q1.group_spectra(lab, law, q1.grid_n(lab["rt"]), through_g=through_g, b230=b230)
    v = np.load(vpath).astype(np.float64)
    out, keep = [], {}
    for c in caps:
        yb = _band(c["file"])
        n = min(len(yb), G["N"])                       # the model's grid; compensated windows end well before
        yb = yb[:n]
        J = q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=c["amp"])])])
        mdl = J.render(v, G, n, c["amp"])
        if c["fade2"] >= 0.15:
            _, (e, sh, gain) = q1._clock_job((yb, mdl, c["fade2"], c["eps"], c["amp"] == 0))
        else:                                          # (almost) the whole file is inside the fade: fit against model x ramp
            g0 = ramp_gain(n, ends[c["file"]])
            _, (e, sh, gain) = q1._clock_job((yb, mdl * g0, comp_end(yb, g0), c["eps"], c["amp"] == 0))
        y = q1.apply_clock(yb, e, sh)[:n] / gain
        end = [PRE + (x + sh - PRE) / (1 + e * 1e-6) for x in ends[c["file"]]]      # the ramp moves with the data
        g = ramp_gain(n, end)
        te = comp_end(y, g)
        keep[c["amp"]] = (y, g, te)
        out.append(dict(key=key, **lab, amp=c["amp"], file=c["file"], held_out=c.get("held_out"), eps=e, gain_db=float(20 * np.log10(gain)),
                        t_cut=c["fade2"], t_comp=te,
                        null_cut=q1.nulls(y, mdl, c["fade2"], W_COMP, narrow=True),
                        null_comp=q1.nulls(y, mdl * g, te, W_COMP, narrow=True)))
    amps = sorted(keep)
    for r in out:
        a = r["amp"]
        b = amps[(amps.index(a) + 1) % len(amps)]
        (ya, ga, ta), (yb_, gb, tb) = keep[a], keep[b]
        m = min(len(ya), len(yb_))
        r["partner"] = b
        r["repeat_cut"] = q1.nulls(ya[:m], yb_[:m], min(r["t_cut"], next(x["t_cut"] for x in out if x["amp"] == b)), W_COMP[1:], narrow=True)
        r["repeat_comp"] = q1.nulls(ya[:m] * gb[:m], yb_[:m] * ga[:m], min(ta, tb), W_COMP[1:], narrow=True)
    print(f"{key}: comp " + " ".join(f"{w}:{'--' if x is None else '%.1f/%.1f' % tuple(x)}" for w, x in out[0]["null_comp"].items() if x)
          + f"  (t_cut {out[0]['t_cut']:.2f} -> t_comp {out[0]['t_comp']:.2f} s)", flush=True)
    return out


def nulls_stage(source, workers, only=None):
    """source 'free' (per-setting models, laws from all files) or a joint split name (e.g. 'rt_x0+Gx1')."""
    from concurrent.futures import ProcessPoolExecutor
    info, fade, law_all, ap, eps0 = q1_tables()
    ends = json.load(open(os.path.join(OUT, "ramp.json")))
    jobs = []
    if source == "free":
        for r in json.load(open(os.path.join(OUT, "free.json"))):
            k = r["key"]
            caps = [dict(amp=a, file=nm, eps=r["eps"][nm], fade2=fade[nm]) for a, nm in files_of(info[k])]
            jobs.append((k, dict(lf=r["lf"], rt=r["rt"], hf=r["hf"]), r["b230"], json.loads(json.dumps(law_all)), False,
                         os.path.join(OUT, "v", k + ".npy"), caps, ends))
    else:
        R = json.load(open(os.path.join(OUT, f"split_{source}.json")))
        by = {}
        for row in R["rows"]:
            by.setdefault(row["key"], []).append(dict(amp=row["amp"], file=row["file"], eps=row["eps_ppm"], fade2=fade[row["file"]],
                                                      held_out=row["held_out"]))
        for k, caps in by.items():
            g = info[k]
            jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], R["law"], True,
                         os.path.join(OUT, f"v_{source}.npy"), caps, ends))
    if only:
        jobs = [j for j in jobs if j[0] in only]
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for rr in ex.map(_nulls_setting, jobs) for r in rr]
    json.dump(rows, open(os.path.join(OUT, f"nulls_{source}" + ("_smoke" if only else "") + ".json"), "w"), indent=0)
    report_nulls(rows, source)
    return rows


def report_nulls(rows, source):
    sel = [r for r in rows if r["held_out"] in (None, True)]
    print(f"\n## {source}: nulls cut at the fade vs fade-compensated ({len(sel)} captures"
          f"{', held out' if any(r['held_out'] for r in sel) else ''}; 100 Hz-10 kHz, clock + gain per capture, fixed fade ramp; median L / R)")
    print("| window (ms) | cut: model | cut: repeat | compensated: model | compensated: repeat | n cut / comp | excess over repeat, comp. (R) |")
    print("|---|---|---|---|---|---|---|")
    for a_, b_ in W_COMP:
        w = f"{a_}-{b_}"
        f = lambda fld: [r[fld][w] for r in sel if r[fld].get(w)]
        cells = []
        for fld in ("null_cut", "repeat_cut", "null_comp", "repeat_comp"):
            x = f(fld)
            cells.append(f"{np.median([v[0] for v in x]):.1f} / {np.median([v[1] for v in x]):.1f}" if x else "--")
        ex = [r["null_comp"][w][1] - r["repeat_comp"][w][1] for r in sel if r["null_comp"].get(w) and r["repeat_comp"].get(w)]
        print(f"| {w} | " + " | ".join(cells) + f" | {len(f('null_cut'))} / {len(f('null_comp'))} | {np.median(ex):+.1f} |" if ex else
              f"| {w} | " + " | ".join(cells) + f" | {len(f('null_cut'))} / {len(f('null_comp'))} | -- |")
    t = np.array([r["t_comp"] - r["t_cut"] for r in sel])
    print(f"- windows extended by median {np.median(t):.2f} s (IQR {np.percentile(t, 25):.2f}..{np.percentile(t, 75):.2f})")
    for rt in sorted({r["rt"] for r in sel}):
        rr = [r for r in sel if r["rt"] == rt]
        line = []
        for w in ("150-300", "300-600", "600-1000", "1000-1500"):
            x = [r["null_comp"][w][1] - r["repeat_comp"][w][1] for r in rr if r["null_comp"].get(w) and r["repeat_comp"].get(w)]
            y = [r["null_comp"][w][1] for r in rr if r["null_comp"].get(w)]
            line.append(f"{w}: {np.median(y):.1f} ({np.median(x):+.1f})" if x else f"{w}: --")
        print(f"- RT {rt}: compensated null R (excess over repeat): " + "; ".join(line))


OUT = os.path.join(ROOT, "analysis", "q2")


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="free", choices=["windows", "ramp", "free", "joint", "nulls", "report"])
    ap.add_argument("--source", default="free", help="nulls: 'free' or a joint split name, e.g. rt_x0+Gx1")
    ap.add_argument("--splits", default="rt,hf,lf")
    ap.add_argument("--rounds", type=int, default=1, help="joint: clock-refit rounds before the final fit (as q1)")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--only", default=None, help="comma-separated setting keys (lf_rt_hf)")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    OUT = a.out
    os.makedirs(OUT, exist_ok=True)
    only = a.only.split(",") if a.only else None
    if a.stage == "windows" or not os.path.exists(os.path.join(OUT, "fade2.json")):
        windows_stage()
    if a.stage == "ramp":
        ramp_stage(a.workers)
    if a.stage == "nulls":
        nulls_stage(a.source, a.workers, only)
    if a.stage == "free":
        free_stage(a.workers, a.iters, only)
    if a.stage == "joint":
        joint_stage(a.splits.split(","), a.workers, a.iters, a.rounds, only)
    if a.stage in ("free", "report"):
        report_free()


if __name__ == "__main__":
    main()
