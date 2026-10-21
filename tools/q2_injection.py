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
  residues per setting: each loop's free excitation (L, R) evaluated at the loop's poles, X_i(z_k) with time from
           cluster tap 0 (the pole residue less the common diffuser and band factors; it does not depend on where the
           excitation window starts or on the first-pass FIR), modes in 100 Hz-10 kHz; with the setting's factors at each
           mode: l = ln G_i(w_k) and its parts from the gain, the LF shelf and the HF lowpass, and the other loops' l_j.
  regress  per mode, the residues across settings against candidate bases, held out by setting (Q1's splits). Each
           residue is referenced to the end of the excitation support, rho = R exp(T l / P) (still analytic in l): with the
           origin at tap 0 a heavily damped mode is extrapolated back by up to e^3 and its estimation noise swamps the fit.
           fixed taps make the residue an analytic function of l alone (sum of exp(-D l / P)); more than that means
           the gain and filters act at different points along the line, coupling, or a non-LTI control dependence.
  branches shared injection with branches through G_i^alpha for a set of exponents (fixed taps on a line whose loss
           is applied in pieces see fractional powers; q1's x0 + G.x1 is alphas 0, 1), fitted on a subset of settings,
           held out by q1's splits (--split rt|hf|lf), fade-compensated nulls. --alphas 0,0.5,1 --subset flat|mild|all --ext.
  filters  per setting, the loop filters refitted with the injection frozen (the flat-only fit, FROZEN): 5 numbers in
           q1.G_loops' parametrisation (loop RT, LF shelf DC and Nyquist RT ratios and pole, HF lowpass ratio; the LF ones
           only where LF != 1.0), clock and gain per capture fitted once before; fade-compensated nulls against the
           repeat-capture null, with the law loops (same frozen injection) for comparison. Second form 'refit+eq': the same
           plus a loop-path EQ (first-order LF shelf and HF one-pole on the loop branches only; 3 numbers), in case the
           LF / HF controls also act outside the recursion. --only for a smoke test.
  tied     the offline 5-tap structure (model_s_results §16, analysis/q2/k1_sparsity.json) on the bench: excitations built
           from the tabled taps (L at their lags, R at lag + 26 or 332, each tap followed by -g460 at +460, which cancels the
           post comb460; g460 = 0.445 is "comb removed"), first-pass FIR 0. --mode tied (8 cluster gains + one magnitude)
           or tapgain (8 + 40 tap gains). Fitted on the flat training settings (RT 0.8 / 1.8 held out), one clock round;
           --score all scores every capture of all 128 settings with the laws. --g460 0.42,0.445 scans. --dense L adds free
           coefficients over the first L samples of every excitation (hybrid: how far the missing dense feed reaches).
  cla      cross-unit check (as q1's cla stage): the CLA IR (another unit and capture chain; HF 0.9, not a library step)
           predicted by library-trained models, fitted on CLA only the clock ratio, a common shift and one gain; nulls in
           100 Hz-10 kHz and 20 Hz-11 kHz -> cla_compare.json. Models: q1's shared injection (reference), the free flat
           injection, the compact 5-tap / -1/sqrt(2) model, and the 5-tap model with q1's -0.70 diffuser.
  (any)    Q2_HP=1 in the environment undoes the measured analogue high-pass on every capture (see hp_undo).
  diffuser the diffuser's all-pass gains (g20, g72, g230) refitted with the 5-tap tied model on the flat training settings
           (Nelder-Mead on the training cost; the 9 linear parameters solved exactly per evaluation); then scored held out
           and, with --score all, on all 128 settings with the laws. --maxfev for a smoke test.
  ap230    the AP230 damping steps re-measured with an all-pass gain --g230 (q1's ap230 stage used -0.70): per setting b
           from the first pass (240-1550 samples after tap 0) of the AMPL-0 capture against the b = 0 reference setting,
           x = ref . AP230(b) / AP230(0) -> ap230_g<g>.json (same layout as q1's ap230.json). --only for a smoke test.
           The diffuser stage takes it with --ap230 <file>.
  loops    the loop refit again, with the compact injection frozen (--frozen, default the 5-tap / -1/sqrt(2) model of
           model_s_results 16.8-16.9, its diffuser and AP230 table): per setting the law loops, 5 shared numbers ('refit',
           as `filters`), and 5 numbers per loop ('perloop', 20 per setting) -> loops_<tag>.json. --only for a smoke test.
  breport  excess of branch runs over the repeat-capture null (from the nulls stage of the same split), by control.
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

# ---------------------------------------------------------------- optional: undo the analogue high-pass on the data
# Q2_HP=1 in the environment: every capture q1.load returns is divided by the measured sub-20 Hz high-pass (the unit's
# coupling corners plus the Tomkin capture chain's 1.8 Hz, analog_results.md §3.3-3.5), above HP_FMIN, so the
# purely digital model and the data no longer differ by it (cluster and reverb alike). Set before workers start.
HP_CORNERS = [3.35, 0.72, 0.23, 0.16, 1.80]
HP_FMIN = 10.0


def hp_undo(h):
    sys.path.insert(0, os.path.join(HERE, "analog"))
    from sample_hold import highpass
    n = len(h)
    N = 1 << int(np.ceil(np.log2(n + 4096)))
    f = np.fft.rfftfreq(N, 1 / FS)
    H = highpass(f, HP_CORNERS)
    inv = np.where(f >= HP_FMIN, 1 / np.where(f > 0, H, 1), 1.0)
    return np.fft.irfft(np.fft.rfft(h, N, axis=0) * inv[:, None], N, axis=0)[:n]


if os.environ.get("Q2_HP") == "1" and not getattr(q1, "_hp_patched", False):
    _load = q1.load

    def _load_hp(out, name, band=False):
        h = hp_undo(_load(out, name))
        return sosfiltfilt(q1.ms.BAND, h, axis=0) if band else h
    q1.load = _load_hp
    q1._hp_patched = True
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


# ---------------------------------------------------------------- residues at the poles, and their setting dependence
F_BAND = (100.0, 10000.0)


def loop_factors(lab, law, w):
    """Per loop: (gain g, LF shelf H_LF(w), HF lowpass H_HF(w)) as in q1.G_law (G = g H_LF H_HF)."""
    import model_s as m
    rt, lf, hf = lab["rt"], lab["lf"], lab["hf"]
    L = law["lf"][lf]
    out = []
    for Pi in q1.P:
        kk = -3 * Pi / FS
        g = 10 ** (kk / (L["hi"] * rt))
        HL = m.shelf(w, 10 ** (kk * (1 / (L["dc"] * rt) - 1 / (L["hi"] * rt))), 1.0, law["a"]) if L["dc"] != L["hi"] else np.ones(len(w), complex)
        HH = np.ones(len(w), complex)
        if hf < 1.0:
            p = q1.onepole_for(np.array(10 ** (kk / rt * L["k"] * (1 / hf - 1))), law["f_ref"])
            HH = (1 - p) / (1 - p * np.exp(-1j * w))
        out.append((g, HL, HH))
    return out


def modes(i, lab, law):
    """Mode frequencies of loop i in F_BAND: w P_i - arg G_i(w) = 2 pi k (fixed point)."""
    Pi = int(q1.P[i])
    k = np.arange(int(np.ceil(F_BAND[0] * Pi / FS)), int(F_BAND[1] * Pi / FS) + 1)
    w = 2 * np.pi * k / Pi
    for _ in range(6):
        g, HL, HH = loop_factors(lab, law, w)[i]
        w = (2 * np.pi * k + np.angle(g * HL * HH)) / Pi
    return k, w


def _residues_setting(a):
    key, lab, law, vpath = a
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    shared = not vpath.startswith(os.path.join(OUT, "v") + os.sep)     # a joint split's x0 + G.x1 (validation)
    G = q1.group_spectra(lab, law, q1.grid_n(lab["rt"]), through_g=shared)
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])        # only to split v
    _, _, xs = J.split(np.load(vpath).astype(np.float64))
    out = {}
    for i in range(4):
        k, w = modes(i, lab, law)
        fac = loop_factors(lab, law, w)
        g, HL, HH = fac[i]
        lG = np.log(g * HL * HH)
        z = np.exp(1j * w + lG.real / q1.P[i])
        t = q1.S0 + np.arange(len(xs[0][i]))
        Z = z[:, None] ** (-t[None, :].astype(float))
        R = np.stack([Z @ xs[c][i] for c in (0, 1)], 1)                  # (modes, 2)
        if shared:                                                         # + G X1 (G = z^P at the pole)
            R = R + (g * HL * HH)[:, None] * np.stack([Z @ xs[c][i + 4] for c in (0, 1)], 1)
        other = np.stack([np.log(f[0] * f[1] * f[2]) for j, f in enumerate(fac) if j != i], 1)
        out[i] = dict(k=k, w=w, R=R, l=lG, lg=np.full(len(w), np.log(g)), lLF=np.log(HL), lHF=np.log(HH), lother=other)
    return key, out


def residues_stage(workers, only=None, source="free"):
    """source 'free' (per-setting injections) or a joint split name (its shared injection and laws: a validation,
    since a fixed-tap model must give residues that are a function of l alone)."""
    from concurrent.futures import ProcessPoolExecutor
    info, _, law, _, _ = q1_tables()
    if source != "free":
        law = json.load(open(os.path.join(OUT, f"split_{source}.json")))["law"]
    vp = lambda k: os.path.join(OUT, "v", k + ".npy") if source == "free" else os.path.join(OUT, f"v_{source}.npy")
    keys = [k for k in info if (only is None or k in only) and os.path.exists(vp(k))]
    jobs = [(k, dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]), json.loads(json.dumps(law)), vp(k)) for k in keys]
    with ProcessPoolExecutor(workers) as ex:
        res = dict(ex.map(_residues_setting, jobs))
    arr = {"keys": np.array(keys), "lab": np.array([[info[k]["lf"], info[k]["rt"], info[k]["hf"]] for k in keys])}
    for i in range(4):
        for f in ("k", "w"):
            arr[f"{f}{i}"] = res[keys[0]][i][f]
        for f in ("w", "R", "l", "lg", "lLF", "lHF", "lother"):
            arr[f"{f}_set{i}"] = np.stack([res[k][i][f] for k in keys])
    np.savez(os.path.join(OUT, "residues" + ("" if source == "free" else "_" + source) + ("_smoke" if only else "") + ".npz"), **arr)
    print(f"residues: {len(keys)} settings, modes per loop " + "/".join(str(len(arr[f"k{i}"])) for i in range(4)))


def _monomials(vs, n):
    """All products of the variables vs (list of (S, K) arrays) with total degree <= n, as a list of (S, K) arrays."""
    import itertools
    out = []
    for d in range(n + 1):
        for c in itertools.combinations_with_replacement(range(len(vs)), d):
            t = np.ones_like(vs[0])
            for j in c:
                t = t * vs[j]
            out.append(t)
    return out


def bases(F, name):
    """Design columns (list of (S, K) complex arrays) for basis `name` from the features F of one loop."""
    l, lg, lF = F["l"], F["lg"].astype(complex), F["lLF"] + F["lHF"]
    kind, n = name[0], int(name[1:])
    if kind == "A":
        return _monomials([l], n)
    if kind == "B":
        return _monomials([lg, lF], n)
    if kind == "T":
        return _monomials([lg, F["lLF"], F["lHF"]], n)
    if kind == "O":
        return _monomials([l], n) + [F["lother"][:, :, j] for j in range(F["lother"].shape[2])]
    if kind == "N":
        return _monomials([l], n) + [np.conj(l) ** d for d in range(1, n)]
    raise ValueError(name)


SPLITS_RES = {"rt": lambda lab: np.isin(lab[:, 1], (0.8, 1.8)), "hf": lambda lab: lab[:, 2] == 0.4,
              "lf": lambda lab: lab[:, 0] == 0.4}
MODELS = ["A0", "A1", "A2", "A3", "A4", "A5", "B2", "B3", "B4", "T2", "T3", "O3", "N3"]
OCT_EDGES = [100, 200, 400, 800, 1600, 3200, 6400, 10001]


def regress_stage(path=None):
    """Held-out error of each basis (dB re the residues' weighted energy), per split, overall and per octave."""
    D = np.load(path or os.path.join(OUT, "residues.npz"))
    lab = D["lab"]
    res = {}
    for sp, test_of in SPLITS_RES.items():
        te = test_of(lab)
        tr = ~te
        for name in MODELS:
            num = np.zeros(len(OCT_EDGES) - 1); den = np.zeros(len(OCT_EDGES) - 1)
            num_in = np.zeros_like(num); den_in = np.zeros_like(den)
            for i in range(4):
                F = {f: D[f"{f}_set{i}"] for f in ("l", "lg", "lLF", "lHF", "lother")}
                cols = np.stack(bases(F, name), 2)                         # (S, K, p)
                wgt = q1.P[i] / np.maximum(-2 * D[f"l_set{i}"].real, 1e-6)  # mode energy 1 / (1 - |z|^2)
                f_hz = D[f"w{i}"] * FS / (2 * np.pi)
                band = np.digitize(f_hz, OCT_EDGES) - 1
                T = q1.S0 + q1.P[i] + q1.EXT                                     # end of the excitation support
                for c in (0, 1):
                    R = D[f"R_set{i}"][:, :, c] * np.exp(T * D[f"l_set{i}"] / q1.P[i])   # mode amplitude at T (analytic in l)
                    for kk in range(R.shape[1]):
                        sw = np.sqrt(wgt[:, kk])
                        A = cols[:, kk, :] * sw[:, None]
                        y = R[:, kk] * sw
                        coef, *_ = np.linalg.lstsq(A[tr], y[tr], rcond=None)
                        e = y - A @ coef
                        b = band[kk]
                        num[b] += np.sum(np.abs(e[te]) ** 2); den[b] += np.sum(np.abs(y[te]) ** 2)
                        num_in[b] += np.sum(np.abs(e[tr]) ** 2); den_in[b] += np.sum(np.abs(y[tr]) ** 2)
            res[(sp, name)] = dict(p=cols.shape[2], held=10 * np.log10(num.sum() / den.sum()), ins=10 * np.log10(num_in.sum() / den_in.sum()),
                                   oct=(10 * np.log10(num / den)).tolist())
    print(f"\n## Residue regression across settings ({len(lab)} settings; weighted by mode energy; dB re the residues)")
    print("| basis | coeffs per mode | " + " | ".join(f"{sp} held out (in-sample)" for sp in SPLITS_RES) + " |")
    print("|---|---|" + "---|" * len(SPLITS_RES))
    for name in MODELS:
        print(f"| {name} | {res[('rt', name)]['p']} | " + " | ".join(f"{res[(sp, name)]['held']:.1f} ({res[(sp, name)]['ins']:.1f})" for sp in SPLITS_RES) + " |")
    for sp in SPLITS_RES:
        print(f"\nheld out per octave, {sp} split (edges {OCT_EDGES[:-1]} Hz)")
        for name in MODELS:
            print(f"  {name:3s} " + " ".join(f"{x:6.1f}" for x in res[(sp, name)]["oct"]))
    json.dump({f"{a}|{b}": v for (a, b), v in res.items()}, open(os.path.join(OUT, "regress.json"), "w"), indent=0)
    return res


# ---------------------------------------------------------------- shared injection with fractional branches
def branch_factor(a, g, HL, HH):
    """One branch's factor: a number alpha -> G^alpha; 'g' gain only, 'H' both filters, 'HL' / 'HH' one filter,
    'gHL' / 'gHH' gain and one filter (a tap on the far side of some of the loop's processing)."""
    if isinstance(a, float):
        return np.exp(a * np.log(g * HL * HH))
    return {"g": g * np.ones_like(HL), "H": HL * HH, "HL": HL, "HH": HH, "gHL": g * HL, "gHH": g * HH}[a]


def diffuser_post(w, b230, diff):
    """comb460 . AP20 . AP72 . AP230 with per-stage gains diff = dict(g20, g72, g230, g460) (q1 uses one gain, -0.7)."""
    z = lambda M: np.exp(-1j * w * M)
    ap = lambda M, g: (-g + z(M)) / (1 - g * z(M))
    return ap(20, diff["g20"]) * ap(72, diff["g72"]) * q1.ap230_resp(b230, w, diff["g230"]) / (1 - diff["g460"] * z(460))


def branch_spectra(lab, law, alphas, b230, diff=None):
    """q1 group spectra with one branch per (alpha, loop): PB / (1 - G z^-P) . factor(alpha); diff overrides the diffuser."""
    G = q1.group_spectra(lab, law, q1.grid_n(lab["rt"]), b230=b230)
    w = 2 * np.pi * np.fft.rfftfreq(G["N"])
    if diff is not None:
        from scipy.signal import sosfreqz
        G["PB"] = diffuser_post(w, b230, diff) * np.abs(sosfreqz(q1.LFBAND, worN=w)[1]) ** 2
        z = lambda M: np.exp(-1j * w * M)
        dt = np.exp(-1j * w * diff.get("dtap", 0.0))       # common fractional delay of the loop paths (re the cluster)
        G["R"] = [G["PB"] * dt / (1 - Gi * z(Pi)) for Gi, Pi in zip(q1.G_law(lab, law, w), q1.P)]
    fac = loop_factors(lab, law, w)
    G["R"] = [r * branch_factor(a, *f) for a in alphas for r, f in zip(G["R"][:4], fac)]
    return G


def score_capture(G, v, amp, nm, fade2, eps, end_dev, ext=None):
    """Clock + gain before the fade, then cut and fade-compensated nulls (as _nulls_setting) for one capture."""
    if ext is not None:
        q1.EXT = ext                                       # excitation length P_i + EXT (worker processes start at q1's 300)
    yb = _band(nm)
    n = min(len(yb), G["N"])
    yb = yb[:n]
    mdl = q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=amp)])]).render(v, G, n, amp)
    if fade2 >= 0.15:
        _, (e, sh, gain) = q1._clock_job((yb, mdl, fade2, eps, amp == 0))
    else:
        g0 = ramp_gain(n, end_dev)
        _, (e, sh, gain) = q1._clock_job((yb, mdl * g0, comp_end(yb, g0), eps, amp == 0))
    y = q1.apply_clock(yb, e, sh)[:n] / gain
    g = ramp_gain(n, [PRE + (x + sh - PRE) / (1 + e * 1e-6) for x in end_dev])
    te = comp_end(y, g)
    return dict(null_cut=q1.nulls(y, mdl, fade2, W_COMP, narrow=True), null_comp=q1.nulls(y, mdl * g, te, W_COMP, narrow=True), t_comp=te)


def _score_job(a):
    return a[0], a[1], score_capture(*a[2:])


def branches_stage(alphas, subset, workers, iters, rounds=1, ext=None, split="rt"):
    from concurrent.futures import ProcessPoolExecutor
    ext = q1.EXT if ext is None else ext
    q1.EXT = ext
    info, fade, _, ap, eps0 = q1_tables()
    ends = json.load(open(os.path.join(OUT, "ramp.json")))
    fits = json.load(open(os.path.join(Q1, "loops.json")))
    keys = [k for k, g in info.items() if subset == "all" or (subset == "flat" and g["lf"] == 1.0 and g["hf"] == 1.0)
            or (subset == "mild" and g["lf"] != 0.1 and g["hf"] != 0.2)]
    test = [k for k in keys if q1.SPLITS[split](info[k])]
    train = [k for k in keys if k not in test]
    law = q1.fit_laws(fits, set(train)) if subset != "flat" else q1.fit_laws(fits)
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    b230 = {k: ap["steps"][ap["settings"][k]["step"]] for k in keys}
    lab = {k: dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]) for k in keys}
    eps = dict(eps0)
    tag = f"branches_{subset}_" + "_".join(f"{a:g}" if isinstance(a, float) else a for a in alphas) + ("" if ext == 300 else f"_ext{ext}") + ("" if split == "rt" else f"_{split}")
    if iters < 100:
        tag += "_smoke"

    def build():
        fl = [(k, amp, nm) for k in train for amp, nm in files_of(info[k]) if amp in q1.TRAIN_AMPL]
        with ProcessPoolExecutor(workers) as ex:
            ys = list(ex.map(q1._prep_job, [(Q1, nm, info[k]["rt"], eps[nm], fade[nm], True, 0.0) for k, _, nm in fl], chunksize=4))
        gs = []
        for k in train:
            G = branch_spectra(lab[k], law, alphas, b230[k])
            G["files"] = [dict(y=y, w=w, amp=amp, name=nm) for (k2, amp, nm), (y, w) in zip(fl, ys) if k2 == k]
            gs.append(G)
        return q1.Joint(gs)
    v = None
    for rd in range(rounds + 1):
        J = build()
        v = J.solve(iters, x0=v)
        if rd < rounds:                                    # clock refit on the training captures
            jobs = []
            for G in J.groups:
                for fl in G["files"]:
                    yb = _band(fl["name"])
                    jobs.append((yb, J.render(v, G, min(len(yb), G["N"]), fl["amp"]), fade[fl["name"]], eps[fl["name"]], False))
            with ProcessPoolExecutor(workers) as ex:
                cl = list(ex.map(q1._clock_job, jobs))
            eps.update({fl["name"]: c[0][0] for G in J.groups for fl in G["files"] for c in [cl.pop(0)]})
    jobs = []
    rng = np.random.default_rng(0)                          # scored training settings: a seeded random 8, spread over the controls
    for k in test + sorted(rng.choice(train, min(8, len(train)), replace=False).tolist()):
        G = branch_spectra(lab[k], law, alphas, b230[k])
        for amp, nm in files_of(info[k]):
            jobs.append((k, amp, G, v, amp, nm, fade[nm], eps0[nm], ends[nm], ext))
    with ProcessPoolExecutor(workers) as ex:
        sc = list(ex.map(_score_job, jobs))
    rows = [dict(key=k, amp=amp, held_out=k in test, **lab[k], **r) for k, amp, r in sc]
    json.dump(dict(alphas=alphas, subset=subset, ext=ext, split=split, law=law, train=train, test=test, params=J.size(), rows=rows),
              open(os.path.join(OUT, tag + ".json"), "w"), indent=0)
    np.save(os.path.join(OUT, tag + ".npy"), v)
    for held in (True, False):
        rr = [r for r in rows if r["held_out"] == held]
        print(f"[{tag}] {'held out' if held else 'training'} ({len(rr)} captures, {J.size()} params), compensated R: " + "  ".join(
            f"{a}-{b}: {np.median([r['null_comp'][f'{a}-{b}'][1] for r in rr if r['null_comp'].get(f'{a}-{b}')]):.1f}"
            for a, b in W_COMP if any(r["null_comp"].get(f"{a}-{b}") for r in rr)), flush=True)


def branch_report(tags, sel=None):
    """Per branch run: held-out compensated null (R) and its excess over the repeat-capture null of the same capture
    (taken from nulls_<split>_x0+Gx1.json), overall and by control; sel optionally restricts the held-out captures."""
    ws = ["50-150", "150-300", "300-600", "600-1000"]
    for tag in tags:
        d = json.load(open(os.path.join(OUT, tag + ".json")))
        rep = {(r["key"], r["amp"]): r for r in json.load(open(os.path.join(OUT, f"nulls_{d.get('split', 'rt')}_x0+Gx1.json")))}
        rows = [r for r in d["rows"] if r["held_out"] and (sel is None or sel(r))]

        def ex(rr, w):
            e = [r["null_comp"][w][1] - rep[(r["key"], r["amp"])]["repeat_comp"][w][1] for r in rr
                 if r["null_comp"].get(w) and (r["key"], r["amp"]) in rep and rep[(r["key"], r["amp"])]["repeat_comp"].get(w)]
            return np.median(e) if e else np.nan
        print(f"{tag} ({len(rows)} held-out captures, {d['params']} params): " + "; ".join(
            f"{w}: {np.median([r['null_comp'][w][1] for r in rows if r['null_comp'].get(w)]):.1f} [{ex(rows, w):+.1f}]" for w in ws))
        for c in ("rt", "hf", "lf"):
            vals = sorted({r[c] for r in rows})
            if len(vals) > 1:
                print(f"   by {c}: " + "; ".join(f"{v}: " + "/".join(f"{ex([r for r in rows if r[c] == v], w):+.1f}" for w in ws[:3]) for v in vals))


# ---------------------------------------------------------------- loop filters refitted per setting, injection frozen
FROZEN = "branches_flat_0_ext1700"         # flat-only shared injection, 600 iterations (at the repeat-capture null there)


def _filters_setting(a):
    """Least-squares refit of one setting's loop parameters (q1.G_loops) against all its captures, injection frozen."""
    from scipy.optimize import least_squares
    key, lab, b230, law, v, ext, caps, ends, v0 = a
    q1.EXT = ext
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    N = q1.grid_n(lab["rt"])
    w = 2 * np.pi * np.fft.rfftfreq(N)
    G0 = q1.group_spectra(lab, law, N, b230=b230)
    J = q1.Joint([dict(G0, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
    cl, f, xs = J.split(v)
    F, X = J._spectra(f, xs, N)                                     # excitation spectra: fixed
    C = J.cluster(cl)
    z = lambda M: np.exp(-1j * w * M)

    import model_s as m

    def eq_of(q):
        """Loop-path EQ (form 'loops+eq'): first-order LF shelf (DC gain q0, pole q1) and HF one-pole (pole q2)."""
        return m.shelf(w, q[0], 1.0, q[1]) * (1 - q[2]) / (1 - q[2] * np.exp(-1j * w))

    def reverb(Gs, E=1.0):
        R = np.array([G0["PB"] * E / (1 - Gi * z(Pi)) for Gi, Pi in zip(Gs, q1.P)])
        return np.fft.irfft(F * G0["PB"][None] + np.einsum("cik,ik->ck", X, R), N, axis=1).T

    # data: clock and gain per capture against the law-loop model, before the fade; then fixed
    data = []
    r_law = reverb(q1.G_law(lab, law, w))
    for c in caps:
        yb = _band(c["file"])
        n = min(len(yb), N, len(C))
        mdl = (c["amp"] / 90) * C[:n] + r_law[:n]
        _, (e, sh, gain) = q1._clock_job((yb[:n], mdl, max(c["fade2"], 0.15), c["eps"], c["amp"] == 0))
        y = q1.apply_clock(yb[:n], e, sh) / gain
        i0, i1 = PRE + int(0.05 * FS), min(PRE + int(max(c["fade2"], 0.15) * FS), n)
        t = (np.arange(i0, i1) - PRE) / FS
        wt = 10 ** (1.5 * np.minimum(t, q1.T_CAP) / lab["rt"])[:, None]
        data.append(dict(c, y=y, n=n, i0=i0, i1=i1, wt=wt / np.sqrt(np.sum((y[i0:i1] * wt) ** 2)), e=e, sh=sh))
    lfree = lab["lf"] != 1.0
    idx = [0, 4] + ([1, 2, 3] if lfree else [])
    lo = np.array([0.05, 0.02, 0.2, -0.5, 0.05, -300]); hi = np.array([20, 20, 5, 0.99999, 1.0, 300])
    vv = np.array(v0, float); vv[5] = 0.0
    vv[4] = min(vv[4], 1.0)

    def res(u, eq=False):
        p = vv.copy(); p[idx] = u[:len(idx)]
        r = reverb(q1.G_loops(p, w), eq_of(u[len(idx):]) if eq else 1.0)
        return np.concatenate([(((d["amp"] / 90) * C[:d["n"]] + r[:d["n"]])[d["i0"]:d["i1"]] - d["y"][d["i0"]:d["i1"]]) * d["wt"] for d in data]).ravel()
    c0 = float(np.sum(res(vv[idx]) ** 2))
    xs_ = np.array([0.001 * lab["rt"], 0.01, 0.01, 0.001, 0.01, 1.0])[idx]
    fit = least_squares(res, vv[idx], bounds=(lo[idx], hi[idx]), x_scale=xs_, diff_step=1e-6, max_nfev=80)
    vv[idx] = fit.x
    u0 = np.r_[fit.x, 1.0, 0.9, 0.0]                              # + loop-path EQ, starting flat
    fe = least_squares(lambda u: res(u, True), u0, bounds=(np.r_[lo[idx], 0.05, 0.0, 0.0], np.r_[hi[idx], 20.0, 0.995, 0.95]),
                       x_scale=np.r_[xs_, 0.01, 0.001, 0.01], diff_step=1e-6, max_nfev=120)
    ve = vv.copy(); ve[idx] = fe.x[:len(idx)]
    E = eq_of(fe.x[len(idx):])
    out = []
    for name, Gs, Ex in (("law", q1.G_law(lab, law, w), None), ("refit", q1.G_loops(vv, w), None), ("refit+eq", q1.G_loops(ve, w), E)):
        G = q1.group_spectra(lab, law, N, Gs=Gs, b230=b230)
        if Ex is not None:
            G["R"] = [r_ * Ex for r_ in G["R"]]
        for c in caps:
            out.append(dict(key=key, **lab, amp=c["amp"], file=c["file"], loops=name,
                            **score_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], ext)))
    print(f"{key}: cost {c0:.3g} -> {fit.cost * 2:.3g} (loops) -> {fe.cost * 2:.3g} (loops+eq); v {np.round(vv, 4).tolist()}; "
          f"eq {np.round(fe.x[len(idx):], 4).tolist()}", flush=True)
    return dict(key=key, **lab, v=vv.tolist(), v_eq=ve.tolist(), eq=fe.x[len(idx):].tolist(), free=idx, cost0=c0,
                cost=float(fit.cost * 2), cost_eq=float(fe.cost * 2), nfev=fit.nfev, rows=out)


def filters_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    info, fade, law, ap, eps0 = q1_tables()
    ends = json.load(open(os.path.join(OUT, "ramp.json")))
    fits = json.load(open(os.path.join(Q1, "loops.json")))
    d = json.load(open(os.path.join(OUT, FROZEN + ".json")))
    v = np.load(os.path.join(OUT, FROZEN + ".npy"))
    flaw = d["law"]
    keys = [k for k in info if only is None or k in only]
    jobs = []
    for k in keys:
        g = info[k]
        caps = [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in files_of(g)]
        v0 = np.median([r["v"] for r in fits if r["file"] in {c["file"] for c in caps}], 0)
        jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(flaw)),
                     v, d["ext"], caps, ends, v0))
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_filters_setting, jobs))
    json.dump(res, open(os.path.join(OUT, "filters" + ("_smoke" if only else "") + ".json"), "w"), indent=0)
    filters_report(res)


def filters_report(res):
    rep = {(r["key"], r["amp"]): r for r in json.load(open(os.path.join(OUT, "nulls_free.json")))}
    ws = ["50-150", "150-300", "300-600", "600-1000"]

    def ex(rows, w):
        e = [r["null_comp"][w][1] - rep[(r["key"], r["amp"])]["repeat_comp"][w][1] for r in rows
             if r["null_comp"].get(w) and rep.get((r["key"], r["amp"])) and rep[(r["key"], r["amp"])]["repeat_comp"].get(w)]
        return np.median(e) if e else np.nan
    rows = [r for s in res for r in s["rows"]]
    print(f"\n## Loop filters refitted per setting, injection frozen ({FROZEN}); {len(res)} settings; R, compensated; null [excess over repeat]")
    for lp in ("law", "refit", "refit+eq"):
        rr = [r for r in rows if r["loops"] == lp]
        print(f"- {lp:5s}: " + "; ".join(f"{w}: {np.median([r['null_comp'][w][1] for r in rr if r['null_comp'].get(w)]):.1f} [{ex(rr, w):+.1f}]" for w in ws))
        for c in ("rt", "hf", "lf"):
            vals = sorted({r[c] for r in rr})
            if len(vals) > 1:
                print(f"    by {c}: " + "; ".join(f"{v_}: " + "/".join(f"{ex([r for r in rr if r[c] == v_], w):+.1f}" for w in ws[:3]) for v_ in vals))
    if len(res) <= 12:                                               # small runs: per setting (group medians mix settings)
        for s_ in res:
            print(f"  {s_['key']:12s} " + "  ".join(f"{lp}: " + "/".join(f"{ex([r for r in s_['rows'] if r['loops'] == lp], w):+.1f}" for w in ws[:3])
                                                for lp in ("law", "refit", "refit+eq")))


# ---------------------------------------------------------------- the 5-tap structure on the bench
def tap_table(path=None):
    """(loop index, channel, excitation index, sign) per tap from k1_sparsity.json (lags from cluster tap 0)."""
    t = json.load(open(path or os.path.join(OUT, "k1_sparsity.json")))["tied"]
    taps = []
    for i, Pi in enumerate(q1.P):
        e = t[str(Pi)]
        for j, ((lag, sg), sh) in enumerate(zip(e["L_taps"], e["R_shift"])):
            taps.append((i, 0, lag - q1.S0, sg, j))
            taps.append((i, 1, lag + sh - q1.S0, sg, j))
    return taps


def tied_map(J, taps, g460, mode, dense=0):
    """Sparse matrix M: full Joint parameter vector = M theta. theta = 8 cluster gains + (one magnitude | one gain per tap)
    [+ with dense > 0: free coefficients for the first `dense` samples of every loop/channel excitation (hybrid)]."""
    from scipy.sparse import lil_matrix
    off = {}
    o = 8 + 2 * q1.LFIR
    for c in (0, 1):
        for i, L in enumerate(J.Lx):
            off[(c, i)] = o
            o += L
    K = int(mode[6:]) if mode.startswith("kernel") else 0       # "kernelK": one kernel at lags -K..K shared by every tap
    ntap = 8 + (1 if mode == "tied" else 2 * K + 1 if mode.startswith("kernel") else len(taps))
    ncol = ntap + dense * 8
    M = lil_matrix((J.size(), ncol))
    if dense:
        assert dense <= min(J.Lx), "dense part longer than the excitation window"
        col = ntap
        for c in (0, 1):
            for i in range(4):
                for t in range(dense):
                    M[off[(c, i)] + t, col] = 1.0
                    col += 1
    for k in range(8):
        M[k, k] = 1.0
    for n, (i, c, t, sg, j) in enumerate(taps):
        assert 0 <= t - K and t + K + 460 < J.Lx[i], "tap outside the excitation window"
        if mode.startswith("kernel"):
            for k in range(-K, K + 1):
                M[off[(c, i)] + t + k, 8 + K + k] += sg
                M[off[(c, i)] + t + k + 460, 8 + K + k] -= g460 * sg
            continue
        col = 8 if mode == "tied" else 8 + n
        M[off[(c, i)] + t, col] += sg
        M[off[(c, i)] + t + 460, col] -= g460 * sg
    return M.tocsr()


def tied_stage(mode, g460s, workers, iters, rounds, ext, score, dense=0, gfix=None, ap230=None):
    from concurrent.futures import ProcessPoolExecutor
    from scipy.sparse.linalg import LinearOperator, lsqr
    q1.EXT = ext
    info, fade, _, ap, eps0 = q1_tables()
    ends = json.load(open(os.path.join(OUT, "ramp.json")))
    fits = json.load(open(os.path.join(Q1, "loops.json")))
    law = q1.fit_laws(fits)
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    flat = [k for k, g in info.items() if g["lf"] == 1.0 and g["hf"] == 1.0]
    test = [k for k in flat if info[k]["rt"] in (0.8, 1.8)]
    train = [k for k in flat if k not in test]
    lab = {k: dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]) for k in info}
    if ap230:
        ap = json.load(open(ap230))
    b230 = {k: ap["steps"][ap["settings"][k]["step"]] for k in info}
    taps = tap_table()
    eps = dict(eps0)
    dif = None if gfix is None else dict(g20=gfix, g72=gfix, g230=gfix, g460=0.445)

    def build():
        fl = [(k, amp, nm) for k in train for amp, nm in files_of(info[k]) if amp in q1.TRAIN_AMPL]
        with ProcessPoolExecutor(workers) as ex:
            ys = list(ex.map(q1._prep_job, [(Q1, nm, info[k]["rt"], eps[nm], fade[nm], True, 0.0) for k, _, nm in fl], chunksize=4))
        gs = []
        for k in train:
            G = branch_spectra(lab[k], law, [0.0], b230[k], dif)
            G["files"] = [dict(y=y, w=w, amp=amp, name=nm) for (k2, amp, nm), (y, w) in zip(fl, ys) if k2 == k]
            gs.append(G)
        return q1.Joint(gs)
    for g460 in g460s:
        tag = (f"tied_{mode}_g{g460:g}" + (f"_dense{dense}_ext{ext}" if dense else "") + ("" if gfix is None else f"_ap{abs(gfix):.4f}")
               + ("_smoke" if iters < 30 else ""))
        th = None
        for rd in range(rounds + 1):
            J = build()
            M = tied_map(J, taps, g460, mode, dense)
            op = LinearOperator((len(J.target()), M.shape[1]), matvec=lambda t: J.forward(M @ t),
                                rmatvec=lambda z: M.T @ J.adjoint(z), dtype=float)
            th = lsqr(op, J.target(), iter_lim=iters, atol=1e-12, btol=1e-12, x0=th)[0]
            v = M @ th
            if rd < rounds:
                jobs, names = [], []
                for G in J.groups:
                    for fl_ in G["files"]:
                        yb = _band(fl_["name"])
                        jobs.append((yb, J.render(v, G, min(len(yb), G["N"]), fl_["amp"]), fade[fl_["name"]], eps[fl_["name"]], False))
                        names.append(fl_["name"])
                with ProcessPoolExecutor(workers) as ex:
                    cl = list(ex.map(q1._clock_job, jobs))
                eps.update({nm: c_[0][0] for nm, c_ in zip(names, cl)})
        keys = test + sorted(np.random.default_rng(0).choice(train, min(8, len(train)), replace=False).tolist())
        if score == "all":
            keys = list(info)
        jobs = [(k, amp, branch_spectra(lab[k], law, [0.0], b230[k], dif), v, amp, nm, fade[nm], eps0[nm], ends[nm], ext)
                for k in keys for amp, nm in files_of(info[k])]
        with ProcessPoolExecutor(workers) as ex:
            sc = list(ex.map(_score_job, jobs))
        rows = [dict(key=k, amp=amp, held_out=k not in train, **lab[k], **r) for k, amp, r in sc]
        nfree = M.shape[1]
        json.dump(dict(mode=mode, g460=g460, ext=ext, split="rt", law=law, train=train, test=test, params=nfree,
                       theta=th.tolist(), rows=rows), open(os.path.join(OUT, tag + ".json"), "w"), indent=0)
        np.save(os.path.join(OUT, tag + ".npy"), v)
        print(f"[{tag}] {nfree} free parameters; tap magnitude / kernel {np.round(th[8:8 + (1 if mode == 'tied' else 17 if mode.startswith('kernel') else len(taps))], 4).tolist()[:17]}",
              flush=True)
        branch_report([tag], sel=lambda r: r["key"] in test)
        if score == "all":
            print("  all 128 settings (flat training settings included):")
            branch_report([tag])


def ap230_stage(g230, only=None):
    """q1.ap230_stage with the all-pass gain g230 (the ratio AP230(b) / AP230(0) depends on it)."""
    from scipy.optimize import minimize_scalar
    info = json.load(open(os.path.join(Q1, "data.json")))
    N = 1 << 14
    w = 2 * np.pi * np.fft.rfftfreq(N)
    ref = sosfiltfilt(q1.LFBAND, q1.load(Q1, info[q1.AP230_REF]["files"]["0"]["file"]), axis=0)
    R = np.fft.rfft(ref[:N], axis=0)
    i0, i1 = PRE + 240, PRE + 1550
    res = {}
    for key, g in info.items():
        if only and key not in only:
            continue
        x = sosfiltfilt(q1.LFBAND, q1.load(Q1, g["files"]["0"]["file"]), axis=0)

        def cost(b):
            mdl = np.fft.irfft(R * (q1.ap230_resp(b, w, g230) / q1.ap230_resp(0.0, w, g230))[:, None], N, axis=0)
            return np.sum((x[i0:i1] - mdl[i0:i1]) ** 2) / np.sum(x[i0:i1] ** 2)
        r = minimize_scalar(cost, bounds=(-0.1, 0.9), method="bounded", options={"xatol": 1e-5})
        res[key] = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"], b=float(r.x), null_db=float(10 * np.log10(r.fun)), null_b0_db=float(10 * np.log10(cost(0.0))))
    bs = np.sort([v["b"] for v in res.values()])
    levels = []
    for b in bs:                                            # the steps: cluster the fitted values (gaps > 0.05), as q1
        if not levels or b - levels[-1][-1] > 0.05:
            levels.append([b])
        else:
            levels[-1].append(b)
    steps = [float(np.median(l)) for l in levels]
    for v in res.values():
        v["step"] = int(np.argmin([abs(v["b"] - st) for st in steps]))
    old = json.load(open(os.path.join(Q1, "ap230.json")))
    path = os.path.join(OUT, f"ap230_g{abs(g230):.4f}" + ("_smoke" if only else "") + ".json")
    json.dump(dict(g230=g230, steps=steps, settings=res), open(path, "w"), indent=0)
    print(f"AP230 steps with g230 = {g230:.5f}: {np.round(steps, 4).tolist()}  (q1, g -0.70: {np.round(old['steps'], 4).tolist()})")
    nb = [v["null_db"] for v in res.values() if v["step"] > 0]
    ch = sum(1 for k, v in res.items() if v["step"] != old["settings"][k]["step"])
    print(f"  {len(res)} settings; first-pass null median {np.median([v['null_db'] for v in res.values()]):.1f} dB "
          f"(step > 0: {np.median(nb) if nb else float('nan'):.1f} dB); settings whose step index changed: {ch}")
    for k, v in sorted(res.items(), key=lambda kv: (kv[1]["rt"], kv[1]["hf"], kv[1]["lf"])):
        if v["step"] > 0 or old["settings"][k]["step"] > 0:
            o = old["settings"][k]
            print(f"  {k}: b {v['b']:.4f} (q1 {o['b']:.4f}), first-pass null {v['null_db']:.1f} dB (q1 {o['null_db']:.1f})")
    return path


def _prep_job_q2(a):
    """q1._prep_job run from this module, so that a Q2_HP patch of q1.load is active in the worker."""
    return q1.prep(*a[:5], weight=a[5], shift=a[6])


def diffuser_stage(workers, score, maxfev, ext=1700, g460=0.445, gfix=None, ap230=None, dtap=False):
    from concurrent.futures import ProcessPoolExecutor
    from scipy.optimize import minimize
    q1.EXT = ext
    info, fade, _, ap, eps0 = q1_tables()
    if ap230:
        ap = json.load(open(ap230))                        # re-measured AP230 steps (ap230 stage)
    ends = json.load(open(os.path.join(OUT, "ramp.json")))
    law = q1.fit_laws(json.load(open(os.path.join(Q1, "loops.json"))))
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    flat = [k for k, g in info.items() if g["lf"] == 1.0 and g["hf"] == 1.0]
    test = [k for k in flat if info[k]["rt"] in (0.8, 1.8)]
    train = [k for k in flat if k not in test]
    lab = {k: dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]) for k in info}
    b230 = {k: ap["steps"][ap["settings"][k]["step"]] for k in info}
    taps = tap_table()
    fl = [(k, amp, nm) for k in train for amp, nm in files_of(info[k]) if amp in q1.TRAIN_AMPL]
    with ProcessPoolExecutor(workers) as ex:
        ys = list(ex.map(_prep_job_q2, [(Q1, nm, info[k]["rt"], eps0[nm], fade[nm], True, 0.0) for k, _, nm in fl], chunksize=4))
    hist = []

    def solve(gs, d=0.0):
        diff = dict(g20=gs[0], g72=gs[1], g230=gs[2], g460=g460, dtap=float(d))
        groups = []
        for k in train:
            G = branch_spectra(lab[k], law, [0.0], b230[k], diff)
            G["files"] = [dict(y=y, w=w, amp=amp, name=nm) for (k2, amp, nm), (y, w) in zip(fl, ys) if k2 == k]
            groups.append(G)
        J = q1.Joint(groups)
        M = tied_map(J, taps, g460, "tied")
        A = np.stack([J.forward(M[:, j].toarray().ravel()) for j in range(M.shape[1])], 1)
        y = J.target()
        th, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ th
        cost = float(r @ r / (y @ y))
        hist.append(dict(g=list(map(float, gs)), dtap=float(d), cost_db=10 * np.log10(cost), a=float(th[8])))
        print(f"  g20 {gs[0]:+.4f} g72 {gs[1]:+.4f} g230 {gs[2]:+.4f} dtap {d:+.4f}: training cost {10 * np.log10(cost):.2f} dB, a {th[8]:.4f}", flush=True)
        return cost, M @ th, diff
    g0 = np.array([-0.70, -0.70, -0.70])
    c0 = solve(g0)[0]
    if gfix is None:
        res = minimize(lambda g: solve(g)[0], g0, method="Nelder-Mead",
                       options=dict(initial_simplex=[g0, g0 + [0.004, 0, 0], g0 + [0, 0.004, 0], g0 + [0, 0, 0.004]],
                                    xatol=5e-4, fatol=1e-6, maxfev=maxfev))
        gbest = res.x
    else:
        gbest = np.full(3, gfix)                             # all three fixed (no gain fitted)
    dbest = 0.0
    if dtap:                                                 # common fractional delay of the loop paths, gains as above
        from scipy.optimize import minimize_scalar
        coarse = {d: solve(gbest, d)[0] for d in (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0)}
        d0 = min(coarse, key=coarse.get)
        dbest = float(minimize_scalar(lambda d: solve(gbest, d)[0], bounds=(d0 - 0.25, d0 + 0.25), method="bounded",
                                      options={"xatol": 0.005, "maxiter": maxfev}).x)
    cost, v, diff = solve(gbest, dbest)
    tag = ("diffuser_tied" + ("" if gfix is None else f"_fix{gfix:g}") + ("_ap230new" if ap230 else "") + ("_dtap" if dtap else "")
           + ("_hp" if os.environ.get("Q2_HP") == "1" else "")
           + ("_smoke" if maxfev < 20 and (gfix is None or dtap) else ""))
    keys = test + sorted(np.random.default_rng(0).choice(train, min(8, len(train)), replace=False).tolist())
    if score == "all":
        keys = list(info)
    rows = []
    base = dict(g20=-0.7, g72=-0.7, g230=-0.7, g460=g460) if not dtap else dict(g20=gbest[0], g72=gbest[1], g230=gbest[2], g460=g460, dtap=0.0)
    for name, dd in (("refit", diff), ("q1 diffuser" if not dtap else "no dtap", base)):
        vv = v if name == "refit" else solve(g0 if not dtap else gbest)[1]
        jobs = [(k, amp, branch_spectra(lab[k], law, [0.0], b230[k], dd), vv, amp, nm, fade[nm], eps0[nm], ends[nm], ext)
                for k in keys for amp, nm in files_of(info[k])]
        with ProcessPoolExecutor(workers) as ex:
            sc = list(ex.map(_score_job, jobs))
        sub = [dict(key=k, amp=amp, held_out=k not in train, diffuser=name, **lab[k], **r) for k, amp, r in sc]
        rows += sub
        json.dump(dict(diff=dd, train=train, test=test, split="rt", params=9 + (3 if name == "refit" and gfix is None else 0) + (1 if name == "refit" and dtap else 0), rows=sub, hist=hist),
                  open(os.path.join(OUT, tag + ("" if name == "refit" else "_q1diff") + ".json"), "w"), indent=0)
    np.save(os.path.join(OUT, tag + ".npy"), v)
    print(f"[{tag}] g20 {diff['g20']:+.4f}  g72 {diff['g72']:+.4f}  g230 {diff['g230']:+.4f}  dtap {diff.get('dtap', 0.0):+.4f}; training cost {10 * np.log10(c0):.2f} -> "
          f"{10 * np.log10(cost):.2f} dB ({len(hist)} evaluations)")
    for t_ in (tag, tag + "_q1diff"):
        branch_report([t_], sel=lambda r: r["key"] in test)
        if score == "all":
            print("  all 128 settings:")
            branch_report([t_])


COMPACT = "diffuser_tied_fix-0.707107_ap230new"     # 5 taps per loop, all-pass gains -1/sqrt(2) (model_s_results 16.8-16.9)


def make_G(lab, law, Gs, b230, diff, w=None):
    """Group spectra with the diffuser override and given loop responses Gs."""
    G = branch_spectra(lab, law, [0.0], b230, diff)
    w = 2 * np.pi * np.fft.rfftfreq(G["N"])
    z = lambda M: np.exp(-1j * w * M)
    dt = np.exp(-1j * w * diff.get("dtap", 0.0))
    G["R"] = [G["PB"] * dt / (1 - Gi * z(Pi)) for Gi, Pi in zip(Gs, q1.P)]
    return G


def _loops_setting(a):
    """Per setting: law loops, 5 shared loop numbers and 5 per loop, fitted with the compact injection frozen."""
    from scipy.optimize import least_squares
    key, lab, b230, law, v, ext, caps, ends, v0, diff = a
    q1.EXT = ext
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    G0 = branch_spectra(lab, law, [0.0], b230, diff)
    N = G0["N"]
    w = 2 * np.pi * np.fft.rfftfreq(N)
    z = lambda M: np.exp(-1j * w * M)
    J = q1.Joint([dict(G0, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
    cl, f, xs = J.split(v)
    F, X = J._spectra(f, xs, N)
    C = J.cluster(cl)

    def reverb(Gs):
        R = np.array([G0["PB"] / (1 - Gi * z(Pi)) for Gi, Pi in zip(Gs, q1.P)])
        return np.fft.irfft(F * G0["PB"][None] + np.einsum("cik,ik->ck", X, R), N, axis=1).T
    data = []
    r_law = reverb(q1.G_law(lab, law, w))
    for c in caps:
        yb = _band(c["file"])
        n = min(len(yb), N, len(C))
        mdl = (c["amp"] / 90) * C[:n] + r_law[:n]
        _, (e, sh, gain) = q1._clock_job((yb[:n], mdl, max(c["fade2"], 0.15), c["eps"], c["amp"] == 0))
        y = q1.apply_clock(yb[:n], e, sh) / gain
        i0, i1 = PRE + int(0.05 * FS), min(PRE + int(max(c["fade2"], 0.15) * FS), n)
        t = (np.arange(i0, i1) - PRE) / FS
        wt = 10 ** (1.5 * np.minimum(t, q1.T_CAP) / lab["rt"])[:, None]
        data.append(dict(c, y=y, n=n, i0=i0, i1=i1, wt=wt / np.sqrt(np.sum((y[i0:i1] * wt) ** 2))))
    idx = [0, 4] + ([1, 2, 3] if lab["lf"] != 1.0 else [])
    lo = np.array([0.05, 0.02, 0.2, -0.5, 0.05, -300]); hi = np.array([20, 20, 5, 0.99999, 1.0, 300])
    xsc = np.array([0.001 * lab["rt"], 0.01, 0.01, 0.001, 0.01, 1.0])[idx]
    vv = np.array(v0, float); vv[5] = 0.0; vv[4] = min(vv[4], 1.0)

    def err(Gs):
        r = reverb(Gs)
        return np.concatenate([(((d["amp"] / 90) * C[:d["n"]] + r[:d["n"]])[d["i0"]:d["i1"]] - d["y"][d["i0"]:d["i1"]]) * d["wt"] for d in data]).ravel()

    def shared(u):
        p = vv.copy(); p[idx] = u
        return q1.G_loops(p, w)

    def perloop(u):
        out = []
        for i in range(4):
            p = vv.copy(); p[idx] = u[i * len(idx):(i + 1) * len(idx)]
            out.append(q1.G_loops(p, w)[i])
        return out
    c_law = float(np.sum(err(q1.G_law(lab, law, w)) ** 2))
    f1 = least_squares(lambda u: err(shared(u)), vv[idx], bounds=(lo[idx], hi[idx]), x_scale=xsc, diff_step=1e-6, max_nfev=80)
    u2 = np.tile(f1.x, 4)
    f2 = least_squares(lambda u: err(perloop(u)), u2, bounds=(np.tile(lo[idx], 4), np.tile(hi[idx], 4)), x_scale=np.tile(xsc, 4),
                       diff_step=1e-6, max_nfev=300)
    out = []
    for name, Gs in (("law", q1.G_law(lab, law, w)), ("refit", shared(f1.x)), ("perloop", perloop(f2.x))):
        G = make_G(lab, law, Gs, b230, diff)
        for c in caps:
            out.append(dict(key=key, **lab, amp=c["amp"], file=c["file"], loops=name,
                            **score_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], ext)))
    per = []
    for i in range(4):
        p = vv.copy(); p[idx] = f2.x[i * len(idx):(i + 1) * len(idx)]
        per.append(p[:5].tolist())
    p_sh = vv.copy(); p_sh[idx] = f1.x
    print(f"{key}: cost law {c_law:.3g} -> shared {2 * f1.cost:.3g} -> per loop {2 * f2.cost:.3g}", flush=True)
    return dict(key=key, **lab, free=idx, v_shared=p_sh[:5].tolist(), v_perloop=per, cost_law=c_law, cost_shared=float(2 * f1.cost),
                cost_perloop=float(2 * f2.cost), nfev=[int(f1.nfev), int(f2.nfev)], rows=out)


def loops_stage(workers, only=None, frozen=COMPACT, ap230="ap230_g0.7071.json"):
    from concurrent.futures import ProcessPoolExecutor
    info, fade, law, ap, eps0 = q1_tables()
    ap = json.load(open(os.path.join(OUT, ap230)))
    ends = json.load(open(os.path.join(OUT, "ramp.json")))
    fits = json.load(open(os.path.join(Q1, "loops.json")))
    d = json.load(open(os.path.join(OUT, frozen + ".json")))
    v = np.load(os.path.join(OUT, frozen + ".npy"))
    law = q1.fit_laws(fits)
    keys = [k for k in info if only is None or k in only]
    jobs = []
    for k in keys:
        g = info[k]
        caps = [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in files_of(g)]
        v0 = np.median([r["v"] for r in fits if r["file"] in {c["file"] for c in caps}], 0)
        jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
                     v, 1700, caps, ends, v0, d["diff"]))
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_loops_setting, jobs))
    path = os.path.join(OUT, "loops_compact" + ("_smoke" if only else "") + ".json")
    json.dump(res, open(path, "w"), indent=0)
    rep = {(r["key"], r["amp"]): r for r in json.load(open(os.path.join(OUT, "nulls_free.json")))}
    ws = ["50-150", "150-300", "300-600", "600-1000"]

    def ex_(rows, w):
        e = [r["null_comp"][w][1] - rep[(r["key"], r["amp"])]["repeat_comp"][w][1] for r in rows
             if r["null_comp"].get(w) and rep.get((r["key"], r["amp"])) and rep[(r["key"], r["amp"])]["repeat_comp"].get(w)]
        return np.median(e) if e else np.nan
    rows = [r for s_ in res for r in s_["rows"]]
    print(f"\n## Loops refitted with the compact injection frozen ({frozen}); {len(res)} settings; R, compensated; null [excess over repeat]")
    for lp in ("law", "refit", "perloop"):
        rr = [r for r in rows if r["loops"] == lp]
        print(f"- {lp:7s}: " + "; ".join(f"{w}: {np.median([r['null_comp'][w][1] for r in rr if r['null_comp'].get(w)]):.1f} [{ex_(rr, w):+.1f}]" for w in ws))
        for c in ("rt", "hf", "lf"):
            vals = sorted({r[c] for r in rr})
            if len(vals) > 1:
                print(f"    by {c}: " + "; ".join(f"{v_}: " + "/".join(f"{ex_([r for r in rr if r[c] == v_], w):+.1f}" for w in ws[:3]) for v_ in vals))
    if len(res) <= 12:
        for s_ in res:
            print(f"  {s_['key']:12s} " + "  ".join(f"{lp}: " + "/".join(f"{ex_([r for r in s_['rows'] if r['loops'] == lp], w):+.1f}" for w in ws[:3])
                                                for lp in ("law", "refit", "perloop")))


CLA_MODELS = {
    "q1 shared x0+Gx1 (reference)": dict(v=os.path.join(Q1, "v_rt_x0+Gx1.npy"), law=os.path.join(Q1, "split_rt_x0+Gx1.json"), ext=300, through_g=True),
    "free flat injection (27,924)": dict(v="branches_flat_0_ext1700.npy", law="branches_flat_0_ext1700.json", ext=1700),
    "5 taps, diffuser -0.70 (9)": dict(v="tied_tied_g0.445.npy", law=None, ext=1700),
    "5 taps, diffuser -1/sqrt(2) (9)": dict(v=COMPACT + ".npy", law=None, ext=1700, diff=COMPACT + ".json"),
}


def cla_compare(names=None):
    import model_s as m
    fits = json.load(open(os.path.join(Q1, "loops.json")))
    h = np.load(q1.CLA_DEV)[m.K0 - PRE:m.K0 + int(1.9 * FS)]
    yb = sosfiltfilt(q1.LFBAND, h, axis=0)
    n = len(yb) - 600
    ws = W_COMP[:6] + [(400, 600)]
    out = {}
    for name, spec in CLA_MODELS.items():
        if names and name not in names:
            continue
        q1.EXT = spec["ext"]
        if spec["law"]:
            law = json.load(open(spec["law"] if os.path.isabs(spec["law"]) else os.path.join(OUT, spec["law"])))["law"]
        else:
            law = q1.fit_laws(fits)
        law["lf"] = {float(k): x for k, x in law["lf"].items()}
        v = np.load(spec["v"] if os.path.isabs(spec["v"]) else os.path.join(OUT, spec["v"]))
        if spec.get("through_g"):
            G = q1.group_spectra(q1.CLA_LAB, law, q1.grid_n(q1.CLA_LAB["rt"]), through_g=True)
        else:
            diff = json.load(open(os.path.join(OUT, spec["diff"])))["diff"] if spec.get("diff") else None
            G = branch_spectra(q1.CLA_LAB, law, [0.0], 0.0, diff)
        J = q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=q1.CLA_LAB["amp"])])])
        mdl = J.render(v, G, n, q1.CLA_LAB["amp"])
        _, (eg, shg, gain) = q1._clock_job((yb, mdl, 1.0, 98.4, True))
        y = q1.apply_clock(yb, eg, shg)[:n] / gain
        out[name] = dict(eps_ppm=eg, shift=shg, gain_db=float(20 * np.log10(gain)),
                         null_100_10k=q1.nulls(y, mdl, 1.5, ws, narrow=True), null_20_11k=q1.nulls(y, mdl, 1.5, ws))
        r = out[name]
        print(f"{name:34s} clock {eg:6.1f} ppm, shift {shg:+.3f}, gain {r['gain_db']:+.2f} dB | 100 Hz-10 kHz L / R: " + "; ".join(
            f"{w}: {x[0]:.1f}/{x[1]:.1f}" for w, x in r["null_100_10k"].items() if x), flush=True)
    path = os.path.join(OUT, "cla_compare.json")
    old = json.load(open(path)) if os.path.exists(path) else {}
    old.update(out)
    json.dump(old, open(path, "w"), indent=1)
    return out


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
    ap.add_argument("--stage", default="free", choices=["windows", "ramp", "free", "joint", "nulls", "residues", "regress", "branches", "breport", "filters", "tied", "diffuser", "ap230", "loops", "cla", "report"])
    ap.add_argument("--g230", type=float, default=-0.7071067812, help="ap230: the all-pass gain for the step measurement")
    ap.add_argument("--ap230", default=None, help="diffuser: an AP230 step table from the ap230 stage (default: q1's)")
    ap.add_argument("--maxfev", type=int, default=120, help="diffuser: Nelder-Mead evaluations")
    ap.add_argument("--dtap", action="store_true", help="diffuser: also fit one common fractional delay of the loop paths")
    ap.add_argument("--gfix", type=float, default=None, help="diffuser: fix all three all-pass gains (no fit), e.g. -0.70710678")
    ap.add_argument("--mode", default="tied", help="tied: one magnitude; tapgain: one gain per tap; kernelK: one kernel at lags -K..K for every tap")
    ap.add_argument("--g460", default="0.445", help="tied: comma-separated g460 values (0.445 = the post comb cancelled)")
    ap.add_argument("--dense", type=int, default=0, help="tied: free dense coefficients over the first N samples of each excitation")
    ap.add_argument("--score", default="test", choices=["test", "all"], help="tied: score the flat test settings, or all 128")
    ap.add_argument("--tags", default="", help="breport: comma-separated branch run names (files in OUT without .json)")
    ap.add_argument("--alphas", default="0,1", help="branches: per branch an exponent of G_i, or g / H / HL / HH / gHL / gHH")
    ap.add_argument("--ext", type=int, default=300, help="branches: excitation length P_i + ext (q1: 300)")
    ap.add_argument("--subset", default="flat", choices=["flat", "mild", "all"],
                    help="branches: LF 1.0 / HF 1.0 only, all but LF 0.1 and HF 0.2, or all settings")
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
    if a.stage == "branches":
        for sp in a.splits.split(","):
            branches_stage([x if x.isalpha() or x.startswith("g") else float(x) for x in a.alphas.split(",")], a.subset, a.workers, a.iters, a.rounds, a.ext, sp)
    if a.stage == "tied":
        tied_stage(a.mode, [float(x) for x in a.g460.split(",")], a.workers, a.iters, a.rounds, a.ext if a.ext != 300 else 1700, a.score, a.dense,
                   a.gfix, a.ap230)
    if a.stage == "diffuser":
        diffuser_stage(a.workers, a.score, a.maxfev, gfix=a.gfix, ap230=a.ap230, dtap=a.dtap)
    if a.stage == "cla":
        cla_compare([a.only] if a.only else None)
    if a.stage == "loops":
        loops_stage(a.workers, only)
    if a.stage == "ap230":
        ap230_stage(a.g230, only)
    if a.stage == "filters":
        filters_stage(a.workers, only)
    if a.stage == "breport":
        branch_report(a.tags.split(","))
    if a.stage == "residues":
        residues_stage(a.workers, only, a.source)
    if a.stage == "regress":
        regress_stage(os.path.join(OUT, "residues" + ("" if a.source == "free" else "_" + a.source) + ("_smoke" if only else "") + ".npz"))
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
