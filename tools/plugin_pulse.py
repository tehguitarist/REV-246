"""The pulse check (tools/pulse_check.py, G1-G6) with the PLUGIN's renders in the research model's place, programs S, M and
L, against the Tomkin 96 kHz captures: plugin vs capture beside model vs capture beside capture vs capture (the tolerance).
Results: analysis/plugin_pulse/, write-up analysis/plugin_check/REPORT_chainfix.md §8.

    .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/plugin_pulse.py --stage s  [--only k,..] [--workers 2]
    .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/plugin_pulse.py --stage ml [--only M_k,..] [--workers 2]
    .venv/bin/python tools/plugin_pulse.py --stage report

The renders are plugin_check's corrected-bridge set (PLUGIN_CHECK_TAG=fix, analysis/plugin_check/renders_fix/96000/,
fixed anchor align_fix.json; M / L at theta 0.13 rad rms, modulation on). Nothing is re-rendered.

  s    per capture: the render through the capture's own path (plugin_check.plugin_q1: chain3 deconvolution, device
       grid, the fixed anchor, per-channel tap-0 shift, q1's LFBAND), fitted as plugin_check's s stage (refit_shift:
       clock, shift and gain, 3 numbers), then pulse_check's prep / compare / nulls exactly as its s stage (capture
       fade-compensated up to comp_end). Every capture of the setting, as pulse_check (DEFECTS flagged in the report).
  ml   per AMPL group: the render in ml_carryover's frame (plugin_check.ml_parts: cluster part = render - AMPL-0 render,
       reverb part = AMPL-0 render), fitted as plugin_check's ml stage (ml_lines.fit_warp: shift, clock and two gains per
       channel on cluster + gap + first pass, the same freedoms as pulse_check's static model), then pulse_check's
       compare (EDCs from R onset + 20 ms, coherence / fine structure over the 250 ms after it) and _ml_nulls.
  report  the plugin's rows joined with pulse_check's own model rows and capture pairs (analysis/pulse/s, analysis/pulse/ml:
       the same metrics, the same captures) -> analysis/plugin_pulse/report.log, summary.json.

Columns in the report: plugin / model [capture vs capture]. For M / L the plugin is time-variant (modulation on) and the
captures are too, so plugin-vs-capture should be compared with capture-vs-capture, not with the static model: metrics of
one IR compared afterwards (decay times, EDT, EDC shape, C50 / C80 / D50 / Ts, tonal balance, echo density, coherence
values, fine-structure sd, mixing time, IACC, L/R level) are meaningful; waveform-matching ones (nulls past the first pass,
fine-structure correlation, narrow-band envelope and EDR per 1/24 octave) measure decorrelation and are read against the
capture pairs' values.
"""
import argparse, glob, json, os, sys, time

os.environ.setdefault("PLUGIN_CHECK_TAG", "fix")
for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import plugin_check as pk  # noqa: E402  (sets EMT_KERNEL=chain3, EMT_PMAX=0.8475)
import pulse_check as pc  # noqa: E402
import numpy as np  # noqa: E402

assert pk.TAG == "fix", "PLUGIN_CHECK_TAG must be 'fix' (renders_fix/, align_fix.json)"
OUT = os.path.join(ROOT, "analysis", "plugin_pulse")
REF = os.path.join(ROOT, "analysis", "pulse")                  # pulse_check's full run: model rows and capture pairs
q2, q1 = pc.q2, pc.q1
FS, PRE = q1.FS, q1.PRE
OCT, THIRD, FS_BANDS = pc.OCT, pc.THIRD, pc.FS_BANDS


# ---------------------------------------------------------------- program S
def _s_job(a):
    key, caps, ends, al, path = a
    t = time.time()
    q1.EXT = 1700
    rows = []
    for c in caps:
        n = len(q2._band(c["file"]))
        p = pk.plugin_q1(pk.rpath(pk.HOST, "S", key, c["amp"]), n, al)
        y, g, te, fit = pk.refit_shift(p, c, ends, n)
        i1 = PRE + int(te * FS)
        Y, P = pc.prep(y, g), pc.prep(p)
        rows.append(dict(amp=c["amp"], file=c["file"], t_comp=te, fit=fit, null=q1.nulls(y, p * g, te, q2.W_COMP, narrow=True),
                         null_oct=q2.octave_nulls(y, p * g, te), m=pc.compare(P, Y, i1)))
        del Y, P
    json.dump(dict(key=key, rows=rows), open(path, "w"))
    nl = rows[-1]["null"]
    print(f"S {key}: {time.time() - t:.0f} s, R null (AMPL {rows[-1]['amp']}) " + " ".join(f"{w}:{x[1]:.1f}" for w, x in nl.items() if x), flush=True)
    return key


def s_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, "s")
    os.makedirs(d, exist_ok=True)
    al = json.load(open(os.path.join(pk.OUT, "align_fix.json")))
    info, fade, _, _, eps0 = q2.q1_tables()
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    jobs = []
    for k in sorted(info):
        path = os.path.join(d, k + ".json")
        if (only and k not in only) or os.path.exists(path):
            continue
        caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(info[k])]
        jobs.append((k, caps, ends, al, path))
    print(f"S: {len(jobs)} settings to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_s_job, jobs):
            pass


# ---------------------------------------------------------------- programs M and L
def _ml_job(a):
    prog, key, path = a
    import ml_carryover as cc
    import ml_lines as mll
    t = time.time()
    K0 = cc.K0
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, _, amps = cc.captures(prog, key, n)
    n = len(caps[amps[-1]])
    C, R, s = pk.ml_parts(prog, key, amps, n, False)
    t0 = (cc.ONSET[prog][1] + 20) / FS
    fw = (t0, t0 + 0.25)
    rows = []
    for amp in amps:
        y = caps[amp]
        i1 = min(mll.usable_end(y, g[amp], prog))
        sl, nf = mll.fit_window(prog, amp > 0)
        mdl, u, gs, _ = mll.fit_warp(y, [C[amp], R] if amp > 0 else [R], sl, nf)
        Y, P = pc.prep(y[:i1], dtype=np.float32), pc.prep(mdl[:i1], dtype=np.float32)
        rows.append(dict(amp=amp, t_end=(i1 - K0) / FS, fit=dict(shift=u[0], ppm=u[1], gains=gs),
                         null=pc._ml_nulls(prog, y, mdl, i1), m=pc.compare(P, Y, i1, K0, t0, fw)))
        del Y, P
    json.dump(dict(key=key, prog=prog, shift=s, rows=rows), open(path, "w"))
    print(f"{prog} {key}: {time.time() - t:.0f} s, recirc1 plugin " + " ".join(f"{r['null']['recirc1']:.1f}" for r in rows if r["null"]["recirc1"] is not None), flush=True)
    return key


def ml_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    import ml_carryover as cc
    d = os.path.join(OUT, "ml")
    os.makedirs(d, exist_ok=True)
    jobs = []
    for prog in ("M", "L"):
        for k in sorted(cc.groups(prog)):
            path = os.path.join(d, f"{prog}_{k}.json")
            if (only and f"{prog}_{k}" not in only) or os.path.exists(path):
                continue
            jobs.append((prog, k, path))
    # longest groups first, so the pool's tail is short ones
    jobs.sort(key=lambda j: -float(j[1].split("_")[1]) * (3 if j[0] == "L" else 1))
    print(f"ML: {len(jobs)} groups to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_ml_job, jobs):
            pass


# ---------------------------------------------------------------- report
def _md(rr, k, f=lambda x: x):
    x = [f(s[k]) for s in rr if s.get(k) is not None and np.isfinite(s[k])]
    return float(np.median(x)) if x else np.nan


def _three(P, M, C, k, f=lambda x: x, sc=1.0, fmt="{:.1f}"):
    """'plugin / model [capture]' for one summarise key; also the numbers."""
    v = [sc * _md(rr, k, f) for rr in (P, M, C)]
    return f"{fmt.format(v[0])} / {fmt.format(v[1])} [{fmt.format(v[2])}]", v


def _within(sets, key, floor=0.03):
    """Share of captures whose |error| (plugin and model) is within its setting's largest capture-pair |difference| (max
    with floor). sets: [(plugin summaries, model summaries, pair summaries)] per setting / group."""
    ok, n = [0, 0], 0
    for P, M, C in sets:
        rep = [abs(c[key]) for c in C if c.get(key) is not None]
        if not rep:
            continue
        tol = max(max(rep), floor)
        for j, rr in enumerate((P, M)):
            for s in rr:
                if s.get(key) is not None:
                    ok[j] += abs(s[key]) <= tol
        n += sum(s.get(key) is not None for s in P)
    return ok, n


def _g6(P, M, C, indent):
    """pulse_check.g6_lines for plugin and for model (each with the capture pairs), item by item."""
    lp, lm = pc.g6_lines(P, C, indent=""), pc.g6_lines(M, C, indent="")
    return [f"{indent}plugin {a}\n{indent[:-2]}    model  {b}" for a, b in zip(lp, lm)]


def _decay_tables(out, groups, summ, key, title, f=abs, sc=100, fmt="{:.1f}"):
    out.append(f"\n### {title}")
    out.append("| Group | " + " | ".join(str(b) for b in OCT) + " |\n|---|" + "---|" * len(OCT))
    for g, (P, M, C) in groups:
        cells = []
        for b in OCT:
            s, v = _three(P, M, C, key.format(b), f, sc, fmt)
            cells.append(s)
            summ.setdefault(g, {})[key.format(b) + ("_abs" if f is abs else "")] = v
        out.append(f"| {g} | " + " | ".join(cells) + " |")


def report_s(summ):
    from defects import DEFECTS
    files = sorted(glob.glob(os.path.join(OUT, "s", "*.json")))
    if not files:
        return []
    S = []                                          # (lab, plugin rows, model rows, cross pairs, in-take pairs)
    for fn in files:
        dp = json.load(open(fn))
        dm = json.load(open(os.path.join(REF, "s", dp["key"] + ".json")))
        S.append((dm, {r["amp"]: r for r in dp["rows"]}, {r["amp"]: r for r in dm["rows"]},
                  [p for p in dm["pairs"] if p["cross"]], [p for p in dm["pairs"] if not p["cross"]]))
    nP = sum(len(x[1]) for x in S)
    out = [f"## Program S: {len(S)} settings, {nP} captures; plugin vs capture / model vs capture [cross-take capture vs capture] "
           "(the plugin: renders_fix, 3-number refit per capture; the model: pulse_check's held-out baseline, analysis/pulse/s)"]
    groups_def = [("All", lambda d, a: True), ("All, clean", lambda d, a: (d["key"], a) not in DEFECTS),
                  ("RT <= 0.6", lambda d, a: d["rt"] <= 0.6), ("RT 0.8-1.4", lambda d, a: 0.8 <= d["rt"] <= 1.4),
                  ("RT >= 1.8", lambda d, a: d["rt"] >= 1.8), ("LF 0.1", lambda d, a: d["lf"] == 0.1), ("LF 2.0", lambda d, a: d["lf"] == 2.0),
                  ("HF 0.2", lambda d, a: d["hf"] == 0.2), ("HF 1.0", lambda d, a: d["hf"] == 1.0)]
    groups = []
    for g, fn in groups_def:
        P = [pc.summarise(r["m"]) for d, rp, rm, cr, _ in S for a, r in rp.items() if fn(d, a)]
        M = [pc.summarise(r["m"]) for d, rp, rm, cr, _ in S for a, r in rm.items() if a in rp and fn(d, a)]
        C = [pc.summarise(p["m"]) for d, rp, rm, cr, _ in S for p in cr if all(fn(d, a) for a in p["pair"])]
        groups.append((g, (P, M, C)))
    G = dict(groups)
    sm = summ.setdefault("S", {})
    _decay_tables(out, groups, sm, "oct_t{}", "G2: octave decay time from the EDC (-5 to -25 dB), |ratio - 1| %, R; median")
    _decay_tables(out, groups, sm, "oct_t{}", "G2: signed octave decay-time ratio - 1 (vs capture) %, R; median (bias)", f=lambda x: x, fmt="{:+.1f}")
    _decay_tables(out, groups, sm, "oct_late{}", "G2: late decay (-25 to -40 dB) |ratio - 1| %, R; median")
    _decay_tables(out, groups, sm, "oct_edt{}", "G2: EDT |ratio - 1| %, R; median")
    _decay_tables(out, groups, sm, "oct_dev{}", "G2: largest EDC difference to -30 dB (dB), R; median", f=lambda x: x, sc=1, fmt="{:.2f}")
    _decay_tables(out, groups, sm, "oct_t{}_L", "G2: octave decay |ratio - 1| %, L; median")
    # share of captures within the setting's largest cross-take difference
    per = []
    for d, rp, rm, cr, _ in S:
        per.append(([pc.summarise(r["m"]) for r in rp.values()], [pc.summarise(rm[a]["m"]) for a in rp if a in rm], [pc.summarise(p["m"]) for p in cr]))
    out.append("\n### G2 / G3: share of captures whose |error| is within the setting's largest cross-take |difference| (floor 3 %), plugin / model")
    for band, fs_ in (("oct", OCT), ("third", THIRD)):
        cells = []
        for f in fs_:
            (op, om), n = _within(per, f"{band}_t{f}")
            cells.append(f"{f}: {op}/{om} of {n}")
            sm.setdefault("within", {})[f"{band}_t{f}"] = [op, om, n]
        out.append(f"- {band} decay: " + "; ".join(cells))
    P, M, C = G["All"]
    out.append("\n### G3: third-octave decay time, signed ratio - 1 (%), R, all: plugin / model [repeat]")
    out.append(" | ".join(f"{f}: {_three(P, M, C, f'third_t{f}', sc=100, fmt='{:+.1f}')[0]}" for f in THIRD))
    out.append("- |ratio - 1| %: " + " | ".join(f"{f}: {_three(P, M, C, f'third_t{f}', abs, 100)[0]}" for f in THIRD))
    for f in THIRD:
        sm.setdefault("All", {})[f"third_t{f}_abs"] = _three(P, M, C, f"third_t{f}", abs, 100)[1]
    out.append("\n### G4: echo density mean |difference| 20-300 ms (R), and L/R coherence (vs capture) per third octave, 50-300 ms")
    for g, (P_, M_, C_) in groups:
        s, v = _three(P_, M_, C_, "ned", fmt="{:.3f}")
        sm.setdefault(g, {})["ned"] = v
        out.append(f"- {g}: NED {s}")
    out.append("- coherence diff, all: " + "; ".join(f"{f}: {_three(P, M, C, f'coh_{f}', fmt='{:+.3f}')[0]} (capture {_md(P, f'cohv_{f}'):.3f})" for f in THIRD))
    out.append("\n### G5: spectral fine-structure correlation 50-300 ms, R (plugin L)")
    for g, (P_, M_, C_) in groups:
        cells = []
        for lo, hi in FS_BANDS:
            w = f"{lo}-{hi}"
            s, v = _three(P_, M_, C_, f"fine_{w}", fmt="{:.3f}")
            sm.setdefault(g, {})[f"fine_{w}"] = v
            cells.append(f"{w}: {s} ({_md(P_, f'fine_{w}_L'):.3f})")
        out.append(f"- {g}: " + "; ".join(cells))
    out.append("\n### G6: reverb-specific metrics, R unless stated: plugin - capture / model - capture [cross-take capture - capture]")
    for g in ("All", "All, clean", "RT <= 0.6", "RT >= 1.8"):
        out.append(f"- {g}:")
        out += _g6(*G[g], indent="  - ")
    for g in ("All",):
        P_, M_, C_ = G[g]
        for k in [f"rv_{x}_{f}" for x in ("c50", "c80", "d50", "ts") for f in OCT] + [f"rv_lev_e_{f}" for f in THIRD] + [f"rv_lev_l_{f}" for f in THIRD] + \
                 ["rv_edr_p90", "rv_edr_max", "rv_nb_p90", "rv_nb_max", "rv_mix_ms", "rv_iacc_e", "rv_iacc_l", "rv_lr_e", "rv_lr_l"] + [f"rv_fsd_{lo}-{hi}" for lo, hi in FS_BANDS]:
            sm[g][k + "_abs"] = [_md(rr, k, abs) for rr in (P_, M_, C_)]
            sm[g][k] = [_md(rr, k) for rr in (P_, M_, C_)]
    out.append("\n### G1: nulls 100 Hz-10 kHz, R, median plugin / model [cross-take repeat]")
    for g, fn in groups_def:
        cells = []
        for w in ("0-50", "50-150", "150-300", "300-600", "600-1000"):
            p_ = [r["null"][w][1] for d, rp, rm, cr, _ in S for a, r in rp.items() if fn(d, a) and r["null"].get(w)]
            m_ = [r["null"][w][1] for d, rp, rm, cr, _ in S for a, r in rm.items() if a in rp and fn(d, a) and r["null"].get(w)]
            c_ = [p["null"][w][1] for d, rp, rm, cr, _ in S for p in cr if all(fn(d, a) for a in p["pair"]) and p["null"].get(w)]
            v = [float(np.median(x)) if x else np.nan for x in (p_, m_, c_)]
            sm.setdefault(g, {})[f"null_{w}"] = v
            cells.append(f"{w}: {v[0]:.1f} / {v[1]:.1f} [{v[2]:.1f}]")
        out.append(f"- {g}: " + "; ".join(cells))
    out.append("\n### G1: octave nulls (R), all, median plugin / model per window 50-150 / 150-300 / 300-600")
    for f in OCT:
        cells = []
        for j in range(3):
            p_ = [r["null_oct"][str(f)][j] for d, rp, rm, cr, _ in S for r in rp.values() if r["null_oct"][str(f)][j] is not None]
            m_ = [r["null_oct"][str(f)][j] for d, rp, rm, cr, _ in S for a, r in rm.items() if a in rp and r["null_oct"][str(f)][j] is not None]
            cells.append(f"{np.median(p_):.1f} / {np.median(m_):.1f}")
        out.append(f"- {f}: " + " | ".join(cells))
    return out


def report_ml(summ):
    from defects import ML_SLIPPED
    files = sorted(glob.glob(os.path.join(OUT, "ml", "*.json")))
    if not files:
        return []
    out = ["\n## Programs M and L: plugin (time-variant, theta 0.13, renders_fix) vs capture / static model (theta 0) vs capture "
           "[capture vs capture, every pair of the group], medians; slipped captures and their pairs left out"]
    D = []
    for fn in files:
        dp = json.load(open(fn))
        dm = json.load(open(os.path.join(REF, "ml", os.path.basename(fn))))
        bad = lambda a: (dp["prog"], dp["key"], a) in ML_SLIPPED
        D.append((dm, {r["amp"]: r for r in dp["rows"] if not bad(r["amp"])}, {r["amp"]: r for r in dm["rows"] if not bad(r["amp"])},
                  [p for p in dm["pairs"] if not any(bad(a) for a in p["pair"])]))
    for prog in ("M", "L"):
        Gp = [x for x in D if x[0]["prog"] == prog]
        if not Gp:
            continue
        sm = summ.setdefault(prog, {})
        sel = lambda fn: ([pc.summarise(r["m"]) for d, rp, rm, cr in Gp if fn(d) for r in rp.values()],
                          [pc.summarise(r["m"]) for d, rp, rm, cr in Gp if fn(d) for a, r in rm.items() if a in rp],
                          [pc.summarise(p["m"]) for d, rp, rm, cr in Gp if fn(d) for p in cr])
        P, M, C = sel(lambda d: True)
        out.append(f"\n### {prog}: {len(Gp)} groups, {len(P)} captures, {len(C)} capture pairs")
        groups = [("All", (P, M, C))] + [(g, sel(fn)) for g, fn in (
            ("RT short", lambda d: d["rt"] <= {"M": 1.0, "L": 3.0}[prog]), ("RT long", lambda d: d["rt"] >= {"M": 3.0, "L": 8.0}[prog]),
            ("LF 0.1", lambda d: d["lf"] == 0.1), ("HF <= 0.2", lambda d: d["hf"] <= 0.2))]
        _decay_tables(out, groups, sm, "oct_t{}", f"{prog} G2: octave decay (-5/-25 dB) |ratio - 1| %, R")
        _decay_tables(out, groups, sm, "oct_t{}", f"{prog} G2: octave decay signed ratio - 1 %, R", f=lambda x: x, fmt="{:+.1f}")
        _decay_tables(out, groups[:1], sm, "oct_late{}", f"{prog} G2: late decay (-25/-40) |ratio - 1| %, R")
        _decay_tables(out, groups[:1], sm, "oct_edt{}", f"{prog} G2: EDT |ratio - 1| %, R")
        _decay_tables(out, groups[:1], sm, "oct_dev{}", f"{prog} G2: largest EDC difference to -30 dB (dB), R", f=lambda x: x, sc=1, fmt="{:.2f}")
        _decay_tables(out, groups[:1], sm, "oct_t{}_L", f"{prog} G2: octave decay |ratio - 1| %, L")
        per = [([pc.summarise(r["m"]) for r in rp.values()], [pc.summarise(rm[a]["m"]) for a in rp if a in rm], [pc.summarise(p["m"]) for p in cr]) for d, rp, rm, cr in Gp]
        cells = []
        for f in OCT:
            (op, om), n = _within(per, f"oct_t{f}")
            cells.append(f"{f}: {op}/{om} of {n}")
            sm.setdefault("within", {})[f"oct_t{f}"] = [op, om, n]
        out.append("- octave decay within the group's largest capture-pair |difference| (floor 3 %), plugin / model: " + "; ".join(cells))
        out.append(f"- G3 third-octave decay signed %, R: " + " ".join(f"{f}: {_three(P, M, C, f'third_t{f}', sc=100, fmt='{:+.1f}')[0]}" for f in THIRD))
        out.append(f"- G3 third-octave decay |ratio - 1| %, R: " + " ".join(f"{f}: {_three(P, M, C, f'third_t{f}', abs, 100)[0]}" for f in THIRD))
        for f in THIRD:
            sm["All"][f"third_t{f}_abs"] = _three(P, M, C, f"third_t{f}", abs, 100)[1]
        s, v = _three(P, M, C, "ned", fmt="{:.3f}")
        sm["All"]["ned"] = v
        out.append(f"- G4 echo density |diff|: {s}")
        out.append("- G4 L/R coherence vs capture (capture value): " + " ".join(f"{f}: {_three(P, M, C, f'coh_{f}', fmt='{:+.3f}')[0]} ({_md(P, f'cohv_{f}'):.2f})" for f in THIRD[::2]))
        for f in THIRD:
            sm["All"][f"coh_{f}_abs"] = _three(P, M, C, f"coh_{f}", abs)[1]
        cells = []
        for lo, hi in FS_BANDS:
            w = f"{lo}-{hi}"
            s, v = _three(P, M, C, f"fine_{w}", fmt="{:.3f}")
            sm["All"][f"fine_{w}"] = v
            cells.append(f"{w}: {s} (plugin L {_md(P, f'fine_{w}_L'):.3f})")
        out.append("- G5 fine-structure correlation R: " + "; ".join(cells))
        out.append(f"- G6, reverb-specific, all groups:")
        out += _g6(P, M, C, indent="  - ")
        for k in [f"rv_{x}_{f}" for x in ("c50", "c80", "d50", "ts") for f in OCT] + [f"rv_lev_e_{f}" for f in THIRD] + [f"rv_lev_l_{f}" for f in THIRD] + \
                 ["rv_edr_p90", "rv_edr_max", "rv_nb_p90", "rv_nb_max", "rv_mix_ms", "rv_iacc_e", "rv_iacc_l", "rv_lr_e", "rv_lr_l"] + [f"rv_fsd_{lo}-{hi}" for lo, hi in FS_BANDS]:
            sm["All"][k + "_abs"] = [_md(rr, k, abs) for rr in (P, M, C)]
            sm["All"][k] = [_md(rr, k) for rr in (P, M, C)]
        ws = ["first_pass", "recirc1", "recirc2", "recirc3"] + [f"{a}-{b}" for a, b in pc.ML_TW]
        cells = []
        for w in ws:
            x = [[r["null"][w] for d, rp, rm, cr in Gp for r in rp.values() if r["null"].get(w) is not None],
                 [r["null"][w] for d, rp, rm, cr in Gp for a, r in rm.items() if a in rp and r["null"].get(w) is not None],
                 [p["null"][w] for d, rp, rm, cr in Gp for p in cr if p["null"].get(w) is not None]]
            v = [float(np.median(z)) if z else np.nan for z in x]
            sm["All"][f"null_{w}"] = v
            cells.append(f"{w}: {v[0]:.1f} / {v[1]:.1f} [{v[2]:.1f}]")
        out.append("- G1 decorrelation: null per window, plugin / static model [capture pair]: " + "; ".join(cells))
    return out


def report_stage():
    summ = {}
    out = report_s(summ) + report_ml(summ)
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(out) + "\n")
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=lambda x: None if not np.isfinite(x) else x)
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["s", "ml", "report"], required=True)
    ap.add_argument("--only", default=None, help="s: setting keys lf_rt_hf; ml: prog_lf_rt_hf; comma-separated")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    only = a.only.split(",") if a.only else None
    if a.stage == "s":
        s_stage(a.workers, only)
    elif a.stage == "ml":
        ml_stage(a.workers, only)
    else:
        report_stage()


if __name__ == "__main__":
    main()
