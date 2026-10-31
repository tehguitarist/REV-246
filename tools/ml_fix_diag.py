"""Why the plugin's M / L renders miss the capture tolerance below 500 Hz (analysis/plugin_check/REPORT_chainfix.md §8):
phase 1, diagnosis only. Results: analysis/ml_fix/ (diag.json, diag.log, per-group JSON).

    .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/ml_fix_diag.py --stage render [--only M_k,..]
    .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/ml_fix_diag.py --stage ml     [--only ..] [--workers 1]
    .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/ml_fix_diag.py --stage port   [--only ..]
    .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/ml_fix_diag.py --stage sweep  [--only ..]
    .venv/bin/python tools/ml_fix_diag.py --stage report

  render  REV246Render (plugin/build-render): per setting, impulse renders at 3 start phases of theta (pre = 0, 1/3, 2/3
          of the 0.83 Hz period: 0, 0.4016, 0.8032 s of silence) with modulation on, and one with modulation off (mod=0),
          every AMPL of the group plus AMPL 0 (host 96 kHz, the whole engine, as plugin_check's ml set); and device-domain
          AMPL-0 renders (the core alone) at the same 3 phases and static, for the port check.
          -> analysis/ml_fix/renders/{96000,dev}/
  ml      per group, per capture: each variant (p0, p1, p2, static) through plugin_pulse's path (plugin_check.ml_parts'
          frame, ml_lines.fit_warp on cluster + gap + first pass, pulse_check.compare / _ml_nulls); plugin vs capture per
          variant, plugin vs plugin (p_i vs p_j, in the capture's frame and window), modulated vs static (p_k vs static).
  port    the Python reference (ml_variance.core, the read-out form) for the device-domain renders: static vs mod=0, and
          theta(t) = A sin(w (n + n0 + 1 + d)) vs each phase render, d (theta's time origin, samples) scanned.
  sweep   the plugin's existing sweep renders (renders_fix/96000/*_sweep0-5, plugin_check's sim set, the 8 SIM_KEYS
          groups) deconvolved as plugin_check._sim_job: decay per octave of the sweep-measured IR against the static
          render's, and against the impulse renders' (does the measurement, not the plugin, move the low-band decay?).
  report  -> analysis/ml_fix/diag.json, diag.log.
"""
import argparse, glob, json, os, subprocess, sys, time

os.environ.setdefault("PLUGIN_CHECK_TAG", "fix")
for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import plugin_check as pk  # noqa: E402
import pulse_check as pc  # noqa: E402
import numpy as np  # noqa: E402
import soundfile as sf  # noqa: E402
from scipy.signal import sosfiltfilt  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_fix")
RND = os.path.join(OUT, "renders")
FS = pc.FS
FSD_FILE = 26374
PERIOD = 1 / 0.83
PRE = [0.0, round(PERIOD / 3, 4), round(2 * PERIOD / 3, 4)]
VARS = ["p0", "p1", "p2", "static"]
KEYS = {"M": ["2.0_0.8_0.9", "0.1_0.8_0.9", "0.1_1.8_0.4", "0.7_3.0_0.2", "0.1_4.0_0.4", "1.2_5.0_0.9"],
        "L": ["2.0_2.6_0.5", "2.0_3.6_0.3", "1.2_4.5_0.8", "0.1_8.0_0.3", "0.1_10.0_0.3", "1.2_20.0_1.0"]}
SIM_KEYS = {"M": ["1.2_5.0_0.9", "2.0_0.8_0.9", "0.7_3.0_0.2", "0.1_1.8_0.4"], "L": ["1.2_20.0_1.0", "1.2_4.5_0.8", "2.0_2.6_0.5", "0.1_10.0_0.3"]}   # = ml_variance's
def tail(key):
    """s of sweep tail rendered and deconvolved: 6 s, or 0.6 RT for the long L settings (RT 20: 12 s)."""
    return max(6.0, min(14.0, 0.6 * float(key.split("_")[1])))


DEV_LEN = 6.0                                   # s, device renders (port check)


def groups(only=None):
    g = [(p, k) for p, ks in KEYS.items() for k in ks]
    return [x for x in g if not only or f"{x[0]}_{x[1]}" in only]


def hpath(prog, key, amp, v):
    return os.path.join(RND, "96000", f"{prog}_{key}_a{amp:g}_{v}.wav")


def spath(prog, key, i):
    """sweep render i (pk.SIM_PRE[i]); i = len(SIM_PRE): modulation off."""
    return os.path.join(RND, "96000", f"{prog}_{key}_a0_sweep{i if i < len(pk.SIM_PRE) else 'st'}.wav")


def dpath(prog, key, v):
    return os.path.join(RND, "dev", f"{prog}_{key}_a0_{v}.wav")


# ---------------------------------------------------------------- renders
def render_stage(only):
    L = []
    for prog, key in groups(only):
        g = pk.ml_groups(prog)[key]
        lf, rt, hf = (float(x) for x in key.split("_"))
        for amp in sorted(set(g["amps"]) | {0}):
            for i, v in enumerate(VARS):
                L.append(pk.render_line(hpath(prog, key, amp, v), pk.HOST, prog, lf, rt, hf, amp, g["len"],
                                        mod=0 if v == "static" else 1, pre=0.0 if v == "static" else PRE[i]))
        if key in SIM_KEYS[prog]:                     # sweeps (plugin_check's sim set, rendered again by this build) + one static
            sw = os.path.join(RND, "sweep_96k.wav")
            if not os.path.exists(sw):
                os.makedirs(RND, exist_ok=True)
                sf.write(sw, (0.25 * pk.sweep_96k()[0]).astype(np.float32), pk.HOST, subtype="FLOAT")
            T = len(pk.sweep_96k()[0]) / pk.HOST + tail(key)
            for i, pre in enumerate(pk.SIM_PRE + [0.0]):
                L.append(pk.render_line(spath(prog, key, i), pk.HOST, prog, lf, rt, hf, 0, T, mod=int(i < len(pk.SIM_PRE)), pre=pre, inp=sw))
        for i, v in enumerate(VARS):
            L.append(pk.render_line(dpath(prog, key, v), pk.HOST, prog, lf, rt, hf, 0, min(g["len"], DEV_LEN),
                                    mod=0 if v == "static" else 1, pre=0.0 if v == "static" else PRE[i], domain="device"))
    todo = [l for l in L if not os.path.exists(l.split()[0][4:])]
    print(f"render: {len(todo)} of {len(L)} to do", flush=True)
    if not todo:
        return
    for l in todo:
        os.makedirs(os.path.dirname(l.split()[0][4:]), exist_ok=True)
    lst = os.path.join(OUT, "render_list.txt")
    open(lst, "w").write("\n".join(todo) + "\n")
    t = time.time()
    r = subprocess.run([pk.RENDER_BIN, lst], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    print(f"render: done in {time.time() - t:.0f} s; " + "; ".join(sorted({l.split('latency=')[1].split()[0] for l in r.stdout.splitlines() if 'latency=' in l})), flush=True)


# ---------------------------------------------------------------- the audit's metrics on IR pairs
def parts(prog, key, amps, n, v):
    """plugin_check.ml_parts for variant v of this tool's renders."""
    import ml_carryover as cc
    import ml_structure as ml
    anchor = pk.onset_anchor(pk.load_render(hpath(prog, key, max(amps), v)))
    h = {}
    for a in sorted(set(amps) | {0}):
        x = sosfiltfilt(ml.BAND, pk.to_device(pk.load_render(hpath(prog, key, a, v)), anchor), axis=0)
        h[a] = np.vstack([x, np.zeros((max(0, n - len(x)), 2))])[:n]
    s = cc.K0 - cc.tap0(h[max(amps)])
    h = {a: cc.fshift(x, s) for a, x in h.items()}
    return {a: h[a] - h[0] for a in amps}, h[0], float(s)


KEEP = [f"oct_t{f}" for f in pc.OCT] + [f"oct_edt{f}" for f in pc.OCT] + [f"oct_late{f}" for f in pc.OCT] + \
       [f"oct_dev{f}" for f in pc.OCT] + [f"oct_t{f}_L" for f in pc.OCT] + [f"third_t{f}" for f in pc.THIRD if f <= 630] + \
       ["ned", "rv_mix_ms", "rvv_mix_ms"]


def slim(m):
    m = json.loads(json.dumps(m))                 # summarise expects the JSON form (string band keys)
    s = pc.summarise(m)
    out = {k: s.get(k) for k in KEEP}
    out.update({k: v for k, v in s.items() if k.startswith(("rv_", "rvv_", "coh_", "fine_")) and k not in out})   # G4-G6 incl. level metrics
    for f in pc.OCT:                                # the absolute decay times (s) of both IRs, R: the bias test needs them
        b = m["oct"][str(f)][1]["t"]
        out[f"T{f}"] = b
    out["ned_def"] = ned_def(m)
    return out


def ned_def(m):
    """A continuous stand-in for the mixing time: the echo-density deficit sum(max(0, 1 - NED)) dt over the first second
    (NED every 20 ms, pulse_check.rv_metrics' ned1), ratio A / B - 1. The mixing time (first 50 ms-smoothed NED >= 0.9)
    flips by 100-200 ms where NED plateaus near 0.9 (L); the deficit moves smoothly."""
    try:
        a, b = (np.array(m["rv"][x]["ned1"], float) for x in ("a", "b"))
    except (KeyError, TypeError):
        return None
    n = min(len(a), len(b))
    if n < 10:
        return None
    da, db = np.sum(np.maximum(0, 1 - a[:n])), np.sum(np.maximum(0, 1 - b[:n]))
    return float(da / db - 1) if db > 0 else None


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
    V = {v: parts(prog, key, amps, n, v) for v in VARS}
    t0 = (cc.ONSET[prog][1] + 20) / FS
    fw = (t0, t0 + 0.25)
    rows = []
    for amp in amps:
        y = caps[amp]
        i1 = min(mll.usable_end(y, g[amp], prog))
        sl, nf = mll.fit_window(prog, amp > 0)
        Y = pc.prep(y[:i1], dtype=np.float32)
        P, mdl, row = {}, {}, dict(amp=amp, t_end=(i1 - K0) / FS, vs_cap={}, null_cap={}, pairs={}, null_pairs={}, vs_static={}, null_static={})
        for v in VARS:
            C, R, _ = V[v]
            mdl[v], u, gs, _ = mll.fit_warp(y, [C[amp], R] if amp > 0 else [R], sl, nf)
            P[v] = pc.prep(mdl[v][:i1], dtype=np.float32)
            row["vs_cap"][v] = slim(pc.compare(P[v], Y, i1, K0, t0, fw))
            row["null_cap"][v] = pc._ml_nulls(prog, y, mdl[v], i1)
        for i, a_ in enumerate(VARS[:3]):
            for b_ in VARS[i + 1:3]:
                row["pairs"][f"{a_}-{b_}"] = slim(pc.compare(P[a_], P[b_], i1, K0, t0, fw))
                row["null_pairs"][f"{a_}-{b_}"] = pc._ml_nulls(prog, mdl[a_], mdl[b_], i1)
            row["vs_static"][a_] = slim(pc.compare(P[a_], P["static"], i1, K0, t0, fw))
            row["null_static"][a_] = pc._ml_nulls(prog, mdl[a_], mdl["static"], i1)
        rows.append(row)
        del P, Y
    lf, rt, hf = (float(x) for x in key.split("_"))
    json.dump(dict(key=key, prog=prog, lf=lf, rt=rt, hf=hf, shift={v: V[v][2] for v in VARS}, rows=rows), open(path, "w"))
    r = rows[-1]
    f = 125 if prog == "M" else 250
    print(f"{prog} {key}: {time.time() - t:.0f} s; {f} Hz decay vs cap " + " ".join(f"{v} {100 * (r['vs_cap'][v][f'oct_t{f}'] or np.nan):+.1f}" for v in VARS)
          + "; pairs " + " ".join(f"{100 * (x[f'oct_t{f}'] or np.nan):+.1f}" for x in r["pairs"].values()), flush=True)
    return key


def ml_stage(workers, only):
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, "ml")
    os.makedirs(d, exist_ok=True)
    jobs = [(p, k, os.path.join(d, f"{p}_{k}.json")) for p, k in groups(only)]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    jobs.sort(key=lambda j: -float(j[1].split("_")[1]) * (3 if j[0] == "L" else 1))
    print(f"ml: {len(jobs)} groups to run", flush=True)
    if workers == 1:
        for j in jobs:
            _ml_job(j)
        return
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_ml_job, jobs):
            pass


# ---------------------------------------------------------------- port check: the Python reference vs the C++ core
def _null(a, b):
    """dB, a against g b (one gain per channel), mean of L and R; windows of the caller's choosing."""
    v = []
    for ch in (0, 1):
        x, y = a[:, ch], b[:, ch]
        g = np.dot(x, y) / np.dot(y, y)
        v.append(10 * np.log10(np.sum((x - g * y) ** 2) / np.sum(x ** 2) + 1e-30))
    return float(np.mean(v))


def py_reverb(M, key, n, theta):
    import ml_carryover as cc
    import ml_variance as mv
    x = np.zeros(n)
    x[0] = 1.0
    for k, Mk in enumerate(cc.DIFF[M.prog]):
        x = cc.ap(x, Mk, cc.G_AP, M.b.get(key, 0.0) if k == len(cc.DIFF[M.prog]) - 1 else 0.0)
    y = mv.core(M, key, x, theta)
    return np.stack([np.convolve(y[:, ch], M.h)[:n] for ch in (0, 1)], 1)


def port_stage(only):
    import ml_lines as mll
    import ml_variance as mv
    FSD = 19.2e6 / 728
    A = 0.13 * np.sqrt(2)
    w = 2 * np.pi * 0.83 / FSD
    res = {}
    for prog, key in groups(only):
        t = time.time()
        M = mll.Model(prog)
        dev = {v: sf.read(dpath(prog, key, v), always_2d=True, dtype="float64")[0] for v in VARS}
        n = len(dev["static"])
        ref0 = py_reverb(M, key, n, np.zeros(n))
        # integer lag of the plugin's core output against the reference (static)
        c = np.fft.irfft(np.fft.rfft(dev["static"][:, 1], 2 * n) * np.conj(np.fft.rfft(ref0[:, 1], 2 * n)), 2 * n)
        lag = int(np.argmax(np.abs(c[:4000])))
        sh = lambda y: np.vstack([np.zeros((lag, 2)), y])[:n]
        W = [(0.0, 0.15), (0.15, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, DEV_LEN - 0.2)]
        wins = lambda a, b: {f"{t0}-{t1}": _null(a[int(t0 * FSD):int(t1 * FSD)], b[int(t0 * FSD):int(t1 * FSD)])
                             for t0, t1 in W if int(t1 * FSD) <= n and np.any(a[int(t0 * FSD):int(t1 * FSD)])}
        r = dict(lag=lag, static=wins(dev["static"], sh(ref0)))
        for i, v in enumerate(VARS[:3]):
            n0 = int(round(PRE[i] * FSD))
            d = dev[v]                              # the file starts at the input sample (pre is cut), theta n0 on
            best = None
            for dd in range(-3, 4):
                th = A * np.sin(w * (np.arange(n) + n0 + 1 + dd + lag))
                ref = py_reverb(M, key, n, th)
                sc = _null(d[int(0.15 * FSD):int(1.0 * FSD)], sh(ref)[int(0.15 * FSD):int(1.0 * FSD)])
                if best is None or sc < best[0]:
                    best = (sc, dd, ref)
            r[v] = dict(d=best[1], null=wins(d, sh(best[2])), vs_static=wins(d, dev["static"]))
        res[f"{prog}_{key}"] = r
        print(f"{prog} {key}: {time.time() - t:.0f} s lag {lag} static " + " ".join(f"{x:.0f}" for x in r["static"].values())
              + " | p1 d=" + str(r["p1"]["d"]) + " " + " ".join(f"{x:.0f}" for x in r["p1"]["null"].values()), flush=True)
        os.makedirs(OUT, exist_ok=True)
        old = json.load(open(os.path.join(OUT, "port.json"))) if os.path.exists(os.path.join(OUT, "port.json")) else {}
        old.update({f"{prog}_{key}": r})
        json.dump(old, open(os.path.join(OUT, "port.json"), "w"), indent=1)


# ---------------------------------------------------------------- sweep vs impulse measurement (the plugin's sim renders)
def sweep_stage(only):
    import ml_carryover as cc
    import ml_lines as mll
    import ml_structure as ml
    import ml_variance as mv
    from scipy.signal import fftconvolve
    s, inv = pk.sweep_96k()
    ref = fftconvolve(s, inv)
    i0, rmax = int(np.argmax(np.abs(ref))), np.abs(ref).max()
    del ref
    out = json.load(open(os.path.join(OUT, "sweep.json"))) if os.path.exists(os.path.join(OUT, "sweep.json")) else {}
    for prog, key in groups(only):
        if key not in SIM_KEYS[prog] or f"{prog}_{key}" in out:
            continue
        t = time.time()
        g = cc.groups(prog)[key]
        n = max(len(cc._IR(f)) for f in g.values()) + 200
        caps, _, amps = cc.captures(prog, key, n)
        n = len(caps[amps[-1]])
        i1 = min(min(mll.usable_end(caps[a], g[a], prog)) for a in amps)
        anchor = pk.onset_anchor(pk.load_render(hpath(prog, key, max(amps), "static")))
        _, _, shv = parts(prog, key, amps, n, "static")
        irs = {}
        for v in VARS:
            _, R, _ = parts(prog, key, amps, n, v)
            irs[v] = R
        for i in range(len(pk.SIM_PRE) + 1):
            x = pk.load_render(spath(prog, key, i)) / 0.25
            h = np.stack([fftconvolve(x[:, c], inv)[i0:i0 + int(tail(key) * pk.HOST)] for c in (0, 1)], 1) / rmax
            del x
            d = sosfiltfilt(ml.BAND, pk.to_device(h, anchor), axis=0)
            d = np.vstack([d, np.zeros((max(0, n - len(d)), 2))])[:n]
            irs[f"sw{i}" if i < len(pk.SIM_PRE) else "swst"] = cc.fshift(d, shv)
        i1 = min(i1, cc.K0 + int((tail(key) - 0.5) * FS))
        t0 = (cc.ONSET[prog][1] + 20) / FS
        fw = (t0, t0 + 0.25)
        Pst = pc.prep(irs["static"][:i1], dtype=np.float32)
        r = dict(i1=(i1 - cc.K0) / FS, vs_static={}, null_static={})
        for v in [x for x in irs if x != "static"]:
            Pv = pc.prep(irs[v][:i1], dtype=np.float32)
            r["vs_static"][v] = slim(pc.compare(Pv, Pst, i1, cc.K0, t0, fw))
            r["null_static"][v] = pc._ml_nulls(prog, irs[v], irs["static"], i1)
            del Pv
        del Pst
        # the sweep-measured plugin against each capture (decay metrics only: they are gain-free and start after the
        # cluster, so the AMPL-0 sweep IR stands for every AMPL), beside the impulse renders' (p0-p2) in the same window
        r["vs_cap"] = {}
        for amp in amps:
            ia = min(min(mll.usable_end(caps[amp], g[amp], prog)), cc.K0 + int((tail(key) - 0.5) * FS))
            Y = pc.prep(caps[amp][:ia], dtype=np.float32)
            r["vs_cap"][str(amp)] = {}
            for v in [x for x in irs if x.startswith("sw") or x in VARS]:
                Pv = pc.prep(irs[v][:ia], dtype=np.float32)
                r["vs_cap"][str(amp)][v] = {k: x for k, x in slim(pc.compare(Pv, Y, ia, cc.K0, t0, fw)).items() if k.startswith(("oct_", "third_", "T"))}
                del Pv
            del Y
        out[f"{prog}_{key}"] = r
        json.dump(out, open(os.path.join(OUT, "sweep.json"), "w"), indent=1)
        f = 125
        print(f"{prog} {key}: {time.time() - t:.0f} s; {f} Hz decay / static - 1: imp " + " ".join(f"{100 * r['vs_static'][v][f'oct_t{f}']:+.1f}" for v in VARS[:3])
              + " sweep " + " ".join(f"{100 * r['vs_static'][f'sw{i}'][f'oct_t{f}']:+.1f}" for i in range(len(pk.SIM_PRE))), flush=True)


# ---------------------------------------------------------------- report
def report_stage():
    print("report: see analysis/ml_fix/diag.log (written by ml_fix_report.py)")


# ---------------------------------------------------------------- item 1 smoke: the audit through the sweep, 3 phases
SW3_KEYS = {"M": KEYS["M"] + ["0.7_1.8_0.2", "0.7_1.2_0.4", "1.2_0.6_0.9", "1.2_2.6_0.7"],
            "L": KEYS["L"] + ["0.8_4.5_1.0", "0.8_6.0_0.8", "1.2_2.2_0.8", "0.8_15.0_0.5"]}


def tail3(key):
    """s of sweep tail: 6 s, or 0.8 RT (L RT 20: 16 s), so the late decay (-25/-40 dB) fits."""
    return max(6.0, min(18.0, 0.8 * float(key.split("_")[1])))


def sw3path(prog, key, i):
    return os.path.join(RND, "96000", f"{prog}_{key}_a0_sw3_{i}.wav")


SW3_DIR = os.environ.get("SW3_DIR", "sw3")         # sw3: the 20-group smoke (31 Oct, pre-AMPL-fix build); sw3_all: all groups


def sw3_groups(only, every=None):
    if every == "all":
        import ml_carryover as cc
        g = [(p, k) for p in ("M", "L") for k in sorted(cc.groups(p))]
    else:
        g = [(p, k) for p, ks in SW3_KEYS.items() for k in ks]
    return [x for x in g if not only or f"{x[0]}_{x[1]}" in only]


def sw3_lines(prog, key):
    sw = os.path.join(RND, "sweep_96k.wav")
    T0 = len(pk.sweep_96k()[0]) / pk.HOST
    g = pk.ml_groups(prog)[key]
    lf, rt, hf = (float(x) for x in key.split("_"))
    L = [pk.render_line(hpath(prog, key, amp, "static"), pk.HOST, prog, lf, rt, hf, amp, g["len"], mod=0)
         for amp in sorted(set(g["amps"]) | {0})]
    L += [pk.render_line(sw3path(prog, key, i), pk.HOST, prog, lf, rt, hf, 0, T0 + tail3(key), pre=pre, inp=sw) for i, pre in enumerate(PRE)]
    return L


def sweep_file():
    sw = os.path.join(RND, "sweep_96k.wav")
    if not os.path.exists(sw):
        os.makedirs(RND, exist_ok=True)
        sf.write(sw, (0.25 * pk.sweep_96k()[0]).astype(np.float32), pk.HOST, subtype="FLOAT")


def render_group(prog, key):
    """Render what this group lacks (in the job: one list per group)."""
    todo = [l for l in sw3_lines(prog, key) if not os.path.exists(l.split()[0][4:])]
    if not todo:
        return
    os.makedirs(os.path.join(RND, "96000"), exist_ok=True)
    lst = os.path.join(RND, f"list_{prog}_{key}.txt")
    open(lst, "w").write("\n".join(todo) + "\n")
    r = subprocess.run([pk.RENDER_BIN, lst], capture_output=True, text=True)
    os.remove(lst)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])


def sw3_render(only):
    sw = os.path.join(RND, "sweep_96k.wav")
    if not os.path.exists(sw):
        os.makedirs(RND, exist_ok=True)
        sf.write(sw, (0.25 * pk.sweep_96k()[0]).astype(np.float32), pk.HOST, subtype="FLOAT")
    T0 = len(pk.sweep_96k()[0]) / pk.HOST
    L = []
    for prog, key in sw3_groups(only):
        g = pk.ml_groups(prog)[key]
        lf, rt, hf = (float(x) for x in key.split("_"))
        for amp in sorted(set(g["amps"]) | {0}):         # the static impulse: the cluster part, anchor and shift
            L.append(pk.render_line(hpath(prog, key, amp, "static"), pk.HOST, prog, lf, rt, hf, amp, g["len"], mod=0))
        for i, pre in enumerate(PRE):
            L.append(pk.render_line(sw3path(prog, key, i), pk.HOST, prog, lf, rt, hf, 0, T0 + tail3(key), pre=pre, inp=sw))
    todo = [l for l in L if not os.path.exists(l.split()[0][4:])]
    print(f"sw3 render: {len(todo)} of {len(L)} to do", flush=True)
    if todo:
        lst = os.path.join(OUT, "render_sw3.txt")
        os.makedirs(os.path.join(RND, "96000"), exist_ok=True)
        open(lst, "w").write("\n".join(todo) + "\n")
        t = time.time()
        r = subprocess.run([pk.RENDER_BIN, lst], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(r.stderr[-2000:])
        print(f"sw3 render: done in {time.time() - t:.0f} s", flush=True)


def _sw3_job(a):
    """Per capture: the sweep-measured plugin (static cluster part + the sweep's reverb part at phase k, fit_warp exactly as
    the audit) vs the capture, the static impulse vs the capture, and the sweep phases against each other."""
    prog, key, path = a[:3]
    inline = len(a) > 3 and a[3]
    import ml_carryover as cc
    import ml_lines as mll
    import ml_structure as ml
    from scipy.signal import fftconvolve
    t = time.time()
    if inline:
        render_group(prog, key)
    tr = time.time() - t
    K0 = cc.K0
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, _, amps = cc.captures(prog, key, n)
    n = len(caps[amps[-1]])
    C, Rst, shv = parts(prog, key, amps, n, "static")
    anchor = pk.onset_anchor(pk.load_render(hpath(prog, key, max(amps), "static")))
    s, inv = pk.sweep_96k()
    ref = fftconvolve(s, inv)
    i0, rmax = int(np.argmax(np.abs(ref))), np.abs(ref).max()
    del ref
    R = {}
    for i in range(len(PRE)):
        x = pk.load_render(sw3path(prog, key, i)) / 0.25
        h = np.stack([fftconvolve(x[:, c], inv)[i0:i0 + int(tail3(key) * pk.HOST)] for c in (0, 1)], 1) / rmax
        del x
        d = sosfiltfilt(ml.BAND, pk.to_device(h, anchor), axis=0)
        del h
        R[f"s{i}"] = cc.fshift(np.vstack([d, np.zeros((max(0, n - len(d)), 2))])[:n], shv)
    R["static"] = Rst
    t0 = (cc.ONSET[prog][1] + 20) / FS
    fw = (t0, t0 + 0.25)
    names = list(R)
    rows = []
    for amp in amps:
        y = caps[amp]
        i1 = min(min(mll.usable_end(y, g[amp], prog)), K0 + int((tail3(key) - 0.5) * FS))
        sl, nf = mll.fit_window(prog, amp > 0)
        Y = pc.prep(y[:i1], dtype=np.float32)
        P, mdl = {}, {}
        row = dict(amp=amp, t_end=(i1 - K0) / FS, vs_cap={}, null_cap={}, pairs={}, null_pairs={})
        for v in names:
            mdl[v], u, gs, _ = mll.fit_warp(y, [C[amp], R[v]] if amp > 0 else [R[v]], sl, nf)
            row.setdefault("fit", {})[v] = dict(shift=u[0], ppm=u[1], gains=[[float(x) for x in gch] for gch in gs])
            P[v] = pc.prep(mdl[v][:i1], dtype=np.float32)
            row["vs_cap"][v] = slim(pc.compare(P[v], Y, i1, K0, t0, fw))
            row["null_cap"][v] = pc._ml_nulls(prog, y, mdl[v], i1)
        sw = [v for v in names if v != "static"]
        for i, a_ in enumerate(sw):
            for b_ in sw[i + 1:]:
                row["pairs"][f"{a_}-{b_}"] = slim(pc.compare(P[a_], P[b_], i1, K0, t0, fw))
                row["null_pairs"][f"{a_}-{b_}"] = pc._ml_nulls(prog, mdl[a_], mdl[b_], i1)
        rows.append(row)
        del P, Y, mdl
    lf, rt, hf = (float(x) for x in key.split("_"))
    json.dump(dict(key=key, prog=prog, lf=lf, rt=rt, hf=hf, render_s=tr, rows=rows), open(path, "w"))
    if inline:                                      # disk: the sweep renders go once scored (the static impulses stay)
        for i in range(len(PRE)):
            os.remove(sw3path(prog, key, i))
    r = rows[-1]
    f = 125 if prog == "M" else 250
    print(f"{prog} {key}: {time.time() - t:.0f} s (render {tr:.0f}); {f} Hz decay vs cap " + " ".join(f"{v} {100 * (r['vs_cap'][v][f'oct_t{f}'] or np.nan):+.1f}" for v in names), flush=True)
    return key


def sw3_stage(workers, only, every=None, chunk=None, rtmin=None, rtmax=None):
    """every='all': every M / L group, rendering inside the job (render, score, delete the sweeps). rtmin / rtmax: L's
    RT range (the long L groups run on one worker), M always below rtmin's cut unless rtmin is set."""
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, SW3_DIR)
    os.makedirs(d, exist_ok=True)
    sweep_file()
    G = sw3_groups(only, every)
    rt = lambda k: float(k.split("_")[1])
    if rtmin is not None:
        G = [(p, k) for p, k in G if p == "L" and rt(k) >= rtmin]
    if rtmax is not None:
        G = [(p, k) for p, k in G if not (p == "L" and rt(k) >= rtmax)]
    if chunk:
        i, m = (int(x) for x in chunk.split("/"))
        G = G[i::m]
    jobs = [(p, k, os.path.join(d, f"{p}_{k}.json"), every == "all") for p, k in G]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    jobs.sort(key=lambda j: -rt(j[1]) * (3 if j[0] == "L" else 1))
    print(f"sw3: {len(jobs)} groups to run", flush=True)
    if workers == 1:
        for j in jobs:
            _sw3_job(j)
        return
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_sw3_job, jobs):
            pass


# ---------------------------------------------------------------- static windows on the (new build's) static impulse renders
def _static_job(a):
    """cluster + gap and first pass, per channel, the static impulse render fitted to each capture exactly as plugin_check's
    ml stage (fit_warp on cluster + gap + first pass): the deterministic part, where the -50 dB bar applies."""
    prog, key = a
    import ml_carryover as cc
    import ml_lines as mll
    K0 = cc.K0
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, _, amps = cc.captures(prog, key, n)
    n = len(caps[amps[-1]])
    C, R, _ = parts(prog, key, amps, n, "static")
    W = mll.windows(prog)
    rows = []
    for amp in amps:
        y = caps[amp]
        end = mll.usable_end(y, g[amp], prog)
        sl, nf = mll.fit_window(prog, amp > 0)
        mdl, u, gs, _ = mll.fit_warp(y, [C[amp], R] if amp > 0 else [R], sl, nf)
        rows.append(dict(amp=amp, gains=[[float(x) for x in gch] for gch in gs],
                         null={w: [mll.null_db(y[:, ch], mdl[:, ch], K0 + W[ch][w][0], min(K0 + W[ch][w][1], end[ch])) for ch in (0, 1)]
                               for w in ("cluster_gap", "first_pass") if amp > 0 or w != "cluster_gap"}))
    return f"{prog}_{key}", rows


def static_stage(workers, only):
    from concurrent.futures import ProcessPoolExecutor
    import ml_carryover as cc
    p = os.path.join(OUT, "static_windows.json")
    res = json.load(open(p)) if os.path.exists(p) else {}
    jobs = [(pr, k) for pr in ("M", "L") for k in sorted(cc.groups(pr)) if f"{pr}_{k}" not in res and (not only or f"{pr}_{k}" in only)
            and os.path.exists(hpath(pr, k, 0, "static"))]
    print(f"static: {len(jobs)} groups", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=4) as ex:
        for i, (k, rows) in enumerate(ex.map(_static_job, jobs)):
            res[k] = rows
            if i % 20 == 0:
                json.dump(res, open(p, "w"))
    json.dump(res, open(p, "w"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["render", "ml", "port", "sweep", "sw3render", "sw3", "static", "report"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--set", default=None, help="sw3: 'all' = every M / L group, rendered inside the job")
    ap.add_argument("--chunk", default=None)
    ap.add_argument("--rtmin", type=float, default=None)
    ap.add_argument("--rtmax", type=float, default=None)
    a = ap.parse_args()
    only = a.only.split(",") if a.only else None
    {"render": lambda: render_stage(only), "ml": lambda: ml_stage(a.workers, only), "port": lambda: port_stage(only),
     "sweep": lambda: sweep_stage(only), "sw3render": lambda: sw3_render(only),
     "static": lambda: static_stage(a.workers, only), "sw3": lambda: sw3_stage(a.workers, only, a.set, a.chunk, a.rtmin, a.rtmax), "report": report_stage}[a.stage]()


if __name__ == "__main__":
    main()
