"""Phase 5 of the plugin engine: the plugin's own output (ReverbEngine: device-rate cores + DeviceBridge with the
analogue chain, at the host rate) against the Tomkin 96 kHz IRs, scored with the research benches' own rules
(CLAUDE.md, docs/NEXT_STEPS.md §0 / §2, plugin/docs/ENGINE.md). Results: analysis/plugin_check/REPORT.md.

    .venv/bin/python tools/memguard.py --limit-gb 6 -- .venv/bin/python tools/plugin_check.py --stage <stage> [...]

The renders come from plugin/build-render's REV246Render (plugin/tools/render/RenderIr.cpp): a unit impulse into the
whole engine at the host rate, DELAY 0, LOW PASS 12 kHz (bypass), FLOAT (16-bit off). Large files go to
analysis/plugin_check/renders/ (gitignored); one small JSON per setting (resumable).

How the plugin enters the benches. A capture reaches the benches as the 96 kHz file deconvolved with the chain3 kernel
and resampled to the device grid (ml_structure._make), each channel's cluster tap 0 then aligned (q1 data stage, S;
ml_carryover, M / L). The plugin's 96 kHz render goes through the same deconvolution and resampling (one fixed anchor:
the engine's timing is the same for every setting) and the same per-channel alignment, so plugin and capture differ
only by what differs between the plugin and the unit; the bench's per-capture fit (clock, shift for AMPL 0, gain) is
then run with the plugin in the model's place.

  render   write the render lists and run REV246Render (--set s|ml|n|e|rates|sim, or all)
  align    the fixed anchor and per-channel tap-0 shifts (S 1.0_1.0_1.0 AMPL 90 render) -> align.json
  s        Program S, every clean capture of the 128 settings: (a) the plugin against the bench's own model (full_model,
           with and without the sub-100 Hz capture transfer), (b) the plugin against the capture: pulse nulls
           (100 Hz-10 kHz, W_COMP windows) and sml_fullband's bands x windows with the cross-take repeat beside them;
           the bench model's nulls for the same captures and fits beside the plugin's (the chain verdict); the plugin
           with the netlists' coupling corners (0.617 / 0.157 Hz) in place of the bridge's (0.72 / 0.16 Hz).
  loff     L's output offset: per capture, a common clock and one shift per channel, plugin render against the 96 kHz
           capture (cluster + first reverb), in the capture domain: the library's L - R offset = 0.5 + dL - dR.
  report   medians -> analysis/plugin_check/summary_<stage>.json and stdout.
"""
import argparse, glob, json, os, subprocess, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("EMT_KERNEL", "chain3")       # the benches' device-domain kernel (NEXT_STEPS §2)
os.environ.setdefault("EMT_PMAX", "0.8475")         # the S baseline (pulse_check, s_spec)
import numpy as np
import soundfile as sf
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "analysis", "plugin_check")
RENDERS = os.path.join(OUT, "renders")
RENDER_BIN = os.path.join(ROOT, "plugin", "build-render", "REV246Render_artefacts", "Release", "REV246Render")
LIB = os.path.join(ROOT, "irs", "EMT 246 captures", "96kHz")
HOST = 96000
LR_DEV = 0.5                                        # DeviceBridge::lLateDeviceSamples
NETLIST_HP = dict(inp=(0.720158113537988, 0.61705364), out=(0.15915494309189537, 0.15671771))   # bridge -> netlist (Hz)


# ---------------------------------------------------------------- renders
def rpath(rate, prog, key, amp, tag=""):
    return os.path.join(RENDERS, str(rate), f"{prog}_{key}_a{amp:g}{tag}.wav")


def render_line(out, rate, prog, lf, rt, hf, amp, length, mod=1, pre=0.0, inp="impulse", domain="host"):
    return (f"out={out} rate={rate} program={prog} lf={lf} rt={rt} hf={hf} ampl={amp} delay=0 lp=12 sixteen=0 "
            f"mod={mod} len={length:.3f} pre={pre:.4f} input={inp} domain={domain}")


def run_renders(lines, name):
    todo = [l for l in lines if not os.path.exists(l.split()[0][4:])]
    print(f"{name}: {len(todo)} of {len(lines)} renders to do", flush=True)
    if not todo:
        return
    for l in todo:
        os.makedirs(os.path.dirname(l.split()[0][4:]), exist_ok=True)
    lst = os.path.join(OUT, f"render_{name}.txt")
    open(lst, "w").write("\n".join(todo) + "\n")
    r = subprocess.run([RENDER_BIN, lst], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    lat = sorted({l.split("latency=")[1].split()[0] for l in r.stdout.splitlines() if "latency=" in l})
    print(f"{name}: done, latency {lat}", flush=True)


def s_settings():
    info = json.load(open(os.path.join(ROOT, "analysis", "q1", "data_chain3.json")))
    out = {}
    for k, g in info.items():
        if k.startswith("_"):
            continue
        n = max(f["n"] for a, f in g["files"].items() if not a.startswith("_"))
        out[k] = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"], len=n / 26373.6 + 0.4)
    return out


def render_stage(which):
    lines = {}
    if which in ("s", "all"):
        L = []
        for k, g in s_settings().items():
            for amp in (0, 30, 60, 90):
                L.append(render_line(rpath(HOST, "S", k, amp), HOST, "S", g["lf"], g["rt"], g["hf"], amp, g["len"]))
        lines["s"] = L
    if which in ("ml", "all"):
        lines["ml"] = ml_render_lines()
    if which in ("sim", "all"):
        lines["sim"] = sim_render_lines()
    if which in ("n", "all"):
        lines["n"] = n_render_lines()
    if which in ("e", "all"):
        lines["e"] = e_render_lines()
    if which in ("rates", "all"):
        lines["rates"] = rates_render_lines()
    for name, L in lines.items():
        run_renders(L, name)


# ---------------------------------------------------------------- the device-domain path (as ml_structure._make)
_H = {}


def _harness():
    if not _H:
        sys.path.insert(0, HERE)
        sys.path.insert(0, os.path.join(HERE, "harness"))
        import ml_structure as ms
        from core import PAD, Kernel, deconvolve, to_device_grid
        _H.update(ms=ms, PAD=PAD, Kernel=Kernel, deconvolve=deconvolve, to_device_grid=to_device_grid,
                  K=Kernel.load(ms.KERNEL))
    return _H


def onset_anchor(x):
    """ml_structure._make's onset (capture samples) of a 96 kHz stereo IR."""
    a = np.abs(x).max(1)
    main = int(np.argmax(a > a.max() * 0.1))
    lo = max(main - int(0.002 * HOST), 0)
    return lo + int(np.argmax(a[lo:main + 1] > a.max() * 1e-2))


def to_device(x, anchor, fs=HOST):
    """A 96 kHz render (n, 2) on the benches' device grid exactly as a capture (_make), with a fixed anchor: device
    index K0 = the anchor; L read 19.19 us later (the cache's L / R offset)."""
    h = _harness()
    ms = h["ms"]
    d = h["deconvolve"](np.vstack([np.zeros((h["PAD"], 2)), x]), fs, h["K"])
    a = float(anchor + h["PAD"])
    k0 = int(np.floor(a / fs * ms.FSD))
    y = np.stack([h["to_device_grid"](d[:, 0], fs, a + ms.LR_US * 1e-6 * fs, k0, ms.CFG),
                  h["to_device_grid"](d[:, 1], fs, a, k0, ms.CFG)], 1)
    return y[max(k0 - ms.K0, 0):]


def load_render(path):
    x, fs = sf.read(path, always_2d=True, dtype="float64")
    assert fs == HOST, (path, fs)
    return x


def align_stage():
    """Fixed anchor (capture samples) and the per-channel shifts that put cluster tap 0 on q1's PRE."""
    sys.path.insert(0, HERE)
    import q1_generalise as q1
    x = load_render(rpath(HOST, "S", "1.0_1.0_1.0", 90))
    anchor = onset_anchor(x)
    h = to_device(x, anchor)
    t0R, t0L = q1._tap0(h)
    # the same tap in the 96 kHz render, for the record (R at 96 kHz, L - R in device samples as the harness measures)
    from core import subsample_peak
    pk = [subsample_peak(x[:, c], int(np.argmax(np.abs(x[:300 + anchor, c])))) for c in (0, 1)]
    out = dict(anchor=anchor, t0=[t0L, t0R], shift=[q1.PRE - t0L, q1.PRE - t0R],
               lr_96k_dev=float((pk[0][0] - pk[1][0]) / HOST * 26373.626), lr_device_after_cache=float(t0L - t0R))
    json.dump(out, open(os.path.join(OUT, "align.json"), "w"), indent=1)
    print(out)


def plugin_q1(path, n, al):
    """A render in q1's frame: device grid, each channel's cluster tap 0 at PRE (al: align.json), n samples."""
    sys.path.insert(0, HERE)
    import q1_generalise as q1
    h = to_device(load_render(path), al["anchor"])
    y = np.stack([q1.fshift(h[:, c], al["shift"][c]) for c in (0, 1)], 1)
    if len(y) < n:
        y = np.vstack([y, np.zeros((n - len(y), 2))])
    return sosfiltfilt(q1.LFBAND, y[:n], axis=0)       # the band every S capture and model has (q2._band, Joint.render)


# ---------------------------------------------------------------- the coupling corners (H item 1)
def corner_ratio(f):
    """Netlist corners over the bridge's (input 0.617 for 0.720 Hz, output 0.157 for 0.159 Hz): an LTI ratio applied to
    a bridge render gives the render with the netlists' corners (first-order high-passes s / (s + w))."""
    s = 2j * np.pi * np.asarray(f, float)
    w = lambda hz: 2 * np.pi * hz
    r = np.ones_like(s)
    for a, b in NETLIST_HP.values():
        r = r * (s + w(a)) / (s + w(b))
    return r


def apply_ratio(x, fs, fun):
    N = 1 << int(np.ceil(np.log2(len(x) + int(4 * fs))))       # 4 s of room for the 0.157 Hz pole's tail
    H = fun(np.fft.rfftfreq(N, 1 / fs))
    X = np.fft.rfft(x, N, axis=0) * (H if H.ndim == 2 else H[:, None])
    return np.fft.irfft(X, N, axis=0)[:len(x)]


# ---------------------------------------------------------------- the analogue chain: bridge (netlists) vs chain3
CHAIN = os.path.join(OUT, "chain_ratio.npz")
_CR = {}


def bridge_chain(f):
    """DeviceBridge's modelled chain (input filter x output filter x S&H lag, coupling corners), unity at 1 kHz (R)."""
    sys.path.insert(0, ROOT)
    from tools.plugin_codegen import bridge as b
    _, _, _, oi = b.side(b.INPUT_FILTER, b.IN_HP_HZ)
    _, _, _, oo = b.side(b.OUTPUT_FILTER, b.OUT_HP_HZ)
    return b.mna_side(b.INPUT_FILTER, oi, b.IN_HP_HZ, f) * b.mna_side(b.OUTPUT_FILTER, oo, b.OUT_HP_HZ, f) * b.lag(f)


def chain_stage():
    """chain3 (the benches' measured per-channel chain) over the bridge's chain on a 2 Hz grid to 13.2 kHz, per channel
    (L with the bridge's 0.5 device samples, chain3's L relative to the cache's 19.19 us read), and the per-band
    mismatch after the best gain and delay (100 Hz-10 kHz, equal weight per octave) -> chain_ratio.npz, chain.json."""
    from scipy.optimize import minimize
    h = _harness()
    f = np.arange(0, 13200, 2.0)
    f[0] = 0.5
    N = 1 << 19
    fg = np.fft.rfftfreq(N, 1 / HOST)
    K3g = h["K"].response(N, HOST)
    K3 = np.stack([np.interp(f, fg, K3g[:, c].real) + 1j * np.interp(f, fg, K3g[:, c].imag) for c in (0, 1)], 1)
    Hb = bridge_chain(f)
    T = 728 / 19.2e6
    w = 2 * np.pi * f
    ratio = np.stack([K3[:, 0] / (Hb * np.exp(-1j * w * LR_DEV * T)), K3[:, 1] / Hb], 1)   # chain3 / bridge
    np.savez(CHAIN, f=f, ratio=ratio)
    edges = [20, 40, 80, 160, 270, 400, 800, 1600, 3200, 6400, 10000, 11000, 12300]
    res = {}
    for c, name in ((0, "L"), (1, "R")):
        r = 1 / ratio[:, c]                          # bridge / chain3
        fit = (f >= 100) & (f <= 10000)
        cost = lambda p: np.mean(np.abs(p[0] * r[fit] * np.exp(1j * w[fit] * p[1] * T) - 1) ** 2 / f[fit])
        x = minimize(cost, [1 / np.median(np.abs(r[fit])), 0.0], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-15)).x
        rr = x[0] * r * np.exp(1j * w * x[1] * T)
        rows = {}
        for a_, b_ in zip(edges[:-1], edges[1:]):
            m = (f >= a_) & (f < b_)
            rows[f"{a_}-{b_}"] = dict(db=[round(float(v), 3) for v in (20 * np.log10(np.abs(rr[m]))).take([0, -1])],
                                      db_range=[round(float(20 * np.log10(np.abs(rr[m])).min()), 3), round(float(20 * np.log10(np.abs(rr[m])).max()), 3)],
                                      deg_range=[round(float(np.degrees(np.angle(rr[m])).min()), 2), round(float(np.degrees(np.angle(rr[m])).max()), 2)],
                                      null_equiv_db=round(float(10 * np.log10(np.mean(np.abs(rr[m] - 1) ** 2))), 1))
        res[name] = dict(gain_db=round(float(20 * np.log10(x[0])), 3), delay_device_samples=round(float(x[1]), 4), bands=rows)
        print(name, res[name]["delay_device_samples"], {k: v["null_equiv_db"] for k, v in rows.items()})
    json.dump(res, open(os.path.join(OUT, "chain.json"), "w"), indent=1)
    # the correction keeps each channel's tap-0 alignment: its best delay (above) taken out
    x = np.array([res["L"]["delay_device_samples"], res["R"]["delay_device_samples"]])
    np.savez(CHAIN, f=f, ratio=ratio * np.exp(-1j * w[:, None] * x[None, :] * T))


def chain3_over_bridge(f):
    """(n, 2): chain3 / the bridge's chain, interpolated from chain_ratio.npz (1 above 13.2 kHz)."""
    if not _CR:
        z = np.load(CHAIN)
        _CR.update(f=z["f"], r=z["ratio"])
    out = np.ones((len(f), 2), complex)
    m = f <= _CR["f"][-1]
    for c in (0, 1):
        out[m, c] = np.interp(f[m], _CR["f"], _CR["r"][:, c].real) + 1j * np.interp(f[m], _CR["f"], _CR["r"][:, c].imag)
    return out


# ---------------------------------------------------------------- program S
def refit_shift(mdl, c, ends, n):
    """sml_fullband.refit_capture with the shift free for every capture (q1._clock_job's AMPL-0 form): y, g, te, fit."""
    import pulse_check as pc
    q2, q1 = pc.q2, pc.q1
    yb = q2._band(c["file"])[:n]
    if c["fade2"] >= 0.15:
        _, (e, sh, gain) = q1._clock_job((yb, mdl, c["fade2"], c["eps"], True))
    else:
        g0 = q2.ramp_gain(n, ends[c["file"]])
        _, (e, sh, gain) = q1._clock_job((yb, mdl * g0, q2.comp_end(yb, g0), c["eps"], True))
    y = q1.apply_clock(yb, e, sh)[:n] / gain
    g = q2.ramp_gain(n, [q1.PRE + (x + sh - q1.PRE) / (1 + e * 1e-6) for x in ends[c["file"]]])
    return y, g, q2.comp_end(y, g), dict(eps=e, shift=sh, gain=gain)


def _s_job(a):
    key, lab, b230, law, fold, fc, v, caps, ends, path, al = a
    import pulse_check as pc
    import sml_fullband as sml
    q2, q1 = pc.q2, pc.q1
    q1.EXT = 1700                                    # fit_capture's excitation layout, as sml_fullband
    G = pc.setting_G(lab, b230, law, fold, fc)
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
    FS, PRE = q1.FS, q1.PRE
    P, out = {}, dict(key=key, **lab, caps={}, pairs={}, pairs_bench={})
    for c in caps:
        # the bench: its model, its fit (clock, gain; shift for AMPL 0), its full-band model with / without the
        # sub-100 Hz capture transfer
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        n = len(y)
        sml.SML_LOWBAND = False
        mf0 = sml.full_model(J, v, G, n, c["amp"], mdl)
        sml.SML_LOWBAND = True
        mf1 = sml.full_model(J, v, G, n, c["amp"], mdl)
        # the plugin through the capture's own path, fitted as the bench fits its model
        p = plugin_q1(rpath(HOST, "S", key, c["amp"]), n, al)
        yp, gp, tep = sml.refit_capture(p, c, ends, n)
        ys, gs, tes, fit_s = refit_shift(p, c, ends, n)   # the same with the shift free for every capture (3 numbers)
        pn = apply_ratio(p, FS, corner_ratio)        # H item 1: the netlists' coupling corners
        p3 = apply_ratio(p, FS, chain3_over_bridge)  # the plugin with chain3 in place of the bridge's chain
        y3, g3, te3, fit_3 = refit_shift(p3, c, ends, n)
        # (a) plugin against the bench model: one gain (both channels) over 0-300 ms, 100 Hz-10 kHz
        i1 = min(PRE + int(0.3 * FS), n)
        pb, mb = sosfiltfilt(q1.ms.BAND, p[:i1], axis=0)[PRE:], sosfiltfilt(q1.ms.BAND, mf0[:i1], axis=0)[PRE:]
        k = float(np.sum(pb * mb) / np.sum(pb * pb))
        te_ = min(te, tep, tes, te3)
        sig = dict(y=yp, e=yp - p * gp, en=yp - pn * gp, ys=ys, es=ys - p * gs, y3=y3, e3=y3 - p3 * g3,
                   yb=y, eb0=y - mf0 * g, eb1=y - mf1 * g, m0=mf0, d0=k * p - mf0, d1=k * p - mf1)
        out["caps"][c["amp"]] = dict(
            file=c["file"], te=te_, te_plugin=tep, te_bench=te, k_plugin_over_bench=k, fit_shift=fit_s, fit_c3=fit_3,
            null=q1.nulls(yp, p * gp, tep, q2.W_COMP, narrow=True),
            null_shift=q1.nulls(ys, p * gs, tes, q2.W_COMP, narrow=True),
            null_netlist_hp=q1.nulls(yp, pn * gp, tep, q2.W_COMP, narrow=True),
            null_c3=q1.nulls(y3, p3 * g3, te3, q2.W_COMP, narrow=True),
            null_bench=q1.nulls(y, mdl * g, te, q2.W_COMP, narrow=True),
            null_oct=q2.octave_nulls(ys, p * gs, tes), bands=sml.band_energies(sig, te_))
        P[c["amp"]] = dict(y=yp, g=gp, te=tep, yb=y, gb=g, teb=te)
    for a_, b_ in pc.CROSS:                          # the cross-take repeat: the bench's own fits, and the plugin's
        if a_ not in P or b_ not in P:
            continue
        A, B = P[a_], P[b_]
        for kk, (yk, gk, tk) in (("pairs_bench", ("yb", "gb", "teb")), ("pairs", ("y", "g", "te"))):
            n = min(len(A[yk]), len(B[yk]))
            te = min(A[tk], B[tk])
            ya, yb = A[yk][:n] * B[gk][:n], B[yk][:n] * A[gk][:n]
            out[kk][f"{a_}-{b_}"] = dict(te=te, null=q1.nulls(ya, yb, te, q2.W_COMP, narrow=True),
                                         bands=sml.band_energies(dict(ya=ya, d=ya - yb), te))
    json.dump(out, open(path, "w"))
    nl = out["caps"][max(out["caps"])]
    print(f"{key}: R null plugin " + " ".join(f"{w}:{x[1]:.1f}" for w, x in nl["null_shift"].items() if x)
          + " | chain3 " + " ".join(f"{x[1]:.1f}" for w, x in nl["null_c3"].items() if x), flush=True)
    return key


def s_stage(workers, only=None, chunk=None):
    from concurrent.futures import ProcessPoolExecutor
    sys.path.insert(0, HERE)
    import pulse_check as pc
    from defects import DEFECTS
    q2 = pc.q2
    d = os.path.join(OUT, "s")
    os.makedirs(d, exist_ok=True)
    al = json.load(open(os.path.join(OUT, "align.json")))
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    folds = {rt: json.load(open(os.path.join(pc.CV, f"fold_{rt}.json"))) for rt in sorted({g["rt"] for g in info.values()})}
    keys = sorted(info)
    if chunk:
        i, m = (int(x) for x in chunk.split("/"))
        keys = keys[i::m]
    jobs = []
    for k in keys:
        g = info[k]
        path = os.path.join(d, k + ".json")
        if (only and k not in only) or os.path.exists(path):
            continue
        caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(g) if (k, a_) not in DEFECTS]
        if caps:
            jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], law,
                         folds[g["rt"]], cells.get(k), v, caps, ends, path, al))
    print(f"{len(jobs)} settings to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_s_job, jobs):
            pass



# ---------------------------------------------------------------- programs M and L
ML_AMPS = {"M": (10, 60, 99), "L": (0, 80)}
SIM_PRE = [0.0, 0.2008, 0.4016, 0.6024, 0.8032, 1.004]       # s of silence before the sweep: theta's phase at 0.83 Hz


def ml_groups(prog):
    sys.path.insert(0, HERE)
    import ml_carryover as cc
    g = cc.groups(prog)
    out = {}
    for k, files in g.items():
        n = max(sf.info(f).frames for f in files.values())
        out[k] = dict(len=n / HOST + 0.5, amps=sorted(files))
    return out


def ml_render_lines():
    L = []
    for prog in ("M", "L"):
        for k, g in ml_groups(prog).items():
            lf, rt, hf = (float(x) for x in k.split("_"))
            for amp in sorted(set(g["amps"]) | {0}):
                L.append(render_line(rpath(HOST, prog, k, amp), HOST, prog, lf, rt, hf, amp, g["len"]))
    return L


def ml_parts(prog, key, amps, n, chain3=False):
    """The plugin in ml_carryover's frame (device grid, 100 Hz-10 kHz, R's cluster tap 0 at K0, one shift for both
    channels from the highest AMPL's render): per AMPL the cluster part (render - AMPL-0 render: the same time variance,
    so exact) and the reverb part (the AMPL-0 render). chain3: the bridge's chain replaced by chain3 first."""
    import ml_carryover as cc
    import ml_structure as ml
    # the anchor as ml_structure._make finds a capture's (its onset: here R's first cluster tap), from the highest AMPL's
    # render, for every render of the group (they share the engine's timing)
    anchor = onset_anchor(load_render(rpath(HOST, prog, key, max(amps))))
    h = {}
    for a in sorted(set(amps) | {0}):
        x = to_device(load_render(rpath(HOST, prog, key, a)), anchor)
        if chain3:
            x = apply_ratio(x, ml.FSD, chain3_over_bridge)
        x = sosfiltfilt(ml.BAND, x, axis=0)
        h[a] = np.vstack([x, np.zeros((max(0, n - len(x)), 2))])[:n]
    s = cc.K0 - cc.tap0(h[max(amps)])
    h = {a: cc.fshift(x, s) for a, x in h.items()}
    return {a: h[a] - h[0] for a in amps}, h[0], float(s)


def _ml_job(a):
    prog, key, path = a
    import ml_carryover as cc
    import ml_lines as mll
    import pulse_check as pc
    K0, FS = cc.K0, cc.FS
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, _, amps = cc.captures(prog, key, n)
    n = len(caps[amps[-1]])
    W = mll.windows(prog)
    bench = {r["amp"]: r for r in json.load(open(mll.group_path(prog, key)))}
    out = dict(key=key, prog=prog, rows=[])
    variants = {}
    for vn, c3 in (("plugin", False), ("chain3", True)):
        C, R, s = ml_parts(prog, key, amps, n, c3)
        variants[vn] = (C, R)
        out[f"shift_{vn}"] = s
    t0 = (cc.ONSET[prog][1] + 20) / FS
    fw = (t0, t0 + 0.25)
    for amp in amps:
        y = caps[amp]
        end = mll.usable_end(y, g[amp], prog)
        sl, nf = mll.fit_window(prog, amp > 0)
        row = dict(amp=amp, end=[(e - K0) / FS for e in end], bench_law=bench[amp]["law"], bench_repeat=bench[amp]["repeat"],
                   bench_bands=bench[amp]["bands"])
        for vn, (C, R) in variants.items():
            cols = [C[amp], R] if amp > 0 else [R]
            mdl, u, gs, _ = mll.fit_warp(y, cols, sl, nf)
            row[f"{vn}_fit"] = dict(shift=u[0], ppm=u[1], gains=gs)
            row[vn] = {w: [mll.null_db(y[:, ch], mdl[:, ch], K0 + W[ch][w][0], min(K0 + W[ch][w][1], end[ch])) for ch in (0, 1)]
                       for w in W[0] if amp > 0 or w != "cluster_gap"}
            if vn == "plugin":
                row["bands"] = mll.band_decay(y, mdl, [K0 + cc.ONSET[prog][ch] for ch in (0, 1)], end)
                i1 = min(end)
                Y, Mp = pc.prep(y[:i1], dtype=np.float32), pc.prep(mdl[:i1], dtype=np.float32)
                row["m"] = pc.compare(Mp, Y, i1, K0, t0, fw)            # pulse_check's G2-G6 metrics, plugin vs capture
                row["tw"] = pc._ml_nulls(prog, y, mdl, i1)
        out["rows"].append(row)
    json.dump(out, open(path, "w"))
    r = out["rows"][-1]
    print(f"{prog} {key}: cluster_gap R plugin {r['plugin'].get('cluster_gap', [None, None])[1]} chain3 {r['chain3'].get('cluster_gap', [None, None])[1]}"
          f" bench {r['bench_law'].get('cluster_gap', [None, None])[1]}", flush=True)
    return key


def ml_stage(workers, only=None, chunk=None, every=1):
    from concurrent.futures import ProcessPoolExecutor
    sys.path.insert(0, HERE)
    import ml_carryover as cc
    from defects import ML_SLIPPED
    d = os.path.join(OUT, "ml")
    os.makedirs(d, exist_ok=True)
    jobs = []
    for prog in ("M", "L"):
        keys = sorted(cc.groups(prog))[::every]
        if chunk:
            i, m = (int(x) for x in chunk.split("/"))
            keys = keys[i::m]
        for k in keys:
            path = os.path.join(d, f"{prog}_{k}.json")
            if (only and f"{prog}_{k}" not in only) or os.path.exists(path):
                continue
            jobs.append((prog, k, path))
    print(f"{len(jobs)} groups to run (slipped captures stay in, flagged in the report: {len(ML_SLIPPED)})", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_ml_job, jobs):
            pass


# ---------------------------------------------------------------- M / L time variance through a sweep
def sweep_96k():
    """ml_variance's exponential sweep (10 Hz-13 kHz, L = 3.4 s, 50 ms sin^2 ends) at 96 kHz, and its inverse."""
    sys.path.insert(0, HERE)
    import ml_variance as mv
    f1, f2 = mv.SWEEP_F
    T = mv.SWEEP_L * np.log(f2 / f1)
    t = np.arange(int(T * HOST)) / HOST
    s = np.sin(2 * np.pi * f1 * mv.SWEEP_L * (np.exp(t / mv.SWEEP_L) - 1))
    m = int(0.05 * HOST)
    s[:m] *= np.sin(np.linspace(0, np.pi / 2, m)) ** 2
    s[-m:] *= np.sin(np.linspace(np.pi / 2, 0, m)) ** 2
    return s, s[::-1] * np.exp(-t / mv.SWEEP_L)


def sim_render_lines():
    sys.path.insert(0, HERE)
    import ml_variance as mv
    s, _ = sweep_96k()
    path = os.path.join(RENDERS, "sweep_96k.wav")
    os.makedirs(RENDERS, exist_ok=True)
    if not os.path.exists(path):
        sf.write(path, (0.25 * s).astype(np.float32), HOST, subtype="FLOAT")
    L = []
    for prog, keys in mv.SIM_KEYS.items():
        for k in keys:
            lf, rt, hf = (float(x) for x in k.split("_"))
            for i, pre in enumerate(SIM_PRE):
                L.append(render_line(rpath(HOST, prog, k, 0, f"_sweep{i}"), HOST, prog, lf, rt, hf, 0, len(s) / HOST + 6.0,
                                     pre=pre, inp=path))
    return L


def _sim_job(a):
    prog, key, path = a
    import ml_carryover as cc
    import ml_lines as mll
    import ml_structure as ml
    import ml_variance as mv
    from scipy.signal import fftconvolve
    al = json.load(open(os.path.join(OUT, "align.json")))
    s, inv = sweep_96k()
    ref = fftconvolve(s, inv)
    i0 = int(np.argmax(np.abs(ref)))
    irs = []
    for i in range(len(SIM_PRE)):
        x = load_render(rpath(HOST, prog, key, 0, f"_sweep{i}")) / 0.25
        h = np.stack([fftconvolve(x[:, c], inv)[i0:i0 + int(6.0 * HOST)] for c in (0, 1)], 1) / np.abs(ref).max()
        d = sosfiltfilt(ml.BAND, to_device(h, onset_anchor(load_render(rpath(HOST, prog, key, ML_AMPS[prog][-1])))), axis=0)
        irs.append(d)
    # one shift for all renders: the impulse render's cluster tap (the sweep renders have no cluster: AMPL 0)
    n = min(len(x) for x in irs)
    C, _, sh = ml_parts(prog, key, [ML_AMPS[prog][-1]], n)
    irs = [cc.fshift(x[:n], sh) for x in irs]
    fit = mll.fitted(prog, key)
    amps = list(fit)
    row = dict(key=key, tomkin=[], plugin=[])
    for i in range(len(amps)):
        for j in range(i + 1, len(amps)):
            e = min(min(fit[amps[i]][2]), min(fit[amps[j]][2]))
            row["tomkin"].append(mv.pair_nulls(prog, fit[amps[i]][0], fit[amps[j]][0], e))
    for i in range(len(irs)):
        for j in range(i + 1, len(irs)):
            row["plugin"].append(mv.pair_nulls(prog, irs[i], irs[j], n - 1000))
    for d_, name in (("sim", "python_sim"), ("sim_0.15", "python_sim_015")):   # theta rms 0.11 / 0.15 (the adopted, §9.7)
        sim = os.path.join(ROOT, "analysis", "ml_variance", d_, f"{prog}_{key}.json")
        if os.path.exists(sim):
            row[name] = json.load(open(sim))["sim"]
    json.dump(row, open(path, "w"))
    print(prog, key, flush=True)
    return key


def sim_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    sys.path.insert(0, HERE)
    import ml_variance as mv
    d = os.path.join(OUT, "sim")
    os.makedirs(d, exist_ok=True)
    jobs = [(p, k, os.path.join(d, f"{p}_{k}.json")) for p, ks in mv.SIM_KEYS.items() for k in ks]
    jobs = [j for j in jobs if not os.path.exists(j[2]) and (not only or f"{j[0]}_{j[1]}" in only)]
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_sim_job, jobs):
            pass


# ---------------------------------------------------------------- L's output offset (H item 2)
def fshift96(x, s):
    N = 1 << int(np.ceil(np.log2(len(x) + 4096)))
    return np.fft.irfft(np.fft.rfft(x, N) * np.exp(-2j * np.pi * np.fft.rfftfreq(N) * s), N)[:len(x)]


def _loff_job(a):
    key, amp, nm, prog = a
    from scipy.optimize import minimize_scalar
    band = butter(4, [100, 10000], btype="band", fs=HOST, output="sos")
    folder = {"S": "Short RT", "M": "Med RT", "L": "Long RT"}[prog]
    y, fs = sf.read(os.path.join(LIB, folder, nm), always_2d=True, dtype="float64")
    p = load_render(rpath(HOST, prog, key, amp))
    n = min(len(y), len(p), int(0.5 * HOST))
    y, p = sosfiltfilt(band, y[:n], axis=0), sosfiltfilt(band, p[:n], axis=0)
    # coarse lag on R over the first 300 ms, then the window: R's cluster tap 0 - 1 ms ... + 60 ms (capture samples)
    N = 1 << int(np.ceil(np.log2(2 * n)))
    c = np.fft.irfft(np.fft.rfft(y[:, 1], N) * np.conj(np.fft.rfft(p[:, 1], N)), N)
    lag = int(np.argmax(np.r_[c[:n], c[-n:]]))
    lag = lag if lag < n else lag - 2 * n
    a0 = onset_anchor(np.roll(p, lag, 0)) - int(0.001 * HOST)
    i0, i1 = max(a0, 0), a0 + int(0.06 * HOST)
    res = []
    for ch in (0, 1):
        def cost(s):
            z = fshift96(p[:, ch], lag + s)[i0:i1]
            g = np.dot(y[i0:i1, ch], z) / np.dot(z, z)
            return np.sum((y[i0:i1, ch] - g * z) ** 2) / np.sum(y[i0:i1, ch] ** 2)
        r = minimize_scalar(cost, bounds=(-3, 3), method="bounded", options=dict(xatol=1e-5))
        res.append((lag + float(r.x), float(10 * np.log10(r.fun))))
    dev = 26373.626373626374 / HOST
    d = (res[0][0] - res[1][0]) * dev                     # device samples by which the capture's L is later than the plugin's
    return dict(key=key, amp=amp, prog=prog, file=nm, shift=[r[0] for r in res], null=[r[1] for r in res],
                dL_minus_dR_dev=d, library_offset_dev=LR_DEV + d)


def loff_stage(workers):
    from concurrent.futures import ProcessPoolExecutor
    sys.path.insert(0, HERE)
    from defects import DEFECTS, ML_SLIPPED
    import ml_carryover as cc
    info = json.load(open(os.path.join(ROOT, "analysis", "q1", "data_chain3.json")))
    jobs = [(k, int(a), f["file"], "S") for k, g in info.items() if not k.startswith("_")
            for a, f in g["files"].items() if not a.startswith("_") and int(a) >= 30 and (k, int(a)) not in DEFECTS]
    for prog in ("M", "L"):
        for k, files in cc.groups(prog).items():
            for a, f in files.items():
                if a >= 60 and (prog, k, a) not in ML_SLIPPED:
                    jobs.append((k, a, os.path.basename(f), prog))
    with ProcessPoolExecutor(workers, max_tasks_per_child=50) as ex:
        rows = list(ex.map(_loff_job, jobs, chunksize=8))
    out = dict(rows=rows)
    for prog in ("S", "M", "L"):
        d = np.array([r["library_offset_dev"] for r in rows if r["prog"] == prog])
        out[prog] = dict(n=len(d), median=float(np.median(d)), iqr=[float(np.percentile(d, 25)), float(np.percentile(d, 75))],
                         us=float(np.median(d) / 26373.626373626374 * 1e6))
        print(prog, out[prog], flush=True)
    # what a 0.5 device-sample offset costs where the library's is the median above, L channel, if all else were exact
    T = 1 / 26373.626373626374
    dm = out["S"]["median"] - LR_DEV
    edges = [20, 160, 800, 1600, 3200, 6400, 10000, 12300]
    out["l_null_limit_db"] = {f"{a_}-{b_}": float(10 * np.log10(np.mean(np.abs(1 - np.exp(-2j * np.pi * np.linspace(a_, b_, 200) * dm * T)) ** 2)))
                              for a_, b_ in zip(edges[:-1], edges[1:])}
    print(out["l_null_limit_db"])
    json.dump(out, open(os.path.join(OUT, "loff.json"), "w"), indent=1)


# ---------------------------------------------------------------- program N (nl_model's bench)
N_EVERY = 4                                         # every 4th capture of the 1,044 (sorted keys), ~260
N_BENCH = os.path.join(ROOT, "analysis", "nl", "bench_warp2_bk_kf_a0.7_lb2_hf_g5")   # the spec's model: branch kernel on


def n_keys():
    sys.path.insert(0, HERE)
    import defects
    import nl_structure as nl
    return sorted(k for k in nl.table() if k not in defects.NL_DEFECTS)[::N_EVERY]


def n_render_lines():
    L = []
    for span, early, late in n_keys():
        k = f"{span}_{early}_{late}"
        L.append(render_line(rpath(HOST, "N", k, span), HOST, "N", early, 1.0, late, span, 3.0))
    return L


def n_device(path):
    """A 96 kHz N render on nl_structure's grid (its wide band window; file sample 0 = the input impulse at K0)."""
    sys.path.insert(0, HERE)
    import nl_structure as nl
    h = _harness()
    x = load_render(path)
    d = h["deconvolve"](np.vstack([np.zeros((h["PAD"], 2)), x]), HOST, h["K"])
    a = float(h["PAD"])
    k0 = int(np.floor(a / HOST * nl.FSD))
    y = np.stack([h["to_device_grid"](d[:, 0], HOST, a + h["ms"].LR_US * 1e-6 * HOST, k0, h["ms"].CFG),
                  h["to_device_grid"](d[:, 1], HOST, a, k0, h["ms"].CFG)], 1)
    return y[k0 - nl.K0:]


def _n_job(key):
    import nl_early as ne
    import nl_model as nm
    import nl_structure as nl
    import ml_lines as mll
    from scipy.optimize import minimize
    FS, K0 = nl.FSD, nl.K0
    span, early, late = key
    path = os.path.join(OUT, "n", f"{span}_{early}_{late}.json")
    if os.path.exists(path):
        return key
    noise = json.load(open(nm.NOISE))
    raw = nl.ir(nl.table()[key])
    h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = ne.front_shift(h)
    end = nm.usable_end(h)
    n = min(end, K0 + int(fr) + int(nm.WINDOWS[-1][1] * FS))
    i0, trimmed = nm.window_start(h, fr)
    y = h[:n]
    bench_f = os.path.join(N_BENCH, "_".join(str(v) for v in key) + ".json")
    out = dict(key=list(key), front_shift=fr, trimmed=trimmed, bench=json.load(open(bench_f)) if os.path.exists(bench_f) else None)
    pr = n_device(rpath(HOST, "N", f"{span}_{early}_{late}", span))
    for vn, c3 in (("plugin", False), ("chain3", True)):
        m0 = apply_ratio(pr, FS, chain3_over_bridge) if c3 else pr
        m0 = np.vstack([m0, np.zeros((max(0, n + 400 - len(m0)), 2))])[:n + 400]
        mb = sosfiltfilt(ne.BAND, m0, axis=0)
        # coarse: the cross-correlation peak of R over the scored span (front taps and early section), both banded
        L_ = min(n, K0 + int(fr) + int(0.6 * FS))
        Nf = 1 << int(np.ceil(np.log2(2 * L_)))
        xc = np.fft.irfft(np.fft.rfft(y[:L_, 1], Nf) * np.conj(np.fft.rfft(mb[:L_, 1], Nf)), Nf)
        lag = int(np.argmax(np.abs(np.r_[xc[:L_], xc[-L_:]])))
        s0 = float(lag if lag < L_ else lag - 2 * L_)       # delay applied to the plugin
        wins = []
        for t0, t1 in nm.WINDOWS:
            a_, b_ = max(i0, K0 + int(fr) + int(t0 * FS)), min(n, K0 + int(fr) + int(t1 * FS))
            if b_ - a_ >= int(0.02 * FS) and np.sum(h[a_:b_, 0] ** 2) >= 10 * noise["ch0"] * (b_ - a_):
                wins.append((a_, b_))

        def model(v):
            m = mll.warp(mb, s0 + v[0], v[1])[:n]
            g = [(y[:, c] @ m[:, c]) / (m[:, c] @ m[:, c]) for c in (0, 1)]
            return m * np.array(g), g

        def cost(v):                                       # the bench's cost: every window with signal counts equally
            m, _ = model(v)
            return float(sum(np.log(np.sum((y[a_:b_] - m[a_:b_]) ** 2) / np.sum(y[a_:b_] ** 2)) for a_, b_ in wins))
        best = min((cost([d, 11.6]), d) for d in np.arange(-2.0, 2.01, 0.1))[1]
        r = minimize(cost, [best, 11.6], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-7, maxiter=300,
                     initial_simplex=[[best, 11.6], [best + 0.05, 11.6], [best, 12.6]]))
        m, g = model(r.x)
        res = dict(shift=float(s0 + r.x[0]), ppm=float(r.x[1]), gain=[float(v) for v in g], windows={})
        for t0, t1 in nm.WINDOWS:
            a_, b_ = max(i0, K0 + int(fr) + int(t0 * FS)), min(n, K0 + int(fr) + int(t1 * FS))
            if b_ - a_ < int(0.02 * FS) or np.sum(y[a_:b_, 0] ** 2) < 10 * noise["ch0"] * (b_ - a_):
                continue
            res["windows"][f"{t0}-{t1}"] = {f"{k}{c}": v for c in (0, 1) for k, v in (
                ("null", float(10 * np.log10(np.sum((y[a_:b_, c] - m[a_:b_, c]) ** 2) / np.sum(y[a_:b_, c] ** 2)))),
                ("floor", float(10 * np.log10(noise[f"ch{c}"] * (b_ - a_) / np.sum(y[a_:b_, c] ** 2)))))}
        # per band over 0.05-0.6 s, full band (the unfiltered capture and plugin, the same fit)
        a_, b_ = max(i0, K0 + int(fr) + int(0.05 * FS)), min(n, K0 + int(fr) + int(0.6 * FS))
        yf = raw[:n] if len(raw) >= n else np.vstack([raw, np.zeros((n - len(raw), 2))])
        mfl = mll.warp(m0, s0 + r.x[0], r.x[1])[:n] * np.array(g)
        res["bands"] = {}
        if b_ - a_ > int(0.05 * FS):
            for lo, hi in nm.BANDS:
                sos = butter(4, [lo, hi], btype="band", fs=FS, output="sos")
                yb_, mb_ = sosfiltfilt(sos, yf, axis=0), sosfiltfilt(sos, mfl, axis=0)
                res["bands"][f"{lo}-{hi}"] = [float(10 * np.log10(np.sum((yb_[a_:b_, c] - mb_[a_:b_, c]) ** 2) / np.sum(yb_[a_:b_, c] ** 2))) for c in (0, 1)]
        out[vn] = res
    json.dump(out, open(path, "w"))
    w = out["plugin"]["windows"]
    print(key, " ".join(f"{t}:{v['null1']:.1f}" for t, v in w.items()), "| chain3", " ".join(f"{v['null1']:.1f}" for v in out["chain3"]["windows"].values()), flush=True)
    return key


def n_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    os.makedirs(os.path.join(OUT, "n"), exist_ok=True)
    keys = n_keys()[:only] if only else n_keys()
    with ProcessPoolExecutor(workers, max_tasks_per_child=10) as ex:
        for _ in ex.map(_n_job, keys, chunksize=2):
            pass


# ---------------------------------------------------------------- program E (tools/e250/e_pulse.py, Altiverb REV-250)
E_LF, E_HF, E_AMPL = 1.49, 0.69, 80.0               # the spec's adopted family values (e_plugin_spec controls)
E_SEED_PRE = 0.37                                    # s of silence per seed: the moving reads' phases at the impulse


def e_render_lines():
    sys.path.insert(0, os.path.join(HERE, "e250"))
    L = []
    for rt in (0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5):
        for seed in range(5):
            L.append(render_line(rpath("dev", "E", f"{rt}", seed), HOST, "E", E_LF, rt, E_HF, E_AMPL, min(3 * rt + 1, 15),
                                 pre=seed * E_SEED_PRE, domain="device"))
    return L


def e_stage(rts):
    """e_pulse's run with the plugin's E core (device-rate render, the plugin's own moving reads, seeds = start phases)
    in place of e_sim's render; 'imp' kind only. -> analysis/plugin_check/e/rt_<rt>.json, e_pulse.json."""
    sys.path.insert(0, os.path.join(HERE, "e250"))
    import e_pulse as ep
    al = 54                                          # ReverbEngine's alignment: E's own origin is not delayed (lateness 54)

    def render(kind, rt, seed, i1, kern):
        """The plugin's device-rate E (its output kernel [0.77, 0.23], e_sim's chain='device') moved to the bench's
        measurement side: that 2-tap kernel divided out, e_pulse's chain kernel (E.kernels: the Altiverb capture chain x
        L's output kernel, tap at index 16) applied, as e_sim.render(chain='altiverb') does."""
        from scipy.signal import lfilter
        x, _ = sf.read(rpath("dev", "E", f"{rt}", seed % 100), always_2d=True, dtype="float64")
        core = lfilter([1.0], [0.77, 0.23], x, axis=0)        # R's direct tap at k0 = 64 = e_pulse's K0
        y = np.stack([np.convolve(core[:, c], kern[c])[16:16 + len(core)] for c in (0, 1)], 1)
        n = i1 + int(0.2 * ep.FS)
        y = np.vstack([y, np.zeros((max(0, n - len(y)), 2))])[:n]
        return y / y[ep.K0, 1]
    ep.render = render
    ep.PARTS = os.path.join(OUT, "e")
    a = argparse.Namespace(rts=rts, kinds=["imp"], seeds=5, pairs=6, no_gain=True, redo=False,
                           out=os.path.join(OUT, "e_pulse_plugin.json"))
    log = lambda m: print(m, flush=True)
    ep.stage_run(a, log)
    ep.stage_report(a, log)


# ---------------------------------------------------------------- host rates
RATE_SET = [("S", "1.0_1.0_1.0", 90), ("S", "0.1_2.5_0.2", 30), ("M", "1.2_5.0_0.9", 60), ("L", "1.2_4.5_0.8", 80),
            ("E", "1.49_2.0_0.69", 80), ("N", "2.0_1.0_1.0", 0.3)]


def rates_render_lines():
    L = []
    for rate in (44100, 48000, 96000):
        for prog, k, amp in RATE_SET:
            lf, rt, hf = (float(x) for x in k.split("_"))
            for mod in (0, 1):
                L.append(render_line(rpath(rate, prog, k, amp, f"_mod{mod}"), rate, prog, lf, rt, hf, amp, 2.0, mod=mod))
    return L


def rates_stage():
    """Each program's 44.1 / 48 kHz render against its 96 kHz render: band-limited resampling to 96 kHz, one gain and a
    fractional shift per channel (the latency's uncompensated fraction), nulls per band 20 Hz-12.8 kHz over 0-1 s."""
    from fractions import Fraction
    from scipy.optimize import minimize_scalar
    from scipy.signal import resample_poly
    edges = [20, 40, 80, 160, 400, 800, 1600, 3200, 6400, 10000, 12000, 12800]
    out = {}
    for prog, k, amp in RATE_SET:
        for mod in (0, 1):
            ref, _ = sf.read(rpath(96000, prog, k, amp, f"_mod{mod}"), always_2d=True, dtype="float64")
            for rate in (44100, 48000):
                x, _ = sf.read(rpath(rate, prog, k, amp, f"_mod{mod}"), always_2d=True, dtype="float64")
                fr = Fraction(96000, rate)
                xr = resample_poly(x, fr.numerator, fr.denominator, axis=0, window=("kaiser", 14.0))
                n = 96000                                   # 0-1 s scored; the 2 s renders leave room for the shift
                m_ = min(len(xr), len(ref))
                xr, r = xr[:m_], ref[:m_]
                lat = {44100: 111, 48000: 117}[rate] * 96000 / rate - 211      # reported latencies, in 96 kHz samples
                row = {"latency_mismatch_96k_samples": lat}
                for ch in (0, 1):
                    bb = butter(4, [20, 12800], btype="band", fs=96000, output="sos")
                    rb, xb = sosfiltfilt(bb, r[:, ch]), sosfiltfilt(bb, xr[:, ch])
                    def cost(s_):
                        z = fshift96(xb, s_)[:n]
                        g = np.dot(rb[:n], z) / np.dot(z, z)
                        return np.sum((rb[:n] - g * z) ** 2) / np.sum(rb[:n] ** 2)
                    sh = minimize_scalar(cost, bounds=(-lat - 3, -lat + 3), method="bounded", options=dict(xatol=1e-6)).x
                    z = fshift96(xr[:, ch], sh)
                    zb = sosfiltfilt(bb, z)[:n]
                    g = np.dot(rb[:n], zb) / np.dot(zb, zb)
                    bands = {}
                    for lo, hi in zip(edges[:-1], edges[1:]):
                        sos = butter(3, [lo, hi], btype="band", fs=96000, output="sos")
                        e_, r_ = sosfiltfilt(sos, r[:, ch] - g * z)[:n], sosfiltfilt(sos, r[:, ch])[:n]
                        bands[f"{lo}-{hi}"] = float(10 * np.log10(np.sum(e_ ** 2) / np.sum(r_ ** 2)))
                    row["LR"[ch]] = dict(shift_96k=float(sh), gain_db=float(20 * np.log10(g)), bands=bands,
                                         broadband=float(10 * np.log10(cost(sh))))
                out[f"{prog}_{k}_mod{mod}_{rate}"] = row
                print(prog, k, mod, rate, "broadband L/R", round(row["L"]["broadband"], 1), round(row["R"]["broadband"], 1),
                      "worst band", round(max(max(row[c]["bands"].values()) for c in "LR"), 1), "shift", round(row["R"]["shift_96k"], 3), flush=True)
    json.dump(out, open(os.path.join(OUT, "rates.json"), "w"), indent=1)


# ---------------------------------------------------------------- report
def _db(e, y):
    return float(10 * np.log10(e / y)) if e > 0 and y > 0 else None


def _med(v):
    v = [x for x in v if x is not None]
    return round(float(np.median(v)), 2) if v else None


def report_s():
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "s", "*.json")))]
    W = ["0-50", "50-150", "150-300", "300-600", "600-1000", "1000-1500", "1500-2500"]
    out = dict(settings=len(rows), captures=sum(len(r["caps"]) for r in rows), pulse={}, fullband={}, corners={})
    for w in W:
        cell = {}
        for kname in ("null", "null_shift", "null_c3", "null_bench", "null_netlist_hp"):
            for ch, cn in ((0, "L"), (1, "R")):
                cell[f"{kname}_{cn}"] = _med([c[kname][w][ch] for r in rows for c in r["caps"].values() if c[kname].get(w)])
        rep = [p["null"][w][1] for r in rows for p in r["pairs_bench"].values() if p["null"].get(w)]
        cell["repeat_R"] = _med(rep)
        for kname in ("null_shift", "null_c3", "null_bench"):
            v = [c[kname][w][1] for r in rows for c in r["caps"].values() if c[kname].get(w)]
            cell[f"{kname}_R_le_-50"] = f"{sum(x <= -50 for x in v)}/{len(v)}"
            if cell["repeat_R"] is not None and v:
                cell[f"{kname}_R_excess"] = round(float(np.median(v)) - cell["repeat_R"], 2)
        out["pulse"][w] = cell
    pairs = [("plugin", "e", "y"), ("plugin_shift", "es", "ys"), ("plugin_chain3", "e3", "y3"), ("bench", "eb0", "yb"),
             ("bench_lowband", "eb1", "yb"), ("plugin_netlist_hp", "en", "y"), ("a_plugin_vs_bench", "d0", "m0"),
             ("a_plugin_vs_bench_lowband", "d1", "m0")]
    bands = list(rows[0]["caps"][next(iter(rows[0]["caps"]))]["bands"])
    for b in bands:
        out["fullband"][b] = {}
        for w in ["0-50", "50-150", "150-300", "300-600", "600-1000"]:
            cell = {}
            for name, e, y in pairs:
                for ch, cn in ((0, "L"), (1, "R")):
                    cell[f"{name}_{cn}"] = _med([_db(c["bands"][b][w][e][ch], c["bands"][b][w][y][ch]) for r in rows
                                                 for c in r["caps"].values() if w in c["bands"][b]])
            cell["repeat_R"] = _med([_db(p["bands"][b][w]["d"][1], p["bands"][b][w]["ya"][1]) for r in rows
                                     for p in r["pairs_bench"].values() if w in p["bands"][b]])
            out["fullband"][b][w] = cell
    for b in ("20-40", "40-80", "80-160", "160-270"):        # H item 1: netlist corners minus bridge corners, per capture
        d = [_db(c["bands"][b][w]["en"][ch], c["bands"][b][w]["y"][ch]) - _db(c["bands"][b][w]["e"][ch], c["bands"][b][w]["y"][ch])
             for r in rows for c in r["caps"].values() for w in ("50-150", "150-300", "300-600") if w in c["bands"][b] for ch in (0, 1)]
        d3 = []
        out["corners"][b] = dict(median_db=_med(d), better=f"{sum(x < 0 for x in d)}/{len(d)}", range=[round(min(d), 3), round(max(d), 3)])
    out["clock_ppm_shift_fit"] = _med([c["fit_shift"]["eps"] for r in rows for c in r["caps"].values()])
    out["shift_fit_device_samples"] = _med([c["fit_shift"]["shift"] for r in rows for c in r["caps"].values()])
    return out


def report_ml():
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "ml", "*.json")))]
    out = {}
    for prog in ("M", "L"):
        rr = [x for r in rows if r["prog"] == prog for x in r["rows"]]
        if not rr:
            continue
        o = dict(groups=sum(r["prog"] == prog for r in rows), captures=len(rr), windows={}, bands={})
        for w in ("cluster_gap", "first_pass", "recirc1", "recirc2", "recirc3"):
            cell = {}
            for kname in ("plugin", "chain3", "bench_law", "bench_repeat"):
                for ch, cn in ((0, "L"), (1, "R")):
                    cell[f"{kname}_{cn}"] = _med([x[kname][w][ch] for x in rr if w in x[kname] and x[kname][w][ch] is not None])
                v = [x[kname][w][1] for x in rr if w in x[kname] and x[kname][w][1] is not None]
                cell[f"{kname}_R_le_-50"] = f"{sum(t <= -50 for t in v)}/{len(v)}"
            o["windows"][w] = cell
        for f in ("125", "250", "500", "1000", "2000", "4000", "8000"):
            for kname, key in (("plugin", "bands"), ("bench", "bench_bands")):
                sl = [x[key][f][ch][0] for x in rr for ch in (0, 1) if x[key].get(f) and x[key][f][ch]]
                lv = [x[key][f][ch][1] for x in rr for ch in (0, 1) if x[key].get(f) and x[key][f][ch]]
                o["bands"].setdefault(f, {})[kname] = dict(slope_db_per_s_median=_med(sl), slope_abs_p90=round(float(np.percentile(np.abs(sl), 90)), 2) if sl else None,
                                                          level_db_median=_med(lv), level_abs_p90=round(float(np.percentile(np.abs(lv), 90)), 2) if lv else None)
        # EDC decay-time ratio (model / capture - 1) per octave, R, plugin vs capture (pulse_check G2)
        o["edc_t_ratio"] = {}
        for f in ("125", "250", "500", "1000", "2000", "4000", "8000"):
            v = [x["m"]["oct"][f][1]["t"][0] / x["m"]["oct"][f][1]["t"][1] - 1 for x in rr if x["m"]["oct"].get(f) and x["m"]["oct"][f][1]["t"]]
            o["edc_t_ratio"][f] = dict(median=_med(v), abs_p90=round(float(np.percentile(np.abs(v), 90)), 4) if v else None)
        out[prog] = o
    sims = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "sim", "*.json")))]
    if sims:
        out["decorrelation"] = {}
        for r in sims:
            row = {}
            for src in ("tomkin", "plugin", "python_sim", "python_sim_015"):
                if src in r:
                    row[src] = {w: _med([p[w] for p in r[src] if p.get(w) is not None]) for w in r["plugin"][0]}
            out["decorrelation"][r["key"] if "prog" not in r else r["key"]] = row
    return out


def report_n():
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "n", "*.json")))]
    o = dict(captures=len(rows), windows={}, bands={})
    for w in ("0.0-0.05", "0.05-0.15", "0.15-0.3", "0.3-0.6", "0.6-1.0", "1.0-2.0"):
        cell = {}
        for vn in ("plugin", "chain3", "bench"):
            src = [(r[vn] if vn != "bench" else r["bench"]) for r in rows if (vn != "bench" or r["bench"])]
            for c in (0, 1):
                cell[f"{vn}_null{c}"] = _med([x["windows"][w][f"null{c}"] for x in src if w in x["windows"]])
            v = [x["windows"][w]["null1"] - x["windows"][w]["floor1"] for x in src if w in x["windows"]]
            cell[f"{vn}_excess1"] = _med(v)
            cell[f"{vn}_R_le_-50"] = f"{sum(x['windows'][w]['null1'] <= -50 for x in src if w in x['windows'])}/{sum(w in x['windows'] for x in src)}"
        cell["floor1"] = _med([r["plugin"]["windows"][w]["floor1"] for r in rows if w in r["plugin"]["windows"]])
        o["windows"][w] = cell
    for b in ("20-100", "100-1000", "1000-5000", "5000-10000", "10000-12000"):
        o["bands"][b] = {vn: _med([(r[vn] if vn != "bench" else r["bench"])["bands"][b][1] for r in rows
                                   if (vn != "bench" or r["bench"]) and b in (r[vn] if vn != "bench" else r["bench"]).get("bands", {})])
                         for vn in ("plugin", "chain3", "bench")}
    return o


def report_stage():
    out = {}
    for name, fn in (("s", report_s), ("ml", report_ml), ("n", report_n)):
        try:
            out[name] = fn()
        except (IndexError, KeyError, ValueError, StopIteration) as e:
            out[name] = dict(error=repr(e))
        json.dump(out[name], open(os.path.join(OUT, f"summary_{name}.json"), "w"), indent=1)
    for f in ("chain.json", "loff.json", "rates.json", "align.json"):
        if os.path.exists(os.path.join(OUT, f)):
            print(f, "present")
    print(json.dumps({k: (v if k != "s" else {kk: vv for kk, vv in v.items() if kk != "fullband"}) for k, v in out.items()}, indent=1)[:20000])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--set", default="all")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--chunk")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--rts", type=float, nargs="*", default=[0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5])
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.stage == "render":
        render_stage(a.set)
    elif a.stage == "align":
        align_stage()
    elif a.stage == "ml":
        ml_stage(a.workers, a.only, a.chunk, a.every)
    elif a.stage == "sim":
        sim_stage(a.workers, a.only)
    elif a.stage == "n":
        n_stage(a.workers, int(a.only[0]) if a.only else None)
    elif a.stage == "e":
        e_stage(a.rts)
    elif a.stage == "rates":
        rates_stage()
    elif a.stage == "report":
        report_stage()
    elif a.stage == "loff":
        loff_stage(a.workers)
    elif a.stage == "chain":
        chain_stage()
    elif a.stage == "s":
        s_stage(a.workers, a.only, a.chunk)
    else:
        raise SystemExit(f"unknown stage {a.stage}")


if __name__ == "__main__":
    main()
