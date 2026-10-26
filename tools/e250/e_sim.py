"""Program E in the time domain with moving reads, and simulated sweep captures of it (docs/EMT250_program.md §12).

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_sim.py --stage data      # statistics of the presets
    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_sim.py --stage sim [--json params]
    .venv/bin/python tools/e250/e_sim.py --stage render --json '{"rt": 2.5, "fs": 48000, "out": "e_2.5.wav"}'   # reference IR

The static network is e_model's (§11.2, the 246 law at LF 1.44 / HF 0.72). The time-variant part (§11.3) is M reads per loop,
each at an integer output delay o(t) (as the 250's addresses: no interpolation) moving between lo and hi, up at v and down at vdn samples per second, a gain +-b (relative to the loop's static read), added to R and, with probability 'share', to L with the opposite sign
(otherwise an independent read for L). Motion laws: 'tri' (bounces between lo and hi), 'sine' (centre in [lo, hi], depth A, rate f), 'walk' (constant speed, direction
re-drawn every 'tdir' seconds). A capture is the response to an exponential sweep (sweep seconds 'tsw', unknown for the Altiverb
captures, a nuisance) deconvolved by its inverse filter, as the library captures were, so the moving reads are smeared the same way.

Statistics (stats()): the random part r = capture - static model, against the static part s:
  e_ratio  random / static energy (dB) in 300-805 / 805-2215 / 2215-4400 / 4400-9000 samples after R's tap 1 (R)
  lobes    the leave-one-out residual's decay-normalised autocorrelation at P_k + (-12..12), mean over the three loops (R)
  lr       correlation of R's and L's leave-one-out residuals at the static alignment (315) and +-1, +-2
  frac     static fraction per octave, 805-9000 samples: rho^2 of the capture against the static model (R)
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import lfilter, butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "tools"))
import e_model as E  # noqa: E402

FS, K0 = E.FS, E.K0
X_STATIC = None                      # e_model parameters (246 law), set by static_params()
DEF = dict(M=6, b=0.45, lo=600, hi=3600, v=60.0, vsd=0.0, vdn=None, f=0.55, A=17.0, fsd=0.0, wave="tri", groups=0, gs=0.17, gw=0.1, Mw=6, dwell=0.1, sets=True, common=False, law="tri", tdir=1.0, share=0.5, tsw=10.0, seed=0, mdel=314)
OCTS = [250, 500, 1000, 2000, 4000]


# The adopted moving-read generator (EMT250_program §16): E's own fixed reads (from its low band), triangle glide, independent phases
ADOPTED = dict(law="fixed", f=0.553, A=32.0, fsd=0.5, common=False, mdel=314, dl_same=1, mlp=0.07)   # §16: E's nine fixed reads; same-sign reads reach L at 315
ADOPTED_13 = dict(law="lfo", wave="tri", f=0.553, A=32.0, fsd=0.3, groups=0, b=0.81, lo=410, hi=3600, share=0.5, M=6, mdel=314)  # §13.4, superseded


LFO_RATE = {"time": 0.553, "samples": 0.553 * 26373.6 / 24000}   # §14.6: the 250's 0.553 Hz kept in seconds (default) or in samples
L_KERNEL = (0.77, 0.23)                    # §14.1: L's output kernel at the device rate (E has it; Tomkin's L p = 0.23)
AMPL_REF = 80.0                            # the Altiverb family's AMPL (H: as the S / M families, library_notes §6)


def lowpass_sos(khz):
    """The LOW PASS lever (2-12 kHz, acting on the whole output, operation card): open at 12 kHz (the analogue band limit is
    the 12.2 kHz input / output filters); below it a 2nd-order Butterworth at the device rate. H: the lever's response has
    never been measured (no LP variation in any existing capture, §14.1)."""
    from scipy.signal import butter
    return None if khz >= 12 else butter(2, khz * 1e3, fs=FS, output="sos")


def render_ir(rt, lf=None, hf=None, seed=0, secs=None, fs_out=None, ampl=AMPL_REF, delay_ms=0.0, lp_khz=12.0, lfo="time",
              chain="altiverb", law=None, moving=True):
    """Reference impulse response of Program E, (n, 2) = (L, R), with every control the 246 panel gives E:
      rt, lf, hf  the 246's loop law (ml_lines.line_coefs); lf / hf default to the Altiverb family's fit (1.38 / 0.74)
      ampl        REFLECTIONS AMPL. %: the direct taps x ampl / 80 (S's linear cluster law, model_s_results §21; H for E)
      delay_ms    REFLECTIONS DELAY: a pure delay of the whole output, rounded to device samples (operation card; H for E)
      lp_khz      LOW PASS: see lowpass_sos (H)
      lfo         'time' (0.553 Hz, default) or 'samples' (0.608 Hz): the moving reads' triangle rate (H)
      chain       'altiverb': the measured chain kernel of the E presets (the capture's chain x L's output kernel; validated
                  against them); 'device': L's output kernel [0.77, 0.23] only, for a host that adds the 246's analogue chain
      law         the moving-read law (None: ADOPTED, §16; 'random': §13.4's random reads; 'plugin': §14.4's plugin-faithful structure)
      moving      False: the static network only (deterministic; the plugin's reference for a sample-exact null)
    Device rate, or resampled to fs_out."""
    global X_STATIC
    static_params(); q = X_STATIC.copy()
    if lf is not None: q[9] = lf
    if hf is not None: q[10] = hf
    kern = E.kernels([r for r in E.RTS if r >= 1.0])
    if chain == "device":
        k = np.zeros(len(kern[0])); k[16:18] = L_KERNEL; kern = [k, k.copy()]
    prm = dict(DEF) | ADOPTED | dict(f=LFO_RATE[lfo], ampl_scale=ampl / AMPL_REF)
    if law == "random":
        prm = dict(DEF) | ADOPTED_13 | dict(f=LFO_RATE[lfo], ampl_scale=ampl / AMPL_REF)
    if law == "plugin":
        prm |= dict(law="plugin", f=LFO_RATE[lfo] + 0.0005, A=35.2, gs=0.17, gw=0.28, Mw=6, lo=300, dwell=0.1, sets=True)
    if not moving:
        prm |= dict(b=0.0, gs=0.0, gw=0.0, mscale=0.0)
    secs = secs or min(3 * rt + 1, 15); n = int(secs * FS) + 4096
    x = np.zeros(n); x[1000] = 1.0
    keep = X_STATIC; X_STATIC = q
    try:
        y, _ = render(x, rt, prm, kern, np.random.default_rng(seed))
    finally:
        X_STATIC = keep
    y = y[1000 + 16 - K0:1000 + 16 - K0 + int(secs * FS)]
    d = int(round(delay_ms * 1e-3 * FS))
    if d:
        y = np.vstack([np.zeros((d, 2)), y])[:len(y)]
    sos = lowpass_sos(lp_khz)
    if sos is not None:
        from scipy.signal import sosfilt
        y = sosfilt(sos, y, axis=0)
    if fs_out:
        from scipy.signal import resample_poly
        from fractions import Fraction
        fr = Fraction(fs_out / FS).limit_denominator(2000); y = resample_poly(y, fr.numerator, fr.denominator, axis=0)
    return y


def static_params():
    """All-data static fit with the 246 law (cached)."""
    global X_STATIC
    f = os.path.join(E.OUT, "static_all_246.json")
    if not os.path.exists(f):
        E.LAWMODE = "246"; E.X0 = np.r_[E.X0[:9], 1.44, 0.72]
        x, kern, _ = E.fitstatic([r for r in E.RTS if r >= 1.0], 1.0)
        json.dump(dict(x=x.tolist()), open(f, "w"))
    E.LAWMODE = "246"
    X_STATIC = np.array(json.load(open(f))["x"])
    return X_STATIC


def loop_iir(rt, Pk, lf, hf):
    import ml_lines as mll
    g, r0, p = mll.line_coefs([Pk], rt, lf, hf)[0]
    a = mll.LAW["a"]; b0, b1 = (r0 * (1 - a) + (1 + a)) / 2, (r0 * (1 - a) - (1 + a)) / 2
    return g * (1 - p) * np.array([b0, b1]), np.convolve([1, -a], [1, -p])


def trajectories(n, prm, rng):
    """Integer output delays (n_reads, n) for each loop's M moving reads, and their gains / L signs."""
    M, lo, hi = prm["M"], prm["lo"], prm["hi"]
    t = np.arange(n) / FS; out = []
    PH = rng.uniform(0, 1, 8)
    for _ in range(3 * M):
        o0 = rng.uniform(lo, hi); d = rng.choice([-1, 1]); v = prm["v"] * rng.uniform(1 - prm["vsd"], 1 + prm["vsd"])
        if prm["law"] == "tri":         # up at v, down at vdn (asymmetric triangle; vdn None: symmetric)
            span = hi - lo; vd = v * (prm["vdn"] / prm["v"]) if prm["vdn"] else v
            Tu, Td = span / v, span / vd; ph = (rng.uniform(0, Tu + Td) + t) % (Tu + Td)
            o = np.where(ph < Tu, lo + v * ph, hi - vd * (ph - Tu))
        elif prm["law"] == "lfo":        # the UAD plugin's law (EMT250_program §13): triangle (or sine) LFO, depth A, rate f,
            #                              phase common to all reads ('groups' phase groups, reads assigned in turn) or independent
            g = (len(out) % prm["groups"]) if prm["groups"] else None
            ph = PH[g] if g is not None else rng.uniform(0, 1)
            x = (prm["f"] * rng.uniform(1 - prm["fsd"], 1 + prm["fsd"]) * t + ph) % 1    # fsd: per-read rate (speed) spread
            wav = (4 * np.abs(x - 0.5) - 1) if prm["wave"] == "tri" else np.cos(2 * np.pi * x)
            o = o0 + prm["A"] * wav
        elif prm["law"] == "sine":       # the UAD 250's law (UAD250_test_plan §2): centre uniform in [lo, hi], depth A, rate f (+- fsd)
            f = prm["f"] * rng.uniform(1 - prm["fsd"], 1 + prm["fsd"])
            o = o0 + prm["A"] * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
        else:                            # walk: constant speed, direction re-drawn every tdir s, reflected at lo / hi
            k = int(np.ceil(t[-1] / prm["tdir"])) + 1
            dirs = rng.choice([-1, 1], k); seg = np.minimum((t / prm["tdir"]).astype(int), k - 1)
            o = o0 + np.cumsum(dirs[seg]) * v / FS
            span = hi - lo; ph = (o - lo) % (2 * span); o = lo + np.where(ph < span, ph, 2 * span - ph)
        out.append(np.round(o).astype(int))
    return np.array(out)


# The UAD plugin's moving reads (EMT250_program §14, `uad_moving.py`), rescaled to E's clock (x 1.0989): (centre, loop, R sign, L sign, set)
PLUGIN_READS = ((635, 0, -1, +1, "Z"), (1199, 2, -1, +1, "Y"), (1412, 1, +1, +1, "X"), (1458, 0, +1, -1, "Y"),
                (1879, 0, -1, -1, "X"), (2078, 1, +1, -1, "Z"), (2391, 2, -1, -1, "X"))
SETS = {"X": 0.62, "Y": 0.29, "Z": 0.09}


def plugin_reads(n, prm, rng, lines, idx):
    """§14's law: every read on one common triangle g(t) (rate f, depth A); the strong reads switch on and off by set (one set
    drawn per dwell segment), gain gs (absolute, direct tap = 1); plus Mw weak reads per loop, always on, gain gw, centres in
    [lo, hi], random signs, L sign shared (opposite) with probability share."""
    t = np.arange(n) / FS; x = (prm["f"] * t + rng.uniform()) % 1; g = prm["A"] * (4 * np.abs(x - 0.5) - 1)
    seg = (t / prm["dwell"]).astype(int); pick = rng.choice(list(SETS), seg.max() + 1, p=list(SETS.values()))[seg]
    movR = np.zeros(n); movL = np.zeros(n)
    def read(line, c):
        j = idx - np.round(c + g).astype(int); return np.where(j >= 0, line[np.maximum(j, 0)], 0)
    for c, k, sr, sl, st in PLUGIN_READS:
        v = read(lines[k], c) * prm["gs"] * (pick == st if prm["sets"] else 1)
        movR += sr * v; movL += sl * v
    for k in range(3):
        for _ in range(prm["Mw"]):
            v = read(lines[k], rng.uniform(prm["lo"], prm["hi"])) * prm["gw"] * rng.choice([-1, 1])
            movR += v
            movL += -v if rng.uniform() < prm["share"] else read(lines[k], rng.uniform(prm["lo"], prm["hi"])) * prm["gw"] * rng.choice([-1, 1])
    return movR, movL


def fixed_reads():
    """§16: E's own moving reads (fixed centres, loops, R and L amplitudes; analysis/e250/fixed_reads.json, from the low band)."""
    return json.load(open(os.path.join(E.OUT, "fixed_reads.json")))["reads"]


def e_reads(n, prm, rng, lines, idx):
    """§16's law: E's fixed reads, each at centre + A tri(f_j t + phase_j); one phase for all reads (common=True) or
    independent ones, rates f (1 +- fsd); L reads the same address with its own signed amplitude."""
    t = np.arange(n) / FS; movR = np.zeros(n); movL = np.zeros(n); ph0 = rng.uniform()
    for r in fixed_reads():
        f = prm["f"] * (rng.uniform(1 - prm["fsd"], 1 + prm["fsd"]) if prm["fsd"] else 1.0)
        x = (f * t + (ph0 if prm["common"] else rng.uniform())) % 1
        wav = np.cos(2 * np.pi * x) if prm.get("wave") == "sine" else (4 * np.abs(x - 0.5) - 1)
        o = np.round(r["centre"] + prm["A"] * wav).astype(int)
        j = idx - o; v = np.where(j >= 0, lines["ABC".index(r["loop"])][np.maximum(j, 0)], 0)
        ms = prm.get("mscale", 1.0)              # 0 = moving reads off (render_ir moving=False)
        movR += ms * r["R"] * v; msL = ms * prm.get("mscale_L", 1.0)
        dl = prm.get("dl_same", 0) if r["L"] * r["R"] > 0 else prm.get("dl_opp", 0)   # extra L delay per sign group (§16.7)
        if dl:
            jl = j - dl; v = np.where(jl >= 0, lines["ABC".index(r["loop"])][np.maximum(jl, 0)], 0)
        movL += msL * r["L"] * v / X_STATIC[2]     # L amplitudes include L's gain (applied again to all of L below)
    if prm.get("mlp"):                            # a one-pole low-pass on the moving reads (§17: low band exact, top trimmed)
        p_ = prm["mlp"]; movR = lfilter([1 - p_], [1, -p_], movR); movL = lfilter([1 - p_], [1, -p_], movL)
    return movR, movL


def render(x, rt, prm, kern, rng):
    """(n, 2) device-rate output of E for a mono input x (R's tap 1 at the input's own time + K0 - 16)."""
    q = E.unpack(X_STATIC); n = len(x)
    # direct path
    direct = x.copy(); direct[429:] += q["a"][0] * x[:-429]; direct[718:] += q["a"][1] * x[:-718]
    direct *= prm.get("ampl_scale", 1.0)                 # REFLECTIONS AMPL. (render_ir)
    # all-pass chain
    u = x.copy()
    for D, g in zip(E.AP_D, q["g"]):
        bb = np.zeros(D + 1); bb[0] = -g; bb[D] = 1; aa = np.zeros(D + 1); aa[0] = 1; aa[D] = -g
        u = lfilter(bb, aa, u)
    # loops: v[t] = u[t] - line[t - P], line = G{v}, block-wise (block < P); static reads read line (after the loop filter),
    # E's fixed moving reads read v (before it, §17)
    lines = []; vlines = []
    for Pk, (lf_, hf_) in zip(E.P, [tuple(q["m"])] * 3):
        bG, aG = loop_iir(rt, Pk, lf_, hf_); line = np.zeros(n); v = np.zeros(n); zi = np.zeros(max(len(aG), len(bG)) - 1); B = 2048
        for t0 in range(0, n, B):
            t1 = min(n, t0 + B); fb = np.zeros(t1 - t0)
            s0 = t0 - Pk
            if t1 - Pk > 0:
                a0 = max(s0, 0); fb[a0 - s0:] = line[a0:t1 - Pk]
            v[t0:t1] = u[t0:t1] - fb
            line[t0:t1], zi = lfilter(bG, aG, v[t0:t1], zi=zi)
        lines.append(line); vlines.append(v)
    stat = np.zeros(n)
    for line, dk, ck in zip(lines, E.DK, q["c"]):
        stat[dk:] += ck * line[:-dk]
    # moving reads
    if prm["law"] == "fixed":
        idx = np.arange(n); movR, movL = e_reads(n, prm, rng, vlines if prm.get("pre", True) else lines, idx)
    elif prm["law"] == "plugin":
        idx = np.arange(n); movR, movL = plugin_reads(n, prm, rng, lines, idx)
    else:
        movR, movL = generic_reads(n, prm, rng, lines, q)
    R = direct + stat + movR
    L = np.zeros(n)
    md = prm["mdel"]
    L[E.LDEL:] = -direct[:-E.LDEL]; L[E.LREV:] += stat[:-E.LREV]; L[md:] += movL[:-md]; L *= q["gL"]
    out = np.stack([np.convolve(L, kern[0])[:n], np.convolve(R, kern[1])[:n]], 1)
    return out, (direct, stat)


def generic_reads(n, prm, rng, lines, q):
    trR = trajectories(n, prm, rng); trL = trajectories(n, prm, rng)
    sg = rng.choice([-1, 1], len(trR)); shared = rng.uniform(size=len(trR)) < prm["share"]
    movR = np.zeros(n); movL = np.zeros(n); idx = np.arange(n)
    for i in range(len(trR)):
        k = i // prm["M"]; ck = abs(q["c"][k]); line = lines[k]
        jR = idx - trR[i]; movR += np.where(jR >= 0, line[np.maximum(jR, 0)], 0) * sg[i] * prm["b"] * ck
        if shared[i]:
            movL -= np.where(jR >= 0, line[np.maximum(jR, 0)], 0) * sg[i] * prm["b"] * ck
        else:
            jL = idx - trL[i]; movL += np.where(jL >= 0, line[np.maximum(jL, 0)], 0) * sg[i] * prm["b"] * ck
    return movR, movL


def sweep(tsw, f1=20.0, f2=12500.0):
    t = np.arange(int(tsw * FS)) / FS; Lr = tsw / np.log(f2 / f1)
    s = np.sin(2 * np.pi * f1 * Lr * (np.exp(t / Lr) - 1))
    fade = int(0.01 * FS); s[:fade] *= np.hanning(2 * fade)[:fade]; s[-fade:] *= np.hanning(2 * fade)[fade:]
    return s


def capture(rt, prm, kern, tail=2.5):
    """Simulated sweep capture: (n, 2) IR with R's tap 1 at K0, and the static model rendered on the same grid."""
    rng = np.random.default_rng(prm["seed"])
    s = sweep(prm["tsw"]); lead = 4096; n = lead + len(s) + int(tail * FS)
    x = np.zeros(n); x[lead:lead + len(s)] = s
    y, _ = render(x, rt, prm, kern, rng)
    N = 1 << int(np.ceil(np.log2(n + len(s))))
    S = np.fft.rfft(s, N); inv = np.conj(S) / (np.abs(S) ** 2 + 1e-6 * np.max(np.abs(S)) ** 2)
    h = np.fft.irfft(np.fft.rfft(y, N, axis=0) * inv[:, None], N, axis=0)
    # the chain kernel puts tap 1's peak at index 16 of the kernel; the input starts at lead
    i0 = lead + 16 - K0
    return h[i0:i0 + K0 + int(tail * FS)]


def stats(Ys, Ss):
    """Statistics of a set of captures Ys (list of (n, 2)) against their static parts Ss (same shapes)."""
    out = {"e_ratio": [], "lobes": [], "lr": [], "frac": []}
    Ya = np.array(Ys)
    for i, (y, s) in enumerate(zip(Ys, Ss)):
        r = y - s
        # lobes and L / R from the leave-one-out residual: a static-model error recirculates at exactly P and would fill
        # the lobes' centre in the data but not in a simulation (whose static model is exact); the LOO mean removes it in both
        rl = y - np.delete(Ya, i, 0).mean(0); R = rl[:, 1]
        out["e_ratio"].append([10 * np.log10(np.sum(r[K0 + a:K0 + b, 1] ** 2) / np.sum(s[K0 + a:K0 + b, 1] ** 2))
                               for a, b in ((300, 805), (805, 2215), (2215, 4400), (4400, 9000))])
        env = np.sqrt(np.convolve(R ** 2, np.ones(801) / 801, "same")) + 1e-15; rn = R / env
        a = rn[K0 + 400:K0 + 6000]
        out["lobes"].append(np.mean([[np.dot(a, rn[K0 + 400 + Pk + d:K0 + 6000 + Pk + d]) / len(a) for d in range(-12, 13)] for Pk in E.P], 0))
        rr = rl[K0 + 700:K0 + 6000, 1]
        out["lr"].append([np.dot(rr, rl[K0 + 700 + E.LREV + d:K0 + 6000 + E.LREV + d, 0])
                          / np.sqrt(np.dot(rr, rr) * np.dot(rl[K0 + 700 + E.LREV + d:K0 + 6000 + E.LREV + d, 0] ** 2, np.ones(5300))) for d in (-2, -1, 0, 1, 2)])
        fr = []
        for f in OCTS:
            sos = butter(4, [f / 1.414, f * 1.414], btype="band", fs=FS, output="sos")
            yb = sosfiltfilt(sos, y[:, 1])[K0 + 805:K0 + 9000]; sb = sosfiltfilt(sos, s[:, 1])[K0 + 805:K0 + 9000]
            fr.append(np.dot(yb, sb) ** 2 / (np.dot(yb, yb) * np.dot(sb, sb)))
        out["frac"].append(fr)
    return {k: np.mean(v, 0).tolist() for k, v in out.items()} | {k + "_sd": np.std(v, 0).tolist() for k, v in out.items()}


def show(st, tag):
    f = lambda v: " ".join(f"{x:+.2f}" for x in v)
    print(f"{tag}: e_ratio {f(st['e_ratio'])} | lr {f(st['lr'])} | frac {f(st['frac'])}")
    lo = st["lobes"]; print(f"   lobes -12..12: {' '.join(f'{x:+.3f}' for x in lo)}", flush=True)


def data_stats():
    static_params(); rts = [r for r in E.RTS if r >= 1.0]; kern = E.kernels(rts)
    n = K0 + int(2.5 * FS)
    Ys = [E.dev(r)[:n] for r in rts]; Ss = [E.render(X_STATIC, r, kern, n) for r in rts]
    st = stats(Ys, Ss); json.dump(st, open(os.path.join(E.OUT, "rand_stats_data.json"), "w"), indent=1)
    return st


def sim_stats(prm, rts=None):
    static_params(); rts = rts or [r for r in E.RTS if r >= 1.0]; kern = E.kernels([r for r in E.RTS if r >= 1.0])
    Ys, Ss = [], []
    for i, rt in enumerate(rts):
        p = dict(prm); p["seed"] = prm["seed"] * 100 + i
        y = capture(rt, p, kern); Ys.append(y); Ss.append(E.render(X_STATIC, rt, kern, len(y)))
    return stats(Ys, Ss)


def cost(st, ref):
    """Mismatch of simulated statistics st against the data's ref, each term in its own unit."""
    d = lambda k: np.array(st[k]) - np.array(ref[k])
    return float(np.sum((d("e_ratio")[1:] / 1.0) ** 2) + (d("e_ratio")[0] / 3.0) ** 2 + np.sum((d("lobes") / 0.01) ** 2)
                 + np.sum((d("lr") / 0.05) ** 2) + np.sum((d("frac") / 0.03) ** 2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--json", default="{}")
    ap.add_argument("--rts", default="")
    a = ap.parse_args()
    if a.stage == "data":
        show(data_stats(), "data")
    elif a.stage == "render":   # --json '{"rt": 2.5, "lf": 1.4, "hf": 0.7, "ampl": 80, "delay_ms": 0, "lp_khz": 12, "lfo": "time", "chain": "altiverb", "seed": 0, "fs": 48000, "out": "e.wav"}'
        o = json.loads(a.json)
        y = render_ir(o["rt"], o.get("lf"), o.get("hf"), o.get("seed", 0), o.get("secs"), o.get("fs"), o.get("ampl", AMPL_REF),
                      o.get("delay_ms", 0.0), o.get("lp_khz", 12.0), o.get("lfo", "time"), o.get("chain", "altiverb"), o.get("law"))
        import soundfile as sf
        sf.write(o.get("out", "e_render.wav"), (y / np.abs(y).max() * 0.5).astype(np.float32), int(o.get("fs") or round(FS)))
        print(f"wrote {o.get('out', 'e_render.wav')}: {y.shape}")
    elif a.stage == "sim":
        prm = DEF | json.loads(a.json)
        show(sim_stats(prm, [float(v) for v in a.rts.split(",")] if a.rts else None), "sim " + json.dumps(json.loads(a.json)))
