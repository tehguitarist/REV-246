"""The out-of-band check (NEXT_STEPS §5 item 5): every null so far is 100 Hz-10 kHz. Do the existing S, M and L models
also hold at 25-100 Hz and 10-12 kHz, where the LF laws (DC floor, beta), the HF pole limit, the AP damping form and
the M/L output kernel's outer taps act most? Nothing is refitted. Results: docs/ml_structure.md §11.

    .venv/bin/python tools/outofband.py --stage score [--prog SML] [--workers 4] [--only N]   # one file per group
    .venv/bin/python tools/outofband.py --stage report
    .venv/bin/python tools/outofband.py --stage sweep [--workers 8]   # the M/L low-band shortfall through a sweep
    .venv/bin/python tools/outofband.py --stage kernel [--workers 8]  # the M/L output kernel to 12 kHz (needs score)

Model: ml_lines.Model (cluster, diffuser with b, 20 taps, output kernel, the four lines with S's laws and limits; M/L
coupled), for S as well as M and L. Per capture the in-band fit of the ml_lines law stage, unchanged: shift, clock and
two gains per channel on [0, onset + shortest line), 100 Hz-10 kHz. The model is then rendered full band (no
band-pass) with the same clock and gains, and multiplied by the pipeline's own response P(f) (the deconvolution's
band window, flat to 11 kHz with a cosine taper to 12.3 kHz, and the resampling to the device rate; measured by
passing an impulse through both), which is 1 in band.

Per band (zero-phase 3rd-order Butterworth band-passes; 100-200 Hz and 5-10 kHz are in-band controls):
  nulls   per window (cluster and gap, first pass, recirculations 1-3, ml_lines.windows), up to the capture's fade:
          strict (in-band gains) and gain-free (the best scalar per window, so level and shape are separated); the
          same for the repeat-capture model (another capture of the group, clock and gains fitted in band as in
          ml_lines.fit_repeat). Low bands also with the unit's measured sub-20 Hz high-passes on the model
          (analog_results.md §3.5: 3.35 / 0.72 / 0.23 / 0.16 Hz and the Tomkin chain's 1.80 Hz), which the 5.2 ms
          deconvolution kernel cannot remove.
  decay   as ml_lines.band_decay: slope (dB/s) and level of 10 log10(E_capture / E_model) per 39 ms frame from the
          onset to the fade, and RT_capture / RT_model = 1 / (1 - slope RT_law(f) / 60); the same slope for capture
          against repeat capture as the reference; for the low bands also the RT x LF grid.
  sweep   M/L only: the model rendered through ml_variance's simulated exponential sweep, static and with the coupling
          rotation theta(t) of ml_structure §9 (SWEEP_N random phases), and the same decay measure, rotating against
          static: does a sweep of the time-variant unit under-read the low-band RT? -> sweep.json
  kernel  M/L only: the output kernel's response to 12 kHz, measured in the R channel's cluster-only window (4 known
          taps, P(f), stored clocks; every third group), against the stored kernel and symmetric kernels. -> kernel.json
"""
import argparse, glob, json, os, sys

for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.optimize import minimize
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
sys.path.insert(0, os.path.join(HERE, "analog"))
import ml_carryover as cc  # noqa: E402
import ml_lines as mlines  # noqa: E402
import ml_structure as ml  # noqa: E402
from core import _band_window, to_device_grid  # noqa: E402
from sample_hold import highpass  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "outofband")
GDIR = os.path.join(OUT, "groups")
FS, K0 = cc.FS, cc.K0
BANDS = {"25-50": (25, 50), "50-100": (50, 100), "100-200": (100, 200), "5k-10k": (5000, 10000),
         "10k-11k": (10000, 11000), "11k-12k": (11000, 12000)}
LOW = ("25-50", "50-100")
SOS = {b: butter(3, lohi, btype="band", fs=FS, output="sos") for b, lohi in BANDS.items()}
HP_CORNERS = [3.35, 0.72, 0.23, 0.16, 1.80]                 # as q2_injection.HP_CORNERS


def pipeline_response(N):
    """P(f) on an N-point device-rate grid: a unit impulse at 96 kHz through the deconvolution's band window and
    to_device_grid, divided by the ideal device-rate impulse at the same place."""
    fs = 96000.0
    n = 1 << 16
    x = np.zeros(n)
    i0 = 20000
    x[i0] = 1.0
    f96 = np.fft.rfftfreq(n, 1 / fs)
    d = np.fft.irfft(np.fft.rfft(x) * _band_window(f96), n)
    k0 = int(np.floor(i0 / fs * FS))
    h = to_device_grid(d, fs, float(i0), k0, ml.CFG)
    w = 2 * np.pi * np.fft.rfftfreq(N)
    P = np.fft.rfft(h, N) * np.exp(1j * w * k0)
    return P / np.abs(P[int(round(1000 / FS * N))])          # unit gain at 1 kHz (resampling scales by FS / fs)


def captures4(prog, key, n):
    """As ml_carryover.captures (same alignment, estimated on the band-passed captures), returning per AMPL the
    band-passed capture and the full-band capture with the same shifts."""
    g = cc.groups(prog)[key]
    amps = sorted(g)
    hi = amps[-1]
    pad = lambda y: np.vstack([y, np.zeros((max(0, n + 4000 * cc.SC[prog] - len(y)), y.shape[1]))])
    both = lambda f: np.hstack([cc._IR(f, band=True), cc._IR(f)])        # columns: band L, band R, raw L, raw R
    yh = pad(both(g[hi]))
    t0 = cc.tap0(yh[:, :2])
    out = {hi: cc.fshift(yh, K0 - t0)[:n]}
    i0, i1 = K0 + cc.CLUSTER_END[prog] + 50, K0 + cc.first_recirc(prog)
    for a in amps[:-1]:
        y = pad(both(g[a]))
        y = pad(align4(out[hi], y, i0, i1 + 2000 * cc.SC[prog], maxlag=cc.ONSET[prog][1] + 200))
        s = cc.sub_align(out[hi][:, :2], y[:, :2], i0, i1)
        out[a] = cc.fshift(y, s)[:n]
    return {a: (v[:, :2], v[:, 2:]) for a, v in out.items()}, amps


def align4(a, b, s0, s1, maxlag):
    """ml_structure.align for any number of columns: the integer lag from column 1 (band-passed R), applied to all."""
    n = 1 << int(np.ceil(np.log2(len(b) + s1)))
    c = np.fft.irfft(np.fft.rfft(b[:s1 + maxlag, 1], n) * np.conj(np.fft.rfft(np.r_[np.zeros(s0), a[s0:s1, 1]], n)), n)
    lags = np.r_[np.arange(maxlag), np.arange(-maxlag, 0)]
    lag = int(lags[np.argmax(np.r_[c[:maxlag], c[-maxlag:]])])
    return b[lag:] if lag >= 0 else np.vstack([np.zeros((-lag, b.shape[1])), b])


def apply_f(x, H):
    """x (n, c) filtered by the response H on its rfft grid (zero-padded to the grid's length)."""
    N = 2 * (len(H) - 1)
    return np.fft.irfft(np.fft.rfft(x, N, axis=0) * H[:, None], N, axis=0)[:len(x)]


def gain_free(y, m):
    """Null with the best scalar on m."""
    yy, mm, ym = np.sum(y * y), np.sum(m * m), np.sum(y * m)
    return float(10 * np.log10(max(1 - ym * ym / (yy * mm), 1e-12))) if yy > 0 and mm > 0 else None


def decay(y, mdl, t0, t1, sos):
    """ml_lines.band_decay for one band: per channel [slope dB/s, level dB, last frame time] or None."""
    yb, mb = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, mdl, axis=0)
    row = []
    for ch in (0, 1):
        a, b = t0[ch], t1[ch]
        k = (b - a) // mlines.FRAME
        if k < 4:
            row.append(None)
            continue
        ey = (yb[a:a + k * mlines.FRAME, ch] ** 2).reshape(k, mlines.FRAME).mean(1)
        em = (mb[a:a + k * mlines.FRAME, ch] ** 2).reshape(k, mlines.FRAME).mean(1)
        ok = (ey > ey.max() * 1e-7) & (em > em.max() * 1e-7)
        if ok.sum() < 4:
            row.append(None)
            continue
        t = (np.arange(k) + 0.5) * mlines.FRAME / FS
        d = 10 * np.log10(ey[ok] / em[ok])
        wgt = np.sqrt(em[ok] / em[ok].max())
        A = np.stack([np.ones(ok.sum()), t[ok]], 1) * wgt[:, None]
        c = np.linalg.lstsq(A, d * wgt, rcond=None)[0]
        row.append([float(c[1]), float(c[0]), float(t[ok][-1])])
    return row


def _job(a):
    prog, key, path = a
    M = mlines.Model(prog)
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, amps = captures4(prog, key, n)
    n = len(caps[amps[-1]][0])
    W = mlines.windows(prog)
    lf, rt, hf = (float(x) for x in key.split("_"))
    Craw, Rraw = M.parts(n, key, band=False)
    C, R = (sosfiltfilt(ml.BAND, x, axis=0) for x in (Craw, Rraw))      # = M.parts(n, key)
    N = 1 << int(np.ceil(np.log2(n + 4096)))
    P = pipeline_response(N)
    HP = highpass(np.fft.rfftfreq(N, 1 / FS), HP_CORNERS)
    rt_law = {b: float(mlines.rt_of(mlines.line_coefs(M.P, rt, lf, hf), M.P, np.sqrt(lo * hi))[0]) for b, (lo, hi) in BANDS.items()}
    rows = []
    for amp in amps:
        yb, yr = caps[amp]
        end = mlines.usable_end(yb)
        sl, nf = mlines.fit_window(prog, amp > 0)
        _, u, gs, (Cw, _) = mlines.fit_warp(yb, [C, R], sl, nf)                    # the law stage's in-band fit
        Cwr, Rwr = mlines.warp(Craw, *u), mlines.warp(Rraw, *u)
        cl_raw = apply_f(Cwr * np.array([gs[ch][0] for ch in (0, 1)])[None], P)
        mraw = cl_raw + apply_f(Rwr * np.array([gs[ch][1] for ch in (0, 1)])[None], P)
        # repeat capture: another capture of the group, clock and gains fitted in band (ml_lines.fit_repeat)
        other = max((x for x in amps if x != amp), key=lambda x: min(mlines.usable_end(caps[x][0])))
        cl_b = Cw * np.array([gs[ch][0] for ch in (0, 1)])[None]
        ob = caps[other][0]

        def cost(v):
            mdl, _ = mlines.gains(yb[:nf], [cl_b[:nf], mlines.warp(ob[:nf], *v)], sl)
            return sum(np.sum((yb[sl[ch], ch] - mdl[sl[ch], ch]) ** 2) for ch in (0, 1))
        v = minimize(cost, (0.0, 0.0), method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-14, initial_simplex=[(0, 0), (0.2, 0), (0, 10)])).x
        _, gr = mlines.gains(yb[:nf], [cl_b[:nf], mlines.warp(ob[:nf], *v)[:nf]], sl)
        orw = mlines.warp(caps[other][1], *v)
        rep = np.stack([gr[ch][0] * cl_raw[:, ch] + gr[ch][1] * orw[:, ch] for ch in (0, 1)], 1)
        oe = mlines.usable_end(ob)
        variants = {"model": mraw, "repeat": rep, "model_hp": apply_f(mraw, HP)}
        row = dict(key=key, amp=int(amp), clock=list(u), rt_law=rt_law, nulls={}, decay={})
        t0 = [K0 + cc.ONSET[prog][ch] for ch in (0, 1)]
        for b, sos in SOS.items():
            yf = sosfiltfilt(sos, yr, axis=0)
            row["nulls"][b] = {}
            for vn, mv in variants.items():
                if vn == "model_hp" and b not in LOW:
                    continue
                mf = sosfiltfilt(sos, mv, axis=0)
                e = [min(end[ch], oe[ch]) if vn == "repeat" else end[ch] for ch in (0, 1)]
                d = {}
                for w in W[0]:
                    if w == "cluster_gap" and amp == 0:
                        continue
                    cells = []
                    for ch in (0, 1):
                        a_, b_ = K0 + W[ch][w][0], min(K0 + W[ch][w][1], e[ch])
                        if b_ - a_ <= 200:
                            cells.append(None)
                            continue
                        cells.append([mlines.null_db(yf[:, ch], mf[:, ch], a_, b_), gain_free(yf[a_:b_, ch], mf[a_:b_, ch])])
                    d[w] = cells
                row["nulls"][b][vn] = d
            row["decay"][b] = {"model": decay(yr, mraw, t0, end, sos),
                               "repeat": decay(yr, rep, t0, [min(end[ch], oe[ch]) for ch in (0, 1)], sos)}
        rows.append(row)
    json.dump(rows, open(path, "w"))
    return prog, key


def score_stage(progs, workers, only=None):
    tag = "_smoke" if only else ""
    gd = GDIR + tag
    os.makedirs(gd, exist_ok=True)
    jobs = []
    for p in progs:
        ks = sorted(cc.groups(p))
        if only:
            ks = ks[::max(1, len(ks) // only)][:only]
        jobs += [(p, k, os.path.join(gd, f"{p}_{k}.json")) for k in ks]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    print(f"score: {len(jobs)} groups to do", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:     # a fresh process per group returns its memory
        fut = [ex.submit(_job, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            p, k = f.result()
            print(f"  {i + 1}/{len(jobs)} {p} {k}", flush=True)


# ---------------------------------------------------------------- report
def report_stage(tag=""):
    out = {}
    for p in "SML":
        rows = [r for f in sorted(glob.glob(os.path.join(GDIR + tag, f"{p}_*.json"))) for r in json.load(open(f))]
        if not rows:
            continue
        W = list(mlines.windows(p)[0])
        print(f"\n## {p}: {len({r['key'] for r in rows})} groups, {len(rows)} captures; median null L / R, dB "
              f"(strict = in-band gains; free = best scalar per window) [repeat capture, same]")
        res = {}
        for b in BANDS:
            print(f"- **{b} Hz**")
            for w in W:
                cells = {}
                for vn in ("model", "model_hp", "repeat"):
                    x = [r["nulls"][b].get(vn, {}).get(w) for r in rows]
                    x = [c for c in x if c and c[0] is not None and c[1] is not None and None not in c[0] + c[1]]
                    if not x:
                        continue
                    a = np.array(x)                                  # (n, ch, strict/free)
                    cells[vn] = np.median(a, 0).tolist()
                if not cells:
                    continue
                res.setdefault(b, {})[w] = cells
                m, q = cells["model"], cells["repeat"]
                s = (f"  {w:12s}: strict {m[0][0]:.1f} / {m[1][0]:.1f}, free {m[0][1]:.1f} / {m[1][1]:.1f} "
                     f"[repeat strict {q[0][0]:.1f} / {q[1][0]:.1f}, free {q[0][1]:.1f} / {q[1][1]:.1f}]")
                if "model_hp" in cells:
                    h = cells["model_hp"]
                    s += f"; with the high-passes strict {h[0][0]:.1f} / {h[1][0]:.1f}, free {h[0][1]:.1f} / {h[1][1]:.1f}"
                print(s)
            # decay: RT ratio capture / model, and the repeat-pair slope as the reference
            rat, sl_m, sl_r, lev = [], [], [], []
            for r in rows:
                dm, dr = r["decay"][b]["model"], r["decay"][b]["repeat"]
                for ch in (0, 1):
                    if dm[ch] is None:
                        continue
                    T = r["rt_law"][b]
                    rat.append(1 / (1 - dm[ch][0] * T / 60))
                    sl_m.append(dm[ch][0])
                    lev.append(dm[ch][1])
                    if dr[ch] is not None:
                        sl_r.append(dr[ch][0])
            if rat:
                rat = np.array(rat)
                print(f"  decay: RT capture / model median {np.median(rat):.3f} (IQR {np.percentile(rat, 25):.3f}-{np.percentile(rat, 75):.3f}); "
                      f"slope capture vs model {np.median(sl_m):+.2f} dB/s, vs repeat {np.median(sl_r):+.2f} dB/s "
                      f"(|slope| median {np.median(np.abs(sl_m)):.2f} / {np.median(np.abs(sl_r)):.2f}); level at onset {np.median(lev):+.2f} dB")
                res.setdefault(b, {})["decay"] = dict(rt_ratio=float(np.median(rat)), iqr=[float(np.percentile(rat, 25)), float(np.percentile(rat, 75))],
                                                      slope_model=float(np.median(sl_m)), slope_repeat=float(np.median(sl_r)),
                                                      abs_slope=[float(np.median(np.abs(sl_m))), float(np.median(np.abs(sl_r)))],
                                                      level=float(np.median(lev)))
                # by the control that acts in this band
                c, ci = ("lf", 0) if b in LOW + ("100-200",) else ("hf", 2) if b != "5k-10k" else ("rt", 1)
                byc = {}
                for r in rows:
                    dm = r["decay"][b]["model"]
                    for ch in (0, 1):
                        if dm[ch] is not None:
                            byc.setdefault(float(r["key"].split("_")[ci]), []).append(1 / (1 - dm[ch][0] * r["rt_law"][b] / 60))
                print(f"  RT capture / model by {c}: " + "; ".join(f"{k}: {np.median(v):.3f}" for k, v in sorted(byc.items())))
                res[b]["decay"]["by_" + c] = {str(k): float(np.median(v)) for k, v in byc.items()}
                if b in LOW + ("100-200",):                          # RT x LF grid (median over HF and channels)
                    grid = {}
                    for r in rows:
                        lf_, rt_ = (float(x) for x in r["key"].split("_")[:2])
                        for d_ in r["decay"][b]["model"]:
                            if d_ is not None:
                                grid.setdefault((rt_, lf_), []).append(1 / (1 - d_[0] * r["rt_law"][b] / 60))
                    lfs = sorted({k[1] for k in grid})
                    print(f"  RT capture / model, rows RT, columns LF {lfs}:")
                    for rt_ in sorted({k[0] for k in grid}):
                        print(f"    {rt_:>5}: " + " ".join(f"{np.median(grid[(rt_, l)]):6.3f}" if (rt_, l) in grid else "     -" for l in lfs))
                    res[b]["decay"]["grid"] = {f"{k[0]}_{k[1]}": float(np.median(v)) for k, v in grid.items()}
        out[p] = res
    json.dump(out, open(os.path.join(OUT, f"report{tag}.json"), "w"), indent=1)


# ---------------------------------------------------------------- sweep: the low-band shortfall of a time-variant unit
SWEEP_KEYS = {"M": ["2.0_5.0_0.9", "0.1_5.0_0.9"], "L": ["2.0_20.0_1.0", "2.0_8.0_1.0", "0.1_20.0_1.0"]}
SWEEP_TAIL, SWEEP_N = 12.0, 3


def _sweep_job(a):
    """One render of the M/L model through ml_variance's simulated sweep (its final band widened to 15 Hz-12 kHz):
    static (ph None) or with theta = THETA_A sin(2 pi THETA_FM t + ph)."""
    import ml_variance as mv
    ml.BAND = butter(4, [15, 12000], btype="band", fs=FS, output="sos")     # module global read by sweep_measure
    prog, key, ph = a
    fn = (lambda t: 0 * t) if ph is None else (lambda t: mv.THETA_A * np.sin(2 * np.pi * mv.THETA_FM * t + ph))
    return prog, key, ph, mv.sweep_measure(mlines.Model(prog), key, fn, tail_s=SWEEP_TAIL)


def sweep_stage(workers):
    rng = np.random.default_rng(7)
    jobs = [(p, k, ph) for p, ks in SWEEP_KEYS.items() for k in ks for ph in [None] + list(rng.uniform(0, 2 * np.pi, SWEEP_N))]
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        res = list(ex.map(_sweep_job, jobs))
    out = {}
    print(f"\n## sweep: RT rotating / static, the M/L model through a simulated sweep ({SWEEP_N} phases, {SWEEP_TAIL} s tail)")
    for p, ks in SWEEP_KEYS.items():
        M = mlines.Model(p)
        for k in ks:
            lf, rt, hf = (float(x) for x in k.split("_"))
            st = [y for q, kk, ph, y in res if q == p and kk == k and ph is None][0]
            rots = [y for q, kk, ph, y in res if q == p and kk == k and ph is not None]
            t0 = [K0 + cc.ONSET[p][0], K0 + cc.ONSET[p][1]]
            t1 = [len(st) - int(0.5 * FS)] * 2
            cells = {}
            for b, sos in SOS.items():
                T = float(mlines.rt_of(mlines.line_coefs(M.P, rt, lf, hf), M.P, np.sqrt(np.prod(BANDS[b])))[0])
                r = [1 / (1 - c[0] * T / 60) for y in rots for c in decay(y, st, t0, t1, sos) if c]
                cells[b] = dict(ratio=float(np.median(r)), rt_law=T)
            out[f"{p}_{k}"] = cells
            print(f"- {p} {k}: " + "; ".join(f"{b} {c['ratio']:.3f} (RT_law {c['rt_law']:.1f} s)" for b, c in cells.items()))
    json.dump(out, open(os.path.join(OUT, "sweep.json"), "w"), indent=1)


# ---------------------------------------------------------------- kernel: the M/L output kernel up to 12 kHz
def _kernel_job(a):
    """Cross and auto spectra of capture against the R channel's 4 cluster taps (no kernel, pipeline P(f), stored clock)
    over the cluster-only window [-30, onset - 3), per capture with a cluster."""
    prog, key = a
    rows = {r["amp"]: r for r in json.load(open(os.path.join(GDIR, f"{prog}_{key}.json")))}
    n = K0 + cc.ONSET[prog][1] + 3000
    caps, amps = captures4(prog, key, n + 4000)
    N = 1 << 14
    w = 2 * np.pi * np.fft.rfftfreq(N)
    C0 = sum(s_ * np.exp(-1j * w * (p_ + K0)) for p_, s_ in cc.cluster_taps(prog)[1])
    c0 = np.fft.irfft(C0 * pipeline_response(N), N)[:n + 2000]
    sl = slice(K0 - 30, K0 + cc.ONSET[prog][1] - 3)
    win = np.hanning(sl.stop - sl.start)
    out = []
    for amp in amps:
        if amp > 0:
            cw = mlines.warp(c0[:, None], *rows[amp]["clock"])[:, 0]
            Y, X = np.fft.rfft(caps[amp][1][sl, 1] * win, N), np.fft.rfft(cw[sl] * win, N)
            out.append((Y * np.conj(X), np.abs(X) ** 2))
    return out


def kernel_stage(workers):
    """Measured kernel response (every third group), gain matched to the stored kernel at 2-6 kHz; against the stored
    11-tap kernel and symmetric 3- and 5-tap kernels fitted to the measured magnitude in band (0.2-10 kHz)."""
    fp = json.load(open(os.path.join(cc.OUT, "firstpass.json")))
    out = {}
    for prog in "ML":
        with ProcessPoolExecutor(workers) as ex:
            res = [x for rr in ex.map(_kernel_job, [(prog, k) for k in sorted(cc.groups(prog))[::3]]) for x in rr]
        N = 1 << 14
        f = np.fft.rfftfreq(N, 1 / FS)
        w = 2 * np.pi * f / FS
        K = sum(r[0] for r in res) / sum(r[1] for r in res)
        Hk = sum(hk * np.exp(-1j * w * (k - cc.KH)) for k, hk in enumerate(fp[prog]["kernel"]))
        i, j = (f > 2000) & (f < 6000), (f > 200) & (f < 10000)
        K = K / np.mean(np.abs(K[i]) / np.abs(Hk[i]))
        sym = {}
        for nt in (1, 2):
            A = np.stack([np.cos(m * w[j]) - 1 for m in range(1, nt + 1)], 1)
            c = np.linalg.lstsq(A, np.abs(K[j]) - 1, rcond=None)[0]
            sym[2 * nt + 1] = 1 + sum(ci * (np.cos(m * w) - 1) for m, ci in enumerate(c, 1))
        print(f"\n## {prog}: output kernel from the R cluster window ({len(res)} captures): |K| dB and phase, "
              f"against the stored 11-tap / symmetric 3-tap / 5-tap (fitted in band) |H| dB")
        out[prog] = {}
        for fq in (1000, 4000, 8000, 9000, 10000, 10500, 11000, 11500, 12000):
            k = int(np.argmin(abs(f - fq)))
            v = dict(mag_db=20 * np.log10(abs(K[k])), phase_deg=float(np.degrees(np.angle(K[k]))), stored_db=20 * np.log10(abs(Hk[k])),
                     sym3_db=20 * np.log10(abs(sym[3][k])), sym5_db=20 * np.log10(abs(sym[5][k])))
            out[prog][fq] = {a: float(b) for a, b in v.items()}
            print(f"  {fq:>6} Hz: {v['mag_db']:+.2f} dB, {v['phase_deg']:+.1f} deg | {v['stored_db']:+.2f} / {v['sym3_db']:+.2f} / {v['sym5_db']:+.2f}")
    json.dump(out, open(os.path.join(OUT, "kernel.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="report", choices=["score", "report", "sweep", "kernel"])
    ap.add_argument("--prog", default="SML")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", type=int)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.stage == "sweep":
        sweep_stage(a.workers)
        return
    if a.stage == "kernel":
        kernel_stage(a.workers)
        return
    if a.stage == "score":
        score_stage(a.prog, a.workers, a.only)
    report_stage("_smoke" if a.only else "")


if __name__ == "__main__":
    main()
