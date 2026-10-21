"""Q2: offline tests of the period-practice candidates (docs/period_dsp_practice.md §4) on the frozen flat injection.
Results: docs/model_s_results.md §16, docs/period_dsp_practice.md §7.

    .venv/bin/python tools/q2_candidates.py [--stage k1,k2,structure] [--src branches_flat_0_ext1700] [--ext 1700]
                                             [--loops 0,1,2,3] [--kmax 64] [--workers 8] [--mstep 1] [--out analysis/q2]

Input: the shared single-branch injection saved by `q2_injection.py --stage branches --alphas 0 --ext 1700 --subset flat`
(layout as q1.Joint.split: 8 cluster gains, the 2 x 1,400 first-pass FIR F, then per channel the four x_i of P_i + ext
samples starting 230 samples after cluster tap 0). Each x_i sits in front of `post` (comb460 . AP20 . AP72 . AP230) and
the loop 1 / (1 - G_i z^-P_i); F in front of `post` only. Reads only that file; nothing is fitted to the captures, so
every number here is in-sample: the held-out test is the bench.

Every score is in band: signals are band-limited with the bench's fit band B = |LFBAND|^2 (zero phase) and fitted by
greedy OMP with band-limited atoms (integer taps, optionally adjacent pairs (1 - f) d[t] + f d[t + 1], f on a 1/8 grid);
the score is the in-band energy explained at K atoms.

Stages:
  k1         K1 (tapped read-out). sparsity: OMP curves per loop and channel, with controls (a) 4 known taps (one past
             P_i, two fractional, equal magnitude) + in-band noise at 1 %: the taps must be found; (b) phase-randomised
             x_i (same magnitude spectrum): must stay dense. comb: the taps come in pairs t, t + 460 with ratio about
             -0.44, i.e. x_i carries (1 - g z^-460), cancelling the bench's comb460: divide it out, scan g by the energy 5
             taps explain, list the taps. tied: x_ic = a . sum_k s_k d[t_k + c . u_k] . (1 - g z^-460) with one magnitude
             a, L taps t_k and signs s_k, R taps = L taps + u_k, u_k in {26, 332} (matched): energy explained.
  k2         K2 / K5 screen (a recursive stage per loop, or a different point on the diffuser chain): for every delay
             M = 2..1700 and gain g on a grid, four re-expressions of x_i: t = x (1 - g z^-M) (an extra comb stage),
             t = x / (1 - g z^-M) (fed before a comb), t = x . AP(M, g) (fed before an all-pass), t = x . AP(M, g)
             time-reversed (an extra all-pass stage). Score: a fast 5-tap proxy (the 5 largest separated peaks of the
             normalised band-limited correlation), then exact OMP for the best few. Run on the raw x_i (does the scan find
             the comb460 factor blind?) and on the comb-divided x_i (is any further stage needed?).
  structure  K3 (tapped injection x tapped read-out): pairwise differences of each loop's 5 taps, coincidences within
             and across loops. K6 (one per-channel distributor for all loops): R_i = (c0 + c26 z^-26 + c332 z^-332) L_i
             with c common to all loops, against the per-tap 26/332 labels. First-pass FIR F: sparsity raw and comb-
             divided. Residual of the tied model: where it sits (near the taps or not) and its spectrum by octave.
  k5         K5 (nested all-pass) and anything else recursive, on what the tied model leaves: normalised in-band
             autocorrelation of the tied residual at lags 10..1700, against Gaussian noise with the residual's smoothed
             spectrum (a nested all-pass leaves peaks at M1, M2, M1 +- M2; lags that are loop lengths or their
             differences point at the free fit trading content between loops instead).
  settings   K7 and K4 on the per-setting free injections (`q2_injection.py --stage free`, analysis/q2/v/, window
             P_i + 300, all 128 settings): K7, are the flat taps there in every setting (5-atom OMP on the comb-divided
             excitation: positions within +-1 of the flat taps; energy the flat tap set explains with free gains) and do
             their gains stay put; K4, does any tap carry the loop filter once more than the others: per-tap gains fitted
             in a low band (40-250 Hz) and a high band (4-10 kHz) against the full band, spread across the 5 taps in dB,
             against the split K4 predicts there (|H_i| of the setting's loop filter in that band re 1 kHz). The decisive
             form: per-tap band-gain profiles at each LF / HF step minus the profile at LF 1.0 / HF 1.0, averaged over
             settings with RT >= 1.4 s (short RT is poorly determined); K4 needs the same taps lower by |H_i| at every
             strong step, in L and R alike.
Output: k1_sparsity.json / .log, k2_stages.json / .log, structure.json / .log, k5_resid.json / .log,
settings.json / .log in --out.
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import lfilter, sosfreqz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import q1_generalise as q1  # noqa: E402

P = q1.P
N = 16384                                      # circular buffer (x_i <= 3,300 samples; room for all-pass tails in k2)
PAD = 6144                                     # x_i sits at [PAD, PAD + L)
FRACS = np.arange(8) / 8                       # adjacent-pair atoms
KS = [1, 2, 4, 8, 16, 32, 64]
NTAP = 5
GRID = [-0.8, -0.7, -0.6, -0.45, -0.3, 0.3, 0.45, 0.6, 0.7, 0.8]
FAMILIES = ["comb_extra", "comb_before", "ap_before", "ap_extra"]
B = np.abs(sosfreqz(q1.LFBAND, worN=2 * np.pi * np.fft.rfftfreq(N))[1]) ** 2
OUT = os.path.join(ROOT, "analysis", "q2")


def band(x):
    b = np.zeros(N); b[PAD:PAD + len(x)] = x
    return np.fft.irfft(np.fft.rfft(b) * B, N)


def kernels(fracs):
    ks = []
    for f in fracs:
        d = np.zeros(N); d[0] = 1 - f; d[1] = f
        ks.append(np.fft.irfft(np.fft.rfft(d) * B, N))
    K = np.array(ks)
    return K, np.fft.rfft(K, axis=1), np.linalg.norm(K, axis=1)


K0 = kernels(np.array([0.0]))


def omp(y, L, fracs, kmax, ks=KS, lo=-8, hi=8):
    """Greedy OMP on band-limited y; atoms at t in [lo, L + hi) relative to PAD. Returns explained fraction at each
    K in ks, the selected (frac, t) and the final gains."""
    Kt, KF, kn = K0 if len(fracs) == 1 and fracs[0] == 0 else kernels(fracs)
    lo, hi = PAD + lo, PAD + L + hi
    e0, r, sel, A, out = y @ y, y.copy(), [], [], {}
    for k in range(1, kmax + 1):
        c = np.fft.irfft(np.fft.rfft(r)[None] * np.conj(KF), N)[:, lo:hi] / kn[:, None]
        j, t = np.unravel_index(np.argmax(np.abs(c)), c.shape)
        sel.append((float(fracs[j]), int(t + lo - PAD)))
        A.append(np.roll(Kt[j], t + lo))
        M = np.array(A).T
        g, *_ = np.linalg.lstsq(M, y, rcond=None)
        r = y - M @ g
        if k in ks:
            out[k] = float(1 - (r @ r) / e0)
    return out, sel, g


def load(path, ext):
    v = np.load(path)
    o, xs = 8 + 2 * q1.LFIR, {}
    assert v.size == o + 2 * int(np.sum(P + ext)), "layout: single-branch injection with this --ext expected"
    F = v[8:o].reshape(2, q1.LFIR)
    for c in (0, 1):
        for i in range(4):
            xs[c, i] = v[o:o + P[i] + ext]; o += P[i] + ext
    return xs, F


def rec(y, M, g):
    """y / (1 - g z^-M), causal, on the whole buffer."""
    if M < 64:
        a = np.zeros(M + 1); a[0] = 1; a[M] = -g
        return lfilter([1.0], a, y)
    out = y.copy()
    for k in range(M, len(y), M):
        out[k:k + M] += g * out[k - M:k][:len(out[k:k + M])]
    return out


def allpass(y, M, g):
    """(-g + z^-M) / (1 - g z^-M), causal."""
    v = rec(y, M, g)
    out = -g * v
    out[M:] += v[:-M]
    return out


def decomb(y, g):
    return rec(y, 460, g)


def restage(y, fam, M, g):
    if fam == "comb_extra":
        out = y.copy(); out[M:] -= g * y[:-M]; return out
    if fam == "comb_before":
        return rec(y, M, g)
    if fam == "ap_before":
        return allpass(y, M, g)
    return allpass(y[::-1], M, g)[::-1]                      # ap_extra: the inverse all-pass is the time-reversed one


def proxy5(t, L):
    """Fast 5-tap score: energy of the 5 largest separated peaks of the normalised correlation with the integer kernel,
    relative to the whole buffer's energy, so energy a stage pushes out of the window counts as unexplained (atoms more
    than a few samples apart are nearly orthogonal)."""
    c = np.fft.irfft(np.fft.rfft(t) * np.conj(K0[1][0]), N)[PAD - 8:PAD + L + 8] / K0[2][0]
    a = np.abs(c)
    s = 0.0
    for _ in range(NTAP):
        j = int(np.argmax(a)); s += a[j] ** 2; a[max(j - 4, 0):j + 5] = 0
    return s / (t @ t)


def _k2_job(args):
    key, y, L, fam, Ms, gs = args
    best = []
    for M in Ms:
        for g in gs:
            best.append((proxy5(restage(y, fam, M, g), L), M, g))
    best.sort(reverse=True)
    return key, fam, best[:10]


def stage_k1(a, xs, ys, loops, say, rng):
    res = dict(src=a.src, ext=a.ext, band="|LFBAND|^2 (20 Hz-11 kHz, zero phase)", sparsity={}, comb={}, tied={})
    fmt = lambda d: " ".join(f"{k}:{100 * v:5.1f}" for k, v in d.items())
    say(f"[sparsity] {a.src}: in-band energy explained (%) at K atoms")
    for i in loops:
        for c in (0, 1):
            x, y = xs[c, i], ys[c, i]; L, e = len(x), ys[c, i] @ ys[c, i]
            row = dict(overhang=float(np.sum(y[PAD + P[i]:PAD + L] ** 2) / e))
            row["int"], _, _ = omp(y, L, np.array([0.0]), a.kmax)
            row["pair"], _, _ = omp(y, L, FRACS, a.kmax)
            X = np.fft.rfft(x); ph = np.exp(2j * np.pi * rng.random(len(X))); ph[0] = 1
            row["ctrl_phase"], _, _ = omp(band(np.fft.irfft(np.abs(X) * ph, L)), L, FRACS, a.kmax)
            ts = np.r_[rng.integers(0, P[i], 3), rng.integers(P[i], L - 2)]
            fs = np.r_[0.0, 0.0, 0.375, 0.625]
            d = np.zeros(L)
            for t, f in zip(ts, fs):
                s = rng.choice([-1, 1]); d[t] += s * (1 - f); d[t + 1] += s * f
            yk = band(d); yk *= np.sqrt(e / (yk @ yk))
            nz = band(rng.standard_normal(L)); nz *= np.sqrt(0.01 * e / (nz @ nz))
            row["ctrl_known"], selk, _ = omp(yk + nz, L, FRACS, 16, [4, 8, 16])
            row["ctrl_known_true"] = sorted((int(t), float(f)) for t, f in zip(ts, fs))
            row["ctrl_known_found"] = sorted((t, f) for f, t in selk[:4])
            res["sparsity"][f"{P[i]}{'LR'[c]}"] = row
            say(f"  {P[i]} {'LR'[c]}  overhang past P_i {100 * row['overhang']:4.1f} %")
            say(f"     integer      {fmt(row['int'])}")
            say(f"     + pairs      {fmt(row['pair'])}")
            say(f"     phase ctrl   {fmt(row['ctrl_phase'])}")
            say(f"     known ctrl   {fmt(row['ctrl_known'])}  taps found: {row['ctrl_known_true'] == row['ctrl_known_found']}")

    say("\n[comb] x / (1 - g z^-460): energy explained by 5 integer taps, mean (min) over the excitations")
    scan = {}
    for g in np.round(np.arange(0.40, 0.4851, 0.005), 3):
        ex = [omp(decomb(ys[c, i], g), len(xs[c, i]), np.array([0.0]), NTAP, [NTAP])[0][NTAP] for i in loops for c in (0, 1)]
        scan[float(g)] = (float(np.mean(ex)), float(np.min(ex)))
        say(f"  g {g:.3f}: {100 * scan[g][0]:.2f} ({100 * scan[g][1]:.2f})")
    g460 = max(scan, key=lambda k: scan[k][0])
    say(f"  best g460 = {g460:.3f}")
    taps = {}
    for i in loops:
        for c in (0, 1):
            out, sel, cf = omp(decomb(ys[c, i], g460), len(xs[c, i]), np.array([0.0]), 8, [4, 5, 8])
            o = np.argsort(-np.abs(cf))
            taps[c, i] = sorted((sel[j][1], float(cf[j])) for j in o[:NTAP])
            res["comb"][f"{P[i]}{'LR'[c]}"] = dict(explained=out, taps=taps[c, i], next_gain=float(np.abs(cf[o[NTAP]])))
            say(f"  {P[i]} {'LR'[c]}  K=4/5/8: {100 * out[4]:.2f} / {100 * out[5]:.2f} / {100 * out[8]:.2f} %   taps (t after 230 : gain) "
                + "  ".join(f"{t}:{v:+.4f}" for t, v in taps[c, i]) + f"   next |gain| {np.abs(cf[o[NTAP]]):.4f}")
    res["comb"]["g460"] = g460
    res["comb"]["scan"] = scan

    say("\n[tied] R taps = L taps + 26 or 332, one magnitude for all, (1 - g z^-460) factor")
    shift = {}
    for i in loops:
        for t, v in taps[1, i]:
            hit = [(tl, s, np.sign(vl) == np.sign(v)) for tl, vl in taps[0, i] for s in (26, 332) if tl + s == t]
            if not hit:
                say(f"  {P[i]}: R tap {t} has no L tap at -26 / -332: tied model not applicable"); break
            shift[i, hit[0][0]] = hit[0][1]
            if not hit[0][2]:
                say(f"  {P[i]}: R tap {t} has the opposite sign of L tap {hit[0][0]}")
    A = {(c, i): tied_atoms(xs, taps, shift, g460, i, c) for i in loops for c in (0, 1)}
    amp = sum(A[q] @ ys[q] for q in A) / sum(A[q] @ A[q] for q in A)
    tr = ty = 0.0
    for i in loops:
        res["tied"][str(P[i])] = dict(L_taps=[(int(t + q1.S0), int(np.sign(v))) for t, v in taps[0, i]],
                                     R_shift=[int(shift.get((i, t), -1)) for t, _ in taps[0, i]])
        for c in (0, 1):
            r = ys[c, i] - amp * A[c, i]; yy = ys[c, i] @ ys[c, i]
            tr += r @ r; ty += yy
            res["tied"][str(P[i])]["explained_" + "LR"[c]] = float(1 - (r @ r) / yy)
            say(f"  {P[i]} {'LR'[c]}: explained {100 * (1 - (r @ r) / yy):.2f} %  (residual {10 * np.log10((r @ r) / yy):.1f} dB)")
        say(f"     L taps (lag after cluster tap 0, sign): {res['tied'][str(P[i])]['L_taps']}   R shifts: {res['tied'][str(P[i])]['R_shift']}")
    res["tied"].update(a=float(amp), explained_all=float(1 - tr / ty),
                       params="20 positions, 20 signs, 20 shift labels (26/332), magnitude a, g460")
    say(f"  all: explained {100 * (1 - tr / ty):.2f} %  (residual {10 * np.log10(tr / ty):.1f} dB), a = {amp:.4f}, g460 = {g460:.3f}")
    json.dump(res, open(os.path.join(OUT, "k1_sparsity.json"), "w"), indent=1)
    return taps, shift, g460, amp, A


def tied_atoms(xs, taps, shift, g460, i, c):
    d = np.zeros(len(xs[c, i]) + 470)
    for t, v in taps[0, i]:
        tt = t + (shift.get((i, t), 0) if c else 0)
        d[tt] += np.sign(v); d[tt + 460] -= g460 * np.sign(v)
    return band(d[:len(xs[c, i])])


def stage_k2(a, xs, ys, loops, say, g460):
    from concurrent.futures import ProcessPoolExecutor
    Ms = list(range(2, 1701, a.mstep))
    res = {}
    for src in ("raw", "comb_divided"):
        say(f"[k2] {src} x_i: best re-expressions by the fast 5-tap score (identity in brackets), then exact OMP")
        jobs = []
        for i in loops:
            for c in (0, 1):
                y = ys[c, i] if src == "raw" else decomb(ys[c, i], g460)
                jobs += [((c, i), y, len(xs[c, i]), fam, Ms, GRID) for fam in FAMILIES]
        with ProcessPoolExecutor(a.workers) as ex:
            out = list(ex.map(_k2_job, jobs))
        res[src] = {}
        for i in loops:
            for c in (0, 1):
                y = ys[c, i] if src == "raw" else decomb(ys[c, i], g460)
                L = len(xs[c, i])
                ident = proxy5(y, L)
                exact_id = omp(y, L, np.array([0.0]), NTAP, [NTAP])[0][NTAP]
                rows = sorted(((s, M, g, fam) for (q, fam, best) in out if q == (c, i) for s, M, g in best), reverse=True)[:5]
                ex_rows = [(float(omp(restage(y, fam, M, g), L, np.array([0.0]), NTAP, [NTAP])[0][NTAP]), fam, M, g) for s, M, g, fam in rows]
                res[src][f"{P[i]}{'LR'[c]}"] = dict(identity_proxy=ident, identity_exact=exact_id,
                                                    best=[dict(proxy=s, exact=e, family=f, M=M, g=g) for (s, M, g, f), (e, *_) in zip(rows, ex_rows)])
                say(f"  {P[i]} {'LR'[c]}  identity {100 * exact_id:.2f} %  [proxy {100 * ident:.2f}]   best: "
                    + "  ".join(f"{f} M={M} g={g:+.2f}: {100 * e:.2f} %" for e, f, M, g in ex_rows[:3]))
    json.dump(res, open(os.path.join(OUT, "k2_stages.json"), "w"), indent=1)


def stage_structure(a, xs, ys, F, loops, say, taps, shift, g460, amp, A):
    res = {}
    say("[k3] pairwise differences of each loop's 5 L taps")
    diffs = {}
    for i in loops:
        ts = [t for t, _ in taps[0, i]]
        diffs[i] = sorted(ts[k] - ts[j] for j in range(5) for k in range(j + 1, 5))
        say(f"  {P[i]}: taps {ts}  differences {diffs[i]}")
    alld = [d for i in loops for d in diffs[i]]
    within = sum(len(v) - len(set(v)) for v in diffs.values())
    across = [(d, [int(P[i]) for i in loops if d in diffs[i]]) for d in sorted(set(alld)) if sum(d in diffs[i] for i in loops) > 1]
    near = [(d, e) for d in sorted(set(alld)) for e in sorted(set(alld)) if 0 < e - d <= 1]
    rng = np.random.default_rng(7)
    sim = []
    for _ in range(2000):                                     # chance level: 5 random taps per loop within the same spans
        ds = []
        for i in loops:
            ts = np.sort(rng.choice(np.arange(min(t for t, _ in taps[0, i]), max(t for t, _ in taps[0, i]) + 1), 5, replace=False))
            ds.append(set(ts[k] - ts[j] for j in range(5) for k in range(j + 1, 5)))
        sim.append(sum(1 for d in set().union(*ds) if sum(d in s for s in ds) > 1))
    say(f"  repeated differences within a loop: {within}; shared across loops: {across} (chance: mean {np.mean(sim):.2f}, "
        f"95th pct {np.percentile(sim, 95):.0f}); differences 1 apart: {near}")
    res["k3"] = dict(differences={int(P[i]): diffs[i] for i in loops}, within=within, across=across,
                     chance_mean=float(np.mean(sim)), chance_p95=float(np.percentile(sim, 95)))

    say("\n[k6] one distributor for all loops: R_i = (c0 + c26 z^-26 + c332 z^-332) L_i, c common, vs the per-tap labels")
    def shifted(y, s):
        return np.roll(y, s)
    cols, tgt = [], []
    for i in loops:
        yl = amp * A[0, i]
        cols.append(np.array([shifted(yl, s) for s in (0, 26, 332)]))
        tgt.append(ys[1, i])
    M = np.concatenate([c.T for c in cols]); T = np.concatenate(tgt)
    cc, *_ = np.linalg.lstsq(M, T, rcond=None)
    res["k6"] = dict(c=[float(v) for v in cc], per_loop={})
    say(f"  common c0 / c26 / c332 = {cc[0]:+.3f} / {cc[1]:+.3f} / {cc[2]:+.3f}")
    for k, i in enumerate(loops):
        r = tgt[k] - cols[k].T @ cc; free, *_ = np.linalg.lstsq(cols[k].T, tgt[k], rcond=None)
        rf = tgt[k] - cols[k].T @ free
        rt = ys[1, i] - amp * A[1, i]
        e = tgt[k] @ tgt[k]
        res["k6"]["per_loop"][int(P[i])] = dict(common=float(1 - r @ r / e), per_loop_c=float(1 - rf @ rf / e), per_tap_labels=float(1 - rt @ rt / e),
                                                c_loop=[float(v) for v in free])
        say(f"  {P[i]} R explained: common distributor {100 * (1 - r @ r / e):.1f} %, per-loop c {100 * (1 - rf @ rf / e):.1f} % "
            f"(c = {free[0]:+.2f} / {free[1]:+.2f} / {free[2]:+.2f}), per-tap 26/332 labels {100 * (1 - rt @ rt / e):.1f} %")

    say("\n[fir] first-pass FIR F (1,400 taps from 230 after tap 0, through post only): in-band energy explained at K")
    res["fir"] = {}
    for c in (0, 1):
        y = band(F[c]); L = F.shape[1]
        raw = omp(y, L, np.array([0.0]), 32, [4, 8, 16, 32])[0]
        dec = omp(decomb(y, g460), L, np.array([0.0]), 32, [4, 8, 16, 32])[0]
        X = np.fft.rfft(F[c]); ph = np.exp(2j * np.pi * np.random.default_rng(3 + c).random(len(X))); ph[0] = 1
        ctl = omp(band(np.fft.irfft(np.abs(X) * ph, L)), L, np.array([0.0]), 32, [4, 8, 16, 32])[0]
        res["fir"]["LR"[c]] = dict(raw=raw, comb_divided=dec, phase_ctrl=ctl, energy_rel_x=float((y @ y) / np.mean([ys[c, i] @ ys[c, i] for i in loops])))
        f = lambda d: " ".join(f"{k}:{100 * v:5.1f}" for k, v in d.items())
        say(f"  {'LR'[c]}  energy re a mean x_i {res['fir']['LR'[c]]['energy_rel_x']:.3f}   raw {f(raw)}   comb-divided {f(dec)}   phase ctrl {f(ctl)}")

    say("\n[resid] tied-model residual: share within +-8 samples of a tap or its +460 partner, and by octave (re the excitation's own octave)")
    res["resid"] = {}
    edges = [20, 177, 354, 707, 1414, 2828, 5657, 11314]
    w = np.fft.rfftfreq(N, 1 / q1.FS)
    for i in loops:
        for c in (0, 1):
            y = ys[c, i]; r = y - amp * A[c, i]
            mask = np.zeros(N, bool)
            for t, _ in taps[0, i]:
                tt = PAD + t + (shift.get((i, t), 0) if c else 0)
                for u in (tt, tt + 460):
                    mask[u - 8:u + 9] = True
            near = float((r[mask] @ r[mask]) / (r @ r))
            Rf, Yf = np.abs(np.fft.rfft(r)) ** 2, np.abs(np.fft.rfft(y)) ** 2
            octs = [float(10 * np.log10(Rf[(w >= lo) & (w < hi)].sum() / Yf[(w >= lo) & (w < hi)].sum())) for lo, hi in zip(edges[:-1], edges[1:])]
            res["resid"][f"{P[i]}{'LR'[c]}"] = dict(near_taps=near, mask_share=float(mask.sum() / len(xs[c, i])), octaves_db=octs)
            say(f"  {P[i]} {'LR'[c]}  near taps {100 * near:4.1f} % of the residual (window {100 * mask.sum() / len(xs[c, i]):.1f} % of samples)"
                "   octaves " + " ".join(f"{v:6.1f}" for v in octs))
    say(f"  octave bands (Hz): {list(zip(edges[:-1], edges[1:]))}")
    json.dump(res, open(os.path.join(OUT, "structure.json"), "w"), indent=1)


def stage_k5(xs, ys, loops, say, taps, shift, g460, amp, A):
    say("[k5] tied residual: largest normalised in-band autocorrelation at lags 10..1700 (control: Gaussian noise with its smoothed spectrum, 95th pct of the max)")
    res = {}
    rng = np.random.default_rng(11)
    for i in loops:
        for c in (0, 1):
            r = (ys[c, i] - amp * A[c, i])
            def acf(z):
                Z = np.abs(np.fft.rfft(z)) ** 2
                ac = np.fft.irfft(Z, N)
                return ac[10:1701] / ac[0]
            ac = acf(r)
            top = np.argsort(-np.abs(ac))[:6]
            env = np.convolve(np.abs(np.fft.rfft(r)) ** 2, np.ones(64) / 64, mode="same")    # smoothed spectrum (~100 Hz)
            ctl = []
            for _ in range(50):
                z = np.fft.rfft(rng.standard_normal(N)) * np.sqrt(env)
                ctl.append(np.max(np.abs(acf(np.fft.irfft(z, N)))))
            p95 = float(np.percentile(ctl, 95))
            loopish = sorted(set([int(v) for v in P] + [int(abs(u - v)) for u in P for v in P if u != v]))
            tag = lambda t: "*" if any(t % d == 0 and t // d <= 6 for d in loopish if d > 0) or any(abs(t - d) <= 1 for d in loopish) else ""
            res[f"{P[i]}{'LR'[c]}"] = dict(top=[(int(t + 10), float(ac[t])) for t in top], ctrl_max_p95=p95)
            say(f"  {P[i]} {'LR'[c]}  top lags " + "  ".join(f"{t + 10}{tag(t + 10)}:{ac[t]:+.3f}" for t in top) + f"   control max (95th pct) {p95:.3f}")
    say(f"  * = a loop length, a difference of two, or a small multiple of a difference ({sorted(set(int(abs(u - v)) for u in P for v in P if u != v))})")
    json.dump(res, open(os.path.join(OUT, "k5_resid.json"), "w"), indent=1)


def stage_settings(a, loops, say, taps, shift, g460, amp, A):
    from scipy.signal import butter
    info = json.load(open(os.path.join(ROOT, "analysis", "q1", "data.json")))
    law = json.load(open(os.path.join(ROOT, "analysis", "q1", "laws.json")))
    law["lf"] = {float(k): v for k, v in law["lf"].items()}
    ext = q1.EXT
    bands = dict(lf=(40, 250), hf=(4000, 10000))
    Bb = {k: B * np.abs(sosfreqz(butter(4, f, btype="band", fs=q1.FS, output="sos"), worN=2 * np.pi * np.fft.rfftfreq(N))[1]) ** 2 for k, f in bands.items()}
    Bb["full"] = B
    fc = dict(lf=np.sqrt(40 * 250), hf=np.sqrt(4000 * 10000))

    def bandk(x, k):
        b = np.zeros(N); b[PAD:PAD + len(x)] = x
        return np.fft.irfft(np.fft.rfft(b) * Bb[k], N)

    def tap_gains(x, i, c, k):
        L = len(x); y = bandk(x, k); cols = []
        for t, _ in taps[0, i]:
            tt = t + (shift.get((i, t), 0) if c else 0)
            for u in (tt, tt + 460):
                d = np.zeros(L)
                if u < L:
                    d[u] = 1
                cols.append(bandk(d, k))
        M = np.array(cols).T
        g, *_ = np.linalg.lstsq(M, y, rcond=None)
        r = y - M @ g
        return g[0::2], float(1 - (r @ r) / (y @ y))
    keys = sorted(os.path.splitext(f)[0] for f in os.listdir(os.path.join(OUT, "v")) if f.endswith(".npy"))
    rows = []
    for key in keys:
        xs, _ = load(os.path.join(OUT, "v", key + ".npy"), ext)
        lab = dict(lf=info[key]["lf"], rt=info[key]["rt"], hf=info[key]["hf"])
        Hs = q1.G_law(lab, law, 2 * np.pi * np.array([fc["lf"], 1000.0, fc["hf"]]) / q1.FS)
        for i in loops:
            pred = {k: float(20 * np.log10(np.abs(Hs[i][j]) / np.abs(Hs[i][1]))) for j, k in ((0, "lf"), (2, "hf"))}
            for c in (0, 1):
                x = xs[c, i]; L = len(x)
                y = decomb(band(x), g460)
                _, sel, _ = omp(y, L, np.array([0.0]), NTAP, [NTAP])
                flat = [t + (shift.get((i, t), 0) if c else 0) for t, _ in taps[0, i]]
                found = sum(any(abs(t - f) <= 1 for f in flat) for _, t in sel)
                gf, ex_flat = tap_gains(x, i, c, "full")
                gl, _ = tap_gains(x, i, c, "lf")
                gh, _ = tap_gains(x, i, c, "hf")
                sg = [float(np.sign(v)) for _, v in taps[0, i]]
                rows.append(dict(key=key, **lab, loop=int(P[i]), ch="LR"[c], found=int(found), explained_flat=ex_flat,
                                 gain_db=[float(20 * np.log10(abs(g) / amp)) for g in gf],
                                 sign_ok=bool(all(np.sign(g) == s for g, s in zip(gf, sg))),
                                 lf_db=[float(20 * np.log10(abs(l / f))) for l, f in zip(gl, gf)],
                                 hf_db=[float(20 * np.log10(abs(h / f))) for h, f in zip(gh, gf)], pred=pred))
    say(f"[settings] K7: the flat taps in the per-setting free injections ({len(keys)} settings x 4 loops x 2 channels = {len(rows)})")
    found = np.array([r["found"] for r in rows]); ex = np.array([r["explained_flat"] for r in rows])
    say(f"  5-atom OMP positions within +-1 of the flat taps: 5/5 in {np.mean(found == 5) * 100:.1f} %, >= 4/5 in {np.mean(found >= 4) * 100:.1f} % of excitations")
    say(f"  energy the flat tap set explains (main + partner atoms, free gains): median {100 * np.median(ex):.1f} %, 10th pct {100 * np.percentile(ex, 10):.1f} %")
    say(f"  signs as in the flat injection: {np.mean([r['sign_ok'] for r in rows]) * 100:.1f} % of excitations")
    gd = np.array([r["gain_db"] for r in rows])
    say(f"  tap gains re the flat magnitude: median {np.median(gd):+.2f} dB, IQR {np.percentile(gd, 25):+.2f}..{np.percentile(gd, 75):+.2f} dB, "
        f"spread across the 5 taps of one excitation (std) median {np.median(np.std(gd, axis=1)):.2f} dB")
    for name, f in (("LF", "lf"), ("HF", "hf"), ("RT", "rt")):
        vals = sorted(set(r[f] for r in rows))
        say(f"  by {name}: " + "  ".join(f"{v}: found {np.mean([r['found'] for r in rows if r[f] == v]):.2f}/5, gain {np.median([np.median(r['gain_db']) for r in rows if r[f] == v]):+.2f} dB, "
                                         f"flat set {100 * np.median([r['explained_flat'] for r in rows if r[f] == v]):.1f} %" for v in vals))
    say("  long RT only (>= 1.4 s), where the per-setting fits are well determined:")
    lr = [r for r in rows if r["rt"] >= 1.4]
    for name, f in (("LF", "lf"), ("HF", "hf")):
        vals = sorted(set(r[f] for r in lr))
        say(f"    by {name}: " + "  ".join(f"{v}: found {np.mean([r['found'] for r in lr if r[f] == v]):.2f}/5, gain {np.median([np.median(r['gain_db']) for r in lr if r[f] == v]):+.2f} dB, "
                                           f"flat set {100 * np.median([r['explained_flat'] for r in lr if r[f] == v]):.1f} %" for v in vals))
    say("\n[settings] K4: per-tap band gain re full band; spread across the 5 taps (max - min, dB) against the split K4 predicts (|H_i| in band re 1 kHz)")
    for f, b in (("lf", "LF"), ("hf", "HF")):
        vals = sorted(set(r[f] for r in rows))
        for v in vals:
            rr = [r for r in rows if r[f] == v]
            spread = [max(r[f + "_db"]) - min(r[f + "_db"]) for r in rr]
            pred = [abs(r["pred"][f]) for r in rr]
            say(f"  {b} {v}: measured spread median {np.median(spread):.2f} dB (90th pct {np.percentile(spread, 90):.2f}); "
                f"K4 predicts {np.median(pred):.2f} dB (median; max {np.max(pred):.2f})")
    say("  per-tap profiles, RT >= 1.4 s: mean band gain of each tap (dB) at a step minus at LF 1.0 / HF 1.0; K4 needs a two-level pattern "
        "with a gap of about the predicted |H|, the same taps in L and R")
    prof = {}
    for f, ref in (("lf", 1.0), ("hf", 1.0)):
        for v in sorted(set(r[f] for r in lr)):
            if v == ref:
                continue
            for lp in sorted(set(r["loop"] for r in lr)):
                line = []
                for ch in "LR":
                    sel_v = [r for r in lr if r[f] == v and r["loop"] == lp and r["ch"] == ch and r["lf" if f == "hf" else "hf"] == 1.0]
                    sel_0 = [r for r in lr if r[f] == ref and r["loop"] == lp and r["ch"] == ch and r["lf" if f == "hf" else "hf"] == 1.0]
                    if not sel_v or not sel_0:
                        continue
                    d = np.mean([r[f + "_db"] for r in sel_v], axis=0) - np.mean([r[f + "_db"] for r in sel_0], axis=0)
                    sd = np.std([r[f + "_db"] for r in sel_v], axis=0)
                    pr = np.median([r["pred"][f] for r in sel_v])
                    prof[f"{f}{v}_{lp}{ch}"] = dict(diff_db=[float(u) for u in d], sd_db=[float(u) for u in sd], pred_db=float(pr), n=len(sel_v))
                    line.append(f"{ch}: " + " ".join(f"{u:+5.1f}" for u in d) + f" (sd {np.median(sd):.1f}, n {len(sel_v)})")
                if line:
                    say(f"    {f.upper()} {v}, loop {lp}: K4 gap {prof[f'{f}{v}_{lp}L']['pred_db']:+5.1f} dB   " + "   ".join(line))
    json.dump(dict(rows=rows, profiles=prof), open(os.path.join(OUT, "settings.json"), "w"), indent=0)


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="k1,k2,structure,k5,settings")
    ap.add_argument("--src", default="branches_flat_0_ext1700")
    ap.add_argument("--ext", type=int, default=1700)
    ap.add_argument("--loops", default="0,1,2,3")
    ap.add_argument("--kmax", type=int, default=64)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--mstep", type=int, default=1, help="k2: delay step of the scan (1 = every M)")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    OUT = a.out
    stages = a.stage.split(",")
    loops = [int(s) for s in a.loops.split(",")]
    xs, F = load(os.path.join(OUT, a.src + ".npy"), a.ext)
    ys = {q: band(x) for q, x in xs.items()}
    logs = {}

    def sayer(name):
        logs[name] = open(os.path.join(OUT, name + ".log"), "w")

        def say(s=""):
            print(s, flush=True); logs[name].write(s + "\n")
        return say
    k1 = stage_k1(a, xs, ys, loops, sayer("k1_sparsity") if "k1" in stages else (lambda s="": None), np.random.default_rng(a.seed))
    if "k2" in stages:
        stage_k2(a, xs, ys, loops, sayer("k2_stages"), k1[2])
    if "structure" in stages:
        stage_structure(a, xs, ys, F, loops, sayer("structure"), *k1)
    if "k5" in stages:
        stage_k5(xs, ys, loops, sayer("k5_resid"), *k1)
    if "settings" in stages:
        stage_settings(a, loops, sayer("settings"), *k1)


if __name__ == "__main__":
    main()
