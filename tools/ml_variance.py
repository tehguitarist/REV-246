"""Q3b: which part of programs M and L varies between captures (NEXT_STEPS §5 item 5, docs/ml_structure.md §9).

    .venv/bin/python tools/ml_variance.py [--stage rot|report|all] [--prog ML] [--workers 4] [--only N] [--chunk i/n]

Works on the residual capture - model of tools/ml_lines.py (`fitted`: stored clock and gains, nothing refitted), which
past the first pass is mostly the capture's own random part (ml_structure §8.1). Each candidate for what varies is a
derivative of the model; if that parameter wanders during the sine sweep, the sweep turns the wander into a kernel on
the derivative's response (a slow modulation becomes a pair of ghosts at lags that shrink with frequency). So per
capture and window the residual is fitted as kernel * derivative (ridge, lags +-LK), on one channel, and the fit is
scored on the *other* channel: L and R read the same four lines through different taps, so a line-side variation
predicts across channels and nothing else does.

The candidate found: one rotation angle theta of both 2 x 2 blocks of the coupling, F(theta) = F R(theta) per block
(dF/dtheta = F_b J, J = [[0, -1], [1, 0]]); the coupling stays orthogonal at every instant. Controls: a common gain
change of the coupling (dF = F), the two blocks rotated in opposite senses, and the two blocks with separate angles.

The model is 1^T (I - Lambda F)^-1 T (ml_lines). A k-pass path has k coupling events; event e (model order, e = 1
next to T) is, in time, q = k - e + 1 line lengths after the input (the data put the rotation at the line outputs and
the taps after it, ml_structure §9.3). recirc k is fitted with one kernel per q; each kernel's correlation with
recirc1's (q = 1) estimates theta's autocorrelation at (q - 1) line lengths, in seconds that differ for M and L.

  rot     per AMPL group (one checkpoint file each in analysis/ml_variance/rot/; finished groups are skipped): per
          capture, recirc1-3 (ml_lines windows, up to the fade): base null, rotation kernel held out across channels,
          the controls (recirc1), one kernel vs one per q (recirc2-3), the kernel's rms (theta rms, 1-9 kHz), its lag
          spread and symmetry per octave (recirc1), the q-kernels' correlation with recirc1's, and recirc1's kernel.
  report  medians by program and control, theta's autocorrelation in seconds, lag x frequency per octave.
  sim     the accuracy test for the M/L reverb (NEXT_STEPS §2, metric 3): the model in the time domain with
          theta(t) = a sin(2 pi f_m t + phi) (read-out form, below), measured through a simulated exponential sweep
          (10 Hz-13 kHz at the device rate, time constant SWEEP_L; the sweep is unknown, only f_m x SWEEP_L is set by the
          data), rendered N_SIM times with random phi; pair nulls per window against Tomkin's AMPL pairs of the same
          groups, and the estimator of `rot` run on the renders (where theta is known) as its calibration.

Read-out form (the transpose of ml_lines' 1^T (I - Lambda F)^-1 T, which the data prefer, §9.3): the diffused input
goes into all four lines, u_i = x + sum_j F(theta)_ji v_j with v_j = G_j * u_j delayed P_j, and each channel's taps
read u_i. With theta = 0 it reproduces ml_lines' model to about -100 dB.
"""
import argparse, json, os, sys
for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                          # one thread per worker: the pool supplies the parallelism
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import ml_carryover as cc  # noqa: E402
import ml_lines as mll  # noqa: E402
import ml_structure as ml  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_variance")
FS, K0 = cc.FS, cc.K0
J = np.array([[0, -1], [1, 0]])
BLOCKS = (([0, 2], [1, 3]), ([1, 3], [0, 2]))               # (lines, loops) of the two blocks of F (ml_structure §7.3)
LK = {"M": 300, "L": 600}                                   # kernel half-length (lags): the ghosts reach ~ +-200 / 400 at 1 kHz
LAM = 1e-2                                                  # ridge, relative to the mean diagonal of X^T X
WIN = ("recirc1", "recirc2", "recirc3")
BAND_K = butter(4, [1000, 9000], btype="band", fs=FS, output="sos")    # where kernels are compared (lags well resolved)
OCT_K = [1000, 2000, 4000, 8000]


# ---------------------------------------------------------------- derivatives
def dF_of(F, kind):
    """The coupling's derivative for one candidate: 'rot' (one angle, both blocks), 'rot1' / 'rot2' (one block),
    'anti' (blocks in opposite senses), 'gain' (dF = F)."""
    if kind == "gain":
        return F.copy()
    D = np.zeros((4, 4))
    for b, (r, c) in enumerate(BLOCKS):
        s = {"rot": 1, "anti": (1, -1)[b], "rot1": (1, 0)[b], "rot2": (0, 1)[b]}[kind]
        D[np.ix_(r, c)] = s * F[np.ix_(r, c)] @ J
    return D


def paths(M, n, key, kinds=("rot", "gain", "anti", "rot1", "rot2"), kmax=3):
    """{(kind, k, e): (n, 2)} the k-pass part of the model's derivative with the coupling perturbed at event e only
    (unit gain, band-limited, nominal clock, as ml_lines.Model.parts). 'rot' gets k = 1..kmax, the controls k = 1."""
    lf, rt, hf = (float(x) for x in key.split("_"))
    coefs = mll.line_coefs(M.P, rt, lf, hf)
    N = 1 << int(np.ceil(np.log2(n + 4096)))
    r = 10 ** (mll.ALIAS_DB / 20 / N)
    w = 2 * np.pi * np.fft.rfftfreq(N)
    zi = np.exp(-1j * w) / r
    zd = lambda d: np.exp(-1j * w * d) / r ** d
    D = np.ones(len(w), complex)
    Ms = cc.DIFF[M.prog]
    for k, Mk in enumerate(Ms):
        bb = M.b.get(key, 0.0) if k == len(Ms) - 1 else 0.0
        Fz = zd(Mk) * ((1 - bb) + bb * zi)
        D *= (-cc.G_AP + Fz) / (1 - cc.G_AP * Fz)
    Hk = sum(hk * zd(K0 + k - cc.KH) for k, hk in enumerate(M.h))
    Lam = mll.line_resp(coefs, zi) * np.stack([zd(Pj) for Pj in M.P], 1)
    T = np.zeros((len(w), 4, 2), complex)
    for ch in (0, 1):
        for p_, sg, i in zip(np.round(M.pos[ch]), M.sgn[ch], M.loop[ch]):
            T[:, i, ch] += sg * zd(p_)
    LF = Lam[:, :, None] * M.F[None]
    out = {}
    for kind in kinds:
        LdF = Lam[:, :, None] * dF_of(M.F, kind)[None]
        for k in range(1, (kmax if kind == "rot" else 1) + 1):
            for e in range(1, k + 1):
                X = T
                for ev in range(1, k + 1):
                    X = np.einsum("bij,bjc->bic", LdF if ev == e else LF, X)
                out[(kind, k, e)] = X.sum(1)
    grow = r ** np.arange(n)[:, None]
    return {kk: sosfiltfilt(ml.BAND, np.fft.irfft(X * D[:, None] * Hk[:, None], N, axis=0)[:n] * grow, axis=0)
            for kk, X in out.items()}


# ---------------------------------------------------------------- kernel fits
def lagged(col, s, lk):
    """(len(s), 2 lk + 1): col shifted by -lk..lk over the slice s (lag l: col[s - l])."""
    v = sliding_window_view(col[s.start - lk:s.stop + lk], s.stop - s.start)
    return np.ascontiguousarray(v[::-1].T)


def ridge(X, y):
    A = X.T @ X
    return np.linalg.solve(A + LAM * np.trace(A) / len(A) * np.eye(len(A)), X.T @ y)


def design(groups, sl, lk):
    """Per channel, the kernel design: for each group of columns (summed), its lags."""
    return [np.hstack([lagged(sum(c[:, ch] for c in g), sl[ch], lk) for g in groups]) for ch in (0, 1)]


def heldout(X, R, E):
    """Fit on one channel, score on the other (both ways); mean null re the capture (dB)."""
    s = []
    for a, b in ((0, 1), (1, 0)):
        k = ridge(X[a], R[a])
        s.append(10 * np.log10(np.sum((R[b] - X[b] @ k) ** 2) / E[b]))
    return float(np.mean(s))


def kcorr(a, b):
    a, b = sosfiltfilt(BAND_K, a), sosfiltfilt(BAND_K, b)
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def _rot_job(a):
    """One AMPL group: every capture's residual against the rotation and the controls (module docstring)."""
    prog, key, path = a
    M = mll.Model(prog)
    fit = mll.fitted(prog, key)
    rows = {r["amp"]: r for r in json.load(open(mll.group_path(prog, key)))}
    n = len(next(iter(fit.values()))[0])
    W = mll.windows(prog)
    lk = LK[prog]
    lags = np.arange(-lk, lk + 1)
    PC = paths(M, n, key)
    out = []
    for amp, (y, mdl, end) in fit.items():
        r_ = rows[amp]
        gw = np.array([r_["gains"][ch][1] for ch in (0, 1)])[None]
        C = {kk: mll.warp(v, *r_["clock"]) * gw for kk, v in PC.items()}
        res = y - mdl
        row = dict(key=key, amp=amp)
        K1 = None
        for k, wn in enumerate(WIN, 1):
            sl = [slice(K0 + W[ch][wn][0], min(K0 + W[ch][wn][1], end[ch])) for ch in (0, 1)]
            if min(s.stop - s.start for s in sl) < 0.8 * min(M.P):
                break
            R = [res[sl[ch], ch] for ch in (0, 1)]
            E = [float(np.sum(y[sl[ch], ch] ** 2)) for ch in (0, 1)]
            w_ = dict(base=float(np.mean([10 * np.log10(np.sum(R[ch] ** 2) / E[ch]) for ch in (0, 1)])))
            byq = [[C[("rot", kk, kk - q + 1)] for kk in range(q, k + 1)] for q in range(1, k + 1)]   # physical event q
            Xq = design(byq, sl, lk)
            w_["rot"] = heldout(Xq, R, E)
            Kq = ridge(np.vstack(Xq), np.r_[R[0], R[1]]).reshape(k, -1)
            if k == 1:
                K1 = Kq[0]
                for ctl, g in (("gain", [[C[("gain", 1, 1)]]]), ("anti", [[C[("anti", 1, 1)]]]),
                               ("two", [[C[("rot1", 1, 1)]], [C[("rot2", 1, 1)]]])):
                    w_[ctl] = heldout(design(g, sl, lk), R, E)
                w_["theta_rms"] = float(np.sqrt(np.sum(sosfiltfilt(BAND_K, K1) ** 2)))
                w_["kernel"] = [round(float(x), 6) for x in K1]
                w_["oct"] = {}
                for f in OCT_K:                               # lag spread and symmetry per octave (joint fit, both channels)
                    sos = mll.OCT_SOS[f]
                    cb = sosfiltfilt(sos, C[("rot", 1, 1)], axis=0)
                    rb = sosfiltfilt(sos, res, axis=0)
                    Xb = np.vstack([lagged(cb[:, ch], sl[ch], lk) for ch in (0, 1)])
                    kb = ridge(Xb, np.r_[rb[sl[0], 0], rb[sl[1], 1]])
                    e = kb ** 2 / np.sum(kb ** 2)
                    al = np.abs(lags)
                    o = np.argsort(al, kind="stable")
                    ce = np.cumsum(e[o])
                    w_["oct"][str(f)] = dict(lag50=int(al[o][np.searchsorted(ce, 0.5)]), lag90=int(al[o][np.searchsorted(ce, 0.9)]),
                                             sym=float(np.sum(((kb + kb[::-1]) / 2) ** 2) / np.sum(kb ** 2)))
            else:
                w_["one"] = heldout(design([sum(byq, [])], sl, lk), R, E)
                w_["corr_r1"] = [kcorr(K1, Kq[q]) for q in range(k)]
                # the opposite time order (rotation at the line inputs, taps before it): event e is e line lengths in
                byq_alt = [[C[("rot", kk, q)] for kk in range(q, k + 1)] for q in range(1, k + 1)]
                w_["rot_alt"] = heldout(design(byq_alt, sl, lk), R, E)
            row[wn] = w_
        out.append(row)
    json.dump(out, open(path, "w"))
    return prog, key


def group_dir(tag=""):
    return os.path.join(OUT, f"rot{tag}")


def rot_stage(progs, workers, only=None, chunk=None):
    """Every AMPL group of the programs (smoke: `only` per program), one checkpoint file each, resumable; chunk (i, n)
    as in ml_lines."""
    tag = "_smoke" if only else ""
    os.makedirs(group_dir(tag), exist_ok=True)
    jobs = []
    for p in progs:
        ks = sorted(cc.groups(p))
        if only:
            ks = ks[::max(1, len(ks) // only)][:only]
        jobs += [(p, k, os.path.join(group_dir(tag), f"{p}_{k}.json")) for k in ks if os.path.exists(mll.group_path(p, k))]
    if chunk:
        jobs = jobs[chunk[0]::chunk[1]]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    print(f"rot: {len(jobs)} groups to do", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        fut = [ex.submit(_rot_job, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            p, k = f.result()
            print(f"  {i + 1}/{len(jobs)} {p} {k}", flush=True)


# ---------------------------------------------------------------- report
def load_rows(prog, tag=""):
    import glob
    return [r for f in sorted(glob.glob(os.path.join(group_dir(tag), f"{prog}_*.json"))) for r in json.load(open(f))]


def report_stage(progs, tag=""):
    for p in progs:
        rows = load_rows(p, tag)
        if not rows:
            continue
        P = np.array(cc.LINES[p]) / FS
        lab = lambda r: dict(zip(("lf", "rt", "hf"), (float(x) for x in r["key"].split("_"))))
        print(f"\n## {p}: {len({r['key'] for r in rows})} groups, {len(rows)} captures; held out across channels, "
              f"mean of L->R and R->L, dB re the capture; median [10-90 %]")
        q = lambda x: f"{np.median(x):6.1f} [{np.percentile(x, 10):6.1f} … {np.percentile(x, 90):6.1f}]"
        for wn in WIN:
            sel = [r[wn] for r in rows if wn in r]
            if not sel:
                continue
            cells = [f"base {q([s['base'] for s in sel])}", f"rotation {q([s['rot'] for s in sel])}"]
            if wn == "recirc1":
                cells += [f"{c} {np.median([s[c] for s in sel]):6.1f}" for c in ("two", "anti", "gain")]
            else:
                cells += [f"one kernel {np.median([s['one'] for s in sel]):6.1f}",
                          f"opposite time order {np.median([s['rot_alt'] for s in sel]):6.1f}"]
            print(f"- {wn} (n {len(sel)}): " + "; ".join(cells))
        rows = [r for r in rows if "recirc1" in r]             # short RT: some captures fade before recirc1
        r1 = [r["recirc1"] for r in rows]
        print("- recirc1 by control: base / rotation / gain control (held out, median), theta rms")
        for ctl in ("rt", "lf", "hf"):
            for x in sorted({lab(r)[ctl] for r in rows}):
                s = [r["recirc1"] for r in rows if lab(r)[ctl] == x]
                print(f"  {ctl.upper()} {x:<5g} {np.median([v['base'] for v in s]):6.1f} / {np.median([v['rot'] for v in s]):6.1f} / "
                      f"{np.median([v['gain'] for v in s]):6.1f}; {np.median([v['theta_rms'] for v in s]):.4f}  ({len(s)})")
        th = np.array([s["theta_rms"] for s in r1])
        print(f"- theta rms (1-9 kHz kernel norm, rad): {np.median(th):.4f} [{np.percentile(th, 10):.4f} … {np.percentile(th, 90):.4f}]")
        print("- recirc1 kernel per octave: median |lag| 50 % / 90 % of its energy (samples), lag50 x f (s Hz), symmetric fraction")
        for f in OCT_K:
            l5 = np.array([s["oct"][str(f)]["lag50"] for s in r1])
            l9 = np.array([s["oct"][str(f)]["lag90"] for s in r1])
            sy = np.array([s["oct"][str(f)]["sym"] for s in r1])
            print(f"  {f:5d} Hz: {np.median(l5):5.0f} / {np.median(l9):5.0f}; {np.median(l5) / FS * f:5.2f}; {np.median(sy):.3f}")
        print(f"- theta's autocorrelation (kernel of event q vs recirc1's, 1-9 kHz; q = 1 is the same instant, so rho = "
              f"corr(q) / corr(1)); line lengths {P.min():.3f}-{P.max():.3f} s:")
        for wn, k in (("recirc2", 2), ("recirc3", 3)):
            sel = np.array([r[wn]["corr_r1"] for r in rows if wn in r])
            if not len(sel):
                continue
            med = np.median(sel, 0)
            print(f"  {wn} (n {len(sel)}): corr " + " ".join(f"q{i + 1} {x:+.2f}" for i, x in enumerate(med))
                  + " | rho " + " ".join(f"{i * P.mean():.2f} s {x / med[0]:+.2f}" for i, x in enumerate(med)))
        by = {}
        for r in rows:
            by.setdefault(r["key"], []).append(np.array(r["recirc1"]["kernel"]))
        cc_ = [kcorr(v[i], v[j]) for v in by.values() for i in range(len(v)) for j in range(i + 1, len(v))]
        print(f"- recirc1 kernels of different captures of one group: corr median {np.median(cc_):+.2f} "
              f"[{np.percentile(cc_, 10):+.2f} … {np.percentile(cc_, 90):+.2f}] (n {len(cc_)})")
        sine_test(p)


def sine_test(p):
    """Periodic theta with a random phase per capture: every capture is the same distance from the theta = 0 model
    (narrow spread), but two captures differ by 2 (1 - cos dphi) times that, dphi uniform (wide spread, deep tail).
    From ml_lines' law files: model and repeat nulls (R) per window, percentiles, and the sine's prediction for the
    repeat nulls from the model median (median +3.0, 90 % +5.9, 10 % -10.1, 1 % -30 dB)."""
    rows = mll.load_rows(p)
    pct = (1, 10, 50, 90)
    pred = {k: 10 * np.log10(2 * (1 - np.cos(np.pi * k / 100))) for k in pct}     # quantile k of the pair null: |dphi| = pi k / 100
    print("- periodic theta test (ml_lines law files, R): percentiles 1 / 10 / 50 / 90, dB; IQR")
    for w in ("recirc1", "recirc2", "recirc3"):
        ok = [r for r in rows if w in r["law"] and r["law"][w][1] is not None and r["repeat"][w][1] is not None]
        m = np.array([r["law"][w][1] for r in ok])
        q = np.array([r["repeat"][w][1] for r in ok])
        iqr = lambda x: np.subtract(*np.percentile(x, [75, 25]))
        print(f"  {w}: model {' '.join(f'{np.percentile(m, k):6.1f}' for k in pct)} ({iqr(m):.1f}); repeat "
              f"{' '.join(f'{np.percentile(q, k):6.1f}' for k in pct)} ({iqr(q):.1f}); sine prediction "
              f"{' '.join(f'{np.median(m) + pred[k]:6.1f}' for k in pct)}  (n {len(ok)})")


# ---------------------------------------------------------------- sim stage
THETA_A, THETA_FM = 0.155, 0.83                             # rad (rms 0.11), Hz: from the rot stage (ml_structure §9.4)
SWEEP_L, SWEEP_F = 3.4, (10.0, 13000.0)                     # s, Hz: f_m x SWEEP_L ~ 2.85 (lag x f of the ghosts, §9.4)
N_SIM = 6
SIM_KEYS = {"M": ["1.2_5.0_0.9", "2.0_0.8_0.9", "0.7_3.0_0.2", "0.1_1.8_0.4"], "L": ["1.2_20.0_1.0", "1.2_4.5_0.8", "2.0_2.6_0.5", "0.1_10.0_0.3"]}
SIM_TW = [(0.15, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 2.0), (2.0, 4.0)]


def F_theta(F, th):
    """(n, 4, 4): F with both blocks rotated by th[n] (F_b R(th))."""
    c, s = np.cos(th), np.sin(th)
    R = np.stack([np.stack([c, -s], -1), np.stack([s, c], -1)], -2)
    out = np.zeros((len(th), 4, 4))
    for r, cl in BLOCKS:
        ix = np.ix_(r, cl)
        out[:, ix[0], ix[1]] = np.einsum("ab,nbc->nac", F[ix], R)
    return out


def line_ba(coefs):
    """Per line, (b, a) of G_j (ml_lines.line_resp) for lfilter."""
    a_ = mll.LAW["a"]
    out = []
    for g, r0, p in coefs:
        b0, b1 = (r0 * (1 - a_) + (1 + a_)) / 2, (r0 * (1 - a_) - (1 + a_)) / 2
        out.append((g * (1 - p) * np.array([b0, b1]), np.convolve([1, -a_], [1, -p])))
    return out


def core(M, key, x, theta):
    """Read-out form, in blocks of the shortest line (every line output in a block was written before it)."""
    from scipy.signal import lfilter, fftconvolve
    lf, rt, hf = (float(v) for v in key.split("_"))
    ba = line_ba(mll.line_coefs(M.P, rt, lf, hf))
    n, P = len(x), M.P
    u = np.zeros((n, 4))
    zi = [np.zeros(2) for _ in range(4)]
    for b0 in range(0, n, min(P)):
        b1 = min(b0 + min(P), n)
        v = np.zeros((b1 - b0, 4))
        for j in range(4):
            src = np.zeros(b1 - b0)
            lo, hi = b0 - P[j], b1 - P[j]
            if hi > 0:
                src[max(lo, 0) - lo:] = u[max(lo, 0):hi, j]
            v[:, j], zi[j] = lfilter(*ba[j], src, zi=zi[j])
        u[b0:b1] = x[b0:b1, None] + np.einsum("mji,mj->mi", F_theta(M.F, theta[b0:b1]), v)
    y = np.zeros((n, 2))
    for ch in (0, 1):
        for i in range(4):
            sel = M.loop[ch] == i
            pos = np.round(M.pos[ch][sel]).astype(int)
            h = np.zeros(pos.max() + 1)
            np.add.at(h, pos, M.sgn[ch][sel])
            y[:, ch] += fftconvolve(u[:, i], h)[:n]
    return y


def sweep_measure(M, key, theta_fn, tail_s=6.0):
    """The reverb part (taps through lines, no cluster) as an exponential-sweep measurement sees it: K0 grid, band-limited."""
    from scipy.signal import fftconvolve
    f1, f2 = SWEEP_F
    T = SWEEP_L * np.log(f2 / f1)
    t = np.arange(int(T * FS)) / FS
    s = np.sin(2 * np.pi * f1 * SWEEP_L * (np.exp(t / SWEEP_L) - 1))
    m = int(0.05 * FS)
    s[:m] *= np.sin(np.linspace(0, np.pi / 2, m)) ** 2
    s[-m:] *= np.sin(np.linspace(np.pi / 2, 0, m)) ** 2
    inv = s[::-1] * np.exp(-t / SWEEP_L)
    tail = int(tail_s * FS)
    x = np.r_[s, np.zeros(tail)]
    for k, Mk in enumerate(cc.DIFF[M.prog]):
        x = cc.ap(x, Mk, cc.G_AP, M.b.get(key, 0.0) if k == len(cc.DIFF[M.prog]) - 1 else 0.0)
    y = core(M, key, x, theta_fn(np.arange(len(x)) / FS))
    y = np.stack([np.convolve(y[:, ch], M.h)[cc.KH:cc.KH + len(y)] for ch in (0, 1)], 1)
    ref = fftconvolve(s, inv)
    i0 = int(np.argmax(np.abs(ref)))
    out = np.zeros((tail + K0, 2))
    out[K0:] = np.stack([fftconvolve(y[:, ch], inv)[i0:i0 + tail] for ch in (0, 1)], 1) / np.abs(ref).max()
    return sosfiltfilt(ml.BAND, out, axis=0)


def pair_nulls(prog, y1, y2, end):
    """Pair null per window, one gain per window and channel: recirc1-3 (R) and SIM_TW (mean of L and R), dB."""
    W = mll.windows(prog)
    def nl(a0, a1, chs):
        v = []
        for ch in chs:
            x1, x2 = y1[a0:a1, ch], y2[a0:a1, ch]
            v.append(10 * np.log10(np.sum((x1 - np.dot(x1, x2) / np.dot(x2, x2) * x2) ** 2) / np.sum(x1 ** 2)))
        return float(np.mean(v))
    out = {}
    for wn in WIN:
        a0, a1 = K0 + W[1][wn][0], min(K0 + W[1][wn][1], end)
        out[wn] = nl(a0, a1, (1,)) if a1 - a0 > 1000 else None
    for t0, t1 in SIM_TW:
        a0, a1 = K0 + int(t0 * FS), min(K0 + int(t1 * FS), end)
        out[f"{t0}-{t1}"] = nl(a0, a1, (0, 1)) if a1 - a0 > 2000 else None
    return out


def _sim_job(a):
    prog, key, path, seed = a
    M = mll.Model(prog)
    fit = mll.fitted(prog, key)
    amps = list(fit)
    row = dict(key=key, tomkin=[], sim=[], calib=[])
    for i in range(len(amps)):
        for j in range(i + 1, len(amps)):
            e = min(min(fit[amps[i]][2]), min(fit[amps[j]][2]))
            row["tomkin"].append(pair_nulls(prog, fit[amps[i]][0], fit[amps[j]][0], e))
    rng = np.random.default_rng(seed)
    irs = [sweep_measure(M, key, lambda t, ph=rng.uniform(0, 2 * np.pi): THETA_A * np.sin(2 * np.pi * THETA_FM * t + ph))
           for _ in range(N_SIM)]
    for i in range(N_SIM):
        for j in range(i + 1, N_SIM):
            row["sim"].append(pair_nulls(prog, irs[i], irs[j], len(irs[i]) - 1000))
    # calibration: the rot stage's estimator on two renders (model at theta = 0 as the reference, one gain per channel)
    n = len(irs[0])
    _, R0 = M.parts(n, key)
    W = mll.windows(prog)
    lk = LK[prog]
    PC = paths(M, n, key, kinds=("rot",))
    for y in irs[:2]:
        g = np.array([np.dot(y[K0:, ch], R0[K0:, ch]) / np.dot(R0[K0:, ch], R0[K0:, ch]) for ch in (0, 1)])[None]
        C = {kk: v * g for kk, v in PC.items()}
        res = y - R0 * g
        cal, K1 = {}, None
        for k, wn in enumerate(WIN, 1):
            sl = [slice(K0 + W[ch][wn][0], K0 + W[ch][wn][1]) for ch in (0, 1)]
            Rr = [res[sl[ch], ch] for ch in (0, 1)]
            X = design([[C[("rot", kk, kk - q + 1)] for kk in range(q, k + 1)] for q in range(1, k + 1)], sl, lk)
            Kq = ridge(np.vstack(X), np.r_[Rr[0], Rr[1]]).reshape(k, -1)
            if k == 1:
                K1 = Kq[0]
                cal["theta_rms"] = float(np.sqrt(np.sum(sosfiltfilt(BAND_K, K1) ** 2)))
            else:
                cal[f"corr_{wn}"] = [kcorr(K1, Kq[q]) for q in range(k)]
        row["calib"].append(cal)
    json.dump(row, open(path, "w"))
    return prog, key


def sim_stage(progs, workers):
    os.makedirs(os.path.join(OUT, "sim"), exist_ok=True)
    jobs = [(p, k, os.path.join(OUT, "sim", f"{p}_{k}.json"), 1000 * i + 1) for p in progs for i, k in enumerate(SIM_KEYS[p])]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for f in as_completed([ex.submit(_sim_job, j) for j in jobs]):
            print("  sim", *f.result(), flush=True)
    print(f"\n## sim: theta = {THETA_A} sin(2 pi {THETA_FM} t + phi), sweep L = {SWEEP_L} s; pair nulls, dB, median [min … max] "
          f"over pairs (Tomkin: its AMPL pairs; sim: {N_SIM} renders)")
    for p in progs:
        P = np.mean(cc.LINES[p]) / FS
        for k in SIM_KEYS[p]:
            row = json.load(open(os.path.join(OUT, "sim", f"{p}_{k}.json")))
            print(f"- {p} {k}:")
            for src in ("tomkin", "sim"):
                cells = []
                for w in list(WIN) + [f"{a}-{b}" for a, b in SIM_TW]:
                    v = [r[w] for r in row[src] if r[w] is not None]
                    cells.append(f"{w} {np.median(v):6.1f} [{min(v):6.1f} … {max(v):6.1f}]" if v else f"{w}   -")
                print(f"  {src:6s} " + "; ".join(cells))
            th = [c["theta_rms"] for c in row["calib"]]
            c2 = np.mean([c["corr_recirc2"] for c in row["calib"]], 0)
            c3 = np.mean([c["corr_recirc3"] for c in row["calib"]], 0)
            print(f"  calibration (true theta rms {THETA_A / np.sqrt(2):.3f}, rho {np.cos(2 * np.pi * THETA_FM * P):+.2f} / "
                  f"{np.cos(4 * np.pi * THETA_FM * P):+.2f} at {P:.2f} / {2 * P:.2f} s): estimated theta rms {np.mean(th):.3f}; "
                  f"recirc2 corr {c2.round(2).tolist()}, recirc3 corr {c3.round(2).tolist()}")


# ---------------------------------------------------------------- main
def main():
    a = argparse.ArgumentParser()
    a.add_argument("--stage", default="rot", help="rot | report | sim | all")
    a.add_argument("--prog", default="ML")
    a.add_argument("--workers", type=int, default=4)
    a.add_argument("--only", type=int, default=None, help="smoke test: this many groups per program")
    a.add_argument("--chunk", default=None, help="i/n: run only every n-th group from the i-th")
    o = a.parse_args()
    os.makedirs(OUT, exist_ok=True)
    tag = "_smoke" if o.only else ""
    if o.stage in ("rot", "all"):
        rot_stage(o.prog, o.workers, o.only, tuple(int(x) for x in o.chunk.split("/")) if o.chunk else None)
    if o.stage in ("report", "all"):
        report_stage(o.prog, tag)
    if o.stage in ("sim", "all"):
        sim_stage(o.prog, o.workers)


if __name__ == "__main__":
    main()
