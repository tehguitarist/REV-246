"""NEXT_STEPS §5 item 5: a per-channel chain kernel for the Tomkin 96 kHz library, replacing the 5.2 ms R-only
deconvolution kernel (docs/ml_structure.md §10-11). Results: docs/ml_structure.md §12.

    .venv/bin/python tools/chain_kernel.py --stage data     # capture windows -> analysis/chain_kernel/data_*.npz (not committed)
    .venv/bin/python tools/chain_kernel.py --stage align    # per-capture tap-0 time and clock -> shifts.json, windows rebuilt
    .venv/bin/python tools/chain_kernel.py --stage fit      # held-out scores of the kernel and its variants -> fit.{json,log}
    .venv/bin/python tools/chain_kernel.py --stage write [--name chain2]
        -> analysis/library/tomkin96/kernels/<name>.npz and analysis/chain_kernel/<name>.json; use with EMT_KERNEL=<name>
           (tools/ml_structure.py: caches and results get the suffix _<name>)

Model. Before the reverb starts a capture holds only the cluster: per channel 4 sign-only taps at known device offsets
(analysis/ml_structure/summary.json), for M and L through a per-program output kernel (taps -8..8, fitted alongside;
the stored ones of analysis/ml_carryover/firstpass.json as the variant), all of one gain A per capture (both
channels), through the analogue and capture chain. The chain is c_ch = h_ch * HP: a free FIR h_ch per channel (LK device
samples after the tap, LEAD before) times the unit's analytic sub-20 Hz high-pass (the schematic's corners, 3.35 / 0.72
/ 0.23 / 0.16 Hz, analog_results.md §3.5; nothing fitted). The captures are divided by HP first (causal, from the file
start), so the FIR only holds the fast part and the high-pass's slow tail, far past any window, is known rather than
truncated. Since every tap is on the device grid, the capture is linear in h_ch sampled on a 4x device grid (4
polyphase problems with one Gram matrix). Each window runs from LEAD device samples before tap 0 to 6 before the reverb
onset (ml_carryover.ONSET): no reverb sample enters any fit. S, M and L (gaps 9 / 18 / 36 ms in L, 22 / 43 / 86 ms in R)
are fitted jointly, S weighted W_S. Per capture: tap 0's time and the capture clock (align stage, from 13 ppm; L =
R + 19.19 us, the cache convention), the gain, and a slow additive term {1, t, t^2} per channel, which is projected out
of model and data alike (it is the same shape in every capture but does not scale with AMPL: capture-side, not the
chain). Gains, output kernels and FIR alternate (ROUNDS).

Scores: nulls over the whole window per capture (gain, shift, clock and the slow term refitted), 100 Hz-10 kHz,
25-100 Hz, 10-12 kHz and 20 Hz-12.3 kHz, medians over the AMPL >= 60 captures of the settings not used in the fit
(odd-indexed per program); S also from M and L only. Reference: the repeat-capture null (AMPL 60 predicted by another
capture of its group with the same freedoms; S and M), and the library kernel ("old").
"""
import argparse, json, os, sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import soundfile as sf
from scipy.signal import bilinear, butter, lfilter, resample_poly, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import ml_carryover as cc  # noqa: E402
import ml_lines as mlines  # noqa: E402   (the Kaiser-sinc table of its clock warp)
import ml_structure as ml  # noqa: E402
from core import Kernel, subsample_peak  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "chain_kernel")
KDIR = os.path.join(ROOT, "analysis", "library", "tomkin96", "kernels")
OLD = os.path.join(KDIR, "kernel_96000.npz")
OLD_CACHE = os.path.join(ROOT, "analysis", "ml_structure", "cache")   # tap 0 and the initial gain come from it
FS96 = 96000
R4 = 91 / 100                   # capture samples per 4x-device sample (fs / fs_dev = 3.64 exactly)
LEAD = 25                       # device samples of kernel before the tap (25 -> 100 on the 4x grid -> 91 at 96 kHz)
END_MARGIN = 6                  # device samples kept clear of the reverb onset
SCHEMATIC_HZ = [3.35, 0.72, 0.23, 0.16]                 # the unit's coupling corners (analog_results.md §3.5)
HP_SETS = {"schematic": SCHEMATIC_HZ,                     # the high-pass divided out before the FIR fit ("var")
           "schematic+1.8": SCHEMATIC_HZ + [1.80],        # + Tomkin's chain corner as analog_results.md §3.5 measured it
           "half+1.8": [v / 2 for v in SCHEMATIC_HZ] + [1.80],
           "none": []}
VAR, LK, W_S, ROUNDS, FREE_HK = "schematic", 500, 3.0, 8, "shared"     # the written kernel
PPM0 = 13.0                     # starting capture clock (ppm; ml_lines' law-stage clocks are 11.6-13.6 in M and L)
SOS4 = {b: butter(4, lohi, btype="band", fs=4 * ml.FSD, output="sos")
        for b, lohi in {"band": (100, 10000), "low": (25, 100), "high": (10000, 12000)}.items()}
SOS4["full"] = butter(4, [20, 12300], btype="band", fs=4 * ml.FSD, output="sos")


def ends(prog):
    """Per channel (L, R): the window end in device samples after tap 0."""
    return [int(np.floor(cc.ONSET[prog][c])) - END_MARGIN for c in (0, 1)]


def hp_inverse(x, corners):
    """x divided by the high-pass cascade (bilinear first-order sections, inverted: causal, from x[0])."""
    for fc in corners:
        b, a = bilinear([1.0, 0.0], [1.0, 2 * np.pi * fc], FS96)
        x = lfilter(a, b, x, axis=0)
    return x


def read_at(y, t):
    """y read at the (fractional) positions t: the Kaiser-windowed sinc of ml_lines.warp."""
    k = np.arange(-mlines.HALF + 1, mlines.HALF + 1)
    i = np.floor(t).astype(int)
    ph = (t - i) * mlines.PHASES
    j = np.minimum(ph.astype(int), mlines.PHASES - 1)
    fr = (ph - j)[:, None]
    w = mlines.WTAB[j] * (1 - fr) + mlines.WTAB[j + 1] * fr
    return np.sum(w * y[i[:, None] + k[None]], 1)


def grid4(x, T, n_before, n, ppm=0.0):
    """x (96 kHz, one channel) on a 4x device grid of a capture whose clock is ppm slow: y[k] = x(T + (k - n_before)
    R4 (1 + ppm 1e-6)), k < n (band-limited)."""
    g, e = 160, 64
    nb, nn = n_before + e, n + 2 * e
    i0 = int(np.floor(T - nb * R4)) - g
    seg = x[i0:i0 + int(np.ceil(nn * R4)) + 2 * g]
    w = np.ones(len(seg)); w[-g // 2:] = np.hanning(g)[g // 2:]          # quiet end, beyond the window and filter
    seg = seg * w
    Tl = T - i0
    J = int(np.floor(Tl / R4))
    s = Tl - J * R4                                                       # advance by s: T lands on grid point J
    N = 1 << int(np.ceil(np.log2(len(seg) + 256)))
    seg = np.fft.irfft(np.fft.rfft(seg, N) * np.exp(2j * np.pi * np.fft.rfftfreq(N) * s), N)[:len(seg)]
    y = resample_poly(seg, 100, 91)[J - nb:J - nb + nn]                   # nominal clock, tap 0 at index nb
    return read_at(y, nb + (np.arange(n) - n_before) * (1 + ppm * 1e-6))


def captures():
    """(prog, key, amp, path) for every capture with a cluster (AMPL > 0) of S, M and L, without the novel clusters."""
    import pandas as pd
    df = pd.read_csv(os.path.join(ROOT, "analysis", "library", "tomkin96", "files.csv"))
    bad = {os.path.abspath(p) for p in df.loc[df["cluster_status"] == "novel_cluster", "path"]}
    out = []
    for prog in "SML":
        for key, g in sorted(cc.groups(prog).items()):
            for amp, path in sorted(g.items()):
                if amp > 0 and os.path.abspath(path) not in bad:
                    out.append((prog, key, amp, path))
    return out


def _data_job(a):
    prog, key, amp, path, (shift, ppm) = a
    x, fs = sf.read(path, always_2d=True)
    am = np.abs(x).max(1)                                 # onset exactly as ml_structure._make
    main = int(np.argmax(am > am.max() * 0.1))
    lo = max(main - int(0.002 * fs), 0)
    onset = lo + int(np.argmax(am[lo:main + 1] > am.max() * 1e-2))
    h = np.load(ml._cache_path(path, OLD_CACHE)).astype(np.float64)
    hb = sosfiltfilt(ml.BAND, h[:4000], axis=0)
    p = cc.tap0(hb)
    A0 = subsample_peak(h[:, 1], int(round(p)))[1]
    PADX = 4096
    p = p + shift                                         # refined tap 0 (align stage), device samples
    T = [onset + PADX + (p - ml.K0) * fs / ml.FSD + ml.LR_US * 1e-6 * fs, onset + PADX + (p - ml.K0) * fs / ml.FSD]
    E = ends(prog)
    n_need = int(max(T) + (max(E) + 50) * fs / ml.FSD) + 1024
    xp = np.vstack([np.zeros((PADX, 2)), x])[:n_need]
    out = dict(A0=A0, p=p, ppm=ppm)
    for var, corners in HP_SETS.items():
        xv = hp_inverse(xp, corners)
        for ch in (0, 1):
            out[f"{var}_{ch}"] = grid4(xv[:, ch], T[ch], 4 * LEAD, 4 * (LEAD + E[ch]), ppm).astype(np.float32)
    return out


def data_stage(workers):
    from concurrent.futures import ProcessPoolExecutor
    os.makedirs(OUT, exist_ok=True)
    caps = captures()
    sp = os.path.join(OUT, "shifts.json")
    shifts = json.load(open(sp)) if os.path.exists(sp) else {}
    with ProcessPoolExecutor(workers, max_tasks_per_child=50) as ex:
        res = list(ex.map(_data_job, [c + (tuple(shifts.get(os.path.relpath(c[3], ROOT), (0.0, PPM0))),) for c in caps],
                          chunksize=8))
    for prog in "SML":
        idx = [i for i, c in enumerate(caps) if c[0] == prog]
        d = dict(key=np.array([caps[i][1] for i in idx]), amp=np.array([caps[i][2] for i in idx]),
                 path=np.array([caps[i][3] for i in idx]), A0=np.array([res[i]["A0"] for i in idx]),
                 p=np.array([res[i]["p"] for i in idx]), ppm=np.array([res[i]["ppm"] for i in idx]))
        for var in HP_SETS:
            for ch in (0, 1):
                d[f"{var}_{ch}"] = np.stack([res[i][f"{var}_{ch}"] for i in idx])
        np.savez_compressed(os.path.join(OUT, f"data_{prog}.npz"), **d)
        print(f"{prog}: {len(idx)} captures, windows L / R {d[VAR + '_0'].shape[1]} / {d[VAR + '_1'].shape[1]} (4x grid)")


# ---------------------------------------------------------------- fit
KHF = 8                         # free output kernel: taps -KHF..KHF (the stored one is -5..5)
AMP_FIT, AMP_SCORE = 30, 60     # captures used in fits (M's AMPL 10 clusters are 19 dB down, some absent) and in scores


def stored_hk():
    """Output kernels per program, taps -KHF..KHF: S none (a unit tap), M and L the stored ones (firstpass.json)."""
    fp = json.load(open(os.path.join(cc.OUT, "firstpass.json")))
    out = {"S": np.r_[np.zeros(KHF), 1.0, np.zeros(KHF)]}
    for p in "ML":
        out[p] = np.pad(np.array(fp[p]["kernel"], float), KHF - cc.KH)
    return out


HK0 = stored_hk()


def taps(prog, ch, hk=None):
    """Device offsets and weights of one channel's cluster, through the program's output kernel (stored by default)."""
    hk = hk or HK0
    h = hk.get((prog, ch), hk[prog])                      # per channel if fitted so
    return [(p_ + k - KHF, s * w) for p_, s in cc.cluster_taps(prog)[ch] for k, w in enumerate(h) if w != 0]


def design(prog, ch, lk, hk=None):
    """U (rows m = -LEAD .. end-1, columns j = -LEAD .. lk-1, device samples): y[m] = sum_j U[m, j] c[j]. Rows stop at
    lk, so every lag is constrained through tap 0 (the output kernel's outer taps past lk are dropped: weights <= 0.04
    times the kernel's tail); returns U and the row count."""
    E = min(ends(prog)[ch], lk)
    U = np.zeros((E + LEAD, lk + LEAD))
    for q, w in taps(prog, ch, hk):
        m = np.arange(max(-LEAD, q - LEAD), min(E, lk + q))      # rows whose lag m - q lies in [-LEAD, lk)
        U[m + LEAD, m - q + LEAD] += w
    return U, E + LEAD


def load_data():
    return {p: dict(np.load(os.path.join(OUT, f"data_{p}.npz"))) for p in "SML"}


NUIS = "quad"                   # per-capture slow term: "quad" {1, t, t^2}; "dc" a DC offset only (through the high-pass inversion)


def nuis_cols(var, n4):
    """(n4, k) nuisance columns on the 4x grid from the window start: a capture's DC offset as it appears in the data
    (constant; after the high-pass inversion its step response sum_k e_k t^k / k!), or a generic quadratic."""
    if NUIS == "quad":
        t = np.linspace(-1, 1, n4)
        return np.stack([np.ones(n4), t, t ** 2], 1)
    if not HP_SETS[var]:
        return np.ones((n4, 1))
    t = np.arange(n4) / (4 * ml.FSD)
    w = 2 * np.pi * np.array(HP_SETS[var])
    e = np.poly(-w)                                       # elementary symmetric polynomials e_0..e_n of the corners
    return np.sum([e[k] * t ** k / np.prod(np.arange(1, k + 1)) for k in range(len(e))], 0)[:, None]


def nuis(nr, var):
    """Orthonormal basis of the nuisance over nr device rows (the 4 polyphase columns share it: smooth)."""
    return np.linalg.qr(nuis_cols(var, 4 * nr)[::4])[0]


def proj(Q, X):
    """X (n, nr, 4) projected off the nuisance basis Q (nr, 3) along its rows."""
    return X - np.einsum("ik,nkr->nir", Q, np.einsum("ik,nir->nkr", Q, X))


def stacked(D, p, var, ch, nr, rows=None):
    seg = D[p][f"{var}_{ch}"] if rows is None else D[p][f"{var}_{ch}"][rows]
    return seg[:, :4 * nr].astype(np.float64).reshape(len(seg), nr, 4)


def solve(D, var=VAR, lk=LK, sel=lambda p, k: True, rounds=ROUNDS, ridge=1e-7, free_hk=FREE_HK, w_s=W_S):
    """Per channel the 4x-grid kernel (4 (lk + LEAD),) from the captures with sel(prog, key) True; gains alternate, and
    with free_hk the M and L output kernels too ("shared" by the channels or "perch"; scaled to the stored ones' sums).
    w_s weights program S's captures. Every capture's window is projected off its nuisance {1, t, t^2}, model and data
    alike. Returns the kernels, the gains and the output kernels."""
    hk = {p: v.copy() for p, v in HK0.items()}
    A = {p: D[p]["A0"].copy() for p in "SML"}
    use = {p: np.array([sel(p, k) for k in D[p]["key"]]) & (D[p]["amp"] >= AMP_FIT) for p in "SML"}
    rows = {(p, ch): design(p, ch, lk)[1] for p in "SML" for ch in (0, 1)}
    Qs = {pc: nuis(nr, var) for pc, nr in rows.items()}
    Y = {pc: proj(Qs[pc], stacked(D, pc[0], var, pc[1], rows[pc])) for pc in rows}
    for it in range(rounds):
        PU = {}
        for (p, ch), Q in Qs.items():
            U = design(p, ch, lk, hk)[0]
            PU[(p, ch)] = U - Q @ (Q.T @ U)
        c = []
        for ch in (0, 1):
            G = 0
            v = 0
            for p in "SML":
                if not use[p].any():
                    continue
                a = A[p][use[p]]
                wp = w_s if p == "S" else 1.0
                G = G + wp * (a @ a) * (PU[(p, ch)].T @ PU[(p, ch)])
                v = v + wp * PU[(p, ch)].T @ np.tensordot(a, Y[(p, ch)][use[p]], 1)
            G = G + ridge * np.trace(G) / len(G) * np.eye(len(G))
            c.append(np.linalg.solve(G, v))                       # (lk + LEAD, 4)
        for p in "SML":                                           # gains: one per capture, both channels
            num = den = 0
            for ch in (0, 1):
                m = PU[(p, ch)] @ c[ch]
                num = num + np.einsum("nij,ij->n", Y[(p, ch)], m)
                den = den + np.sum(m ** 2)
            A[p] = num / den
        if free_hk and it < rounds - 1:
            for p in "ML":
                if not use[p].any():
                    continue
                a = A[p][use[p]]
                M, r = [0, 0], [0, 0]
                for ch in (0, 1):
                    nr = rows[(p, ch)]
                    X = np.stack([render(c[ch].reshape(-1), p, ch, nr, {p: np.eye(2 * KHF + 1)[k]}).reshape(nr, 4)
                                  for k in range(2 * KHF + 1)])                     # (17, nr, 4)
                    X = proj(Qs[(p, ch)], X).reshape(2 * KHF + 1, -1)
                    M[ch] = (a @ a) * (X @ X.T)
                    r[ch] = X @ np.tensordot(a, Y[(p, ch)][use[p]], 1).reshape(-1)
                if free_hk == "perch":
                    for ch in (0, 1):
                        h = np.linalg.solve(M[ch], r[ch])
                        hk[(p, ch)] = h * HK0[p].sum() / h.sum()
                else:
                    h = np.linalg.solve(M[0] + M[1], r[0] + r[1])
                    hk[p] = h * HK0[p].sum() / h.sum()
    return [x.reshape(-1) for x in c], A, hk                      # interleaved: 4x-grid index 4 (j + LEAD) + r


def old_kernel4(lk):
    g = Kernel.load(OLD).g
    g4 = resample_poly(g, 100, 91)
    out = np.zeros(4 * (lk + LEAD))
    out[4 * LEAD:4 * LEAD + min(len(g4), 4 * lk)] = g4[:4 * lk]
    return [out, out]


def render(c4, prog, ch, nrow, hk=None):
    """Model on the 4x grid (unit gain): sum_q w_q c(t - 4 q)."""
    y = np.zeros(4 * nrow + 4 * 4000)
    for q, w in taps(prog, ch, hk):
        a = 4 * q
        if a >= 0:
            y[a:a + len(c4)] += w * c4
        else:
            y[:len(c4) + a] += w * c4[-a:]
    return y[:4 * nrow]


def deriv(x):
    N = 1 << int(np.ceil(np.log2(len(x) + 512)))
    return np.fft.irfft(np.fft.rfft(x, N) * 2j * np.pi * np.fft.rfftfreq(N), N)[:len(x)]


def stack_cols(Xs):
    """Both channels' columns as one system: the first 3 (gain, shift, clock) shared, the nuisances per channel."""
    n0, k = Xs[0].shape
    X = np.zeros((n0 + len(Xs[1]), 3 + 2 * (k - 3)))
    X[:n0, :3], X[n0:, :3] = Xs[0][:, :3], Xs[1][:, :3]
    X[:n0, 3:k], X[n0:, k:] = Xs[0][:, 3:], Xs[1][:, 3:]
    return X


def fit_captures(D, var, c4, p, rows=None, hk=None):
    """Per capture (both channels): gain a, shift and clock (linearised: y = a m + b m' + e tau m', tau from tap 0) and
    the nuisance (nuis_cols) per channel, by least squares over the whole window. Returns the data, the fitted models,
    the shifts (4x-grid samples) and the clock errors (ppm)."""
    Ys, Xs = [], []
    for ch in (0, 1):
        nr = ends(p)[ch] + LEAD
        Ys.append(stacked(D, p, var, ch, nr, rows).reshape(-1, 4 * nr))
        m = render(c4[ch], p, ch, nr, hk)
        dm = deriv(m)
        Xs.append(np.column_stack([m, dm, (np.arange(4 * nr) - 4 * LEAD) * dm, nuis_cols(var, 4 * nr)]))
    X = stack_cols(Xs)
    coef = np.linalg.lstsq(X, np.hstack(Ys).T, rcond=None)[0]
    fit = (X @ coef).T
    n0 = Ys[0].shape[1]
    return Ys, [fit[:, :n0], fit[:, n0:]], -coef[1] / coef[0], -coef[2] / coef[0] * 1e6


def score(D, var, c4, sel, hk=None):
    """Median nulls per program and channel over the captures with sel True and AMPL >= AMP_SCORE (gain, shift and
    nuisance fitted per capture, on the whole window), per band."""
    out = {}
    for p in "SML":
        keep = np.where([sel(p, k) and a >= AMP_SCORE for k, a in zip(D[p]["key"], D[p]["amp"])])[0]
        if not len(keep):
            continue
        Ys, Fs, _, _ = fit_captures(D, var, c4, p, keep, hk)
        out[p] = dict(n=int(len(keep)), **band_nulls(Ys, Fs))
    return out


def band_nulls(Ys, Fs):
    res = {}
    for b, sos in SOS4.items():
        row = []
        for ch in (0, 1):
            z = np.zeros((len(Ys[ch]), 400))
            y = sosfiltfilt(sos, np.hstack([z, Ys[ch], z]), axis=1)
            r = sosfiltfilt(sos, np.hstack([z, Ys[ch] - Fs[ch], z]), axis=1)
            row.append(float(np.median(10 * np.log10(np.sum(r ** 2, 1) / np.sum(y ** 2, 1)))))
        res[b] = row
    return res


def repeat_floor(D, var, sel):
    """The same nulls with another capture as the model: per group, AMPL 60 predicted by the highest AMPL (S 90, M 99;
    L has no second cluster capture), same freedoms. Both captures carry noise, so a noise-free model sits up to 3 dB
    deeper."""
    out = {}
    for p in "SM":
        Ys, Fs = [[], []], [[], []]
        for k in sorted(set(D[p]["key"])):
            if not sel(p, k):
                continue
            idx = np.where(D[p]["key"] == k)[0]
            lo = [i for i in idx if D[p]["amp"][i] == 60]
            hi = idx[np.argmax(D[p]["amp"][idx])]
            for i in lo:
                ys, xs = [], []
                for ch in (0, 1):
                    nr = ends(p)[ch] + LEAD
                    y = D[p][f"{var}_{ch}"][i, :4 * nr].astype(np.float64)
                    m = D[p][f"{var}_{ch}"][hi, :4 * nr].astype(np.float64)
                    ys.append(y)
                    dm = deriv(m)
                    xs.append(np.column_stack([m, dm, (np.arange(4 * nr) - 4 * LEAD) * dm, nuis_cols(var, 4 * nr)]))
                n0 = len(ys[0])
                X = stack_cols(xs)
                f = X @ np.linalg.lstsq(X, np.r_[ys[0], ys[1]], rcond=None)[0]
                for ch, (y, fc) in enumerate(((ys[0], f[:n0]), (ys[1], f[n0:]))):
                    Ys[ch].append(y)
                    Fs[ch].append(fc)
        out[p] = dict(n=len(Ys[0]), **band_nulls([np.array(v) for v in Ys], [np.array(v) for v in Fs]))
    return out


def align_stage(rounds):
    """Refine every capture's tap-0 time and clock: fit the kernel on all captures, measure each capture's shift and
    clock error, rebuild the windows with them (shifts.json: [device samples, ppm], cumulative)."""
    sp = os.path.join(OUT, "shifts.json")
    for it in range(rounds):
        D = load_data()
        c, _, hk = solve(D)
        shifts = json.load(open(sp)) if os.path.exists(sp) else {}
        for p in "SML":
            _, _, d4, dppm = fit_captures(D, VAR, c, p, hk=hk)
            ok = D[p]["amp"] >= AMP_FIT
            for path, d, e, o in zip(D[p]["path"], d4, dppm, ok):
                key = os.path.relpath(str(path), ROOT)
                s0, e0 = shifts.get(key, (0.0, PPM0))
                shifts[key] = (s0 + float(d) / 4, e0 + float(e)) if o else (s0, e0)
            print(f"round {it + 1} {p}: shift (device samples) median {np.median(d4[ok]) / 4:+.4f}, p5..p95 "
                  f"{np.percentile(d4[ok], 5) / 4:+.4f}..{np.percentile(d4[ok], 95) / 4:+.4f}; clock correction (ppm) median "
                  f"{np.median(dppm[ok]):+.2f}, p5..p95 {np.percentile(dppm[ok], 5):+.2f}..{np.percentile(dppm[ok], 95):+.2f}",
                  flush=True)
        json.dump(shifts, open(sp, "w"), indent=0)
        data_stage(4)


def fir_numbers(lk):
    return 2 * 4 * (lk + LEAD)


def tail_db(c4):
    """Per channel, the FIR's energy per lag range (dB re its total)."""
    out = []
    for c in c4:
        e = c ** 2
        lag = np.arange(len(e)) // 4 - LEAD
        out.append({f"{a}-{b}": round(float(10 * np.log10(e[(lag >= a) & (lag < b)].sum() / e.sum() + 1e-30)), 1)
                    for a, b in ((-LEAD, 0), (0, 50), (50, 150), (150, 300), (300, 600), (600, 1200), (1200, 2300))
                    if a < len(e) // 4 - LEAD})
    return out


def fit_stage():
    global NUIS
    D = load_data()
    keys = {p: sorted(set(D[p]["key"])) for p in "SML"}
    even = lambda p, k: keys[p].index(k) % 2 == 0
    odd = lambda p, k: not even(p, k)
    res = {}
    lines = ["# Chain kernel: nulls over the cluster-only windows, median over AMPL >= 60 captures, L / R, dB "
             "(per capture a gain, a shift, a clock and the slow term refitted). Fitted on even settings, odd scored, "
             "unless stated."]

    def show(name, sc, params, extra=None):
        res[name] = dict(sc, params=params, **(extra or {}))
        lines.append(f"\n## {name} ({params} numbers)")
        for p, s_ in sc.items():
            lines.append(f"- {p} ({s_['n']}): " + "; ".join(f"{b} {v[0]:.1f} / {v[1]:.1f}"
                                                         for b, v in s_.items() if b in SOS4))
        if extra and "tail_db" in extra:
            lines.append(f"- FIR energy by lag, L: {extra['tail_db'][0]}; R: {extra['tail_db'][1]}")
        print("\n".join(lines[-1 - len(sc) - (1 if extra and 'tail_db' in extra else 0):]), flush=True)
    show("repeat capture (AMPL 60 from the highest AMPL of its group)", repeat_floor(D, VAR, odd), 0)
    show("library kernel (kernel_96000, 5.2 ms, R only)", score(D, "none", old_kernel4(LK), odd), 0)
    hk_n = 2 * (2 * KHF + 1)
    runs = [("chosen: schematic high-pass, FIR 500, output kernels free, S weight 3", {}),
            ("output kernels stored", dict(free_hk=False)),
            ("S weight 1", dict(w_s=1.0)),
            ("no high-pass (FIR alone)", dict(var="none")),
            ("high-pass schematic + 1.8 Hz (analog_results.md §3.5)", dict(var="schematic+1.8")),
            ("high-pass half the schematic + 1.8 Hz", dict(var="half+1.8")),
            ("FIR 300", dict(lk=300)), ("FIR 800", dict(lk=800)), ("FIR 1200", dict(lk=1200)), ("FIR 2300", dict(lk=2300)),
            ("FIR 2300, no high-pass", dict(lk=2300, var="none"))]
    for name, kw in runs:
        c, _, hk = solve(D, sel=even, **kw)
        lk = kw.get("lk", LK)
        show(name, score(D, kw.get("var", VAR), c, odd, hk), fir_numbers(lk) + (hk_n if kw.get("free_hk", FREE_HK) else 0),
             dict(tail_db=tail_db(c), fir_sum_over_peak=[round(float(x.sum() / np.abs(x).max()), 2) for x in c]))
    NUIS = "dc"
    c, _, hk = solve(D, sel=even)
    show("slow term: a DC offset only (through the high-pass inversion)", score(D, VAR, c, odd, hk), fir_numbers(LK) + hk_n,
         dict(tail_db=tail_db(c)))
    NUIS = "quad"
    c, _, _ = solve(D, sel=lambda p, k: p != "S", free_hk=False)
    show("S held out: fitted on M and L only (stored output kernels), every S capture scored",
         score(D, VAR, c, lambda p, k: p == "S"), fir_numbers(LK))
    json.dump(res, open(os.path.join(OUT, "fit.json"), "w"), indent=1)
    open(os.path.join(OUT, "fit.log"), "w").write("\n".join(lines) + "\n")


def write_stage(name):
    """Fit on all captures, convert to 96 kHz, pin R's mid-band gain and delay to the old kernel (the cache's level
    and time convention; L keeps its own relative to R), save."""
    D = load_data()
    c4, A, hk = solve(D)
    lead96 = int(round(4 * LEAD * R4))
    g = np.stack([resample_poly(c, 91, 100) for c in c4], 1)
    hp = HP_SETS[VAR]
    N = 1 << 16
    f = np.fft.rfftfreq(N, 1 / FS96)
    Gn = Kernel(g, FS96, lead=lead96, hp_hz=hp).response(N)
    Go = np.fft.rfft(Kernel.load(OLD).g, N)
    m = (f > 300) & (f < 3000)
    H = Gn[m, 1] / Go[m]
    ph = np.polyfit(f[m], np.unwrap(np.angle(H)), 1)
    gain = np.exp(np.mean(np.log(np.abs(H))))
    tau = -ph[0] / (2 * np.pi) * FS96                                    # samples the new R kernel is late by
    Nk = 1 << int(np.ceil(np.log2(len(g) + 512)))
    shift = np.exp(2j * np.pi * np.fft.rfftfreq(Nk) * tau)[:, None]
    g = np.fft.irfft(np.fft.rfft(g, Nk, axis=0) * shift, Nk, axis=0)[:len(g)] / gain
    src = (f"chain_kernel.py: per-channel chain kernel, FIR {LK} + {LEAD} device samples x high-pass {hp} Hz, fitted on "
           f"the cluster-only windows of {sum(int((D[p]['amp'] >= AMP_FIT).sum()) for p in 'SML')} S/M/L captures; "
           f"R pinned to kernel_96000 at 300-3000 Hz (gain {20 * np.log10(gain):+.3f} dB, delay {tau:+.4f} samples removed)")
    k = Kernel(g, FS96, source=src, lead=lead96, hp_hz=hp)
    k.save(os.path.join(KDIR, name + ".npz"))
    Gn = k.response(N)
    rel = Gn / Go[:, None]
    LR = Gn[m, 0] / Gn[m, 1]
    pl = np.polyfit(f[m], np.unwrap(np.angle(LR)), 1)
    info = dict(name=name, source=src, hp_hz=hp, fir_device_samples=LK, lead96=lead96, numbers=fir_numbers(LK),
                pinned_gain_db=float(20 * np.log10(gain)), pinned_delay_samples=float(tau),
                L_vs_R_300_3000=dict(gain_db=float(np.mean(20 * np.log10(np.abs(LR)))),
                                     delay_samples=float(-pl[0] / (2 * np.pi) * FS96)),
                response_vs_old_dB_deg={str(fr): [[round(float(20 * np.log10(abs(rel[i, ch]))), 3),
                                                   round(float(np.degrees(np.angle(rel[i, ch]))), 2)] for ch in (0, 1)]
                                        for fr in (10, 25, 50, 100, 160, 280, 400, 1000, 4000, 8000, 10000, 11000, 11500, 12000)
                                        for i in [int(round(fr / f[1]))]},
                tail_db=tail_db(c4), output_kernels={p: [round(float(v), 5) for v in hk[p]] for p in "ML"},
                stored_output_kernels={p: [round(float(v), 5) for v in HK0[p]] for p in "ML"})
    json.dump(info, open(os.path.join(OUT, f"{name}.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", required=True, choices=["data", "align", "fit", "write"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--rounds", type=int, default=2, help="align: rounds")
    ap.add_argument("--name", default="chain2", help="write: kernel name")
    a = ap.parse_args()
    if a.stage == "data":
        data_stage(a.workers)
    elif a.stage == "align":
        align_stage(a.rounds)
    elif a.stage == "fit":
        fit_stage()
    else:
        write_stage(a.name)


if __name__ == "__main__":
    main()
