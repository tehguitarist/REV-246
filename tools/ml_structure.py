"""Programs M and L from the Tomkin 96 kHz library (NEXT_STEPS Q3): cluster patterns, reverb onset,
loop recursion and time variance, with program S as the control. Results: docs/ml_structure.md.

    .venv/bin/python tools/ml_structure.py [--stage all|clusters|loops|lags|pairs] [--out analysis/ml_structure]

Writes <out>/summary.json (and caches device-rate IRs as <out>/cache/*.npy, not committed). EMT_KERNEL=<name> in the
environment deconvolves with analysis/library/tomkin96/kernels/<name>.npz instead (cache_<name>/; chain2 is the
per-channel chain kernel of tools/chain_kernel.py, docs/ml_structure.md §12).

Methods (all on device-rate IRs: deconvolved with the library kernel, resampled to 19.2 MHz/728):

  clusters  AMPL pairs (same LF/RT/HF, different AMPL), aligned on the reverb (0.3-0.8 s) and
            subtracted; the four largest peaks per channel in the first 2,000 samples are the cluster.
            Reverb onset = first sample of the low-AMPL file above -30 dB of its 0-0.3 s peak (cluster masked).
  loops     Tail predictor: y[n] = sum_k s_k y[n-k] (k = 1..300: spectral colour and band-limit)
            + sum_d sum_j w_dj y[n-d-j] (j = -3..3) over a candidate lag set. Fitted on one window,
            null reported on the next (held out). Once the input has ended, any LTI state-space
            output obeys the characteristic polynomial, so the right lag set predicts the tail to the
            noise; random lag sets of the same size are the control.
  lags      Sub-sample position of each lag (peak of the low-passed lag FIR) across settings;
            S is the calibration (its lags are known integers).
  pairs     AMPL pairs as repeat captures of the reverb: AMPL changes only the cluster, so an LTI
            unit gives identical tails. Null per time window and per octave.
  floors    Capture-domain noise floor against the CLA IR: in-band level per 25 ms frame relative to the
            50-300 ms reverb; how deep the natural decay holds before a floor or Tomkin's fade.
  clock     S AMPL pairs in the capture domain with a per-capture clock ratio fitted (on R, 0.05-0.6 s):
            the repeat-capture null a perfect model could reach against one library capture.
"""
import argparse, collections, glob, itertools, json, os, re, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import soundfile as sf
from scipy.signal import butter, find_peaks, firwin, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
from core import PAD, Kernel, deconvolve, fs_dev, load_config, subsample_peak, to_device_grid  # noqa: E402

CFG = load_config()
FSD = fs_dev(CFG)
LIB = os.path.join(ROOT, "irs", "EMT 246 captures", "96kHz")
# EMT_KERNEL=<name> in the environment deconvolves with analysis/library/tomkin96/kernels/<name>.npz instead of the
# library kernel (kernel_96000); every cache and result file that depends on it gets the suffix "_<name>" (tagged()).
KTAG = os.environ.get("EMT_KERNEL", "")
KERNEL = os.path.join(ROOT, "analysis", "library", "tomkin96", "kernels", (KTAG or "kernel_96000") + ".npz")


def tagged(path):
    """path with the kernel's suffix (unchanged for the library kernel)."""
    return f"{path}_{KTAG}" if KTAG else path


CACHE = tagged(os.path.join(ROOT, "analysis", "ml_structure", "cache"))
PROG = {"Short RT": "S", "Med RT": "M", "Long RT": "L"}
LR_US = 19.19                  # library L/R offset (library_notes.md §3)
K0 = 50                        # device index of the onset in the cached IRs
SHORT, HW = 300, 3             # predictor: short block length, half-width of each lag FIR
BAND = butter(4, [100, 10000], btype="band", fs=FSD, output="sos")
OCT = [125, 250, 500, 1000, 2000, 4000, 8000]
WIN = [(0.05, 0.15), (0.15, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 2.0), (2.0, 4.0)]

S_LOOPS = [1321, 1389, 1461, 1587]
ML_TAPS = {"M": [5416, 5694, 5812, 6090, 11506], "L": [10826, 11384, 11618, 12176, 23002]}
# loop test files and windows (fit [ta, tb), held out [tb, tc)): flattest filters with a long tail
LOOP_TEST = {"S": (dict(lf=1.0, hf=1.0, amp=0, rt=2.5), 0.4, 1.0, 1.4),
             "M": (dict(lf=1.2, hf=0.9, amp=10, rt=5.0), 1.0, 2.3, 3.3),
             "L": (dict(lf=1.2, hf=1.0, amp=0, rt=20.0), 2.0, 4.5, 6.5)}
# lag precision: (file filter, windows). S is the calibration, so flat loop filters only (LF 1.0, HF 1.0)
LAG_TEST = {"S": (lambda l: l["lf"] == 1.0 and l["hf"] == 1.0 and l["rt"] >= 1.8, 0.35, 0.9, 1.2),
            "M": (lambda l: l["hf"] >= 0.9, 0.9, 2.0, 2.6),
            "L": (lambda l: l["hf"] >= 1.0, 1.8, 4.0, 5.0)}
LATE_TEST = {"L": (8.0, 10.5, 12.5)}   # a later window of the same file: does the plateau stay put?


# ---------------------------------------------------------------- data
def labels(path):
    name = os.path.basename(path).replace("_dc.wav", "")
    lp, delay, amp, lf, rt, hf = (int(v) for v in name.split("_"))
    return dict(prog=PROG[os.path.basename(os.path.dirname(path))], amp=amp, lf=lf / 10, rt=rt / 10, hf=hf / 10)


def files(prog):
    return [f for f in sorted(glob.glob(os.path.join(LIB, "*", "*.wav")))
            if os.path.basename(os.path.dirname(f)) in PROG and labels(f)["prog"] == prog]


def pick(prog, **want):
    return [f for f in files(prog) if all(labels(f)[k] == v for k, v in want.items())]


def _cache_path(path, cache):
    return os.path.join(cache, os.path.basename(os.path.dirname(path)).replace(" ", "") + "__" + os.path.basename(path)[:-4] + ".npy")


def _make(args):
    path, cache = args
    out = _cache_path(path, cache)
    if os.path.exists(out):
        return
    x, fs = sf.read(path, always_2d=True)
    a = np.abs(x).max(1)
    main = int(np.argmax(a > a.max() * 0.1))
    lo = max(main - int(0.002 * fs), 0)
    onset = lo + int(np.argmax(a[lo:main + 1] > a.max() * 1e-2))
    d = deconvolve(np.vstack([np.zeros((PAD, 2)), x]), fs, Kernel.load(KERNEL))
    anchor = float(onset + PAD)
    k0 = int(np.floor(anchor / fs * FSD))
    h = np.stack([to_device_grid(d[:, 0], fs, anchor + LR_US * 1e-6 * fs, k0, CFG),
                  to_device_grid(d[:, 1], fs, anchor, k0, CFG)], 1)
    np.save(out, h[max(k0 - K0, 0):].astype(np.float32))


class IRs:
    """Device-rate IRs (n, 2) with the onset at index K0, cached per file."""

    def __init__(self, cache):
        self.cache = cache
        os.makedirs(cache, exist_ok=True)

    def build(self, paths, workers=10):
        todo = [(p, self.cache) for p in paths if not os.path.exists(_cache_path(p, self.cache))]
        if todo:
            with ProcessPoolExecutor(workers) as ex:
                list(ex.map(_make, todo, chunksize=4))

    def __call__(self, path, band=False):
        h = np.load(_cache_path(path, self.cache)).astype(np.float64)
        return sosfiltfilt(BAND, h, axis=0) if band else h


def ampl_groups(prog):
    g = collections.defaultdict(dict)
    for f in files(prog):
        l = labels(f)
        g[(l["lf"], l["rt"], l["hf"])][l["amp"]] = f
    return {k: v for k, v in sorted(g.items(), key=lambda kv: (-kv[0][1], kv[0])) if len(v) > 1}


def align(a, b, s0, s1, maxlag=2000):
    """Integer lag of b against a, from the cross-correlation over a[s0:s1]; returns b shifted onto a."""
    n = 1 << int(np.ceil(np.log2(len(b) + s1)))
    c = np.fft.irfft(np.fft.rfft(b[:s1 + maxlag, 1], n) * np.conj(np.fft.rfft(np.r_[np.zeros(s0), a[s0:s1, 1]], n)), n)
    lags = np.r_[np.arange(maxlag), np.arange(-maxlag, 0)]
    lag = int(lags[np.argmax(np.r_[c[:maxlag], c[-maxlag:]])])
    b = b[lag:] if lag >= 0 else np.vstack([np.zeros((-lag, 2)), b])
    m = min(len(a), len(b))
    return a[:m], b[:m]


# ---------------------------------------------------------------- clusters and onsets
def cluster_stage(ir, n_pairs=40):
    out = {}
    for prog in "SML":
        offs, vals, onsets = [], [], []
        for key, d in list(ampl_groups(prog).items())[:n_pairs]:
            amps = sorted(d)
            a, b = align(ir(d[amps[-1]])[:int(1.0 * FSD)], ir(d[amps[0]])[:int(1.0 * FSD)], int(0.3 * FSD), int(0.8 * FSD))
            a, b = a[:int(0.4 * FSD)], b[:int(0.4 * FSD)]
            g = [np.dot(a[2000:, c], b[2000:, c]) / np.dot(b[2000:, c], b[2000:, c]) for c in range(2)]
            diff = a - np.array(g) * b
            o, v, on = [], [], []
            for c in range(2):
                e = np.abs(diff[:2000, c])
                pk, _ = find_peaks(e, height=0.3 * e.max(), distance=20)
                t = [subsample_peak(diff[:, c], p) for p in sorted(sorted(pk, key=lambda p: -e[p])[:4])]
                o.append([p - t[0][0] for p, _ in t]); v.append([x for _, x in t])
                x = np.abs(b[:, c]).copy()
                for p, _ in t:
                    x[max(int(round(p)) - 25, 0):int(round(p)) + 26] = 0
                on.append(int(np.argmax(x > 10 ** (-30 / 20) * x.max())) - t[0][0])
            offs.append(o); vals.append(v); onsets.append(on)
        O, V, N = np.array(offs), np.array(vals), np.array(onsets)
        out[prog] = {"pairs": len(O),
                     "offsets": {ch: np.round(O[:, c].mean(0), 3).tolist() for c, ch in enumerate("LR")},
                     "offset_std": {ch: np.round(O[:, c].std(0), 3).tolist() for c, ch in enumerate("LR")},
                     "signs": {ch: np.sign(np.median(V[:, c], 0)).astype(int).tolist() for c, ch in enumerate("LR")},
                     "tap_spread_db": float(20 * np.log10(np.abs(np.median(V, 0)).max() / np.abs(np.median(V, 0)).min())),
                     "reverb_onset_after_tap0": {ch: float(np.median(N[:, c])) for c, ch in enumerate("LR")},
                     "reverb_onset_iqr": {ch: np.percentile(N[:, c], [25, 75]).tolist() for c, ch in enumerate("LR")}}
        print(prog, json.dumps(out[prog]))
    return out


# ---------------------------------------------------------------- tail predictor
def _cols(y, n0, n1, lags, short=True):
    c = [y[n0 - k:n1 - k] for k in range(1, SHORT + 1)] if short else []
    c += [y[n0 - d - j:n1 - d - j] for d in lags for j in range(-HW, HW + 1)]
    return np.asarray(c, np.float32).T


def predict(y, lags, ta, tb, tc, short=True):
    """Fit on [ta, tb), return (held-out null dB, in-sample null dB, weights, held-out residual, held-out target)."""
    n0 = max(int(ta * FSD), max(list(lags) + [SHORT]) + HW + 1)
    n1, n2 = int(tb * FSD), min(int(tc * FSD), len(y))
    A = _cols(y, n0, n1, lags, short)
    G = (A.T @ A).astype(np.float64)
    G += 1e-6 * np.trace(G) / len(G) * np.eye(len(G))
    w = np.linalg.solve(G, (A.T @ y[n0:n1]).astype(np.float64))
    ins = 10 * np.log10(np.sum((y[n0:n1] - A @ w) ** 2) / np.sum(y[n0:n1] ** 2))
    r = y[n1:n2] - _cols(y, n1, n2, lags, short) @ w
    return 10 * np.log10(np.sum(r ** 2) / np.sum(y[n1:n2] ** 2)), ins, w, r, y[n1:n2]


def subset_sums(D, kmax=None):
    return sorted({sum(c) for k in range(1, (kmax or len(D)) + 1) for c in itertools.combinations(D, k)})


def octave_ratio(r, y):
    out = []
    for fc in OCT:
        s = butter(4, [fc / 2 ** 0.5, min(fc * 2 ** 0.5, 12000)], btype="band", fs=FSD, output="sos")
        out.append(round(float(10 * np.log10(np.sum(sosfiltfilt(s, r) ** 2) / np.sum(sosfiltfilt(s, y) ** 2))), 1))
    return out


def loop_stage(ir):
    rng = np.random.default_rng(1)
    out = {}
    for prog, (want, ta, tb, tc) in LOOP_TEST.items():
        f = pick(prog, **want)[0]
        h = ir(f, band=True)
        if prog == "S":
            models = {"short block only": [], "4 single lags": S_LOOPS, "single lags + pair sums": subset_sums(S_LOOPS, 2),
                      "4 parallel loops (15 subset sums)": subset_sums(S_LOOPS)}
            base = subset_sums(S_LOOPS)
        else:
            T = ML_TAPS[prog]
            models = {"short block only": [], "4 mid taps": T[:4], "5 taps": T,
                      "4 mid taps as parallel loops (15 subset sums)": subset_sums(T[:4])}
            base = T
        for i in range(3):
            models[f"control {i}: random lags, same count"] = sorted(int(v) for v in rng.uniform(0.92, 1.08, len(base)) * np.array(base))
        res = {}
        for name, lags in models.items():
            r = [predict(h[:, c], lags, ta, tb, tc) for c in range(2)]
            res[name] = {"params": SHORT + len(lags) * (2 * HW + 1), "held_out_db": [round(float(x[0]), 2) for x in r],
                         "in_sample_db": [round(float(x[1]), 2) for x in r]}
        r = [predict(h[:, c], base, ta, tb, tc) for c in range(2)]
        out[prog] = {"file": os.path.basename(f), "labels": labels(f), "fit_s": [ta, tb], "held_out_s": [tb, tc],
                     "models": res, "residual_per_octave_db": {ch: octave_ratio(x[3], x[4]) for ch, x in zip("LR", r)}}
        if prog in LATE_TEST:
            t = LATE_TEST[prog]
            out[prog]["later_window"] = {"fit_s": t[:2], "held_out_s": t[1:],
                                         "held_out_db": [round(float(predict(h[:, c], base, *t)[0]), 2) for c in range(2)]}
        print(prog, json.dumps(out[prog]))
    return out


def lag_stage(ir, n=12):
    lp = firwin(15, 9000, fs=FSD)
    out = {}
    for prog, (ok, ta, tb, tc) in LAG_TEST.items():
        lags = subset_sums(S_LOOPS) if prog == "S" else ML_TAPS[prog]
        show = [lags.index(d) for d in (S_LOOPS + [sum(S_LOOPS)] if prog == "S" else lags)]
        fl = [f for f in files(prog) if ok(labels(f)) and labels(f)["rt"] >= tc]
        fl = [f for f in fl if len(ir(f)) > tc * FSD * 1.05]
        fl = fl[::max(len(fl) // n, 1)][:n]
        pos = []
        for f in fl:
            h = ir(f, band=True)
            p = []
            for c in range(2):
                W = predict(h[:, c], lags, ta, tb, tc)[2][SHORT:].reshape(len(lags), -1)
                for j in show:
                    s = np.convolve(np.r_[np.zeros(20), W[j], np.zeros(20)], lp, "same")
                    i = int(np.argmax(np.abs(s)))
                    a, b, cc = np.abs(s[i - 1:i + 2])
                    p.append(i - 20 - HW + 0.5 * (a - cc) / (a - 2 * b + cc))
            pos.append(p)
        P = np.array(pos).reshape(len(fl), 2, len(show))
        out[prog] = {"files": len(fl), "lags": [lags[j] for j in show],
                     "offset_mean": {ch: np.round(P[:, c].mean(0), 3).tolist() for c, ch in enumerate("LR")},
                     "offset_std": {ch: np.round(P[:, c].std(0), 3).tolist() for c, ch in enumerate("LR")}}
        print(prog, json.dumps(out[prog]))
    return out


# ---------------------------------------------------------------- AMPL pairs as repeat captures
def _shift(x, s):
    N = 1 << int(np.ceil(np.log2(len(x) + 64)))
    return np.fft.irfft(np.fft.rfft(x, N) * np.exp(-2j * np.pi * np.fft.rfftfreq(N) * s), N)[:len(x)]


def pair_nulls(ir, fa, fb):
    a, b = align(ir(fa, band=True), ir(fb, band=True), int(0.3 * FSD), int(0.8 * FSD))
    s0, s1 = int(0.3 * FSD), int(0.8 * FSD)

    def cost(s):
        x = _shift(b[:, 1], s)[s0:s1]
        return np.sum((a[s0:s1, 1] - x * np.dot(a[s0:s1, 1], x) / np.dot(x, x)) ** 2)
    s = min(np.linspace(-0.5, 0.5, 21), key=cost)
    b = np.stack([_shift(b[:, c], s) for c in range(2)], 1)
    win, octs = [], {}
    for t0, t1 in WIN:
        i0, i1 = int(t0 * FSD), int(t1 * FSD)
        if i1 + 2000 > len(a):
            win.append(None); continue
        nl = []
        for c in range(2):
            g = np.dot(a[i0:i1, c], b[i0:i1, c]) / np.dot(b[i0:i1, c], b[i0:i1, c])
            nl.append(float(10 * np.log10(np.sum((a[i0:i1, c] - g * b[i0:i1, c]) ** 2) / np.sum(a[i0:i1, c] ** 2))))
        win.append(nl)
        if 0.15 <= t0 < 2.0:
            row = []
            for fc in OCT:
                sos = butter(4, [fc / 2 ** 0.5, min(fc * 2 ** 0.5, 12000)], btype="band", fs=FSD, output="sos")
                x, y = (sosfiltfilt(sos, z[i0 - 2000:i1 + 2000, 1])[2000:-2000] for z in (a, b))
                row.append(float(10 * np.log10(np.sum((x - y * np.dot(x, y) / np.dot(y, y)) ** 2) / np.sum(x ** 2))))
            octs[(t0, t1)] = row
    return win, octs


def pairs_stage(ir, n_pairs=40):
    out = {}
    for prog in "SML":
        W, O = [], collections.defaultdict(list)
        for key, d in list(ampl_groups(prog).items())[:n_pairs]:
            amps = sorted(d)
            w, o = pair_nulls(ir, d[amps[-1]], d[amps[-2]])
            W.append(w)
            for k, v in o.items():
                O[k].append(v)
        res = {"pairs": len(W), "windows": {}, "octaves_R": {}}
        for j, (t0, t1) in enumerate(WIN):
            v = np.array([w[j] for w in W if w[j] is not None])
            if len(v):
                res["windows"][f"{t0}-{t1}"] = {"n": len(v), "median_null_db_LR": np.round(np.median(v, 0), 1).tolist()}
        for (t0, t1), v in O.items():
            res["octaves_R"][f"{t0}-{t1}"] = dict(zip(map(str, OCT), np.round(np.median(v, 0), 1).tolist()))
        out[prog] = res
        print(prog, json.dumps(res))
    # M has three AMPL values: per-capture randomness gives the same null for every pairing
    trip = collections.defaultdict(list)
    for key, d in list(ampl_groups("M").items())[:25]:
        if len(d) == 3:
            for hi, lo in [(60, 10), (99, 60), (99, 10)]:
                trip[f"{hi}-{lo}"].append(pair_nulls(ir, d[hi], d[lo])[0])
    out["M_ampl_pairings"] = {k: {f"{t0}-{t1}": round(float(np.median([w[j][1] for w in v if w[j] is not None])), 1)
                                  for j, (t0, t1) in enumerate(WIN[:5])} for k, v in trip.items()}
    print("M pairings (R)", json.dumps(out["M_ampl_pairings"]))
    return out


# ---------------------------------------------------------------- noise floors (capture domain, 96 kHz WAVs)
CLA = os.path.join(ROOT, "irs", "eugtone_EMT246_CLA_shortRT.wav")
FLOOR_BAND = (20, 11300)
FRAME = 0.025


def _floor_file(path):
    """In-band level per 25 ms frame relative to the 50-300 ms reverb level; how far the natural decay holds and how it ends."""
    x, fs = sf.read(path, always_2d=True)
    a = np.abs(x).max(1)
    on = int(np.argmax(a > a.max() * 0.01))
    end = np.nonzero(np.any(x != 0, 1))[0][-1]
    y = sosfiltfilt(butter(4, FLOOR_BAND, btype="band", fs=fs, output="sos"), x[:end], axis=0)
    ref = 10 * np.log10(np.mean(y[on + int(0.05 * fs):on + int(0.3 * fs)] ** 2, 0))
    n = int(FRAME * fs)
    E = np.array([10 * np.log10(np.mean(y[i:i + n] ** 2, 0) + 1e-30) for i in range(on, len(y) - n, n)])
    out = []
    for c in range(2):
        e = E[:, c] - ref[c]
        t = np.arange(len(e)) * FRAME
        w = 6                                                # 150 ms local slope
        sl = np.array([np.polyfit(t[i:i + w], e[i:i + w], 1)[0] for i in range(len(e) - w)])
        q = int(np.argmin(e))                                # quietest frame, and the slope over 250 ms before it
        last = e[max(q - 10, 0):q + 1]
        r = dict(quietest_rel=float(e[q]), quietest_dbfs=float(E[q, c]),
                 slope_before_quietest=float(np.polyfit(np.arange(len(last)) * FRAME, last, 1)[0]) if len(last) > 2 else None)
        i0, i1 = int(0.3 / FRAME), (int(np.argmax(e < -30)) if np.any(e < -30) else 0)
        if i0 + w < i1 < len(sl):
            s0 = float(np.median(sl[i0:i1]))                 # natural decay rate (dB/s)
            k = i1
            while k < len(sl) and 0.7 * s0 >= sl[k] >= 1.4 * s0:
                k += 1
            how = "end" if k >= len(sl) else ("flattens" if sl[k] > 0.7 * s0 else "fade")
            r.update(decay_db_s=s0, natural_to=float(e[min(k + w // 2, len(e) - 1)]), ends_by=how)
            if how == "flattens" and e[k] < -70:
                seg = slice(k + w // 2, min(k + w // 2 + 6, len(e)))
                r.update(plateau_rel=float(np.median(e[seg])), plateau_dbfs=float(np.median(E[seg, c])))
        out.append(r)
    return out


def floors_stage(ir=None):
    """Is the library's noise floor like the CLA capture's? Natural-decay depth, plateaus and fades, relative to 50-300 ms."""
    fl = {"CLA": [CLA], **{p: files(p) for p in "SML"}}
    with ProcessPoolExecutor(10) as ex:
        res = {k: list(ex.map(_floor_file, v, chunksize=8)) for k, v in fl.items()}
    out = {"band_hz": FLOOR_BAND, "reference": "50-300 ms after onset, in band",
           "format_floor_dbfs": {"24-bit rounding, full band": round(20 * np.log10(2.0 ** -23 / 12 ** 0.5), 1),
                                 "24-bit rounding, in band at 96 kHz": round(20 * np.log10(2.0 ** -23 / 12 ** 0.5) + 10 * np.log10((FLOOR_BAND[1] - FLOOR_BAND[0]) / 48000), 1)}}
    for k, rows in res.items():
        o = {}
        for c, ch in enumerate("LR"):
            v = [r[c] for r in rows]
            nat = [r for r in v if "natural_to" in r]
            d = np.array([r["natural_to"] for r in nat]) if nat else np.array([np.nan])
            deep = [r for r in nat if r["natural_to"] < -70]
            pl = np.array([[r["plateau_rel"], r["plateau_dbfs"]] for r in v if "plateau_rel" in r])
            ratio = [r["slope_before_quietest"] / r["decay_db_s"] for r in nat if r["slope_before_quietest"] is not None]
            o[ch] = {"files": len(v), "with_natural_decay": len(nat),
                     "natural_decay_holds_to_db": {"deepest": round(float(np.nanmin(d)), 1), "p5": round(float(np.nanpercentile(d, 5)), 1),
                                                   "median": round(float(np.nanmedian(d)), 1)},
                     "below_-70_db": len(deep), "below_-70_ends_by": dict(collections.Counter(r["ends_by"] for r in deep)),
                     "plateau_rel_db_median": round(float(np.median(pl[:, 0])), 1) if len(pl) else None,
                     "plateau_dbfs_median": round(float(np.median(pl[:, 1])), 1) if len(pl) else None,
                     "plateau_files": len(pl),
                     "quietest_frame_rel_db_median": round(float(np.median([r["quietest_rel"] for r in v])), 1),
                     "quietest_frame_dbfs_median": round(float(np.median([r["quietest_dbfs"] for r in v])), 1),
                     "slope_before_quietest_over_decay_median": round(float(np.median(ratio)), 2) if ratio else None}
        out[k] = o
        print(k, json.dumps(o))
    return out


# ---------------------------------------------------------------- S repeat captures with a per-capture clock ratio
CLOCK_BAND = butter(4, [100, 10000], btype="band", fs=96000, output="sos")


def _warp(x, shift, ppm, block=4096):
    """x read at n*(1 + ppm*1e-6) + shift (samples), band-limited; the delay is constant within each block."""
    n = len(x)
    N = 1 << int(np.ceil(np.log2(n + 256)))
    X, f = np.fft.rfft(x, N), np.fft.rfftfreq(N)
    out = np.empty(n)
    for i in range(0, n, block):
        j = min(i + block, n)
        d = shift + (i + j) / 2 * ppm * 1e-6
        out[i:j] = np.fft.irfft(X * np.exp(2j * np.pi * f * d), N)[i:j]
    return out


def _clock_pair(args):
    from scipy.optimize import minimize
    fa, fb = args
    A, fs = sf.read(fa, always_2d=True)
    B, _ = sf.read(fb, always_2d=True)
    a, b = sosfiltfilt(CLOCK_BAND, A, axis=0), sosfiltfilt(CLOCK_BAND, B, axis=0)
    i0, i1, m = int(0.05 * fs), int(0.6 * fs), 3000
    c = np.correlate(b[i0 - m:i1 + m, 1], a[i0:i1, 1], "valid")        # c[k]: b[n + k - m] ~ a[n]
    lag = int(np.argmax(c)) - m
    b = np.roll(b, -lag, axis=0)
    n = min(len(a), len(b), int(1.0 * fs))
    a, b = a[:n], b[:n]

    def cost(p):
        y = _warp(b[:, 1], p[0], p[1])[i0:i1]
        x = a[i0:i1, 1]
        return np.sum((x - y * np.dot(x, y) / np.dot(y, y)) ** 2) / np.sum(x ** 2)
    s0 = min(np.linspace(-1, 1, 21), key=lambda s: cost([s, 0.0]))
    p_shift = minimize(lambda p: cost([p[0], 0.0]), [s0], method="Nelder-Mead").x[0]
    p = minimize(cost, [p_shift, 0.0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-12)).x

    def nulls(sh, ppm):
        r = []
        for t0, t1 in WIN[:4]:
            j0, j1 = int(t0 * fs), int(t1 * fs)
            row = []
            for ch in range(2):
                x, y = a[j0:j1, ch], _warp(b[:, ch], sh, ppm)[j0:j1]
                row.append(float(10 * np.log10(np.sum((x - y * np.dot(x, y) / np.dot(y, y)) ** 2) / np.sum(x ** 2))))
            r.append(row)
        return r
    return dict(pair=[os.path.basename(fa), os.path.basename(fb)], lag=lag, shift=float(p[0]), ppm=float(p[1]),
                shift_only=nulls(p_shift, 0.0), clock=nulls(p[0], p[1]))


def clock_stage(ir=None, n_pairs=40):
    """S AMPL pairs (90 vs 60) in the capture domain: null with a sub-sample shift only, then with shift + clock ratio
    (fitted on R over 0.05-0.6 s, applied to both channels)."""
    jobs = [(d[90], d[60]) for d in list(ampl_groups("S").values())[:n_pairs] if 90 in d and 60 in d]
    with ProcessPoolExecutor(10) as ex:
        res = list(ex.map(_clock_pair, jobs))
    ppm = np.array([r["ppm"] for r in res])
    out = {"pairs": len(res), "band_hz": [100, 10000], "fit_window_s": [0.05, 0.6],
           "ppm": {"median": round(float(np.median(ppm)), 2), "min": round(float(ppm.min()), 2), "max": round(float(ppm.max()), 2)},
           "ppm_by_lf": {str(lf): sorted({round(r["ppm"], 1) for r, (fa, _) in zip(res, jobs) if labels(fa)["lf"] == lf})
                         for lf in sorted({labels(fa)["lf"] for fa, _ in jobs})},
           "windows": {}}
    for j, (t0, t1) in enumerate(WIN[:4]):
        s = np.array([r["shift_only"][j] for r in res])
        k = np.array([r["clock"][j] for r in res])
        out["windows"][f"{t0}-{t1}"] = {"shift_only_LR": np.round(np.median(s, 0), 1).tolist(),
                                        "clock_LR": np.round(np.median(k, 0), 1).tolist(),
                                        "clock_p10_p90_R": np.round(np.percentile(k[:, 1], [10, 90]), 1).tolist()}
    print("S clock", json.dumps(out))
    return out


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=["all", "clusters", "loops", "lags", "pairs", "floors", "clock"])
    ap.add_argument("--out", default=os.path.join(ROOT, "analysis", "ml_structure"))
    a = ap.parse_args()
    ir = IRs(tagged(os.path.join(a.out, "cache")))
    print("device-rate IRs ...", flush=True)
    ir.build(files("S") + files("M") + files("L"))
    path = os.path.join(a.out, "summary.json")
    summary = json.load(open(path)) if os.path.exists(path) else {}
    stages = {"clusters": cluster_stage, "loops": loop_stage, "lags": lag_stage, "pairs": pairs_stage,
              "floors": floors_stage, "clock": clock_stage}
    for name, fn in stages.items():
        if a.stage in ("all", name):
            print(f"== {name}", flush=True)
            summary[name] = fn(ir)
            with open(path, "w") as f:
                json.dump(summary, f, indent=1)
    print("wrote", path)


if __name__ == "__main__":
    main()
