"""Programs M and L from the Tomkin 96 kHz library (NEXT_STEPS Q3): cluster patterns, reverb onset,
loop recursion and time variance, with program S as the control. Results: docs/ml_structure.md.

    .venv/bin/python tools/ml_structure.py [--stage all|clusters|loops|lags|pairs] [--out analysis/ml_structure]

Writes <out>/summary.json (and caches device-rate IRs as <out>/cache/*.npy, not committed).

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
KERNEL = os.path.join(ROOT, "analysis", "library", "tomkin96", "kernels", "kernel_96000.npz")
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


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=["all", "clusters", "loops", "lags", "pairs"])
    ap.add_argument("--out", default=os.path.join(ROOT, "analysis", "ml_structure"))
    a = ap.parse_args()
    ir = IRs(os.path.join(a.out, "cache"))
    print("device-rate IRs ...", flush=True)
    ir.build(files("S") + files("M") + files("L"))
    path = os.path.join(a.out, "summary.json")
    summary = json.load(open(path)) if os.path.exists(path) else {}
    stages = {"clusters": cluster_stage, "loops": loop_stage, "lags": lag_stage, "pairs": pairs_stage}
    for name, fn in stages.items():
        if a.stage in ("all", name):
            print(f"== {name}", flush=True)
            summary[name] = fn(ir)
            with open(path, "w") as f:
                json.dump(summary, f, indent=1)
    print("wrote", path)


if __name__ == "__main__":
    main()
