"""Program N (Non-Lin): the whole program as a static model, and its bench against the Tomkin 96 kHz library
(docs/NonLin_program.md §10-11).

    .venv/bin/python tools/nl_model.py --stage noise                 # the capture noise floor (LATE 0 tails)
    .venv/bin/python tools/nl_model.py --stage bench [--workers 6] [--only N]
    .venv/bin/python tools/nl_model.py --stage report

Model (device samples, mono in, stereo out; everything from nl_early / nl_late):
  front taps    L 41 / 410 / 463, R 41 / 278 / 463, 1 : 2 : 1 main-tap units, gains on the EARLY law (front ticks);
  early         AP(124) -> 27 main taps in 14 stages (SPAN runs N_STAGES of them, EARLY sets their law gains) ->
                AP(180) L / AP(188) R; R taps 221 earlier, R's first-tap direct term at 465 in the front frame;
  late          AP(124) -> the main-tap line at FEED[SPAN] (R 220 earlier) -> AP(314) -> three flat loops 1108 / 1254 /
                1386 at the Jot gain for RT = LATE -> four equal signed read taps per loop and channel, 0.5 main-tap
                units, times the late gain (EARLY > 0: the law at stage min(ns, 13); else 0 dB); LATE 0: no late.
Per capture (bench): the front taps' trim shift (measured), the branch offset o (AP(124) branch re the front taps), the
clock ratio, one gain per channel; nothing else is fitted. Nulls 100 Hz-10 kHz per window (s from the front frame's
zero): 0-0.05, 0.05-0.15, 0.15-0.3, 0.3-0.6, 0.6-1.0, 1.0-2.0, each up to FADE_GUARD before the file's end, with the
noise-limited null of the window (the capture noise floor of --stage noise over the window's data energy) for excess.
"""
import argparse, functools, json, os, sys
from concurrent.futures import ProcessPoolExecutor
for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("NL_KERNEL", "0")
import numpy as np
from scipy.optimize import minimize
from scipy.signal import lfilter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import defects  # noqa: E402
import nl_early as ne  # noqa: E402
import nl_late as nlt  # noqa: E402
import nl_structure as nl  # noqa: E402

FS, K0 = nl.FSD, nl.K0
LATE_AMP = 0.5                                  # read taps in main-tap units (0.2494-0.2500 of the 410 front tap)
WINDOWS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 2.0)]
FADE_GUARD = 0.75                               # s before the file's end (S's 0.70 s fade ramp + a guard)
NOISE = os.path.join(nl.OUT, "noise.json")
BENCH = os.path.join(nl.OUT, "bench")


def late_db(span, early):
    """The late section's gain (dB, main-tap units re EARLY 10): the EARLY law at stage min(ns, 13) (EARLY > 0)."""
    if early >= 10 or early < 0:
        return 0.0
    j = min(ne.N_STAGES[span], len(ne.STAGES) - 1)
    return -60 * (4 + 9 * j) * ne.TICK / early


@functools.lru_cache(maxsize=8)
def spectra(key, ch, N):
    """Frequency-domain model of one channel at o = 0 (front frame zero at index 0): (front part, branch part)."""
    span, early, late = key
    ns = ne.N_STAGES[span]
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    fdb, sdb = ne.law_db(span, early, ch)
    front = sum(ne.FRONT_AMP[ch][i] * 10 ** (fdb[i] / 20) * np.exp(ph * f) for i, f in enumerate(ne.FRONT[ch]))
    d = ne.diffuser(ch, N // 2)
    Df = np.fft.rfft(d, N)
    branch = np.zeros(len(ph), complex)
    for t, s, i in ne.taps(ch, ns):
        g = s * 10 ** (sdb[[j for j, st in enumerate(ne.STAGES) if i in st][0]] / 20)
        if ch == 1 and i == 0:                  # R's first tap: direct term in the front frame
            branch += g * (Df - d[0]) * np.exp(ph * t)
            front = front + g * d[0] * np.exp(ph * ne.R_FIRST_DIRECT)
        else:
            branch += g * Df * np.exp(ph * t)
    if late > 0:
        u = np.zeros(N)
        u[0] = 1.0
        U = np.fft.rfft(ne.ap(ne.ap(u, ne.AP_IN), nlt.AP_LATE))
        gl = LATE_AMP * 10 ** (late_db(span, early) / 20)
        base = nlt.FEED[span] + (nlt.R_LATE if ch else 0)
        for P, qs in nlt.READS.items():
            v = np.zeros(N)
            v[0] = 1.0
            V = np.fft.rfft(lfilter([1.0], np.r_[1.0, np.zeros(P - 1), -nlt.loop_gain(P, late)], v))
            taps = sum(s * np.exp(ph * (base + q)) for s, q in zip(nlt.READ_SIGNS[ch][P], qs))
            branch += gl * U * V * taps
    return front, branch


def render(key, n, o=0.0, fr=0.0, band=True):
    """The model (n, 2) with the front frame's zero at K0 + fr and the branch o samples later; band-limited to
    100 Hz-10 kHz unless band=False."""
    N = 1 << int(np.ceil(np.log2(n + 2 * ne.TAPS[-1] + 4096)))
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    out = []
    for ch in (0, 1):
        F, B = spectra(key, ch, N)
        out.append(np.fft.irfft((F + B * np.exp(ph * o)) * np.exp(ph * (K0 + fr)), N)[:n])
    y = np.array(out).T
    return sosfiltfilt(ne.BAND, y, axis=0) if band else y


def usable_end(h):
    a = np.abs(h)
    last = max(int(np.nonzero(a[:, ch] > a[:, ch].max() * 1e-6)[0][-1]) for ch in (0, 1))
    return last - int(min(FADE_GUARD * FS, 0.3 * last))       # short files: at most 30 % of the file


# ---------------------------------------------------------------- noise floor
def noise_stage():
    """Absolute noise power per sample (100 Hz-10 kHz) in the LATE 0 captures, from 0.25 s after the last active main
    tap to the fade (the early section is then 110 dB down); median over files, per channel."""
    vals = {0: [], 1: []}
    for k, p in sorted(nl.table().items()):
        if k[2] != 0.0 or k in defects.NL_DEFECTS:
            continue
        h = sosfiltfilt(ne.BAND, nl.ir(p), axis=0)
        ns = ne.N_STAGES[k[0]]
        i0 = K0 + (ne.taps(0, ns)[-1][0] if ns else 600) + int(0.25 * FS)
        i1 = usable_end(h)
        if i1 - i0 < int(0.1 * FS):
            continue
        for ch in (0, 1):
            vals[ch].append(float(np.mean(h[i0:i1, ch] ** 2)))
    out = {f"ch{ch}": float(np.median(v)) for ch, v in vals.items()}
    out["n_files"] = len(vals[0])
    out["spread_db"] = [float(10 * np.log10(np.percentile(vals[0], 90) / np.percentile(vals[0], 10)))]
    json.dump(out, open(NOISE, "w"), indent=1)
    print(out)


# ---------------------------------------------------------------- bench
def bench_file(key):
    out_f = os.path.join(BENCH, "_".join(str(v) for v in key) + ".json")
    if os.path.exists(out_f):
        return key, json.load(open(out_f))
    noise = json.load(open(NOISE))
    raw = nl.ir(nl.table()[key])
    h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = ne.front_shift(h)
    end = usable_end(h)
    n = min(end, K0 + int(fr) + int(WINDOWS[-1][1] * FS))
    span, early, late = key
    ns = ne.N_STAGES[span]
    # branch offset candidates: the early section's template match (SPAN > 0) and the common small offsets
    cands = [0, -64, -84, -87, -163] + (ne.coarse_offsets(h, key, fr) if ns else [])
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])

    def resid(o, ppm):
        y = nlt.unwarp(rawp, ppm)[:n] if ppm else rawp[:n]
        y = sosfiltfilt(ne.BAND, y, axis=0)
        m = render(key, n, o, fr)
        g = [(y[:, c] @ m[:, c]) / (m[:, c] @ m[:, c]) for c in (0, 1)]
        return y, m * np.array(g), g
    i0 = K0 + int(fr) + ne.FRONT[0][0] - 20

    noise = json.load(open(NOISE))
    wins = []                                           # the windows with signal: each counts equally in the fit
    for t0, t1 in WINDOWS:
        a, b = max(i0, K0 + int(fr) + int(t0 * FS)), min(n, K0 + int(fr) + int(t1 * FS))
        if b - a >= int(0.02 * FS) and np.sum(h[a:b, 0] ** 2) >= 10 * noise["ch0"] * (b - a):
            wins.append((a, b))

    def cost(v):
        if abs(v[1] - nlt.CLOCK_PPM) > 3:                   # the session's clock is 11.6 +- 0.5 ppm
            return 1e3
        y, m, _ = resid(v[0], v[1])
        return float(sum(np.log(np.sum((y[a:b] - m[a:b]) ** 2) / np.sum(y[a:b] ** 2)) for a, b in wins))
    best = min(((cost([float(o), nlt.CLOCK_PPM]), o) for o in sorted({c + d for c in cands for d in (-1, 0, 1)})))
    # the valley in o is ~0.1 sample wide: a fine scan around the best integer, then a small simplex
    best = min(((cost([float(o), nlt.CLOCK_PPM]), o) for o in best[1] + np.arange(-1.0, 1.001, 0.02)))
    r = minimize(cost, [best[1], nlt.CLOCK_PPM], method="Nelder-Mead",
                 options={"xatol": 1e-4, "fatol": 1e-7, "maxiter": 200,
                          "initial_simplex": [[best[1], 11.6], [best[1] + 0.02, 11.6], [best[1], 12.1]]})
    o, ppm = float(r.x[0]), float(r.x[1])
    y, m, g = resid(o, ppm)
    res = {"o": o, "ppm": ppm, "gain": [float(v) for v in g], "front_shift": fr, "windows": {}}
    for t0, t1 in WINDOWS:
        a = max(i0, K0 + int(fr) + int(t0 * FS))
        b = min(n, K0 + int(fr) + int(t1 * FS))
        if b - a < int(0.02 * FS):
            continue
        w = {}
        if np.sum(y[a:b, 0] ** 2) < 10 * noise["ch0"] * (b - a):   # within 10 dB of the noise: nothing to score
            continue
        for c in (0, 1):
            e = float(np.sum(y[a:b, c] ** 2))
            w[f"null{c}"] = float(10 * np.log10(np.sum((y[a:b, c] - m[a:b, c]) ** 2) / e))
            w[f"floor{c}"] = float(10 * np.log10(noise[f"ch{c}"] * (b - a) / e))
        res["windows"][f"{t0}-{t1}"] = w
    os.makedirs(BENCH, exist_ok=True)
    json.dump(res, open(out_f, "w"))
    return key, res


def bench_stage(workers, only):
    keys = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS)
    if only:
        keys = keys[:: max(1, len(keys) // only)][:only]
    with ProcessPoolExecutor(workers) as ex:
        for k, res in ex.map(bench_file, keys, chunksize=2):
            w = res["windows"]
            print(k, f"o {res['o']:+.1f} ppm {res['ppm']:+.1f} | " + " ".join(
                f"{t}: {v['null0']:.1f}/{v['null1']:.1f} ({v['floor0']:.0f})" for t, v in w.items()), flush=True)


def report_stage():
    rows = []
    for f in sorted(os.listdir(BENCH)):
        k = tuple(float(v) for v in f[:-5].split("_"))
        rows.append((k, json.load(open(os.path.join(BENCH, f)))))
    print(f"{len(rows)} captures (NL_DEFECTS left out)")
    print("window      n    null L / R (median)   floor L     excess L / R (median, windows above the floor)")
    for t0, t1 in WINDOWS:
        t = f"{t0}-{t1}"
        ws = [r["windows"][t] for _, r in rows if t in r["windows"]]
        if not ws:
            continue
        nl0 = np.median([w["null0"] for w in ws])
        nl1 = np.median([w["null1"] for w in ws])
        fl = np.median([w["floor0"] for w in ws])
        ex = [max(w[f"null{c}"] - w[f"floor{c}"], 0) for w in ws for c in (0, 1)]
        print(f"{t:10s} {len(ws):4d}   {nl0:6.1f} / {nl1:6.1f}        {fl:6.1f}     {np.median(ex):5.1f}")
    for name, sel in [("EARLY > 0", lambda k: 0 < k[1] < 10), ("EARLY < 0", lambda k: k[1] < 0),
                      ("EARLY 10", lambda k: k[1] >= 10), ("LATE 0", lambda k: k[2] == 0), ("LATE > 0", lambda k: k[2] > 0)]:
        sub = [r for k, r in rows if sel(k)]
        line = []
        for t0, t1 in WINDOWS:
            t = f"{t0}-{t1}"
            ws = [r["windows"][t] for r in sub if t in r["windows"]]
            if ws:
                line.append(f"{t}: {np.median([w['null0'] for w in ws]):.1f}")
        print(f"  {name:10s} ({len(sub)}): " + "  ".join(line))


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--stage", default="bench")
    ap_.add_argument("--workers", type=int, default=6)
    ap_.add_argument("--only", type=int, default=0)
    a = ap_.parse_args()
    {"noise": noise_stage, "report": report_stage}.get(a.stage, lambda: bench_stage(a.workers, a.only))()


if __name__ == "__main__":
    main()
