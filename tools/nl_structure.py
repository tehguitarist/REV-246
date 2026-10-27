"""Program N (Non-Lin) from the Tomkin 96 kHz library: device-rate IRs and the structure tests of docs/NonLin_program.md.

    .venv/bin/python tools/nl_structure.py --stage cache [--workers 6]

The Non-Lin files are trimmed to their first sample (DELAY 0 throughout), so no onset is searched: file sample 0 is put
at device index K0 and the files are aligned on their content afterwards. Deconvolved with the chain3 kernel (the M / L
default, ml_structure.md 12, 18) unless EMT_KERNEL says otherwise, resampled to 19.2 MHz / 728, L half a sample late
as in the other programs (library_notes.md 3). Cache: analysis/nl/cache_<kernel>/<name>.npy, (n, 2) float32, not
committed.
"""
import argparse, functools, glob, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import soundfile as sf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
os.environ.setdefault("EMT_KERNEL", "chain3")
import core  # noqa: E402
from core import PAD, Kernel, deconvolve, fs_dev, load_config, to_device_grid  # noqa: E402

# A wider analysis band than the harness's (11.0-12.3 kHz raised cosine), the default: the unit's chain is only 4 dB
# down at 12 kHz and 26 dB at 12.5 kHz, so 12.0-12.8 kHz keeps the content that is there (NonLin_program 11.3). Own
# cache. NL_WIDE=0: the harness's band (cache_chain3).
WIDE = os.environ.get("NL_WIDE", "1") == "1"
if WIDE:
    core.F_ROLL, core.F_HI = 12000.0, 12800.0

CFG = load_config()
FSD = fs_dev(CFG)
LIB = os.path.join(ROOT, "irs", "EMT 246 captures", "96kHz", "Non-Lin")
KTAG = os.environ["EMT_KERNEL"]
KERNEL = os.path.join(ROOT, "analysis", "library", "tomkin96", "kernels", KTAG + ".npz")
OUT = os.path.join(ROOT, "analysis", "nl")
CACHE = os.path.join(OUT, f"cache_{KTAG}" + ("_wide" if WIDE else ""))
EPS_REL = 1e-3                 # core.deconvolve's regularisation
K0 = 50                        # device index of file sample 0
LR_US = 19.19                  # library L/R offset (library_notes.md 3)


def labels(path):
    """LP_DELAY_SPAN_EARLY_LATE, one implied decimal on the last three (library_notes.md 2)."""
    lp, delay, span, early, late = (int(v) for v in os.path.basename(path).replace("_dc.wav", "").split("_"))
    return dict(span=span / 10, early=early / 10, late=late / 10)


def key(l):
    return (l["span"], l["early"], l["late"])


def files():
    return sorted(glob.glob(os.path.join(LIB, "*.wav")))


def table():
    """{(span, early, late): path}"""
    return {key(labels(f)): f for f in files()}


def _cache_path(path):
    return os.path.join(CACHE, os.path.basename(path)[:-7] + ".npy")


def _make(path):
    out = _cache_path(path)
    if os.path.exists(out):
        return
    x, fs = sf.read(path, always_2d=True)
    d = deconvolve(np.vstack([np.zeros((PAD, 2)), x]), fs, Kernel.load(KERNEL))
    anchor = float(PAD)
    k0 = int(np.floor(anchor / fs * FSD))
    h = np.stack([to_device_grid(d[:, 0], fs, anchor + LR_US * 1e-6 * fs, k0, CFG),
                  to_device_grid(d[:, 1], fs, anchor, k0, CFG)], 1)
    np.save(out, h[k0 - K0:].astype(np.float32))


def ir(path_or_key):
    """Device-rate IR (n, 2), file sample 0 at K0."""
    p = table()[path_or_key] if isinstance(path_or_key, tuple) else path_or_key
    return np.load(_cache_path(p)).astype(np.float64)


def processing_response(f):
    """The deconvolution's own response to a flat input, per channel (len(f), 2), zero phase: the regularised inverse
    of the chain kernel times the analysis band window (what every cached IR carries and the model must carry too)."""
    k = Kernel.load(KERNEL)
    fs = k.fs
    N = 1 << 18
    ff = np.fft.rfftfreq(N, 1 / fs)
    G = k.response(N, fs)
    G = G if G.ndim == 2 else np.stack([G, G], 1)
    g2 = np.abs(G) ** 2
    win = core._band_window(ff)
    if k.hp_hz:                                       # the low band edge core.deconvolve applies with high-pass corners
        win = win * np.sin(np.pi / 2 * np.clip((ff - core.F_LO[0]) / (core.F_LO[1] - core.F_LO[0]), 0, 1)) ** 2
    R = g2 / (g2 + EPS_REL * np.max(np.abs(G), 0) ** 2) * win[:, None]
    R = R / R[np.argmin(np.abs(ff - 1000))]          # unit at 1 kHz (the normalisation is the tap gain's)
    out = np.stack([np.interp(f, ff, R[:, c]) for c in (0, 1)], 1)
    lb = int(os.environ.get("NL_LOWBAND", "2"))      # the measured low-band transfer of the captures (tools/lowband_joint.py; complex, so no longer zero phase); 2 = adopted, 0 = off
    if lb:
        import lowband_joint
        out = out * lowband_joint.response(f, lb)
    return out


@functools.lru_cache(maxsize=1)
def _resampler_response():
    """(frequency in Hz at the capture rate, gain) of core.to_device_grid's resample_poly FIR (the 96 kHz -> device rate
    step every cached IR went through), DC-normalised."""
    from scipy.signal import firwin
    fs = 96000
    up, down = core.dev_ratio(fs, CFG)
    mr = max(up, down)
    h = firwin(2 * 10 * mr + 1, 1.0 / mr, window=("kaiser", 5.0))
    N = 1 << 22
    H = np.abs(np.fft.rfft(h, N))
    return np.fft.rfftfreq(N) * fs * up, H / H[0]


def resampler_gain(f):
    """Gain of to_device_grid's polyphase low-pass at f (Hz), zero phase. Every cached IR carries it (NonLin_program
    11.6): it is the 'top-end EQ' of 11.3 (-0.22 / -0.90 / -1.50 / -2.33 / -3.65 dB at 11.5 / 12 / 12.25 / 12.5 / 12.8 kHz),
    not a property of the unit."""
    ff, H = _resampler_response()
    return np.interp(f, ff, H)


def cache_stage(workers):
    os.makedirs(CACHE, exist_ok=True)
    todo = [p for p in files() if not os.path.exists(_cache_path(p))]
    print(f"{len(todo)} of {len(files())} to build", flush=True)
    with ProcessPoolExecutor(workers) as ex:
        for i, _ in enumerate(ex.map(_make, todo, chunksize=4)):
            if i % 100 == 0:
                print(i, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="cache")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", type=int, default=0)
    a = ap.parse_args()
    if a.stage == "cache":
        if a.only:
            global files
            fl = files()[:a.only]
            files = lambda: fl  # noqa: E731
        cache_stage(a.workers)


if __name__ == "__main__":
    main()
