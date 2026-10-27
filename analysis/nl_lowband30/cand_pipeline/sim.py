"""Forward simulation of the Program N capture pipeline (cand_pipeline): model device-rate response -> 96 kHz -> chain kernel
-> (file length, fade) -> core.deconvolve -> to_device_grid, to compare with the cached IR. Read-only on tools/, analysis/nl/."""
import os, sys, json
os.environ["NL_BK"] = "1"; os.environ.setdefault("NL_LOWBAND", "0")
import numpy as np, soundfile as sf
from scipy.signal import resample_poly
ROOT = "/Users/leighpierce/Git/246"
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "tools", "harness"))
import nl_model as M, nl_structure as nl, nl_late as nlt, nl_early as ne          # noqa: E402
import core                                                                         # noqa: E402
BENCH = os.path.join(ROOT, "analysis", "nl", "bench_warp2_bk_kf_a0.7" + ("_lb" + os.environ["NL_LOWBAND"] if os.environ["NL_LOWBAND"] != "0" else ""))
FS96 = 96000
K0, FSD, PAD, CFG = nl.K0, nl.FSD, core.PAD, nl.CFG
KERN = core.Kernel.load(nl.KERNEL)
UP, DOWN = core.dev_ratio(FS96, CFG)               # 25, 91
RATIO = FS96 / FSD
KK0 = int(np.floor(PAD / FS96 * FSD))              # to_device_grid's k0 as in nl_structure._make


def align(key):
    return json.load(open(os.path.join(BENCH, "_".join(str(v) for v in key) + ".json")))


def file_path(key):
    return nl.table()[key]


def model_cache(key, n, r):
    """Device-rate model (n, 2) on the cache grid (file sample 0 at K0), no processing, no band limit, bench gain."""
    return M.render(key, n, r["o"], r["front_shift"], band=False, proc=False) * np.array(r["gain"])


def forward_warp(m, ppm):
    """raw[i] = m(i / (1 + ppm 1e-6)): inverse of nlt.unwarp (the cache clock vs the model's)."""
    up = 32
    u = nlt.upsample(m, up)
    t = np.arange(len(m)) / (1 + ppm * 1e-6) * up
    return np.stack([np.interp(t, np.arange(len(u)), u[:, c]) for c in range(m.shape[1])], 1)


def to_capture_rate(raw):
    """Inverse of nl._make's to_device_grid chain (per channel; L anchor is LR_US later): the 96 kHz signal dc whose
    to_device_grid(dc, anchor_c, KK0) gives the cache raw. Returns (n96, 2), index = file sample + PAD."""
    out = []
    for c in (0, 1):
        anchor = PAD + (nl.LR_US * 1e-6 * FS96 if c == 0 else 0.0)
        full = np.concatenate([np.zeros(KK0 - K0), raw[:, c]])        # device samples k = 0 ...
        u = resample_poly(full, DOWN, UP)                              # u[j] = signal at capture index j (device k = j UP/DOWN)
        s = anchor - KK0 * RATIO                                       # dc[t] = u[t - s]
        N = 1 << int(np.ceil(np.log2(2 * len(u))))
        fr = np.fft.rfftfreq(N)
        dc = np.fft.irfft(np.fft.rfft(u, N) * np.exp(-2j * np.pi * fr * s), N)[:len(u)]
        out.append(dc)
    return np.array(out).T[PAD:]                                       # index = file sample (the part before the file start is lost)


def chain(dc, nout, N=None):
    """dc (n, 2) through the chain kernel (FIR + lead + analytic high-pass corners), circular FFT of size N (large)."""
    N = N or 1 << int(np.ceil(np.log2(len(dc) + 6 * FS96)))
    G = KERN.response(N, FS96)
    X = np.fft.rfft(dc, N, axis=0) * G
    return np.fft.irfft(X, N, axis=0)[:nout]


def pipeline(x96, key_unused=None):
    """core.deconvolve + to_device_grid exactly as nl_structure._make, returns the cache-grid IR (n,2) with sample 0 at K0."""
    d = core.deconvolve(np.vstack([np.zeros((PAD, 2)), x96]), FS96, KERN)
    anchor = float(PAD)
    h = np.stack([core.to_device_grid(d[:, 0], FS96, anchor + nl.LR_US * 1e-6 * FS96, KK0, CFG),
                  core.to_device_grid(d[:, 1], FS96, anchor, KK0, CFG)], 1)
    return h[KK0 - K0:], d


from scipy.signal import butter, sosfiltfilt
_SOS = butter(4, [300, 8000], btype="band", fs=FS96, output="sos")


def fit_scale(x, xs, t1=0.4):
    """per-channel scale of x_sim to the real file (300 Hz-8 kHz, first t1 s)"""
    n = int(t1 * FS96)
    xb, sb = sosfiltfilt(_SOS, x[:n], axis=0), sosfiltfilt(_SOS, xs[:n], axis=0)
    return np.array([(xb[:, c] @ sb[:, c]) / (sb[:, c] @ sb[:, c]) for c in (0, 1)])


def fade_curve(x, xs, A, blk=2400):
    """block rms ratio real / sim (300 Hz-8 kHz) -> estimated linear fade: returns (t_start, gain array per sample) or (None, ones)"""
    xb, sb = sosfiltfilt(_SOS, x, axis=0), sosfiltfilt(_SOS, xs * A, axis=0)
    nf = len(x); nb = nf // blk
    rat = np.array([np.sqrt((xb[i * blk:(i + 1) * blk, 0] ** 2).mean() / max((sb[i * blk:(i + 1) * blk, 0] ** 2).mean(), 1e-30)) for i in range(nb)])
    act = np.array([np.sqrt((sb[i * blk:(i + 1) * blk, 0] ** 2).mean()) for i in range(nb)])
    return rat, act


def make_x(key, r, nf, A, tail_s=0.0, fade=None):
    """simulated file samples (nf + tail, 2): model -> clock -> 96 kHz -> chain kernel, scale A (per ch); fade = (t0, t1) seconds of
    a linear amplitude ramp 1 -> 0 applied to the first nf samples."""
    nsamp = nf + int(tail_s * FS96)
    nn = int(nsamp / RATIO) + 40
    m = model_cache(key, nn, r)
    raw = forward_warp(m, r["ppm"])
    xs = chain(to_capture_rate(raw), nsamp) * A
    if fade is not None:
        t = np.arange(nsamp) / FS96
        g = np.clip((fade[1] - t) / (fade[1] - fade[0]), 0, 1)
        xs = xs * g[:, None]
    return xs
