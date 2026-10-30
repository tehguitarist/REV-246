"""Model-free: two raw 96 kHz captures of one (SPAN, EARLY) differing in LATE. Lag (96 kHz samples, sub-sample) of the
front-tap segment and of the main-tap segment (before the late section starts), channel R and L, band 300 Hz-8 kHz."""
import sys, json, os, glob
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
import nl_structure as nl
FS = 96000; STEP = FS / nl.FSD
BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7_lb2_hf_g5"
SOS = butter(4, [300, 8000], btype="band", fs=FS, output="sos")

def lag(a, b):
    """sub-sample lag of b re a (b later = positive), by cross-correlation and parabolic peak, upsampled x8 by FFT"""
    n = 1 << int(np.ceil(np.log2(2 * len(a)))); U = 8
    X = np.fft.rfft(b, n) * np.conj(np.fft.rfft(a, n))
    Xp = np.zeros(n * U // 2 + 1, complex); Xp[:len(X)] = X
    xc = np.fft.irfft(Xp, n * U)
    i = int(np.argmax(xc)); y0, y1, y2 = xc[i - 1], xc[i], xc[(i + 1) % len(xc)]
    d = i + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    if d > n * U / 2: d -= n * U
    return d / U

def load(key):
    x, fs = sf.read(nl.table()[key], always_2d=True); assert fs == FS
    return sosfiltfilt(SOS, x, axis=0)

def front_idx(x):
    """the L 410 tap peak at 96k (largest peak in the first 30 ms)"""
    return int(np.argmax(np.abs(x[: int(0.03 * FS), 0])))

span, early = float(sys.argv[1]), float(sys.argv[2])
keys = sorted(k for k in nl.table() if k[0] == span and k[1] == early and k[2] > 0)
ref = None
for k in keys:
    b = json.load(open(os.path.join(BENCH, "_".join(str(v) for v in k) + ".json")))
    nb = (b["o"] - 1) * STEP
    x = load(k); i = front_idx(x)
    if ref is None:
        ref, xr, ir_ = k, x, i; nref = nb
    # front segment: from 410-tap - 380 dev to + 60 dev (L: 41..463), R similar; main: 26 ms .. min(SPAN, 0.12) s after
    t0 = lambda dev: int(dev * STEP)
    seg = {"front": (t0(-380), t0(60)), "main": (t0(685 - 410 - 10), t0(685 - 410) + int(0.08 * FS))}
    out = []
    for name, (a, bb) in seg.items():
        for c in (0, 1):
            out.append(lag(xr[ir_ + a: ir_ + bb, c], x[i + a: i + bb, c]) + (i - ir_))
    print(k, f"bench n={nb:+.2f} (rel {nb - nref:+.2f}) | front L {out[0]:+.2f} R {out[1]:+.2f} | main L {out[2]:+.2f} R {out[3]:+.2f} | main-front L {out[2]-out[0]:+.2f} R {out[3]-out[1]:+.2f}")
