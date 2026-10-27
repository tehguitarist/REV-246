"""Does the late reverb run on a slightly different time scale than the model? Per capture, capture and model
band-limited to 177 Hz - 2.8 kHz, their relative lag tracked in 40 ms windows (hop 20 ms, cross-correlation peak with
parabolic interpolation, channel R), and a straight line fitted to lag against time (energy-weighted): slope in ppm
(lag samples per sample of time). A per-pass advance d on a loop of length P gives slope -d / P, the same in every
program if the mechanism is a fraction of the line length (so M / L would show d x 2 / x 4). model_s_results §25.16.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/warp.py S [keys]          # dense_mid prep (the §25.7 baseline, clocks per capture)
    .venv/bin/python tools/loop_read/warp.py M|L [keys]                          # ml_lines.fitted (law model, stored clocks and gains)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt

BAND, WIN, HOP, MAXLAG = (177, 2800), 0.040, 0.020, 12


def lag_track(y, m, fs, i0, i1):
    sos = butter(4, BAND, btype="band", fs=fs, output="sos")
    yb, mb = sosfiltfilt(sos, y), sosfiltfilt(sos, m)
    W, H = int(WIN * fs), int(HOP * fs)
    T, L, E = [], [], []
    for a in range(i0, i1 - W, H):
        ys, ms = yb[a:a + W], mb[a - MAXLAG:a + W + MAXLAG]
        c = np.array([ys @ ms[MAXLAG + k:MAXLAG + k + W] for k in range(-MAXLAG, MAXLAG + 1)])
        j = int(np.argmax(c))
        if 0 < j < len(c) - 1:
            d = 0.5 * (c[j - 1] - c[j + 1]) / (c[j - 1] - 2 * c[j] + c[j + 1])
            T.append(a + W / 2 - i0); L.append(j - MAXLAG + d); E.append(ys @ ys)
    return np.array(T), np.array(L), np.array(E)


def slope(T, L, E):
    if len(T) < 4:
        return np.nan, np.nan
    w = E / E.sum()
    A = np.stack([np.ones_like(T), T], 1) * np.sqrt(w)[:, None]
    th, *_ = np.linalg.lstsq(A, L * np.sqrt(w), rcond=None)
    return th[1] * 1e6, th[0]


def run_S(keys):
    import dense_mid as dm
    out = {}
    for k in keys:
        z = np.load(os.path.join(dm.OUT, "prep", k + ".npz"))
        res = []
        for j in sorted(int(n[3:]) for n in z.files if n.startswith("res")):
            y, m = z[f"y{j}"][:, 1].astype(float), z[f"m{j}"][:, 1].astype(float)
            i0, i1 = dm.PRE + int(0.06 * dm.FS), min(dm.PRE + int(min(float(z[f"te{j}"]), 0.6) * dm.FS), len(y) - MAXLAG)
            res.append(slope(*lag_track(y, m, dm.FS, i0, i1))[0])
        out[k] = res
    return out


def run_ML(prog, keys):
    import ml_lines as ml
    out = {}
    for k in keys:
        res = []
        for amp, (y, m, end) in ml.fitted(prog, k).items():
            i0, i1 = ml.K0 + int(0.06 * ml.FS), min(end[1], ml.K0 + int(0.6 * ml.FS), len(y) - MAXLAG)
            res.append(slope(*lag_track(y[:, 1], m[:, 1], ml.FS, i0, i1))[0])
        out[k] = res
    return out


if __name__ == "__main__":
    prog = sys.argv[1]
    keys = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    if prog == "S":
        import dense_mid as dm
        keys = keys or sorted(p[:-4] for p in os.listdir(os.path.join(dm.OUT, "prep")))
        out = run_S(keys)
    else:
        import ml_lines as ml
        keys = keys or sorted(ml.cc.groups(prog))
        out = run_ML(prog, keys)
    for k in keys:
        v = np.array(out[k])
        print(f"{prog} {k}: slope ppm per capture " + " ".join(f"{x:+.0f}" for x in v) + f" | median {np.nanmedian(v):+.0f}", flush=True)
