"""Where does the sub-20 Hz high-pass come from? Fit high-pass candidates to the low end of the cluster region.

Before the reverb starts (R: 570 device samples after tap 0, L: 240), an S IR is the 8 cluster taps
(known positions and signs, unit gain) through the analogue chain and the capture chain. Below
about 1 kHz the band-limiting part is flat, so after a common low-pass the only unknowns are the
high-pass(es), a gain, a delay and a DC offset. A high-pass with corners w_i leaves a slow negative
tail after each tap whose start level is -sum(w_i) and whose curvature sets the order.

Candidates (corners in Hz):
  schematic   the transcribed coupling corners only (unit, fixed; tools/analog/sample_hold.py)
  sch+1       schematic + one free first-order corner
  sch+2       schematic + two free first-order corners
  free2       two free corners alone (the model_s.py 15 + 15 Hz form)
Fit on R (to 20 ms after tap 0), held out on L (corners fixed, gain/delay/DC refitted, to 8 ms).

    .venv/bin/python tools/analog/lf_highpass.py [--n-tomkin 24] [--lp 1000] [--out lf_highpass.json]
Writes analysis/analog/lf_highpass.json.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import soundfile as sf
from scipy.optimize import least_squares
from scipy.signal import butter, sosfiltfilt, sosfreqz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from sample_hold import FS_DEV, IN_COUPLING_HZ, IN_AMP_COUPLING_HZ, highpass  # noqa: E402

# Coupling corners in the transcribed schematics (Hz): input board C106/R117 and C101/R105,
# input filter C3109/R3134, output filter C3009/R3031. The output 100 uF capacitors depend on the
# load and sit far lower; left out.
SCHEMATIC_HZ = [IN_COUPLING_HZ, IN_AMP_COUPLING_HZ, 1 / (2 * np.pi * 221e3 * 1e-6), 1 / (2 * np.pi * 1e6 * 1e-6)]
CLUSTER = {1: ([0, 143, 194, 353], [1, 1, -1, 1]), 0: ([0, 54, 194], [1, 1, 1])}  # R, L (L's 288 is after onset)
WIN_MS = {1: 20.0, 0: 8.0}
LP_HZ = 1000.0
FS = 96000


def lp_sos():
    return butter(4, LP_HZ, fs=FS, output="sos")


def model_wave(n, t0, A, corners, sos, ch):
    """Cluster taps (fractional positions, device-rate spacing) through the high-passes and the low-pass."""
    N = 1 << int(np.ceil(np.log2(n + 8192)))
    f = np.fft.rfftfreq(N, 1 / FS)
    taps, signs = CLUSTER[ch]
    X = sum(s * np.exp(-2j * np.pi * f * (t0 + d * FS / FS_DEV) / FS) for d, s in zip(taps, signs))
    _, Hlp = sosfreqz(sos, worN=f, fs=FS)
    H = X * highpass(f, corners) * np.abs(Hlp) ** 2          # filtfilt = |H|^2, zero phase
    return A * np.fft.irfft(H, N)[:n]


def fit_one(y, t0_guess, ch, n_free, use_sch, fixed=None, sos=None):
    """y: low-passed channel. Returns parameters and window rms error re window rms."""
    i0 = max(0, int(t0_guess - 0.001 * FS)); i1 = int(t0_guess + WIN_MS[ch] / 1000 * FS)
    n = i1
    base = list(SCHEMATIC_HZ) if use_sch else []

    def corners(q):
        if fixed is not None: return base + list(fixed)
        return base + list(np.exp(q[3:3 + n_free]))

    def res(q):
        m = model_wave(n, q[1], q[0], corners(q), sos, ch) + q[2]
        return (m - y[:n])[i0:i1]
    x0 = [float(np.max(y[i0:i1])), t0_guess, 0.0] + ([np.log(10.0)] * n_free if fixed is None else [])
    if fixed is None and n_free == 2: x0[-1] = np.log(3.0)
    nq = len(x0) - 3
    lo = [-np.inf, t0_guess - 40, -np.inf] + [np.log(0.05)] * nq      # delay kept within 0.4 ms of tap 0
    hi = [np.inf, t0_guess + 40, np.inf] + [np.log(200.0)] * nq
    r = least_squares(res, x0, x_scale="jac", bounds=(lo, hi))
    e = r.fun
    sig = y[i0:i1] - np.mean(y[i0:i1])
    return {"A": float(r.x[0]), "t0": float(r.x[1]), "dc": float(r.x[2]), "corners_hz": [float(c) for c in corners(r.x)],
            "free_hz": [float(c) for c in (np.exp(r.x[3:]) if fixed is None else fixed)],
            "err_db": float(10 * np.log10(np.sum(e ** 2) / np.sum(sig ** 2))), "n_params": 3 + (n_free if fixed is None else 0)}


def noise_db(y, t0, ch, sos):
    """Low-passed pre-onset level (before tap 0, guard 1 ms) against the window signal, as a floor."""
    pre = y[:max(0, int(t0 - 0.0015 * FS))]
    if len(pre) < 48: return None
    i0 = int(t0 - 0.001 * FS); i1 = int(t0 + WIN_MS[ch] / 1000 * FS)
    sig = y[i0:i1] - np.mean(y[i0:i1])
    return float(10 * np.log10(np.mean((pre - pre.mean()) ** 2) * (i1 - i0) / np.sum(sig ** 2)))


CANDIDATES = {"schematic": (0, True), "sch+1": (1, True), "sch+2": (2, True), "free2": (2, False)}


def analyse(path, tap0_s, sos):
    x, fs = sf.read(path, always_2d=True)
    assert fs == FS
    out = {}
    yR = sosfiltfilt(sos, x[:, 1]); yL = sosfiltfilt(sos, x[:, 0])
    t0R = tap0_s * FS
    for name, (nf, sch) in CANDIDATES.items():
        r = fit_one(yR, t0R, 1, nf, sch, sos=sos)
        # held out on L: corners fixed from R
        free = r["free_hz"]
        l = fit_one(yL, t0R + 0.49 * FS / FS_DEV, 0, nf, sch, fixed=free, sos=sos)
        out[name] = {"free_hz": [round(c, 2) for c in free], "sum_fc_hz": round(sum(r["corners_hz"]), 2),
                     "R_err_db": round(r["err_db"], 1), "L_heldout_err_db": round(l["err_db"], 1), "n_params_R": r["n_params"]}
    out["noise_R_db"] = noise_db(yR, t0R, 1, sos)
    return out


def tap0_cla(path):
    x, _ = sf.read(path, always_2d=True)
    r = x[:, 1]
    a = np.abs(r[: int(0.05 * FS)])
    i = int(np.argmax(a > 0.5 * a.max()))                    # first cluster tap (tap 0 is positive)
    i += int(np.argmax(a[i:i + 8]))
    return (i - 19) / FS                                     # kernels put the peak 19 samples after tap 0


def main():
    global LP_HZ
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-tomkin", type=int, default=24)
    ap.add_argument("--lp", type=float, default=LP_HZ, help="analysis low-pass (Hz)")
    ap.add_argument("--out", default="lf_highpass.json")
    a = ap.parse_args()
    LP_HZ = a.lp
    sos = lp_sos()
    res = {"schematic_hz": [round(c, 3) for c in SCHEMATIC_HZ], "lp_hz": LP_HZ, "win_ms": WIN_MS, "files": {}}
    cla = os.path.join(ROOT, "irs", "eugtone_EMT246_CLA_shortRT.wav")
    res["files"]["CLA"] = analyse(cla, tap0_cla(cla), sos)
    print("CLA", json.dumps(res["files"]["CLA"]), flush=True)
    d = pd.read_csv(os.path.join(ROOT, "analysis", "library", "tomkin96", "files.csv"))
    s = d[(d.fs == 96000) & (d.label_program == "S") & (d.cluster_status == "ok")]
    s = s[~s["flags"].fillna("").str.contains("clipped")]
    s = s.sort_values("cluster_to_reverb_db", ascending=False).drop_duplicates("relpath").head(a.n_tomkin)
    for _, r in s.iterrows():
        res["files"][r.relpath] = analyse(os.path.join(ROOT, r.path), r.tap0_ms / 1000, sos)
        print(r.relpath, json.dumps(res["files"][r.relpath]), flush=True)
    summ = {}
    for name in CANDIDATES:
        for grp, keys in (("CLA", ["CLA"]), ("Tomkin", [k for k in res["files"] if k != "CLA"])):
            v = [res["files"][k][name] for k in keys]
            summ[f"{grp} {name}"] = {"sum_fc_hz_median": round(float(np.median([q["sum_fc_hz"] for q in v])), 2),
                                    "sum_fc_hz_iqr": [round(float(np.percentile([q["sum_fc_hz"] for q in v], p)), 2) for p in (25, 75)],
                                    "free_hz_median": [round(float(np.median([q["free_hz"][i] for q in v])), 2) for i in range(len(v[0]["free_hz"]))],
                                    "R_err_db_median": round(float(np.median([q["R_err_db"] for q in v])), 1),
                                    "L_heldout_err_db_median": round(float(np.median([q["L_heldout_err_db"] for q in v])), 1),
                                    "n_params_R": v[0]["n_params_R"], "n_files": len(v)}
    nz = [res["files"][k]["noise_R_db"] for k in res["files"] if res["files"][k]["noise_R_db"] is not None]
    summ["noise_R_db_median"] = round(float(np.median(nz)), 1) if nz else None
    res["summary"] = summ
    print(json.dumps(summ, indent=1))
    json.dump(res, open(os.path.join(ROOT, "analysis", "analog", a.out), "w"), indent=1)


if __name__ == "__main__":
    main()
