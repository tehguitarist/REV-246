"""Compare the simulated analogue chain (input filter x output filter x one-pole lag) with the
reference kernel (ref/kernel_eugtone_96k.npz) and the manual's THROUGH figures.

    .venv/bin/python tools/analog/compare_kernel.py [--kernel <npz> --tag <name>] [--input-sh] [--fit-ti] [--skirt-shift]
Writes analysis/analog/{summary,chain}[_<tag>].{json,png}. Default kernel: the CLA reference kernel.
Options (all off by default, so the default chain is unchanged):
  --input-sh     add the input integrate/hold stage (sample_hold.input_sh, 9.1 us window from p. 39)
  --fit-ti       fit the input window length freely (implies --input-sh)
  --skirt-shift  fit a common frequency scale on both filters (component tolerance)
  --ron R        input S&H switch on-resistance in ohms (default 0)
  --zoh          add the output hold's zero-order-hold sinc (staircase between S&H windows)
  --fit-hi F     fit the S&H (and options) on 0.5 kHz..F instead of 0.5-11 kHz (rms fields keep their names)
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from mna import response  # noqa: E402
from netlists import INPUT_FILTER, OUTPUT_FILTER  # noqa: E402
from sample_hold import input_sh, zoh_db  # noqa: E402

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--kernel", default=os.path.join("ref", "kernel_eugtone_96k.npz"))
ap.add_argument("--tag", default="")
ap.add_argument("--input-sh", action="store_true")
ap.add_argument("--fit-ti", action="store_true")
ap.add_argument("--skirt-shift", action="store_true")
ap.add_argument("--ron", type=float, default=0.0, help="input S&H switch on-resistance (ohm)")
ap.add_argument("--zoh", action="store_true")
ap.add_argument("--fit-hi", type=float, default=11000.0, help="upper end (Hz) of the S&H/option fit; 12300 includes the skirt")
args = ap.parse_args()
args.input_sh = args.input_sh or args.fit_ti
TAG = f"_{args.tag}" if args.tag else ""
z = np.load(os.path.join(ROOT, args.kernel))
fs = float(z["fs"]); N = 1 << 16
mag = 20 * np.log10(np.abs(np.fft.rfft(z["g"], N)) + 1e-30)
fr = np.fft.rfftfreq(N, 1 / fs)
kernel = lambda f: np.interp(f, fr, mag) - np.interp(1e3, fr, mag)


def chain(f, gbw, scale=1.0):
    """Both filters in dB re 1 kHz; scale > 1 moves their response up in frequency by that factor."""
    f = np.asarray(f, float) / scale
    out = 0
    for net in (INPUT_FILTER, OUTPUT_FILTER):
        h = 20 * np.log10(np.abs(response(net, np.r_[1e3 / scale, f], gbw)))
        out = out + (h[1:] - h[0])
    return out


lag = lambda f, fc: -10 * np.log10(1 + (np.asarray(f) / fc) ** 2) + 10 * np.log10(1 + (1e3 / fc) ** 2)
FS_DEV = 19.2e6 / 728


def iir_lag(f, pole):
    """Output S&H charging through R202 + P202 into C205 (service manual p. 48) acts as y[n] = y[n-1] + (1-p)(x[n] - y[n-1]) at the device rate."""
    w = 2 * np.pi * np.asarray(f, float) / FS_DEV
    h = (1 - pole) / np.abs(1 - pole * np.exp(-1j * w))
    h1 = (1 - pole) / abs(1 - pole * np.exp(-2j * np.pi * 1e3 / FS_DEV))
    return 20 * np.log10(h / h1)


def in_sh(f, ti=None):
    h = 20 * np.log10(np.abs(input_sh(np.r_[1e3, np.asarray(f, float)], ti, args.ron)[0]))
    return h[1:] - h[0]


BANDS = ((500, 1000), (1000, 2000), (2000, 4000), (4000, 8000), (8000, 11000), (11000, 12300))


def bands(err_f, f):
    """rms dB error per band; the last band is the 11-12.3 kHz skirt (not in the fit range)."""
    return {f"{a / 1e3:g}-{b / 1e3:g} kHz": round(float(np.sqrt(np.mean(err_f[(f >= a) & (f < b)] ** 2))), 3) for a, b in BANDS}


Ffit = np.arange(500, 11000, 250.0)
Fopt = Ffit if args.fit_hi <= 11000 else np.r_[Ffit, np.arange(11000, args.fit_hi + 1, 50.0)]
Fband = np.r_[np.arange(500, 11000, 100.0), np.arange(11000, 12301, 25.0)]
Fall = np.r_[np.arange(500, 12000, 250.0), np.arange(12000, 14001, 100.0)]
res, curves = {}, {}
for name, gbw in (("ideal", None), ("GBW 3 MHz", 3e6), ("GBW 5 MHz", 5e6)):
    r = least_squares(lambda p: chain(Ffit, gbw) + lag(Ffit, p[0]) - kernel(Ffit), [6000.0])
    fc = float(r.x[0])
    ri = least_squares(lambda p: chain(Fopt, gbw) + iir_lag(Fopt, p[0]) - kernel(Fopt), [0.3], bounds=(0, 0.95))
    if args.input_sh or args.skirt_shift or args.zoh:
        # pole p, then [input window (us)], then [filter frequency scale]
        def unpack(q):
            q = list(q); pole = q.pop(0)
            ti = q.pop(0) * 1e-6 if args.fit_ti else None
            sc = q.pop(0) if args.skirt_shift else 1.0
            return pole, ti, sc
        def model(f, q):
            pole, ti, sc = unpack(q)
            return (chain(f, gbw, sc) + iir_lag(f, pole) + (in_sh(f, ti) if args.input_sh else 0)
                    + (zoh_db(f) if args.zoh else 0))
        x0, lo, hi = [0.3], [0], [0.95]
        if args.fit_ti: x0, lo, hi = x0 + [9.1], lo + [1.0], hi + [30.0]
        if args.skirt_shift: x0, lo, hi = x0 + [1.0], lo + [0.97], hi + [1.03]
        ri = least_squares(lambda q: model(Fopt, q) - kernel(Fopt), x0, bounds=(lo, hi))
    else:
        model = lambda f, q: chain(f, gbw) + iir_lag(f, q[0])
    curves[name] = chain(Fall, gbw) + lag(Fall, fc)
    thr = {f"{f / 1e3:g} kHz": round(float(chain([f], gbw)[0] + lag([f], fc)[0]), 2) for f in (10e3, 12e3)}
    res[name] = {"lag_corner_hz": round(fc), "rms_db_0p5_to_11k": round(float(np.sqrt(np.mean(r.fun ** 2))), 3), "through_like": thr,
                 "sandh_pole": round(float(ri.x[0]), 3),
                 "sandh_rms_db_0p5_to_11k": round(float(np.sqrt(np.mean((model(Ffit, ri.x) - kernel(Ffit)) ** 2))), 3),
                 "sandh_bands_rms_db": bands(model(Fband, ri.x) - kernel(Fband), Fband)}
    if args.fit_ti: res[name]["input_window_us"] = round(float(ri.x[1]), 2)
    if args.skirt_shift: res[name]["skirt_scale"] = round(float(ri.x[-1]), 5)
    if args.input_sh or args.skirt_shift or args.zoh: curves[name + " + S&H fit"] = model(Fall, ri.x)
res["kernel"] = {f"{f / 1e3:g} kHz": round(float(kernel([f])[0]), 2) for f in (10e3, 12e3)}
res["manual_through_spec_db"] = {"10 kHz": "-0.5 to -1", "12 kHz": -2.5}
res["options"] = {"input_sh": args.input_sh, "fit_ti": args.fit_ti, "skirt_shift": args.skirt_shift, "zoh": args.zoh, "ron_ohm": args.ron, "fit_hi_hz": args.fit_hi}
os.makedirs(os.path.join(ROOT, "analysis", "analog"), exist_ok=True)
json.dump(res, open(os.path.join(ROOT, "analysis", "analog", f"summary{TAG}.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for a, lim in zip(ax, ((0, 11.5), (11.5, 14))):
    m = (Fall / 1e3 >= lim[0]) & (Fall / 1e3 <= lim[1])
    a.plot(Fall[m] / 1e3, kernel(Fall)[m], "k", lw=2, label="reference kernel")
    for n, c in curves.items(): a.plot(Fall[m] / 1e3, c[m], label=n)
    a.set_xlabel("kHz"); a.set_ylabel("dB re 1 kHz"); a.grid(alpha=.3)
ax[0].set_ylim(-8, 4); ax[1].set_ylim(-90, 5); ax[0].legend()
plt.tight_layout(); plt.savefig(os.path.join(ROOT, "analysis", "analog", f"chain{TAG}.png"), dpi=110)
