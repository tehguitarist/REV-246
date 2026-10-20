"""Compare the simulated analogue chain (input filter x output filter x one-pole lag) with the
reference kernel (ref/kernel_eugtone_96k.npz) and the manual's THROUGH figures.

    .venv/bin/python tools/analog/compare_kernel.py
Writes analysis/analog/{summary.json,chain.png}.
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

z = np.load(os.path.join(ROOT, "ref", "kernel_eugtone_96k.npz"))
fs = float(z["fs"]); N = 1 << 16
mag = 20 * np.log10(np.abs(np.fft.rfft(z["g"], N)) + 1e-30)
fr = np.fft.rfftfreq(N, 1 / fs)
kernel = lambda f: np.interp(f, fr, mag) - np.interp(1e3, fr, mag)


def chain(f, gbw):
    f = np.asarray(f, float)
    out = 0
    for net in (INPUT_FILTER, OUTPUT_FILTER):
        h = 20 * np.log10(np.abs(response(net, np.r_[1e3, f], gbw)))
        out = out + (h[1:] - h[0])
    return out


lag = lambda f, fc: -10 * np.log10(1 + (np.asarray(f) / fc) ** 2) + 10 * np.log10(1 + (1e3 / fc) ** 2)
Ffit = np.arange(500, 11000, 250.0)
Fall = np.r_[np.arange(500, 12000, 250.0), np.arange(12000, 14001, 100.0)]
res, curves = {}, {}
for name, gbw in (("ideal", None), ("GBW 3 MHz", 3e6), ("GBW 5 MHz", 5e6)):
    r = least_squares(lambda p: chain(Ffit, gbw) + lag(Ffit, p[0]) - kernel(Ffit), [6000.0])
    fc = float(r.x[0])
    curves[name] = chain(Fall, gbw) + lag(Fall, fc)
    thr = {f"{f / 1e3:g} kHz": round(float(chain([f], gbw)[0] + lag([f], fc)[0]), 2) for f in (10e3, 12e3)}
    res[name] = {"lag_corner_hz": round(fc), "rms_db_0p5_to_11k": round(float(np.sqrt(np.mean(r.fun ** 2))), 3), "through_like": thr}
res["kernel"] = {f"{f / 1e3:g} kHz": round(float(kernel([f])[0]), 2) for f in (10e3, 12e3)}
res["manual_through_spec_db"] = {"10 kHz": "-0.5 to -1", "12 kHz": -2.5}
os.makedirs(os.path.join(ROOT, "analysis", "analog"), exist_ok=True)
json.dump(res, open(os.path.join(ROOT, "analysis", "analog", "summary.json"), "w"), indent=1)
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
plt.tight_layout(); plt.savefig(os.path.join(ROOT, "analysis", "analog", "chain.png"), dpi=110)
