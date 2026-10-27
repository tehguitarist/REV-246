"""Per-capture HF chain-correction table for the bench (NL_HF=1): z(key) = z0 + s u at 9 knots (10.5 ... 12.8 kHz), a complex log-gain
applied to the model above 10 kHz (unity below 10 kHz and above 12.8 kHz). Shape (z0, u) = rank-1 fit over the capture runs
(fork D, d_perrun.py), fitted on the OTHER runs for a capture inside a run (leave-one-run-out); the scalar s from the R-channel front
taps of the 15 captures nearest in time (the capture's own front taps excluded). Captures without a timestamp: the all-run mean.
Writes analysis/nl/hf_corr.json."""
import os, sys, json, pickle
import numpy as np
OUT = "/Users/leighpierce/Git/246/analysis/nl_hfdrift"
KI = [1, 2, 3, 4, 5, 6, 7, 8, 9]
KN = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
rows = [r for r in pickle.load(open(OUT + "/d_front_all.pkl", "rb")) if r["t"] is not None]
rows.sort(key=lambda r: r["t"]); t = np.array([r["t"] for r in rows])
runs = [np.array(ix) for ix in np.split(np.arange(len(rows)), np.where(np.diff(t) > 20 * 60)[0] + 1)]
runs = [ix for ix in runs if len(ix) >= 8]
times = json.load(open(OUT + "/d_times.json"))
def zvec(idx, ch=1):
    c = sum(rows[i][ch]["num"] for i in idx) / sum(rows[i][ch]["den"] for i in idx)
    return np.log(c[KI])
def split(z): return np.r_[z.real, z.imag]
Z = np.array([zvec(ix) for ix in runs]); nrun = np.array([len(ix) for ix in runs], float)
def fit_shape(Zs, w):
    V = np.array([split(z) for z in Zs]); z0 = (w[:, None] * V).sum(0) / w.sum()
    U, S, Vt = np.linalg.svd((V - z0) * np.sqrt(w)[:, None], full_matrices=False)
    u = Vt[0]
    # sign convention: gain part of u at 11.5 kHz negative (as in fork D's shape)
    if u[4] > 0: u = -u
    return z0, u
z0_all, u_all = fit_shape(Z, nrun)
span = [(t[ix].min() - 600, t[ix].max() + 600) for ix in runs]
out = {}
for k, tk in times.items():
    key = tuple(float(v) for v in k.split("_"))
    own = [j for j, (a, b) in enumerate(span) if a <= tk <= b]
    others = [j for j in range(len(runs)) if j not in own]
    z0, u = fit_shape(Z[others], nrun[others])
    order = sorted((i for i in range(len(rows)) if rows[i]["key"] != key), key=lambda i: abs(rows[i]["t"] - tk))[:15]
    zl = zvec(order)
    s = float(np.dot(split(zl) - z0, u) / np.dot(u, u))
    v = z0 + s * u
    out[k] = dict(s=s, run=own[0] if own else -1, z_re=v[:9].tolist(), z_im=v[9:].tolist())
mean = (nrun[:, None] * np.array([split(z) for z in Z])).sum(0) / nrun.sum()
json.dump(dict(note="HF chain correction (NonLin_program 11.6, fork D/F): complex log-gain at knots_hz per capture, exp(z) applied to the model above 10 kHz (unity <= 10 kHz and >= 12.8 kHz). z = z0 + s u; shape fitted on the other capture runs, s from the R front taps of the 15 nearest-in-time captures. Captures without a timestamp use `default`.",
               knots_hz=KN[KI].tolist(), default=dict(z_re=mean[:9].tolist(), z_im=mean[9:].tolist()), captures=out),
          open("/Users/leighpierce/Git/246/analysis/nl/hf_corr.json", "w"))
s = np.array([v["s"] for v in out.values()]); print(len(out), "captures; s mean %.2f sd %.2f range %.2f..%.2f" % (s.mean(), s.std(), s.min(), s.max()))
print("shape sign check: u gain part", np.round(u_all[:9], 3), " phase part mrad", np.round(u_all[9:] * 1e3, 1))
