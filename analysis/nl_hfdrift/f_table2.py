"""Per-capture HF chain-correction table v2 (fork H's best of both, analysis/nl_hfdrift/h2_est.log): the scalar from the R front taps of the
captures weighted by a Gaussian in time (tau 5 min, the capture's own front taps excluded) instead of the 15 nearest, and NO correction (unity)
where the 15 nearest front-tap estimates scatter above the 90th percentile of the library's scatter (a front-tap-only test: the capture's reverb
is never used). Shape fitted on the other capture runs, as f_table.py. Writes analysis/nl/hf_corr_g5.json (use: NL_HF=1 NL_HF_FILE=...)."""
import os, sys, json
import numpy as np
OUT = "/Users/leighpierce/Git/246/analysis/nl_hfdrift"
sys.path.insert(0, OUT)
import h_est as H
times = json.load(open(OUT + "/d_times.json"))
def own_runs(tk): return [j for j, (a, b) in enumerate(H.span) if a <= tk <= b]
def scatter(idx, w, z0, u):
    si = np.array([H.row_s(i, z0, u) for i in idx]); m = np.average(si, weights=w); return float(np.sqrt(np.average((si - m) ** 2, weights=w)))
z0a, ua = H.fit_shape(H.Z, H.nrun)
sc_all = []
for i in range(0, len(H.rows), 2):
    s_n, inf_n, idx, w = H.s_nb(H.KEYS[i], H.T[i], [], z0a, ua); sc_all.append(scatter(idx, w, z0a, ua))
X90 = float(np.percentile(sc_all, 90)); print("scatter p50/p75/p90/p95", np.round(np.percentile(sc_all, [50, 75, 90, 95]), 4), "-> threshold", round(X90, 4))
out = {}; nskip = 0
for k, tk in times.items():
    key = tuple(float(v) for v in k.split("_")); own = own_runs(tk); z0, u = H.shape_for(own)
    z, idx = H.est("gauss", dict(tau=5), tk, key, own, z0, u)
    s_n, inf_n, idx15, w15 = H.s_nb(key, tk, own, z0, u, k=15); sc = scatter(idx15, w15, z0, u)
    skip = bool(sc > X90); nskip += skip
    out[k] = dict(skip=skip, scatter=sc, z_re=z.real.tolist(), z_im=z.imag.tolist())
mean = (H.nrun[:, None] * np.array([H.split(z) for z in H.Z])).sum(0) / H.nrun.sum()
json.dump(dict(note="HF chain correction v2 (NonLin_program 11.6): complex log-gain at knots_hz per capture, exp(z) applied to the model above 10 kHz (unity <= 10 kHz and >= 12.8 kHz); `skip`: no correction. z = z0 + s u, shape from the other capture runs, s from the R front taps weighted by a 5-minute Gaussian in time; skip where the 15 nearest front-tap estimates scatter above the library's 90th percentile (%.4f). Captures without a timestamp: `default` (the all-run mean)." % X90,
               knots_hz=H.KH.tolist(), default=dict(z_re=mean[:9].tolist(), z_im=mean[9:].tolist()), captures=out),
          open("/Users/leighpierce/Git/246/analysis/nl/hf_corr_g5.json", "w"))
print(len(out), "captures,", nskip, "skipped (%.1f %%)" % (100 * nskip / len(out)))
