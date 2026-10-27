"""Does the scoring path (EMT_XCORR in q2_injection.final_G) apply the per-loop correction the way xfit / xprobe's columns
assume? One setting, baseline clock and gain kept: residual energy per window and band from (a) the model with the
corrected G and (b) the prediction e - A h from the shifted loop contributions.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/xprobe_apply_check.py 1.0_0.8_0.4
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sml_fullband as sf
from sml_fullband import q1, q2, pc, np as _np, sosfiltfilt

key = sys.argv[1]
xc = json.load(open(os.path.join(sf.OUT, "xcorr.json")))
info, fade, _, _, eps0 = q2.q1_tables()
ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
v, law, _ = q2.final_injection(compact=pc.COMPACT)
cells = json.load(open(os.path.join(pc.CV, "cells.json")))
g_ = info[key]
fold_f = json.load(open(os.path.join(pc.CV, f"fold_{g_['rt']}.json")))
lab = dict(lf=g_["lf"], rt=g_["rt"], hf=g_["hf"])
b230 = ap["steps"][ap["settings"][key]["step"]]
caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(g_) if (key, a_) not in sf.DEFECTS]
q1.EXT = 1700
assert q2.XCORR is None
G0 = pc.setting_G(lab, b230, law, fold_f, cells.get(key))
q2.XCORR = xc
G1 = pc.setting_G(lab, b230, law, fold_f, cells.get(key))
q2.XCORR = None
fold = next((f for f, d in xc["folds"].items() if lab["rt"] not in d["train_rt"]), "all")
h = np.array(xc["folds"][fold]["h"])
lags = xc["lags"]
print(key, "fold used", fold, "|G1.R/G0.R - 1| max per loop", [float(np.max(np.abs(b / np.where(a == 0, 1, a) - 1))) for a, b in zip(G0["R"], G1["R"])])
J = q1.Joint([dict(G0, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
cl, f_, xs = J.split(v)
F_, X_ = J._spectra(f_, xs, G0["N"])
c = caps[0]
y, mdl, g, te = q2.fit_capture(G0, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
n = len(y)
N = 1 << int(np.ceil(np.log2(n + 8192)))
Pc, Pr = sf.pipeline(N, False)[:, None] * np.ones((1, 2)), sf.pipeline(N, True)[:, None] * np.ones((1, 2))
r0, r1 = J.reverb(G0, F_, X_)[:n], J.reverb(G1, F_, X_)[:n]
mf0 = sf.apply_f(mdl - r0, Pc) + sf.apply_f(r0, Pr)
mf1 = sf.apply_f(mdl - r0, Pc) + sf.apply_f(r1, Pr)
print("baseline full_model matches sf.full_model:", float(np.max(np.abs(mf0 - sf.full_model(J, v, G0, n, c["amp"], mdl)))))
e0, e1 = y - mf0 * g, y - mf1 * g
P = sf.pipeline(N, True)
for ch in (0, 1):
    ri = [sf.apply_f(np.fft.irfft(X_[ch, i] * G0["R"][i], G0["N"])[:n, None] * np.ones((1, 2)), Pr)[:, 0] for i in range(4)]
    A = np.stack([np.roll(r, lag) * g[:, ch] for r in ri for lag in lags], 1)
    pred = e0[:, ch] - A @ h.reshape(-1)
    print(f"ch {ch}: sum of per-loop contributions vs reverb (rms ratio) {np.std(sum(ri) - sf.apply_f(r0, Pr)[:, ch]) / np.std(sf.apply_f(r0, Pr)[:, ch]):.2e}; "
          f"actual change of reverb (rms of (mf1-mf0)*g) {np.std((mf1 - mf0)[:, ch] * g[:, ch]):.3e} vs A h {np.std(A @ h.reshape(-1)):.3e}; "
          f"rms diff {np.std((mf1 - mf0)[:, ch] * g[:, ch] + A @ h.reshape(-1)):.3e}")
    for a_ms, b_ms in ((0, 50), (50, 150), (150, 300)):
        ix = sf.win_idx(a_ms, b_ms, te, n)
        if ix is None:
            continue
        i0, i1 = ix
        row = []
        for name, sos in (("full", None), ("<400", sf.XP_SOS[0]), ("400-3200", sf.XP_SOS[1]), (">3200", sf.XP_SOS[2])):
            f = (lambda x: x) if sos is None else (lambda x: sosfiltfilt(sos, x))
            E0, E1, Ep = (np.sum(f(z)[i0:i1] ** 2) for z in (e0[:, ch], e1[:, ch], pred))
            row.append(f"{name} actual {10 * np.log10(E1 / E0):+.2f} / predicted {10 * np.log10(Ep / E0):+.2f}")
        print(f"   {a_ms}-{b_ms} ms: " + "; ".join(row))
