"""HF filter forms against the loop_fir targets (law x free 9-tap correction) at the strong HF cells.
model_s_results §25.14.

    EMT_KERNEL=chain2 .venv/bin/python tools/loop_read/hfform.py   # old loop_fir corrections (chain2, §22 law)"""
import json, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.optimize import least_squares
import q2_injection as q2
q1 = q2.q1
law = q2.final_injection(compact="diffuser_tied_fix-0.707107_ap230new_chain2")[1]
FS = q1.FS
f = np.geomspace(100, 10000, 200); w = 2 * np.pi * f / FS
z1 = np.exp(-1j * w)
def op(p): return (1 - p) / (1 - p * z1)
FORMS = {
 "pole":      (lambda t: t[0] * op(t[1]), [1.0, 0.8]),
 "2 poles":   (lambda t: t[0] * op(t[1]) * op(t[2]), [1.0, 0.8, 0.1]),
 "pole+zero": (lambda t: t[0] * op(t[1]) * (1 - t[2] * z1) / (1 - t[2]), [1.0, 0.8, 0.05]),
 "pole+F(b)": (lambda t: t[0] * op(t[1]) * ((1 - t[2]) + t[2] * z1), [1.0, 0.8, 0.1]),
 "biquad":    (lambda t: t[0] * (1 + t[2] * z1 + t[3] * z1 ** 2) / (1 - t[1] * z1 + t[4] * z1 ** 2), [0.16, 0.84, 0.0, 0.0, 0.0]),
}
cells = [(0.2, 0.6), (0.2, 0.8), (0.4, 0.3)]
for p in sorted(glob.glob(os.path.join(ROOT, "analysis/loop_fir/chain2/*.json"))):
    d = json.load(open(p))
    if (d["hf"], d["rt"]) not in cells or d["lf"] not in (1.0, 2.0):
        continue
    lab = dict(lf=d["lf"], rt=d["rt"], hf=d["hf"])
    Gl = q2.G_capped(lab, law, w, beta=q2.FINAL_BETA)
    lab1 = dict(lab, hf=1.0); Gs = q2.G_capped(lab1, law, w, beta=q2.FINAL_BETA)   # without the HF part
    lags = np.array(d["lags"])
    out = [f"{d['key']} b {d['b230']:.2f}:"]
    for i in range(4):
        C = np.exp(-1j * np.outer(w, lags)) @ np.array(d["fir"][i])
        T = Gl[i] * C / Gs[i]                                        # target HF part (shelf and gain divided out)
        H0 = Gl[i] / Gs[i]                                           # the law's capped one-pole
        e0 = 10 * np.log10(np.mean(np.abs(T - H0) ** 2) / np.mean(np.abs(T - 1) ** 2 + 1e-30))
        res = {}
        for nm, (fn, x0) in FORMS.items():
            r = least_squares(lambda t: np.r_[(fn(t) - T).real, (fn(t) - T).imag], x0, bounds=([0.0 if len(x0) == 5 else 0.5] + [-1.99 if len(x0) == 5 else -0.99] * (len(x0) - 1), [1.5] + [1.99 if len(x0) == 5 else 0.999] * (len(x0) - 1)))
            res[nm] = (20 * np.log10(np.sqrt(np.mean(np.abs(fn(r.x) - T) ** 2)) / np.sqrt(np.mean(np.abs(T - H0) ** 2))), r.x)
        if i == 0:
            out.append(f"   law pole {(-1 if False else 0):d}")
        out.append(f"   loop {i}: law misfit re target {20*np.log10(np.sqrt(np.mean(np.abs(T-H0)**2))):.1f} dB; " +
                   "; ".join(f"{nm} {v[0]:+.1f}" for nm, v in res.items()) + "; biquad poles " + str(np.round(np.roots([1, -res["biquad"][1][1], res["biquad"][1][4]]), 3).tolist()) + " zeros " + str(np.round(np.roots([1, res["biquad"][1][2], res["biquad"][1][3]]), 3).tolist()))
    print("\n".join(o for o in out if not o.startswith("   law pole")), flush=True)
