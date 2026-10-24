"""Capped HF pole: does the unit move the missing HF loss into the broadband gain? Per loop of each short-RT, HF < 1
setting: the law's uncapped pole p_law, the baseline's pole p_b, and the shortfall factor k = |op(p_law)| / |op(p_b)| at
the law's f_ref, against the broadband loss read from the free correction (poles.json: a / a_b).
model_s_results §25.14.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/capcomp.py"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import dense_mid as dm
q2, q1, ls = dm.q2, dm.q1, dm.ls
law = q2.final_injection(compact=ls.lb.lf_.COMPACT)[1]
law["lf"] = {float(k): x for k, x in law["lf"].items()}
fref = law["f_ref"]; FS = q1.FS
rows = json.load(open(os.path.join(ROOT, "analysis/loop_fir/chain3_current_kc8_read/poles.json")))
wr = 2 * np.pi * np.atleast_1d(fref) / FS if np.ndim(fref) == 0 else None
print("f_ref", fref)
def mag(p, f):
    w = 2 * np.pi * f / FS
    return abs((1 - p) / (1 - p * np.exp(-1j * w)))
out = []
for r in rows:
    L = law["lf"][r["lf"]]
    kk = -3 * r["P"] / FS
    p_law = float(q1.onepole_for(np.array(10 ** (kk / r["rt"] * L["k"] * (1 / r["hf"] - 1))), law["f_ref"]))
    k = mag(p_law, fref) / mag(r["pb"], fref) if p_law > r["pb"] else 1.0
    out.append((r, p_law, k))
print("setting loop | p_law  p_base  p_read | predicted broadband factor k-1 | read a/a_b-1 | one-pole fit dB")
for r, pl, k in out:
    if pl > 0.80 or r["a"] / r["ab"] - 1 < -0.01:
        print(f"{r['key']} L{r['loop']} | {pl:.4f} {r['pb']:.4f} {r['p']:.4f} | {k - 1:+.4f} | {r['a'] / r['ab'] - 1:+.4f} | {r['fit_db']:+.0f}")
