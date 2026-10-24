"""Each loop's HF pole read from the current-baseline free corrections: the baseline's HF one-pole H_b (its loop response
over the same with HF 1.0), the target H_b C, a one-pole a (1 - p) / (1 - p z^-1) fitted over 177 Hz - 5.6 kHz.
model_s_results §25.14.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/poles.py    # needs loop_fir.py --current --kc 8 fits; -> analysis/loop_fir/chain3_current_kc8_read/poles.json"""
import json, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.optimize import least_squares
import dense_mid as dm
q2, q1, ls = dm.q2, dm.q1, dm.ls
FS = q1.FS
f = np.geomspace(177, 5600, 120); w = 2 * np.pi * f / FS; z1 = np.exp(-1j * w)
op = lambda p: (1 - p) / (1 - p * z1)
law = q2.final_injection(compact=ls.lb.lf_.COMPACT)[1]
info = q2.q1_tables()[0]
rows = []
for p in sorted(glob.glob(os.path.join(ROOT, "analysis/loop_fir/chain3_current_kc8/*.json"))):
    d = json.load(open(p)); key = d["key"]
    if d["hf"] >= 1.0:
        continue
    t, lfdc, kw = ls.base_kw(key, ls.HF_CELLS_JSON)
    L = json.loads(json.dumps(law)); L["lf"] = {float(k): x for k, x in L["lf"].items()}
    for k_, x in lfdc.items():
        L["lf"][float(k_)]["dc"] = float(x)
    lab = dict(lf=d["lf"], rt=d["rt"], hf=d["hf"])
    Gb = q2.final_Gs(lab, L, w, lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=0.0, **kw)
    G1 = q2.final_Gs(dict(lab, hf=1.0), L, w, lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=0.0, **kw)
    lags = np.array(d["lags"])
    for i in range(4):
        Hb = Gb[i] / G1[i]
        pb = least_squares(lambda u: np.r_[(u[0] * op(u[1]) - Hb).real, (u[0] * op(u[1]) - Hb).imag], [1, 0.7], bounds=([0.5, 0], [1.5, 0.999])).x
        C = np.exp(-1j * np.outer(w, lags)) @ np.array(d["fir"][i])
        T = Hb * C
        r = least_squares(lambda u: np.r_[(u[0] * op(u[1]) - T).real, (u[0] * op(u[1]) - T).imag], pb, bounds=([0.5, 0], [1.5, 0.999]))
        mis = np.sqrt(np.mean(np.abs(r.x[0] * op(r.x[1]) - T) ** 2)) / np.sqrt(np.mean(np.abs(T - Hb) ** 2) + 1e-30)
        rows.append(dict(key=key, lf=d["lf"], rt=d["rt"], hf=d["hf"], loop=i, P=int(q1.P[i]), pb=float(pb[1]), ab=float(pb[0]), p=float(r.x[1]), a=float(r.x[0]), fit_db=float(20 * np.log10(mis))))
json.dump(rows, open(os.path.join(ROOT, "analysis/loop_fir/chain3_current_kc8_read/poles.json"), "w"))
import collections
print("cell (HF, RT): per LF, baseline pole -> read pole (median over loops), dp, one-pole residual re the correction (dB)")
C_ = collections.defaultdict(list)
for r in rows: C_[(r["hf"], r["rt"], r["lf"])].append(r)
for k in sorted(C_):
    rr = C_[k]
    print(f"HF {k[0]} RT {k[1]} LF {k[2]:<3}: " + "  ".join(f"L{r['loop']} {r['pb']:.4f}->{r['p']:.4f} ({r['p']-r['pb']:+.4f}, a {r['a']/r['ab']-1:+.4f}, {r['fit_db']:+.0f})" for r in rr))
