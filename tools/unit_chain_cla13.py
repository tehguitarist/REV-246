"""§37.x: a 13-number CLA-only cascade fitted to CLA's measured plain transfer (the method of unit_transfer.fit_bq with
unit_chain's third peak; 40 Hz-1.6 kHz sixth octaves, held to 1 at 2-10 kHz), started from bq.json's 10 numbers plus a
third peak, and its fullband score on CLA alone beside the 10 numbers. -> analysis/unit_chain/cla13.json, cla13_run/.
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/unit_chain_cla13.py"""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "harness"))
import unit_transfer as ut  # noqa: E402
import unit_chain as uc  # noqa: E402
import unit_chain_fullband as fb  # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = fb.OUT
R = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "measure.json")))
f6 = np.array(R["sixth"]); sel = (f6 >= 40) & (f6 <= 1600)
T = ut.plain_T(R, "CLA"); T = T / T[-1]
fa = np.array([2000, 3000, 4000, 6000, 8000, 10000.0])
x10 = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))["CLA"]["x"]


def res(v):
    r = uc.lb_resp(v[1:], f6[sel]) * np.exp(-2j * np.pi * f6[sel] * v[0] / ut.FS) - T[sel]
    ra = uc.lb_resp(v[1:], fa) - 1
    return np.r_[r.real, r.imag, ra.real, ra.imag]


best = None
for x3 in ([420.0, -0.2, 2.0], [390.0, -0.5, 3.0], [500.0, 0.3, 1.0]):
    s_ = least_squares(res, np.r_[0.0, x10, x3], bounds=([-2] + uc.LO3, [2] + uc.HI3) if False else ([-2] + ut.LO + [350.0, -3.0, 0.5], [2] + ut.HI + [700.0, 3.0, 8.0]), x_scale="jac", max_nfev=2000)
    rms = float(np.sqrt(np.mean(s_.fun[:2 * sel.sum()] ** 2)))
    if best is None or rms < best[1]:
        best = (s_.x[1:].tolist(), rms)
print("CLA 13 numbers fitted to CLA's transfer: rms", round(best[1], 4), "(10 numbers:", round(json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))["CLA"]["rms"], 4), ")")
print([round(v, 4) for v in best[0]])
json.dump(dict(CLA13T=dict(x=best[0], rms=best[1])), open(os.path.join(OUT, "cla13.json"), "w"), indent=1)

import altiverb_predict as av  # noqa: E402
av.PRESETS = {}
_prev = ut.variant_H
ut.variant_H = lambda vn: (lambda f: np.repeat(uc.lb_resp(best[0], f)[:, None], 2, 1)) if vn == "c13T" else _prev(vn)
SUB = os.path.join(OUT, "cla13_run"); os.makedirs(SUB, exist_ok=True)
ut.OUT = fb.OUT = SUB
V = ["plain", "bq_CLA", "c13T"]
ut.fullband_stage(V)
os.replace(os.path.join(SUB, "fullband.json"), os.path.join(SUB, "cla13.json"))
r = json.load(open(os.path.join(SUB, "cla13.json")))["irs"]["CLA"]
print("CLA, null R per band (dB): window 50-150 / 150-300 / 300-600 / 600-1000, bq_CLA then c13T, and delta vs plain")
for b in json.load(open(os.path.join(SUB, "cla13.json")))["bands"]:
    W = ["50-150", "150-300", "300-600", "600-1000"]
    print(f"{b:12s}", " ".join(f"{r['plain']['null'][b][w][1]:6.1f}" for w in W), "|", " ".join(f"{r['bq_CLA']['null'][b][w][1]:6.1f}" for w in W), "|", " ".join(f"{r['c13T']['null'][b][w][1]:6.1f}" for w in W))
