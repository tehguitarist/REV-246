"""§37.x: ONE low-band cascade for CLA and the Altiverb (equal weight per unit), in the 10-number (bq.json) and 13-number
(unit_chain.py) forms. Writes analysis/unit_chain/joint.json.

  J_T  fitted to the two units' measured plain transfers (analysis/unit_transfer/measure.json: CLA and Altiverb RT 1.6,
       mean of L and R, sixth octaves 40 Hz-1.6 kHz, normalised at 1.6 kHz, held to 1 at 2-10 kHz), as unit_transfer.fit_bq
       does for one unit. Data window 50-600 ms; the other four Altiverb presets are held out.
  J_0  fitted to the two units' 13-number cascades fitted on 0-50 ms only (unit_chain s.json chain@0-50), complex response on
       20 Hz-2 kHz, log spaced: the data windows from 50 ms are held out (a surrogate, the cascades being the data's summary).

    .venv/bin/python tools/unit_chain_joint.py
"""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "harness"))
import unit_transfer as ut  # noqa: E402
import unit_chain as uc  # noqa: E402

ROOT = os.path.dirname(HERE)
FS = ut.FS
R = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "measure.json")))
S = json.load(open(os.path.join(ROOT, "analysis", "unit_chain", "s.json")))["irs"]
BQ = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))
f6 = np.array(R["sixth"])
sel = (f6 >= 40) & (f6 <= 1600)
fa = np.array([2000, 3000, 4000, 6000, 8000, 10000.0])
U = ["CLA", "Altiverb RT 1.6"]
T = {u: ut.plain_T(R, u) / ut.plain_T(R, u)[-1] for u in U}


def resp(x, f):
    return ut.bq_resp(x, f) if len(x) == 10 else uc.lb_resp(x, f)


def fit_T(n, x0):
    lo, hi = (ut.LO, ut.HI) if n == 10 else (uc.LO3, uc.HI3)

    def res(v):
        out = []
        for u in U:                                    # a free delay per unit (the comparison refits the shift), equal weight
            d = np.exp(-2j * np.pi * f6[sel] * v[len(x0) + U.index(u)] / FS)
            r = resp(v[:len(x0)], f6[sel]) * d - T[u][sel]
            ra = resp(v[:len(x0)], fa) - 1
            out += [r.real, r.imag, ra.real, ra.imag]
        return np.concatenate(out)
    s_ = least_squares(res, np.r_[np.clip(x0, lo, hi), 0, 0], bounds=(lo + [-2, -2], hi + [2, 2]), x_scale="jac", max_nfev=2000)
    nin = 2 * sel.sum() * 2
    return s_.x[:len(x0)].tolist(), float(np.sqrt(np.mean(s_.fun ** 2)))


def fit_0(n, x0):
    lo, hi = (ut.LO, ut.HI) if n == 10 else (uc.LO3, uc.HI3)
    f = np.geomspace(20, 2000, 120)
    tg = [uc.lb_resp(S[u]["chain@0-50"]["chain"]["bq"], f) for u in U]

    def res(v):
        h = resp(v, f)
        return np.concatenate([np.r_[(h - t).real, (h - t).imag] for t in tg])
    s_ = least_squares(res, np.clip(x0, lo, hi), bounds=(lo, hi), x_scale="jac", max_nfev=2000)
    return s_.x.tolist(), float(np.sqrt(np.mean(s_.fun ** 2)))


if __name__ == "__main__":
    out = {}
    xc10 = BQ["CLA"]["x"]
    xc13 = S["CLA"]["chain@150-300"]["chain"]["bq"]
    for tag, fn in (("J_T", fit_T), ("J_0", fit_0)):
        for n, x0 in ((10, xc10), (13, xc13)):
            best = None
            for start in (x0, BQ["Altiverb RT 1.6"]["x"] if n == 10 else S["Altiverb RT 1.6"]["chain@150-300"]["chain"]["bq"]):
                x, r = fn(n, start)
                if best is None or r < best[1]:
                    best = (x, r)
            out[f"{tag}_{n}"] = dict(x=best[0], rms=best[1])
            print(tag, n, round(best[1], 4), [round(v, 3) for v in best[0]])
    json.dump(out, open(os.path.join(ROOT, "analysis", "unit_chain", "joint.json"), "w"), indent=1)
