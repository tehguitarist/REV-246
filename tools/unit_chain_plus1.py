"""§37.2: "CLA 10 + 1": CLA's 10 bq.json numbers kept exactly plus one RBJ peaking section (3 numbers) fitted to the units'
measured plain transfers (analysis/unit_transfer/measure.json) on 140-600 Hz only, a free delay per unit (the comparison
refits the shift). Candidates: P_all (CLA and the five Altiverb presets, CLA weighted 1/2 and the Altiverb set 1/2),
P_alt (five Altiverb presets), P_cla (CLA), P_ho (CLA + Altiverb RT 1.6, 1/2 each; the other four presets held out).
    .venv/bin/python tools/unit_chain_plus1.py       # -> analysis/unit_chain/plus1.json and the footprint table
"""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "harness"))
import unit_transfer as ut  # noqa: E402
import unit_chain as uc  # noqa: E402

ROOT = os.path.dirname(HERE)
R = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "measure.json")))
X10 = json.load(open(os.path.join(ROOT, "analysis", "unit_transfer", "bq.json")))["CLA"]["x"]
f6 = np.array(R["sixth"])
sel = (f6 >= 140) & (f6 <= 600)
ALT = [u for u in R["irs"] if u.startswith("Altiverb")]
T = {u: ut.plain_T(R, u) / ut.plain_T(R, u)[-1] for u in R["irs"]}
LO, HI = [200.0, -6.0, 0.5], [700.0, 6.0, 8.0]


def fit(units, weights):
    def res(v):
        out = []
        for u, w in zip(units, weights):
            d = np.exp(-2j * np.pi * f6[sel] * v[3 + units.index(u)] / ut.FS)
            r = (ut.bq_resp(X10, f6[sel]) * peak(v[:3], f6[sel]) * d - T[u][sel]) * np.sqrt(w)
            out += [r.real, r.imag]
        return np.concatenate(out)
    best = None
    for x0 in ([380.0, -0.6, 2.9], [300.0, 0.5, 1.0], [450.0, -0.3, 4.0], [250.0, 0.5, 1.0]):
        s_ = least_squares(res, np.r_[x0, np.zeros(len(units))], bounds=(LO + [-2] * len(units), HI + [2] * len(units)), x_scale="jac", max_nfev=2000)
        if best is None or s_.cost < best.cost:
            best = s_
    return best.x[:3].tolist(), float(np.sqrt(np.mean(best.fun ** 2)))


def peak(p, f):
    s_ = ut.peak_sos(*p)
    z = np.exp(-2j * np.pi * np.asarray(f) / ut.FS)
    return (s_[0] + s_[1] * z + s_[2] * z ** 2) / (s_[3] + s_[4] * z + s_[5] * z ** 2)


def resp(p, f):
    return ut.bq_resp(X10, f) * peak(p, f)


if __name__ == "__main__":
    cand = {"P_all": (["CLA"] + ALT, [0.5] + [0.5 / len(ALT)] * len(ALT)),
            "P_alt": (ALT, [1 / len(ALT)] * len(ALT)),
            "P_cla": (["CLA"], [1.0]),
            "P_ho": (["CLA", "Altiverb RT 1.6"], [0.5, 0.5])}
    out = {}
    for k, (us, ws) in cand.items():
        p, r = fit(us, ws)
        base = np.sqrt(np.mean([w * np.mean(np.abs(ut.bq_resp(X10, f6[sel]) - T[u][sel]) ** 2) for u, w in zip(us, ws)]))
        out[k] = dict(x=p, rms=r)
        print(k, [round(v, 4) for v in p], "rms 140-600", round(r, 4), "(10 numbers, no delay:", round(float(base), 4), ")")
    json.dump(out, open(os.path.join(ROOT, "analysis", "unit_chain", "plus1.json"), "w"), indent=1)
    f = np.geomspace(20, 13000, 2000)
    h0 = ut.bq_resp(X10, f)
    print("\nfootprint vs shipped CLA 10: max |dB| / max phase deg per band")
    for k in out:
        print(k)
        for lo, hi in [(20, 40), (40, 80), (80, 160), (160, 270), (270, 400), (400, 800), (800, 1600), (1600, 3200), (3200, 6400), (6400, 13000)]:
            m = (f >= lo) & (f < hi)
            d = peak(out[k]["x"], f[m])
            print(f"  {lo:5d}-{hi:5d}: {np.abs(20 * np.log10(abs(d))).max():.3f} dB  {np.abs(np.degrees(np.angle(d))).max():.2f} deg")
