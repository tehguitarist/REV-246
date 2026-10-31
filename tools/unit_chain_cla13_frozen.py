"""§37.3: the frozen scoring of unit_chain_plus1_score.py (each capture's clock, shift, gain, 14-number FIR as fitted for the shipped
CLA 10) for §37.1's CLA-only 13-number cascade (analysis/unit_chain/cla13.json), on CLA and the five Altiverb presets; footprint.
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/unit_chain_cla13_frozen.py
    .venv/bin/python tools/unit_chain_cla13_frozen.py --report"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import unit_chain_plus1_score as sc  # noqa: E402
import unit_chain as uc  # noqa: E402
import unit_transfer as ut  # noqa: E402

X13 = json.load(open(os.path.join(sc.OUT, "cla13.json")))["CLA13T"]["x"]
sc.P = {"C13T": {"x": X13}}
sc.VARS = ["plain", "bq_CLA", "C13T"]
_h = sc.H_of
sc.H_of = lambda vn: (lambda ff: np.repeat(uc.lb_resp(X13, ff)[:, None], 2, 1)) if vn == "C13T" else _h(vn)
_dump = json.dump
sc.OUTNAME = "cla13_score.json"

if __name__ == "__main__":
    if "--report" in sys.argv:
        f = np.geomspace(20, 13000, 3000)
        d = uc.lb_resp(X13, f) / ut.bq_resp(sc.p1.X10, f)
        print("footprint vs shipped CLA 10: max |dB| / max |phase deg| per band")
        for lo, hi in [(20, 40), (40, 80), (80, 160), (160, 270), (270, 400), (400, 800), (800, 1600), (1600, 3200), (3200, 6400), (6400, 13000)]:
            m = (f >= lo) & (f < hi)
            print(f"  {lo:5d}-{hi:5d}: {np.abs(20 * np.log10(abs(d[m]))).max():.3f} dB  {np.abs(np.degrees(np.angle(d[m]))).max():.2f} deg")
        os.replace(os.path.join(sc.OUT, "cla13_score.json"), os.path.join(sc.OUT, "plus1_score.json")) if False else None
    sc.report() if "--report" in sys.argv else sc.run()
