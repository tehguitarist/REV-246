"""§37.x: score the 13-number low-band cascade of tools/unit_chain.py (§37) against §27's 10-number one (bq.json) over the
whole band, 20 Hz-12.3 kHz, with the same fullband method (unit_transfer.fullband_stage, nothing else changed).

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 3 -- \
        .venv/bin/python tools/unit_chain_fullband.py            # -> analysis/unit_chain/fullband13.{json,log}

Variants: plain; bq_CLA / bq_Altiverb RT 1.6 (the 10 numbers, bq.json); b13_<unit>@<window> (the 13 numbers of
analysis/unit_chain/s.json "chain@<window>": 0-50 held out of every window from 50 ms, 150-300 in-sample at 150-300 ms, all in-sample).
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import unit_transfer as ut  # noqa: E402
import unit_chain as uc  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "analysis", "unit_chain")
S = json.load(open(os.path.join(OUT, "s.json")))["irs"]
VARS = ["plain", "bq_CLA", "bq_Altiverb RT 1.6"]
for u, short in (("CLA", "CLA"), ("Altiverb RT 1.6", "Alt")):
    for w in ("0-50", "150-300", "all"):
        VARS.append(f"b13_{short}@{w}")
ROOT_UT = os.path.join(os.path.dirname(HERE), "analysis", "unit_transfer")
_orig = ut.variant_H


def variant_H(vn):
    if vn.startswith("b13_"):
        short, w = vn[4:].split("@")
        x = S["CLA" if short == "CLA" else "Altiverb RT 1.6"][f"chain@{w}"]["chain"]["bq"]
        return lambda f: np.repeat(uc.lb_resp(x, f)[:, None], 2, 1)
    if vn.startswith("bq_"):
        x = json.load(open(os.path.join(ROOT_UT, "bq.json")))[vn[3:]]["x"]
        return lambda f: np.repeat(ut.bq_resp(x, f)[:, None], 2, 1)
    return _orig(vn)


ut.variant_H = variant_H
ut.OUT = OUT
ut.FB_VARIANTS = VARS


def report(tag="fullband13", VARS=VARS):
    R = json.load(open(os.path.join(OUT, "fullband.json")))
    os.replace(os.path.join(OUT, "fullband.json"), os.path.join(OUT, tag + ".json"))
    sets = {"CLA": ["CLA"], "Altiverb": [u for u in R["irs"] if u.startswith("Altiverb")]}
    lines = ["# 10-number (bq.json) vs 13-number (unit_chain s.json) low-band cascade, full band; null R (L) dB, variant minus plain",
             "# negative = better. Each cell: R median over the unit set (Altiverb: five S presets)", ""]
    for sn, us in sets.items():
        for vn in VARS:
            lines.append(f"## {sn}: {vn}  (absolute null R, then delta vs plain)")
            lines.append("band         " + " ".join(f"{w:>9}" for w in R["wins"]) + "  |  " + " ".join(f"{w:>9}" for w in R["wins"]))
            for b in R["bands"]:
                a, d = [], []
                for w in R["wins"]:
                    va = [R["irs"][u][vn]["null"][b][w][1] for u in us if R["irs"][u][vn]["null"][b][w]]
                    vp = [R["irs"][u]["plain"]["null"][b][w][1] for u in us if R["irs"][u]["plain"]["null"][b][w]]
                    a.append(np.median(va) if va else np.nan)
                    d.append(np.median(np.array(va) - np.array(vp)) if va else np.nan)
                lines.append(f"{b:12s} " + " ".join(f"{x:9.1f}" for x in a) + "  |  " + " ".join(f"{x:+9.1f}" for x in d))
            lines.append("")
    open(os.path.join(OUT, tag + ".log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    ut.fullband_stage(VARS)
    report()
