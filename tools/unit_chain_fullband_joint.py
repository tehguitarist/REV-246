"""§37.x: the joint (CLA + Altiverb) cascades of unit_chain_joint.py through the same fullband scoring as
unit_chain_fullband.py (imported; plain, CLA's and the Altiverb's own 10- and 13-number cascades as comparison).
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 3 -- \
        .venv/bin/python tools/unit_chain_fullband_joint.py        # -> analysis/unit_chain/fullband_joint.{json,log}
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import unit_chain_fullband as fb  # noqa: E402
import unit_transfer as ut  # noqa: E402
import unit_chain as uc  # noqa: E402

J = json.load(open(os.path.join(fb.OUT, "joint.json")))
VARS = ["plain", "bq_CLA", "bq_Altiverb RT 1.6", "b13_CLA@150-300", "b13_Alt@150-300", "b13_CLA@0-50", "b13_Alt@0-50"] + ["j_" + k for k in J]
_prev = ut.variant_H


def variant_H(vn):
    if vn.startswith("j_"):
        x = J[vn[2:]]["x"]
        return lambda f: np.repeat(uc.lb_resp(x, f)[:, None], 2, 1)
    return _prev(vn)


ut.variant_H = variant_H
SUB = os.path.join(fb.OUT, "joint_run")           # own folder: runs beside unit_chain_fullband.py
os.makedirs(SUB, exist_ok=True)
ut.OUT = fb.OUT = SUB
fb.VARS = VARS
if __name__ == "__main__":
    ut.fullband_stage(VARS)
    fb.report("fullband_joint", VARS)
