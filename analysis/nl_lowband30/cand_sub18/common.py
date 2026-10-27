"""Helpers for the sub-18 Hz / estimation-design study (cand_sub18). Read-only on tools/ and the existing analysis folders."""
import os, sys
os.environ.setdefault("NL_BK", "1"); os.environ.setdefault("NL_BK_ALPHA", "0.7")
os.environ.setdefault("LB_BENCH", "bench_warp2_bk_kf_a0.7")
os.environ["NL_LOWBAND"] = os.environ.get("NL_LOWBAND", "0")
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30")
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
from lb_common import *     # noqa
import lowband_joint as LJ  # noqa
HERE = os.path.dirname(os.path.abspath(__file__))
FEED = {0.0: None, 0.1: None, 0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}

def sess(r):
    return int(round(r["front_shift"]))

def catalog(span_lo=0.3, span_hi=0.5, late_max=0.4, early_min=0.3):
    B = bench_all()
    out = []
    for k, r in sorted(B.items()):
        if r["trimmed"] or not (span_lo <= k[0] <= span_hi) or k[2] > late_max or k[1] < early_min: continue
        out.append(k)
    return B, out
