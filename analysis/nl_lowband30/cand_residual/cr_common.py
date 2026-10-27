"""Shared helpers for the anatomy of the 100 Hz-1 kHz (and 20-100 Hz) residual left after the shipped low-band transfer (NL_LOWBAND=2).
Read-only on tools/, docs/, analysis/ outside this directory. Model: NL_BK=1 (alpha 0.7 branch kernel), bench alignment bench_warp2_bk_kf_a0.7_lb2."""
import os, sys, json, time
os.environ["NL_BK"] = "1"
os.environ["LB_BENCH"] = "bench_warp2_bk_kf_a0.7_lb2"
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30")
from lb_common import *        # noqa  (np, M, ne, nlt, nl, FS, K0, bands, load, bench_all, ...)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
os.makedirs(CACHE, exist_ok=True)
NCACHE = 12600
T0 = 0.05
T1_MAX = 0.45
EARLIES = (0.7, 1.0, 1.5, 2.0, 3.0, 10.0)
# last audible tap per SPAN, seconds (early section only: after it only the AP tails and noise are left)
TLAST = {s: ne.taps(0, n)[-1][0] / FS for s, n in ne.N_STAGES.items() if n}


def tag(k):
    return "_".join(str(v) for v in k)


def window(k, r):
    """(a, b) sample window in cache indices: 0.05 s ... min(0.45 s, last tap + 0.02 s) from the front frame zero."""
    t1 = min(T1_MAX, TLAST.get(k[0], 9) + 0.02)
    z = K0 + int(r["front_shift"])
    return z + int(T0 * FS), z + int(t1 * FS)


def cached(k, r):
    f = os.path.join(CACHE, tag(k) + ".npz")
    if os.path.exists(f):
        d = np.load(f)
        return d["y"], d["m"]
    y, m, _ = load(k, r, NCACHE, band=False)
    np.savez(f, y=y.astype(np.float64), m=m.astype(np.float64))
    return y, m


def band_sos(lo, hi):
    return bands(lo, hi)


def fpow(x):
    return float(np.sum(x * x))
