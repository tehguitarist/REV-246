"""LATE 0.2: the late loops' effective RT (NonLin_program 10.4). The adopted bench (bench_warp2_bk_kf_a0.7_lb2_hf_g5) with
its own per-capture o, clock and front shift fixed; only the loop gains change (gains per channel by least squares as the
bench). Stages:
  scan  [keys]   per capture: the bench's window nulls with the loops at RT_eff = LATE (1 + e), e on a grid, and per-loop
                 (3 free e) -> late02/scan/<key>.json
  report         per-step medians, held-out split
"""
import glob, json, os, sys
os.environ.setdefault("NL_BK", "1")
os.environ.setdefault("EMT_KERNEL", "chain3")
for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
import numpy as np
from scipy.signal import sosfiltfilt
from scipy.optimize import minimize
import defects
import nl_structure as nl
import nl_early as ne
import nl_late as nlt
import nl_model as nm

FS, K0 = nl.FSD, nl.K0
OUT = os.path.join(nl.OUT, "late02", "scan")
EPS = np.round(np.arange(-0.02, 0.0401, 0.0025), 4)
_ORIG = nlt.loop_gain
_E = {}


def _gain(P, late):
    return 10 ** (-3 * P / (FS * late * (1 + _E.get(P, 0.0))))


def set_e(e):
    """e: scalar (all loops) or {P: e}."""
    _E.clear()
    _E.update({P: e for P in nlt.LOOPS} if np.isscalar(e) else e)
    nlt.loop_gain = _gain
    nm.spectra.cache_clear()


def kstr(k):
    return "_".join(str(v) for v in k)


def keys_for(late, every=1):
    ks = []
    for k in sorted(nl.table()):
        if k[2] != late or k in defects.NL_DEFECTS:
            continue
        b = json.load(open(os.path.join(nm.BENCH, kstr(k) + ".json")))
        n = int(round((b["o"] - 1) * nm.OSTEP))
        if b["trimmed"] or not (n == 0 or k in defects.NL_SLIPPED):
            continue
        ks.append(k)
    return ks[::every]


class Cap:
    def __init__(self, key):
        self.key = key
        b = json.load(open(os.path.join(nm.BENCH, kstr(key) + ".json")))
        self.b = b
        raw = nl.ir(nl.table()[key])
        h = sosfiltfilt(ne.BAND, raw, axis=0)
        fr = b["front_shift"]
        self.n = n = min(nm.usable_end(h), K0 + int(fr) + int(nm.WINDOWS[-1][1] * FS))
        rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
        y = nlt.unwarp_from(nlt.upsample(rawp), n, b["ppm"])[:n]
        self.y = sosfiltfilt(ne.BAND, y, axis=0)
        i0, _ = nm.window_start(h, fr)
        noise = json.load(open(nm.NOISE))
        self.wins = {}
        for t0, t1 in nm.WINDOWS:
            a, bb = max(i0, K0 + int(fr) + int(t0 * FS)), min(n, K0 + int(fr) + int(t1 * FS))
            if bb - a >= int(0.02 * FS) and np.sum(self.y[a:bb, 0] ** 2) >= 10 * noise["ch0"] * (bb - a):
                self.wins[f"{t0}-{t1}"] = (a, bb, [float(10 * np.log10(noise[f"ch{c}"] * (bb - a) / np.sum(self.y[a:bb, c] ** 2))) for c in (0, 1)])

    def nulls(self, e):
        set_e(e)
        m = nm.render(self.key, self.n, self.b["o"], self.b["front_shift"])
        g = [(self.y[:, c] @ m[:, c]) / (m[:, c] @ m[:, c]) for c in (0, 1)]
        m = m * np.array(g)
        return {w: [float(10 * np.log10(np.sum((self.y[a:b, c] - m[a:b, c]) ** 2) / np.sum(self.y[a:b, c] ** 2))) for c in (0, 1)]
                for w, (a, b, _) in self.wins.items()}

    def cost(self, e):                                  # the bench's cost: every window with signal counts equally, both channels
        v = self.nulls(e)
        return float(sum(x for w in v.values() for x in w))


def scan(key):
    path = os.path.join(OUT, kstr(key) + ".json")
    if os.path.exists(path):
        return
    c = Cap(key)
    grid = {float(e): c.nulls(float(e)) for e in EPS}
    costs = np.array([sum(x for w in grid[float(e)].values() for x in w) for e in EPS])
    j = int(np.argmin(costs))
    e1 = float(EPS[j])
    if 0 < j < len(EPS) - 1:                            # parabola through the 3 best
        y0, y1, y2 = costs[j - 1:j + 2]
        d = (y0 - y2) / (2 * (y0 - 2 * y1 + y2)) if (y0 - 2 * y1 + y2) > 0 else 0.0
        e1 = float(EPS[j] + d * (EPS[1] - EPS[0]))
    r = minimize(lambda v: c.cost(dict(zip(nlt.LOOPS, v))), [e1] * 3, method="Nelder-Mead",
                 options=dict(xatol=2e-4, fatol=1e-3, maxiter=120, initial_simplex=[[e1] * 3, [e1 + 0.005, e1, e1], [e1, e1 + 0.005, e1], [e1, e1, e1 + 0.005]]))
    out = dict(key=list(key), floor={w: v[2] for w, v in c.wins.items()}, grid={str(k): v for k, v in grid.items()},
               e_best=e1, null_best=c.nulls(e1), e_loops=[float(x) for x in r.x], null_loops=c.nulls(dict(zip(nlt.LOOPS, r.x))),
               null_nominal=grid[0.0])
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(path, "w"))
    print(kstr(key), f"e_best {e1:+.4f} loops", " ".join(f"{x:+.4f}" for x in r.x), "| R nominal",
          " ".join(f"{v[1]:.1f}" for v in grid[0.0].values()), "| best", " ".join(f"{v[1]:.1f}" for v in out["null_best"].values()), flush=True)


if __name__ == "__main__":
    st = sys.argv[1]
    if st == "scan":
        late = float(sys.argv[2])
        every = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        lim = int(sys.argv[4]) if len(sys.argv) > 4 else None
        ks = keys_for(late, every)[:lim]
        print(len(ks), "captures", flush=True)
        for k in ks:
            scan(k)
