"""Fork G: why does NL_HF fail in some clusters? Per capture, on its own reverb windows (0.05-0.6 s), the 10-12 kHz null with
  none            : no HF correction (= the lb2 bench's model)
  used            : the correction in hf_corr.json (z = z0_LORO + s u_LORO, scalar from the 15 nearest-in-time front-tap captures)
  s_best          : the ORACLE scalar s on the all-run shape (z = z0 + s u), 1-D minimum of the 10-12 kHz null
  free9 / free3   : free per-capture complex gains on 9 / 3 hat knots (oracle upper bound of any per-capture complex scheme)
plus diagnostics (own front-tap scalar, HF energy, noise-limited null, response time of the window).
usage: g_oracle.py N_PER_CELL [workers]   -> analysis/nl_hfdrift/g_oracle_<N>.pkl"""
import os, sys, json, pickle, random
os.environ["NL_BK"] = "1"
os.environ.pop("NL_HF", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter, sosfiltfilt
from scipy.optimize import minimize_scalar
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7_lb2"
HFB = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7_lb2_hf"
import defects
KN9 = np.array([10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800.0])
KN3 = np.array([11000, 11750, 12500.0])
shape = json.load(open(OUT + "/d_hf_shape.json"))
Z0 = np.array(shape["z0_re"]) + 1j * np.array(shape["z0_im"]); U = np.array(shape["u_re"]) + 1j * np.array(shape["u_im"])
def split(z): return np.r_[z.real, z.imag]
UV = split(U)
SOS = butter(4, [10000, 12000], btype="band", fs=FS, output="sos")
NOISE = json.load(open(OUT + "/a_noise.json"))["noise"]["b10_12"]            # per-sample noise power L, R (tail of the LATE 0 files)
HFC = json.load(open("/Users/leighpierce/Git/246/analysis/nl/hf_corr.json"))
FRONT = {}; SOWN = {}; FROWS = []
TIMES = {tuple(float(v) for v in k.split("_")): t for k, t in json.load(open(OUT + "/d_times.json")).items()}
_tt = sorted(TIMES.items(), key=lambda kv: kv[1]); RUNID = {}; _rid = 0
for _i, (_k, _t) in enumerate(_tt):
    if _i and _t - _tt[_i - 1][1] > 1200: _rid += 1
    RUNID[_k] = _rid
def load_front():
    rows = pickle.load(open(OUT + "/d_front_all.pkl", "rb"))
    KI = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    for r in rows:
        c = r[1]["num"] / r[1]["den"]; k = tuple(r["key"])
        FRONT[k] = np.log(c[KI])
        if r["t"] is not None and k in RUNID:
            SOWN[k] = float(np.dot(split(FRONT[k]) - split(Z0), UV) / np.dot(UV, UV)); FROWS.append((r["t"], k))
    FROWS.sort()


def E_scalar(f, z):
    kx = np.r_[9500.0, 10000.0, KN9, 13200.0]; g = np.r_[1.0, 1.0, np.exp(z), 1.0]
    return np.interp(f, kx, g.real) + 1j * np.interp(f, kx, g.imag)


def hats(f, knots):
    """Hat functions at `knots`, neighbours (and the anchors 10 kHz / 13.2 kHz) as the foot points."""
    kx = np.r_[10000.0, knots, 13200.0]; H = []
    for i in range(1, len(kx) - 1):
        H.append(np.clip(np.minimum((f - kx[i - 1]) / (kx[i] - kx[i - 1]), (kx[i + 1] - f) / (kx[i + 1] - kx[i])), 0, 1))
    return np.array(H)


def apply_E(m, E):
    return np.fft.irfft(np.fft.rfft(m, axis=0) * E[:, None], len(m), axis=0)


def null_b(yb, mm):
    mb = sosfiltfilt(SOS, mm, axis=0)
    return np.array([10 * np.log10(np.sum((yb[:, c] - mb[:, c]) ** 2) / np.sum(yb[:, c] ** 2)) for c in (0, 1)])


def free_fit(y, m, knots):
    L = len(m); w = np.hanning(L); f = np.fft.rfftfreq(L) * FS
    sel = (f >= 10000) & (f <= 13200)
    Y = np.fft.rfft(y * w[:, None], axis=0)[sel]; Mm = np.fft.rfft(m * w[:, None], axis=0)[sel]
    H = hats(f, knots)[:, sel]
    A = np.concatenate([(Mm[:, c][None, :] * H).T for c in (0, 1)], 0)            # (2 nbins, K)
    t = np.concatenate([Y[:, c] - Mm[:, c] for c in (0, 1)])
    d = np.linalg.lstsq(A, t, rcond=None)[0]
    return d, f


def process(key):
    b = common.bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
    a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if b2 - a < int(0.05 * FS):
        return None
    n = b2 + 200
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = nlt.unwarp(rawp, b["ppm"])[a:b2]
    m = (M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"]))[a:b2]
    L = len(m); f = np.fft.rfftfreq(L) * FS
    yb = sosfiltfilt(SOS, y, axis=0)
    r = dict(key=key, trimmed=bool(tr), a=a, b=b2, o=b["o"], fr=fr, ppm=b["ppm"],
             t_resp0=(a - M.K0 - fr - b["o"]) / FS, t=TIMES.get(key))
    r["none"] = null_b(yb, m)
    ck = HFC["captures"].get("_".join(str(v) for v in key), HFC["default"])
    zu = np.array(ck["z_re"]) + 1j * np.array(ck["z_im"])
    r["used"] = null_b(yb, apply_E(m, E_scalar(f, zu)))
    r["s_used"] = float(np.dot(split(zu) - split(Z0), UV) / np.dot(UV, UV))               # projection on the all-run shape
    r["s_used_json"] = ck.get("s")
    cost = lambda s: float(np.mean(null_b(yb, apply_E(m, E_scalar(f, Z0 + s * U)))))
    grid = np.arange(-0.5, 1.01, 0.1); cs = [cost(s) for s in grid]; s0 = grid[int(np.argmin(cs))]
    res = minimize_scalar(cost, bounds=(s0 - 0.1, s0 + 0.1), method="bounded", options={"xatol": 1e-3})
    r["s_best"] = float(res.x); r["oracle"] = null_b(yb, apply_E(m, E_scalar(f, Z0 + res.x * U)))
    r["null_s0"] = null_b(yb, apply_E(m, E_scalar(f, Z0)))                                 # shape mean only (s = 0)
    r["cost_curve"] = (grid.tolist(), cs)
    for tag, kn in (("free9", KN9), ("free3", KN3)):
        d, ff = free_fit(y, m, kn)
        E = 1 + (d[:, None] * hats(f, kn)).sum(0)
        r[tag] = null_b(yb, apply_E(m, E)); r[tag + "_d"] = d
    if key in FRONT:
        zo = FRONT[key]; r["s_own"] = float(np.dot(split(zo) - split(Z0), UV) / np.dot(UV, UV))
    # estimator variants: own-run mean of the front-tap scalars (others in the same run), own front-tap scalar, closest point of the line to z = 0
    rid = RUNID.get(key); oth = [v for k2, v in SOWN.items() if RUNID.get(k2) == rid and k2 != key]
    r["run_id"] = rid; r["n_run_front"] = len(oth)
    if len(oth) >= 3:
        r["s_run"] = float(np.mean(oth)); r["run_null"] = null_b(yb, apply_E(m, E_scalar(f, Z0 + r["s_run"] * U)))
    if "s_own" in r:
        r["own_null"] = null_b(yb, apply_E(m, E_scalar(f, Z0 + r["s_own"] * U)))
    r["s_zero"] = float(-np.dot(split(Z0), UV) / np.dot(UV, UV))
    r["zline_null"] = null_b(yb, apply_E(m, E_scalar(f, Z0 + r["s_zero"] * U)))
    r["zline_dist"] = float(np.linalg.norm(split(Z0) + r["s_zero"] * UV))
    if r["t"] is not None and FROWS:
        near = sorted((x for x in FROWS if x[1] != key), key=lambda x: abs(x[0] - r["t"]))[:15]
        r["neigh_out_frac"] = float(np.mean([RUNID.get(k2) != rid for _, k2 in near])); r["neigh_span_s"] = float(max(abs(x[0] - r["t"]) for x in near))
    # energy / noise
    Eb = float(np.sum(yb ** 2, 0).mean())
    r["hf_level_db"] = 10 * np.log10(Eb / L / 1.0 + 1e-30)
    r["noise_null"] = [10 * np.log10(NOISE[c] * L / np.sum(yb[:, c] ** 2)) for c in (0, 1)]
    r["frac_hf_energy_db"] = 10 * np.log10(np.sum(yb ** 2) / np.sum(y ** 2))
    return r


def cluster(key, trimmed):
    if key[0] == 0.2 and key[1] == -0.2: return "fam(0.2,-0.2,*)"
    if trimmed: return "trimmed"
    if key[0] == 0.3: return "SPAN0.3"
    if key[0] == 0.1: return "SPAN0.1"
    return "rest"


def select(npc):
    random.seed(21)
    keys = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS)
    cells = {}
    for k in keys:
        n = "_".join(str(v) for v in k) + ".json"
        try:
            A = json.load(open(f"{HFB}/{n}")); B = json.load(open(f"{common.BENCH}/{n}"))
        except FileNotFoundError:
            continue
        if "10000-12000" not in A["bands"]: continue
        worse = bool(max(np.array(A["bands"]["10000-12000"]) - np.array(B["bands"]["10000-12000"])) > 0.5)
        cells.setdefault((cluster(k, A["trimmed"]), worse), []).append(k)
    sel = []
    for c, v in sorted(cells.items()):
        sel += random.sample(v, min(npc, len(v)))
    return sel, {c: len(v) for c, v in cells.items()}


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    npc = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    sel, sizes = select(npc); print(len(sel), "captures; cell sizes", sizes, flush=True)
    load_front()
    out = []   # (workers load their own copy)
    if nw == 1:
        for i, k in enumerate(sel):
            r = process(k)
            if r: out.append(r)
            if i % 20 == 0: print(i, flush=True)
    else:
        with ProcessPoolExecutor(nw, initializer=load_front) as ex:
            for i, r in enumerate(ex.map(process, sel, chunksize=2)):
                if r: out.append(r)
                if i % 20 == 0: print(i, flush=True)
    pickle.dump(out, open(OUT + f"/g_oracle_{npc}.pkl", "wb"))
    print("done", len(out))
