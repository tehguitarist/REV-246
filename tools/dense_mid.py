"""Is the §25.7 baseline's remaining mid-band excess (250 Hz-1 kHz, +8...+13 dB over repeat at 50-300 ms in every
group, incoherent with the model: model_s_results §25.9) a setting-independent error of the loop excitation? A diagnostic,
not a model: one dense correction per loop and channel added to the excitations (q1.Joint's X, lags S0 .. S0 + P_i + EXT),
fitted by LSQR on the 80-1500 Hz residual of the training settings' captures, and scored on settings it never saw.

    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/dense_mid.py --stage prep [--only keys]
    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/dense_mid.py --stage fit --fold A|B [--iters 200]
    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/dense_mid.py --stage score --fold A|B
    EMT_KERNEL=chain3 .venv/bin/python tools/dense_mid.py --stage report

prep    per setting (lowband_shelf.bb_prep: the §25.7 model, clocks, fade): the 80-1500 Hz residual of each capture, its
        weight (fade ramp, 50 ms to the compensated end), the loop spectra R_i, the baseline model -> analysis/dense_mid/<kernel>/prep/
fold    A trains on RT 0.3 / 0.8 / 1.4 / 2.2 and scores RT 0.6 / 1.0 / 1.8 / 2.5; B the other way round.
score   per test setting: R nulls per octave and window, baseline and baseline + correction, against the next AMPL capture.
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.signal import butter, sosfiltfilt, sosfreqz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lowband_shelf as ls  # noqa: E402

q2, q1 = ls.q2, ls.q1
FS, PRE = ls.FS, ls.PRE
BAND = (80, 1500)
SOSM = butter(3, BAND, btype="band", fs=FS, output="sos")
FOLDS = {"A": (0.3, 0.8, 1.4, 2.2), "B": (0.6, 1.0, 1.8, 2.5)}
T_FROM = 0.05
OUT = os.path.join(ROOT, "analysis", "dense_mid", q1.ms.KTAG or "library")


def prep(key):
    path = os.path.join(OUT, "prep", key + ".npz")
    if os.path.exists(path):
        return
    d = ls.bb_prep(key)
    r = ls.reverb(d, d["Gs0"])
    N = d["N"]
    w = d["w"]
    Bm = np.abs(sosfreqz(SOSM, worN=w)[1]) ** 2                     # filtfilt = |H|^2, zero phase
    R = np.array([d["G"]["PB"] / (1 - Gi * np.exp(-1j * w * Pi)) for Gi, Pi in zip(d["Gs0"], q1.P)])
    arrs = dict(R=R.astype(np.complex64), Bm=Bm.astype(np.float32), N=N, lab=np.array([d["lab"]["lf"], d["lab"]["rt"], d["lab"]["hf"]]))
    for j, c in enumerate(d["caps"]):
        n = len(c["y"])
        m = ls.model(d, c, r)[:n] * c["g"][:n]
        wt = np.zeros(n)
        i0, i1 = PRE + int(T_FROM * FS), min(PRE + int(c["te"] * FS), n)
        wt[i0:i1] = c["g"][i0:i1, 1]
        arrs[f"res{j}"] = sosfiltfilt(SOSM, c["y"] - m, axis=0).astype(np.float32)
        arrs[f"w{j}"] = wt.astype(np.float32)
        arrs[f"y{j}"] = c["y"].astype(np.float32)
        arrs[f"m{j}"] = m.astype(np.float32)
        arrs[f"g{j}"] = c["g"].astype(np.float32)
        arrs[f"amp{j}"] = np.array(c["amp"])
        arrs[f"te{j}"] = np.array(c["te"])
    np.savez(path, **arrs)


THROUGH_G = False                                                   # --through-g: a second correction per loop through one pass of G_i
_LAW = {}


def loop_G(key, N):
    """The §25.7 baseline's loop responses G_i for one setting on an N grid (no captures needed)."""
    if "law" not in _LAW:
        _LAW["law"] = q2.final_injection(compact=ls.lb.lf_.COMPACT)[1]
        _LAW["ap"] = json.load(open(os.path.join(q2.OUT, ls.lb.lf_.AP230)))
    info = q2.q1_tables()[0]
    t, lfdc, kw = ls.base_kw(key, ls.HF_CELLS_JSON)
    law = json.loads(json.dumps(_LAW["law"]))
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    for lf_k, x in lfdc.items():
        law["lf"][float(lf_k)]["dc"] = float(x)
    lab = dict(lf=info[key]["lf"], rt=info[key]["rt"], hf=info[key]["hf"])
    b230 = _LAW["ap"]["steps"][_LAW["ap"]["settings"][key]["step"]]
    w = 2 * np.pi * np.fft.rfftfreq(N)
    return q2.final_Gs(lab, law, w, lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=b230, **kw)


def load(key, full=False):
    z = np.load(os.path.join(OUT, "prep", key + ".npz"))
    k = sorted(int(n[3:]) for n in z.files if n.startswith("res"))
    files = []
    for j in k:
        w = z[f"w{j}"].astype(np.float64)
        e = int(np.nonzero(w)[0][-1]) + 1                               # no samples past the weighted window
        res = z[f"res{j}"] if FIT_BAND == BAND else sosfiltfilt(SOS_FIT, (z[f"y{j}"] - z[f"m{j}"]).astype(np.float64), axis=0)
        files.append(dict(y=res[:e].astype(np.float64), w=w[:e], amp=0))
    R = list(z["R"].astype(np.complex128))
    if THROUGH_G:
        R = R + [r * Gi for r, Gi in zip(R, loop_G(key, int(z["N"])))]
    Bm = z["Bm"] if FIT_BAND == BAND else np.abs(sosfreqz(SOS_FIT, worN=2 * np.pi * np.fft.rfftfreq(int(z["N"])))[1]) ** 2
    g = dict(N=int(z["N"]), R=[r * Bm for r in R], PB=np.zeros(int(z["N"]) // 2 + 1), files=files)
    if full:
        g["R_full"] = R
        g["caps"] = [dict(y=z[f"y{j}"].astype(np.float64), m=z[f"m{j}"].astype(np.float64), g=z[f"g{j}"].astype(np.float64),
                          amp=int(z[f"amp{j}"]), te=float(z[f"te{j}"])) for j in k]
    g["lab"] = z["lab"].tolist()
    g["key"] = key
    return g


DIFF_COLS = ("g20", "g72", "g230", "M20", "M72", "M230", "level")
DIFF2_COLS = ("ff20", "ff72", "ff230", "fb20", "fb72", "fb230", "level")   # --basis diffuser2: feedforward and feedback gains apart


def diff_D(key, N):
    """--basis diffuser: the relative change of the diffuser's response per unit change of each all-pass gain (d log AP / dg),
    and for each stage's delay one sample longer (AP(M + 1) / AP(M) - 1); AP230 with this setting's two-tap damping; and
    the loop paths' level re the cluster (a reference column: the same as --basis gain with L = R)."""
    w = 2 * np.pi * np.fft.rfftfreq(N)
    ap = json.load(open(os.path.join(q2.OUT, ls.lb.lf_.AP230)))
    b = ap["steps"][ap["settings"][key]["step"]]
    g0 = q2.FINAL_DIFF["g20"]
    F = lambda M: ((1 - b) + b * np.exp(-1j * w)) if M >= 200 else 1.0
    A = lambda M, g: (-g + np.exp(-1j * w * M) * F(M)) / (1 - g * np.exp(-1j * w * M) * F(M))
    h = 1e-5
    if BASIS == "diffuser2":                                        # A = (-g1 + z^-M F) / (1 - g2 z^-M F): d log A / d g1 and d g2
        A2 = lambda M, g1, g2: (-g1 + np.exp(-1j * w * M) * F(M)) / (1 - g2 * np.exp(-1j * w * M) * F(M))
        return [(A2(M, g0 + h, g0) - A2(M, g0 - h, g0)) / (2 * h) / A(M, g0) for M in (20, 72, 230)] + \
            [(A2(M, g0, g0 + h) - A2(M, g0, g0 - h)) / (2 * h) / A(M, g0) for M in (20, 72, 230)] + [np.ones_like(w)]
    return [(A(M, g0 + h) - A(M, g0 - h)) / (2 * h) / A(M, g0) for M in (20, 72, 230)] + \
        [A(M + 1, g0) / A(M, g0) - 1 for M in (20, 72, 230)] + [np.ones_like(w)]


def x_base():
    """The baseline's own excitation (the compact model's 5 signed taps per loop and channel), Joint layout, no cluster / F."""
    v = q2.final_injection(compact=ls.lb.lf_.COMPACT)[0]
    set_win()
    return v[8 + 2 * q1.LFIR:]


def fit_diffuser(fold):
    """Least squares for the 6 diffuser numbers on the training settings' fit-band residual."""
    x0 = x_base()
    groups = [load(k) for k in keys_for(FOLDS[fold])]
    off = 8 + 2 * q1.LFIR
    tgt = q1.Joint(groups).target()
    cols = []
    cols_ = DIFF2_COLS if BASIS == "diffuser2" else DIFF_COLS
    for k in range(len(cols_)):
        Gk = [dict(g, R=[r * D for r in g["R"]]) for g in groups for D in [diff_D(g["key"], g["N"])[k]]]
        cols.append(q1.Joint(Gk).forward(np.r_[np.zeros(off), x0]))
        del Gk
    A = np.array(cols)
    th = np.linalg.solve(A @ A.T, A @ tgt)
    res = tgt - th @ A
    print(f"fold {fold}: {len(groups)} settings, {len(tgt)} samples; diffuser " + ", ".join(f"{n} {v:+.6f}" for n, v in zip(cols_, th)) +
          f"; training residual {10 * np.log10(np.sum(res ** 2) / np.sum(tgt ** 2)):+.2f} dB re the baseline's", flush=True)
    np.save(os.path.join(OUT, f"theta_{tag(fold)}.npy"), th)


SUBSET = {}                                                         # --subset lf=1.0,hf=1.0: restrict train and test


def keys_for(rts, subset=True):
    """The prepared settings at those RT steps (within SUBSET unless subset=False)."""
    info = q2.q1_tables()[0]
    ap = json.load(open(os.path.join(q2.OUT, ls.lb.lf_.AP230)))
    b0 = lambda k: ap["steps"][ap["settings"][k]["step"]] < 0.01                     # AP230 without damping
    return sorted(k for k in info if info[k]["rt"] in rts and os.path.exists(os.path.join(OUT, "prep", k + ".npz"))
                  and (not subset or all(info[k][c] == v for c, v in SUBSET.items() if c not in ("b0", "rtmin"))
                  and (not SUBSET.get("b0") or b0(k)) and info[k]["rt"] >= SUBSET.get("rtmin", 0)))


def tag(fold):
    return fold + ("_" + "_".join(f"{c}{v:g}" for c, v in sorted(SUBSET.items())) if SUBSET else "") + ("_G" if THROUGH_G else "") + \
        ("" if BASIS == "dense" else "_" + BASIS) + ("" if FIT_BAND == BAND else f"_band{FIT_BAND[0]}-{FIT_BAND[1]}") + \
        ("" if WIN is None else f"_win{WIN[0]}-{WIN[1]}") + ("" if BASE is None else f"_on-{BASE}")


def stag(fold):
    """Score tag: the fit's tag, plus the projection band if any."""
    return tag(fold) + ("" if PROJECT is None else f"_proj{PROJECT[0]}-{PROJECT[1]}") + ("_all" if TEST_ALL else "")


PROJECT = None                                                      # --project lo,hi: score the correction band-limited to lo-hi Hz
TEST_ALL = False                                                    # --test-all: score every setting at the test RT steps, outside SUBSET too


BASIS = "dense"                                                     # --basis dense | taps | tapsall | gain | gainx | kernelK | ckernelK | lkernelK | common | diffuser | diffuser2
FIT_BAND, SOS_FIT = BAND, SOSM                                      # --band lo,hi: the fitting band (prep's residual is 80-1500 Hz)


def basis_map(J):
    """Sparse map from the basis coefficients to the dense excitation vector (None = dense: identity)."""
    if BASIS == "dense":
        return None
    from scipy.sparse import lil_matrix
    off, o = {}, 0
    for c in (0, 1):
        for i, L in enumerate(J.Lx):
            off[(c, i)] = o
            o += L
    if BASIS == "common":                                           # one correction per channel shared by the four loops
        n = min(J.Lx)
        M = lil_matrix((o, 2 * n))
        for c in (0, 1):
            for i in range(4):
                for t in range(n):
                    M[off[(c, i)] + t, c * n + t] = 1.0
        return M.tocsr()
    taps = q2.tap_table()
    if BASIS in ("tapsall", "gain", "gainx"):                       # each tap's correction in its own loop and/or the other three
        ncol = dict(tapsall=len(taps), gain=2, gainx=4)[BASIS]
        M = lil_matrix((o, ncol))
        for n, (i, c, t, sg, j) in enumerate(taps):
            for k in range(4):
                if BASIS == "tapsall":
                    M[off[(c, k)] + t, n] = sg
                elif BASIS == "gainx":
                    M[off[(c, k)] + t, 2 * c + (k != i)] = sg
                elif k == i:
                    M[off[(c, k)] + t, c] = sg
        return M.tocsr()
    K = int(BASIS.split("kernel")[1]) if "kernel" in BASIS else 0     # kernelK: one kernel at lags -K..K for every tap; ckernelK: one per channel;
    nc = 2 if BASIS.startswith("ckernel") else 8 if BASIS.startswith("lkernel") else 1   # lkernelK: one per loop and channel
    M = lil_matrix((o, len(taps) if BASIS == "taps" else nc * (2 * K + 1)))
    for n, (i, c, t, sg, j) in enumerate(taps):
        if BASIS == "taps":
            M[off[(c, i)] + t, n] = sg
        else:
            col = (c if nc == 2 else 4 * c + i if nc == 8 else 0) * (2 * K + 1)
            for k in range(max(-K, -t), min(K, J.Lx[i] - 1 - t) + 1):   # within this loop's excitation
                M[off[(c, i)] + t + k, col + k + K] += sg
    return M.tocsr()


WIN = None                                                          # --win s0,ext: the correction's support (default q1's 230, 1700)
BASE = None                                                         # --base TAG: a fitted correction held fixed underneath (x_<fold>_<subset>_<TAG>.npy)


def base_x(fold):
    """The --base correction of this fold (dense layout), or None."""
    if BASE is None:
        return None
    sub = "_".join(f"{c}{v:g}" for c, v in sorted(SUBSET.items()))
    return np.load(os.path.join(OUT, f"x_{fold}{'_' + sub if sub else ''}_{BASE}.npy"))


def set_win():
    q1.EXT = 1700 if WIN is None else WIN[1]
    if WIN is not None:
        q1.S0 = WIN[0]


def fit(fold, iters):
    from scipy.sparse.linalg import LinearOperator, lsqr
    set_win()
    if BASIS.startswith("diffuser"):
        return fit_diffuser(fold)
    groups = [load(k) for k in keys_for(FOLDS[fold])]
    J = q1.Joint(groups)
    off = 8 + 2 * q1.LFIR
    M = basis_map(J)
    nx = (J.size() - off) if M is None else M.shape[1]
    fwd = lambda x: J.forward(np.r_[np.zeros(off), x if M is None else M @ x])
    adj = lambda z: J.adjoint(z)[off:] if M is None else M.T @ J.adjoint(z)[off:]
    tgt = J.target()
    xb = base_x(fold)
    if xb is not None:                                              # fit what the base correction leaves
        tgt = tgt - J.forward(np.r_[np.zeros(off), xb])
    op = LinearOperator((len(tgt), nx), matvec=fwd, rmatvec=adj, dtype=float)
    x = lsqr(op, tgt, iter_lim=iters, atol=1e-12, btol=1e-12)[0]
    res = tgt - fwd(x)
    print(f"fold {fold}: {len(groups)} settings, {len(tgt)} samples, {nx} unknowns; training residual "
          f"{10 * np.log10(np.sum(res ** 2) / np.sum(tgt ** 2)):+.2f} dB re the baseline's ({FIT_BAND[0]}-{FIT_BAND[1]} Hz, weighted)", flush=True)
    if M is not None:
        np.save(os.path.join(OUT, f"theta_{tag(fold)}.npy"), x)
        x = M @ x
    elif xb is not None:
        np.save(os.path.join(OUT, f"theta_{tag(fold)}.npy"), x)    # the increment over the base
    np.save(os.path.join(OUT, f"x_{tag(fold)}.npy"), x if xb is None else x + xb)   # x: the total


def score(fold):
    set_win()
    if BASIS.startswith("diffuser"):
        th, x = np.load(os.path.join(OUT, f"theta_{tag(fold)}.npy")), x_base()
    else:
        x = np.load(os.path.join(OUT, f"x_{tag(fold)}.npy"))
    test = [r for r in sorted({rt for f_ in FOLDS.values() for rt in f_}) if r not in FOLDS[fold]]
    os.makedirs(os.path.join(OUT, f"score_{stag(fold)}"), exist_ok=True)
    for k in keys_for(test, subset=not TEST_ALL):
        path = os.path.join(OUT, f"score_{stag(fold)}", k + ".json")
        if os.path.exists(path):
            continue
        g = load(k, full=True)
        J = q1.Joint([dict(g, R=g["R_full"])])
        xb = base_x(fold) if not BASIS.startswith("diffuser") else None
        _, _, xs = J.split(np.r_[np.zeros(8 + 2 * q1.LFIR), x if xb is None else x - xb])
        if PROJECT is not None:                                     # keep only the correction's PROJECT band (zero-phase, zero-padded; not the base)
            sos = butter(4, PROJECT, btype="band", fs=FS, output="sos")
            xs = [[sosfiltfilt(sos, np.r_[np.zeros(4000), v, np.zeros(4000)], padtype=None)[4000:4000 + len(v)] for v in row] for row in xs]
        if xb is not None:
            xs = [[u + v for u, v in zip(r1, r2)] for r1, r2 in zip(xs, J.split(np.r_[np.zeros(8 + 2 * q1.LFIR), xb])[2])]
        F, X = J._spectra(np.zeros((2, q1.LFIR)), xs, g["N"])
        R = np.array(g["R_full"])
        if BASIS.startswith("diffuser"):                                     # the change of the baseline's own loop paths
            R = R * sum(t * D for t, D in zip(th, diff_D(k, g["N"])))[None, :]
        add = np.fft.irfft(np.einsum("cik,ik->ck", X, R), g["N"], axis=1).T
        rows = []
        caps = g["caps"]
        for j, c in enumerate(caps):
            cb = caps[(j + 1) % len(caps)]
            n = min(len(c["y"]), len(cb["y"]))
            m1 = c["m"][:n] + add[:n] * c["g"][:n]
            te = min(c["te"], cb["te"])
            res = {}
            for f, sos in ls.SOS_OCT.items():
                y, a0, a1, yb = (sosfiltfilt(sos, v[:n], axis=0)[:, 1] for v in (c["y"], c["m"], m1, cb["y"]))
                for a, b in ls.W_OCT:
                    i0, i1 = PRE + int(a / 1000 * FS), min(PRE + int(b / 1000 * FS), PRE + int(te * FS), n)
                    if i1 - i0 >= int(0.05 * FS):
                        e = lambda u, v: float(10 * np.log10(np.sum((u - v) ** 2) / np.sum(u ** 2)))
                        res[f"{f}_{a}-{b}"] = [e(y[i0:i1], a0[i0:i1]), e(y[i0:i1], a1[i0:i1]), e(y[i0:i1], yb[i0:i1])]
            bb = ls.nulls_broad(c["y"][:n], c["m"][:n], te), ls.nulls_broad(c["y"][:n], m1, te), ls.nulls_broad(c["y"][:n], cb["y"][:n], te)
            rows.append(dict(oct=res, broad=bb))
        json.dump(dict(key=k, lab=g["lab"], rows=rows), open(path, "w"))
        print(k, flush=True)


def report():
    S = [json.load(open(p)) for f_ in FOLDS for p in glob.glob(os.path.join(OUT, f"score_{stag(f_)}", "*.json"))]
    out = [f"## Dense mid-band excitation correction, held out by RT step ({len(S)} test settings): R, median over captures; "
           "octave nulls baseline -> + correction [excess over repeat]"]
    for a, b in ls.W_OCT:
        w = f"{a}-{b}"
        out.append(f"\n### {w} ms")
        for f in ls.OCT:
            v = np.array([r["oct"][f"{f}_{w}"] for s in S for r in s["rows"] if f"{f}_{w}" in r["oct"]])
            if len(v):
                out.append(f"  {f:5d} Hz: {np.median(v[:, 0]):6.1f} -> {np.median(v[:, 1]):6.1f} [{np.median(v[:, 0] - v[:, 2]):+5.1f} -> {np.median(v[:, 1] - v[:, 2]):+5.1f}]")
        v = np.array([[r["broad"][0][w][1], r["broad"][1][w][1]] for s in S for r in s["rows"] if w in r["broad"][0]])
        if len(v):
            out.append(f"  100 Hz-10 kHz: {np.median(v[:, 0]):.1f} -> {np.median(v[:, 1]):.1f}")
    if TEST_ALL:                                                    # 100 Hz-10 kHz by group: inside SUBSET's kind of setting or not
        ap = json.load(open(os.path.join(q2.OUT, ls.lb.lf_.AP230)))
        inside = lambda s: ap["steps"][ap["settings"][s["key"]]["step"]] < 0.01 and s["lab"][1] >= SUBSET.get("rtmin", 0)
        out.append("\n### 100 Hz-10 kHz by group, R, median over captures: baseline -> + correction [excess over repeat]")
        for nm, sel in (("all", lambda s: True), ("undamped, RT >= rtmin", inside), ("damped or RT < rtmin", lambda s: not inside(s))):
            row = []
            for a, b in ls.W_OCT:
                w = f"{a}-{b}"
                v = np.array([[r["broad"][k][w][1] for k in range(len(r["broad"]))] for s in S if sel(s) for r in s["rows"] if w in r["broad"][0]])
                if len(v):
                    ex = f" [{np.median(v[:, 0] - v[:, 2]):+.1f} -> {np.median(v[:, 1] - v[:, 2]):+.1f}]" if v.shape[1] > 2 else ""
                    row.append(f"{w}: {np.median(v[:, 0]):.1f} -> {np.median(v[:, 1]):.1f}{ex}")
            out.append(f"  {nm} ({sum(map(sel, S))}): " + "; ".join(row))
    open(os.path.join(OUT, f"report{stag('')}.log"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def decode(lo=200, hi=1500, H=30):
    """The common kernel (--basis common) read off in the lo-hi band: consistency between folds and fit bands, L/R lag,
    lag profile, and its shape around every read tap (projection on the band-limited impulse and its derivative)."""
    sos = butter(4, (lo, hi), btype="band", fs=FS, output="sos")
    bp = lambda v: sosfiltfilt(sos, np.r_[np.zeros(4000), v, np.zeros(4000)], padtype=None)[4000:4000 + len(v)]
    cc = lambda u, v: float(u @ v / np.sqrt(u @ u * v @ v))
    th = {f: np.load(os.path.join(OUT, f"theta_{tag(f)}.npy")) for f in FOLDS}
    n = len(th["A"]) // 2
    k = {f: [bp(th[f][c * n:(c + 1) * n]) for c in (0, 1)] for f in FOLDS}
    out = [f"## Common kernel {tag('')}, read in {lo}-{hi} Hz (lags from {q1.S0})"]
    for c, nm in ((0, "L"), (1, "R")):
        raw = [th[f][c * n:(c + 1) * n] for f in FOLDS]
        out.append(f"{nm}: fold A vs B {cc(k['A'][c], k['B'][c]):.3f}; "
                   f"in-band / total energy {np.mean([(k[f][c] ** 2).sum() / (r ** 2).sum() for f, r in zip(FOLDS, raw)]):.2f}")
        if BASIS == "common" and FIT_BAND != BAND:
            t0 = tag("A").replace(f"_band{FIT_BAND[0]}-{FIT_BAND[1]}", "")
            p0 = os.path.join(OUT, f"theta_{t0}.npy")
            if os.path.exists(p0):
                out.append(f"   vs the {BAND[0]}-{BAND[1]} Hz fit (fold A): {cc(k['A'][c], bp(np.load(p0)[c * n:(c + 1) * n])):.3f}")
    for f in FOLDS:
        a, b = k[f]
        r = np.array([a[max(0, -s):n - max(0, s)] @ b[max(0, s):n - max(0, -s)] for s in range(-500, 501)]) / np.sqrt(a @ a * b @ b)
        i = np.argsort(-np.abs(r))[:3]
        out.append(f"fold {f}: R vs L cross-correlation, best lags " + ", ".join(f"{j - 500:+d} ({r[j]:+.2f})" for j in i) +
                   f"; at +332 {r[832]:+.2f}")
    imp = np.zeros(n); imp[n // 2] = 1
    bi = bp(imp)[n // 2 - H:n // 2 + H + 1]
    di = np.gradient(bi)
    bi, di = bi / np.linalg.norm(bi), di / np.linalg.norm(di)
    out.append(f"\nAround each read tap (+-{H} lags, sign removed, mean of folds): gain / delay projections x 1e5, unexplained fraction")
    segs = []
    for (i, c, t, sg, j) in sorted(q2.tap_table(), key=lambda z: (z[1], z[2])):
        m = np.mean([sg * np.r_[np.zeros(H), k[f][c], np.zeros(H)][t:t + 2 * H + 1] for f in FOLDS], axis=0)
        segs.append(m)
        res = m - (m @ bi) * bi - (m @ di) * di
        out.append(f"  loop {i} {'LR'[c]} lag {q1.S0 + t:5d} sign {sg:+d}: gain {1e5 * m @ bi:+6.2f}  delay {1e5 * m @ di:+6.2f}  rest {res @ res / (m @ m):.2f}")
    m = np.mean(segs, axis=0)
    out.append(f"mean over taps: gain {1e5 * m @ bi:+.2f}, delay {1e5 * m @ di:+.2f}, unexplained {1 - ((m @ bi) ** 2 + (m @ di) ** 2) / (m @ m):.3f}")
    for c, nm in ((0, "L"), (1, "R")):
        e = np.mean([k[f][c] ** 2 for f in FOLDS], axis=0)
        near = np.zeros(n, bool)
        for (i, cc_, t, sg, j) in q2.tap_table():
            if cc_ == c:
                near[max(0, t - 15):t + 16] = True
        out.append(f"{nm}: energy within +-15 lags of a tap {e[near].sum() / e.sum():.2f} (those lags are {near.mean():.2f} of the window)")
    open(os.path.join(OUT, f"decode{tag('')}.log"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def decode_loops(lo=177, hi=354):
    """The --base fit's increment (theta_*: dense, per loop and channel) read in the lo-hi band: fold agreement, the part
    common to the four loops, the lag profile (first-pass taps, one loop length on, the rest), and its ratio to the loop's
    own excitation in frequency (a filter on the loop's input or a time pattern)."""
    set_win()
    sos = butter(4, (lo, hi), btype="band", fs=FS, output="sos")
    bp = lambda v: sosfiltfilt(sos, np.r_[np.zeros(4000), v, np.zeros(4000)], padtype=None)[4000:4000 + len(v)]
    cc = lambda u, v: float(u @ v / np.sqrt(u @ u * v @ v))
    Lx = [int(P_) + q1.EXT for P_ in q1.P]
    cut = lambda x: [[x[sum(Lx) * c + sum(Lx[:i]):sum(Lx) * c + sum(Lx[:i + 1])] for i in range(4)] for c in (0, 1)]
    d = {f: [[bp(v) for v in row] for row in cut(np.load(os.path.join(OUT, f"theta_{tag(f)}.npy")))] for f in FOLDS}
    x0 = cut(x_base())
    taps = q2.tap_table()
    out = [f"## Increment {tag('')} over the base, per loop and channel, read in {lo}-{hi} Hz (lags from {q1.S0}; loop lengths {q1.P.tolist()})"]
    tot = sum((d[f][c][i] ** 2).sum() for f in FOLDS for c in (0, 1) for i in range(4))
    out.append("loop ch | share | fold A vs B | vs mean of the other loops | energy: taps +-15 / taps + P +-15 / lags < P / lags >= P")
    for c in (0, 1):
        for i in range(4):
            m = np.mean([d[f][c][i] for f in FOLDS], axis=0)
            n = min(Lx)
            oth = np.mean([np.mean([d[f][c][k][:n] for f in FOLDS], axis=0) for k in range(4) if k != i], axis=0)
            e = m ** 2
            near, nearP = np.zeros(len(m), bool), np.zeros(len(m), bool)
            for (ii, cc_, t, sg, j) in taps:
                if ii == i and cc_ == c:
                    near[max(0, t - 15):t + 16] = True
                    tp = t + int(q1.P[i])
                    if tp < len(m):
                        nearP[max(0, tp - 15):tp + 16] = True
            Pi = int(q1.P[i]) - q1.S0                                # excitation index of lag P_i
            sh = sum((d[f][c][i] ** 2).sum() for f in FOLDS) / tot
            out.append(f"  {i} {'LR'[c]} | {sh:.2f} | {cc(d['A'][c][i], d['B'][c][i]):+.2f} | {cc(m[:n], oth):+.2f} | "
                       f"{e[near].sum() / e.sum():.2f} ({near.mean():.2f}) / {e[nearP].sum() / e.sum():.2f} ({nearP.mean():.2f}) / "
                       f"{e[:Pi].sum() / e.sum():.2f} / {e[Pi:].sum() / e.sum():.2f}")
    out.append("\nRatio to the loop's own excitation (5 signed taps), mean of folds, at 180 / 220 / 260 / 300 / 340 Hz: |dX/X0| x 1e3 and phase (deg)")
    N = 1 << 15
    fr = np.fft.rfftfreq(N, 1 / FS)
    ks = [int(np.argmin(np.abs(fr - f0))) for f0 in (180, 220, 260, 300, 340)]
    for c in (0, 1):
        for i in range(4):
            X0 = np.fft.rfft(x0[c][i], N)
            r = [np.fft.rfft(np.mean([d[f][c][i] for f in FOLDS], axis=0), N)[k] / X0[k] for k in ks]
            rA = [np.fft.rfft(d["A"][c][i], N)[k] / X0[k] for k in ks]
            rB = [np.fft.rfft(d["B"][c][i], N)[k] / X0[k] for k in ks]
            out.append(f"  {i} {'LR'[c]}: " + "  ".join(f"{1e3 * abs(z):5.2f} {np.degrees(np.angle(z)):+5.0f}" for z in r) +
                       "   | A/B: " + " ".join(f"{1e3 * abs(a_):.1f}/{1e3 * abs(b_):.1f}" for a_, b_ in zip(rA, rB)))
    for c in (0, 1):                                                # where in lag the mean increment's energy sits (100-lag bins)
        e = np.mean([np.mean([d[f][c][i][:min(Lx)] for f in FOLDS], axis=0) ** 2 for i in range(4)], axis=0)
        b = e[:len(e) // 100 * 100].reshape(-1, 100).sum(1)
        out.append(f"{'LR'[c]} per-loop energy by 100-lag bin from lag {q1.S0} (% of total): " + " ".join(f"{100 * x / b.sum():.0f}" for x in b))
    open(os.path.join(OUT, f"decodeloops{tag('')}_{lo}-{hi}.log"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["prep", "fit", "score", "report", "decode", "decodeloops"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--fold", default="A", choices=sorted(FOLDS))
    ap.add_argument("--iters", type=int, default=200)
    ap.add_argument("--subset", default="")
    ap.add_argument("--through-g", action="store_true")
    ap.add_argument("--basis", default="dense")
    ap.add_argument("--band", default=None)
    ap.add_argument("--win", default=None)
    ap.add_argument("--project", default=None)
    ap.add_argument("--test-all", action="store_true")
    ap.add_argument("--base", default=None, help="a fitted correction's tag (after the subset) held fixed underneath, e.g. gain_band200-1500")
    a = ap.parse_args()
    global THROUGH_G, BASIS, FIT_BAND, SOS_FIT, WIN, PROJECT, TEST_ALL, BASE
    TEST_ALL, BASE = a.test_all, a.base
    PROJECT = tuple(int(x) for x in a.project.split(",")) if a.project else None
    WIN = tuple(int(x) for x in a.win.split(",")) if a.win else None
    THROUGH_G, BASIS = a.through_g, a.basis
    if a.band:
        FIT_BAND = tuple(int(x) for x in a.band.split(","))
        SOS_FIT = butter(3, FIT_BAND, btype="band", fs=FS, output="sos")
    SUBSET.update({c: float(v) for c, v in (x.split("=") for x in a.subset.split(",") if x)})
    os.makedirs(os.path.join(OUT, "prep"), exist_ok=True)
    if a.stage == "prep":
        info = q2.q1_tables()[0]
        keys = sorted(k for k in info if a.only is None or k in a.only.split(","))
        for i, k in enumerate(keys):
            prep(k)
            print(f"  {i + 1}/{len(keys)} {k}", flush=True)
    elif a.stage == "fit":
        fit(a.fold, a.iters)
    elif a.stage == "score":
        score(a.fold)
    elif a.stage == "decode":
        decode()
    elif a.stage == "decodeloops":
        decode_loops(*(int(x) for x in (a.project or "177,354").split(",")))
    else:
        report()


if __name__ == "__main__":
    main()
