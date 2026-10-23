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
        files.append(dict(y=z[f"res{j}"][:e].astype(np.float64), w=w[:e], amp=0))
    R = list(z["R"].astype(np.complex128))
    if THROUGH_G:
        R = R + [r * Gi for r, Gi in zip(R, loop_G(key, int(z["N"])))]
    g = dict(N=int(z["N"]), R=[r * z["Bm"] for r in R], PB=np.zeros(int(z["N"]) // 2 + 1), files=files)
    if full:
        g["R_full"] = R
        g["caps"] = [dict(y=z[f"y{j}"].astype(np.float64), m=z[f"m{j}"].astype(np.float64), g=z[f"g{j}"].astype(np.float64),
                          amp=int(z[f"amp{j}"]), te=float(z[f"te{j}"])) for j in k]
    g["lab"] = z["lab"].tolist()
    return g


SUBSET = {}                                                         # --subset lf=1.0,hf=1.0: restrict train and test


def keys_for(rts):
    """The prepared settings at those RT steps (within SUBSET)."""
    info = q2.q1_tables()[0]
    ap = json.load(open(os.path.join(q2.OUT, ls.lb.lf_.AP230)))
    b0 = lambda k: ap["steps"][ap["settings"][k]["step"]] < 0.01                     # AP230 without damping
    return sorted(k for k in info if info[k]["rt"] in rts and os.path.exists(os.path.join(OUT, "prep", k + ".npz"))
                  and all(info[k][c] == v for c, v in SUBSET.items() if c not in ("b0", "rtmin"))
                  and (not SUBSET.get("b0") or b0(k)) and info[k]["rt"] >= SUBSET.get("rtmin", 0))


def tag(fold):
    return fold + ("_" + "_".join(f"{c}{v:g}" for c, v in sorted(SUBSET.items())) if SUBSET else "") + ("_G" if THROUGH_G else "")


def fit(fold, iters):
    from scipy.sparse.linalg import LinearOperator, lsqr
    q1.EXT = 1700
    groups = [load(k) for k in keys_for(FOLDS[fold])]
    J = q1.Joint(groups)
    off = 8 + 2 * q1.LFIR
    nx = J.size() - off
    fwd = lambda x: J.forward(np.r_[np.zeros(off), x])
    adj = lambda z: J.adjoint(z)[off:]
    tgt = J.target()
    op = LinearOperator((len(tgt), nx), matvec=fwd, rmatvec=adj, dtype=float)
    x = lsqr(op, tgt, iter_lim=iters, atol=1e-12, btol=1e-12)[0]
    res = tgt - fwd(x)
    print(f"fold {fold}: {len(groups)} settings, {len(tgt)} samples, {nx} unknowns; training residual "
          f"{10 * np.log10(np.sum(res ** 2) / np.sum(tgt ** 2)):+.2f} dB re the baseline's (80-1500 Hz, weighted)", flush=True)
    np.save(os.path.join(OUT, f"x_{tag(fold)}.npy"), x)


def score(fold):
    q1.EXT = 1700
    x = np.load(os.path.join(OUT, f"x_{tag(fold)}.npy"))
    test = [r for r in sorted({rt for f_ in FOLDS.values() for rt in f_}) if r not in FOLDS[fold]]
    os.makedirs(os.path.join(OUT, f"score_{tag(fold)}"), exist_ok=True)
    for k in keys_for(test):
        path = os.path.join(OUT, f"score_{tag(fold)}", k + ".json")
        if os.path.exists(path):
            continue
        g = load(k, full=True)
        J = q1.Joint([dict(g, R=g["R_full"])])
        _, _, xs = J.split(np.r_[np.zeros(8 + 2 * q1.LFIR), x])
        F, X = J._spectra(np.zeros((2, q1.LFIR)), xs, g["N"])
        add = np.fft.irfft(np.einsum("cik,ik->ck", X, np.array(g["R_full"])), g["N"], axis=1).T
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
            bb = ls.nulls_broad(c["y"][:n], c["m"][:n], te), ls.nulls_broad(c["y"][:n], m1, te)
            rows.append(dict(oct=res, broad=bb))
        json.dump(dict(key=k, lab=g["lab"], rows=rows), open(path, "w"))
        print(k, flush=True)


def report():
    S = [json.load(open(p)) for f_ in FOLDS for p in glob.glob(os.path.join(OUT, f"score_{tag(f_)}", "*.json"))]
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
    open(os.path.join(OUT, f"report{tag('')}.log"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["prep", "fit", "score", "report"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--fold", default="A", choices=sorted(FOLDS))
    ap.add_argument("--iters", type=int, default=200)
    ap.add_argument("--subset", default="")
    ap.add_argument("--through-g", action="store_true")
    a = ap.parse_args()
    global THROUGH_G
    THROUGH_G = a.through_g
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
    else:
        report()


if __name__ == "__main__":
    main()
