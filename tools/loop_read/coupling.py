"""Is there a small coupling between S's four loops? (NEXT_STEPS §0 item 6; model_s_results §31.)

On the §29 baseline (fold terms of analysis/loop_cv/cap0.8475_folds, the hfrt cells, EMT_PMAX 0.8475, EMT_HFRT_MIN
0.12), per setting and clean capture (tools/defects.py), the residual capture - model (fade-compensated: model x ramp)
is regressed on the first-order terms a coupling would add. A coupling eps_ji (loop i's output fed into loop j's input)
adds, to first order, loop j's read-out path driven by loop i's recirculating signal:

    b_ij = X_j R_j . G_i Z_i / (1 - G_i Z_i)     (X_j: loop j's excitation, R_j = PB / (1 - G_j Z_j), Z_i = z^-P_i)

12 cross terms; the 4 self terms (i = j: a per-pass gain error of loop i) as a control. Per setting the normal
equations per capture, channel and window (50-150, 150-300, 300-600 ms; 177 Hz - 4 kHz) are stored (prep, one file per
setting, resumable); fit chooses the coefficients leaving one RT step out (one set for all settings, L and R shared)
and scores the held-out settings' nulls from the stored equations.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds \\
        SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_hfrt/cells.json \\
        .venv/bin/python tools/loop_read/coupling.py --stage prep [--only keys] [--workers 8]
    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/coupling.py --stage fit       # -> analysis/lowband/chain3_coupling.{json,log}"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "lowband", "chain3_coupling")
WIN = [(0.05, 0.15), (0.15, 0.30), (0.30, 0.60)]
PAIRS = [(i, j) for j in range(4) for i in range(4)]        # (i, j): loop i's output into loop j; i == j: self terms


def prep_one(a):
    key, path = a
    import lowband_shelf as ls
    from defects import DEFECTS
    q2, q1 = ls.q2, ls.q1
    FS, PRE = q1.FS, q1.PRE
    d = ls.bb_prep(key)
    w, N, PB = d["w"], d["N"], d["G"]["PB"]
    Z = [np.exp(-1j * w * Pi) for Pi in q1.P]
    R = [PB / (1 - Gi * Zi) for Gi, Zi in zip(d["Gs0"], Z)]
    sos = butter(4, [177, 4000], btype="band", fs=FS, output="sos")
    basis = []
    for i, j in PAIRS:
        S = d["X"][:, j] * (R[j] * d["Gs0"][i] * Z[i] / (1 - d["Gs0"][i] * Z[i]))[None]    # R_i / PB, without dividing by PB
        basis.append(np.fft.irfft(S, N, axis=1).T)        # (N, 2)
    r0 = ls.reverb(d, d["Gs0"])
    out = dict(key=key, **d["lab"], caps=[])
    for c in d["caps"]:
        if (key, c["amp"]) in DEFECTS:
            continue
        n = len(c["y"])
        g = c["g"][:n]
        res = sosfiltfilt(sos, c["y"] - ls.model(d, c, r0) * g, axis=0)
        yb = sosfiltfilt(sos, c["y"], axis=0)
        B = np.stack([sosfiltfilt(sos, b[:n] * g, axis=0) for b in basis])        # (16, n, 2)
        cap = dict(amp=c["amp"], win={})
        for t0, t1 in WIN:
            i0, i1 = PRE + int(t0 * FS), min(PRE + int(t1 * FS), PRE + int(c["te"] * FS), n)
            if i1 - i0 < int(0.04 * FS):
                continue
            cap["win"][f"{int(t0 * 1000)}-{int(t1 * 1000)}"] = [dict(
                A=(B[:, i0:i1, ch] @ B[:, i0:i1, ch].T).tolist(), c=(B[:, i0:i1, ch] @ res[i0:i1, ch]).tolist(),
                rr=float(res[i0:i1, ch] @ res[i0:i1, ch]), yy=float(yb[i0:i1, ch] @ yb[i0:i1, ch])) for ch in (0, 1)]
        out["caps"].append(cap)
    json.dump(out, open(path, "w"))
    print(f"{key}: {len(out['caps'])} clean captures", flush=True)
    return key


def prep_stage(only, workers):
    from concurrent.futures import ProcessPoolExecutor
    import q2_injection as q2
    os.makedirs(OUT, exist_ok=True)
    keys = sorted(q2.q1_tables()[0]) if not only else only
    jobs = [(k, os.path.join(OUT, k + ".json")) for k in keys if not os.path.exists(os.path.join(OUT, k + ".json"))]
    print(f"{len(jobs)} settings to prep", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(prep_one, jobs):
            pass


def fit_stage():
    D = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(OUT, "*.json")))]
    rts = sorted({d["rt"] for d in D})
    sets = {"cross (12)": [k for k, (i, j) in enumerate(PAIRS) if i != j], "self (4)": [k for k, (i, j) in enumerate(PAIRS) if i == j],
            "cross + self (16)": list(range(16))}
    lines = [f"# Loop coupling on the §29 baseline ({len(D)} settings, clean captures; 177 Hz - 4 kHz; coefficients chosen leaving one "
             "RT step out, L and R shared, every capture weighted by its window energy); R null median, baseline -> with the terms"]
    res = {}
    for name, idx in sets.items():
        per = {w: ([], []) for w in ("50-150", "150-300", "300-600")}
        coefs = {}
        for rt in rts:
            A = np.zeros((len(idx), len(idx))); cc = np.zeros(len(idx))
            for d in D:
                if d["rt"] == rt:
                    continue
                for cap in d["caps"]:
                    for wn, chs in cap["win"].items():
                        if wn != "150-300":                    # trained on 150-300 ms (the window with the excess)
                            continue
                        for e in chs:
                            s = 1.0 / e["yy"]
                            A += s * np.array(e["A"])[np.ix_(idx, idx)]; cc += s * np.array(e["c"])[idx]
            x = np.linalg.solve(A + 1e-12 * np.trace(A) * np.eye(len(idx)), cc)
            coefs[str(rt)] = x.tolist()
            for d in D:
                if d["rt"] != rt:
                    continue
                for cap in d["caps"]:
                    for wn, chs in cap["win"].items():
                        e = chs[1]
                        A_, c_ = np.array(e["A"])[np.ix_(idx, idx)], np.array(e["c"])[idx]
                        after = e["rr"] - 2 * x @ c_ + x @ A_ @ x
                        per[wn][0].append(10 * np.log10(e["rr"] / e["yy"])); per[wn][1].append(10 * np.log10(max(after, 1e-30) / e["yy"]))
        res[name] = dict(coefs=coefs, nulls={w: [float(np.median(a)), float(np.median(b))] for w, (a, b) in per.items() if a})
        lines.append(f"- {name}: " + "; ".join(f"{w}: {a:.2f} -> {b:.2f}" for w, (a, b) in res[name]["nulls"].items())
                     + " | coefficients (median over folds) " + " ".join(f"{np.median([coefs[r][k] for r in coefs]):+.4f}" for k in range(len(idx))))
    json.dump(res, open(OUT + ".json", "w"), indent=1)
    open(OUT + ".log", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["prep", "fit"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    if a.stage == "prep":
        prep_stage(a.only.split(",") if a.only else None, a.workers)
    else:
        fit_stage()


if __name__ == "__main__":
    main()
