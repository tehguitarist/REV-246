"""Is the baseline's common residual an output-side transfer error? Per channel a FIR h at lags -K..K with
residual ~ h * model, least squares over the training RT fold's captures (every setting, damped too; t0 .. min(te, 600 ms)),
scored on the other fold: octave nulls per window, baseline -> + h, against the next AMPL capture (excess).
model_s_results §25.12. h per fold -> analysis/dense_mid/<kernel>/resid/h_<fold>_K<K>_t<t0>.npy. t0 must leave K samples
before the window (t0 >= 12 ms at K = 300).

    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 6 -- .venv/bin/python tools/resid_common/outfir.py <K> <t0 ms>"""
import collections, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt
import dense_mid as dm
FS, PRE = dm.FS, dm.PRE
K, t0 = int(sys.argv[1]), float(sys.argv[2])            # half-length, fit start (ms)
band = butter(4, (100, 10000), btype="band", fs=FS, output="sos")
info = dm.q2.q1_tables()[0]
ap = json.load(open(os.path.join(dm.q2.OUT, dm.ls.lb.lf_.AP230)))
bstep = lambda k: ap["steps"][ap["settings"][k]["step"]]
FOLDS = dm.FOLDS
def caps(k):
    z = np.load(os.path.join(dm.OUT, "prep", k + ".npz"))
    out = []
    for j in sorted(int(n[3:]) for n in z.files if n.startswith("res")):
        out.append(dict(y=z[f"y{j}"].astype(float), m=z[f"m{j}"].astype(float), te=float(z[f"te{j}"]), amp=int(z[f"amp{j}"])))
    return out
def lagmat(x, i0, i1):
    return np.stack([x[i0 - s:i1 - s] for s in range(-K, K + 1)], 1)
keys = sorted(p[:-4] for p in os.listdir(os.path.join(dm.OUT, "prep")))
H = {}
for f, rts in FOLDS.items():
    AtA, Atb = np.zeros((2, 2 * K + 1, 2 * K + 1)), np.zeros((2, 2 * K + 1))
    for k in keys:
        if info[k]["rt"] not in rts:
            continue
        for c_ in caps(k):
            n = len(c_["y"])
            i0, i1 = PRE + int(t0 / 1000 * FS), min(PRE + int(min(c_["te"], 0.6) * FS), n - K)
            if c_["amp"] == 0 and t0 < 50:
                continue                                     # AMPL 0 has no cluster
            for c in (0, 1):
                m = sosfiltfilt(band, c_["m"][:, c]); r = sosfiltfilt(band, c_["y"][:, c] - c_["m"][:, c])
                A = lagmat(m, i0, i1)
                AtA[c] += A.T @ A; Atb[c] += A.T @ r[i0:i1]
    H[f] = np.array([np.linalg.solve(AtA[c] + 1e-9 * np.trace(AtA[c]) / len(AtA[c]) * np.eye(len(AtA[c])), Atb[c]) for c in (0, 1)])
    np.save(os.path.join(dm.OUT, "resid", f"h_{f}_K{K}_t{t0:g}.npy"), H[f])
cc = lambda u, v: float(u.ravel() @ v.ravel() / np.sqrt(u.ravel() @ u.ravel() * v.ravel() @ v.ravel()))
print(f"K {K}, fit from {t0:g} ms: h fold A vs B {cc(H['A'], H['B']):+.3f}; L vs R {cc(H['A'][0], H['A'][1]):+.3f}")
w = np.linspace(0, np.pi, 4097); fr = w / np.pi * FS / 2
for f in "AB":
    for c in (0, 1):
        Hw = np.exp(-1j * np.outer(w, np.arange(-K, K + 1))) @ H[f][c]
        print(f"  {f} {'LR'[c]}: 1+H at " + "  ".join(f"{x}Hz {100 * (abs(1 + Hw[np.argmin(abs(fr - x))]) - 1):+.3f}% {np.degrees(np.angle(1 + Hw[np.argmin(abs(fr - x))])):+.3f}deg" for x in (125, 250, 500, 1000, 2000, 4000, 8000, 10000)))
# score held out
OCT = [125, 250, 500, 1000, 2000, 4000]
W = [(0, 50), (50, 150), (150, 300), (300, 600), (600, 1000)]
osos = {fc: butter(3, (fc / 1.41, fc * 1.41), btype="band", fs=FS, output="sos") for fc in OCT}
rows = []
grp = lambda k: "undamped" if bstep(k) < 0.01 else "damped"
for f, rts in FOLDS.items():
    h = H["B" if f == "A" else "A"]                          # scored with the fold that did not train on these RT steps
    for k in keys:
        if info[k]["rt"] not in rts:
            continue
        cs = caps(k)
        for j, c_ in enumerate(cs):
            cb = cs[(j + 1) % len(cs)]
            n = min(len(c_["y"]), len(cb["y"]))
            te = min(c_["te"], cb["te"])
            y, m, yb = c_["y"][:n, 1], c_["m"][:n, 1], cb["y"][:n, 1]
            m1 = m + sum(h[1][s + K] * np.r_[np.zeros(max(s, 0)), m[:n - s] if s >= 0 else m[-s:], np.zeros(max(-s, 0))][:n] for s in range(-K, K + 1))   # m + h * m
            for (a, b) in W:
                i0, i1 = PRE + int(a / 1000 * FS), min(PRE + int(b / 1000 * FS), PRE + int(te * FS), n)
                if i1 - i0 < int(0.03 * FS) or (c_["amp"] == 0 and a == 0):
                    continue
                for fc in OCT + [0]:
                    flt = (lambda v: sosfiltfilt(osos[fc], v)) if fc else (lambda v: sosfiltfilt(band, v))
                    Y, M0, M1, YB = (flt(v)[i0:i1] for v in (y, m, m1, yb))
                    e = lambda u, v: 10 * np.log10(np.sum((u - v) ** 2) / np.sum(u ** 2))
                    rows.append((grp(k), a, fc, e(Y, M0), e(Y, M1), e(Y, YB)))
D = collections.defaultdict(list)
for g, a, fc, x0, x1, xr in rows:
    D[(g, a, fc)].append((x0, x1, xr)); D[("all", a, fc)].append((x0, x1, xr))
for g in ("all", "undamped", "damped"):
    print(f"\n{g}: R, median null baseline -> + h [excess over next AMPL]")
    for (a, b) in W:
        cells = []
        for fc in OCT + [0]:
            v = np.array(D.get((g, a, fc), []))
            if len(v):
                cells.append(f"{fc or 'bb'}: {np.median(v[:, 0]):.1f}->{np.median(v[:, 1]):.1f} [{np.median(v[:, 0] - v[:, 2]):+.1f}->{np.median(v[:, 1] - v[:, 2]):+.1f}]")
        print(f"  {a}-{b}: " + "; ".join(cells))
