"""Is the §25.11 baseline's 250 Hz residual a fixed waveform across settings? Per capture (prep: y - model, clocks and
gains fixed), 250 Hz octave, R; correlation of each capture's residual with the mean residual of the OTHER settings
(leave-one-setting-out), per window; repeat captures (other AMPL) give the ceiling. Groups: undamped / damped.
model_s_results §25.12.

    EMT_KERNEL=chain3 .venv/bin/python tools/resid_common/fixed.py [lo,hi Hz]"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt
import dense_mid as dm
FS, PRE = dm.FS, dm.PRE
lo, hi = (int(x) for x in sys.argv[1].split(",")) if len(sys.argv) > 1 else (177, 354)
sos = butter(4, (lo, hi), btype="band", fs=FS, output="sos")
info = dm.q2.q1_tables()[0]
ap = json.load(open(os.path.join(dm.q2.OUT, dm.ls.lb.lf_.AP230)))
bstep = lambda k: ap["steps"][ap["settings"][k]["step"]]
W = [(0, 30), (30, 50), (50, 75), (75, 100), (100, 125), (125, 150), (150, 200), (200, 300)]
n_end = PRE + int(0.3 * FS)
R, Y, lab = {}, {}, {}
for p in sorted(os.listdir(os.path.join(dm.OUT, "prep"))):
    k = p[:-4]
    z = np.load(os.path.join(dm.OUT, "prep", p))
    caps = sorted(int(n[3:]) for n in z.files if n.startswith("res"))
    rr, yy = [], []
    for j in caps:
        y, m = z[f"y{j}"][:, 1].astype(float), z[f"m{j}"][:, 1].astype(float)
        if len(y) < n_end:
            continue
        rr.append(sosfiltfilt(sos, (y - m)[:n_end + 2000])[:n_end]); yy.append(sosfiltfilt(sos, y[:n_end + 2000])[:n_end])
    if rr:
        R[k], Y[k], lab[k] = rr, yy, (info[k]["lf"], info[k]["rt"], info[k]["hf"], bstep(k))
und = [k for k in R if lab[k][3] < 0.01]
dam = [k for k in R if lab[k][3] >= 0.01]
cc = lambda u, v: float(u @ v / np.sqrt(u @ u * v @ v))
print(f"{lo}-{hi} Hz, R; {len(und)} undamped / {len(dam)} damped settings; median over captures")
print("window ms | resid vs mean resid of other undamped settings (undamped caps / damped caps) | resid vs other AMPL resid (same setting) | resid / signal energy (dB)")
for a, b in W:
    i0, i1 = PRE + int(a / 1000 * FS), PRE + int(b / 1000 * FS)
    tot = {k: np.mean([r[i0:i1] for r in R[k]], axis=0) for k in R}
    S = sum(tot[k] for k in und)
    cu = [cc(r[i0:i1], (S - tot[k]) / (len(und) - 1)) for k in und for r in R[k]]
    cd = [cc(r[i0:i1], S / len(und)) for k in dam for r in R[k]]
    cr = [cc(R[k][j][i0:i1], R[k][(j + 1) % len(R[k])][i0:i1]) for k in und for j in range(len(R[k])) if len(R[k]) > 1]
    e = [10 * np.log10((r[i0:i1] @ r[i0:i1]) / (y[i0:i1] @ y[i0:i1])) for k in und for r, y in zip(R[k], Y[k])]
    print(f"  {a:4d}-{b:<4d} | {np.median(cu):+.2f} / {np.median(cd):+.2f} | {np.median(cr):+.2f} | {np.median(e):6.1f}")

print("\namplitude of each setting's residual on the leave-one-out mean waveform (undamped mean), median by RT / by AP230 step b, per window")
rts = sorted({lab[k][1] for k in R}); bs = sorted({round(lab[k][3], 2) for k in R})
for a, b in W[1:]:
    i0, i1 = PRE + int(a / 1000 * FS), PRE + int(b / 1000 * FS)
    tot = {k: np.mean([r[i0:i1] for r in R[k]], axis=0) for k in R}
    S = sum(tot[k] for k in und)
    amp = {}
    for k in R:
        ref = (S - tot[k]) / (len(und) - 1) if k in und else S / len(und)
        amp[k] = float(tot[k] @ ref / (ref @ ref))
    byrt = " ".join(f"{rt:g}:{np.median([amp[k] for k in und if lab[k][1] == rt]):.2f}" for rt in rts if any(lab[k][1] == rt for k in und))
    byb = " ".join(f"{bb:g}:{np.median([amp[k] for k in R if round(lab[k][3], 2) == bb]):.2f}" for bb in bs)
    bylf = " ".join(f"{lf:g}:{np.median([amp[k] for k in und if lab[k][0] == lf]):.2f}" for lf in (0.1, 0.4, 1.0, 2.0))
    byhf = " ".join(f"{hf:g}:{np.median([amp[k] for k in und if lab[k][2] == hf]):.2f}" for hf in (0.2, 0.4, 0.6, 1.0))
    print(f"  {a:4d}-{b:<4d} undamped by RT {byrt} | by b {byb} | LF {bylf} | HF {byhf}")
