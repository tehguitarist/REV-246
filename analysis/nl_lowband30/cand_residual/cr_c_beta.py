"""Fitted values of the front-section perturbations (all captures, per channel) and the per-EARLY-cell values, 20-100 Hz and 100-1000 Hz."""
from cr_c_fit import *
G = CAND["front taps gain+dt (6)"]; nm = ["fr0_gain", "fr1_gain", "fr2_gain", "fr0_dt", "fr1_dt", "fr2_dt"]
G3 = CAND["front gain+dt+rdir_dt (3)"]
for bd in (20, 100):
    print(f"\nband {bd}: front taps gain (relative) / dt (samples), train = all SPAN 0.5 captures of the EARLY; L then R")
    for e_ in [None] + EL:
        tr = E5 if e_ is None else [k for k in E5 if k[1] == e_]
        out = []
        for ch in (0, 1):
            Aw = np.zeros((6, 6)); bw = np.zeros(6)
            for k in tr:
                A, b, c = reduce(D[k], ch, bd, G, False); w = 1 / D[k][f"yty_{ch}_{bd}"]; Aw += w * A; bw += w * b
            beta = np.linalg.solve(Aw + 1e-7 * np.trace(Aw) / 6 * np.eye(6), bw); out.append(beta)
        lab = "all " if e_ is None else f"E{e_:<3g}"
        print(f"  {lab}: " + " | ".join("g " + " ".join(f"{1e3*v:+6.2f}" for v in o[:3]) + "e-3  dt " + " ".join(f"{1e4*v:+6.2f}" for v in o[3:]) + "e-4" for o in out))
