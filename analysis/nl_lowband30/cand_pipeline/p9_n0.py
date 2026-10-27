import sys, json; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=220, precision=2, suppress=True)
P4 = json.load(open("/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline/p4_res.json"))
N0S = (16, 30, 60, 100, 150)
rows = {n: [] for n in N0S}; cm = []
for ks, v in P4.items():
    key = tuple(float(a) for a in ks.split("_")); r = align(key)
    x, _ = sf.read(file_path(key), always_2d=True); nf = len(x); cache = nl.ir(key); nlti = len(cache)
    mp = M.render(key, nlti, r["o"], r["front_shift"], band=False, proc=True) * np.array(r["gain"])
    z = int(round(K0 + r["front_shift"])); u = 1e4 / np.abs(mp[z:z + 800, 0]).max()
    k, _ = pipeline(np.ones((nf, 2))); kpre = k[8:40].mean(0) * u
    cm.append(v["c_meas"])
    for n in N0S: rows[n].append(x[:n].mean(0) * kpre)
cm = np.array(cm).ravel()
for n in N0S:
    c = np.array(rows[n]).ravel(); print(f"n0 {n:4d}: corr {np.corrcoef(cm, c)[0,1]:.4f} slope {(cm@c)/(c@c):.3f} rms resid {np.sqrt(((cm-c)**2).mean()):.3f} (rms c_meas {np.sqrt((cm**2).mean()):.2f})")
