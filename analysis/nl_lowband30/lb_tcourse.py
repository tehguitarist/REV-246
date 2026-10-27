import pickle, collections, numpy as np
from scipy.signal import butter, sosfiltfilt
FS = 26373.626373626373
S = pickle.load(open("low_sel.pkl", "rb"))
bd = {"20-100": butter(4, [20, 100], btype="band", fs=FS, output="sos"), "100-1k": butter(4, [100, 1000], btype="band", fs=FS, output="sos")}
edges = [0, 0.01, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
print("local null per time bin (dB, e vs y energy in the bin; t re front zero), median over captures (LATE<=1, SPAN 0.5), L/R mean")
for name, sos in bd.items():
    print(name, " bins:", edges)
    for grp in ("fs0", "fs23"):
        for early in sorted({k[1] for k in S}):
            ks = [k for k, d in S.items() if k[1] == early and (abs(d["fs"]) < .5) == (grp == "fs0")]
            if not ks: continue
            row = []
            for t0, t1 in zip(edges[:-1], edges[1:]):
                v = []
                for k in ks:
                    d = S[k]; a, b = int(d["z"] + t0 * FS), int(d["z"] + t1 * FS)
                    nn = 0.0
                    for c in (0, 1):
                        yb, mb = sosfiltfilt(sos, d["y"][:, c]), sosfiltfilt(sos, d["m"][:, c])
                        v.append(10 * np.log10(np.sum((yb[a:b] - mb[a:b]) ** 2) / np.sum(yb[a:b] ** 2)))
                row.append(np.median(v))
            print(f"  {grp} E{early:<4} n{len(ks):2d} " + " ".join(f"{x:6.1f}" for x in row))
