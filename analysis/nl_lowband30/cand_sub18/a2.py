from solve import *
D = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb")); lam = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
rng = np.random.default_rng(0)
def est(ks, v, ch, lam):
    G, b, yy, n = pool(D, ks, v, ch); return Tof(solve(G, b, lam))
idx = [1, 2, 3, 4, 5, 6, 7, 8, 9]   # 4..20 Hz
print("knots", KNOTS[idx])
for v in ("w45", "span", "pre_feed", "w60", "to_end"):
    print("##### variant", v)
    for ch in (0, 1):
        for s in (0, 20, 23):
            ks = sel(D, sess=[s])
            ks = [k for k in ks if v in D[k]["g"]]
            if len(ks) < 4: continue
            T = est(ks, v, ch, lam)
            bs = np.array([est([ks[i] for i in rng.integers(0, len(ks), len(ks))], v, ch, lam) for _ in range(30)])
            sd = np.sqrt(np.mean(abs(bs - bs.mean(0)) ** 2, 0))
            print(f"ch{ch} fs{s:3d} n{len(ks):3d} T: " + " ".join(f"{20*np.log10(abs(T[i])):+5.1f}/{np.degrees(np.angle(T[i])):+4.0f}" for i in idx))
            print(f"{'':16s}sd(|dT|,dB re 1): " + " ".join(f"{20*np.log10(max(sd[i],1e-9)+1)-0:5.2f}   " for i in idx))
