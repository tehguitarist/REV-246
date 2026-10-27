from solve import *
import collections
D = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb")); v = sys.argv[2]
print(len(D), collections.Counter((d["sess"], k[0]) for k, d in D.items()))
np.set_printoptions(linewidth=200)
for ch in (0, 1):
    print("=== ch", ch, "variant", v)
    for s in sorted(set(d["sess"] for d in D.values())):
        ks = sel(D, sess=[s]); G, b, yy, n = pool(D, ks, v, ch)
        c = solve(G, b); T = Tof(c)
        print(f"fs {s:3d} n{len(ks):3d}  ", fmt(T, range(0, 11)), "  gain %.2f dB" % gain_db(G, b, yy, c))
    ks = sel(D); G, b, yy, n = pool(D, ks, v, ch); c = solve(G, b); print("all      ", fmt(Tof(c), range(0, 11)))
