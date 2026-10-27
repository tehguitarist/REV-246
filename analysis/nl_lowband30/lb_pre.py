import pickle, numpy as np
np.set_printoptions(linewidth=200, precision=5, suppress=True)
D = pickle.load(open("time.pkl", "rb"))
seen = set()
for (sp, fs, k), rec in D.items():
    if (sp, fs) in seen: continue
    seen.add((sp, fs))
    r = rec[0]
    print(sp, fs, k)
    print("  y[0:90:6]*1e4 ", (r["y"][0:90:6] * 1e4))
    print("  m[0:90:6]*1e4 ", (r["m"][0:90:6] * 1e4))
    print("  e[0:90:6]*1e4 ", (r["e"][0:90:6] * 1e4))
    print("  mean y[0:30]*1e4 %.2f  mean m[0:30] %.2f" % (r["y"][:30].mean() * 1e4, r["m"][:30].mean() * 1e4))
