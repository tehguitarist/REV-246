import sys, pickle
from fitlib2 import *
out, slo, shi, lmax, emin = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
stride = int(sys.argv[6]) if len(sys.argv) > 6 else 1; off = int(sys.argv[7]) if len(sys.argv) > 7 else 0
B, ks = catalog(slo, shi, lmax, emin); ks = ks[off::stride]
path = os.path.join(HERE, out)
D = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
t0 = time.time()
for i, k in enumerate(ks):
    if k in D: continue
    D[k] = dict(sess=sess(B[k]), r=B[k], g=build2(k, B[k]))
    if i % 10 == 0:
        pickle.dump(D, open(path, "wb")); print(i, len(ks), "%.0f s" % (time.time() - t0), flush=True)
pickle.dump(D, open(path, "wb")); print("done", len(D))
