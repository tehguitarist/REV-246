from cr_common import *
import resource
B = bench_all()
keys = [(0.5, e, l) for e in (1.0, 2.0) for l in (0.0, 0.2, 0.4)]
t = time.time()
S = {b: band_sos(*b) for b in ((20, 100), (100, 1000), (1000, 5000))}
for k in keys:
    r = B[k]; y, m = cached(k, r); a, b = window(k, r)
    out = []
    for bd, sos in S.items():
        yb = sosfiltfilt(sos, y, axis=0)[a:b]; mb = sosfiltfilt(sos, m, axis=0)[a:b]
        out.append(" ".join(f"{10*np.log10(fpow(yb[:,c]-mb[:,c])/fpow(yb[:,c])):6.1f}" for c in (0, 1)))
    print(k, f"o {r['o']:.3f} fs {r['front_shift']:.1f}", "| " + " | ".join(out), "| bench bands:", [round(v[0], 1) for v in (r["bands"]["20-100"], r["bands"]["100-1000"], r["bands"]["1000-5000"])], f"{time.time()-t:.1f}s")
print("maxrss MB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
