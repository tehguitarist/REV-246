import sys; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=200, precision=3)
key = tuple(float(v) for v in sys.argv[1:4])
r = align(key); print(key, {k: r[k] for k in ("o", "ppm", "gain", "front_shift", "trimmed")})
x, fs = sf.read(file_path(key), always_2d=True); nf = len(x)
cache = nl.ir(key); print("cache", cache.shape, "file", x.shape)
nn = int(nf / RATIO) + 40
m = model_cache(key, nn, r)
raw = forward_warp(m, r["ppm"])
dc = to_capture_rate(raw)
xs = chain(dc, nf)
print("x_sim shape", xs.shape)
# scale-free comparison of x_real and x_sim: least squares gain per channel over active region
for c in (0, 1):
    g = (x[:, c] @ xs[:, c]) / (xs[:, c] @ xs[:, c]); e = x[:, c] - g * xs[:, c]
    print("ch", c, "gain", g, "null", 10 * np.log10((e @ e) / (x[:, c] @ x[:, c])))
    # windows of 0.1 s
    seg = [(int(a), int(a + 4800)) for a in range(0, nf - 4800, 9600)]
    print("  windowed null dB:", np.array([10 * np.log10(((x[a:b, c] - g * xs[a:b, c]) ** 2).sum() / (x[a:b, c] ** 2).sum()) for a, b in seg]))
    print("  first 8 real", x[:8, c], " sim", g * xs[:8, c])
    print("  mean first 150: real %.3e sim %.3e" % (x[:150, c].mean(), g * xs[:150, c].mean()))
    print("  last 8 real", x[-8:, c], "sim", g * xs[-8:, c])
