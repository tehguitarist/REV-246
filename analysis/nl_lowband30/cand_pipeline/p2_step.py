import sys; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=220, precision=2, suppress=True)
key = tuple(float(v) for v in sys.argv[1:4]); r = align(key)
fade = None if len(sys.argv) < 5 else tuple(float(v) for v in sys.argv[4:6])
x, fs = sf.read(file_path(key), always_2d=True); nf = len(x)
cache = nl.ir(key)
print(key, "file s", nf / FS96, "o %.3f fr %.2f ppm %.2f" % (r["o"], r["front_shift"], r["ppm"]), "fade", fade)
xs0 = make_x(key, r, nf, 1.0)
A = fit_scale(x, xs0); print("A", A)
x_tr = make_x(key, r, nf, A, 0.0, fade)
x_long = make_x(key, r, nf, A, 6.0)
h_tr, _ = pipeline(x_tr); h_lg, _ = pipeline(x_long)
nlti = len(cache)
mp = M.render(key, nlti, r["o"], r["front_shift"], band=False, proc=True) * np.array(r["gain"])
lti = forward_warp(mp, r["ppm"])
z = int(round(K0 + r["front_shift"])); pk = np.abs(mp[z:z + 800, 0]).max()
u = 1e4 / pk
idx = np.arange(34, 70, 2)
print("index", idx)
sl = slice(8, 40)
for name, a, b in (("real-LTI", cache, lti), ("sim_trunc-LTI", h_tr[:nlti], lti), ("sim_long-LTI", h_lg[:nlti], lti), ("sim_trunc-sim_long", h_tr[:nlti], h_lg[:nlti])):
    e = (a[:nlti] - b[:nlti]) * u
    print(f"{name:20s} c(8:40) L %.2f R %.2f   (1e-4 of front-tap peak)" % tuple(e[sl].mean(0)), " | e after idx 60-200 rms: %.2f %.2f" % tuple(np.sqrt((e[60:200] ** 2).mean(0))))
    for c in (0, 1): print("   ch", c, e[idx, c])
# 96k domain: real pre-roll in d
print("---- residual file x_real - x_sim through the pipeline")
xr = x - x_tr
h_r, d_r = pipeline(xr)
e = h_r[:nlti] * u
print("c(8:40) L %.2f R %.2f" % tuple(e[sl].mean(0)))
for c in (0, 1): print("  ch", c, e[idx, c])
print("pipeline(real)-pipeline(sim_trunc) == real - sim_trunc  -> cross-check with  real - sim_trunc:")
e2 = (cache - h_tr[:nlti]) * u
print("c(8:40) L %.2f R %.2f" % tuple(e2[sl].mean(0)))
for c in (0, 1): print("  ch", c, e2[idx, c])
print("residual file stats: mean first150 %s, mean all %s, rms all %s" % (xr[:150].mean(0), xr.mean(0), np.sqrt((xr**2).mean(0))))
# spectrum of residual file: DC / sub-20 Hz and Nyquist tone
Xr = np.fft.rfft(xr, axis=0); f = np.fft.rfftfreq(len(xr), 1 / FS96)
for lo, hi in ((0, 1), (1, 3), (3, 12), (12, 50), (50, 300), (300, 3000), (3000, 12800), (12800, 47000), (47990, 48001)):
    m_ = (f >= lo) & (f < hi); print(f"  band {lo}-{hi} Hz  rms amp per ch: %s" % (np.sqrt((np.abs(Xr[m_]) ** 2).sum(0)) * np.sqrt(2) / len(xr)))
# time domain: first 400 samples of the residual, mean of even/odd
print("xr[0:6]", xr[:6, 0], " pre-roll (0:170) mean %.3e, alt-amp %.3e" % (xr[:170, 0].mean(), (xr[0:170:2, 0].mean() - xr[1:170:2, 0].mean()) / 2))
