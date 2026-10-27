import sys; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=220, precision=3, suppress=False)
key = tuple(float(v) for v in sys.argv[1:4]); r = align(key)
x, fs = sf.read(file_path(key), always_2d=True); nf = len(x)
cache = nl.ir(key); nlti = len(cache)
xs0 = make_x(key, r, nf, 1.0); A = fit_scale(x, xs0)
x_tr = make_x(key, r, nf, A)
xr = x - x_tr
mp = M.render(key, nlti, r["o"], r["front_shift"], band=False, proc=True) * np.array(r["gain"])
z = int(round(K0 + r["front_shift"])); pk = np.abs(mp[z:z + 800, 0]).max(); u = 1e4 / pk
sl = slice(8, 40)
def run(tag, xx):
    h, _ = pipeline(xx); e = h[:nlti] * u
    print(f"{tag:42s} c L %6.2f R %6.2f | idx44 %5.2f %5.2f idx48 %5.2f %5.2f  idx60 %5.2f %5.2f" % (*e[sl].mean(0), *e[44], *e[48], *e[60]))
    return h
run("full residual", xr)
# low-pass (< 100 Hz) / high part of the residual
X = np.fft.rfft(xr, axis=0); f = np.fft.rfftfreq(len(xr), 1 / FS96)
for fc in (3, 12, 30, 100, 1000, 12000):
    lo = np.fft.irfft(X * (f < fc)[:, None], len(xr), axis=0)
    run(f"residual lowpassed < {fc} Hz", lo)
    run(f"residual high-passed > {fc} Hz", xr - lo)
# pre-roll only, and constant
for t in (100, 175, 260, 400, 1000, 5000):
    q = np.zeros_like(xr); q[:t] = xr[:t]; run(f"residual first {t} samples only", q)
dcm = xr.mean(0); run("constant = mean(residual) whole file", np.ones_like(xr) * dcm)
dcp = xr[:170].mean(0); run("constant = mean(first 170) whole file", np.ones_like(xr) * dcp)
print("mean resid whole", dcm, "first170", dcp)
# Nyquist tone amplitude through the file, per 0.1 s
alt = xr * (-1.0) ** np.arange(len(xr))[:, None]
print("alt (Nyquist) amplitude per 0.1 s ch0:", [float(alt[i:i+9600, 0].mean()) for i in range(0, len(xr) - 9600, 9600)])
print("DC-ish: block means per 0.1s ch0:", [float(xr[i:i+9600, 0].mean()) for i in range(0, len(xr) - 9600, 9600)])
