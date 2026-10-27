"""Per-capture clock ratio fitted in a HIGH band (A) and scored in other bands (B): is the 11-12.8 kHz residual a clock-ratio
precision problem? EQ = complex pooled EQ from the train half (heldout.json); score on the test half only.
Outputs analysis/nl_hfdrift/ppmfit.json"""
import sys, json, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
from scipy.signal import butter, sosfiltfilt
from scipy.optimize import minimize_scalar
M.EQ = False; M._proc.cache_clear()
H = json.load(open(OUT + "/heldout.json")); knots = np.array(H["knots"])
C = np.array([complex(a, b) for a, b in H["c_new_tr"]])
keys = [tuple(k) for k in H["keys"]]
test = [k for i, k in enumerate(keys) if i % 2 == 1]
if len(sys.argv) > 1: test = test[:int(sys.argv[1])]
SOS = {nm: butter(4, [lo, hi], btype="band", fs=FS, output="sos") for nm, (lo, hi) in
       dict(A="8000-11000".split("-") and (8000, 11000), B=(11000, 12800), C=(10000, 12000), LO=(1000, 5000), MID=(5000, 10000), HI1=(9000, 11000), HI2=(11000, 12000), HI3=(11800, 12800)).items()}
def nullb(y, m, nm):
    yb, mb = sosfiltfilt(SOS[nm], y, axis=0), sosfiltfilt(SOS[nm], m, axis=0)
    return np.array([10 * np.log10(np.sum((yb[:, c] - mb[:, c]) ** 2) / np.sum(yb[:, c] ** 2)) for c in (0, 1)])
rows = []
for i, key in enumerate(test):
    b = bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
    a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if b2 - a < int(0.05 * FS): continue
    n = b2 + 200
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    m = M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"])
    ms = m[a:b2]; L = len(ms); f = np.fft.rfftfreq(L) * FS
    E = np.interp(f, knots, C.real) + 1j * np.interp(f, knots, C.imag)
    me = np.fft.irfft(np.fft.rfft(ms, axis=0) * E[:, None], L, axis=0)
    cache = {}
    def Y(ppm):
        k = round(ppm, 4)
        if k not in cache: cache[k] = unwarp2(rawp, k)[a:b2]
        return cache[k]
    p0 = b["ppm"]
    out = dict(key=key, ppm_bench=p0)
    out["bench_ppm"] = {nm: nullb(Y(p0), me, nm).tolist() for nm in ("LO", "MID", "A", "B", "C", "HI1", "HI2", "HI3")}
    for fitband in ("A", "HI1", "B"):
        r = minimize_scalar(lambda p: float(np.mean(nullb(Y(p), me, fitband))), bounds=(p0 - 3, p0 + 3), method="bounded", options={"xatol": 0.01})
        out[f"ppm_fit_{fitband}"] = float(r.x)
        out[f"fit_{fitband}"] = {nm: nullb(Y(r.x), me, nm).tolist() for nm in ("LO", "MID", "A", "B", "C", "HI1", "HI2", "HI3")}
    rows.append(out)
    print(i, key, "bench %.2f  A %.2f  HI1 %.2f  B %.2f" % (p0, out["ppm_fit_A"], out["ppm_fit_HI1"], out["ppm_fit_B"]), flush=True)
json.dump(rows, open(OUT + "/ppmfit.json", "w"), default=float)
def med(sel, nm): return np.round(np.median(np.array([r[sel][nm] for r in rows]), 0), 1)
print("median over", len(rows), "test captures (L, R)")
for nm in ("LO", "MID", "A", "HI1", "HI2", "HI3", "B", "C"):
    print(f"{nm:4s} bench-ppm {med('bench_ppm', nm)}  fitA {med('fit_A', nm)}  fitHI1(9-11k) {med('fit_HI1', nm)}  fitB(11-12.8k, in-sample) {med('fit_B', nm)}")
d = np.array([[r["ppm_fit_A"] - r["ppm_bench"], r["ppm_fit_HI1"] - r["ppm_bench"], r["ppm_fit_B"] - r["ppm_bench"]] for r in rows])
print("ppm shift fitA / HI1 / B vs bench: mean", d.mean(0).round(3), "sd", d.std(0).round(3))
print("corr(fitA shift, fitB shift) = %.2f" % np.corrcoef(d[:, 0], d[:, 2])[0, 1])
