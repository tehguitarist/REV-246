"""Fork A part 1+3a: repeat-capture ceiling per sub-band (NEW default warp) and the model's null on the SAME captures/windows.
Pairs = captures with the same SPAN and EARLY, different LATE (the early section is identical until the late section arrives: SPAN 0.5 window
0.05-0.45 s, SPAN 0.4 window 0.05-0.30 s). N has no true same-setting repeats (one file per setting).
Alignment fitted on 8-9 kHz ONLY (held out: 9-11, 11-11.6, 11.6-12.3, 12.3-12.8, 10-12 kHz):
  V1 delay + gain;  V2 delay + gain + relative clock (ppm);  V3 = V2 plus one complex gain per scored band (pure repeatability).
Model null = bench fits (o, ppm, gain), exact resampler gain, branch kernel on (bench_warp2_bk_kf_a0.7), 'exact gain only', same interior window.
usage: a_repeat.py [n_groups]"""
import os, sys, json, itertools
os.environ["NL_WARP2"] = "1"; os.environ["NL_BK"] = "1"
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter, sosfiltfilt, hilbert
from scipy.optimize import minimize
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2_bk_kf_a0.7"
assert nlt.WARP2 and M.BK is not None
PAD = 200
BANDS = dict(fit=(8000, 9000), b9_11=(9000, 11000), b11_11p6=(11000, 11600), b11p6_12p3=(11600, 12300), b12p3_12p8=(12300, 12800), b10_12=(10000, 12000))
SOS = {k: butter(4, v, btype="band", fs=FS, output="sos") for k, v in BANDS.items()}
def bp(x, k): return sosfiltfilt(SOS[k], x, axis=0)
def seg(key, a_ref, tend):
    b = common.bench(key); raw = nl.ir(nl.table()[key]); fr = b["front_shift"]
    a = M.K0 + int(a_ref) + int(0.05 * FS); e = M.K0 + int(a_ref) + int(tend * FS)
    rawp = np.vstack([raw[:e + PAD + 80], np.zeros((max(0, e + PAD + 80 - len(raw)), 2))])
    up = nlt.upsample(rawp)
    return dict(key=key, b=b, up=up, a=a - PAD, e=e + PAD, n=e + PAD + 1)
def warped(c, dppm, tau):
    """segment [a-PAD, e+PAD) of the capture, extra clock dppm (stretch about the segment start) and delay tau (samples), (n, 2)"""
    full = nlt.unwarp_from(c["up"], c["n"], c["b"]["ppm"] + dppm)[c["a"]:c["e"]]
    if tau == 0: return full
    f = np.fft.rfftfreq(len(full))
    return np.fft.irfft(np.fft.rfft(full, axis=0) * np.exp(-2j * np.pi * f * tau)[:, None], len(full), axis=0)
def nulls(y1, y2, bands, cg=False):
    """per band, per channel null (dB) of y1 - y2 (interior), y2 optionally with one complex gain per band and channel"""
    out = {}
    for k in bands:
        a, b = bp(y1, k)[PAD:-PAD], bp(y2, k)[PAD:-PAD]
        r = []
        for c in range(2):
            if cg:
                h = np.imag(hilbert(b[:, c])); X = np.vstack([b[:, c], h]).T
                co = np.linalg.lstsq(X, a[:, c], rcond=None)[0]; d = a[:, c] - X @ co
            else:
                d = a[:, c] - b[:, c]
            r.append(10 * np.log10(np.sum(d ** 2) / np.sum(a[:, c] ** 2)))
        out[k] = r
    return out
def model_seg(c):
    key = c["key"]; b = c["b"]; fr = b["front_shift"]
    n = c["e"] + 5
    m = M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"])
    return m[c["a"]:c["e"]]
if __name__ == "__main__":
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 999
    groups = {}
    for k in nl.table():
        if k[0] in (0.5, 0.4): groups.setdefault((k[0], k[1]), []).append(k)
    groups = {g: sorted(v, key=lambda k: k[2]) for g, v in groups.items() if len(v) >= 4}
    glist = sorted(groups)
    glist = glist[:: max(1, len(glist) // ng)][:ng]
    res = []
    for g in glist:
        keys = groups[g]; tend = 0.45 if g[0] == 0.5 else 0.30
        a_ref = float(np.median([common.bench(k)["front_shift"] for k in keys]))
        cs = [seg(k, a_ref, tend) for k in keys]
        ms = [model_seg(c) for c in cs]
        ys = [warped(c, 0.0, 0.0) for c in cs]
        pairs = sorted({(i, j) for i in range(len(cs)) for j in (i + 1, i + 3) if j < len(cs)})
        gres = dict(group=g, n=len(cs), pairs=[], model=[])
        for (i, j) in pairs:
            y1 = ys[i]
            def fitcost(v, use_ppm):
                y2 = warped(cs[j], v[2] if use_ppm else 0.0, v[0])
                a, b = bp(y1, "fit")[PAD:-PAD], bp(y2, "fit")[PAD:-PAD]
                gn = np.sum(a * b) / np.sum(b * b)
                return float(np.log(np.sum((a - gn * b) ** 2) / np.sum(a ** 2)))
            # V1: delay by scan then refine; V2: + relative clock
            sc = [(fitcost([t, 0, 0], False), t) for t in np.arange(-1.0, 1.001, 0.05)]
            t0 = min(sc)[1]
            r1 = minimize(lambda v: fitcost([v[0], 0, 0], False), [t0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-8))
            r2 = minimize(lambda v: fitcost([v[0], 0, v[1]], True), [r1.x[0], 0.0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-8, initial_simplex=[[r1.x[0], 0], [r1.x[0] + 0.02, 0], [r1.x[0], 0.3]]))
            out = {}
            for nm, (tau, dp) in (("V1", (r1.x[0], 0.0)), ("V2", (r2.x[0], r2.x[1]))):
                y2 = warped(cs[j], dp, tau)
                a, b = bp(y1, "fit")[PAD:-PAD], bp(y2, "fit")[PAD:-PAD]
                gn = np.sum(a * b) / np.sum(b * b)
                out[nm] = nulls(y1, gn * y2, list(BANDS))
                if nm == "V2": out["V3"] = nulls(y1, y2, list(BANDS), cg=True); out["tau_ppm"] = [float(tau), float(dp)]
            gres["pairs"].append(dict(i=keys[i], j=keys[j], **out))
        for c, m, y in zip(cs, ms, ys):
            gres["model"].append(dict(key=c["key"], **{"model": nulls(y, m, list(BANDS)), "model_cg": nulls(y, m, list(BANDS), cg=True)}))
        res.append(gres); print(g, "pairs", len(pairs), flush=True)
    json.dump(res, open(OUT + "/a_repeat.json", "w"), default=float)
    def med(sel, k, ch=None):
        v = np.array([p[sel][k] for gr in res for p in gr["pairs"]]) if sel in ("V1", "V2", "V3") else np.array([m[sel][k] for gr in res for m in gr["model"]])
        return np.round(np.median(v, 0), 1)
    print("median over", sum(len(g["pairs"]) for g in res), "pairs,", sum(len(g["model"]) for g in res), "captures (L, R)")
    for k in BANDS:
        print(f"{k:11s} repeat V1 {med('V1', k)} V2 {med('V2', k)} V3(cplx gain) {med('V3', k)} | model {med('model', k)} model+cplx gain {med('model_cg', k)}")
