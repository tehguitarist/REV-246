"""Bench-style scoring (tools/nl_model.bench_file logic, re-implemented) of the LTI model against the capture, with and without the
file-DC step response added to the model. DC = mean of the file's first N0 samples (before any arrival), per channel: measured, not fitted."""
import sys, json, time; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
from scipy.signal import butter, sosfiltfilt
np.set_printoptions(linewidth=220, precision=2, suppress=True)
N0 = int(os.environ.get("N0", "30"))
B20, B100 = (butter(4, [lo, hi], btype="band", fs=FSD, output="sos") for lo, hi in ((20, 100), (100, 1000)))
WIN = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.6)]
def nul(y, m, a, b, sos=None):
    return np.array([10 * np.log10(np.sum((y[a:b, c] - m[a:b, c]) ** 2) / np.sum(y[a:b, c] ** 2)) for c in (0, 1)])
keys = [tuple(float(v) for v in a.split(",")) for a in sys.argv[2:]]
res = {}
for key in keys:
    t0 = time.time()
    r = align(key); ppm = r["ppm"]; fr = r["front_shift"]
    x, fs_ = sf.read(file_path(key), always_2d=True); nf = len(x)
    raw = nl.ir(key); h = sosfiltfilt(ne.BAND, raw, axis=0)
    end = M.usable_end(h); n = min(end, K0 + int(fr) + int(2.0 * FSD)); i0, trimmed = M.window_start(h, fr)
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    yf = nlt.unwarp(rawp, ppm)[:n]
    mf = M.render(key, n, r["o"], fr, band=False) * np.array(r["gain"])
    dc = x[:N0].mean(0)
    hdc, _ = pipeline(np.ones((nf, 2)) * dc)
    hp = np.vstack([hdc[:n + 64], np.zeros((max(0, n + 64 - len(hdc)), 2))])
    corr = nlt.unwarp(hp, ppm)[:n]
    out = {}
    tt = (np.arange(n) - (K0 + int(fr))) / FSD
    wt = np.clip((0.08 - tt) / 0.04, 0, 1); wt = 0.5 - 0.5 * np.cos(np.pi * wt)            # 1 until 40 ms after the front frame, 0 from 80 ms
    for tag, mm in (("base", mf), ("dcfix", mf + corr), ("dctaper", mf + corr * wt[:, None])):
        y, m = sosfiltfilt(ne.BAND, yf, axis=0), sosfiltfilt(ne.BAND, mm, axis=0)
        w = {}
        for t0_, t1_ in WIN:
            a, b = max(i0, K0 + int(fr) + int(t0_ * FSD)), min(n, K0 + int(fr) + int(t1_ * FSD))
            if b - a >= int(0.02 * FSD): w[f"{t0_}-{t1_}"] = nul(y, m, a, b).tolist()
        for nm, sos in (("20-100", B20), ("100-1000", B100)):
            yb, mb = sosfiltfilt(sos, yf, axis=0), sosfiltfilt(sos, mm, axis=0)
            for t0_, t1_ in ((0.0, 0.05), (0.05, 0.6)):
                a, b = max(i0, K0 + int(fr) + int(t0_ * FSD)), min(n, K0 + int(fr) + int(t1_ * FSD))
                if b - a > int(0.03 * FSD): w[f"{nm}@{t0_}-{t1_}"] = nul(yb, mb, a, b).tolist()
        out[tag] = w
    bj = json.load(open(os.path.join(BENCH, "_".join(map(str, key)) + ".json")))
    out["bench_json_0.0-0.05"] = [bj["windows"].get("0.0-0.05", {}).get("null0"), bj["windows"].get("0.0-0.05", {}).get("null1")]
    out["bench_json_20-100"] = bj["bands"]["20-100"]
    out["dc"] = dc.tolist(); out["trimmed"] = bool(trimmed); out["nf_s"] = nf / FS96
    res["_".join(map(str, key))] = out
    print(key, f"trimmed {trimmed} {time.time()-t0:.0f}s dc {dc*1e6}")
    for k in out["base"]:
        print(f"   {k:16s} base {np.array(out['base'][k])}  dcfix {np.array(out['dcfix'][k])}  delta {np.array(out['dcfix'][k]) - np.array(out['base'][k])}")
    print("   bench json 0.0-0.05", out["bench_json_0.0-0.05"], " 20-100", out["bench_json_20-100"], flush=True)
json.dump(res, open(sys.argv[1], "w"))
