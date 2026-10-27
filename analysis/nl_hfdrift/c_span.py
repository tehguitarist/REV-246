"""Fork C: is the residual per capture SESSION (SPAN) in both the front taps and the reverb? Per SPAN, R channel (clean three-tap-free ch1) and L up to 12.25 kHz."""
import os, sys, random, pickle
os.environ["NL_WARP2"] = "1"; os.environ.pop("NL_BK", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter, sosfiltfilt
import c_front as CF
import defects
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
front = pickle.load(open(OUT + "/c_front_150.pkl", "rb"))["all"]
fkeys = {tuple(r["key"]) for r in front}
ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and k not in fkeys and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
random.seed(29); sample = random.sample(ks, 150)
sel = [0, 3, 5, 6, 7, 8, 9]     # knots 10, 11.25, 11.75, 12, 12.25, 12.5, 12.8 kHz
rv = {}
for key in sample:
    b = common.bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
    a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if b2 - a < int(0.05 * FS): continue
    n = b2 + 200; rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = nlt.unwarp(rawp, b["ppm"])[a:b2]; m = (M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"]))[a:b2]
    L = len(m); w = np.hanning(L); f = np.fft.rfftfreq(L) * FS; W = CF.hats(f)
    Y = np.fft.rfft(y * w[:, None], axis=0); Mm = np.fft.rfft(m * w[:, None], axis=0)
    rv[key] = [((W * (Y[:, c] * np.conj(Mm[:, c]))[None, :]).sum(1), (W * (np.abs(Mm[:, c]) ** 2)[None, :]).sum(1)) for c in (0, 1)] + [b["ppm"]]
spans = sorted({k[0] for k in rv})
print("phase mrad at knots", [CF.KNOTS[i] / 1000 for i in sel], "\n(front taps | reverb), per SPAN, n front / n reverb, median ppm of the reverb set")
for ch in (0, 1):
    for sp in spans:
        fr_rows = [r for r in front if r["key"][0] == sp]; rv_keys = [k for k in rv if k[0] == sp]
        if not fr_rows or not rv_keys: continue
        cf = sum(r[ch]["num"] for r in fr_rows) / sum(r[ch]["den"] for r in fr_rows)
        cr = sum(rv[k][ch][0] for k in rv_keys) / sum(rv[k][ch][1] for k in rv_keys)
        print(f"ch{ch} SPAN {sp:.1f} n {len(fr_rows):2d}/{len(rv_keys):2d} ppm {np.median([rv[k][2] for k in rv_keys]):5.2f} | front " + " ".join("%+6.1f" % (np.angle(cf[i]) * 1e3) for i in sel) + " | reverb " + " ".join("%+6.1f" % (np.angle(cr[i]) * 1e3) for i in sel))
