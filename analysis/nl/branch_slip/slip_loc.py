"""Where does the branch's shift happen? Per window null of the capture against the bench model at o=+1 and at the
bench's fitted o (same ppm, per-window gain). Windows: front taps, then 10 ms windows along the branch."""
import json, os, sys
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
import nl_model as nm, nl_early as ne, nl_late as nlt, nl_structure as nl
from scipy.signal import sosfiltfilt
FS, K0 = nl.FSD, nl.K0
STEP = 96000 / FS

def run(key, tmax=0.4, wlen=0.01):
    b = json.load(open(os.path.join(nm.BENCH, "_".join(str(v) for v in key) + ".json")))
    raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; o = b["o"]; ppm = b["ppm"]
    n = min(nm.usable_end(h), K0 + int(fr) + int(tmax * FS) + 10)
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = sosfiltfilt(ne.BAND, nlt.unwarp_from(nlt.upsample(rawp), n, ppm)[:n], axis=0)
    m1 = nm.render(key, n, 1.0, fr); mo = nm.render(key, n, o, fr)
    i0, _ = nm.window_start(h, fr)
    wins = [(max(i0, K0 + int(fr) + 20), K0 + int(fr) + 470)]
    t = 0.026
    while t < tmax:
        wins.append((K0 + int(fr) + int(t * FS), K0 + int(fr) + int((t + wlen) * FS))); t += wlen
    out = []
    for a, bb in wins:
        a = max(a, i0); bb = min(bb, n)
        if bb - a < 50: continue
        r = []
        for m in (m1, mo):
            c = 1  # R
            e = []
            for c in (0, 1):
                ym, mm = y[a:bb, c], m[a:bb, c]
                if mm @ mm == 0: e.append(np.nan); continue
                g = (ym @ mm) / (mm @ mm)
                e.append(10 * np.log10(np.sum((ym - g * mm) ** 2) / np.sum(ym ** 2)))
            r.append(np.nanmean(e))
        out.append(((a - K0 - fr) / FS * 1e3, r[0], r[1]))
    return o, (o - 1) * STEP, out

if __name__ == "__main__":
    for s in sys.argv[1:]:
        key = tuple(float(v) for v in s.split("_"))
        o, nh, out = run(key)
        print(f"{key} o={o:+.3f} n_host={nh:+.2f}")
        print("  " + " ".join(f"{t:.0f}:{a:.0f}/{b:.0f}" for t, a, b in out))
