"""Normal-equation builder: time-domain LS of the capture residual on frequency-hat columns applied to the FULL model record (no window leak).
r = y - m_base (window a..b);  columns B_k m_base, B_k = hat_k(f) (re) and 1j*hat_k(f) (im).  m_base = model through NL_LOWBAND=BASE (0: raw, 2: shipped)."""
from common import *
import time
BASE = int(os.environ.get("SUB_BASE", "0"))
KNOTS = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 34, 40], float)    # hat knots (Hz), 17 -> 34 real columns
NC = 2 * len(KNOTS)
SQ = 1 << 20

def hats(f):
    H = np.zeros((len(KNOTS), len(f)))
    for i, k in enumerate(KNOTS):
        lo = KNOTS[i - 1] if i else 0.0; hi = KNOTS[i + 1] if i + 1 < len(KNOTS) else k + 8.0
        H[i] = np.interp(f, [lo, k, hi], [0, 1, 0], left=0, right=0)
    return H

def window_variants(key, r, end, raw_len):
    """name -> (a, b) sample indices."""
    z = K0 + r["front_shift"]; a = int(z + 0.05 * FS)
    v = {"w45": (a, int(min(z + 0.45 * FS, end)))}
    b2 = z + {0.0: 0.1, 0.1: 0.2}.get(key[0], key[0] + 0.15) * FS
    if key[2] > 0 and FEED.get(key[0]):
        v["pre_feed"] = (a, int(min(z + FEED[key[0]] + r["o"] - 0.01 * FS, end)))
    else:
        v["pre_feed"] = v["w45"]
    v["w60"] = (a, int(min(z + 0.60 * FS, end)))
    v["to_end"] = (a, int(end))                      # to the end of usable data
    v["span"] = (a, int(min(b2, end)))               # SPAN + 0.15 s (LATE 0: all the signal there is)
    return {k: ab for k, ab in v.items() if ab[1] - ab[0] > 0.1 * FS}

def build(key, r, variants=("w45", "pre_feed", "w60", "to_end", "span")):
    """Returns dict variant -> {ch: (G (NC,NC), b (NC,), yy, n, scale)} and the model energy; also needs nothing stored but sums."""
    raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h)
    wv = {k: ab for k, ab in window_variants(key, r, end, len(raw)).items() if k in variants}
    n = max(b for a, b in wv.values()) + 64
    rawp = np.vstack([raw[:n], np.zeros((max(0, n - len(raw)), 2))])
    y = nlt.unwarp(rawp, r["ppm"])[:n]
    N = 1 << int(np.ceil(np.log2(n + 2 * ne.TAPS[-1] + 4096)))
    f = np.fft.rfftfreq(N, 1 / FS); ph = -2j * np.pi * np.fft.rfftfreq(N)
    M._proc.cache_clear(); os.environ["NL_LOWBAND"] = str(BASE); R = M._proc(N)
    Hh = hats(f)
    out = {v: {} for v in wv}
    for ch in (0, 1):
        F, Br = M.spectra(key, ch, N)
        S = (F + Br * np.exp(ph * r["o"])) * np.exp(ph * (K0 + r["front_shift"])) * R[:, ch] * r["gain"][ch]
        m = np.fft.irfft(S, N)[:n]
        cols = np.empty((n, NC))
        for i in range(len(KNOTS)):
            cols[:, 2 * i] = np.fft.irfft(S * Hh[i], N)[:n]
            cols[:, 2 * i + 1] = np.fft.irfft(S * 1j * Hh[i], N)[:n]
        for v, (a, b) in wv.items():
            sc = np.sqrt(np.mean(m[a:b] ** 2)) + 1e-30
            res = (y[a:b, ch] - m[a:b]) / sc; Cc = cols[a:b] / sc
            out[v][ch] = (Cc.T @ Cc, Cc.T @ res, float(res @ res), b - a, sc)
    return out

if __name__ == "__main__":
    B, ks = catalog(0.3, 0.5, 0.4, 0.3)
    import random; random.seed(1); sel = random.sample(ks, 3)
    for k in sel:
        t0 = time.time(); o = build(k, B[k]); print(k, sess(B[k]), "time %.1f s" % (time.time() - t0), list(o), flush=True)
