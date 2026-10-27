"""If every main tap carries an independent relative gain error d (rms), the model error / signal ratio at frequency f is
d * sqrt(sum g_i^2) / |sum s_i g_i exp(-j w t_i)|.  Compare the amplification across bands with the observed band errors."""
from lb_common import *
for early in (10.0, 1.0, 0.5):
  for ch in (0, 1):
    ns = ne.N_STAGES[0.5]
    fdb, sdb = ne.law_db(0.5, early, ch)
    tl = ne.taps(ch, ns)
    g = np.array([s * 10 ** (sdb[[j for j, st in enumerate(ne.STAGES) if i in st][0]] / 20) for t, s, i in tl]); t = np.array([t for t, s, i in tl], float)
    # window 0.05-0.45 s: taps with t in that range (device samples) -- the window sees taps up to 0.45 s
    sel = (t > 0.05 * FS * 0) & (t < 0.45 * FS)
    out = []
    for lo, hi in ((20, 100), (100, 1000), (1000, 5000), (5000, 10000)):
        f = np.geomspace(lo, hi, 400)
        H = np.array([abs(np.sum(g[sel] * np.exp(-2j * np.pi * ff / FS * t[sel]))) for ff in f])
        amp = np.sqrt(np.sum(g[sel] ** 2)) / np.sqrt(np.mean(H ** 2))
        out.append(20 * np.log10(amp))
    print(f"EARLY {early:5.1f} ch{ch}: amplification of a per-tap relative gain error (dB re the error itself): 20-100 {out[0]:5.1f} | 100-1k {out[1]:5.1f} | 1-5k {out[2]:5.1f} | 5-10k {out[3]:5.1f}   (n taps {sel.sum()})")
print("observed model error (SPAN 0.5, fixed bench, EARLY 1-2): 20-100 -39.5 | 100-1k -57 | 1-5k -61 | 5-10k -58 dB;  repeat floor -58..-62 | -67..-69 | -60 | -59")
