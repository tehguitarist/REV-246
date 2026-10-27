"""build2: as fitlib.build plus the 20-100 Hz bench band (4th-order band-pass, sosfiltfilt, window 0.05-0.6 s) in normal-equation form, raw and with the pre-roll step subtracted,
so that ANY correction c (held out any way) can be scored exactly:  null(c) = ee - 2 c.bb + c.Gb.c  over yyb."""
from fitlib import *
BAND = [("20-100", bands(20, 100)), ("100-300", bands(100, 300))]
T_ = np.arange(4096); STEP = 0.5 * (1 - np.tanh((T_ - 45.5) / 1.5))

def build2(key, r, variants=("w45", "pre_feed", "w60", "to_end", "span")):
    raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h)
    wv = {k: ab for k, ab in window_variants(key, r, end, len(raw)).items() if k in variants}
    zi = K0 + int(r["front_shift"]); ba = zi + int(0.05 * FS); bb_ = int(min(end, zi + 0.6 * FS))
    n = max(max(b for a, b in wv.values()), bb_) + 64
    rawp = np.vstack([raw[:n], np.zeros((max(0, n - len(raw)), 2))])
    y = nlt.unwarp(rawp, r["ppm"])[:n]
    N = 1 << int(np.ceil(np.log2(n + 2 * ne.TAPS[-1] + 4096)))
    f = np.fft.rfftfreq(N, 1 / FS); ph = -2j * np.pi * np.fft.rfftfreq(N)
    M._proc.cache_clear(); os.environ["NL_LOWBAND"] = str(BASE); R = M._proc(N)
    Hh = hats(f)
    out = {"ls": {v: {} for v in wv}, "band": {}, "c": []}
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
            out["ls"][v][ch] = (Cc.T @ Cc, Cc.T @ res, float(res @ res), b - a, sc)
        e = y[:, ch] - m; cc = float(e[8:40].mean()); out["c"].append(cc / np.abs(m[zi:zi + 800]).max())
        for nm, sos in BAND:
            P = lambda x: sosfiltfilt(sos, x, axis=0)
            yb = P(y[:, ch])[ba:bb_]; Pc = P(cols)[ba:bb_]
            Pm = P(m)[ba:bb_]; Pu = P(STEP[:n] if n <= len(STEP) else np.r_[STEP, np.zeros(n - len(STEP))])[ba:bb_]
            er = yb - Pm; ec = er - cc * Pu
            out["band"][(nm, ch)] = dict(yy=float(yb @ yb), ee_raw=float(er @ er), ee_cor=float(ec @ ec), Gb=Pc.T @ Pc,
                                         bb_raw=Pc.T @ er, bb_cor=Pc.T @ ec, n=bb_ - ba)
    return out
