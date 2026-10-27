"""(b) Is the deterministic error concentrated at the taps? Undo the diffusers on the LATE-averaged residual (all-pass inverse = time-reversed all-pass), then compare the
error energy in +-24 samples around each tap (and R's) with the error between taps, per tap relative to the tap's amplitude."""
from cr_common import *
B = bench_all()
BD = [(20, 100), (100, 1000), (1000, 5000)]
SOS = {b: band_sos(*b) for b in BD}
def db(x): return 10 * np.log10(x)
def undiff(x, ch):
    z = x[::-1].copy()
    for Mm in (ne.AP_IN, ne.AP_OUT[ch]): z = ne.ap(z, Mm)
    return z[::-1]
H = 24
for (e_, fsel) in ((2.0, 2), (3.0, 0), (10.0, 0), (1.0, 2)):
    ks = [k for k in sorted(B) if k[0] == 0.5 and k[1] == e_ and not B[k]["trimmed"] and round(B[k]["front_shift"] / 10) == fsel]
    ref = ks[0]; ks = [k for k in ks if abs(B[k]["o"] - B[ref]["o"]) < 0.02]
    r = B[ref]; ns = 14; fdb, sdb = ne.law_db(0.5, e_, 0)
    print(f"\nSPAN 0.5 EARLY {e_} session {fsel}: {len(ks)} captures averaged, o {r['o']:.3f}")
    for bd in BD:
        for ch in (0, 1):
            ea = np.mean([sosfiltfilt(SOS[bd], (lambda y, m: y - m)(*cached(k, B[k])), axis=0)[:, ch] for k in ks], axis=0)
            ya = np.mean([sosfiltfilt(SOS[bd], cached(k, B[k])[0], axis=0)[:, ch] for k in ks], axis=0)
            edd = undiff(ea, ch); ydd = undiff(ya, ch)
            z = K0 + r["front_shift"] + r["o"]
            tp = [(t, i) for t, s, i in ne.taps(ch, ns)]
            ntap = np.zeros(len(edd), bool); idxs = []
            rows = []
            tot_e = 0
            for t, i in tp:
                c = int(round(z + t)) 
                lo, hi = c - H, c + H + 1
                if hi > len(edd): continue
                e_loc = fpow(edd[lo:hi]); y_loc = fpow(ydd[lo:hi]); ntap[lo:hi] = True
                j = [jj for jj, st in enumerate(ne.STAGES) if i in st][0]
                rows.append((i, c, db(e_loc / y_loc)))
            a, b = window(ref, r)
            win = np.zeros(len(edd), bool); win[a:b] = True
            inside = fpow(edd[win & ntap]); outside = fpow(edd[win & ~ntap]); yin = fpow(ydd[win & ntap]); yout = fpow(ydd[win & ~ntap])
            print(f"  {bd[0]:>4}-{bd[1]:<4} {'LR'[ch]}: in-window de-diffused error: at taps (+-{H}) {db(inside/yin):6.1f} dB re signal there; between taps {db(outside/max(yout,1e-30)):6.1f} dB re signal there "
                  f"(signal outside taps is {db(yout/yin):5.1f} dB re at taps; error outside/inside {db(outside/inside):5.1f} dB, tap share of samples {np.mean(ntap[win]):.2f})")
