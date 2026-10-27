"""Is the 20-100 Hz (and 100-400 Hz) error in the 0.05-0.45 s window the band-pass response to the file-start step (pre-roll baseline offset, cache index < 46)?
For each capture and channel: s = mean(e[8:41]); e' = e - s * 1[idx<46]; band-pass and score in the window. Also the same with a free step per capture fitted in-band (LS on the filtered box column)."""
from cr_common import *
B = bench_all()
EAR = (0.7, 1.0, 1.5, 2.0, 3.0, 10.0)
sel = [k for k in sorted(B) if k[1] in EAR and not B[k]["trimmed"] and (k[0] == 0.5 or (k[0] in (0.2, 0.3, 0.4) and k[2] == 0.0))]
BD = [(20, 100), (100, 400), (400, 1000), (100, 1000), (1000, 5000)]
SOS = {b: band_sos(*b) for b in BD}
def db(x): return 10 * np.log10(x)
rows = []
for k in sel:
    r = B[k]; y, m = cached(k, r); a, b = window(k, r)
    e = y - m
    box = np.zeros(NCACHE); box[:46] = 1.0
    row = {}
    for bd in BD:
        out = []
        for c in (0, 1):
            s = e[8:41, c].mean()
            yb = sosfiltfilt(SOS[bd], y[:, c])[a:b]; eb = sosfiltfilt(SOS[bd], e[:, c])[a:b]
            e2 = sosfiltfilt(SOS[bd], e[:, c] - s * box)[a:b]
            xb = sosfiltfilt(SOS[bd], box)[a:b]; beta = (xb @ eb) / (xb @ xb); e3 = eb - beta * xb
            out.append((db(fpow(eb) / fpow(yb)), db(fpow(e2) / fpow(yb)), db(fpow(e3) / fpow(yb)), s, beta))
        row[bd] = out
    rows.append((k, row, e[8:41, :].mean(0), e[:, :][60:, :].mean(0)))
print("median over captures x channels: null (bench alignment) | step removed from the pre-roll mean | step amplitude LS-fitted in band")
for bd in BD:
    v = np.array([[x[:3] for x in row[bd]] for _, row, _, _ in rows])           # (cap, ch, 3)
    print(f"  {bd[0]:>5}-{bd[1]:<5}: {np.median(v[:,:,0]):6.1f} | {np.median(v[:,:,1]):6.1f} | {np.median(v[:,:,2]):6.1f}")
print("by EARLY (20-100 and 100-1000):")
for e_ in EAR:
    for bd in ((20, 100), (100, 1000)):
        v = np.array([[x[:3] for x in row[bd]] for k, row, _, _ in rows if k[1] == e_ and k[0] == 0.5])
        print(f"  EARLY {e_:5.1f} {bd}: {np.median(v[:,:,0]):6.1f} | {np.median(v[:,:,1]):6.1f} | {np.median(v[:,:,2]):6.1f}")
# step amplitude vs LATE in SPAN 0.5, EARLY 2 (units of full-band model peak)
print("pre-roll mean of e (units of 1e-4 x peak of |m|) for SPAN 0.5 EARLY 2 by LATE, L / R:")
for k, row, pm, _ in rows:
    if k[0] == 0.5 and k[1] == 2.0:
        pk = np.abs(cached(k, B[k])[1]).max(0)
        print("  ", k[2], np.round(1e4 * pm / pk, 2))
