"""Within one session: does the first-25 ms residual scale with the front-tap level (a transfer / processing effect on the content) or stay fixed
(additive artefact)? e (20-300 Hz, front frame, absolute model units) of every EARLY value vs the EARLY-10 e of the same session."""
import sys
from lb_common import *
B = bench_all()
sos = bands(20, 300)
ia, ib = K0 - 40, K0 + int(0.025 * FS)
for span, lates in ((0.3, (0.0, 0.2, 0.4)), (0.5, (0.0, 0.2, 0.4))):
    sel = {}
    for k, r in sorted(B.items()):
        if k[0] != span or k[2] not in lates or abs(r["front_shift"] + 0.13) > 0.2 or abs(r["o"] - 0.5) > 1.2: continue
        y, m, _ = load(k, r, n=3000)
        yf, mf = frame(y, r["front_shift"]), frame(m, r["front_shift"])
        sel.setdefault(k[1], []).append((k, yf, mf))
    if 10.0 not in sel: print("SPAN", span, "no EARLY 10 in fs0 / o~1 set; EARLY present", sorted(sel)); 
    ref_early = 10.0 if 10.0 in sel else (3.0 if 3.0 in sel else None)
    if ref_early is None: continue
    ref = sel[ref_early][0]
    print(f"\nSPAN {span} session fs0: reference {ref[0]}")
    for c in (0, 1):
        er = sosfiltfilt(sos, (ref[1] - ref[2])[:, c])[ia:ib]
        print(f"  ch{c}: EARLY  n  | front-tap peak (model, abs units) | e rms (abs) | corr with ref e | LS gain e/e_ref  | gain of 410-tap vs ref")
        for early in sorted(sel, reverse=True):
            row = []
            for k, yf, mf in sel[early]:
                e = sosfiltfilt(sos, (yf - mf)[:, c])[ia:ib]
                pk = np.abs(mf[K0 + 300:K0 + 500, c]).max()
                row.append((pk, np.sqrt(np.mean(e ** 2)), np.corrcoef(e, er)[0, 1], (e @ er) / (er @ er)))
            row = np.array(row).mean(0)
            pk_ref = np.abs(ref[2][K0 + 300:K0 + 500, c]).max()
            print(f"        {early:5.1f} {len(sel[early]):2d} | {20*np.log10(row[0]):6.1f} dB | {20*np.log10(row[1]):6.1f} dB | {row[2]:6.2f} | {row[3]:6.2f} | {row[0]/pk_ref:6.3f}")
