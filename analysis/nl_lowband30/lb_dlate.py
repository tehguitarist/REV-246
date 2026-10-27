"""Captures that differ only in LATE share the same device output for t < FEED (causality): any difference in the first 25 ms is a
capture/processing artefact. Measure dY = y(LATE) - y(LATE 0.0) on [file+0, z+25ms] and compare to the model's dM."""
import sys
from lb_common import *
B = bench_all()
sosA, sosB = bands(20, 100), bands(100, 1000)
for sp, ea in ((0.4, 10.0), (0.5, 10.0), (0.3, 10.0), (0.3, 1.0), (0.2, 10.0), (0.1, 10.0), (0.0, 10.0)):
    ref = (sp, ea, 0.0)
    if ref not in B: print("no ref", ref); continue
    n = 2400
    r0 = B[ref]; y0, m0, _ = load(ref, r0, n=n)
    z = int(round(K0 + r0["front_shift"]))
    print(f"\nSPAN {sp} EARLY {ea}: dY = y(LATE)-y(0) in [file+0, z+25 ms], rms dB re 410-peak; (dM = same for the model); fs0={r0['front_shift']:.2f}")
    print("  LATE   fs   | full-band dY  dM  | 20-100: dY  dM  | 100-1k: dY dM | mean dY[10:70]e-4 ")
    for late in (0.2, 0.4, 0.8, 1.2, 2.0):
        k = (sp, ea, late)
        if k not in B: continue
        r = B[k]
        if abs(r["front_shift"] - r0["front_shift"]) > 0.2: print("  session differs", late, r["front_shift"]); continue
        y1, m1, _ = load(k, r, n=n)
        pk = np.abs(m0[z:z + 800, 0]).max()
        ib = int(z + 0.025 * FS)
        q = lambda v: 20 * np.log10(np.sqrt(np.mean(v ** 2)) / pk + 1e-12)
        out = []
        for c in (0,):
            dy, dm = y1[:, c] - y0[:, c], m1[:, c] - m0[:, c]
            out.append(f"{q(dy[:ib]):6.1f} {q(dm[:ib]):6.1f} | {q(sosfiltfilt(sosA, dy)[:ib]):6.1f} {q(sosfiltfilt(sosA, dm)[:ib]):6.1f} | {q(sosfiltfilt(sosB, dy)[:ib]):6.1f} {q(sosfiltfilt(sosB, dm)[:ib]):6.1f} | {dy[10:70].mean() / pk * 1e4:6.2f} (R {(y1[10:70,1]-y0[10:70,1]).mean()/np.abs(m0[z:z+800,1]).max()*1e4:6.2f})")
        print(f"  {late:4.1f} {r['front_shift']:6.2f} | " + out[0])
