import sys, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
key = (0.0, 10.0, 2.0)
y, m, info = pair(key)
ppm = 12.0
print("kernel gain default at 11,11.5,12,12.3,12.6,12.8,13.2 kHz:", np.round(20*np.log10(unwarp_kernel_response(f_hz=[11000,11500,12000,12300,12600,12800,13190])), 2), "dB")
for name, fn in (("default", nlt.unwarp), ("long", unwarp2)):
    z = fn(np.vstack([m, np.zeros((64, 2))]), ppm)[:len(m)]
    print(name)
    for lo, hi in [(6000, 8000), (11000, 11600), (11600, 12300), (12300, 12800)]:
        r = local_delay(z[:, 0], m[:, 0], lo, hi, win=4096, hop=4096, start=info["i0"])
        sel = r[r[:, 4] > r[:, 4].max() - 40][:6]
        print("  %5d-%5d " % (lo, hi) + " ".join("%.2f:%+.3f" % (t, tp) for t, tp, tg, c, e in sel))
    # amplitude ratio at the band edge, 0.4 s
    a = int(0.4 * FS); w = np.hanning(4096); f = np.fft.rfftfreq(4096) * FS
    Z = np.fft.rfft(z[a:a + 4096, 0] * w); Mm = np.fft.rfft(m[a:a + 4096, 0] * w)
    print("  |Z/M| at 11.6/12.1/12.6k:", [round(float(np.sqrt(np.sum(abs(Z[(f>=c-125)&(f<c+125)])**2)/np.sum(abs(Mm[(f>=c-125)&(f<c+125)])**2))),3) for c in (11625,12125,12625)])
