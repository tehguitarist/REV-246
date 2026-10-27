"""Predict the drift without any capture: apply the harness's own unwarp to the MODEL at the capture's ppm and measure the
delay of the 11-12.8 kHz bands relative to 6-8 kHz. Compare with the measured capture-vs-model delay (old unwarp)."""
import sys, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
bands = [(11000, 11600), (11600, 12300), (12300, 12800)]
for key in [(0.0, 10.0, 2.0), (0.5, 10.0, 2.0), (0.3, 1.0, 1.0)]:
    y, m, info = pair(key); ppm = info["b"]["ppm"]; a = info["i0"]
    z = nlt.unwarp(np.vstack([m, np.zeros((64, 2))]), ppm)[:len(m)]
    print(key, "ppm %.2f" % ppm, " windows centred at", end=" ")
    ref = local_delay(z[:, 0], m[:, 0], 6000, 8000, win=4096, hop=4096, start=a)
    tt = ref[:5, 0]; print(np.round(tt, 2))
    for lo, hi in bands:
        pr = local_delay(z[:, 0], m[:, 0], lo, hi, win=4096, hop=4096, start=a); me = local_delay(y[:, 0], m[:, 0], lo, hi, win=4096, hop=4096, start=a)
        k = min(len(pr), len(me), len(ref))
        sel = np.arange(k)[:5]
        pred = (pr[:k, 1] - ref[:k, 1])[sel]; meas = me[:k, 1][sel]
        print("  %5d-%5d predicted (unwarp artefact) %s | measured %s | measured-predicted %s" % (lo, hi, np.round(pred, 3), np.round(meas, 3), np.round(meas - pred, 3)))
