"""What fills 0.3-1.0 s of LATE 0.2 captures: capture, model (table) and residual levels per 50 ms, re the capture's
0.05-0.1 s level, R; the capture noise floor beside."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from late02 import Cap, nlt, nm, FS, K0
import numpy as np
T = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "heldout.json")))["table"]
_o = nlt.loop_gain
nlt.loop_gain = lambda P, late: T[str(P)] if late == 0.2 else _o(P, late)
noise = json.load(open(nm.NOISE))
for k in sys.argv[1:]:
    c = Cap(tuple(float(v) for v in k.split("_")))
    nm.spectra.cache_clear()
    m = nm.render(c.key, c.n, c.b["o"], c.b["front_shift"])
    g = (c.y[:, 1] @ m[:, 1]) / (m[:, 1] @ m[:, 1]); m = m[:, 1] * g; y = c.y[:, 1]
    z = K0 + int(c.b["front_shift"]); w = int(0.05 * FS)
    ref = np.mean(y[z + w:z + 2 * w] ** 2)
    row = []
    for t in np.arange(0.1, 1.0, 0.1):
        a = z + int(t * FS)
        if a + w > c.n: break
        L = lambda x: 10 * np.log10(np.mean(x[a:a + w] ** 2) / ref + 1e-30)
        row.append(f"{t:.1f}s cap {L(y):.0f} mdl {L(m):.0f} res {L(y - m):.0f}")
    print(k, f"noise {10*np.log10(noise['ch1']/ref):.0f} |", " | ".join(row))
