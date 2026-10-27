"""Does the first-30 ms LF residual look like the acausal leak of the branch (main taps + late) through the processing response?"""
import sys
from lb_common import *
B = bench_all()
sosA, sosB = bands(20, 100), bands(100, 1000)
want = [(0.0, 10.0, 0.0), (0.1, 10.0, 0.0), (0.2, 10.0, 0.4), (0.3, 10.0, 0.4), (0.4, 10.0, 0.4), (0.5, 10.0, 0.4), (0.3, 1.0, 0.4), (0.3, 0.3, 0.4), (0.3, 3.0, 1.0), (0.5, 3.0, 1.0)]
print("key | fs | window [file+20, z+25ms]: rms dB re 410 peak:  y  leak(branch in model)  e(res) | by band 20-100 / 100-1k: e/leak (dB)  e vs y-front")
for k in want:
    if k not in B: print("missing", k); continue
    r = B[k]; n = 2400
    y, m, _ = load(k, r, n=n)
    fm, bm = render_parts(k, n, r["o"], r["front_shift"])
    g = np.array(r["gain"]); fm, bm = fm * g, bm * g
    z = int(round(K0 + r["front_shift"]))
    ia, ib = K0 + 20, int(z + 0.025 * FS)
    line = f"{k} fs {r['front_shift']:6.1f} |"
    for c in (0,):
        pk = np.abs(m[z:z + 800, c]).max()
        e = y[:, c] - m[:, c]
        line += f" ch{c}"
        for nm, sos in (("20-100", sosA), ("100-1k", sosB)):
            lk = sosfiltfilt(sos, bm[:, c])[ia:ib]; ee = sosfiltfilt(sos, e)[ia:ib]; yy = sosfiltfilt(sos, y[:, c])[ia:ib]
            ref = sosfiltfilt(sos, y[:, c] - fm[:, c])[ia:ib]
            q = lambda v: 20 * np.log10(np.sqrt(np.mean(v ** 2)) / pk + 1e-12)
            cc = np.dot(lk, ee) / np.sqrt(np.dot(lk, lk) * np.dot(ee, ee))
            line += f" | {nm}: y {q(yy):6.1f} leak {q(lk):6.1f} e {q(ee):6.1f} (e-leak {q(ee)-q(lk):5.1f} dB, corr {cc:5.2f})"
    print(line)
