"""Residual peaks vs each loop's modes. Settings: one key (argv), all its captures, R; residual and model spectra over
50-600 ms (Hann), E/M per bin; each loop's mode frequencies from the baseline's G_i (phase of G_i e^{-jwP_i} = 0 mod 2 pi).
model_s_results §25.12.

    EMT_KERNEL=chain3 .venv/bin/python tools/resid_common/modes.py <key>"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import dense_mid as dm
FS, PRE = dm.FS, dm.PRE
key = sys.argv[1]; lo, hi = 180, 345
z = np.load(os.path.join(dm.OUT, "prep", key + ".npz"))
caps = sorted(int(n[3:]) for n in z.files if n.startswith("res"))
te = min(float(z[f"te{j}"]) for j in caps)
i0, i1 = PRE + int(0.05 * FS), PRE + int(min(0.6, te) * FS)
win = np.hanning(i1 - i0); Nf = 1 << 18
E = np.array([np.fft.rfft((z[f"y{j}"][:, 1] - z[f"m{j}"][:, 1])[i0:i1] * win, Nf) for j in caps if int(z[f"amp{j}"]) > 0])
M = np.array([np.fft.rfft(z[f"m{j}"][:, 1][i0:i1] * win, Nf) for j in caps if int(z[f"amp{j}"]) > 0])
f = np.fft.rfftfreq(Nf, 1 / FS)
r = 10 * np.log10((np.abs(E) ** 2).mean(0) / (np.abs(M) ** 2).mean(0))
s = (f >= lo) & (f <= hi)
idx = np.nonzero(s)[0]
pk = [i for i in idx[1:-1] if r[i] >= r[i - 1] and r[i] >= r[i + 1] and r[i] > np.median(r[s]) + 4]
G = dm.loop_G(key, Nf)
w = 2 * np.pi * f / FS
modes = []
for i_, (Gi, Pi) in enumerate(zip(G, dm.q1.P)):
    ph = np.unwrap(np.angle(Gi * np.exp(-1j * w * Pi))[s])
    n = np.floor(ph / (2 * np.pi))
    for k in np.nonzero(np.diff(n))[0]:
        modes.append((f[idx[k]], i_))
print(f"{key}: {len(E)} captures, 50-{1000 * min(0.6, te):.0f} ms; median E/M {np.median(r[s]):.1f} dB in {lo}-{hi} Hz; resolution {FS / (i1 - i0):.1f} Hz")
for i in pk:
    near = sorted(modes, key=lambda m: abs(m[0] - f[i]))[:2]
    print(f"  peak {f[i]:6.1f} Hz  E/M {r[i]:6.1f} dB   nearest modes: " + ", ".join(f"loop {l} {fm:.1f} ({fm - f[i]:+.1f})" for fm, l in near))

# every bin: distance to the nearest mode, and the gap between the two nearest modes of different loops around it
mf = np.array([m[0] for m in modes]); ml = np.array([m[1] for m in modes])
rows = []
for i in idx:
    d = np.abs(mf - f[i]); o = np.argsort(d)
    j1 = o[0]; j2 = next(j for j in o[1:] if ml[j] != ml[j1])
    rows.append((r[i], d[j1], abs(mf[j1] - mf[j2]), d[j2]))
R_ = np.array(rows)
on = R_[:, 1] < 0.3                                           # on a mode
pair = on & (R_[:, 2] < 0.4)                                  # on a mode with another loop's mode within 0.4 Hz
single = on & (R_[:, 2] > 1.5)                                # on a mode, the other loops' modes > 1.5 Hz away
off = R_[:, 1] > 1.0
print(f"  E/M median (dB): on a coincident pair {np.median(R_[pair, 0]):.1f} ({pair.sum()} bins); on a lone mode {np.median(R_[single, 0]):.1f} ({single.sum()}); "
      f"off modes (> 1 Hz) {np.median(R_[off, 0]):.1f} ({off.sum()})")
for g in (0.2, 0.5, 1.0, 2.0, 4.0):
    sel = on & (R_[:, 2] < g) & (R_[:, 2] >= g / 2.5 if g > 0.2 else True)
    if sel.sum():
        print(f"    on a mode, other loop's mode gap < {g:g} Hz: E/M {np.median(R_[sel, 0]):.1f} dB ({sel.sum()} bins)")
