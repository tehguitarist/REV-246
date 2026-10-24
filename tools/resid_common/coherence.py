"""Coherence of the baseline's residual E with its own model M across captures (window a-b ms, R), per frequency bin:
gamma^2 = |sum M* E|^2 / (sum |M|^2 sum |E|^2): the fraction of E that one fixed transfer error (E = H M) explains; and
the same against the mean residual waveform (a fixed additive waveform: E = c W). Undamped, RT >= 1.0, AMPL > 0.
model_s_results §25.12.

    EMT_KERNEL=chain3 .venv/bin/python tools/resid_common/coherence.py <a,b ms>"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import dense_mid as dm
FS, PRE = dm.FS, dm.PRE
a, b = (float(x) for x in sys.argv[1].split(","))
ap = json.load(open(os.path.join(dm.q2.OUT, dm.ls.lb.lf_.AP230)))
info = dm.q2.q1_tables()[0]
i0, i1 = PRE + int(a / 1000 * FS), PRE + int(b / 1000 * FS)
win = np.hanning(i1 - i0)
Ms, Es = [], []
for p in sorted(os.listdir(os.path.join(dm.OUT, "prep"))):
    k = p[:-4]
    if ap["steps"][ap["settings"][k]["step"]] >= 0.01 or info[k]["rt"] < 1.0:
        continue
    z = np.load(os.path.join(dm.OUT, "prep", p))
    for j in sorted(int(n[3:]) for n in z.files if n.startswith("res")):
        y, m = z[f"y{j}"][:, 1].astype(float), z[f"m{j}"][:, 1].astype(float)
        if len(y) < i1 or int(z[f"amp{j}"]) == 0:
            continue
        Ms.append(np.fft.rfft(m[i0:i1] * win)); Es.append(np.fft.rfft((y - m)[i0:i1] * win))
M, E = np.array(Ms), np.array(Es)
f = np.fft.rfftfreq(i1 - i0, 1 / FS)
g_t = np.abs((M.conj() * E).sum(0)) ** 2 / ((np.abs(M) ** 2).sum(0) * (np.abs(E) ** 2).sum(0))
W = E.mean(0)
g_w = np.abs((W.conj()[None] * E).sum(0)) ** 2 / (len(E) * np.abs(W) ** 2 * (np.abs(E) ** 2).sum(0))
print(f"{a:g}-{b:g} ms, {len(E)} captures; per band: median over bins of coherence with the model (transfer) / with the mean residual (fixed waveform); E/M level dB")
for lo, hi in ((150, 200), (200, 250), (250, 300), (300, 350), (350, 500), (500, 700), (700, 1000), (1000, 1400), (2000, 3000), (3500, 5000)):
    s = (f >= lo) & (f < hi)
    print(f"  {lo}-{hi} Hz: transfer {np.median(g_t[s]):.2f}  fixed waveform {np.median(g_w[s]):.2f}  E/M {10 * np.log10((np.abs(E[:, s]) ** 2).sum() / (np.abs(M[:, s]) ** 2).sum()):.1f} dB")
print("fine, 180-340 Hz: bin Hz | E/M dB | transfer coh | |H| % / phase deg of sum M*E / sum |M|^2")
Hh = (M.conj() * E).sum(0) / (np.abs(M) ** 2).sum(0)
for i in np.nonzero((f >= 180) & (f <= 340))[0][::max(1, int(4 / (f[1] - f[0])))]:
    print(f"  {f[i]:6.1f} | {10 * np.log10((np.abs(E[:, i]) ** 2).sum() / (np.abs(M[:, i]) ** 2).sum()):6.1f} | {g_t[i]:.2f} | {100 * abs(Hh[i]):.3f} {np.degrees(np.angle(Hh[i])):+5.0f}")
