"""(b) spectral structure: 1/6-octave profile of the model error (total), its coherent part across LATE-only repeats, and the incoherent (noise) part.
SPAN 0.5, EARLY 1..10, window 0.05-0.45 s, Hann-tapered segment (suppresses window-edge leakage), series pre-band-limited 12-1600 Hz (order 6)."""
from cr_common import *
import itertools
B = bench_all()
EAR = (1.0, 1.5, 2.0, 3.0, 10.0)
sel = [k for k in sorted(B) if k[0] == 0.5 and k[1] in EAR and not B[k]["trimmed"]]
NF = 16384
fr = np.fft.rfftfreq(NF, 1 / FS)
PRE = butter(6, [12, 3000], btype='band', fs=FS, output="sos")
edges = 16 * 2 ** (np.arange(0, 60) / 6.0); edges = edges[edges <= 2600]
def db(x): return 10 * np.log10(x)
S = {}
for k in sel:
    r = B[k]; y, m = cached(k, r); a, b = window(k, r)
    w = np.hanning(b - a)[:, None]
    yf = sosfiltfilt(PRE, y, axis=0)[a:b] * w; mf = sosfiltfilt(PRE, m, axis=0)[a:b] * w
    S[k] = (np.fft.rfft(yf, NF, axis=0), np.fft.rfft(mf, NF, axis=0))
gi = np.digitize(fr, edges)
nb = len(edges) - 1
def bsum(x):               # (NF/2+1,) -> per band sums
    return np.array([x[gi == i + 1].sum() for i in range(nb)])
tot = {e: [] for e in EAR}; coh = {e: [] for e in EAR}; inc = {e: [] for e in EAR}; rep = {e: [] for e in EAR}
for e_ in EAR:
    ks = [k for k in sel if k[1] == e_]
    for k in ks:
        Y, M = S[k]
        for c in (0, 1):
            tot[e_].append(bsum(np.abs(Y[:, c] - M[:, c]) ** 2) / bsum(np.abs(Y[:, c]) ** 2))
    pairs = [(p, q) for p, q in itertools.combinations(ks, 2) if abs(B[p]["front_shift"] - B[q]["front_shift"]) < 0.3 and abs(B[p]["o"] - B[q]["o"]) < 0.02]
    for p, q in pairs:
        for c in (0, 1):
            Yp, Mp = S[p]; Yq, Mq = S[q]
            ep, eq = Yp[:, c] - Mp[:, c], Yq[:, c] - Mq[:, c]
            den = 0.5 * (bsum(np.abs(Yp[:, c]) ** 2) + bsum(np.abs(Yq[:, c]) ** 2))
            coh[e_].append(bsum((ep * np.conj(eq)).real) / den)
            # repeat: capture A against capture B (complex gain per 1/6-oct band would hide the error; use one real gain over the whole band)
            g = np.real(np.vdot(Yq[:, c], Yp[:, c])) / np.real(np.vdot(Yq[:, c], Yq[:, c]))
            rep[e_].append(bsum(np.abs(Yp[:, c] - g * Yq[:, c]) ** 2) / bsum(np.abs(Yp[:, c]) ** 2))
cf = np.sqrt(edges[:-1] * edges[1:])
print("centre Hz | per EARLY: total model error dB / coherent dB / pair-repeat dB (raw, 2 captures' noise)")
print("   Hz   " + "".join(f"  E={e_:<4}              " for e_ in EAR))
for i in range(nb):
    if cf[i] < 17: continue
    line = f"{cf[i]:7.1f} "
    for e_ in EAR:
        t = np.median([x[i] for x in tot[e_]]); c = np.mean([x[i] for x in coh[e_]]); r_ = np.median([x[i] for x in rep[e_]])
        line += f" {db(t):6.1f} {db(max(c, 1e-9)):6.1f} {db(r_):6.1f} |"
    print(line)
np.savez(os.path.join(HERE, "b_profile.npz"), cf=cf, **{f"tot{e}": np.array(tot[e]) for e in EAR}, **{f"coh{e}": np.array(coh[e]) for e in EAR})
