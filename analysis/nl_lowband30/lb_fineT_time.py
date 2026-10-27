import pickle, numpy as np
FS = 26373.626373626373
D = pickle.load(open("fineT.pkl", "rb")); fr = D["fr"]; NF = (len(fr) - 1) * 2
for c in (0, 1):
    h = np.fft.irfft(D["T"][c] - 1, NF)
    t = (np.arange(NF) + NF // 2) % NF - NF // 2
    mag = np.abs(h); order = np.argsort(-mag)[:40]
    print(f"ch{c}: |h| total energy {np.sum(h**2):.3e}; largest |h| samples (lag samples = ms):")
    seen = []
    for i in order:
        if any(abs(t[i] - s) < 3 for s in seen): continue
        seen.append(t[i]); print(f"   lag {t[i]:6d} samples = {t[i]/FS*1e3:8.2f} ms   h = {h[i]:+.2e}")
        if len(seen) >= 12: break
    # energy distribution over lag bins
    for a, b in ((-10000, -2000), (-2000, -300), (-300, 0), (0, 300), (300, 2000), (2000, 10000), (10000, 32768)):
        sel = (t >= a) & (t < b)
        print(f"   lag {a:6d}..{b:<6d} energy fraction {np.sum(h[sel]**2)/np.sum(h**2):.3f}")
