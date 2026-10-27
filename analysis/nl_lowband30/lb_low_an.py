import pickle, itertools, collections, numpy as np
from scipy.signal import butter, sosfiltfilt
FS = 26373.626373626373; K0 = 50
S = pickle.load(open("low_sel.pkl", "rb"))
def win(x, z, t0, t1):
    a, b = int(z + t0 * FS), int(z + t1 * FS)
    return x[a:b] * np.hanning(b - a)[:, None]
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
def spec(d, t0=0.05, t1=0.45):
    return np.fft.rfft(win(d["y"], d["z"], t0, t1), NF, axis=0), np.fft.rfft(win(d["m"], d["z"], t0, t1), NF, axis=0)
sp = {k: spec(d) for k, d in S.items()}
# repeat pairs: same EARLY, same fs, same o
pairs = [(a, b) for a, b in itertools.combinations(sorted(S), 2) if a[1] == b[1] and abs(S[a]["o"] - S[b]["o"]) < 0.02 and abs(S[a]["fs"] - S[b]["fs"]) < 0.3]
print(len(S), "captures,", len(pairs), "repeat pairs")
edges = [10, 15, 20, 25, 32, 40, 50, 63, 80, 100, 130, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
print("band(Hz)    model err |E|/|M| dB: rms over caps | coherent mean (|.|, phase deg) | repeat floor |dY|/(sqrt2 |M|) dB  | gap rms-floor")
for lo, hi in zip(edges[:-1], edges[1:]):
    mk = (fr >= lo) & (fr < hi)
    rat = []; fl = []
    for k, (Y, Mo) in sp.items():
        for c in (0, 1):
            den = np.sum(np.abs(Mo[mk, c]) ** 2)
            rat.append(np.sum((Y[mk, c] - Mo[mk, c]) * Mo[mk, c].conj()) / den)      # complex relative error per capture/channel
    rat = np.array(rat)
    for a, b in pairs:
        for c in (0, 1):
            Ya, Ma = sp[a][0][mk, c], sp[a][1][mk, c]; Yb = sp[b][0][mk, c]
            # relative repeat difference after per-band LS gain (the repeat fit freedom)
            g = np.vdot(Yb, Ya) / np.vdot(Yb, Yb)
            fl.append(np.sum(np.abs(Ya - g * Yb) ** 2) / np.sum(np.abs(Ya) ** 2) / 2)
    # model error energy relative to Y energy
    errE = np.array([np.sum(np.abs(sp[k][0][mk, c] - sp[k][1][mk, c]) ** 2) / np.sum(np.abs(sp[k][0][mk, c]) ** 2) for k in sp for c in (0, 1)])
    coh = rat.mean()
    print(f"{lo:5d}-{hi:<5d} err {10*np.log10(np.median(errE)):6.1f} dB | coherent |E/M| {20*np.log10(abs(coh)):6.1f} / {np.degrees(np.angle(coh)):5.0f} deg  (rms {20*np.log10(np.sqrt(np.mean(abs(rat)**2))):6.1f}) | repeat floor {10*np.log10(np.median(fl)):6.1f} | gap {10*np.log10(np.median(errE))-10*np.log10(np.median(fl)):5.1f}")
