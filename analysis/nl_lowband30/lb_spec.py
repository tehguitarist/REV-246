"""Relative error E/M = Y/M - 1 (front section, [0, z+25 ms]) per group: coherent mean (fixed part) and rms (incl. noise/unfixed),
in 1/3-octave-ish bands; phase of the coherent mean."""
import pickle, collections, numpy as np
D = pickle.load(open("detail.pkl", "rb")); freqs = D["freqs"]; data = D["data"]
edges = [20, 31.5, 50, 80, 125, 200, 315, 500, 800, 1250, 2000, 3150]
def bandmean(vals, lo, hi):
    m = (freqs >= lo) & (freqs < hi)
    return vals[:, m]
for c in (0, 1):
    print(f"\n=== channel {'LR'[c]}:  band: coherent-mean |E/M| dB / phase deg  | rms |E/M| dB (rms over captures, power-avg over band, weighted by |M|^2)")
    hdr = "band(Hz)    " + "".join(f"{g:>26s}" for g in ("fs0", "fs16", "fs19.6", "fs22.9", "fs36"))
    print(hdr)
    rows = {}
    for g in ("fs0", "fs16", "fs19.6", "fs22.9", "fs36"):
        ks = [kk for kk in data if kk[0] == g]
        if not ks: continue
        Y = np.array([data[kk]["out"][c][0] for kk in ks]); M = np.array([data[kk]["out"][c][1] for kk in ks])
        rows[g] = (Y, M)
    for lo, hi in zip(edges[:-1], edges[1:]):
        line = f"{lo:6.1f}-{hi:6.0f} "
        for g in ("fs0", "fs16", "fs19.6", "fs22.9", "fs36"):
            if g not in rows: line += " " * 26; continue
            Y, M = rows[g]; mk = (freqs >= lo) & (freqs < hi)
            E = Y[:, mk] - M[:, mk]
            # per-capture band-weighted complex ratio  <E M*>/<M M*>
            ratio = (E * M[:, mk].conj()).sum(1) / (np.abs(M[:, mk]) ** 2).sum(1)
            coh = ratio.mean()
            rms = np.sqrt(np.mean(np.abs(ratio) ** 2))
            line += f"   {20*np.log10(abs(coh)+1e-12):6.1f}/{np.degrees(np.angle(coh)):5.0f}  | {20*np.log10(rms):6.1f}"
        print(line)
