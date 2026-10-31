"""Per-pass amplitudes of the three late loops (NonLin_program 10.4). Per capture: the bench's alignment (o, clock, front
shift); the model's front part and main-tap branch (late left out) each with one free gain per channel; every loop pass
k = 0..K-1 a free amplitude (its 4 signed read taps through AP(124) + AP(314) as the model). The ratio a_k+1 / a_k is the
loop's gain per pass, measured model-free of the loop law.
    passes.py <LATE> [K] [every]  -> analysis/nl/late02/passes_<LATE>.json"""
import json, os, sys
os.environ.setdefault("NL_BK", "1"); os.environ.setdefault("EMT_KERNEL", "chain3")
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.signal import sosfiltfilt
import nl_structure as nl, nl_early as ne, nl_late as nlt, nl_model as nm
from late02 import keys_for, kstr
FS, K0 = nl.FSD, nl.K0


def passes(key, K):
    b = json.load(open(os.path.join(nm.BENCH, kstr(key) + ".json")))
    raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr, o = b["front_shift"], b["o"]
    n = min(nm.usable_end(h), K0 + int(fr) + int(0.6 * FS))
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = sosfiltfilt(ne.BAND, nlt.unwarp_from(nlt.upsample(rawp), n, b["ppm"])[:n], axis=0)
    N = 1 << 17
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    u = np.zeros(N); u[0] = 1
    U = np.fft.rfft(ne.ap(ne.ap(u, ne.AP_IN), nlt.AP_LATE))
    out = {}
    for ch in (0, 1):
        F, B = nm.spectra((key[0], key[1], 0.0), ch, N)
        cols = [np.fft.irfft(F * np.exp(ph * (K0 + fr)), N)[:n], np.fft.irfft(B * np.exp(ph * (K0 + fr + o)), N)[:n]]
        base = nlt.FEED[key[0]] + (nlt.R_LATE if ch else 0)
        for P, qs in nlt.READS.items():
            for k in range(K):
                T = sum(s * np.exp(ph * (K0 + fr + o + base + q + k * P)) for s, q in zip(nlt.READ_SIGNS[ch][P], qs))
                cols.append(np.fft.irfft(U * T, N)[:n])
        A = sosfiltfilt(ne.BAND, np.array(cols).T, axis=0)
        x, *_ = np.linalg.lstsq(A, y[:, ch], rcond=None)
        a = x[2:].reshape(3, K)
        out[ch] = dict(a=a.tolist(), ratio=[float(np.median(r[1:3] / r[0:2])) for r in a], a0=a[:, 0].tolist())
    return out


if __name__ == "__main__":
    late = float(sys.argv[1]); K = int(sys.argv[2]) if len(sys.argv) > 2 else 6; every = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    res = {}
    for k in keys_for(late, every):
        r = passes(k, K)
        res[kstr(k)] = r
        print(kstr(k), " ".join(f"{v:.4f}" for v in r[0]["ratio"]), "|", " ".join(f"{v:.4f}" for v in r[1]["ratio"]), flush=True)
    json.dump(dict(late=late, K=K, nominal=[nlt.loop_gain(P, late) for P in nlt.LOOPS], captures=res),
              open(os.path.join(nl.OUT, "late02", f"passes_{late}.json"), "w"))
    R = np.array([[v[c]["ratio"] for c in (0, 1)] for v in res.values()])
    print("nominal", [round(nlt.loop_gain(P, late), 4) for P in nlt.LOOPS], "median L", np.median(R[:, 0], 0).round(4), "R", np.median(R[:, 1], 0).round(4),
          "IQR L", (np.percentile(R[:, 0], 75, 0) - np.percentile(R[:, 0], 25, 0)).round(4))
