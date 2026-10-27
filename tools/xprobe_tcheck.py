"""Score the free fixed transfer T (model_s_results 41.3, fold_T) truncated to the probe's short lags inside the probe's
normal equations, to compare the time-domain proxy with the bench: does the proxy see T's real gain?

    EMT_KERNEL=chain3 .venv/bin/python tools/xprobe_tcheck.py analysis/fullband/xprobe/groups_ch.npz
T of fold A is evaluated on fold B's blocks and vice versa (as the bench), common to all four loops, per channel."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import sml_fullband as sf

Z = np.load(sys.argv[1])
GR, LAGS = list(Z["groups"]), list(Z["lags"])
ATA, ATb, BTB = Z["ATA"], Z["ATb"], Z["BTB"]
NLAG = len(LAGS)
SHORT = [l for l in LAGS if abs(l) <= 10 or l in [s * (e + d) for e in (230, 460) for d in range(-2, 3) for s in (1, -1)]]
TB, NB_ = sf.TBIN, 6150
print("TBIN", TB, "FS", sf.FS)
for kind, keep in (("short (|lag|<=10 + AP230 echoes: the scored 'short_ap')", SHORT), ("|lag| <= 10 only", [l for l in LAGS if abs(l) <= 10]), ("all probe lags", LAGS)):
    print("\nT truncated to", kind)
    for tr, te in (("A", "B"), ("B", "A")):
        T = sf.fold_T(tr)                       # (2, TNB) complex, fold `tr` LF >= 1.0
        out = np.zeros((4, 3))
        for c in (0, 1):
            nfft = int(round(sf.FS / TB))
            full = np.ones(nfft // 2 + 1, complex)
            full[:T.shape[1]] = T[c]
            ir = np.fft.irfft(full - 1, nfft)
            h = np.zeros(4 * NLAG)
            for i in range(4):
                for j, l in enumerate(LAGS):
                    if l in keep:
                        h[i * NLAG + j] = ir[l % nfft]
            gs = [i for i, g in enumerate(GR) if g.split("|")[0] == te and int(g.split("|")[2]) == c]
            for w in range(4):
                for b in range(3):
                    A, bb, t = ATA[gs, w, b].sum(0), ATb[gs, w, b].sum(0), BTB[gs, w, b].sum()
                    out[w, b] += 10 * np.log10((t - 2 * h @ bb + h @ A @ h) / t) / 2
        print(f"  T from fold {tr} on fold {te} (mean of L and R), dB per window x band(<400 / 400-3200 / >3200 Hz):")
        for w, wn in enumerate(("0-50", "50-150", "150-300", "300-600")):
            print(f"    {wn:8s}" + "  ".join(f"{out[w, b]:+7.2f}" for b in range(3)))
