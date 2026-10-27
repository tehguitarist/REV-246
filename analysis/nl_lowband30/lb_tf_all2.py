"""Per-capture, per-channel complex data/model ratio on a log grid 12-200 Hz for all untrimmed LATE 0 captures (EARLY >= 0.3, plus ramps with enough signal),
Hann window z .. min(z+1 s, usable end).  Saves tf_all.npz (F, num, den per capture/channel) for pooling by SPAN / session / EARLY.
env: NL_BK=1 NL_BK_ALPHA=0.65 LB_BENCH=bench_bk_a0.65"""
import sys, json
from lb_common import *
B = bench_all()
F = np.geomspace(12, 200, 36)
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
rows = []
keys = [k for k in sorted(B) if k[2] == 0.0 and not B[k]["trimmed"]]
for k in keys:
    r = B[k]; fs_ = r["front_shift"]
    raw = nl.ir(nl.table()[k]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    end = M.usable_end(h); z = K0 + fs_
    b = int(min(z + 1.0 * FS, end)); a = int(z + 0.05 * FS)
    if b - a < int(0.3 * FS) or k[0] < 0.2: continue
    n = b + 64
    y, m, _ = load(k, r, n=n, band=False)
    w = np.hanning(b - a)[:, None]
    Y = np.fft.rfft(y[a:b] * w, NF, axis=0); Mo = np.fft.rfft(m[a:b] * w, NF, axis=0)
    num = np.zeros((len(F), 2), complex); den = np.zeros((len(F), 2))
    for i, f in enumerate(F):
        mk = (fr >= f / 1.12) & (fr < f * 1.12)
        for c in (0, 1):
            num[i, c] = np.sum(Y[mk, c] * Mo[mk, c].conj()); den[i, c] = np.sum(abs(Mo[mk, c]) ** 2)
    rows.append(dict(key=k, fs=fs_, num=num, den=den, win=(b - a) / FS))
    print(k, round((b - a) / FS, 2), flush=True)
import pickle; pickle.dump(dict(F=F, rows=rows), open("tf_all2.pkl", "wb"))
