"""Per-capture transfer sums of Program N on S's 2 Hz grid (sml_fullband convention: TN = 2^16 rfft at the device rate, 2 Hz bins up to 12.3 kHz), on the CURRENT N model
(warp2 + adopted branch kernel 0.7): window front-zero + 0.05 s ... min(0.45 s, late feed - 10 ms, usable end), Hann; SPAN >= 0.2, LATE <= 0.2, EARLY >= 0.3, untrimmed.
Saves joint_N_sums.pkl: {key: (span, sym (2, 6150) complex, smm (2, 6150), syy)}.
env: NL_BK=1 NL_BK_ALPHA=0.7 LB_BENCH=bench_warp2_bk_kf_a0.7"""
import pickle, collections
from lb_common import *
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS); TNB = int(12300 / 2.0); idx = (fr / 2.0).astype(int); ok = idx < TNB
FEED = {0.2: 3834, 0.3: 6696, 0.4: 9270, 0.5: 12254}
out = {}
for k in sorted(B):
    r = B[k]
    if k[0] < 0.2 or k[1] < 0.3 or k[2] > 0.2 or r["trimmed"] or abs(r["front_shift"] - 26.2) < 1 or abs(r["front_shift"] - 29.5) < 1: continue
    raw = nl.ir(nl.table()[k]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h); z = K0 + r["front_shift"]
    b = min(z + 0.45 * FS, end)
    if k[2] > 0: b = min(b, z + FEED[k[0]] + r["o"] - 0.01 * FS)
    a = z + 0.05 * FS
    if b - a < 0.1 * FS: continue
    a, b = int(a), int(b)
    y, m, _ = load(k, r, n=b + 64, band=False); w = np.hanning(b - a)[:, None]
    Y = np.fft.rfft(y[a:b] * w, NF, axis=0); Mo = np.fft.rfft(m[a:b] * w, NF, axis=0)
    ym = Y * Mo.conj()
    acc = lambda zz: np.stack([np.bincount(idx[ok], zz[ok, c], TNB) for c in (0, 1)])
    out[k] = (k[0], acc(ym.real) + 1j * acc(ym.imag), acc(abs(Mo) ** 2), acc(abs(Y) ** 2))
print(len(out), "captures", collections.Counter(v[0] for v in out.values()))
pickle.dump(out, open("joint_N_sums.pkl", "wb"))
