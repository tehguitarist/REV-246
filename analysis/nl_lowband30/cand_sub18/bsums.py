"""N cross sums per capture and window type, 2 Hz bins on S's grid (as lb_joint_sums, same convention), NL_LOWBAND=0 model, bench_warp2_bk_kf_a0.7 alignment.
Window types: w45f/w60f = front-zero+0.05 s ... min(0.45|0.6 s, late feed - 10 ms, usable end); w45n/w60n = same without the feed limit. Resumable."""
import pickle, json, time
from common import *
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS); TNB = int(12300 / 2.0); idx = (fr / 2.0).astype(int); ok = idx < TNB
path = os.path.join(HERE, "bsums.pkl"); out = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
ks = [k for k in sorted(B) if k[0] >= 0.1 and not B[k]["trimmed"] and abs(B[k]["front_shift"] - 26.2) >= 1 and abs(B[k]["front_shift"] - 29.5) >= 1]
print(len(ks), "candidates", flush=True); t0 = time.time()
for i, k in enumerate(ks):
    if k in out: continue
    r = B[k]; raw = nl.ir(nl.table()[k]); h = sosfiltfilt(ne.BAND, raw, axis=0); end = M.usable_end(h); z = K0 + r["front_shift"]; a = int(z + 0.05 * FS)
    feed = z + nlt.FEED[k[0]] + r["o"] - 0.01 * FS
    W = {}
    for nm, tl, lim in (("w45f", 0.45, True), ("w60f", 0.6, True), ("w45n", 0.45, False), ("w60n", 0.6, False)):
        b = min(z + tl * FS, end)
        if lim and k[2] > 0: b = min(b, feed)
        if b - a >= 0.1 * FS: W[nm] = (a, int(b))
    if not W: out[k] = None; continue
    nmax = max(b for a, b in W.values()) + 64
    y, m, _ = load(k, r, n=nmax, band=False)
    res = {}
    for nm, (a_, b_) in W.items():
        if nm in res: continue
        w = np.hanning(b_ - a_)[:, None]; Y = np.fft.rfft(y[a_:b_] * w, NF, axis=0); Mo = np.fft.rfft(m[a_:b_] * w, NF, axis=0); ym = Y * Mo.conj()
        acc = lambda zz: np.stack([np.bincount(idx[ok], zz[ok, c], TNB) for c in (0, 1)])
        res[nm] = (acc(ym.real) + 1j * acc(ym.imag), acc(abs(Mo) ** 2), acc(abs(Y) ** 2), b_ - a_)
    fl = {w: [B[k]["windows"][w].get(f"floor{c}") for c in (0, 1)] for w in B[k]["windows"]}
    out[k] = dict(span=k[0], early=k[1], late=k[2], sess=sess(r), win=res, floor=fl, same=(W["w45f"] == W["w45n"]) if "w45f" in W and "w45n" in W else None)
    if i % 25 == 0:
        pickle.dump(out, open(path, "wb")); print(i, "%.0f s" % (time.time() - t0), flush=True)
pickle.dump(out, open(path, "wb")); print("done", len(out))
