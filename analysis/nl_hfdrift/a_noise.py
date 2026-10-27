"""Fork A part 2+3b: noise floor per sub-band (LATE 0 tails, new warp) vs the signal in the 0.05-0.6 s windows of ~N captures, and the model's
null on the same windows ('exact gain only', bench_warp2_bk_kf_a0.7 fits). usage: a_noise.py [n_captures]"""
import os, sys, json, random
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from a_repeat import *          # sets NL_WARP2 / NL_BK, common.BENCH, BANDS, bp, nulls, ...
import defects
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
NB = [k for k in BANDS if k != "fit"]
# --- noise: LATE 0 tails
nz = {k: [[], []] for k in NB}; nf = 0
for key, p in sorted(nl.table().items()):
    if key[2] != 0.0 or key in defects.NL_DEFECTS: continue
    try: b = common.bench(key)
    except FileNotFoundError: continue
    raw = nl.ir(p); h = M.sosfiltfilt(ne.BAND, raw, axis=0)
    ns = ne.N_STAGES[key[0]]; i0 = M.K0 + (ne.taps(0, ns)[-1][0] if ns else 600) + int(0.25 * FS); i1 = M.usable_end(h)
    if i1 - i0 < int(0.1 * FS): continue
    y = nlt.unwarp(np.vstack([raw[:i1 + 64], np.zeros((max(0, i1 + 64 - len(raw)), 2))]), b["ppm"])[:i1]
    for k in NB:
        z = bp(y[i0 - 200:i1], k)[200:]
        for c in (0, 1): nz[k][c].append(float(np.mean(z[:, c] ** 2)))
    nf += 1
noise = {k: [float(np.median(nz[k][c])) for c in (0, 1)] for k in NB}
print("noise tails from", nf, "LATE-0 files; per-sample noise power (dB re 1) L/R:", {k: [round(10 * np.log10(v), 1) for v in noise[k]] for k in NB})
print("spread (p10-p90, dB) of tail noise per band, L:", {k: round(10 * np.log10(np.percentile(nz[k][0], 90) / np.percentile(nz[k][0], 10)), 1) for k in NB})
# --- captures
random.seed(5)
ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
rows = []
for key in random.sample(ks, N):
    b = common.bench(key); raw = nl.ir(nl.table()[key]); h = M.sosfiltfilt(ne.BAND, raw, axis=0); fr = b["front_shift"]
    i0, tr = M.window_start(h, fr)
    a, e = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if e - a < int(0.05 * FS): continue
    n = e + PAD + 80
    rawp = np.vstack([raw[:n], np.zeros((max(0, n - len(raw)), 2))])
    y = nlt.unwarp(rawp, b["ppm"])[a - PAD:e + PAD]
    m = (M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"]))[a - PAD:e + PAD]
    mn = nulls(y, m, NB); mcg = nulls(y, m, NB, cg=True)
    r = dict(key=key, L=e - a)
    for k in NB:
        sig = [float(np.sum(bp(y, k)[PAD:-PAD][:, c] ** 2)) for c in (0, 1)]
        r[k] = dict(model=mn[k], model_cg=mcg[k], snr=[10 * np.log10(sig[c] / (noise[k][c] * (e - a))) for c in (0, 1)],
                    noise_null=[10 * np.log10(noise[k][c] * (e - a) / sig[c]) for c in (0, 1)])
    rows.append(r)
json.dump(dict(noise=noise, rows=rows), open(OUT + "/a_noise.json", "w"), default=float)
print(len(rows), "captures. Medians over captures whose band is >= 10 dB above the noise (L, R); model null, noise-limited null (single capture vs model), excess of model over noise")
for k in NB:
    sel = [r for r in rows if min(r[k]["snr"]) >= 10]
    if not sel: print(k, "none above noise"); continue
    mo = np.array([r[k]["model"] for r in sel]); no = np.array([r[k]["noise_null"] for r in sel]); cg = np.array([r[k]["model_cg"] for r in sel])
    print(f"{k:11s} n {len(sel):3d}/{len(rows)} model {np.round(np.median(mo, 0), 1)}  model+cplx gain {np.round(np.median(cg, 0), 1)}  noise-limited {np.round(np.median(no, 0), 1)}  excess {np.round(np.median(mo - no, 0), 1)}")
