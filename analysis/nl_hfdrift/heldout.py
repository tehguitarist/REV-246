"""Held-out test: does a long-kernel clock warp remove the 11-12.3 kHz drift and the 10-12 kHz null excess?
Per capture: bench's own o, ppm, gain, front shift (read only). Windows as the bench's per-band check: 0.05-0.6 s.
Variants of the capture's clock correction: old = nlt.unwarp (bench), new = unwarp2. The model is the bench's model with
its EQ off (render with EQ off, then either the bench's EQ or a pooled EQ refit on the TRAIN half applied in this script).
Outputs analysis/nl_hfdrift/heldout.json"""
import sys, json, random, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
from scipy.signal import butter, sosfiltfilt
M.EQ = False; M._proc.cache_clear()
eqj = json.load(open(M.EQ_FILE))
random.seed(7)
allk = sorted(nl.table().keys())
import defects
bad = set(getattr(defects, "NL_DEFECTS", []))
ks = [k for k in allk if k not in bad]
ks = [k for k in ks if (nl.table() and True)]
sample = random.sample(ks, int(sys.argv[1]) if len(sys.argv) > 1 else 80)
BANDS2 = [(8000, 10000), (10000, 11000), (11000, 12000), (11000, 11600), (11600, 12300), (12300, 12800), (10000, 12000)]
rec = []
def tau_mean(y, m, lo, hi, a, b2):
    q = local_delay(y, m, lo, hi, win=1024, hop=1024, start=a, stop=b2)
    return float(np.mean(q[:, 1])) if q.ndim == 2 and len(q) else float("nan")
for i, key in enumerate(sample):
    try:
        b = bench(key)
    except FileNotFoundError:
        continue
    raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, trimmed = M.window_start(h, fr)
    a, b2 = max(i0, M.K0 + int(fr) + int(0.05 * FS)), min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    if b2 - a < int(0.05 * FS):
        continue
    n = b2 + 200
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    yo = nlt.unwarp(rawp, b["ppm"])[:n]
    yn = unwarp2(rawp, b["ppm"])[:n]
    m = M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"])      # EQ off here
    r = dict(key=key, a=a, b=b2, trimmed=bool(trimmed), bench_10_12=b["bands"].get("10000-12000"))
    # per-band delay of old/new vs model in the 3 HF bands, over the whole window, per channel
    for nm, y in (("old", yo), ("new", yn)):
        for lo, hi in ((11000, 11600), (11600, 12300), (12300, 12800)):
            r[f"tau_{nm}_{lo}"] = [tau_mean(y[:, c], m[:, c], lo, hi, a, b2) for c in (0, 1)]
    # spectra of the window for pooled EQ fit and for scoring
    seg = slice(a, b2); L = b2 - a; w = np.hanning(L)
    Y = {nm: np.fft.rfft(y[seg] * w[:, None], axis=0) for nm, y in (("old", yo), ("new", yn))}
    Mm = np.fft.rfft(m[seg] * w[:, None], axis=0)
    f = np.fft.rfftfreq(L) * FS
    r["spec"] = dict(f=f, M=Mm, Y=Y)
    r["y"] = {"old": yo[seg], "new": yn[seg]}; r["m"] = m[seg]
    rec.append(r)
    print(i, key, "n", L, flush=True)
print("captures", len(rec))
# --- pooled EQ fit on TRAIN (even index) and TEST (odd), real zero-phase gain on knots
knots = np.array([9000, 10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0])
def fit_eq(recs, which, CPLX=False):
    num = np.zeros(len(knots), complex if CPLX else float); den = np.zeros(len(knots))
    for r in recs:
        f = r["spec"]["f"]; Mm = r["spec"]["M"]; Y = r["spec"]["Y"][which]
        for ci in (0, 1):
            x = (Y[:, ci] * np.conj(Mm[:, ci])) if CPLX else np.real(Y[:, ci] * np.conj(Mm[:, ci])); d = np.abs(Mm[:, ci]) ** 2
            # hat-function weights onto knots
            for kk in range(len(knots)):
                lo = knots[kk - 1] if kk else knots[0] - 1; hi = knots[kk + 1] if kk < len(knots) - 1 else knots[-1] + 1
                wgt = np.clip(np.minimum((f - lo) / (knots[kk] - lo), (hi - f) / (hi - knots[kk])), 0, 1)
                num[kk] += np.sum(wgt * x); den[kk] += np.sum(wgt * d)
    return num / den
def apply_eq(r, which, g):
    f = r["spec"]["f"]; E = np.interp(f, knots, g)
    return np.fft.irfft(r["spec"]["M"] * E[:, None], len(r["m"]) , axis=0)  # windowed by hanning: not used for scoring
def score(r, y, mm, lo, hi):
    sos = butter(4, [lo, hi], btype="band", fs=FS, output="sos")
    yb, mb = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, mm, axis=0)
    return [float(10 * np.log10(np.sum((yb[:, c] - mb[:, c]) ** 2) / np.sum(yb[:, c] ** 2))) for c in (0, 1)]
def eq_time(r, g):
    """EQ applied to the (unwindowed) model segment in the frequency domain, zero phase."""
    m = r["m"]; L = len(m); f = np.fft.rfftfreq(L) * FS; E = np.interp(f, knots, g.real) + 1j * np.interp(f, knots, g.imag)
    return np.fft.irfft(np.fft.rfft(m, axis=0) * E[:, None], L, axis=0)
train = [r for i, r in enumerate(rec) if i % 2 == 0]; test = [r for i, r in enumerate(rec) if i % 2 == 1]
g_old_tr, g_new_tr = fit_eq(train, "old") + 0j, fit_eq(train, "new") + 0j
g_old_te, g_new_te = fit_eq(test, "old") + 0j, fit_eq(test, "new") + 0j
c_old_tr, c_new_tr = fit_eq(train, "old", True), fit_eq(train, "new", True)
c_old_te, c_new_te = fit_eq(test, "old", True), fit_eq(test, "new", True)
print('complex new-unwarp EQ (train) phase mrad:', np.round(np.angle(c_new_tr)*1e3, 1)); print('complex new-unwarp EQ (test)  phase mrad:', np.round(np.angle(c_new_te)*1e3, 1)); print('complex old-unwarp EQ (train) phase mrad:', np.round(np.angle(c_old_tr)*1e3, 1))
print("EQ knots (kHz):", knots / 1000)
print("old-unwarp EQ gain dB (train/test):", np.round(20*np.log10(np.abs(g_old_tr)), 2), np.round(20*np.log10(np.abs(g_old_te)), 2))
print("new-unwarp EQ gain dB (train/test):", np.round(20*np.log10(np.abs(g_new_tr)), 2), np.round(20*np.log10(np.abs(g_new_te)), 2))
print("bench's out_eq.json (dB at 11.5/12/12.3/12.6k):", np.round(np.interp([11500, 12000, 12300, 12600], eqj["f_hz"], eqj["db"]), 2))
res = {}
flat = np.ones(len(knots), complex)
for tag, rs, g_o, g_n, c_o, c_n in (("test", test, g_old_tr, g_new_tr, c_old_tr, c_new_tr), ("train", train, g_old_tr, g_new_tr, c_old_tr, c_new_tr)):
    for lo, hi in BANDS2:
        row = {}
        for nm, y_key, g in (("old/no-EQ", "old", flat), ("old/EQ(train)", "old", g_o), ("new/no-EQ", "new", flat), ("new/EQ(train)", "new", g_n), ("old/cEQ(train)", "old", c_o), ("new/cEQ(train)", "new", c_n)):
            v = [score(r, r["y"][y_key], eq_time(r, g), lo, hi) for r in rs]
            row[nm] = (float(np.median([x[0] for x in v])), float(np.median([x[1] for x in v])))
        res[f"{tag} {lo}-{hi}"] = row
        print(tag, f"{lo}-{hi}", {k: tuple(round(x, 1) for x in v) for k, v in row.items()})
# the bench's reported 10-12 kHz (EQ from out_eq.json, old unwarp), same captures
bn = [r["bench_10_12"] for r in test if r["bench_10_12"]]
print("bench json 10-12k median on test half:", np.median([x[0] for x in bn]), np.median([x[1] for x in bn]))
# delays
for nm in ("old", "new"):
    for lo in (11000, 11600, 12300):
        v = np.array([np.mean(r[f"tau_{nm}_{lo}"]) for r in rec])
        print("mean delay y-m %s band %d: median %.4f  (p10 %.4f p90 %.4f) samples" % (nm, lo, np.median(v), np.percentile(v, 10), np.percentile(v, 90)))
json.dump({"res": res, "knots": knots.tolist(), "g_old_tr": g_old_tr.real.tolist(), "g_new_tr": g_new_tr.real.tolist(), "g_new_te": g_new_te.real.tolist(), "c_new_tr": [[z.real, z.imag] for z in c_new_tr], "c_new_te": [[z.real, z.imag] for z in c_new_te], "c_old_tr": [[z.real, z.imag] for z in c_old_tr],
           "keys": [list(r["key"]) for r in rec]}, open("/Users/leighpierce/Git/246/analysis/nl_hfdrift/heldout.json", "w"))
