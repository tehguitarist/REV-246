"""Fork B: is the 11-12.8 kHz residual (after the warp fix and the exact resampler gain) structure independent?
Complex pooled residual E(f) fitted separately per group (time window, LATE, SPAN, EARLY class, channel, trimmed), bootstrap spread
over captures, and a held-out (capture halves) comparison of group-specific vs pooled fits.
usage: b_groups.py N_CAPTURES [out_prefix]   (writes analysis/nl_hfdrift/b_groups_<prefix>.json / .txt)"""
import os, sys, json, random, pickle
os.environ["NL_WARP2"] = "1"
os.environ.pop("NL_BK", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
import defects
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
TAG = sys.argv[2] if len(sys.argv) > 2 else "run"
WINS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.6)]
knots = np.array([10000, 10500, 11000, 11250, 11500, 11750, 12000, 12250, 12500, 12800, 13200.0]); K = len(knots)
SHOW = [4, 5, 6, 7, 8, 9]                       # 11.5 ... 12.8 kHz
BANDS = dict(S1=(11000, 11600), S2=(11600, 12300), S3=(12300, 12800), C=(10000, 12000))
bk = lambda k: common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"
ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and os.path.exists(bk(k)))
random.seed(5)
meta = {k: common.bench(k) for k in ks}
base = set(random.sample(ks, int(N * 0.65)))
for name, pred, n_min in (("LATE0", lambda k: k[2] == 0, 0.2), ("EARLY10", lambda k: k[1] >= 10, 0.2), ("trim", lambda k: meta[k]["trimmed"], 0.2), ("SPAN0", lambda k: k[0] == 0, 0.2)):
    pool = [k for k in ks if pred(k) and k not in base]
    have = sum(1 for k in base if pred(k))
    need = max(0, int(N * n_min) - have)
    base |= set(random.sample(pool, min(need, len(pool))))
sample = sorted(base)
print(len(sample), "captures", flush=True)

def hat(f):
    W = np.zeros((K, len(f)))
    for kk in range(K):
        lo = knots[kk - 1] if kk else knots[0] - 1; hi = knots[kk + 1] if kk < K - 1 else knots[-1] + 1
        W[kk] = np.clip(np.minimum((f - lo) / (knots[kk] - lo), (hi - f) / (hi - knots[kk])), 0, 1)
    return W

caps = []
for i, key in enumerate(sample):
    b = meta[key]; raw = nl.ir(nl.table()[key]); h = M.sosfiltfilt(ne.BAND, raw, axis=0)
    fr = b["front_shift"]; i0, tr = M.window_start(h, fr)
    end = min(M.usable_end(h), M.K0 + int(fr) + int(0.6 * FS))
    n = end + 200
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = nlt.unwarp(rawp, b["ppm"])[:n]
    m = M.render(key, n, b["o"], fr, band=False) * np.array(b["gain"])
    rec = dict(key=key, trimmed=bool(b["trimmed"]), win={})
    for wi, (t0, t1) in enumerate(WINS):
        a = max(i0, M.K0 + int(fr) + int(t0 * FS)); e = min(end, M.K0 + int(fr) + int(t1 * FS))
        if e - a < 400: continue
        L = e - a; w = np.hanning(L); f = np.fft.rfftfreq(L) * FS
        sel = (f >= 9500) & (f <= 13200)
        Y = np.fft.rfft(y[a:e] * w[:, None], axis=0)[sel]; Mm = np.fft.rfft(m[a:e] * w[:, None], axis=0)[sel]; ff = f[sel]
        W = hat(ff)
        num = np.stack([(W * (Y[:, c] * np.conj(Mm[:, c]))).sum(1) for c in (0, 1)]); den = np.stack([(W * np.abs(Mm[:, c]) ** 2).sum(1) for c in (0, 1)])
        rec["win"][wi] = dict(f=ff, Y=Y, M=Mm, num=num, den=den)
    caps.append(rec)
    if i % 40 == 0: print(i, key, flush=True)
nc = len(caps)
NUM = np.zeros((nc, 4, 2, K), complex); DEN = np.zeros((nc, 4, 2, K)); HAS = np.zeros((nc, 4), bool)
for ci, r in enumerate(caps):
    for wi, d in r["win"].items():
        NUM[ci, wi] = d["num"]; DEN[ci, wi] = d["den"]; HAS[ci, wi] = True
key_of = [r["key"] for r in caps]
attr = dict(
    LATE=np.array([0 if k[2] == 0 else 1 for k in key_of]), SPAN=np.array([0 if k[0] == 0 else 1 for k in key_of]),
    EARLY=np.array([0 if k[1] >= 10 else (1 if k[1] > 0 else 2) for k in key_of]), TRIM=np.array([int(r["trimmed"]) for r in caps]))
LEVELS = dict(LATE=("LATE0", "LATE>0"), SPAN=("SPAN0", "SPAN>0"), EARLY=("EARLY10", "EARLY>0", "EARLY<0"), TRIM=("untrimmed", "trimmed"))

def fit(capidx, wins=(0, 1, 2, 3), chans=(0, 1)):
    n = NUM[np.ix_(capidx, wins, chans)].sum((0, 1, 2)); d = DEN[np.ix_(capidx, wins, chans)].sum((0, 1, 2))
    return n / np.maximum(d, 1e-300)
def gp(g): return np.angle(g) * 1e3, 20 * np.log10(np.abs(g))
rng = np.random.default_rng(1)
def boot(capidx, wins=(0, 1, 2, 3), chans=(0, 1), B=300):
    out = []
    for _ in range(B):
        out.append(fit(rng.choice(capidx, len(capidx)), wins, chans))
    return np.array(out)
lines = []
def P(s=""):
    print(s, flush=True); lines.append(s)
hdr = "knot kHz        " + "  ".join(f"{knots[k]/1000:>12.2f}" for k in SHOW)
groups = []      # (label, capidx, wins, chans)
allc = np.arange(nc)
groups.append(("pooled", allc, (0, 1, 2, 3), (0, 1)))
for wi, wn in enumerate(("t 0-0.05", "t 0.05-0.15", "t 0.15-0.3", "t 0.3-0.6")):
    groups.append((wn, allc[HAS[:, wi]], (wi,), (0, 1)))
for an, lv in LEVELS.items():
    for j, nm in enumerate(lv):
        groups.append((nm, allc[attr[an] == j], (0, 1, 2, 3), (0, 1)))
groups.append(("chan L", allc, (0, 1, 2, 3), (0,))); groups.append(("chan R", allc, (0, 1, 2, 3), (1,)))
P(f"{nc} captures; knots 11.5 ... 12.8 kHz; each cell = phase mrad (bootstrap sd) / gain dB (sd)")
P(hdr)
res = {}
for lab, idx, wins, chans in groups:
    if len(idx) < 5:
        P(f"{lab:14s} (n {len(idx)}) too few"); continue
    g = fit(idx, wins, chans); bt = boot(idx, wins, chans)
    ph, gd = gp(g); pb, gb = gp(bt); res[lab] = (g, bt, len(idx))
    P(f"{lab:14s} n{len(idx):3d} " + "  ".join(f"{ph[k]:+6.1f}({pb[:, k].std():3.1f})/{gd[k]:+5.2f}({gb[:, k].std():.2f})" for k in SHOW))
# significance: each level vs the pooled fit of the others (independent bootstraps)
P(); P("level vs the rest (difference of complex fits): phase diff mrad / gain diff dB, * = |diff| > 3 sd of the bootstrap difference")
pairs = [("t 0-0.05", "t 0.3-0.6"), ("t 0.05-0.15", "t 0.3-0.6"), ("t 0-0.05", "t 0.05-0.15"), ("t 0.15-0.3", "t 0.3-0.6"), ("LATE0", "LATE>0"), ("SPAN0", "SPAN>0"), ("EARLY10", "EARLY>0"), ("EARLY>0", "EARLY<0"), ("EARLY10", "EARLY<0"), ("untrimmed", "trimmed"), ("chan L", "chan R")]
sigs = []
for a_, b_ in pairs:
    if a_ not in res or b_ not in res: continue
    ga, ba, _ = res[a_]; gb_, bb, _ = res[b_]
    dd = ga / gb_; bd = ba / bb[rng.permutation(len(bb))]
    ph = np.angle(dd) * 1e3; gd = 20 * np.log10(np.abs(dd)); pb = np.angle(bd) * 1e3; gbb = 20 * np.log10(np.abs(bd))
    cells = []
    for k in SHOW:
        sp, sg = pb[:, k].std(), gbb[:, k].std()
        f1, f2 = abs(ph[k]) > 3 * sp, abs(gd[k]) > 3 * sg
        cells.append(f"{ph[k]:+5.1f}{'*' if f1 else ' '}/{gd[k]:+5.2f}{'*' if f2 else ' '}")
        if f1 or f2: sigs.append((a_, b_, knots[k], ph[k], gd[k]))
    P(f"{a_:11s} vs {b_:11s} " + "  ".join(cells))
P(); P("significant (3 sd) differences: " + ("; ".join(f"{a} vs {b} @{kz/1000:.2f}k ph {p:+.1f} mrad gain {g:+.2f} dB" for a, b, kz, p, g in sigs) if sigs else "none"))

# --- held out: group-specific vs pooled, capture halves, scored in the band bins
def E_of(g):
    kx = np.r_[9500.0, knots]; return lambda f: np.interp(f, kx, np.r_[1.0, g.real]) + 1j * np.interp(f, kx, np.r_[0.0, g.imag])
def score(ci, wi, c, Efn, bands):
    d = caps[ci]["win"][wi]; f = d["f"]; Yc = d["Y"][:, c]; Mc = Efn(f) * d["M"][:, c]
    out = {}
    for bn, (lo, hi) in bands.items():
        s = (f >= lo) & (f < hi); out[bn] = (np.sum(np.abs(Yc[s] - Mc[s]) ** 2), np.sum(np.abs(Yc[s]) ** 2))
    return out
VARIANTS = ["none", "pooled", "time", "LATE", "SPAN", "EARLY", "TRIM", "chan", "time+LATE", "time+chan"]
def gkey(var, ci, wi, c):
    return {"none": None, "pooled": (0,), "time": (wi,), "LATE": (attr["LATE"][ci],), "SPAN": (attr["SPAN"][ci],), "EARLY": (attr["EARLY"][ci],), "TRIM": (attr["TRIM"][ci],),
            "chan": (c,), "time+LATE": (wi, attr["LATE"][ci]), "time+chan": (wi, c)}[var]
def fit_variant(var, train):
    gs = {}
    keys = {gkey(var, ci, wi, c) for ci in train for wi in range(4) if HAS[ci, wi] for c in (0, 1)} if var != "none" else set()
    for gk in keys:
        n = 0; d = 0
        for ci in train:
            for wi in range(4):
                if not HAS[ci, wi]: continue
                for c in (0, 1):
                    if gkey(var, ci, wi, c) == gk: n = n + NUM[ci, wi, c]; d = d + DEN[ci, wi, c]
        gs[gk] = n / np.maximum(d, 1e-300)
    return gs
folds = [(allc[0::2], allc[1::2]), (allc[1::2], allc[0::2])]
acc = {v: {bn: [] for bn in BANDS} for v in VARIANTS}      # per test capture-channel nulls
who = {v: [] for v in VARIANTS}
for train, test in folds:
    fits = {v: fit_variant(v, train) for v in VARIANTS}
    for ci in test:
        for c in (0, 1):
            tot = {v: {bn: [0.0, 0.0] for bn in BANDS} for v in VARIANTS}
            for wi in range(4):
                if not HAS[ci, wi]: continue
                for v in VARIANTS:
                    gk = gkey(v, ci, wi, c)
                    Efn = (lambda f: np.ones(len(f), complex)) if v == "none" else E_of(fits[v].get(gk, fits["pooled"][(0,)]))
                    s = score(ci, wi, c, Efn, BANDS)
                    for bn in BANDS: tot[v][bn][0] += s[bn][0]; tot[v][bn][1] += s[bn][1]
            for v in VARIANTS:
                for bn in BANDS: acc[v][bn].append(10 * np.log10(max(tot[v][bn][0], 1e-300) / max(tot[v][bn][1], 1e-300)))
P(); P("held out by capture halves (2 folds), median null dB per capture-channel (error energy over the Y energy in the band, windows summed), "
       "params per group: 22 (11 complex knots)")
P(f"{'variant':10s} {'groups':>6s}  " + "  ".join(f"{bn:>6s} med" for bn in BANDS) + "   | paired vs pooled (C band): mean dB, worse>0.5 dB, better>0.5 dB")
ng = {"none": 0, "pooled": 1, "time": 4, "LATE": 2, "SPAN": 2, "EARLY": 3, "TRIM": 2, "chan": 2, "time+LATE": 8, "time+chan": 8}
for v in VARIANTS:
    d = np.array(acc[v]["C"]) - np.array(acc["pooled"]["C"])
    P(f"{v:10s} {ng[v]:6d}  " + "  ".join(f"{np.median(acc[v][bn]):10.1f}" for bn in BANDS) + f"   | {d.mean():+.2f}, {int((d > 0.5).sum())}, {int((d < -0.5).sum())} of {len(d)}")
open(f"/Users/leighpierce/Git/246/analysis/nl_hfdrift/b_groups_{TAG}.txt", "w").write("\n".join(lines))

# --- per capture-channel: is the 0-0.05 s deviation at 12.5-12.8 kHz systematic or a few captures?
P(); P("per capture-channel complex gain at 12.25 / 12.5 / 12.8 kHz, window 0-0.05 s vs 0.05-0.6 s (median [25%-75%] phase mrad; gain dB); only captures with HF energy >= 1% of the median")
def pc(sel_caps, wins, k):
    ph, gd = [], []
    for ci in sel_caps:
        for c in (0, 1):
            n = NUM[ci, wins, c, k].sum(); d = DEN[ci, wins, c, k].sum()
            if d > 0: ph.append((np.angle(n / d) * 1e3, d)); gd.append(20 * np.log10(abs(n / d)))
    dd = np.array([x[1] for x in ph]); keep = dd >= 0.01 * np.median(dd)
    p = np.array([x[0] for x in ph])[keep]; g = np.array(gd)[keep]
    return p, g
for gname, gidx in (("all", allc), ("EARLY10", allc[attr["EARLY"] == 0]), ("EARLY>0", allc[attr["EARLY"] == 1]), ("EARLY<0", allc[attr["EARLY"] == 2])):
    for wl, wins in (("0-0.05", [0]), ("0.05-0.6", [1, 2, 3])):
        cells = []
        for k in (7, 8, 9):
            p, g = pc(gidx, wins, k)
            cells.append(f"{np.median(p):+6.1f} [{np.percentile(p, 25):+6.1f},{np.percentile(p, 75):+6.1f}] / {np.median(g):+5.2f}")
        P(f"{gname:8s} {wl:9s} n{len(p):4d}  " + "   ".join(cells))
# fraction of the 12.5-12.8 kHz energy that sits in window 0
e0 = DEN[:, 0, :, 8].sum((1)); eall = DEN[:, :, :, 8].sum((1, 2))
for j, nm in enumerate(LEVELS["EARLY"]):
    s = attr["EARLY"] == j; P(f"{nm}: share of 12.5 kHz model energy in 0-0.05 s: median {np.median(e0[s] / eall[s]):.2f}")
open(f"/Users/leighpierce/Git/246/analysis/nl_hfdrift/b_groups_{TAG}.txt", "w").write("\n".join(lines))

# --- extra: is what is left after the pooled E a per-capture alignment error? (in-sample, windows 1-3, one fractional delay per capture-channel)
P(); P("after the pooled E: best per-capture-channel fractional delay (samples) on 11.6-12.8 kHz, windows 0.05-0.6 s, in-sample")
gp_ = fit(allc); Ep = E_of(gp_)
taus = np.arange(-0.03, 0.0301, 0.0005)
rows = []
for ci in allc:
    for c in (0, 1):
        errs = np.zeros(len(taus)); ref = 0.0
        for wi in (1, 2, 3):
            if not HAS[ci, wi]: continue
            d = caps[ci]["win"][wi]; f = d["f"]; s = (f >= 11600) & (f < 12800)
            Mc = (Ep(f) * d["M"][:, c])[s]; Yc = d["Y"][:, c][s]; fs_ = f[s] / FS
            for j, t in enumerate(taus):
                errs[j] += np.sum(np.abs(Yc - Mc * np.exp(-2j * np.pi * fs_ * t)) ** 2)
            ref += np.sum(np.abs(Yc) ** 2)
        if ref > 0:
            j0 = int(np.argmin(np.abs(taus)))
            rows.append((ci, c, taus[int(np.argmin(errs))], 10 * np.log10(errs[j0] / ref), 10 * np.log10(errs.min() / ref)))
R_ = np.array([(r[2], r[3], r[4]) for r in rows])
P(f"{len(rows)} capture-channels: median null before {np.median(R_[:,1]):.1f} dB, after the per-capture delay {np.median(R_[:,2]):.1f} dB; median |tau| {np.median(np.abs(R_[:,0])):.4f}, IQR of tau [{np.percentile(R_[:,0],25):+.4f}, {np.percentile(R_[:,0],75):+.4f}]")
tl = {(r[0]): {} for r in rows}
for r in rows: tl[r[0]][r[1]] = r[2]
both = np.array([(v[0], v[1]) for v in tl.values() if 0 in v and 1 in v])
P(f"correlation of tau between L and R of the same capture: {np.corrcoef(both[:,0], both[:,1])[0,1]:+.2f} (n {len(both)})")
for j, nm in enumerate(LEVELS["EARLY"]):
    sel = np.array([attr["EARLY"][r[0]] == j for r in rows]); P(f"  {nm}: median tau {np.median(R_[sel,0]):+.4f}, null before {np.median(R_[sel,1]):.1f} after {np.median(R_[sel,2]):.1f}")
open(f"/Users/leighpierce/Git/246/analysis/nl_hfdrift/b_groups_{TAG}.txt", "w").write("\n".join(lines))
