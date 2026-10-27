"""Fork E: is the bench's per-capture alignment (branch offset o, front shift, clock ratio) too coarse at 12 kHz?
Held-out refits of (o, fr, ppm) on one band, scored on others and on the headline 100 Hz-10 kHz windows (strict rule),
with the pooled HF correction (R front-tap estimate from fork C, on both channels) OFF and ON; and a front-part vs
branch-part delay test per band.  usage: e_align.py N [seed]   (no shared file touched; bench_warp2 fits read-only)"""
import os, sys, json, random, pickle, time
os.environ["NL_WARP2"] = "1"; os.environ.pop("NL_BK", None)
import numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
import common
from common import *
from scipy.signal import butter, sosfiltfilt
from scipy.optimize import minimize
import c_front as CF
import defects
common.BENCH = "/Users/leighpierce/Git/246/analysis/nl/bench_warp2"
N_CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 10
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 41
KNOTS = CF.KNOTS
front = pickle.load(open(OUT + "/c_front_150.pkl", "rb"))["all"]
fkeys = {tuple(r["key"]) for r in front}
corr = CF.pooled(front, 1)                       # the R front-tap estimate, applied to both channels
noise = json.load(open(M.NOISE))
WINDOWS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.6)]
BAND_NAMES = dict(B20=(20, 100), B100=(100, 1000), B1k=(1000, 5000), B5k=(5000, 10000), A=(9000, 11000), S1=(11000, 11600),
                  S2=(11600, 12300), S3=(12300, 12800), B=(11000, 12800), C=(10000, 12000))
SOS = {nm: butter(4, [lo, hi], btype="band", fs=FS, output="sos") for nm, (lo, hi) in BAND_NAMES.items()}
FIT_BANDS = dict(A=(9000, 11000), B=(11000, 12800), LO=(5000, 10000), AB=(9000, 12800))


def corr_gain(f):
    kx = np.r_[9500.0, KNOTS[:-1]]; g = corr[:-1]
    return np.interp(f, kx, np.r_[1.0, g.real], right=1.0) + 1j * np.interp(f, kx, np.r_[0.0, g.imag], right=0.0)


class Cap:
    def __init__(self, key, on):
        b = common.bench(key); raw = nl.ir(nl.table()[key]); h = sosfiltfilt(ne.BAND, raw, axis=0)
        self.key, self.b, self.on = key, b, on
        self.fr0, self.o0, self.ppm0 = b["front_shift"], b["o"], b["ppm"]
        i0, self.trimmed = M.window_start(h, self.fr0)
        self.a = max(i0, M.K0 + int(self.fr0) + int(0.05 * FS)); self.b2 = min(M.usable_end(h), M.K0 + int(self.fr0) + int(0.6 * FS))
        self.i0 = i0
        self.ok = self.b2 - self.a >= int(0.05 * FS)
        if not self.ok: return
        self.n = self.b2 + 200
        self.rawp = np.vstack([raw[:self.n + 64], np.zeros((max(0, self.n + 64 - len(raw)), 2))])
        self.up32 = nlt.upsample(self.rawp)
        self.NN = 1 << int(np.ceil(np.log2(self.n + 2 * ne.TAPS[-1] + 4096)))
        self.ph = -2j * np.pi * np.fft.rfftfreq(self.NN)
        Rp = M._proc(self.NN)
        if on: Rp = Rp * corr_gain(np.fft.rfftfreq(self.NN) * FS)[:, None]
        self.spec = [(M.spectra(key, ch, self.NN), Rp[:, ch]) for ch in (0, 1)]
        self.af = max(self.i0, 0)                        # the fit window starts at the content start, so the front taps pin the front shift (o and fr are degenerate on the branch alone)
        self.L = self.b2 - self.af; self.w = np.hanning(self.L); self.f = np.fft.rfftfreq(self.L) * FS

    def model(self, do=0.0, dfr=0.0, parts=("F", "B")):
        out = []
        for (F, B), R in self.spec:
            S = (F if "F" in parts else 0) + (B * np.exp(self.ph * (self.o0 + do)) if "B" in parts else 0)
            out.append(np.fft.irfft(S * np.exp(self.ph * (M.K0 + self.fr0 + dfr)) * R, self.NN)[:self.n])
        return np.array(out).T

    def data(self, dppm=0.0):
        return nlt.unwarp_from(self.up32, self.n, self.ppm0 + dppm)[:self.n]

    def cost(self, p, band):
        do, dfr, dppm = p
        y = self.data(dppm)[self.af:self.b2]; m = self.model(do, dfr)[self.af:self.b2]
        sel = (self.f >= band[0]) & (self.f < band[1])
        tot = 0.0
        for c in (0, 1):
            Y = np.fft.rfft(y[:, c] * self.w)[sel]; Mm = np.fft.rfft(m[:, c] * self.w)[sel]
            yy = np.sum(np.abs(Y) ** 2); mm = np.sum(np.abs(Mm) ** 2); ym = np.sum((Y * np.conj(Mm)).real)
            tot += 10 * np.log10(max(yy - ym * ym / mm, 1e-30) / yy)
        return tot / 2

    def fit(self, band, free):
        """Nelder-Mead on the free subset of (do, dfr, dppm) ; returns the full 3-vector."""
        idx = [i for i, v in enumerate(("o", "fr", "ppm")) if v in free]
        steps = np.array([0.01, 0.01, 0.3])[idx]
        def f(x):
            p = np.zeros(3); p[idx] = x
            return self.cost(p, band)
        x0 = np.zeros(len(idx)); sim = np.vstack([x0] + [x0 + np.eye(len(idx))[i] * steps[i] for i in range(len(idx))])
        r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=2e-4, fatol=1e-3, maxiter=120, initial_simplex=sim))
        p = np.zeros(3); p[idx] = r.x
        return p, r.fun

    def score(self, p):
        """bands (C-style: whole arrays filtered, then cut) and headline windows (bench style), gain refit on 100 Hz-10 kHz."""
        do, dfr, dppm = p
        y = self.data(dppm); m = self.model(do, dfr)
        yb, mb = sosfiltfilt(ne.BAND, y, axis=0), sosfiltfilt(ne.BAND, m, axis=0)
        g = np.array([(yb[self.a:self.b2, c] @ mb[self.a:self.b2, c]) / (mb[self.a:self.b2, c] @ mb[self.a:self.b2, c]) for c in (0, 1)])
        out = dict(g=g.tolist(), bands={}, win={})
        for nm in BAND_NAMES:
            ya, ma = sosfiltfilt(SOS[nm], y, axis=0)[self.a:self.b2], sosfiltfilt(SOS[nm], m * g, axis=0)[self.a:self.b2]
            out["bands"][nm] = [float(10 * np.log10(np.sum((ya[:, c] - ma[:, c]) ** 2) / np.sum(ya[:, c] ** 2))) for c in (0, 1)]
        for t0, t1 in WINDOWS:
            a = max(self.i0, M.K0 + int(self.fr0) + int(t0 * FS)); b = min(self.b2, M.K0 + int(self.fr0) + int(t1 * FS))
            if b - a < int(0.02 * FS) or np.sum(yb[a:b, 0] ** 2) < 10 * noise["ch0"] * (b - a): continue
            out["win"][f"{t0}-{t1}"] = [float(10 * np.log10(np.sum((yb[a:b, c] - mb[a:b, c] * g[c]) ** 2) / np.sum(yb[a:b, c] ** 2))) for c in (0, 1)]
        return out


def med_iqr(v):
    v = np.asarray(v, float); return "%+.4f [%+.4f, %+.4f]" % (np.median(v), np.percentile(v, 25), np.percentile(v, 75))


def part_delays(caps):
    """front-part window (front taps only) and branch-part window (0.05-0.6 s): pooled cross-spectrum phase delay per band, L / R."""
    NF = 4096; fq = np.fft.rfftfreq(NF) * FS
    bands = dict(LO=(5000, 10000), MID=(10000, 11600), HI=(11600, 12800))
    acc = {(w, c): np.zeros(len(fq), complex) for w in ("front", "branch") for c in (0, 1)}
    per = {(w, c, bn): [] for w in ("front", "branch") for c in (0, 1) for bn in bands}
    for cap in caps:
        y = cap.data(); m = cap.model()
        for c in (0, 1):
            s0 = int(M.K0 + cap.fr0 + CF.B0.get(c, 685) * 0 + 41 - 30); s1 = int(M.K0 + cap.fr0 + CF.B0[c] - CF.SPAN_BEFORE_BRANCH)
            for w, (a, b) in (("front", (s0, s1)), ("branch", (cap.a, cap.b2))):
                hw = np.hanning(b - a)
                X = np.fft.rfft(y[a:b, c] * hw, NF) * np.conj(np.fft.rfft(m[a:b, c] * hw, NF))
                acc[(w, c)] += X
                for bn, (lo, hi) in bands.items():
                    k = (fq >= lo) & (fq < hi); ph = np.angle(X[k]); wt = np.abs(X[k])
                    if wt.sum() > 0: per[(w, c, bn)].append(float(-np.sum(wt * ph * fq[k]) / (2 * np.pi * np.sum(wt * fq[k] ** 2)) * FS))
    out = {}
    for (w, c), X in acc.items():
        for bn, (lo, hi) in bands.items():
            k = (fq >= lo) & (fq < hi); ph = np.angle(X[k]); wt = np.abs(X[k])
            out[(w, c, bn)] = (float(-np.sum(wt * ph * fq[k]) / (2 * np.pi * np.sum(wt * fq[k] ** 2)) * FS), per[(w, c, bn)])
    return out


if __name__ == "__main__":
    ks = sorted(k for k in nl.table() if k not in defects.NL_DEFECTS and k not in fkeys and os.path.exists(common.BENCH + "/" + "_".join(str(v) for v in k) + ".json"))
    random.seed(SEED); sample = random.sample(ks, N_CAP)
    res = {}
    for on in (False, True):
        t0 = time.time(); tag = "ON" if on else "OFF"
        caps = [c for c in (Cap(k, on) for k in sample) if c.ok]
        print(f"=== pooled HF correction {tag}: {len(caps)} captures, sample seed {SEED}, built in {time.time()-t0:.0f} s", flush=True)
        variants = {"base": lambda c: (np.zeros(3), None),
                    "fit A(9-11) o,fr,ppm": lambda c: c.fit(FIT_BANDS["A"], ("o", "fr", "ppm")),
                    "fit B(11-12.8) o,fr,ppm": lambda c: c.fit(FIT_BANDS["B"], ("o", "fr", "ppm")),
                    "fit AB(9-12.8) o,fr,ppm": lambda c: c.fit(FIT_BANDS["AB"], ("o", "fr", "ppm")),
                    "fit LO(5-10) o,fr,ppm": lambda c: c.fit(FIT_BANDS["LO"], ("o", "fr", "ppm")),
                    "fit LO(5-10) o,fr": lambda c: c.fit(FIT_BANDS["LO"], ("o", "fr")),
                    "fit LO(5-10) o": lambda c: c.fit(FIT_BANDS["LO"], ("o",)),
                    "fit A(9-11) o,fr": lambda c: c.fit(FIT_BANDS["A"], ("o", "fr")),
                    "fit B(11-12.8) o,fr": lambda c: c.fit(FIT_BANDS["B"], ("o", "fr"))}
        sc = {v: [] for v in variants}; shifts = {v: [] for v in variants}
        for i, cap in enumerate(caps):
            for v, fn in variants.items():
                p, _ = fn(cap); shifts[v].append(p); sc[v].append(cap.score(p))
            if i % 10 == 0: print(f"  capture {i+1}/{len(caps)} {cap.key} ({time.time()-t0:.0f} s)", flush=True)
        base = sc["base"]
        print(f"--- shifts (samples; ppm in ppm), median [IQR]  [{tag}]")
        for v in variants:
            if v == "base": continue
            s = np.array(shifts[v]); print(f"   {v:28s} do {med_iqr(s[:,0])}   dfr {med_iqr(s[:,1])}   dppm {med_iqr(s[:,2])}")
        print(f"--- held-out scoring [{tag}]: median null L/R (dB) per band; diff vs base median, count worse >0.5 dB of 2x{len(caps)}")
        for v in variants:
            line = []
            for nm in ("B5k", "A", "S1", "S2", "S3", "B", "C"):
                x = np.array([s["bands"][nm] for s in sc[v]]); b = np.array([s["bands"][nm] for s in base]); d = x - b
                line.append(f"{nm} {np.median(x[:,0]):.1f}/{np.median(x[:,1]):.1f} ({np.median(d):+.2f}; w{int((d>0.5).sum())})")
            print(f"   {v:28s} " + " | ".join(line))
        print(f"--- strict rule, headline windows 100 Hz-10 kHz [{tag}]: channel-windows worse than base by >0.5 dB; worst; median diff")
        for v in variants:
            if v == "base": continue
            line = []
            for wn in ("0.0-0.05", "0.05-0.15", "0.15-0.3", "0.3-0.6"):
                d = [s["win"][wn][c] - b["win"][wn][c] for s, b in zip(sc[v], base) if wn in s["win"] and wn in b["win"] for c in (0, 1)]
                if d: line.append(f"{wn}: {int((np.array(d)>0.5).sum())}/{len(d)} worst {max(d):+.2f} med {np.median(d):+.2f}")
            print(f"   {v:28s} " + " | ".join(line))
        pd = part_delays(caps)
        print(f"--- front part vs branch part delay of capture re model (samples; positive = capture later), pooled [per-capture median, IQR] [{tag}]")
        for bn in ("LO", "MID", "HI"):
            for c in (0, 1):
                f_, fl = pd[("front", c, bn)]; b_, bl = pd[("branch", c, bn)]
                print(f"   {bn:3s} ch{c}: front {f_:+.4f} [{med_iqr(fl)}]   branch {b_:+.4f} [{med_iqr(bl)}]   branch-front {b_-f_:+.4f}")
        res[tag] = dict(shifts={v: np.array(s).tolist() for v, s in shifts.items()}, sc=sc, keys=[list(c.key) for c in caps])
    json.dump(res, open(OUT + f"/e_align_{N_CAP}.json", "w"), default=float)
