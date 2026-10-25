"""Plugin (UAD EMT 250) versus the Altiverb REV-250 presets: which plugin LF/HF settings match each preset's decay, and how the
early cluster, spectrum and loops compare.  .venv/bin/python tools/uad250/compare_altiverb.py

Altiverb: analysis/altiverb/cache/REV-250 <RT>s.npy (deconvolved, 48 kHz, 2 x N). Plugin: A1 impulse grid (Noise off).
Both are single realisations of a time-varying reverb, so only statistics and the static parts (early cluster, decay, spectrum)
are compared; tails cannot be nulled (correlation of the diffuse tails is ~0 by construction, see UAD250_test_plan.md 2.2).
"""
import csv, os, sys
import numpy as np, soundfile as sf
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "tools", "harness"))
from harness.core import band_rt, find_loops  # noqa: E402
r = 48000; OCT = [125, 250, 500, 1000, 2000, 4000, 8000]
rows = list(csv.DictReader(open(os.path.join(ROOT, "analysis", "uad250", "a1", "a1.csv"))))
RTS = [0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5]

def first(x, th=0.02): return int(np.argmax(np.abs(x) > th * np.abs(x).max()))
def t60(h): return [band_rt(h, r, f / 1.414, f * 1.414)["T"] for f in OCT]

print("RT   best plugin (LF, HF)  rms ln-error  | Altiverb T60 L (125..8k) | plugin T60 L at best | first arrival ms (Alt L/R | plug L/R)")
best = {}
for RT in RTS:
    h = np.load(os.path.join(ROOT, "analysis", "altiverb", "cache", f"REV-250 {RT}s.npy"))[:, :int(min(3 * RT + 2, 12) * r)]
    f0, f1 = first(h[0]), first(h[1]); R_i = 0 if f0 < f1 else 1  # which channel leads (R)
    tl = [t60(h[1 - R_i]), t60(h[R_i])]  # [L, R]
    cand = []
    for d in rows:
        if float(d["T"]) != RT: continue
        pl = [float(d[f"T60_{f}"]) for f in OCT]; pr = [float(d[f"T60R_{f}"]) for f in OCT]
        err = np.sqrt(np.mean([np.log(np.array(tl[0]) / np.array(pl)) ** 2, np.log(np.array(tl[1]) / np.array(pr)) ** 2]))
        cand.append((err, d))
    err, d = min(cand, key=lambda c: c[0]); best[RT] = (d, h, R_i, tl)
    print(f"{RT:3.1f}  LF x{d['LF']}, HF {d['HF']:>4s}   {err:.3f}  | " + " ".join(f"{v:4.2f}" for v in tl[0]) + " | " + " ".join(f"{float(d[f'T60_{f}']):4.2f}" for f in OCT) +
          f" | {min(f0, f1) / r * 1e3:5.2f}/{max(f0, f1) / r * 1e3:5.2f} | {float(d['first_R_ms']) + 11.7:.2f}/{float(d['first_L_ms']) + 11.7:.2f}")

print("\nEarly cluster (0-60 ms after the R onset), best-matching plugin setting: correlation, level step and spectral difference")
for RT in RTS:
    d, h, R_i, tl = best[RT]
    p, _ = sf.read(os.path.join(ROOT, "analysis", "uad250", "captures", "A_static", f"imp_T{RT}_LF{d['LF']}_HF{d['HF']}.wav"))
    out = []
    for nm, ca, cp in (("R", R_i, 1), ("L", 1 - R_i, 0)):
        a = h[ca]; b = p[:, cp].astype(float); ia, ib = first(a), first(b)
        n = int(0.060 * r); x = a[ia - 5:ia - 5 + n]; y = b[ib - 5:ib - 5 + n]
        c = (x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum()); g = np.sqrt((y * y).sum() / (x * x).sum())
        Fx = np.abs(np.fft.rfft(x, 8192)) ** 2; Fy = np.abs(np.fft.rfft(y, 8192)) ** 2; fr = np.fft.rfftfreq(8192, 1 / r)
        sp = [10 * np.log10(Fy[(fr >= f / 1.414) & (fr < f * 1.414)].sum() / Fx[(fr >= f / 1.414) & (fr < f * 1.414)].sum()) for f in (250, 1000, 4000, 8000)]
        out.append(f"{nm}: corr {c:.3f}, level {20 * np.log10(g):+.1f} dB, spectrum plug-Alt at 250/1k/4k/8k Hz " + "/".join(f"{v:+.1f}" for v in sp))
    print(f"RT {RT}: " + " | ".join(out))

print("\nLoops (decay-flattened ACF, 8 sigma, 48 kHz) Altiverb versus plugin at RT >= 3.0:")
for RT in (3.0, 3.6, 4.0, 4.5):
    d, h, R_i, tl = best[RT]
    L, _ = find_loops(h.T[:int(10 * r)], int(0.3 * r), int(min(len(h[0]), (0.3 + 1.2 * RT) * r)), RT, float(r), n_sigma=6.0, max_lag=int(0.4 * r))
    print(f"  RT {RT}: Altiverb loops ms " + ", ".join(f"{p / r * 1e3:.2f}" for p in L) + f" | plugin {d['loops_ms'] or '-'}")
