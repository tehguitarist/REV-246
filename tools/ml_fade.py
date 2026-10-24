"""NEXT_STEPS §0 item 1: Tomkin's fade-out in programs M and L, for ml_lines.usable_end. Results: docs/ml_structure.md §15.

    EMT_KERNEL=chain3 .venv/bin/python tools/ml_fade.py --stage measure [--prog ML] [--only N] [--keys P_key,..] [--workers 8]
    EMT_KERNEL=chain3 .venv/bin/python tools/ml_fade.py --stage report
    EMT_KERNEL=chain3 .venv/bin/python tools/ml_fade.py --stage rt [--workers 6]     # then --stage rtreport

  measure  per AMPL group (ml_lines.fitted: the static model with each capture's stored clock and gains, nothing refitted):
           octave-band energies (250 Hz-4 kHz) of capture and model, per channel, in FRAME-sample frames on the K0 grid
           from the R onset + 0.1 s to the capture's end (the last sample above 1e-6 of its peak: the raw file's end to
           one device sample). One JSON per group in analysis/ml_fade/curves/ (resumable).
  report   per file and channel, the capture / model level difference d(t) (mean over bands, frames above the floor)
           fitted as a + b t + 10 log10 <g(t)^2>_frame, g a fade ending at the file's end (a and b absorb the gain and any
           small decay error of the static model, which the M / L time variance makes up to ~1 %); g's shape (linear,
           quarter-sine, raised cosine, dB-linear) and length T by grid search, the end offset too in a check. Then T
           against program, file length, RT and AMPL; a check against partner captures (energy ratio of a file to a
           longer AMPL partner, no model) -> analysis/ml_fade/fade.json (read by ml_lines.usable_end), report.log.
  rt       ml_structure §8.3's per-octave RT / law re-read: ml_lines.band_decay from the onset to the old usable end
           (0.72 s before the file's end) and to the new one, per group in analysis/ml_fade/rt/ -> rt.log.
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.signal import sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import ml_lines as mll  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_fade")
FS, K0 = mll.FS, mll.K0
FRAME = 660                                                 # 25 ms
BANDS = [250, 500, 1000, 2000, 4000]
T0 = 0.1                                                    # frames start this long after the R onset (s)


# ---------------------------------------------------------------- measure
def end_of(y):
    """Per channel, one past the last sample above 1e-6 of the peak (the raw file's end, ml_structure's grid)."""
    a = np.abs(y)
    return [int(np.nonzero(a[:, ch] > a[:, ch].max() * 1e-6)[0][-1]) + 1 for ch in (0, 1)]


def band_frames(x, k0, k1):
    """dB energy per FRAME frame [k0 + j FRAME, ...) up to k1, per band and channel: (bands, 2, frames)."""
    nf = (k1 - k0) // FRAME
    out = np.empty((len(BANDS), 2, nf))
    for i, f in enumerate(BANDS):
        xb = sosfiltfilt(mll.OCT_SOS[f], x[:k1], axis=0)[k0:k0 + nf * FRAME]
        out[i] = 10 * np.log10(np.sum(xb.reshape(nf, FRAME, 2) ** 2, 1).T + 1e-30)
    return out


def _measure_job(a):
    prog, key, path = a
    fit = mll.fitted(prog, key)
    k0 = K0 + mll.cc.ONSET[prog][1] + int(T0 * FS)
    rows = []
    for amp, (y, mdl, _) in fit.items():
        e = end_of(y)
        k1 = max(e)
        Y, M = band_frames(y, k0, k1), band_frames(mdl, k0, k1)
        rows.append(dict(amp=amp, end=e, Y=np.round(Y, 2).tolist(), M=np.round(M, 2).tolist()))
    lf, rt, hf = (float(x) for x in key.split("_"))
    json.dump(dict(prog=prog, key=key, lf=lf, rt=rt, hf=hf, k0=k0, frame=FRAME, bands=BANDS, rows=rows), open(path, "w"))
    print(f"{prog} {key}: {len(rows)} captures, ends " + " ".join(f"{(min(r['end']) - K0) / FS:.2f}" for r in rows) + " s", flush=True)
    return key


def measure_stage(progs, workers, only=None, keys=None):
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, "curves")
    os.makedirs(d, exist_ok=True)
    jobs = []
    for p in progs:
        ks = sorted(mll.cc.groups(p))
        if keys:
            ks = [k for q, k in keys if q == p]
        elif only:
            ks = ks[::max(1, len(ks) // only)][:only]
        jobs += [(p, k, os.path.join(d, f"{p}_{k}.json")) for k in ks if not os.path.exists(os.path.join(d, f"{p}_{k}.json"))]
    print(f"{len(jobs)} groups to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_measure_job, jobs):
            pass


# ---------------------------------------------------------------- report
SHAPES = {"linear": lambda u: u, "quarter-sine": lambda u: np.sin(np.pi / 2 * u),
          "raised-cosine": lambda u: 0.5 - 0.5 * np.cos(np.pi * u), "dB-linear 60": lambda u: 10 ** (-3 * (1 - u)) * (u > 0)}
FLOOR_DB = 90.0                 # a band's frame counts while the model is within this of its first frame
DROP_DB = 20.0                  # ... and while d is within this of the file's early level (S's RAMP_MIN, 0.1)
SUB = 16                        # sub-points per frame for <g^2>
TS = np.round(np.concatenate([np.arange(0.2, 2.0, 0.02), np.arange(2.0, 8.0, 0.05), np.arange(8.0, 30.0, 0.25)]), 3)


def fade_db(tl, T, g, off=0.0):
    """10 log10 <g(u)^2> over each frame; tl (frames, SUB) = time from each sub-point to the file's end (s)."""
    u = np.clip((tl - off) / T, 0.0, 1.0)
    return 10 * np.log10(np.mean(g(u) ** 2, 1) + 1e-12)


def curves(prog=None):
    for p in sorted(glob.glob(os.path.join(OUT, "curves", "*.json"))):
        D = json.load(open(p))
        if prog is None or D["prog"] == prog:
            yield D


def file_d(D, r, ch, ref=None):
    """d(t) per frame for one capture and channel against the model (or ref, another capture's row, on the same frames):
    sub-point times to the file's end tl (frames, SUB) and d (dB, mean over bands), frames above the floor only."""
    Y, M = np.array(r["Y"])[:, ch], np.array(r["M"])[:, ch]
    e = r["end"][ch]
    nf = (e - D["k0"]) // D["frame"]
    if ref is not None:
        nf = min(nf, (ref["end"][ch] - D["k0"]) // D["frame"])
    Y, M = Y[:, :nf], M[:, :nf]
    Z = M if ref is None else np.array(ref["Y"])[:, ch, :nf]
    keep = ((M > M[:, :1] - FLOOR_DB) & (Y > M[:, :1] - FLOOR_DB - 20)).all(0)
    d = np.mean(Y - Z, 0)
    early = np.median(d[:max(4, nf // 4)])
    keep &= d > early - DROP_DB
    starts = D["k0"] + D["frame"] * np.arange(nf)[keep]
    tl = (e - (starts[:, None] + (np.arange(SUB) + 0.5) / SUB * D["frame"])) / FS
    return tl, d[keep]


def fit_fade(tl, d, g, Ts=TS, off=0.0, sb=None):
    """Grid over T: d - fade = a (sb None) or a + b t with b ~ N(0, sb^2) against 0.5 dB frame noise; (T, rms, a, b)."""
    t = -tl.mean(1)
    best = None
    for T in Ts:
        r = d - fade_db(tl, T, g, off)
        if sb is None:
            c = np.array([r.mean(), 0.0])
        else:
            X = np.vstack([np.stack([np.ones_like(t), t], 1), [0.0, 0.5 / sb]])
            c = np.linalg.lstsq(X, np.append(r, 0.0), rcond=None)[0]
        rms = float(np.sqrt(np.mean((r - c[0] - c[1] * t) ** 2)))
        if best is None or rms < best[1]:
            best = (float(T), rms, float(c[0]), float(c[1]))
    return best


def fit_power(tl, d, T0):
    """g = u^p: (T, p, rms) from the linear fit's T."""
    from scipy.optimize import minimize

    def cost(v):
        r = d - fade_db(tl, np.exp(v[0]), lambda u: u ** v[1])
        return float(np.sqrt(np.mean((r - r.mean()) ** 2)))
    s = minimize(cost, [np.log(T0), 1.0], method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5))
    return float(np.exp(s.x[0])), float(s.x[1]), float(s.fun)


def fit_start(tl, d, T0, length):
    """Linear ramp with its start (s before the file's end, inside the file) and its length T (<= 50 s) both free; the
    ramp's zero point off = start - T (< 0: the file was cut inside the ramp). On a partial ramp T and off trade, the
    start does not. Returns (T, off, rms)."""
    from scipy.optimize import minimize

    def cost(v):
        r = d - fade_db(tl, v[1], SHAPES["linear"], v[0] - v[1])
        return float(np.sqrt(np.mean((r - r.mean()) ** 2)))
    bounds = [(0.05, length), (0.05, 50.0)]
    s = min((minimize(cost, [float(np.clip(st, 0.06, length * 0.99)), float(np.clip(t0, 0.06, 49.0))], method="Powell", bounds=bounds, options=dict(xtol=1e-4, ftol=1e-7))
             for st, t0 in ((T0, T0), (1.5 * T0, 1.5 * T0), (1.0, 1.0), (1.5 * T0, 3 * T0))), key=lambda s: s.fun)
    return float(s.x[1]), float(s.x[0] - s.x[1]), float(s.fun)


def _fit_group(D):
    import ml_carryover as cc
    names = cc.groups(D["prog"])[D["key"]]
    rows = []
    for r in D["rows"]:
        for ch in (0, 1):
            tl, d = file_d(D, r, ch)
            if len(d) < 20:
                continue
            f = {s: fit_fade(tl, d, g) for s, g in SHAPES.items()}
            T1 = f["linear"][0]
            length = (r["end"][ch] - K0) / FS
            T, off, rms = fit_start(tl, d, T1, length)
            row = dict(prog=D["prog"], key=D["key"], amp=r["amp"], ch=ch, rt=D["rt"], lf=D["lf"], hf=D["hf"],
                       file=os.path.basename(names[r["amp"]]), length=length, fits=f,
                       ridge=fit_fade(tl, d, SHAPES["linear"], sb=0.01 * 60 / D["rt"]), power=fit_power(tl, d, T1),
                       T=T, off=off, rms=rms, start=T + off, tmin=float(tl.min()), rms0=float(np.std(d)),
                       # the fade is seen when its start lies inside the frames above the floor and it lowers the
                       # residual; otherwise the reverb met the floor first, and the file is fade-free up to tmin
                       observable=bool(T1 > 0.2 and T + off > tl.min() + 0.1 and np.std(d) - rms > 0.05))
            row["partners"] = {}                    # model-free check: the same fit against every longer AMPL partner
            for q in D["rows"]:
                if q["amp"] != r["amp"] and q["end"][ch] > r["end"][ch] + int(0.2 * FS):
                    tq, dq = file_d(D, r, ch, ref=q)
                    if len(dq) >= 20:
                        Tq, oq, _ = fit_start(tq, dq, T, length)
                        row["partners"][str(q["amp"])] = (Tq + oq, (q["end"][ch] - r["end"][ch]) / FS)
            rows.append(row)
    return rows


def pct(x, qs=(25, 50, 75)):
    return "/".join(f"{np.percentile(x, q):.3f}" for q in qs)


FALLBACK = {"M": 1.10, "L": 0.50}   # fade start for files with no frames above the floor: M 1.1 s, L 0.5 x the file length
GUARD = 0.05                        # s before the fade start (90 % of L / R start differences are under 0.03-0.04 s)


def report_stage(workers):
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for rr in ex.map(_fit_group, curves()) for r in rr]
    json.dump(rows, open(os.path.join(OUT, "fits.json"), "w"))
    out = ["# Tomkin's fade in M and L (tools/ml_fade.py; per file and channel, frames 25 ms, 250 Hz-4 kHz, mean over octaves)"]
    for p in ("M", "L"):
        A = [r for r in rows if r["prog"] == p]
        R = [r for r in A if r["observable"]]
        if not R:
            continue
        ln = np.array([r["length"] for r in R])
        T, off, st = (np.array([r[k] for r in R]) for k in ("T", "off", "start"))
        out.append(f"\n## {p}: {len({r['key'] for r in A})} groups, {len(A)} file-channels, fade seen in {len(R)} (elsewhere the reverb "
                   f"meets the floor first); file length {pct(ln)} s (quartiles)")
        out.append("- shape, ramp ending at the file's end (gain only): T quartiles, median rms (dB):")
        for s in SHAPES:
            out.append(f"  {s:14s} T {pct([r['fits'][s][0] for r in R])} s; rms {np.median([r['fits'][s][1] for r in R]):.3f}")
        pw = np.array([r["power"] for r in R])
        out.append(f"  g = u^p free: T {pct(pw[:, 0])} s, p {pct(pw[:, 1])}, rms {np.median(pw[:, 2]):.3f}")
        out.append(f"  linear, slope free within 1 % of the decay rate: T {pct([r['ridge'][0] for r in R])} s, "
                   f"b {pct([r['ridge'][3] for r in R])} dB/s, rms {np.median([r['ridge'][1] for r in R]):.3f}")
        out.append(f"- linear, length and end free: T {pct(T)} s (10/90 % {np.percentile(T, 10):.2f}/{np.percentile(T, 90):.2f}); "
                   f"end {pct(off)} s before the file's end (< 0: cut inside the ramp); rms {np.median([r['rms'] for r in R]):.3f}")
        h = np.histogram(T, bins=[0, 0.25, 0.35, 0.5, 0.7, 0.9, 0.97, 1.03, 1.1, 1.5, 3, 30])
        out.append("  T histogram: " + " ".join(f"{a:g}-{b:g}: {n}" for a, b, n in zip(h[1][:-1], h[1][1:], h[0])))
        out.append(f"  fade start {pct(st)} s before the end; start / length {pct(st / ln)} (10/90 % {np.percentile(st / ln, 10):.2f}/"
                   f"{np.percentile(st / ln, 90):.2f}); T against file length: Spearman {spearman(T, ln):+.2f}, RT "
                   f"{spearman(T, np.array([r['rt'] for r in R])):+.2f}, AMPL {spearman(T, np.array([r['amp'] for r in R])):+.2f}")
        for lo, hi in ((0, 1.5), (1.5, 3), (3, 6), (6, 12), (12, 50)):
            s = (ln >= lo) & (ln < hi)
            if s.sum():
                out.append(f"  length {lo}-{hi} s ({s.sum()}): T {pct(T[s])}, start / length {pct(st[s] / ln[s])}")
        by = {(r["key"], r["amp"], r["ch"]): r for r in R}
        lr = [by[k]["start"] - by[(k[0], k[1], 1)]["start"] for k in by if k[2] == 0 and (k[0], k[1], 1) in by]
        out.append(f"- fade start, L - R: {pct(lr)} s; |L - R| 90 % {np.percentile(np.abs(lr), 90):.3f} s ({len(lr)} files)")
        pr = [(s_ - r["start"], r["start"]) for r in R for a, (s_, dl) in r["partners"].items()
              if (r["key"], int(a), r["ch"]) in by and dl > by[(r["key"], int(a), r["ch"])]["start"] + 0.1]
        if pr:
            pr = np.array(pr)
            out.append(f"- partner check (no model; {len(pr)} file-channel / partner pairs, the partner unfaded up to this file's end): "
                       f"start against the partner - against the model {pct(pr[:, 0])} s")
    # the table ml_lines.usable_end reads: per file the fade start (s before the end; mean over the channels that see it)
    # (seen: the mean over the channels that see it; not seen: the latest floor point over the channels, before which
    # the file shows no fade)
    tab = {"_fallback": FALLBACK, "_guard": GUARD}
    for r in rows:
        e = tab.setdefault(r["file"], {"prog": r["prog"], "seen": [], "floor": []})
        (e["seen"] if r["observable"] else e["floor"]).append(round(r["start"] if r["observable"] else r["tmin"], 4))
    for k, v in tab.items():
        if not k.startswith("_"):
            v["start"], v["seen"] = (float(np.mean(v["seen"])), True) if v["seen"] else (float(max(v["floor"])), False)
            del v["floor"]
    json.dump(tab, open(os.path.join(OUT, "fade.json"), "w"), indent=0)
    nseen = sum(1 for k, v in tab.items() if not k.startswith("_") and v["seen"])
    out.append(f"\nfade.json: {len(tab) - 2} files, fade start seen in {nseen}, floor point otherwise; fallback (no frames) {FALLBACK}, guard {GUARD} s")
    print("\n".join(out))
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(out) + "\n")


def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b)[0])


# ---------------------------------------------------------------- rt
def _rt_job(a):
    """One group: ml_lines.band_decay (per octave, onset to the usable end) with the old end (0.72 s before the file's
    end, every program) and the new one (usable_end with the per-file fade)."""
    prog, key, path = a
    import ml_carryover as cc
    g = cc.groups(prog)[key]
    rows = []
    for amp, (y, mdl, new) in mll.fitted(prog, key).items():
        old = mll.usable_end(y)
        t0 = [K0 + cc.ONSET[prog][ch] for ch in (0, 1)]
        rows.append(dict(amp=amp, file=os.path.basename(g[amp]), end_old=old, end_new=new,
                         old=mll.band_decay(y, mdl, t0, old), new=mll.band_decay(y, mdl, t0, new)))
    json.dump(dict(prog=prog, key=key, rows=rows), open(path, "w"))
    print(f"{prog} {key}: ends moved " + " ".join(f"{(o - n) / FS:+.2f}" for r in rows for o, n in zip(r['end_old'][1:], r['end_new'][1:])) + " s", flush=True)
    return key


def rt_stage(progs, workers):
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, "rt")
    os.makedirs(d, exist_ok=True)
    jobs = [(p, k, os.path.join(d, f"{p}_{k}.json")) for p in progs for k in sorted(mll.cc.groups(p))
            if not os.path.exists(os.path.join(d, f"{p}_{k}.json"))]
    print(f"{len(jobs)} groups to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_rt_job, jobs):
            pass


def rt_report():
    """ml_structure §8.3's table, measured RT / law RT per octave (median over captures and channels) by LF step, with the
    old and the new usable end."""
    out = ["\n# Per-octave RT / law (ml_structure §8.3), usable end old (0.72 s) -> new (per-file fade)"]
    for p in ("M", "L"):
        G = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "rt", f"{p}_*.json")))]
        if not G:
            continue
        mv = [(o - n) / FS for D in G for r in D["rows"] for o, n in zip(r["end_old"], r["end_new"])]
        out.append(f"\n## {p}: {len(G)} groups; end moved earlier by {pct(mv)} s (quartiles)")
        out.append("  " + " " * 16 + " ".join(f"{f:>13d}" for f in mll.OCT))
        for lf in sorted({float(D["key"].split("_")[0]) for D in G}):
            cells = []
            for i, f in enumerate(mll.OCT):
                rr = {"old": [], "new": []}
                for D in G:
                    if float(D["key"].split("_")[0]) != lf:
                        continue
                    rtl = mll.rt_law_bands(p, D["key"])[i]
                    for r in D["rows"]:
                        for v in rr:
                            for ch in (0, 1):
                                b = r[v][str(f)][ch]
                                if b is not None and b[2] >= mll.T_MIN and 60 / rtl - b[0] > 0:
                                    rr[v].append(60 / (60 / rtl - b[0]) / rtl)
                cells.append(f"{np.median(rr['old']):.3f}->{np.median(rr['new']):.3f}" if rr["new"] else "-")
            out.append(f"  LF {lf:<5}         " + " ".join(f"{c:>13s}" for c in cells))
    print("\n".join(out))
    open(os.path.join(OUT, "rt.log"), "w").write("\n".join(out) + "\n")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["measure", "report", "rt", "rtreport"], required=True)
    ap.add_argument("--prog", default="ML")
    ap.add_argument("--only", type=int, default=None, help="measure: smoke test, this many groups per program")
    ap.add_argument("--keys", default=None, help="measure: comma-separated prog_key, e.g. L_0.1_10.0_0.1")
    ap.add_argument("--workers", type=int, default=8)
    o = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if o.stage == "measure":
        keys = [(k[0], k[2:]) for k in o.keys.split(",")] if o.keys else None
        measure_stage(o.prog, o.workers, o.only, keys)
    elif o.stage == "report":
        report_stage(o.workers)
    elif o.stage == "rt":
        rt_stage(o.prog, o.workers)
    if o.stage in ("rt", "rtreport"):
        rt_report()


if __name__ == "__main__":
    main()
