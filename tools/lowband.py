"""NEXT_STEPS §1 next candidate: program S's low band (25-100 Hz) in the recirculations. ml_structure §12.5 found S's
recirculations 15-20 dB above the repeat-capture null at 25-100 Hz (ml_lines' law model), a band no null so far scores.
This scores the current S baseline (model_s_results §22.9: the final model with the cross-validated loop terms, each
setting with the constants of the fold that held out its RT step) there. Results: docs/model_s_results.md §23.

    EMT_KERNEL=chain2 .venv/bin/python tools/lowband.py --stage score [--only keys] [--workers 4] [--variant name]
    EMT_KERNEL=chain2 .venv/bin/python tools/lowband.py --stage report [--variant name]

Per capture as q2_injection.final (clock, shift and gain fitted before the fade, 20 Hz-11 kHz data, fade-compensated
model), then per band (zero-phase 3rd-order Butterworth band-passes, as outofband.py; 100 Hz-10 kHz is the scored band
of every earlier null):
  nulls   per window of q2's W_COMP up to the compensated end: strict (the in-band gain) and gain-free (best scalar per
          window and channel, so a level error and a shape error separate);
  repeat  the same against the next AMPL capture of the setting (as q2's nulls stage: y_a g_b against y_b g_a);
  decay   the slope (dB/s) of 10 log10(E_capture / E_model) over 39 ms frames from 50 ms to the compensated end
          (energy-weighted, as outofband.decay), and the same capture against repeat capture.

--variant applies a candidate loop change on top of the baseline (VARIANTS); outputs go to analysis/lowband/<kernel>_<variant>/.
One result file per setting, skipped on rerun.
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import loop_fir as lf_  # noqa: E402

q2 = lf_.q2
q1 = q2.q1
FS, PRE = q2.FS, q2.PRE
CV = os.path.join(ROOT, "analysis", "loop_cv")
BANDS = {"25-50": (25, 50), "50-100": (50, 100), "100-200": (100, 200), "100-10k": (100, 10000)}
SOS = {b: butter(3, lohi, btype="band", fs=FS, output="sos") for b, lohi in BANDS.items()}
FRAME = 1024                                    # 39 ms, as ml_lines.FRAME
VARIANTS = {"base": None, "tfix": None}       # name -> Gmod (a loop change, q2.final_G) or None
POST = {"tfix"}                                 # variants that filter the whole model by a fixed per-channel T(f) (tfit stage)
TFOLDS = {"A": (0.3, 0.8, 1.4, 2.2), "B": (0.6, 1.0, 1.8, 2.5)}   # tfit: trained on these RT steps, applied to the others
TF_SMOOTH, TF_TAPER = 1.5, (300.0, 400.0)      # tfit: running mean (+/- Hz) and the taper to 1


def fold_terms(rt):
    """The loop terms of the fold that held out RT step rt (loop_cv.py), as final_G arguments."""
    f = json.load(open(os.path.join(CV, f"fold_{rt}.json")))
    return dict(gamma=f["gamma"], lfrt=f["lfrt"], hfdp=f["hfdp"], pramp=tuple(f["pramp"]))


def outdir(variant):
    return os.path.join(ROOT, "analysis", "lowband", (q1.ms.KTAG or "library") + ("" if variant == "base" else "_" + variant))


def band_nulls(y, m, te):
    """Per band and window: [strict, gain-free] per channel (None where the window is under 50 ms)."""
    out = {}
    for b, sos in SOS.items():
        yf, mf = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, m, axis=0)
        d = {}
        for a_ms, b_ms in q2.W_COMP:
            i0, i1 = PRE + int(a_ms / 1000 * FS), min(PRE + int(b_ms / 1000 * FS), PRE + int(te * FS), len(y))
            if i1 - i0 < int(0.05 * FS):
                continue
            cells = []
            for c in (0, 1):
                yy, mm = yf[i0:i1, c], mf[i0:i1, c]
                e2, y2, m2, ym = np.sum((yy - mm) ** 2), np.sum(yy * yy), np.sum(mm * mm), np.sum(yy * mm)
                cells.append([float(10 * np.log10(e2 / y2)), float(10 * np.log10(max(1 - ym * ym / (y2 * m2), 1e-12)))])
            d[f"{a_ms}-{b_ms}"] = cells
        out[b] = d
    return out


def band_decay(y, m, te, bands=None):
    """Per band (default SOS), per channel: [slope dB/s, level dB] of 10 log10(E_y / E_m) from 50 ms to te, or None."""
    out = {}
    i0, i1 = PRE + int(0.05 * FS), min(PRE + int(te * FS), len(y))
    k = (i1 - i0) // FRAME
    for b, sos in (bands or SOS).items():
        yf, mf = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, m, axis=0)
        row = []
        for c in (0, 1):
            if k < 4:
                row.append(None)
                continue
            ey = (yf[i0:i0 + k * FRAME, c] ** 2).reshape(k, FRAME).mean(1)
            em = (mf[i0:i0 + k * FRAME, c] ** 2).reshape(k, FRAME).mean(1)
            ok = (ey > ey.max() * 1e-7) & (em > em.max() * 1e-7)
            if ok.sum() < 4:
                row.append(None)
                continue
            t = (np.arange(k) + 0.5) * FRAME / FS
            wgt = np.sqrt(em[ok] / em[ok].max())
            A = np.stack([np.ones(ok.sum()), t[ok]], 1) * wgt[:, None]
            cf = np.linalg.lstsq(A, 10 * np.log10(ey[ok] / em[ok]) * wgt, rcond=None)[0]
            row.append([float(cf[1]), float(cf[0])])
        out[b] = row
    return out


def tfix_path(fold):
    return os.path.join(ROOT, "analysis", "lowband", (q1.ms.KTAG or "library") + f"_tfix_{fold}.json")


def post_filter(mdl, rt):
    """The whole model (cluster and reverb) times the fixed T(f) of the fold that did not train on RT step rt."""
    fold = next(k for k, rts in TFOLDS.items() if rt not in rts)
    d = json.load(open(tfix_path(fold)))
    N = 1 << int(np.ceil(np.log2(len(mdl) + 8192)))
    f = np.fft.rfftfreq(N, 1 / FS)
    out = np.empty_like(mdl)
    for c in (0, 1):
        T = np.interp(f, d["f"], d["T_re"][c], right=1.0) + 1j * np.interp(f, d["f"], d["T_im"][c], right=0.0)
        out[:, c] = np.fft.irfft(np.fft.rfft(mdl[:, c], N) * T, N)[:len(mdl)]
    return out


def tfit_stage():
    """Per fold: the energy-weighted capture / model T(f) of the transfer stage's LF 1.0 / 2.0, HF 1.0 settings at the
    fold's RT steps, smoothed (+/- TF_SMOOTH Hz), held at its 15 Hz value below 15 Hz, tapered to 1 over TF_TAPER."""
    D = [d for d in (json.load(open(p)) for p in glob.glob(os.path.join(outdir("base") + "_transfer", "*.json"))) if d]
    for fold, rts in TFOLDS.items():
        ds = [d for d in D if d["lf"] in (1.0, 2.0) and d["hf"] == 1.0 and d["rt"] in rts]
        f = np.array(ds[0]["f"])
        E = sum(np.array(d["e_model"]) for d in ds)
        T = sum((np.array(d["T_re"]) + 1j * np.array(d["T_im"])) * np.array(d["e_model"]) for d in ds) / np.maximum(E, 1e-30)
        k = int(round(TF_SMOOTH / (f[1] - f[0])))
        ker = np.ones(2 * k + 1) / (2 * k + 1)
        Ts = np.stack([np.convolve(T[:, c], ker, mode="same") for c in (0, 1)], 1)
        Ts[f < 15] = Ts[np.searchsorted(f, 15)]
        a = np.clip((f - TF_TAPER[0]) / (TF_TAPER[1] - TF_TAPER[0]), 0, 1)[:, None]
        Ts = (1 - a) * Ts + a
        json.dump(dict(fold=fold, rt=list(rts), settings=sorted(d["key"] for d in ds), f=f.tolist(),
                       T_re=Ts.real.T.tolist(), T_im=Ts.imag.T.tolist()), open(tfix_path(fold), "w"))
        print(f"fold {fold}: {len(ds)} settings; " + "; ".join(
            f"{fb} Hz " + " / ".join(f"{20 * np.log10(abs(Ts[i, c])):+.2f} dB {np.angle(Ts[i, c], deg=True):+.1f}" for c in (0, 1))
            for fb in (25, 40, 60, 100, 150) for i in [int(np.searchsorted(f, fb))]))


def _setting(a_):
    key, lab, b230, law, v, caps, ends, variant, path = a_
    q1.EXT = 1700
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    G = q2.final_G(lab, b230, law, **fold_terms(lab["rt"]), Gmod=VARIANTS[variant])
    keep, rows = {}, []
    for c in caps:
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        if variant in POST:
            mdl = post_filter(mdl, lab["rt"])
        keep[c["amp"]] = (y, g, te)
        rows.append(dict(key=key, **lab, amp=c["amp"], file=c["file"], t_comp=te, nulls=band_nulls(y, mdl * g, te),
                         decay=band_decay(y, mdl * g, te)))
    amps = sorted(keep)
    for r in rows:
        b = amps[(amps.index(r["amp"]) + 1) % len(amps)]
        (ya, ga, ta), (yb, gb, tb) = keep[r["amp"]], keep[b]
        m = min(len(ya), len(yb))
        r["partner"] = b
        r["repeat"] = band_nulls(ya[:m] * gb[:m], yb[:m] * ga[:m], min(ta, tb))
        r["repeat_decay"] = band_decay(ya[:m] * gb[:m], yb[:m] * ga[:m], min(ta, tb))
    json.dump(rows, open(path, "w"))
    return key


def _transfer(a_):
    """Per setting: the capture / model transfer T(f) = sum Y M* / sum |M|^2 over the reverb (T_FROM to each capture's
    compensated end; Hann-tapered ends), all captures pooled, per channel, on a TF_DF grid to TF_FMAX."""
    key, lab, b230, law, v, caps, ends, variant, path = a_
    q1.EXT = 1700
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    G = q2.final_G(lab, b230, law, **fold_terms(lab["rt"]), Gmod=VARIANTS[variant])
    N = int(round(FS / TF_DF))
    f = np.fft.rfftfreq(N, 1 / FS)
    sel = f <= TF_FMAX
    num, den, yy = np.zeros((sel.sum(), 2), complex), np.zeros((sel.sum(), 2)), np.zeros((sel.sum(), 2))
    for c in caps:
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        i0, i1 = PRE + int(T_FROM * FS), min(PRE + int(te * FS), len(y), PRE + N)
        if i1 - i0 < int(0.3 * FS):
            continue
        w = np.ones(i1 - i0)
        nt = int(0.02 * FS)
        w[:nt] = w[-nt:] = 0
        w[:nt], w[-nt:] = np.hanning(2 * nt)[:nt], np.hanning(2 * nt)[nt:]
        Y = np.fft.rfft(y[i0:i1] * w[:, None], N, axis=0)[sel]
        M = np.fft.rfft((mdl * g)[i0:i1] * w[:, None], N, axis=0)[sel]
        num += Y * np.conj(M)
        den += np.abs(M) ** 2
        yy += np.abs(Y) ** 2
    if not den.any():
        json.dump(None, open(path, "w"))
        return key
    T = num / np.maximum(den, 1e-30)
    coh = np.abs(num) ** 2 / np.maximum(den * yy, 1e-30)
    json.dump(dict(key=key, **lab, f=f[sel].tolist(), T_re=T.real.tolist(), T_im=T.imag.tolist(), coh=coh.tolist(),
                   e_model=den.tolist()), open(path, "w"))
    return key


T_FROM, TF_DF, TF_FMAX = 0.15, 0.5, 400.0
THIRDS = 20.0 * 2 ** (np.arange(16) / 3)        # slopes: third-octave centres 20-640 Hz


def _slopes(a_):
    """Per setting and capture: the decay-rate difference capture - model (dB/s; band_decay's slope) in third-octave
    bands, from 50 ms to the compensated end, and the same for the repeat capture."""
    key, lab, b230, law, v, caps, ends, variant, path = a_
    q1.EXT = 1700
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    G = q2.final_G(lab, b230, law, **fold_terms(lab["rt"]), Gmod=VARIANTS[variant])
    bands = {f"{fc:.0f}": butter(3, (fc * 2 ** (-1 / 6), fc * 2 ** (1 / 6)), btype="band", fs=FS, output="sos") for fc in THIRDS}
    keep, rows = {}, []
    for c in caps:
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        keep[c["amp"]] = (y, mdl * g, te)
    amps = sorted(keep)
    for a in amps:
        y, m, te = keep[a]
        yb, gb_, tb = keep[amps[(amps.index(a) + 1) % len(amps)]]
        n = min(len(y), len(yb))
        rows.append(dict(key=key, **lab, amp=a, f=THIRDS.tolist(), model=list(band_decay(y, m, te, bands).values()),
                         repeat=list(band_decay(y[:n], yb[:n], min(te, tb), bands).values())))
    json.dump(rows, open(path, "w"))
    return key


def _ml_transfer(a_):
    """Cross-program check: the same T(f) for one S / M / L AMPL group with ml_lines' law model (outofband's full-band
    render: in-band clock and gains of the law stage, the pipeline response P(f)), over the deterministic part, from
    25 samples before tap 0 to the end of the first pass (onset + shortest line), Hann-tapered at the end."""
    import outofband as ob
    mlines, cc = ob.mlines, ob.cc
    prog, key, path = a_
    M = mlines.Model(prog)
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, amps = ob.captures4(prog, key, n)
    n = len(caps[amps[-1]][0])
    Craw, Rraw = M.parts(n, key, band=False)
    C, R = (sosfiltfilt(ob.ml.BAND, x, axis=0) for x in (Craw, Rraw))
    Nr = 1 << int(np.ceil(np.log2(n + 4096)))
    P = ob.pipeline_response(Nr)
    N = int(round(FS / TF_DF))
    f = np.fft.rfftfreq(N, 1 / FS)
    sel = f <= TF_FMAX
    num, den, yy = np.zeros((sel.sum(), 2), complex), np.zeros((sel.sum(), 2)), np.zeros((sel.sum(), 2))
    i0 = ob.K0 - 25
    for amp in amps:
        if amp < 30:
            continue
        yb, yr = caps[amp]
        sl, nf = mlines.fit_window(prog, True)
        _, u, gs, _ = mlines.fit_warp(yb, [C, R], sl, nf)
        m = ob.apply_f(mlines.warp(Craw, *u) * np.array([gs[ch][0] for ch in (0, 1)])[None], P) + \
            ob.apply_f(mlines.warp(Rraw, *u) * np.array([gs[ch][1] for ch in (0, 1)])[None], P)
        for ch in (0, 1):
            i1 = ob.K0 + cc.ONSET[prog][ch] + min(cc.LINES[prog])
            w = np.ones(i1 - i0)
            nt = min(int(0.02 * FS), (i1 - i0) // 4)
            w[-nt:] = np.hanning(2 * nt)[nt:]
            Y = np.fft.rfft(yr[i0:i1, ch] * w, N)[sel]
            Mm = np.fft.rfft(m[i0:i1, ch] * w, N)[sel]
            num[:, ch] += Y * np.conj(Mm)
            den[:, ch] += np.abs(Mm) ** 2
            yy[:, ch] += np.abs(Y) ** 2
    T = num / np.maximum(den, 1e-30)
    json.dump(dict(prog=prog, key=key, f=f[sel].tolist(), T_re=T.real.tolist(), T_im=T.imag.tolist(),
                   coh=(np.abs(num) ** 2 / np.maximum(den * yy, 1e-30)).tolist(), e_model=den.tolist()), open(path, "w"))
    return prog, key


def _ml_check(a_):
    """Cross-program test of the fixed T (fitted on S's reverb only, fold A): per S / M / L group, nulls of ml_lines'
    law model with and without T in the cluster, first-pass and first-recirculation windows, per band, AMPL >= 30."""
    import outofband as ob
    mlines, cc = ob.mlines, ob.cc
    prog, key, path = a_
    M = mlines.Model(prog)
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, amps = ob.captures4(prog, key, n)
    n = len(caps[amps[-1]][0])
    Craw, Rraw = M.parts(n, key, band=False)
    C, R = (sosfiltfilt(ob.ml.BAND, x, axis=0) for x in (Craw, Rraw))
    P = ob.pipeline_response(1 << int(np.ceil(np.log2(n + 4096))))
    W = mlines.windows(prog)
    rows = []
    for amp in amps:
        if amp < 30:
            continue
        yb, yr = caps[amp]
        sl, nf = mlines.fit_window(prog, True)
        _, u, gs, _ = mlines.fit_warp(yb, [C, R], sl, nf)
        m = ob.apply_f(mlines.warp(Craw, *u) * np.array([gs[ch][0] for ch in (0, 1)])[None], P) + \
            ob.apply_f(mlines.warp(Rraw, *u) * np.array([gs[ch][1] for ch in (0, 1)])[None], P)
        end = mlines.usable_end(yb)
        row = dict(prog=prog, key=key, amp=amp, nulls={})
        for vn, mv in (("model", m), ("tfix", post_filter(m, 0.6))):          # rt 0.6 -> fold A's T
            row["nulls"][vn] = {}
            for b, sos in SOS.items():
                yf, mf = sosfiltfilt(sos, yr, axis=0), sosfiltfilt(sos, mv, axis=0)
                row["nulls"][vn][b] = {w: [mlines.null_db(yf[:, ch], mf[:, ch], ob.K0 + W[ch][w][0], min(ob.K0 + W[ch][w][1], end[ch]))
                                           for ch in (0, 1)] for w in ("cluster_gap", "first_pass", "recirc1")}
        rows.append(row)
    json.dump(rows, open(path, "w"))
    return prog, key


def ml_transfer_stage(workers, every=4, stage="mltransfer"):
    import outofband as ob
    od = os.path.join(ROOT, "analysis", "lowband", (q1.ms.KTAG or "library") + ob.HK + "_" + stage)
    os.makedirs(od, exist_ok=True)
    jobs = [(p, k, os.path.join(od, f"{p}_{k}.json")) for p in "SML" for k in sorted(ob.cc.groups(p))[::every]]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    print(f"{stage}: {len(jobs)} groups to do", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        fut = [ex.submit(_ml_check if stage == "mlcheck" else _ml_transfer, j) for j in jobs]
        for i, fu in enumerate(as_completed(fut)):
            print(f"  {i + 1}/{len(jobs)} {fu.result()}", flush=True)


def score_stage(workers, only=None, variant="base", stage="score"):
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=lf_.COMPACT)
    od = outdir(variant) + ("" if stage == "score" else "_" + stage)
    os.makedirs(od, exist_ok=True)
    keys = sorted(k for k in info if only is None or k in only)
    jobs = [(k, dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
             v, [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in q2.files_of(info[k])], ends, variant,
             os.path.join(od, k + ".json")) for k in keys]
    jobs = [j for j in jobs if not os.path.exists(j[-1])]
    print(f"{stage} ({variant}): {len(jobs)} settings to do", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        fut = [ex.submit({"transfer": _transfer, "slopes": _slopes}.get(stage, _setting), j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            print(f"  {i + 1}/{len(jobs)} {f.result()}", flush=True)


def load(variant="base"):
    return [r for p in sorted(glob.glob(os.path.join(outdir(variant), "*.json"))) for r in json.load(open(p))]


def report_stage(variant="base"):
    rows = load(variant)
    ws = ["0-50", "50-150", "150-300", "300-600", "600-1000"]
    out = [f"## S baseline (§22.9, CV constants) by band ({variant}; {len({r['key'] for r in rows})} settings, {len(rows)} captures); "
           "strict null L / R [repeat L / R]; (gain-free R)"]

    def med(rr, fld, b, w, c, j):
        x = [r[fld][b][w][c][j] for r in rr if w in r[fld][b]]
        return np.median(x) if x else np.nan

    def ex(rr, b, w, c=1):
        x = [r["nulls"][b][w][c][0] - r["repeat"][b][w][c][0] for r in rr if w in r["nulls"][b] and w in r["repeat"][b]]
        return np.median(x) if x else np.nan
    for b in BANDS:
        out.append(f"\n### {b} Hz\n| window | model L / R | repeat L / R | gain-free R | excess L / R |\n|---|---|---|---|---|")
        for w in ws:
            out.append(f"| {w} | {med(rows, 'nulls', b, w, 0, 0):.1f} / {med(rows, 'nulls', b, w, 1, 0):.1f} | "
                       f"{med(rows, 'repeat', b, w, 0, 0):.1f} / {med(rows, 'repeat', b, w, 1, 0):.1f} | {med(rows, 'nulls', b, w, 1, 1):.1f} | "
                       f"{ex(rows, b, w, 0):+.1f} / {ex(rows, b, w, 1):+.1f} |")
    out.append("\n### Excess R over repeat by control (25-50 / 50-100 / 100-200 Hz), windows 50-150 | 150-300 | 300-600")
    for c in ("lf", "rt", "hf"):
        for val in sorted({r[c] for r in rows}):
            rr = [r for r in rows if r[c] == val]
            out.append(f"- {c} {val}: " + " | ".join("/".join(f"{ex(rr, b, w):+.1f}" for b in ("25-50", "50-100", "100-200")) for w in ws[1:4]))
    out.append("\n### Decay: RT capture / model from the slope (RT_law at the band from the label x LF; slope dB/s, median R) [capture / repeat]")
    for c in ("lf",):
        for val in sorted({r[c] for r in rows}):
            rr = [r for r in rows if r[c] == val]
            cells = []
            for b in ("25-50", "50-100", "100-200", "100-10k"):
                s = [r["decay"][b][1][0] for r in rr if r["decay"][b][1]]
                sr = [r["repeat_decay"][b][1][0] for r in rr if r["repeat_decay"][b][1]]
                cells.append(f"{b}: {np.median(s):+.2f} dB/s [{np.median(sr):+.2f}; IQR {np.percentile(s, 25):+.2f}..{np.percentile(s, 75):+.2f}]")
            out.append(f"- lf {val}: " + "; ".join(cells))
    open(os.path.join(os.path.dirname(outdir(variant)), (q1.ms.KTAG or "library") + f"_{variant}_report.log"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["score", "transfer", "slopes", "mltransfer", "mlcheck", "tfit", "report"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--variant", default="base", choices=sorted(VARIANTS))
    a = ap.parse_args()
    if a.stage == "tfit":
        tfit_stage()
    elif a.stage in ("mltransfer", "mlcheck"):
        ml_transfer_stage(a.workers, stage=a.stage)
    elif a.stage in ("score", "transfer", "slopes"):
        score_stage(a.workers, a.only.split(",") if a.only else None, a.variant, a.stage)
    else:
        report_stage(a.variant)


if __name__ == "__main__":
    main()
