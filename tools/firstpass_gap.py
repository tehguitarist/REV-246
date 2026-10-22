"""NEXT_STEPS §5 item 5: the first-pass gap shared by S's compact model and M/L (docs/ml_structure.md §7.6, §8.1;
model_s_results.md §16.11, §17). The first pass is the cluster + diffuser + 20 equal sign-only taps per channel (+ the M/L
output kernel), with per capture one shift, one cluster gain and one tap gain (tools/ml_carryover.py, firstpass stage).
Its held-out nulls sit 10-15 dB above the repeat-capture null; this tool asks what the residual is.

    .venv/bin/python tools/firstpass_gap.py --stage collect [--model lines|carry] [--only N] [--workers 8]
    .venv/bin/python tools/firstpass_gap.py --stage report [--model lines|carry]
    .venv/bin/python tools/firstpass_gap.py --stage chain                 # needs collect --model lines

Results: docs/ml_structure.md §10, §12. With EMT_KERNEL=<name> (tools/ml_structure.py) the captures come from that
kernel's cache and every output file gets the suffix _<name>. Models: "lines" (default) is ml_lines.Model with a shift and clock ratio fitted per
capture (§8.1); "carry" is ml_carryover's first pass without a clock fit (§7.6).

Stages:
  collect  every AMPL group of S, M and L: the captures on the K0 grid, the model parts, the per-capture fit (shift,
           clock for "lines", cluster gain, tap gain per channel), and per channel the first-pass window of: the capture
           minus its fitted cluster, the fitted tap part, and the same for the repeat model (another capture of the
           group, with the same freedoms); "lines" also keeps the long spans (capture, fitted cluster, fitted reverb) from
           index 0 to the first pass + SPAN. One .npz per group in analysis/firstpass_gap/groups_<model>/ (not committed).
  report   on the cache: (1) the transfer error H(f) = sum Y T* / sum |T|^2 between data and model tap part, pooled over
           b = 0 captures, per program and channel, for the even- and odd-indexed settings; (2) the residual spectrum
           against the repeat-capture residual spectrum per band; (3) the held-out gain of a setting-independent
           correction of the tap part: short FIRs (ridge) and the mean residual template, fitted on even settings,
           scored on odd. -> report_<model>.{json,log}
  chain    one kernel g (lags LAGS) on the whole fitted model, m -> m + g * m, fitted by least squares on: another
           program, the even settings (joint or per channel), or the cluster-only regions only (no reverb sample; all
           programs, or M and L only); scored per capture (two gains per channel refitted) in the cluster-and-gap and
           first-pass windows. -> chain.{json,log}; the kernel "cluster_perch" is the deconvolution correction
           (tools/q2_injection.py: Q2_CHAIN=analysis/firstpass_gap/chain.json:cluster_perch).
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import minimize, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import ml_carryover as mc  # noqa: E402
import ml_lines  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "firstpass_gap")
K0 = mc.K0
PAD = 64                        # extra tap-part samples kept either side of the window (for the FIR test)
HK = ml_lines.SUFFIX                                         # ML_FP / ML_HK (ml_lines): outputs suffixed _fp / _hk
SPAN = 600                      # chain stage: samples kept after the first pass (spans start at index 0, lag -K0)


def fit_parts(y, parts, prog, has_cluster):
    """As ml_carryover.fit_null, but returns the shift, per-channel coefficients and the shifted capture."""
    W = mc.win(prog)
    i_end = K0 + W[1][1][1]

    def solve(s):
        ys = mc.fshift(y[:i_end + 200], -s)
        res, coefs = 0.0, []
        for ch in (0, 1):
            sl = slice(K0 + W[ch][0 if has_cluster else 1][0], K0 + W[ch][1][1])
            X = np.stack([p[sl, ch] for p in parts], 1)
            c = np.linalg.lstsq(X, ys[sl, ch], rcond=None)[0]
            coefs.append(c)
            res += np.sum((ys[sl, ch] - X @ c) ** 2)
        return res, ys, coefs
    s = minimize_scalar(lambda s: solve(s)[0], bounds=(-1.0, 1.0), method="bounded", options={"xatol": 1e-3}).x
    _, ys, coefs = solve(s)
    return float(s), ys, coefs


def _collect_job(a):
    """One group, model "carry" (ml_carryover's first pass: fractional tap positions, no clock fit) or "lines"
    (ml_lines.Model: integer taps, per capture a shift and clock ratio fitted on [0, onset + shortest line), as §8.1)."""
    prog, key, P, sgn, h, b, train, model = a
    path = os.path.join(gdir(model), f"{prog}_{key}.npz")
    if os.path.exists(path):
        return path
    caps, _, amps = mc.captures(prog, key)
    W = mc.win(prog)
    out = dict(amps=np.array(amps), b=b, train=train)
    if model == "carry":
        n = K0 + W[1][1][1] + 400
        C, T = mc.render(prog, P, sgn, mc.cluster_taps(prog), n, b=b, h=h)
    else:
        n = len(caps[amps[-1]])
        M = ml_lines.Model(prog)
        C, T = M.parts(n, key)
    for amp in amps:
        y = caps[amp][:n]
        other = [x for x in amps if x != amp][0]
        if model == "carry":
            s, ys, co = fit_parts(y, (C, T), prog, amp > 0)
            sr, ysr, cor = fit_parts(y, (C, caps[other][:n]), prog, amp > 0)
            Cw = Cr = C
            Tw, Ow = T, caps[other][:n]
            out[f"clock_{amp}"] = (s, 0.0)
        else:
            sl, nf = ml_lines.fit_window(prog, amp > 0)
            _, u, gs, (Cw, Tw) = ml_lines.fit_warp(y, [C, T], sl, nf)
            ys, co = y, [np.array(g) for g in gs]
            cl_part = Cw * np.array([gs[ch][0] for ch in (0, 1)])[None]
            Cr = cl_part

            def cost(v):
                mdl, _ = ml_lines.gains(y[:nf], [cl_part[:nf], ml_lines.warp(caps[other][:nf], *v)], sl)
                return sum(np.sum((y[sl[ch], ch] - mdl[sl[ch], ch]) ** 2) for ch in (0, 1))
            v = minimize(cost, (0.0, 0.0), method="Nelder-Mead",
                         options=dict(xatol=1e-3, fatol=1e-14, initial_simplex=[(0, 0), (0.2, 0), (0, 10)])).x
            Ow = ml_lines.warp(caps[other][:nf + 200], *v)
            _, gr = ml_lines.gains(y[:nf], [cl_part[:nf], Ow[:nf]], sl)
            ysr, cor = y, [np.array(g) for g in gr]
            out[f"clock_{amp}"] = u
        if model == "lines":                                                     # long spans for the chain stage
            fs_ = slice(0, K0 + max(W[ch][1][1] for ch in (0, 1)) + SPAN)                # from index 0 (lag -K0)
            out[f"yf_{amp}"] = y[fs_]
            out[f"cf_{amp}"] = Cw[fs_] * np.array([co[ch][0] for ch in (0, 1)])[None]
            out[f"tf_{amp}"] = Tw[fs_] * np.array([co[ch][1] for ch in (0, 1)])[None]
        for ch in (0, 1):
            sl = slice(K0 + W[ch][1][0], K0 + W[ch][1][1])
            out[f"y_{amp}_{ch}"] = ys[sl, ch] - co[ch][0] * Cw[sl, ch]          # capture minus its cluster
            out[f"t_{amp}_{ch}"] = co[ch][1] * Tw[sl.start - PAD:sl.stop + PAD, ch]   # fitted tap part, padded
            out[f"yr_{amp}_{ch}"] = ysr[sl, ch] - cor[ch][0] * Cr[sl, ch]
            out[f"tr_{amp}_{ch}"] = cor[ch][1] * Ow[sl, ch]
    np.savez_compressed(path, **out)
    return path


def gdir(model):
    return os.path.join(OUT, mc.ml.tagged(f"groups_{model}") + HK)


def collect_stage(workers, model, only=None):
    os.makedirs(gdir(model), exist_ok=True)
    fp = json.load(open(ml_lines.fp_path()))
    jobs = []
    for p in "SML":
        P = np.array(fp[p]["positions"])
        sg = np.array(fp[p]["signs"])
        h = np.array(fp[p]["kernel"])
        bk = {}
        for r in fp[p]["rows"]:
            bk[r["key"]] = r["b"]
        ks = sorted(mc.groups(p))
        if only:
            ks = ks[::max(1, len(ks) // only)][:only]
        jobs += [(p, k, P, sg, h, bk[k], k in fp[p]["train"], model) for k in ks]
    with ProcessPoolExecutor(workers, max_tasks_per_child=4) as ex:
        for i, path in enumerate(ex.map(_collect_job, jobs, chunksize=1)):
            if i % 20 == 0:
                print(f"{i + 1}/{len(jobs)} {os.path.basename(path)}", flush=True)


# ---------------------------------------------------------------- report
def load(prog, model):
    rows = []
    for path in sorted(glob.glob(os.path.join(gdir(model), f"{prog}_*.npz"))):
        z = np.load(path)
        key = os.path.basename(path)[len(prog) + 1:-4]
        for amp in z["amps"]:
            rows.append(dict(key=key, amp=int(amp), b=float(z["b"]), train=bool(z["train"]),
                             y=[z[f"y_{amp}_{c}"] for c in (0, 1)], t=[z[f"t_{amp}_{c}"][PAD:-PAD] for c in (0, 1)],
                             tp=[z[f"t_{amp}_{c}"] for c in (0, 1)],
                             yr=[z[f"yr_{amp}_{c}"] for c in (0, 1)], tr=[z[f"tr_{amp}_{c}"] for c in (0, 1)]))
    return rows


def db(x):
    return 10 * np.log10(x)


def null(y, m):
    return db(np.sum((y - m) ** 2) / np.sum(y ** 2))


BANDS = [(100, 200), (200, 400), (400, 800), (800, 1600), (1600, 3200), (3200, 5000), (5000, 7000), (7000, 10000)]


def spectra(rows, ch, nfft, which=("y", "t")):
    """Pooled cross and auto spectra over rows (Hann-tapered first-pass window)."""
    ky, kt = which
    w = np.hanning(len(rows[0][ky][ch]))
    Syt = Stt = Syy = 0
    for r in rows:
        Y = np.fft.rfft(w * r[ky][ch], nfft)
        T = np.fft.rfft(w * r[kt][ch], nfft)
        Syt, Stt, Syy = Syt + Y * np.conj(T), Stt + np.abs(T) ** 2, Syy + np.abs(Y) ** 2
    return Syt, Stt, Syy


def report_prog(prog, out, model):
    rows = load(prog, model)
    if not rows:
        return
    keys = sorted({r["key"] for r in rows})
    half = {k: i % 2 for i, k in enumerate(keys)}            # even-indexed settings train, odd held out
    flat = [r for r in rows if r["b"] < 0.01]                  # b = 0: the model's first pass is setting-independent
    nfft = 1 << int(np.ceil(np.log2(len(rows[0]["y"][0]) * 2)))
    f = np.fft.rfftfreq(nfft, 1 / mc.FS)
    res = dict(n_groups=len(keys), n_captures=len(rows), n_flat=len(flat))
    lines = [f"\n## {prog}: {len(keys)} groups, {len(rows)} captures ({len(flat)} with b = 0)"]

    # baseline nulls (all rows) and repeat
    for ch in (0, 1):
        m = np.median([null(r["y"][ch], r["t"][ch]) for r in rows])
        q = np.median([null(r["yr"][ch], r["tr"][ch]) for r in rows])
        res[f"null_{ch}"], res[f"repeat_{ch}"] = m, q
    lines.append(f"- first pass (minus fitted cluster), median: model {res['null_0']:.1f} / {res['null_1']:.1f}, "
                 f"repeat {res['repeat_0']:.1f} / {res['repeat_1']:.1f} dB (L / R)")

    # (1) transfer error, b = 0 rows, per half
    for ch in (0, 1):
        H = {}
        for hv in (0, 1):
            sel = [r for r in flat if half[r["key"]] == hv] or flat
            Syt, Stt, _ = spectra(sel, ch, nfft)
            H[hv] = Syt / Stt
        lines.append(f"- {'LR'[ch]}: transfer error H = Y/T (b = 0 rows), per band: |H| dB, angle deg, "
                     f"[train half / held-out half]")
        bandrows = []
        for lo, hi in BANDS:
            i = (f >= lo) & (f < hi)
            a = [20 * np.log10(np.abs(np.mean(H[hv][i]))) for hv in (0, 1)]
            ph = [np.degrees(np.angle(np.mean(H[hv][i]))) for hv in (0, 1)]
            e = [db(np.mean(np.abs(H[hv][i] - 1) ** 2)) for hv in (0, 1)]
            bandrows.append(dict(lo=lo, hi=hi, mag_db=a, phase_deg=ph, err_db=e))
            lines.append(f"  {lo:>5}-{hi:<5} Hz: {a[0]:+.3f} / {a[1]:+.3f} dB, {ph[0]:+.2f} / {ph[1]:+.2f} deg, "
                         f"|H-1|^2 {e[0]:.1f} / {e[1]:.1f} dB")
        res[f"H_bands_{ch}"] = bandrows
        res[f"H_curve_{ch}"] = dict(f=f[::4].tolist(), re=[np.real(H[0][::4]).tolist(), np.real(H[1][::4]).tolist()],
                                    im=[np.imag(H[0][::4]).tolist(), np.imag(H[1][::4]).tolist()])

    # (2) residual spectrum vs repeat residual spectrum, all rows, per band (re the data)
    for ch in (0, 1):
        w = np.hanning(len(rows[0]["y"][ch]))
        Sr = Sq = Sy = 0
        for r in rows:
            Sr = Sr + np.abs(np.fft.rfft(w * (r["y"][ch] - r["t"][ch]), nfft)) ** 2
            Sq = Sq + np.abs(np.fft.rfft(w * (r["yr"][ch] - r["tr"][ch]), nfft)) ** 2
            Sy = Sy + np.abs(np.fft.rfft(w * r["y"][ch], nfft)) ** 2
        lines.append(f"- {'LR'[ch]}: residual re data per band, model / repeat (excess)")
        bb = []
        for lo, hi in BANDS:
            i = (f >= lo) & (f < hi)
            m, q = db(Sr[i].sum() / Sy[i].sum()), db(Sq[i].sum() / Sy[i].sum())
            bb.append(dict(lo=lo, hi=hi, model=m, repeat=q, share=float(Sr[i].sum() / Sr[(f >= 100) & (f < 10000)].sum())))
            lines.append(f"  {lo:>5}-{hi:<5} Hz: {m:.1f} / {q:.1f} ({m - q:+.1f}); share of the model residual "
                         f"{100 * bb[-1]['share']:.0f} %")
        res[f"resid_bands_{ch}"] = bb

    # (3) setting-independent corrections of the tap part, b = 0 rows: train on even settings, score odd
    tr = [r for r in flat if half[r["key"]] == 0]
    te = [r for r in flat if half[r["key"]] == 1]
    for ch in (0, 1):
        base = np.median([null(r["y"][ch], r["t"][ch]) for r in te])
        # mean residual template, normalised by the tap gain (the tap part's scale)
        g = lambda r: np.sqrt(np.sum(r["t"][ch] ** 2))
        tmpl = np.mean([(r["y"][ch] - r["t"][ch]) / g(r) for r in tr], 0)
        tm = np.median([null(r["y"][ch], r["t"][ch] + g(r) * tmpl) for r in te])
        # with the template, the tap gain refit (template scaled together with the taps)
        tm2 = []
        for r in te:
            m0 = r["t"][ch] / g(r) + tmpl
            c = np.sum(r["y"][ch] * m0) / np.sum(m0 * m0)
            tm2.append(null(r["y"][ch], c * m0))
        tm2 = np.median(tm2)
        line = f"- {'LR'[ch]}: held out (b = 0, odd settings, n {len(te)}): model {base:.1f}"
        line += f"; + mean residual template ({len(tmpl)} numbers) {tm:.1f}, with gain refit {tm2:.1f}"
        res[f"corr_{ch}"] = dict(base=base, template=tm, template_refit=tm2)
        for K in (4, 16, 48):
            lags = range(-K, K + 1)
            def design(r):
                n = len(r["t"][ch])
                return np.stack([r["tp"][ch][PAD - k:PAD - k + n] for k in lags], 1)
            A = sum(design(r).T @ design(r) for r in tr)
            v = sum(design(r).T @ r["y"][ch] for r in tr)
            lam = 1e-6 * np.trace(A) / len(lags)
            hf = np.linalg.solve(A + lam * np.eye(len(lags)), v)
            nn = []
            for r in te:
                m0 = design(r) @ hf
                c = np.sum(r["y"][ch] * m0) / np.sum(m0 * m0)
                nn.append(null(r["y"][ch], c * m0))
            res[f"corr_{ch}"][f"fir{K}"] = float(np.median(nn))
            res[f"corr_{ch}"][f"fir{K}_h"] = hf.tolist()
            line += f"; FIR ±{K} {np.median(nn):.1f}"
        lines.append(line)
    out[prog] = res
    print("\n".join(lines), flush=True)


def report_stage(model):
    out = {}
    for p in "SML":
        report_prog(p, out, model)
    json.dump(out, open(os.path.join(OUT, mc.ml.tagged(f"report_{model}") + HK + ".json"), "w"), indent=1)



# ---------------------------------------------------------------- chain: one low-frequency kernel on the whole model
LAGS = np.arange(-32, 401)          # kernel lags (device samples); the transfer error's kernel sits at 50-150


def fit_slices(prog, has_cluster):
    """Per channel, the fit span (cluster and gap + first pass; first pass only without a cluster) as indices into the
    stored long spans (which start at lag -K0), and the two scoring windows."""
    W = mc.win(prog)
    out = []
    for ch in (0, 1):
        a = K0 + (W[ch][0][0] - 30 if has_cluster else W[ch][1][0])
        out.append((slice(a, K0 + W[ch][1][1]), {"cluster_gap": slice(K0 + W[ch][0][0], K0 + W[ch][0][1]),
                                                 "first_pass": slice(K0 + W[ch][1][0], K0 + W[ch][1][1])}))
    return out


def lagged(x, sl):
    """Columns x(t - lag) for t in sl, lag in LAGS (zero outside x)."""
    n = len(x)
    idx = np.arange(sl.start, sl.stop)[:, None] - LAGS[None]
    return np.where((idx >= 0) & (idx < n), x[np.clip(idx, 0, n - 1)], 0.0)


def _gram_job(path):
    """Normal equations of g for y - m = (g * m), per channel, summed over the group's captures: on the fit span of
    the b = 0 captures ("fit"), and on the cluster-only span [lag -30, onset - 3) of every capture with a cluster
    ("cluster": no reverb sample enters it, any b)."""
    z = np.load(path)
    prog = os.path.basename(path)[0]
    out = {"fit": [[0, 0], [0, 0]], "cluster": [[0, 0], [0, 0]]}
    W = mc.win(prog)
    for amp in z["amps"]:
        y, m = z[f"yf_{amp}"], z[f"cf_{amp}"] + z[f"tf_{amp}"]
        for ch, (sl, _) in enumerate(fit_slices(prog, amp > 0)):
            spans = []
            if float(z["b"]) < 0.01:
                spans.append(("fit", sl))
            if amp > 0:
                spans.append(("cluster", slice(K0 + W[ch][0][0] - 30, K0 + W[ch][0][1])))
            for name, sp in spans:
                X = lagged(m[:, ch], sp)
                out[name][0][ch] = out[name][0][ch] + X.T @ X
                out[name][1][ch] = out[name][1][ch] + X.T @ (y[sp, ch] - m[sp, ch])
    return path, out


def conv_lags(x, g):
    """(g * x)(t) = sum_k g_k x(t - LAGS_k), same length as x."""
    out = np.zeros_like(x)
    for k, gk in zip(LAGS, g):
        if k >= 0:
            out[k:] += gk * x[:len(x) - k]
        else:
            out[:k] += gk * x[-k:]
    return out


def _score_job(a):
    """Every capture of one group: nulls in the two windows with each named kernel (gains refitted on the fit span)."""
    path, kernels = a
    z = np.load(path)
    prog = os.path.basename(path)[0]
    key = os.path.basename(path)[2:-4]
    rows = []
    for amp in z["amps"]:
        y = z[f"yf_{amp}"]
        row = dict(prog=prog, key=key, amp=int(amp), b=float(z["b"]))
        for name, gk in kernels.items():
            res = {}
            for ch, (sl, wins) in enumerate(fit_slices(prog, amp > 0)):
                cols = []
                for part in ("cf", "tf"):
                    x = z[f"{part}_{amp}"][:, ch]
                    cols.append(x + conv_lags(x, gk[ch]) if gk is not None else x)
                X = np.stack([c[sl] for c in cols], 1)
                c = np.linalg.lstsq(X, y[sl, ch], rcond=None)[0]
                mdl = c[0] * cols[0] + c[1] * cols[1]
                for w, ws in wins.items():
                    if w == "cluster_gap" and amp == 0:
                        continue
                    res.setdefault(w, [None, None])[ch] = null(y[ws, ch], mdl[ws])
            row[name] = res
        row["repeat_first_pass"] = [null(z[f"yr_{amp}_{ch}"], z[f"tr_{amp}_{ch}"]) for ch in (0, 1)]
        rows.append(row)
    return rows


def chain_stage(workers):
    paths = sorted(glob.glob(os.path.join(gdir("lines"), "*.npz")))
    with ProcessPoolExecutor(workers) as ex:
        grams = dict(ex.map(_gram_job, paths, chunksize=4))
    prog_of = lambda p: os.path.basename(p)[0]
    keys = {q: sorted(os.path.basename(p)[2:-4] for p in paths if prog_of(p) == q) for q in "SML"}
    even = lambda p: keys[prog_of(p)].index(os.path.basename(p)[2:-4]) % 2 == 0

    def solve(sel, joint, span="fit"):
        G = [sum(grams[p][span][0][ch] for p in paths if sel(p) and not np.isscalar(grams[p][span][0][ch])) for ch in (0, 1)]
        v = [sum(grams[p][span][1][ch] for p in paths if sel(p) and not np.isscalar(grams[p][span][1][ch])) for ch in (0, 1)]
        if joint:
            G = [G[0] + G[1]] * 2
            v = [v[0] + v[1]] * 2
        return [np.linalg.solve(Gc + 1e-9 * np.trace(Gc) / len(LAGS) * np.eye(len(LAGS)), vc) for Gc, vc in zip(G, v)]
    K = {"none": None}
    for q in "SML":
        K[f"from_{q}"] = solve(lambda p, q=q: prog_of(p) == q, True)
    K["from_even"] = solve(even, True)
    K["from_even_perch"] = solve(even, False)
    K["cluster_perch"] = solve(lambda p: True, False, "cluster")              # never sees a reverb sample
    K["cluster_ML_perch"] = solve(lambda p: prog_of(p) != "S", False, "cluster")
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for rr in ex.map(_score_job, [(p, K) for p in paths], chunksize=2) for r in rr]
    json.dump(dict(lags=LAGS.tolist(), kernels={k: (None if g is None else [x.tolist() for x in g]) for k, g in K.items()},
                   rows=rows), open(os.path.join(OUT, mc.ml.tagged("chain") + HK + ".json"), "w"))
    names = list(K)
    for q in "SML":
        print(f"\n## {q}: one kernel on the whole model (cluster + reverb), lags {LAGS[0]}..{LAGS[-1]}; median null L / R, dB; "
              f"held out = settings or programs the kernel was not fitted on")
        for sel_name, sel in (("odd settings", lambda r: keys[q].index(r["key"]) % 2 == 1), ("all", lambda r: True)):
            rr = [r for r in rows if r["prog"] == q and sel(r)]
            rep = np.median([r["repeat_first_pass"] for r in rr], 0)
            print(f"- {sel_name} ({len(rr)} captures; first-pass repeat {rep[0]:.1f} / {rep[1]:.1f}):")
            for nm in names:
                cells = []
                for w in ("cluster_gap", "first_pass"):
                    x = np.array([r[nm][w] for r in rr if w in r[nm]])
                    cells.append(f"{w} {np.median(x[:, 0]):.1f} / {np.median(x[:, 1]):.1f}")
                held = ("in-sample" if nm == f"from_{q}" or (nm.startswith("from_even") and sel_name == "all") else
                        "held out first pass" if nm.startswith("cluster") else "held out")
                print(f"  {nm:16s} ({held}): " + "; ".join(cells))
    for k, g in K.items():
        if g is not None:
            e = g[0]
            top = LAGS[np.argsort(np.abs(e))[::-1][:6]]
            print(f"- kernel {k}: energy {db(np.sum(e ** 2)):.1f} dB, largest at lags {sorted(top.tolist())}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="report", choices=["collect", "report", "chain"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", type=int)
    ap.add_argument("--model", default="lines", choices=["carry", "lines"])
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.stage == "collect":
        collect_stage(a.workers, a.model, a.only)
    elif a.stage == "chain":
        chain_stage(a.workers)
    else:
        report_stage(a.model)


if __name__ == "__main__":
    main()
