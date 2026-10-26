"""Program E on the reverb metrics of the pulse check (tools/pulse_check.py, G1-G6), as programs S / M / L were audited.

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_pulse.py [--rts 1.0 1.5] [--seeds 5] [--out file.json]

Data: the nine Altiverb REV-250 presets on the device grid (e_model.dev, R's direct tap at K0 = 64, normalised to it). Model:
e_sim's Program E with the adopted moving reads, in two renders per seed, (a) an impulse render (S.render_ir, no sweep) and
(b) a simulated 10 s sweep capture deconvolved like the library's (S.capture; EMT250_program 13.6: the sweep length changes the
measured high-band T60 of a time-variant system). Every metric function is pulse_check's (compare, rv_metrics, band_decay,
ned_curve, iacc, spectra, coherence, fine, edr_times, nb_env_dev, summarise); only the band lists are extended (20-100 Hz,
31.5 / 63 Hz octaves, thirds 25-80 Hz, 10 kHz third, 10-12.4 kHz, 1/24-octave decays up to 11.5 kHz), plus the L channel
(rv_metrics on the channel-swapped bands, L's tap 0 at K0 + 314), time-resolved IACC, a late coherence / tail-coloration window.

Both IRs are cut at the same index i1 per RT: the capture's file fade (and its ~-95 dB floor) sets it. Per metric the capture value
is compared with the spread of the model's realisations (seeds): inside if |model median - capture| <= max(K_SD x model sd,
an estimator floor per metric family). Pair metrics (EDC difference, 1/24-octave decay and narrow-band envelope deviations,
fine-structure correlation) compare model-vs-capture with model-vs-model (seed pairs).
"""
import argparse, itertools, json, os, re, sys, time

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.signal import butter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, HERE)
import pulse_check as pc  # noqa: E402
import e_model as E  # noqa: E402
import e_sim as S  # noqa: E402

FS, K0, LD = E.FS, E.K0, E.LDEL
K_SD = 2.5                                   # inside if |median - capture| <= max(K_SD x sd, floor)
EXTRA_OCT = [31.5, 63, "20-100", "10k-12.4k"]
EXTRA_THIRD = [25, 32, 40, 50, 63, 80, 10000]
EXT_BANDS = {"31.5", "63", "20-100", "10k-12.4k", "25", "32", "40", "50", "80", "10000"}   # outside the earlier 100 Hz-10 kHz scoring
# "63" and "10000" below: octave 63 / third 63 are both extended


def setup_bands():
    """Extend pulse_check's band lists in place (its functions read these module globals)."""
    sos = {31.5: butter(3, [31.5 / 2 ** .5, 31.5 * 2 ** .5], btype="band", fs=FS, output="sos"),
           63: butter(3, [63 / 2 ** .5, 63 * 2 ** .5], btype="band", fs=FS, output="sos"),
           "20-100": butter(3, [20, 100], btype="band", fs=FS, output="sos"),
           "10k-12.4k": butter(3, [10000, 12400], btype="band", fs=FS, output="sos")}
    for f, s in sos.items():
        pc.SOS[("o", f)] = s
    for f in EXTRA_THIRD:
        pc.SOS[("t", f)] = butter(3, [f / 2 ** (1 / 6), min(f * 2 ** (1 / 6), 11000)], btype="band", fs=FS, output="sos")
    pc.OCT = list(pc.OCT) + EXTRA_OCT
    pc.THIRD_REF = [100, 126, 159, 200, 252, 317, 400, 504, 635, 800, 1008, 1270, 1600, 2016, 2540, 3200, 4032, 5080, 6400]   # the earlier 100 Hz-6.3 kHz thirds
    pc.THIRD = [25, 32, 40, 50, 63, 80] + list(pc.THIRD) + [10000]
    eb = 100 * 2 ** (np.arange(0, 24 * np.log2(115) + 1) / 24)
    pc.EDR_BANDS = eb[eb <= 11500]


# ---------------------------------------------------------------- the end of the usable capture
def i1_of(rt):
    """Common cut index: the capture's file fade (knee, less 0.1 s) or where R's 50 ms level first reaches -90 dB (its floor
    is ~-95 dB), whichever is first; at most 5.2 s."""
    y = E.dev(rt)[:, 1]; w = int(0.05 * FS)
    lv = np.array([10 * np.log10((y[i:i + w] ** 2).mean() + 1e-30) for i in range(0, len(y) - w, w)])
    d = np.diff(lv)
    knee = next((j for j in range(4, len(d)) if d[j] < -5 and d[j] < 2.2 * np.mean(d[j - 4:j])), len(d) - 2)
    floor = next((j for j in range(2, len(lv)) if lv[j] < -90), len(lv))
    j = min(knee - 2, floor, int(5.2 / 0.05))
    return int(j * w)


# ---------------------------------------------------------------- per-IR metrics beyond compare()
def swap(X):
    return {k: v[:, ::-1] for k, v in X.items()}


def ring(t):
    t = np.array([np.nan if x is None else x for x in t], float)
    v = [t[i] / np.nanmedian(t[max(0, i - 4):i + 5]) for i in range(len(t)) if np.isfinite(t[i])]
    rough = [np.log(t[i] / np.nanmedian(t[max(0, i - 4):i + 5])) for i in range(len(t)) if np.isfinite(t[i])]
    return (float(np.nanmax(v)), float(np.std(rough))) if v else (np.nan, np.nan)


def derive_rv(rv, pre=""):
    """Scalars from an rv_metrics dict (+ derived: ringing, 1/24-octave roughness, echo density means)."""
    o = {}
    for k, v in rv.items():
        if k == "edr":
            if v:
                o[pre + "edr_ring"], o[pre + "edr_rough"] = ring(v)
        elif k == "ned1":
            a = np.array(v)
            for tag, lo, hi in (("a", 0.02, 0.1), ("b", 0.1, 0.3), ("c", 0.3, 1.0)):
                s = a[int(lo / 0.02):int(hi / 0.02)]
                if len(s):
                    o[pre + f"ned_{tag}"] = float(np.mean(s))
        elif v is not None:
            o[pre + k] = float(v)
    return o


def extras(X, i1, k0=K0):
    r = {}
    bb = X["bb"].astype(np.float64)
    for tag, fw in (("e", (0.05, 0.30)), ("l", (0.5, 1.0))):                      # tail coloration, both channels, early / late
        if i1 - k0 - int(fw[0] * FS) > int(0.15 * FS):
            Xs, f = pc.spectra(bb, i1, k0, fw); F = pc.fine(Xs, f)
            for lo, hi in pc.FS_BANDS:
                m_ = (f >= lo) & (f < hi)
                for c, ch in ((0, "L"), (1, "R")):
                    r[f"fsd{tag}_{lo}-{hi}_{ch}"] = float(np.std(F[m_, c]))
            if tag == "l":
                for fc, v in pc.coherence(Xs, f).items():
                    r[f"cohl_{fc}"] = v
    for t0, t1 in ((0.0, 0.05), (0.05, 0.15), (0.15, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)):   # IACC over time
        a, b = k0 + int(t0 * FS), min(k0 + int(t1 * FS), i1)
        if b - a > int(0.02 * FS):
            r[f"iacct_{int(t0 * 1000)}-{int(t1 * 1000)}"] = pc.iacc(bb[a:b, 0], bb[a:b, 1])
    nl, ne = int(pc.EARLY * FS), int(0.5 * FS)
    for f in pc.OCT:                                                                  # per octave, extended bands, L/R
        o = X[("o", f)].astype(np.float64)
        for tag, a, b in (("e", k0, k0 + nl), ("l", k0 + nl, min(i1, k0 + ne))):
            if b - a > int(0.02 * FS):
                r[f"iaccb_{tag}_{f}"] = pc.iacc(o[a:b, 0], o[a:b, 1])
                r[f"lrb_{tag}_{f}"] = float(10 * np.log10((np.sum(o[a:b, 0] ** 2) + 1e-30) / (np.sum(o[a:b, 1] ** 2) + 1e-30)))
    rv = pc.rv_metrics(swap(X), i1, k0 + LD)                                          # L as the "R" of rv_metrics
    for k in [k for k in rv if k.startswith(("iacc", "lr_"))]:
        rv.pop(k)
    r.update(derive_rv(rv, "L_"))
    for a_, b_w in ((0.03, 0.08), (0.08, 0.2), (0.2, 0.5), (0.5, 1.0)):              # broadband reverb level by window (dB re R's tap 1)
        a, b = k0 + int(a_ * FS), min(k0 + int(b_w * FS), i1)
        if b - a > int(0.02 * FS):
            for c, ch in ((0, "L"), (1, "R")):
                r[f"ew_{int(a_ * 1000)}-{int(b_w * 1000)}_{ch}"] = float(10 * np.log10(np.sum(bb[a:b, c] ** 2) + 1e-30))
    b_ = pc.band_decay(bb[:, 0], bb[:, 0], k0 + int(pc.T0 * FS), i1)                  # broadband decay (100 Hz-10 kHz), both channels
    r["bb_t_L"], r["bb_edt_L"] = (b_["t"] or [np.nan])[0], (b_["edt"] or [np.nan])[0]
    b_ = pc.band_decay(bb[:, 1], bb[:, 1], k0 + int(pc.T0 * FS), i1)
    r["bb_t_R"], r["bb_edt_R"] = (b_["t"] or [np.nan])[0], (b_["edt"] or [np.nan])[0]
    return r


def add_shape(o):
    """Tonal balance with the level taken out: each third-octave level (early / late) less the 100 Hz-8 kHz total, per kind."""
    for tag in ("e", "l"):
        for pre in ("", "L_"):
            ks = [k for k in o if k.startswith(f"{pre}lev_{tag}_")]
            ref = [k for k in ks if int(float(k.split("_")[-1])) in pc.THIRD_REF]
            if ref:
                tot = 10 * np.log10(sum(10 ** (o[k] / 10) for k in ref))
                for k in ks:
                    o[k.replace("lev_", "levs_", 1)] = o[k] - tot
    return o


def flat(m, side):
    """Absolute metrics of one side (0 = model a, 1 = capture b) from a JSON-roundtripped compare() result."""
    o = {}
    for band in ("oct", "third"):
        for fk, bs in m[band].items():
            for tag, b in (("L", bs[0]), ("R", bs[1])):
                for k in ("t", "edt", "late"):
                    v = b[k]
                    o[f"{band}_{k}_{fk}_{tag}"] = np.nan if v is None else v[side]
    o.update(derive_rv(m["rv"]["ab"[side]]))
    for fk, v in m.get("coh", {}).get("ab"[side], {}).items():
        o[f"coh_{fk}"] = v
    for w, x in m.get("fine_sd", {}).items():
        o[f"fsd_{w}_R"] = x[side]
    return o


PAIR_PREFIX = ("oct_dev", "third_dev", "rv_edr_p90", "rv_edr_max", "rv_nb_p90", "rv_nb_max", "fine_")


def pairkeys(m):
    s = pc.summarise(m)
    return {k: v for k, v in s.items() if k.startswith(PAIR_PREFIX) and v is not None and np.isfinite(v)}


def roundtrip(m):
    return json.loads(json.dumps(m, default=float))


# ---------------------------------------------------------------- tolerance per metric family
FLOORS = [  # (regex on the key without a leading L_, mode, floor)
    (r"^(oct|third)_t_", "ln", 0.03), (r"^(oct|third)_edt_", "ln", 0.05), (r"^(oct|third)_late_", "ln", 0.06),
    (r"^bb_t_", "ln", 0.03), (r"^bb_edt_", "ln", 0.05),
    (r"^rv_ts_", "ln", 0.03), (r"^rv_mix_ms", "ln", 0.10), (r"^mix_ms", "ln", 0.10),
    (r"^(rv_)?c(50|80)_", "lin", 0.3), (r"^(rv_)?d50_", "lin", 0.01), (r"^(rv_)?lev_", "lin", 0.3),
    (r"^fsd", "lin", 0.15), (r"^coh", "lin", 0.05), (r"^iacc", "lin", 0.03), (r"^rv_iacc", "lin", 0.03),
    (r"^lrb_|^rv_lr_", "lin", 0.3), (r"^(rv_)?ned_", "lin", 0.02), (r"edr_ring", "lin", 0.03), (r"edr_rough", "lin", 0.01),
]


def mode_floor(key):
    k = key[2:] if key.startswith("L_") else key
    for rx, mode, fl in FLOORS:
        if re.search(rx, k):
            return mode, fl
    return "lin", 0.0


def family(key):
    k = key[2:] if key.startswith("L_") else key
    ch = "L" if key.startswith("L_") else ""
    m = re.match(r"^(oct|third)_(t|edt|late)_(.+)_([LR])$", k)
    if m:
        band, kind, f, c = m.groups()
        return f"{band} {'T60' if kind == 't' else kind.upper()} {c}", f
    m = re.match(r"^(iaccb|lrb)_([el])_(.+)$", k)
    if m:
        return f"{m.group(1)} {m.group(2)}", m.group(3)
    m = re.match(r"^ew_(.+)_([LR])$", k)
    if m:
        return f"ew {m.group(2)}", m.group(1)
    m = re.match(r"^(rv_)?(c50|c80|d50|ts|lev_e|lev_l|levs_e|levs_l)_(.+)$", k)
    if m:
        return f"{m.group(2)}{' L' if ch else ' R'}", m.group(3)
    m = re.match(r"^(fsd[el])_(.+)_([LR])$", k)
    if m:
        return f"{m.group(1)} {m.group(3)}", m.group(2)
    m = re.match(r"^(cohl|coh)_(.+)$", k)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r"^iacct_(.+)$", k)
    if m:
        return "iacct", m.group(1)
    return k.replace("rv_", "") + (" L" if ch else ""), ""


def tr(x, mode):
    return np.log(x) if mode == "ln" else x


def judge(vals, cap, key):
    """vals: model values over seeds; cap: capture value. Returns dict (None if not computable)."""
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], float)
    if len(v) < 3 or cap is None or not np.isfinite(cap):
        return None
    mode, fl = mode_floor(key)
    if mode == "ln" and (cap <= 0 or np.any(v <= 0)):
        return None
    a, c = tr(v, mode), tr(cap, mode)
    med, sd = float(np.median(a)), float(np.std(a, ddof=1))
    tol = max(K_SD * sd, fl)
    off = med - c
    res = dict(cap=float(cap), med=float(np.median(v)), sd_model=float(np.std(v, ddof=1)), lo=float(v.min()), hi=float(v.max()),
               off=float(np.expm1(off) if mode == "ln" else off), off_mode="ratio-1" if mode == "ln" else "diff",
               tol=float(np.expm1(tol) if mode == "ln" else tol), z=float(off / (sd + 1e-12)), excess=float(abs(off) - tol),
               inside=bool(abs(off) <= tol), n=int(len(v)))
    return res


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


# ---------------------------------------------------------------- runs
def render(kind, rt, seed, i1, kern):
    if kind == "imp":
        y = S.render_ir(rt, seed=seed, secs=(i1 + int(0.2 * FS)) / FS)
    else:                                                    # 'sweep' = 10 s; 'sweep3' / 'sweep5' ...: that many seconds
        prm = dict(S.DEF) | S.ADOPTED | dict(tsw=float(kind[5:] or 10), seed=seed)
        y = S.capture(rt, prm, kern, tail=(i1 - K0 + int(0.2 * FS)) / FS)
    return y / y[K0, 1]


class Gain:
    """Scale the model's reverb (static reads c_k and the moving reads' amplitudes, not the direct path) by g, for a diagnostic
    run in which the reverb level law is taken as right: what is left outside is then not a level error."""
    def __init__(self, g):
        self.g = g

    def __enter__(self):
        self.sp, self.fr, self.xs = S.static_params, S.fixed_reads, S.X_STATIC
        q = S.X_STATIC.copy(); q[3:6] *= self.g
        S.X_STATIC = q; S.static_params = lambda: q
        g = self.g; fr = self.fr
        S.fixed_reads = lambda: [dict(r, R=r["R"] * g, L=r["L"] * g) for r in fr()]
        return self

    def __exit__(self, *a):
        S.static_params, S.fixed_reads, S.X_STATIC = self.sp, self.fr, self.xs


def run_kind(kind, rt, seeds, npairs, i1, n, kern, Xc, capm, g, log):
    """Model seeds of one kind ('imp' / 'sweep', '_g' = reverb scaled by g) against the capture. Returns (result, capm)."""
    base = kind.split("_")[0]
    t0 = time.time(); X, flats, pairs_mc = [], [], []
    for s in seeds:
        if g is not None:
            with Gain(g):
                y = render(base, rt, (0 if base == "imp" else 100) + s, i1, kern)
        else:
            y = render(base, rt, (0 if base == "imp" else 100) + s, i1, kern)
        Xm = pc.prep(y[:n], None, np.float32); X.append(Xm)
        m = roundtrip(pc.compare(Xm, Xc, i1, k0=K0))
        flats.append(add_shape(flat(m, 0) | extras(Xm, i1)))
        pairs_mc.append(pairkeys(m))
        if capm is None:
            capm = add_shape(flat(m, 1) | extras(Xc, i1))
        if s == 0:
            log(f"  RT {rt} {kind} seed 0: {time.time() - t0:.0f} s")
    allp = list(itertools.combinations(range(len(seeds)), 2))
    allp = [allp[i] for i in np.linspace(0, len(allp) - 1, min(npairs, len(allp))).astype(int)] if npairs else []
    pairs_mm = [pairkeys(roundtrip(pc.compare(X[i], X[j], i1, k0=K0))) for i, j in allp]
    log(f"  RT {rt} {kind}: {len(seeds)} seeds, {len(allp)} model pairs, {time.time() - t0:.0f} s")
    return dict(model=flats, pairs_mc=pairs_mc, pairs_mm=pairs_mm), capm


def run_rt(rt, kinds, seeds, npairs, log, gain_kinds=True):
    i1 = i1_of(rt); n = i1 + int(0.1 * FS)
    kern = E.kernels([r for r in E.RTS if r >= 1.0])
    cap = E.dev(rt)[:n]
    Xc = pc.prep(cap, None, np.float32)
    capm = None; res = {}; gs = {}
    a = K0 + int(0.08 * FS)
    for kind in kinds:
        res[kind], capm = run_kind(kind, rt, seeds, npairs, i1, n, kern, Xc, capm, None, log)
        if gain_kinds:
            key = "ew_80-200_R"
            # level correction from the model seeds' own windows (R, 0.08 s to i1): g = sqrt(capture / model energy)
            Em = np.mean([10 ** (f["ew_80-200_R"] / 10) + 10 ** (f.get("ew_200-500_R", -99) / 10) + 10 ** (f.get("ew_500-1000_R", -99) / 10)
                          for f in res[kind]["model"]])
            Ec = 10 ** (capm["ew_80-200_R"] / 10) + 10 ** (capm.get("ew_200-500_R", -99) / 10) + 10 ** (capm.get("ew_500-1000_R", -99) / 10)
            gs[kind] = float(np.sqrt(Ec / Em))
            res[kind + "_g"], _ = run_kind(kind + "_g", rt, seeds, npairs, i1, n, kern, Xc, capm, gs[kind], log)
    return dict(i1=i1, t_end=(i1 - K0) / FS, cap=capm, gains=gs, kinds=res)


def judge_all(R):
    """Per kind, per RT, per metric key."""
    out = {}
    for kind in next(iter(R.values()))["kinds"]:
        out[kind] = {}
        for rt, d in R.items():
            k = d["kinds"][kind]; cap = d["cap"]; J = {}
            keys = sorted(set().union(*[f.keys() for f in k["model"]]) & set(cap.keys()))
            for key in keys:
                j = judge([f.get(key) for f in k["model"]], cap.get(key), key)
                if j:
                    J[key] = j
            # pair metrics: model-vs-capture (median over seeds) against model-vs-model (all seed pairs)
            P = {}
            pk = set().union(*[p.keys() for p in k["pairs_mc"]])
            for key in sorted(pk):
                mc = [p[key] for p in k["pairs_mc"] if key in p]
                mm = [p[key] for p in k["pairs_mm"] if key in p]
                if len(mc) >= 3 and len(mm) >= 3:
                    mcm, mmq = float(np.median(mc)), float(np.max(mm))
                    P[key] = dict(mc_med=mcm, mc_max=float(np.max(mc)), mm_med=float(np.median(mm)), mm_max=mmq, mm_p90=float(np.percentile(mm, 90)),
                                  inside=bool(mcm <= mmq * 1.1 + (0.02 if key.startswith("fine") else 0.3 if "dev" in key or "nb" in key else 0.02)))
            out[kind][str(rt)] = dict(metrics=J, pairs=P)
    return out


def summarise_all(J, rts_main):
    """Per kind and metric key: RTs >= 1.0 inside / n, median offset, sign counts of the outside ones."""
    S_ = {}
    for kind, byrt in J.items():
        keys = sorted(set().union(*[v["metrics"].keys() for v in byrt.values()]))
        rows = {}
        for key in keys:
            js = [(rt, byrt[rt]["metrics"].get(key)) for rt in byrt if float(rt) in rts_main]
            js = [(rt, j) for rt, j in js if j]
            if not js:
                continue
            offs = np.array([j["off"] for _, j in js]); ins = [j["inside"] for _, j in js]
            outs = [j for _, j in js if not j["inside"]]
            fam, band = family(key)
            rows[key] = dict(family=fam, band=band, ext=band in EXT_BANDS, n=len(js), inside=int(sum(ins)), med_off=float(np.median(offs)),
                             n_pos=int(sum(1 for j in outs if j["off"] > 0)), n_neg=int(sum(1 for j in outs if j["off"] < 0)),
                             max_excess=float(max([j["excess"] for _, j in js])), mode=js[0][1]["off_mode"], rts_out=[float(rt) for rt, j in js if not j["inside"]])
        fams = {}
        for key, r in rows.items():
            fn = r["family"] + (" [ext]" if r["ext"] else "")
            f = fams.setdefault(fn, dict(n=0, inside=0, keys_out=[]))
            f["n"] += r["n"]; f["inside"] += r["inside"]
            if r["inside"] < r["n"]:
                f["keys_out"].append((key, r["inside"], r["n"], round(r["med_off"], 4), r["n_pos"], r["n_neg"]))
        S_[kind] = dict(by_key=rows, by_family=fams)
    return S_


PARTS = os.path.join(ROOT, "analysis", "e250", "e_pulse_parts")


def part_path(rt):
    return os.path.join(PARTS, f"rt_{rt}.json")


def stage_run(a, log):
    setup_bands(); S.static_params()
    os.makedirs(PARTS, exist_ok=True)
    for rt in a.rts:
        if os.path.exists(part_path(rt)) and not a.redo:
            log(f"RT {rt}: part exists, skipped"); continue
        t0 = time.time(); r = run_rt(rt, a.kinds, range(a.seeds), a.pairs, log, not a.no_gain)
        json.dump(clean(r), open(part_path(rt), "w"))
        log(f"RT {rt}: i1 {r['t_end']:.2f} s, {time.time() - t0:.0f} s, gains {r['gains']}")


def stage_report(a, log):
    R = {}
    for rt in a.rts:
        if os.path.exists(part_path(rt)):
            R[rt] = json.load(open(part_path(rt)))
    J = judge_all(R)
    SM = summarise_all(J, set(r for r in R if r >= 1.0))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(clean(dict(meta=dict(rts=sorted(R), seeds=a.seeds, k_sd=K_SD, note=__doc__.split("\n")[0],
                                   t_end={str(r): R[r]["t_end"] for r in R}, gains={str(r): R[r]["gains"] for r in R}),
                         capture={str(r): R[r]["cap"] for r in R},
                         model={str(r): {k: v["model"] for k, v in R[r]["kinds"].items()} for r in R},
                         judged=J, summary=SM)), open(a.out, "w"))
    for kind, s in SM.items():
        log(f"\n== {kind}: families (RT >= 1.0), inside / n")
        for fn, f in sorted(s["by_family"].items()):
            log(f"  {fn:28s} {f['inside']:4d} / {f['n']:4d}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="run", choices=["run", "report"])
    ap.add_argument("--rts", type=float, nargs="*", default=E.RTS)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--pairs", type=int, default=6, help="model seed pairs per (RT, kind)")
    ap.add_argument("--kinds", nargs="*", default=["imp", "sweep"])
    ap.add_argument("--no-gain", action="store_true", help="skip the level-matched diagnostic kinds (imp_g, sweep_g)")
    ap.add_argument("--redo", action="store_true", help="recompute RTs whose part file exists")
    ap.add_argument("--out", default=os.path.join(ROOT, "analysis", "e250", "e_pulse.json"))
    ap.add_argument("--parts", default=None, help="folder of the per-RT part files")
    a = ap.parse_args()
    global PARTS
    PARTS = a.parts or PARTS
    log = lambda s: print(s, flush=True)
    if a.stage == "run":
        stage_run(a, log)
    else:
        setup_bands()
        stage_report(a, log)


if __name__ == "__main__":
    main()
