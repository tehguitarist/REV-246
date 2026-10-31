"""The pulse check's per-IR metrics (tools/pulse_check.py G2-G6, plus G1 nulls) for program N: the PLUGIN's N renders
against the Tomkin 96 kHz Non-Lin captures, beside the N bench model against the captures, beside a tolerance.
Results: analysis/plugin_pulse/n/, report analysis/plugin_pulse/report_n.log, write-up analysis/plugin_check/REPORT_chainfix.md §9.

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/plugin_pulse_n.py --stage run [--only k,..] [--limit N] [--workers 1]
    .venv/bin/python tools/plugin_pulse_n.py --stage report

Inputs (nothing re-rendered, nothing refitted against the plugin):
  - renders: plugin_check's corrected-bridge set, analysis/plugin_check/renders_fix/96000/N_*.wav (30 October build),
    put on nl_structure's device grid by plugin_check.n_device (fixed anchor), exactly as plugin_check's n stage;
  - captures: plugin_check.n_keys() (N_EVERY = 4: 261 of the 1,044 clean captures), the library's 96 kHz-sample slip
    undone as analysis/nl/branch_slip/n_deslip.py does (the part after a cut in the front / branch gap moved back by n
    file samples; n from tools/defects.NL_SLIPPED, checked against the bench's o);
  - the plugin's fit per capture (front shift measured, then shift, clock and 2 gains, plugin_check._n_job's fit): read
    from plugin_check's own result files (analysis/plugin_check/n_fix/ for n = 0, analysis/nl/branch_slip/deslip/n_fix/
    for the de-slipped ones) and re-checked here (the window nulls are recomputed and must equal the stored ones);
  - the bench model: nl_model.render with the adopted bench's own per-capture numbers (bench_warp2_bk_kf_a0.7_lb2_hf_g5:
    o, clock, front shift; gains by least squares as the bench), against the same de-slipped capture clock-corrected as
    the bench does it; o moved by the slip (o - n / 3.64), so plugin and bench see the same capture.
Skipped / flagged: SPAN 0 / LATE 0 (no branch: three front taps only, skipped), captures whose bench o is not +1 plus a
whole-number slip in -400..-1 (Tomkin's trims, §11.4: flagged "unresolved", out of the medians), trimmed starts (flagged,
reported apart).

Tolerance (the library has no repeat captures of an N setting: all 1,046 labels are distinct):
  - twin: the capture against itself plus an independent stretch of the library's own capture noise (from the LATE 0
    tails, as nl_model's noise stage), 3 realisations per capture; the tolerance is the largest |twin error|. It is what a
    perfect, noise-free model would score against the capture (a real repeat would differ by sqrt(2) x that noise and by
    whatever else changes between captures, so it is a lower bound on a repeat spread);
  - LATE pair (real capture vs capture, early section only): the sampled capture against the library capture that
    differs only in LATE (identical early section up to the late's first arrival, NonLin §10.1), SPAN >= 0.2, aligned with
    one shift and 2 gains (the plugin's freedoms less the clock); early metrics (0-80 ms third-octave levels, early IACC
    and L/R) and the 0-50 ms null.

Metric definitions are pulse_check.compare's (EDC from the front frame's zero, t0 = 0: N has no cluster to skip), plus two
for N: the capture's EDC linearity per band and fit range (max |dev| from the fitted line, dB: a decay time is reported only
where the capture's EDC is a straight line within 1 dB over that range) and the tail decay (the late section's RT: EDC
from 40 ms after the last main tap, -5/-25 dB, LATE > 0).
"""
import argparse, glob, json, os, sys, time, zlib

os.environ.setdefault("PLUGIN_CHECK_TAG", "fix")
os.environ.setdefault("NL_BK", "1")                  # the adopted bench's branch kernel (bench_warp2_bk_kf_a0.7_...)
for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import plugin_check as pk  # noqa: E402  (sets EMT_KERNEL=chain3)
import pulse_check as pc  # noqa: E402
import numpy as np  # noqa: E402
from scipy.signal import butter, sosfiltfilt  # noqa: E402
import defects  # noqa: E402
import nl_structure as nl  # noqa: E402
import nl_early as ne  # noqa: E402
import nl_model as nm  # noqa: E402
import nl_late as nlt  # noqa: E402
import ml_lines as mll  # noqa: E402

assert pk.TAG == "fix" and nm.BENCH == pk.N_BENCH, (pk.TAG, nm.BENCH)
FS, K0 = nl.FSD, nl.K0
STEP = 96000 / FS
OUT = os.path.join(ROOT, "analysis", "plugin_pulse", "n")
LOG = os.path.join(ROOT, "analysis", "plugin_pulse", "report_n.log")
FIT_N0 = os.path.join(pk.OUT, "n_fix")                                       # plugin_check n stage (PLUGIN_CHECK_TAG=fix)
FIT_DS = os.path.join(ROOT, "analysis", "nl", "branch_slip", "deslip", "n_fix")  # n_deslip.py (slip undone)
OCT, THIRD = pc.OCT, pc.THIRD
RANGES = {"t": (-5, -25), "edt": (0, -10), "late": (-25, -40)}
LIN_MAX = 1.0                                         # dB: capture EDC straight within this over the fit range
TAIL_GUARD = 0.04                                     # s after the last main tap (AP(180) tail 3 dB / 180 samples)
N_TWIN = 3


def kstr(k):
    return "_".join(str(v) for v in k)


def bench(key):
    return json.load(open(os.path.join(pk.N_BENCH, kstr(key) + ".json")))


def slip(key):
    """(n file samples, status): n from defects.NL_SLIPPED, checked against the bench's o."""
    b = bench(key)
    nb = int(round((b["o"] - 1) * STEP))
    if key in defects.NL_SLIPPED:
        assert defects.NL_SLIPPED[key] == nb, (key, defects.NL_SLIPPED[key], nb)
        return nb, "slipped"
    if nb == 0:
        return 0, "ok"
    return nb, "unresolved"                            # trims (§11.4) and the off-grid offsets: not undone


def deslip(y, fr, n):
    """n_deslip.ir_deslip: the part after the cut in the front / branch gap moved later by -n file samples."""
    if n == 0:
        return y
    s = -n / STEP
    cut = K0 + int(round(fr + (466 + 588 - s) / 2))
    post = np.zeros_like(y); post[cut:] = y[cut:]
    out = y.copy(); out[cut:] = 0
    return out + mll.warp(post, s, 0.0)


# ---------------------------------------------------------------- the capture noise pool (twins)
_POOL = []


def noise_pool():
    """The library's capture noise, full band, both channels: LATE 0 tails from 0.25 s after the last main tap to the usable
    end (nl_model.noise_stage's segments), joined with 5 ms crossfades."""
    if _POOL:
        return _POOL[0]
    segs = []
    for k, p in sorted(nl.table().items()):
        if k[2] != 0.0 or k in defects.NL_DEFECTS:
            continue
        raw = nl.ir(p)
        h = sosfiltfilt(ne.BAND, raw, axis=0)
        ns = ne.N_STAGES[k[0]]
        i0 = K0 + (ne.taps(0, ns)[-1][0] if ns else 600) + int(0.25 * FS)
        i1 = nm.usable_end(h)
        if i1 - i0 >= int(0.1 * FS):
            segs.append(raw[i0:i1])
    xf = int(0.005 * FS)
    w = np.sin(np.linspace(0, np.pi / 2, xf)) ** 2
    out = segs[0]
    for s in segs[1:]:
        out = np.vstack([out[:-xf], out[-xf:] * (1 - w)[:, None] + s[:xf] * w[:, None], s[xf:]])
    _POOL.append(out)
    return out


# ---------------------------------------------------------------- N extras
def edc_lin(X, k0, i1):
    """Per band and channel: max |EDC - its fitted line| (dB) over each fit range (points before i1 - 150 ms)."""
    out = {}
    for (kind, f) in pc.SOS:
        for c in (0, 1):
            d = pc.edc_curve(X[(kind, f)][:, c], k0, i1)
            n_ok = max(0, i1 - k0 - int(0.15 * FS))
            for nm_, (hi, lo) in RANGES.items():
                m = (d <= hi) & (d >= lo)
                m[n_ok:] = False
                if m.sum() < int(0.03 * FS):
                    v = None
                else:
                    t = np.nonzero(m)[0] / FS
                    p = np.polyfit(t, d[m], 1)
                    v = float(np.max(np.abs(d[m] - np.polyval(p, t))))
                out[f"{'oct' if kind == 'o' else 'third'}_{nm_}{f}{'' if c else '_L'}"] = v
    return out


def tail_decay(X, ta, i1):
    """The late section's decay (s / 60 dB) from the EDC after ta, -5 to -25 dB, per octave and broadband, R and L."""
    out = {}
    if i1 - ta < int(0.1 * FS):
        return out
    n_ok = max(0, i1 - ta - int(0.15 * FS))
    for name, key in [("bb", "bb")] + [(str(f), ("o", f)) for f in OCT]:
        for c in (0, 1):
            d = pc.edc_curve(X[key][:, c], ta, i1)
            t = pc.edc_time(d, -5, -25, n_ok)
            out[f"tail_{name}{'' if c else '_L'}"] = t
    return out


def g1_nulls(y, m, i0, n, fr, noise, mb=None, yb=None):
    """nl_model's window nulls (100 Hz-10 kHz, from i0) and the excess over the noise floor, both channels. mb: the banded
    model if it was banded before the warp (plugin_check._n_job's order)."""
    yb = sosfiltfilt(ne.BAND, y, axis=0) if yb is None else yb
    mb = sosfiltfilt(ne.BAND, m, axis=0) if mb is None else mb
    out = {}
    for t0, t1 in nm.WINDOWS:
        a_, b_ = max(i0, K0 + int(fr) + int(t0 * FS)), min(n, K0 + int(fr) + int(t1 * FS))
        if b_ - a_ < int(0.02 * FS) or np.sum(yb[a_:b_, 0] ** 2) < 10 * noise["ch0"] * (b_ - a_):
            continue
        out[f"{t0}-{t1}"] = {f"{k}{c}": v for c in (0, 1) for k, v in (
            ("null", float(10 * np.log10(np.sum((yb[a_:b_, c] - mb[a_:b_, c]) ** 2) / np.sum(yb[a_:b_, c] ** 2)))),
            ("floor", float(10 * np.log10(noise[f"ch{c}"] * (b_ - a_) / np.sum(yb[a_:b_, c] ** 2)))))}
    a_, b_ = max(i0, K0 + int(fr) + int(0.05 * FS)), min(n, K0 + int(fr) + int(0.6 * FS))
    bands = {}
    if b_ - a_ > int(0.05 * FS):
        for lo, hi in nm.BANDS:
            sos = butter(4, [lo, hi], btype="band", fs=FS, output="sos")
            Y, M = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, m, axis=0)
            bands[f"{lo}-{hi}"] = [float(10 * np.log10(np.sum((Y[a_:b_, c] - M[a_:b_, c]) ** 2) / np.sum(Y[a_:b_, c] ** 2))) for c in (0, 1)]
    return out, bands


def signal_end(h, noise, n):
    """Last index where the capture's 100 Hz-10 kHz 10 ms energy (either channel) is 10 dB above the noise floor."""
    w = int(0.01 * FS)
    e = np.stack([np.convolve(h[:n, c] ** 2, np.ones(w) / w, "same") for c in (0, 1)], 1)
    ok = (e[:, 0] > 10 * noise["ch0"]) | (e[:, 1] > 10 * noise["ch1"])
    idx = np.nonzero(ok)[0]
    return int(min(n, idx[-1] + w)) if len(idx) else n


def early_only(r):
    return {k: v for k, v in r.items() if k.startswith(("lev_e_", "iacc_e", "lr_e"))}


def late_pair(key, fr1, y1, i0):
    """The capture against the library capture differing only in LATE (SPAN >= 0.2; early metrics and 0-50 ms null)."""
    span, early, late = key
    if span < 0.2:
        return None
    cands = sorted((abs(k[2] - late), k) for k in nl.table() if k[:2] == key[:2] and k != key and k not in defects.NL_DEFECTS)
    for _, k2 in cands:
        n2, st2 = slip(k2)
        b2 = bench(k2)
        if st2 == "unresolved" or b2["trimmed"]:
            continue
        raw2 = nl.ir(k2)
        h2 = sosfiltfilt(ne.BAND, raw2, axis=0)
        fr2 = ne.front_shift(h2)
        y2 = deslip(raw2, fr2, n2)
        nn = min(len(y1), len(y2))
        y2 = y2[:nn]                                   # both in the capture clock (|d ppm| ~ 1-3: < 0.01 sample over 120 ms)
        a_, b_ = max(i0, K0 + int(fr1)), K0 + int(fr1) + int(0.05 * FS)
        y1b = sosfiltfilt(ne.BAND, y1[:nn], axis=0)
        y2b = sosfiltfilt(ne.BAND, y2, axis=0)
        from scipy.optimize import minimize_scalar

        def cost(s):
            z = mll.warp(y2b[:b_ + 200], s, 0.0)[a_:b_]
            g = [(y1b[a_:b_, c] @ z[:, c]) / (z[:, c] @ z[:, c]) for c in (0, 1)]
            return float(sum(np.log(np.sum((y1b[a_:b_, c] - g[c] * z[:, c]) ** 2)) for c in (0, 1)))
        s0 = fr1 - fr2
        r = minimize_scalar(cost, bounds=(s0 - 0.5, s0 + 0.5), method="bounded", options=dict(xatol=1e-4))
        z = mll.warp(y2, r.x, 0.0)
        zb = sosfiltfilt(ne.BAND, z, axis=0)
        g = np.array([(y1b[a_:b_, c] @ zb[a_:b_, c]) / (zb[a_:b_, c] @ zb[a_:b_, c]) for c in (0, 1)])
        null = [float(10 * np.log10(np.sum((y1b[a_:b_, c] - g[c] * zb[a_:b_, c]) ** 2) / np.sum(y1b[a_:b_, c] ** 2))) for c in (0, 1)]
        k0 = K0 + int(round(fr1))
        i1 = min(nn, k0 + int(0.12 * FS))           # the late's first arrival is >= 150 ms after the input at SPAN >= 0.2
        A, B = pc.prep((z * g)[:i1 + 2000]), pc.prep(y1[:i1 + 2000])
        ra, rb = early_only(pc.rv_metrics(A, i1, k0)), early_only(pc.rv_metrics(B, i1, k0))
        return dict(partner=list(k2), shift=float(r.x), gain=g.tolist(), null_0_50=null, a=ra, b=rb)
    return None


# ---------------------------------------------------------------- one capture
def job(key):
    path = os.path.join(OUT, kstr(key) + ".json")
    if os.path.exists(path):
        return key
    t = time.time()
    span, early, late = key
    noise = json.load(open(nm.NOISE))
    nsl, status = slip(key)
    b = bench(key)
    raw0 = nl.ir(nl.table()[key])
    raw = deslip(raw0, b["front_shift"], nsl) if status == "slipped" else raw0   # n_deslip: the cut from the bench's front shift
    h = sosfiltfilt(ne.BAND, raw, axis=0)
    fr = ne.front_shift(h)                                                      # then _n_job's own, on the de-slipped capture
    n = min(nm.usable_end(h), K0 + int(fr) + int(nm.WINDOWS[-1][1] * FS))
    i0, trimmed = nm.window_start(h, fr)
    out = dict(key=list(key), slip=nsl, status=status, trimmed=trimmed, front_shift=fr, n=n, bench_o=b["o"])
    # the plugin in the capture's frame: plugin_check._n_job's fit, from its result files
    ff = os.path.join(FIT_DS if status == "slipped" else FIT_N0, kstr(key) + ".json")
    fit = json.load(open(ff))
    assert abs(fit["front_shift"] - fr) < 1e-5, (key, fit["front_shift"], fr)
    fp = fit["plugin"]
    pr = pk.n_device(pk.rpath(pk.HOST, "N", kstr(key), span))
    m0 = np.vstack([pr, np.zeros((max(0, n + 400 - len(pr)), 2))])[:n + 400]
    P = mll.warp(m0, fp["shift"], fp["ppm"])[:n] * np.array(fp["gain"])
    Y = raw[:n] if len(raw) >= n else np.vstack([raw, np.zeros((n - len(raw), 2))])
    Pb = mll.warp(sosfiltfilt(ne.BAND, m0, axis=0), fp["shift"], fp["ppm"])[:n] * np.array(fp["gain"])
    wp, bp = g1_nulls(Y, P, i0, n, fr, noise, Pb, h[:n])
    chk = max([abs(wp[w][f"null{c}"] - fp["windows"][w][f"null{c}"]) for w in fp["windows"] for c in (0, 1)] or [0.0])
    out["plugin"] = dict(fit=dict(shift=fp["shift"], ppm=fp["ppm"], gain=fp["gain"]), windows=wp, bands=bp, check_db=chk)
    # the bench model against the same capture, clock-corrected as the bench does
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    Yb = nlt.unwarp_from(nlt.upsample(rawp), n, b["ppm"])[:n]
    ob = b["o"] - nsl / STEP if status == "slipped" else b["o"]
    mf = nm.render(key, n, ob, b["front_shift"], band=False)
    hb, mbb = sosfiltfilt(ne.BAND, Yb, axis=0), sosfiltfilt(ne.BAND, mf, axis=0)
    g = np.array([(hb[:, c] @ mbb[:, c]) / (mbb[:, c] @ mbb[:, c]) for c in (0, 1)])
    Bm = mf * g
    wb, bb = g1_nulls(Yb, Bm, i0, n, fr, noise)
    chkb = max([abs(wb[w][f"null{c}"] - b["windows"][w][f"null{c}"]) for w in b["windows"] if w in wb for c in (0, 1)] or [0.0]) if status == "ok" else None
    out["bench"] = dict(o=ob, ppm=b["ppm"], gain=g.tolist(), windows=wb, bands=bb, check_db=chkb)
    # metrics (pulse_check.compare), up to where the capture is 10 dB above the noise
    i1 = signal_end(h, noise, n)
    k0 = K0 + int(round(fr))
    out.update(i1=i1, k0=k0, t_end=(i1 - k0) / FS)
    fw = pc.FS_WIN
    YP = pc.prep(Y)
    out["lin"] = edc_lin(YP, k0, i1)
    ns = ne.N_STAGES[span]
    ta = k0 + (ne.taps(0, ns)[-1][0] + 1 if ns else 600) + int(TAIL_GUARD * FS)
    out["tail_capture"] = tail_decay(YP, ta, i1) if late > 0 else {}
    PP = pc.prep(P)
    out["plugin"]["m"] = pc.compare(PP, YP, i1, k0, 0.0, fw)
    out["plugin"]["tail"] = tail_decay(PP, ta, i1) if late > 0 else {}
    del PP
    YB, BB = pc.prep(Yb), pc.prep(Bm)
    out["bench"]["m"] = pc.compare(BB, YB, i1, k0, 0.0, fw)
    out["bench"]["tail"] = tail_decay(BB, ta, i1) if late > 0 else {}
    out["tail_capture_bench"] = tail_decay(YB, ta, i1) if late > 0 else {}
    del YB, BB
    pool = noise_pool()
    rng = np.random.default_rng(zlib.crc32(kstr(key).encode()))
    out["twin"] = []
    for _ in range(N_TWIN):
        o_ = int(rng.integers(0, len(pool) - n))
        T = Y + pool[o_:o_ + n]
        TP = pc.prep(T)
        out["twin"].append(dict(offset=o_, m=pc.compare(TP, YP, i1, k0, 0.0, fw), tail=tail_decay(TP, ta, i1) if late > 0 else {}))
        del TP
    out["pair"] = late_pair(key, fr, raw, i0) if not trimmed and status != "unresolved" else None
    out["sec"] = time.time() - t
    json.dump(out, open(path, "w"))
    w = out["plugin"]["windows"]
    print(f"{kstr(key)} [{status}{' trimmed' if trimmed else ''}] {out['sec']:.0f} s check {chk:.3f}/{chkb if chkb is None else round(chkb, 3)} | plugin "
          + " ".join(f"{k}:{v['null1']:.1f}" for k, v in w.items()) + " | bench " + " ".join(f"{v['null1']:.1f}" for v in wb.values()), flush=True)
    return key


def keys_all():
    return [k for k in pk.n_keys() if not (k[0] == 0.0 and k[2] == 0.0)]          # SPAN 0 / LATE 0: no branch, skipped


def run_stage(workers, only=None, limit=None):
    os.makedirs(OUT, exist_ok=True)
    keys = keys_all()
    if only:
        keys = [k for k in keys if kstr(k) in only]
    keys = [k for k in keys if not os.path.exists(os.path.join(OUT, kstr(k) + ".json"))]
    if limit:
        keys = keys[:limit]
    print(f"N: {len(keys)} captures to run", flush=True)
    if workers == 1:
        for k in keys:
            job(k)
        return
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(workers, max_tasks_per_child=10) as ex:
        for _ in ex.map(job, keys):
            pass


# ---------------------------------------------------------------- report
def _summ(m, lin, tail, tail_ref):
    """pulse_check.summarise, decay fits dropped where the capture's EDC is not straight (LIN_MAX), plus the tail decay."""
    s = pc.summarise(m)
    for band, fs_ in (("oct", OCT), ("third", THIRD)):
        for f in fs_:
            for tag in ("", "_L"):
                for nm_ in RANGES:
                    k = f"{band}_{nm_}{f}{tag}"
                    lv = lin.get(k)
                    if lv is None or lv > LIN_MAX:
                        s[f"{band}_{nm_}{f}{tag}"] = None
    for k, v in tail.items():
        r = tail_ref.get(k)
        s[k] = None if v is None or r is None else v / r - 1
    for k in [x for x in s if x.startswith("rv_edr") or x.startswith("rv_ring")]:
        s[k] = None                                    # EDR: per 1/24-octave -5/-25 dB fits, not meaningful on N's envelopes
    return s


def _md(rr, k, f=abs):
    x = [f(r[k]) for r in rr if r.get(k) is not None and np.isfinite(r[k])]
    return (float(np.median(x)), len(x)) if x else (np.nan, 0)


def load():
    D = []
    for fn in sorted(glob.glob(os.path.join(OUT, "*.json"))):
        d = json.load(open(fn))
        tc = d["tail_capture"]
        d["S_p"] = _summ(d["plugin"]["m"], d["lin"], d["plugin"].get("tail", {}), tc)
        d["S_b"] = _summ(d["bench"]["m"], d["lin"], d["bench"].get("tail", {}), d.get("tail_capture_bench", tc))
        d["S_t"] = [_summ(t["m"], d["lin"], t.get("tail", {}), tc) for t in d["twin"]]
        D.append(d)
    return D


def report_stage():
    D = load()
    span = lambda d: d["key"][0]
    early = lambda d: d["key"][1]
    late = lambda d: d["key"][2]
    for d in D:                    # the bench itself fails (a capture that is not its label's response): listed, out of the medians
        w = d["bench"]["windows"]
        d["anomalous"] = w.get("0.05-0.15", {}).get("null1", -99) > -20     # gross: not its label's response
    main = [d for d in D if d["status"] != "unresolved" and not d["trimmed"] and not d["anomalous"]]
    groups = [("All untrimmed", main),
              ("EARLY > 0 (decaying ramp, < 10)", [d for d in main if 0 < early(d) < 10]),
              ("EARLY 10 (flat)", [d for d in main if early(d) == 10]),
              ("EARLY < 0 (reverse)", [d for d in main if early(d) < 0]),
              ("LATE 0 (gated)", [d for d in main if late(d) == 0]),
              ("LATE > 0", [d for d in main if late(d) > 0]),
              ("SPAN 0 (front + late)", [d for d in main if span(d) == 0]),
              ("LATE 0.2", [d for d in main if late(d) == 0.2]),
              ("LATE >= 0.4", [d for d in main if late(d) >= 0.4]),
              ("slip undone", [d for d in main if d["status"] == "slipped"]),
              ("trimmed start", [d for d in D if d["trimmed"] and d["status"] != "unresolved" and not d["anomalous"]])]
    out = [f"# Program N: audit metrics on the plugin's renders (tools/plugin_pulse_n.py), {len(D)} captures "
           f"({len(main)} untrimmed with the slip resolved, {sum(d['trimmed'] for d in D)} trimmed, "
           f"{sum(d['status'] == 'unresolved' for d in D)} unresolved slip, out), R channel unless stated",
           "Cells: plugin |err| / bench |err| [twin tolerance |err|] (n), medians; 'in tol': captures whose |err| <= the largest of "
           "their 3 twins' |err| (plugin / bench); 'p>b': captures where the plugin's |err| exceeds the bench's by more than 2 x "
           "the twin tolerance. Decay fits only where the capture's EDC is straight within 1 dB.",
           "Anomalous (the bench itself nulls above -20 dB at 0.05-0.15 s; out): " + ", ".join(
               kstr(d["key"]) for d in D if d["anomalous"]) + "; unresolved slip (out): " + ", ".join(
               kstr(d["key"]) for d in D if d["status"] == "unresolved"),
           f"Fit checks: plugin window nulls recomputed vs plugin_check files, max |diff| {max(d['plugin']['check_db'] for d in D):.4f} dB; "
           f"bench (n = 0 captures) vs bench files, max |diff| {max((d['bench']['check_db'] or 0) for d in D):.4f} dB"]
    summ = {}

    def cell(G, k, sc=1.0, fmt="{:.2f}", f=abs):
        p, n1 = _md([d["S_p"] for d in G], k, f)
        b_, _ = _md([d["S_b"] for d in G], k, f)
        tw = [max((abs(t[k]) for t in d["S_t"] if t.get(k) is not None and np.isfinite(t[k])), default=None) for d in G]
        tv = [x for x in tw if x is not None]
        t_ = float(np.median(tv)) if tv else np.nan
        worse = sum(1 for d, x in zip(G, tw) if x is not None and d["S_p"].get(k) is not None and d["S_b"].get(k) is not None
                    and abs(d["S_p"][k]) > abs(d["S_b"][k]) + 2 * x)
        ip = sum(1 for d, x in zip(G, tw) if x is not None and d["S_p"].get(k) is not None and abs(d["S_p"][k]) <= x)
        ib = sum(1 for d, x in zip(G, tw) if x is not None and d["S_b"].get(k) is not None and abs(d["S_b"][k]) <= x)
        nn = sum(1 for d, x in zip(G, tw) if x is not None and d["S_p"].get(k) is not None)
        summ.setdefault(gname, {})[k] = dict(plugin=p * sc, bench=b_ * sc, twin=t_ * sc, n=n1, in_tol=[ip, ib, nn], plugin_worse=worse)
        return f"{fmt.format(p * sc)} / {fmt.format(b_ * sc)} [{fmt.format(t_ * sc)}] ({n1}; in tol {ip}/{ib}; p>b {worse})"

    for gname, G in groups:
        if not G:
            continue
        out.append(f"\n## {gname}: {len(G)} captures")
        out.append("### G2 octave decay (-5/-25 dB) |ratio - 1| %")
        out += [f"- {f}: " + cell(G, f"oct_t{f}", 100, "{:.2f}") for f in OCT]
        out.append("### G2 EDT (0/-10 dB) |ratio - 1| %")
        out += [f"- {f}: " + cell(G, f"oct_edt{f}", 100, "{:.2f}") for f in OCT]
        out.append("### G2 late decay (-25/-40 dB) |ratio - 1| %")
        out += [f"- {f}: " + cell(G, f"oct_late{f}", 100, "{:.2f}") for f in OCT]
        out.append("### N tail decay (late section, EDC from 40 ms after the last main tap, -5/-25) |ratio - 1| %")
        out += [f"- {f}: " + cell(G, f"tail_{f}", 100, "{:.2f}") for f in ["bb"] + [str(x) for x in OCT]]
        out.append("### G2 EDC shape: largest EDC difference to -30 dB (dB), R / L")
        out += [f"- {f}: " + cell(G, f"oct_dev{f}") + " | L " + cell(G, f"oct_dev{f}_L") for f in OCT]
        out.append("### G3 third-octave decay (-5/-25) |ratio - 1| %")
        out.append("- " + " | ".join(f"{f}: " + cell(G, f"third_t{f}", 100, "{:.2f}") for f in THIRD))
        out.append("### G3 third-octave EDC shape (dB)")
        out.append("- " + " | ".join(f"{f}: " + cell(G, f"third_dev{f}") for f in THIRD))
        out.append("### G4 echo density mean |diff| (50-300 ms), mixing time |ratio - 1| %")
        out.append("- NED " + cell(G, "ned", 1, "{:.4f}") + "; mixing time " + cell(G, "rv_mix_ms", 100, "{:.2f}"))
        out.append("### G4 L/R coherence |diff| per third octave 50-300 ms")
        out.append("- " + " | ".join(f"{f}: " + cell(G, f"coh_{f}", 1, "{:.4f}") for f in THIRD[::2]))
        out.append("### G5 fine-structure correlation (value, not |err|: 1 = identical), R")
        out.append("- " + " | ".join(f"{lo}-{hi}: " + cell(G, f"fine_{lo}-{hi}", 1, "{:.4f}", f=lambda x: x) for lo, hi in pc.FS_BANDS))
        out.append("### G6 C50 / C80 |diff| dB, D50 |diff|, centre time Ts |diff| ms per octave")
        for x, sc, fmt in (("c50", 1, "{:.3f}"), ("c80", 1, "{:.3f}"), ("d50", 1, "{:.4f}"), ("ts", 1, "{:.3f}")):
            out.append(f"- {x}: " + " | ".join(f"{f}: " + cell(G, f"rv_{x}_{f}", sc, fmt) for f in OCT))
        out.append("### G6 tonal balance, third-octave level |diff| dB, early 0-80 ms / late 80 ms on")
        out.append("- early: " + " | ".join(f"{f}: " + cell(G, f"rv_lev_e_{f}", 1, "{:.3f}") for f in THIRD[::2]))
        out.append("- late: " + " | ".join(f"{f}: " + cell(G, f"rv_lev_l_{f}", 1, "{:.3f}") for f in THIRD[::2]))
        out.append("### G6 colouration: fine-structure sd |diff| dB")
        out.append("- " + " | ".join(f"{lo}-{hi}: " + cell(G, f"rv_fsd_{lo}-{hi}", 1, "{:.3f}") for lo, hi in pc.FS_BANDS))
        out.append("### G6 narrow-band envelope (1/24 octave, first 30 dB), largest |level diff| dB: p90 / worst band")
        out.append("- p90 " + cell(G, "rv_nb_p90") + "; worst " + cell(G, "rv_nb_max"))
        out.append("### G6 IACC |diff| early / late, L/R level |diff| dB early / late")
        out.append("- IACC e " + cell(G, "rv_iacc_e", 1, "{:.4f}") + "; l " + cell(G, "rv_iacc_l", 1, "{:.4f}"))
        out.append("- L/R e " + cell(G, "rv_lr_e", 1, "{:.3f}") + "; l " + cell(G, "rv_lr_l", 1, "{:.3f}"))
        # G1 nulls
        out.append("### G1 nulls 100 Hz-10 kHz, R (L), median plugin / bench, excess over noise floor, n at or under -50 dB")
        for t0, t1 in nm.WINDOWS:
            w = f"{t0}-{t1}"
            rows = [d for d in G if w in d["plugin"]["windows"] and w in d["bench"]["windows"]]
            if not rows:
                continue
            v = lambda who, c: [d[who]["windows"][w][f"null{c}"] for d in rows]
            fl = np.median([d["plugin"]["windows"][w]["floor1"] for d in rows])
            pr_, br_ = np.median(v("plugin", 1)), np.median(v("bench", 1))
            summ.setdefault(gname, {})[f"null_{w}"] = dict(plugin=[float(np.median(v("plugin", 1))), float(np.median(v("plugin", 0)))],
                                                           bench=[float(br_), float(np.median(v("bench", 0)))], floor=float(fl), n=len(rows))
            out.append(f"- {w} s ({len(rows)}): plugin {pr_:.1f} ({np.median(v('plugin', 0)):.1f}), excess {pr_ - fl:+.1f}, "
                       f"{sum(x <= -50 for x in v('plugin', 1))} <= -50 | bench {br_:.1f} ({np.median(v('bench', 0)):.1f}), "
                       f"{sum(x <= -50 for x in v('bench', 1))} <= -50 | floor {fl:.1f}")
        cells = []
        for lo, hi in nm.BANDS:
            bk = f"{lo}-{hi}"
            rows = [d for d in G if bk in d["plugin"]["bands"] and bk in d["bench"]["bands"]]
            if rows:
                cells.append(f"{bk}: {np.median([d['plugin']['bands'][bk][1] for d in rows]):.1f} / {np.median([d['bench']['bands'][bk][1] for d in rows]):.1f}")
        out.append("- per band 0.05-0.6 s, R, plugin / bench: " + "; ".join(cells))
    # LATE pairs: real capture-vs-capture of the early section
    Pp = [d for d in main if d.get("pair")]
    if Pp:
        gname = "LATE pairs"
        out.append(f"\n## Real capture vs capture, early section (LATE pairs, SPAN >= 0.2): {len(Pp)} pairs")
        out.append("Cells: plugin |err| / bench |err| on these captures [capture pair |diff|], medians")

        def pcell(k, fmt="{:.3f}"):
            p = _md([d["S_p"] for d in Pp], f"rv_{k}")[0]
            b_ = _md([d["S_b"] for d in Pp], f"rv_{k}")[0]
            c = [abs(d["pair"]["a"][k] - d["pair"]["b"][k]) for d in Pp if k in d["pair"]["a"] and k in d["pair"]["b"]]
            cm = float(np.median(c)) if c else np.nan
            ip = sum(1 for d in Pp if d["S_p"].get(f"rv_{k}") is not None and k in d["pair"]["a"] and abs(d["S_p"][f"rv_{k}"]) <= abs(d["pair"]["a"][k] - d["pair"]["b"][k]))
            summ.setdefault(gname, {})[k] = dict(plugin=p, bench=b_, pair=cm, n=len(c), plugin_within=ip)
            return f"{fmt.format(p)} / {fmt.format(b_)} [{fmt.format(cm)}]"
        out.append("- early third-octave level |diff| dB: " + " | ".join(f"{f}: {pcell(f'lev_e_{f}')}" for f in THIRD[::2]))
        out.append(f"- early IACC {pcell('iacc_e', '{:.4f}')}; early L/R {pcell('lr_e')}")
        np_ = [d["pair"]["null_0_50"][1] for d in Pp]
        pp_ = [d["plugin"]["windows"]["0.0-0.05"]["null1"] for d in Pp if "0.0-0.05" in d["plugin"]["windows"]]
        bp_ = [d["bench"]["windows"]["0.0-0.05"]["null1"] for d in Pp if "0.0-0.05" in d["bench"]["windows"]]
        summ[gname]["null_0-0.05"] = dict(plugin=float(np.median(pp_)), bench=float(np.median(bp_)), pair=float(np.median(np_)))
        out.append(f"- 0-0.05 s null R: plugin {np.median(pp_):.1f} / bench {np.median(bp_):.1f} [pair {np.median(np_):.1f}]")
    open(LOG, "w").write("\n".join(out) + "\n")
    json.dump(summ, open(os.path.join(os.path.dirname(OUT), "summary_n.json"), "w"), indent=1,
              default=lambda x: None if isinstance(x, float) and not np.isfinite(x) else x)
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["run", "report"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    if a.stage == "run":
        run_stage(a.workers, a.only.split(",") if a.only else None, a.limit)
    else:
        report_stage()


if __name__ == "__main__":
    main()
