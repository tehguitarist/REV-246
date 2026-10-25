"""NEXT_STEPS §5 item 6a: the pulse check. Stage G of docs/analysis_harness_spec.md as one validation report: the current
models against every dataset, on the metrics the nulls can hide. Results: docs/model_s_results.md §26.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/pulse_check.py --stage s [--only keys] [--workers 8]
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/pulse_check.py --stage units
    EMT_KERNEL=chain3 .venv/bin/python tools/pulse_check.py --stage ml [--prog ML] [--only N] [--workers 6]
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/pulse_check.py --stage report
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/pulse_check.py --stage gate --tag <name>   # ~3 min sanity gate

  s       Program S, Tomkin 96 kHz library, all 128 settings held out: each setting rendered with the constants of the
          loop_cv fold that left its RT step out (the §32.6 baseline, analysis/loop_cv/hfcells_chain3_cap0.8475_pallf with EMT_HFRT_MIN=0.12: fold
          constants + cells.json), per capture clock, shift and gain as q2_injection.fit_capture (3 numbers). The capture is
          fade-compensated (/ ramp, up to comp_end). One JSON per setting in analysis/pulse/s/ (resumable).
  units   Program S on other units, nothing refitted: the CLA IR (LF 1.0, RT 1.0, HF 0.9, AMPL 80) and the five Altiverb
          S presets (RT 0.6-2.2, same panel otherwise; altiverb_predict.py). Model constants: every RT step (fold_all,
          the ramp as the median of the cap folds'); AP230 step 0. Per IR clock, shift and gain, then the 7-tap channel FIR
          per channel on the cluster-only region (altiverb_predict.channel_fir, 14 numbers), clock refitted. Metrics up
          to where the IR falls 60 dB below its 50-300 ms level (max 1.5 s). Unit against unit: the Altiverb RT 1.0
          against CLA (the same panel setting), both aligned to the model -> analysis/pulse/units.json.
          --chain: the model first through CLA's capture chain relative to Tomkin's (chain_rel: the per-chain high-pass
          and S&H pole of NEXT_STEPS §4, nothing fitted here) -> units_chain.json.
  ml      Programs M and L, every AMPL group of the library: the static model (ml_lines law with S's HF pole curve,
          theta = 0; ml_lines.fitted: per capture clock, shift and two gains per channel fitted on the deterministic part
          before any line, everything after predicted) against each capture, and every capture pair. EDCs from the R
          onset + 20 ms, coherence / fine structure / echo density over the 250 ms after it, up to the capture's fade.
          Nulls (one gain per window and channel, as ml_variance.pair_nulls) per recirculation and time window give the
          capture-to-capture decorrelation curve beside the static model's error. One JSON per group in analysis/pulse/ml/.
  gate    The sanity gate after any model change (~3 min): 16 S settings (two per RT step, every LF and HF step), CLA
          and the Altiverb RT 1.6 preset, 8 M / L groups (ml_variance's SIM_KEYS), each scored as above, then PASS / FAIL
          against the full run on the same keys (nulls not worse by more than 0.5 dB; the largest decay error not more
          than 0.5 % worse, and its 90th percentile within +-3 %; EDC: the largest not 0.1 dB worse, the 90th percentile
          within 0.5 dB; the S decay and EDC tails leave out the bad AMPL-0 and slipped captures, DEFECTS, model_s_results §26, §28-29) -> analysis/pulse/gate_<tag>/gate.log. The model is whatever
          setting_G builds: PULSE_CV=<loop_cv fold folder> and EMT_PMAX select it; a new model term must be wired into
          setting_G (S) or ml_lines (M / L) before the gate can see it.
  report  Medians per group, model against capture beside capture against capture -> analysis/pulse/report.log.

Metrics, each for model vs capture and for capture vs capture (the tolerance; S's four AMPL captures are two takes,
0 + 90 and 30 + 60, so the cross-take pairs 0-30 and 60-90 are the repeat spread and the in-take pairs 0-90, 30-60 a
check), L and R:
  G1 nulls   100 Hz-10 kHz per window (q1.nulls) and per octave (q2.octave_nulls, R).
  G2 EDC     per octave 125 Hz-8 kHz from 50 ms after tap 0: decay time from the EDC between -5 and -25 dB (fit points
             at least 150 ms before the window's end; ratio a / b - 1), the late decay time (-25 to -40 dB, where the
             window reaches it), EDT (0 to -10 dB, ratio), and the largest EDC difference down to -30 dB.
  G3 decays  the same EDC decay time per third octave 100 Hz-8 kHz (the per-band decay the null weights by energy).
  G4         normalised echo density (Abel-Huang, 20 ms Hann, 100 Hz-10 kHz) over 0-300 ms: mean |difference| 20-300 ms;
             L/R coherence |sum L R*|^2 / (sum |L|^2 sum |R|^2) per third octave over 50-300 ms: each IR's value per band.
  G5         spectral fine structure: log magnitude 50-300 ms (Hann), less its 50 Hz moving average, correlated between
             the two IRs in 100 Hz-1 kHz, 1-4 kHz and 4-10 kHz.
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
if not os.environ.get("PULSE_CV"):                 # the baseline's RT x HF floor (model_s_results 29), read by q2 at import
    os.environ.setdefault("EMT_HFRT_MIN", "0.12")
import numpy as np
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import loop_fir as lf_  # noqa: E402

q2 = lf_.q2
q1 = q2.q1
FS, PRE = q1.FS, q1.PRE
CV = os.environ.get("PULSE_CV") or os.path.join(ROOT, "analysis", "loop_cv", "hfcells_chain3_cap0.8475_pallf")   # model_s_results 32.6
REF = os.path.join(ROOT, "analysis", "pulse")                   # the full run: the gate's reference
OUT = os.environ.get("PULSE_OUT") or REF                        # PULSE_OUT: a candidate's per-group files and report beside the full run's
COMPACT = f"diffuser_tied_fix-0.707107_ap230new_{q1.ms.KTAG}"
BB = butter(4, (100, 10000), btype="band", fs=FS, output="sos")
OCT = q2.OCT                                                    # 125 Hz-8 kHz
THIRD = [round(100 * 2 ** (k / 3)) for k in range(19)]           # 100 Hz-6.3 kHz ... 8 kHz
THIRD = [f for f in THIRD if f <= 8100]
SOS = {("o", f): q2.OCT_SOS[f] for f in OCT}
SOS.update({("t", f): butter(3, [f / 2 ** (1 / 6), min(f * 2 ** (1 / 6), 11000)], btype="band", fs=FS, output="sos") for f in THIRD})
T0 = 0.05                                                       # EDC start, s after tap 0 (past the cluster)
FS_WIN = (0.05, 0.30)                                           # coherence / fine-structure window
FS_BANDS = [(100, 1000), (1000, 4000), (4000, 10000)]
CROSS, INTAKE = [(0, 30), (60, 90)], [(0, 90), (30, 60)]


# ---------- per-IR preparation: band signals, compensated after filtering ----------
def prep(x, g=None, dtype=np.float64):
    """x (n, 2) aligned to tap 0; g the fade ramp (capture) or None (model). Band-filtered copies / g."""
    inv = 1.0 / np.maximum(g, 0.1) if g is not None else 1.0
    d = {"bb": (sosfiltfilt(BB, x, axis=0) * inv).astype(dtype)}
    for k, sos in SOS.items():
        d[k] = (sosfiltfilt(sos, x, axis=0) * inv).astype(dtype)
    return d


def edc_curve(s, i0, i1):
    e = s[i0:i1].astype(np.float64) ** 2
    c = np.cumsum(e[::-1])[::-1]
    return 10 * np.log10(c / c[0] + 1e-30)


def edc_time(d, hi, lo, n_ok):
    """Decay time (s per 60 dB) from a line through the EDC between hi and lo dB, points before n_ok only."""
    m = (d <= hi) & (d >= lo)
    m[n_ok:] = False
    if m.sum() < int(0.03 * FS):
        return None
    t = np.nonzero(m)[0] / FS
    return float(-60 / np.polyfit(t, d[m], 1)[0])


def band_decay(sa, sb, i0, i1):
    """EDC metrics of one band and channel, IR a against IR b (b is the reference)."""
    if i1 - i0 < int(0.1 * FS):
        return dict(t=None, edt=None, late=None, edc_dev=None)
    da, db = edc_curve(sa, i0, i1), edc_curve(sb, i0, i1)
    n_ok = max(0, i1 - i0 - int(0.15 * FS))
    ta, tb = edc_time(da, -5, -25, n_ok), edc_time(db, -5, -25, n_ok)
    ea, eb = edc_time(da, 0, -10, n_ok), edc_time(db, 0, -10, n_ok)
    la, lb = edc_time(da, -25, -40, n_ok), edc_time(db, -25, -40, n_ok)
    m = db[:n_ok] >= -30
    return dict(t=None if ta is None or tb is None else [ta, tb], edt=None if ea is None or eb is None else [ea, eb],
                late=None if la is None or lb is None else [la, lb],
                edc_dev=float(np.max(np.abs(da[:n_ok][m] - db[:n_ok][m]))) if m.sum() > int(0.03 * FS) else None)


def ned(s, i1, win=0.02, k0=PRE, t0=0.0):
    """Normalised echo density (Abel-Huang) from t0 after tap 0 (index k0), 20 ms Hann, 5 ms hop, 300 ms (or to i1)."""
    w = int(win * FS); hw = np.hanning(w); hw /= hw.sum()
    out = []
    a = k0 + int(t0 * FS)
    for i in range(a, min(i1, a + int(0.3 * FS)) - w, w // 4):
        seg = s[i:i + w]
        sd = np.sqrt(np.sum(hw * seg ** 2))
        out.append(float(np.sum(hw * (np.abs(seg) > sd)) / 0.3173))
    return np.array(out)


def spectra(x, i1, k0=PRE, fwin=FS_WIN):
    i0, j1 = k0 + int(fwin[0] * FS), min(k0 + int(fwin[1] * FS), i1)
    seg = x[i0:j1] * np.hanning(j1 - i0)[:, None]
    N = 1 << 16
    return np.fft.rfft(seg, N, axis=0), np.fft.rfftfreq(N, 1 / FS)


def coherence(X, f):
    out = {}
    for fc in THIRD:
        m = (f >= fc / 2 ** (1 / 6)) & (f < fc * 2 ** (1 / 6))
        out[fc] = float(np.abs(np.sum(X[m, 0] * np.conj(X[m, 1]))) ** 2 / (np.sum(np.abs(X[m, 0]) ** 2) * np.sum(np.abs(X[m, 1]) ** 2)))
    return out


def fine(X, f):
    L = 20 * np.log10(np.abs(X) + 1e-12)
    k = int(50 / (f[1] - f[0]))
    ker = np.ones(k) / k
    return L - np.stack([np.convolve(L[:, c], ker, "same") for c in (0, 1)], 1)


# ---------- G6: reverb-specific metrics of one IR (model or capture), compared afterwards ----------
EDR_BANDS = 100 * 2 ** (np.arange(0, 24 * np.log2(100) + 1) / 24)   # 1/24-octave centres, 100 Hz - 10 kHz
EDR_BANDS = EDR_BANDS[EDR_BANDS <= 10000]
EARLY = 0.08                                                    # s after tap 0: early / late split (C80's)


def edr_times(x, i0, i1, nfft=4096, hop=512):
    """Decay time (s per 60 dB, EDC -5 to -25 dB) per 1/24-octave band from a short-time spectrum of x[i0:i1] (an energy
    decay relief at fine frequency resolution); None where the band's range is too short. nfft 4096 (6.4 Hz): fine enough
    to separate neighbouring modes of different loops; with 2048 (12.9 Hz) two modes ~7 Hz apart share a bin and their
    beating makes the -5/-25 dB fit flip by 20-30 % for a 0.1 dB envelope difference (model_s_results 33.1)."""
    seg = x[i0:i1]
    if len(seg) < nfft + 8 * hop:
        return None
    win = np.hanning(nfft)
    fr = np.lib.stride_tricks.sliding_window_view(seg, nfft)[::hop] * win
    P = np.abs(np.fft.rfft(fr, axis=1)) ** 2
    f = np.fft.rfftfreq(nfft, 1 / FS)
    out = []
    for fc in EDR_BANDS:
        m = (f >= fc * 2 ** (-1 / 48)) & (f < fc * 2 ** (1 / 48))
        if not m.any():
            m = np.abs(f - fc) == np.min(np.abs(f - fc))
        e = P[:, m].sum(1)
        c = np.cumsum(e[::-1])[::-1]
        d = 10 * np.log10(c / c[0] + 1e-30)
        k = (d <= -5) & (d >= -25)
        k[max(0, len(d) - int(0.15 * FS / hop)):] = False
        out.append(float(-60 / np.polyfit(np.nonzero(k)[0] * hop / FS, d[k], 1)[0]) if k.sum() >= 4 else None)
    return out


def nb_env_dev(a, b, i0, i1, top=30.0):
    """Narrow-band ringing check (model_s_results 33.1): per 1/24-octave band (EDR_BANDS, 4th-order band-pass, at least
    +-1.6 Hz), the largest |level difference| (dB) between the 10 ms energy envelopes of a and b while b is within `top` dB
    of its peak after i0. Unlike a fitted decay time it does not flip when two loops' modes beat inside one band."""
    w = int(0.01 * FS)
    if i1 - i0 < 4 * w:
        return None
    out = []
    for fc in EDR_BANDS:
        bw = fc * (2 ** (1 / 48) - 2 ** (-1 / 48)) / 2 + 1.6
        sos = butter(4, [fc - bw, min(fc + bw, 0.45 * FS)], btype="band", fs=FS, output="sos")
        ea, eb = (np.convolve(sosfiltfilt(sos, x[:i1]) ** 2, np.ones(w) / w, "same")[i0:i1] for x in (a, b))
        k = eb > eb.max() * 10 ** (-top / 10)
        out.append(float(np.max(np.abs(10 * np.log10((ea[k] + 1e-30) / (eb[k] + 1e-30))))) if k.sum() > w else None)
    return out


def ned_curve(s, i1, k0, span=1.0, win=0.02):
    """Normalised echo density from tap 0 over up to `span` s (20 ms Hann, 5 ms hop)."""
    w = int(win * FS); hw = np.hanning(w); hw /= hw.sum()
    out = []
    for i in range(k0, min(i1, k0 + int(span * FS)) - w, w // 4):
        seg = s[i:i + w]
        sd = np.sqrt(np.sum(hw * seg ** 2))
        out.append(float(np.sum(hw * (np.abs(seg) > sd)) / 0.3173))
    return np.array(out)


def iacc(xl, xr, lag=int(0.001 * FS)):
    """Largest normalised L/R cross-correlation within +-1 ms."""
    n = np.sqrt(np.sum(xl ** 2) * np.sum(xr ** 2)) + 1e-30
    c = np.correlate(xl, xr, "full")
    mid = len(xr) - 1
    return float(np.max(np.abs(c[mid - lag:mid + lag + 1])) / n)


def rv_metrics(X, i1, k0):
    """G6 per IR (prepared dict, R unless stated): C50 / C80 (dB), D50, centre time Ts (ms) per octave from tap 0;
    third-octave level (dB) early (0-80 ms) and late (80 ms - i1); EDR decay times per 1/24 octave from 50 ms;
    echo density over 1 s and the mixing time (first time the 50 ms-smoothed NED reaches 0.9, ms); L/R: IACC and level
    difference (dB) early and late, broadband and per octave 500 Hz - 4 kHz."""
    r = {}
    ne, nl = int(0.05 * FS), int(EARLY * FS)
    for f in OCT:
        e = X[("o", f)][k0:i1, 1].astype(np.float64) ** 2
        if len(e) <= nl:
            continue
        t = np.arange(len(e)) / FS
        E = e.sum() + 1e-30
        r[f"c50_{f}"] = float(10 * np.log10(e[:ne].sum() / (e[ne:].sum() + 1e-30)))
        r[f"c80_{f}"] = float(10 * np.log10(e[:nl].sum() / (e[nl:].sum() + 1e-30)))
        r[f"d50_{f}"] = float(e[:ne].sum() / E)
        r[f"ts_{f}"] = float(1000 * np.sum(t * e) / E)
    for f in THIRD:
        e = X[("t", f)][k0:i1, 1].astype(np.float64) ** 2
        if len(e) > nl:
            r[f"lev_e_{f}"] = float(10 * np.log10(e[:nl].sum() + 1e-30))
            r[f"lev_l_{f}"] = float(10 * np.log10(e[nl:].sum() + 1e-30))
    bb = X["bb"].astype(np.float64)
    r["edr"] = edr_times(bb[:, 1], k0 + int(0.05 * FS), i1)
    nc = ned_curve(bb[:, 1], i1, k0)
    r["ned1"] = nc[::4].round(3).tolist()
    sm = np.convolve(nc, np.ones(10) / 10, "valid") if len(nc) >= 10 else nc
    hit = np.nonzero(sm >= 0.9)[0]
    r["mix_ms"] = float((hit[0] + 5) * 5.0) if len(hit) else None
    for tag, a, b in (("e", k0, k0 + nl), ("l", k0 + nl, min(i1, k0 + int(0.5 * FS)))):
        if b - a > int(0.02 * FS):
            r[f"iacc_{tag}"] = iacc(bb[a:b, 0], bb[a:b, 1])
            r[f"lr_{tag}"] = float(10 * np.log10((np.sum(bb[a:b, 0] ** 2) + 1e-30) / (np.sum(bb[a:b, 1] ** 2) + 1e-30)))
            for f in (500, 1000, 2000, 4000):
                o = X[("o", f)][a:b].astype(np.float64)
                r[f"iacc_{tag}_{f}"] = iacc(o[:, 0], o[:, 1])
    return r


def compare(A, B, i1, k0=PRE, t0=T0, fwin=FS_WIN):
    """All G2-G5 metrics of IR A against IR B (prepared dicts) up to device index i1; tap 0 at k0, EDCs from t0,
    coherence and fine structure over fwin (s after tap 0)."""
    i0 = k0 + int(t0 * FS)
    r = dict(oct={}, third={})
    for (kind, f), key in [(k, k) for k in SOS]:
        r["oct" if kind == "o" else "third"][f] = [band_decay(A[key][:, c], B[key][:, c], i0, i1) for c in (0, 1)]
    na, nb = ned(A["bb"][:, 1], i1, k0=k0, t0=0.0 if k0 == PRE else fwin[0]), ned(B["bb"][:, 1], i1, k0=k0, t0=0.0 if k0 == PRE else fwin[0])
    r["ned"] = dict(a=na[::4].round(3).tolist(), b=nb[::4].round(3).tolist(), dev=float(np.mean(np.abs(na[3:] - nb[3:]))) if len(na) > 3 else None)
    if i1 - k0 - int(fwin[0] * FS) > int(0.15 * FS):
        XA, f = spectra(A["bb"], i1, k0, fwin)
        XB, _ = spectra(B["bb"], i1, k0, fwin)
        r["coh"] = dict(a=coherence(XA, f), b=coherence(XB, f))
        FA, FB = fine(XA, f), fine(XB, f)
        r["fine"] = {f"{lo}-{hi}": [float(np.corrcoef(FA[(f >= lo) & (f < hi), c], FB[(f >= lo) & (f < hi), c])[0, 1]) for c in (0, 1)]
                     for lo, hi in FS_BANDS}
        r["fine_sd"] = {f"{lo}-{hi}": [float(np.std(FA[(f >= lo) & (f < hi), 1])), float(np.std(FB[(f >= lo) & (f < hi), 1]))]
                        for lo, hi in FS_BANDS}                 # G6: tail coloration (Gaussian tail ~5.6 dB)
    r["rv"] = dict(a=rv_metrics(A, i1, k0), b=rv_metrics(B, i1, k0))   # G6
    r["nbenv"] = nb_env_dev(A["bb"][:, 1].astype(np.float64), B["bb"][:, 1].astype(np.float64), k0 + int(0.05 * FS), i1)
    return r


# ---------- program S on the library ----------
def setting_G(lab, b230, law, fold, fc):
    """The baseline's group spectra for one setting with its fold's constants (as loop_cv run -> q2 _final_setting)."""
    law = json.loads(json.dumps(law))
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    for lf_k, dc in fold.get("lfdc", {}).items():
        law["lf"][float(lf_k)]["dc"] = float(dc)
    gamma, fl, hfdp = fold["gamma"], {}, fold["hfdp"]
    if fc and "dc_min" in fc:
        law["lf"][0.1]["dc"] = fc["dc"]
        gamma, fl = fc["gamma"], dict(dc_min=fc["dc_min"], beta=fc["beta"], gamma_dc=fc.get("gamma_dc", True))
    elif q2.dcfloor_form(lab):                      # EMT_DCFLOOR: the floor's compact form (model_s_results 34)
        fl = dict(floor_form=q2.dcfloor_form(lab))
    if fc and "hf_dp" in fc:
        hfdp = {str(lab["hf"]): fc["hf_dp"], "p_above": 0.0, "p_below": 1.0, "exclude_lf": [], "slope": fc.get("hf_slope", 0.0)}
    if fc and "hf_poles" in fc:                     # each loop's HF pole outright (model_s_results 32.4)
        hfdp = {"poles": fc["hf_poles"]}
    if fc and "lf_gain" in fc:                      # the (LF, RT) cell's per-loop gain and DC section (model_s_results 32.2)
        fl["loopcell"] = dict(gain=fc["lf_gain"], dc=fc["lf_dc"])
    return q2.final_G(lab, b230, law, 0.0, 0.0, fold["lfrt"], gamma, None, hfdp, tuple(fold["pramp"]), **fl)


def _s_job(a):
    key, lab, b230, law, fold, fc, v, caps, ends, path = a
    G = setting_G(lab, b230, law, fold, fc)
    P, rows = {}, []
    for c in caps:
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        i1 = PRE + int(te * FS)
        Y, M = prep(y, g), prep(mdl)
        P[c["amp"]] = (Y, i1, y, g)
        rows.append(dict(amp=c["amp"], file=c["file"], t_comp=te, null=q1.nulls(y, mdl * g, te, q2.W_COMP, narrow=True),
                         null_oct=q2.octave_nulls(y, mdl * g, te), m=compare(M, Y, i1)))
    pairs = []
    for a_, b_ in CROSS + INTAKE:
        (Ya, ia, ya, ga), (Yb, ib, yb, gb) = P[a_], P[b_]
        i1 = min(ia, ib)
        te = (i1 - PRE) / FS
        n = min(len(ya), len(yb))
        pairs.append(dict(pair=[a_, b_], cross=(a_, b_) in CROSS, t_comp=te,
                          null=q1.nulls(ya[:n] * gb[:n], yb[:n] * ga[:n], te, q2.W_COMP[1:], narrow=True), m=compare(Ya, Yb, i1)))
    out = dict(key=key, **lab, b230=b230, fold=fold["hold"], rows=rows, pairs=pairs)
    json.dump(out, open(path, "w"))
    nl = rows[1]["null"]
    print(f"{key}: R null " + " ".join(f"{w}:{x[1]:.1f}" for w, x in nl.items() if x), flush=True)
    return key


def s_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, "s")
    os.makedirs(d, exist_ok=True)
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=COMPACT)
    cells = json.load(open(os.path.join(CV, "cells.json")))
    folds = {rt: json.load(open(os.path.join(CV, f"fold_{rt}.json"))) for rt in sorted({g["rt"] for g in info.values()})}
    jobs = []
    for k, g in info.items():
        path = os.path.join(d, k + ".json")
        if (only and k not in only) or os.path.exists(path):
            continue
        lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
        caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(g)]
        jobs.append((k, lab, ap["steps"][ap["settings"][k]["step"]], law, folds[g["rt"]], cells.get(k), v, caps, ends, path))
    print(f"{len(jobs)} settings to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_s_job, jobs):
            pass


# ---------- program S on other units (CLA, Altiverb) ----------
def chain_rel(x, hp=(6.65, 6.65), hp_ref=(1.8,), p=0.325, p_ref=0.301):
    """x (n, 2) through the CLA capture chain relative to Tomkin's (NEXT_STEPS §4, analog_results §2-3, nothing fitted
    here): the per-chain high-pass terms (CLA 2 x 6.65 Hz, Tomkin 1.8 Hz; the unit's own corners cancel) and the output
    S&H pole (CLA 0.325, Tomkin 0.301)."""
    sys.path.insert(0, os.path.join(HERE, "analog"))
    from sample_hold import highpass, output_sh
    n = len(x)
    N = 1 << int(np.ceil(np.log2(n + 8192)))
    f = np.fft.rfftfreq(N, 1 / FS)
    H = highpass(f, hp) / np.where(f > 0, highpass(f, hp_ref), 1.0) * output_sh(f, p) / output_sh(f, p_ref)
    H[0] = 0.0
    return np.fft.irfft(np.fft.rfft(x, N, axis=0) * H[:, None], N, axis=0)[:n]


def units_stage(chain=False, only=None):
    import altiverb_predict as av
    q1.EXT = 1700
    v, law, _ = q2.final_injection(compact=COMPACT)
    ap = json.load(open(os.path.join(q2.OUT, lf_.AP230)))
    fold = json.load(open(os.path.join(ROOT, "analysis", "loop_cv", "fold_all.json")))
    fold["pramp"] = [float(np.median([json.load(open(p))["pramp"][0] for p in
                                      glob.glob(os.path.join(ROOT, "analysis", "loop_cv", f"cap{os.environ['EMT_PMAX']}_folds", "fold_*.json"))])), 0.62]
    irs = [("CLA", 1.0, av.cla_h())] + [(f"Altiverb RT {k}", rt, av.altiverb_h(k)[0]) for k, rt in av.PRESETS.items() if not only or k in only]
    res, keep = {}, {}
    for name, rt, h in irs:
        lab = dict(lf=1.0, rt=rt, hf=0.9)
        G = setting_G(lab, ap["steps"][0], law, fold, None)
        yb = sosfiltfilt(q1.LFBAND, h * (np.sign(h[PRE].sum()) or 1.0), axis=0)
        n = min(len(yb) - 600, G["N"])
        mdl = q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=80)])]).render(v, G, n, 80)
        if chain:
            mdl = chain_rel(mdl)
        eg, shg, gain = av.fit_clock(yb, mdl, starts=[98.4] if name == "CLA" else [42.8])
        y = q1.apply_clock(yb, eg, shg)[:n] / gain
        plain = q1.nulls(y, mdl, 1.5, q2.W_COMP[:6], narrow=True)
        mdl2, hs = av.channel_fir(y, mdl)
        eg2, shg2, gain2 = av.fit_clock(yb, mdl2, starts=[eg])
        y2 = q1.apply_clock(yb, eg2, shg2)[:n] / gain2
        te = min(1.5, q2.comp_end(y2, np.ones_like(y2)))
        i1 = PRE + int(te * FS)
        Y = prep(y2)
        keep[name] = (Y, i1)
        res[name] = dict(rt=rt, eps_ppm=eg2, gain_db=float(20 * np.log10(gain2)), t_end=te, null_plain=plain, fir=[np.round(x, 4).tolist() for x in hs],
                         null=q1.nulls(y2, mdl2, te, q2.W_COMP[:6], narrow=True), null_oct=q2.octave_nulls(y2, mdl2, te), m=compare(prep(mdl2), Y, i1))
        print(f"{name}: clock {eg2:+.1f} ppm, t_end {te:.2f} s | null R " + " ".join(f"{w}:{x[1]:.1f}" for w, x in res[name]["null"].items() if x), flush=True)
    if "Altiverb RT 1.0" in keep:
        (YA, ia), (YC, ic) = keep["Altiverb RT 1.0"], keep["CLA"]
        res["unit_vs_unit"] = dict(t_end=(min(ia, ic) - PRE) / FS, m=compare(YA, YC, min(ia, ic)))
    json.dump(res, open(os.path.join(OUT, "units_chain.json" if chain else "units.json"), "w"))


def units_report(chain=False):
    path = os.path.join(OUT, "units_chain.json" if chain else "units.json")
    if not os.path.exists(path):
        return []
    U = json.load(open(path))
    if "unit_vs_unit" not in U:
        return []
    ref = json.load(open(os.path.join(OUT, "s", "1.0_1.0_1.0.json")))
    rp = [summarise(p["m"]) for p in ref["pairs"] if p["cross"]]
    uu = summarise(U["unit_vs_unit"]["m"])
    out = ["\n## Other units (S model, nothing refitted" + ("; through CLA's capture chain relative to Tomkin's, chain_rel" if chain else "") + "; 14-number channel FIR on the cluster-only region): model vs unit; "
           "[Altiverb RT 1.0 vs CLA, unit vs unit]; {Tomkin 1.0_1.0_1.0 cross-take repeat, median}"]
    pr = lambda k, sc=100, fmt="{:+.1f}": ((fmt.format(sc * uu[k]) if uu.get(k) is not None else "--"),
                                         (fmt.format(sc * np.median([r[k] for r in rp if r.get(k) is not None])) if any(r.get(k) is not None for r in rp) else "--"))
    for name, r in U.items():
        if name == "unit_vs_unit":
            continue
        s = summarise(r["m"])
        out.append(f"- **{name}** (t_end {r['t_end']:.2f} s, clock {r['eps_ppm']:+.1f} ppm): null R " + " ".join(f"{w}: {x[1]:.1f}" for w, x in r["null"].items() if x)
                   + " | plain " + " ".join(f"{w}: {x[1]:.1f}" for w, x in r["null_plain"].items() if x))
        out.append("  - octave decay (-5/-25) ratio-1 %, R: " + " ".join(f"{f}: {100 * s[f'oct_t{f}']:+.1f}" if s.get(f"oct_t{f}") is not None else f"{f}: --" for f in OCT))
        out.append("  - octave late decay (-25/-40) %, R: " + " ".join(f"{f}: {100 * s[f'oct_late{f}']:+.1f}" if s.get(f"oct_late{f}") is not None else f"{f}: --" for f in OCT))
        out.append("  - octave EDT %, R: " + " ".join(f"{f}: {100 * s[f'oct_edt{f}']:+.1f}" if s.get(f"oct_edt{f}") is not None else f"{f}: --" for f in OCT))
        out.append("  - octave decay %, L: " + " ".join(f"{f}: {100 * s[f'oct_t{f}_L']:+.1f}" if s.get(f"oct_t{f}_L") is not None else f"{f}: --" for f in OCT))
        out.append(f"  - NED {s['ned']:.3f}; fine R " + " ".join(f"{lo}-{hi}: {s.get(f'fine_{lo}-{hi}', np.nan):.3f}" for lo, hi in FS_BANDS)
                   + "; coherence model-unit " + " ".join(f"{f}: {s[f'coh_{f}']:+.2f}" for f in THIRD[::3] if f"coh_{f}" in s))
        out.append("  - octave nulls R (50-150/150-300/300-600): " + "; ".join(f"{f}: " + "/".join("--" if x is None else f"{x:.1f}" for x in r["null_oct"][str(f)]) for f in OCT))
        out += g6_lines([s], rp, indent="  - G6 ")
    out.append("- **unit vs unit (Altiverb RT 1.0 vs CLA) and Tomkin repeat:** octave decay %: " + " ".join(f"{f}: {pr(f'oct_t{f}')[0]} {{{pr(f'oct_t{f}')[1]}}}" for f in OCT))
    out.append("  - octave late decay %: " + " ".join(f"{f}: {pr(f'oct_late{f}')[0]} {{{pr(f'oct_late{f}')[1]}}}" for f in OCT))
    out.append(f"  - NED {uu['ned']:.3f}; fine R " + " ".join(f"{lo}-{hi}: {pr(f'fine_{lo}-{hi}', 1, '{:.3f}')[0]} {{{pr(f'fine_{lo}-{hi}', 1, '{:.3f}')[1]}}}" for lo, hi in FS_BANDS)
               + "; coherence diff " + " ".join(f"{f}: {pr(f'coh_{f}', 1, '{:+.2f}')[0]}" for f in THIRD[::3]))
    return out


# ---------- programs M and L on the library ----------
ML_TW = [(0.3, 0.6), (0.6, 1.0), (1.0, 2.0), (2.0, 4.0)]


def _ml_nulls(prog, y1, y2, end):
    """Null of y1 against g y2, one gain per window and channel: recirc1-3 (R) and ML_TW (mean of L and R), dB."""
    import ml_lines as mll
    K0 = mll.K0

    def nl(a0, a1, chs):
        v = []
        for ch in chs:
            x1, x2 = y1[a0:a1, ch], y2[a0:a1, ch]
            v.append(10 * np.log10(np.sum((x1 - np.dot(x1, x2) / np.dot(x2, x2) * x2) ** 2) / np.sum(x1 ** 2)))
        return float(np.mean(v))
    W = mll.windows(prog)
    out = {}
    for wn in ("first_pass", "recirc1", "recirc2", "recirc3"):
        a0, a1 = K0 + W[1][wn][0], min(K0 + W[1][wn][1], end)
        out[wn] = nl(a0, a1, (1,)) if a1 - a0 > 1000 else None
    for t0, t1 in ML_TW:
        a0, a1 = K0 + int(t0 * FS), min(K0 + int(t1 * FS), end)
        out[f"{t0}-{t1}"] = nl(a0, a1, (0, 1)) if a1 - a0 > 2000 else None
    return out


def _ml_job(a):
    prog, key, path = a
    import ml_lines as mll
    fit = mll.fitted(prog, key)
    k0 = mll.K0
    t0 = (mll.cc.ONSET[prog][1] + 20) / FS
    fw = (t0, t0 + 0.25)
    P, rows = {}, []
    for amp, (y, mdl, ue) in fit.items():
        i1 = min(ue)
        Y, M = prep(y[:i1], dtype=np.float32), prep(mdl[:i1], dtype=np.float32)
        P[amp] = (Y, i1, y)
        rows.append(dict(amp=amp, t_end=(i1 - k0) / FS, null=_ml_nulls(prog, y, mdl, i1), m=compare(M, Y, i1, k0, t0, fw)))
        del M
    pairs = []
    amps = sorted(P)
    for i in range(len(amps)):
        for j in range(i + 1, len(amps)):
            (Ya, ia, ya), (Yb, ib, yb) = P[amps[i]], P[amps[j]]
            i1 = min(ia, ib)
            pairs.append(dict(pair=[amps[i], amps[j]], cross=True, t_end=(i1 - k0) / FS, null=_ml_nulls(prog, ya, yb, i1), m=compare(Ya, Yb, i1, k0, t0, fw)))
    lf, rt, hf = (float(x) for x in key.split("_"))
    json.dump(dict(key=key, prog=prog, lf=lf, rt=rt, hf=hf, rows=rows, pairs=pairs), open(path, "w"))
    print(f"{prog} {key}: {len(rows)} captures, recirc1 model " + " ".join(f"{r['null']['recirc1']:.1f}" for r in rows if r["null"]["recirc1"] is not None), flush=True)
    return key


def ml_stage(progs, workers, only=None, keys_=None):
    from concurrent.futures import ProcessPoolExecutor
    import ml_carryover as cc
    d = os.path.join(OUT, "ml")
    os.makedirs(d, exist_ok=True)
    jobs = []
    for prog in progs:
        keys = sorted(cc.groups(prog)) if keys_ is None else keys_[prog]
        keys = keys[::max(1, len(keys) // only)][:only] if only else keys
        jobs += [(prog, k, os.path.join(d, f"{prog}_{k}.json")) for k in keys if not os.path.exists(os.path.join(d, f"{prog}_{k}.json"))]
    print(f"{len(jobs)} groups to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_ml_job, jobs):
            pass


def ml_report():
    D = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(OUT, "ml", "*.json")))]
    if not D:
        return []
    out = ["\n## Programs M and L (static model, theta = 0, HF pole curve): model vs capture [capture vs capture], medians"]

    def med(rr, k, f=lambda x: x):
        x = [f(s[k]) for s in rr if s.get(k) is not None and np.isfinite(s[k])]
        return np.median(x) if x else np.nan
    for prog in ("M", "L"):
        G = [d for d in D if d["prog"] == prog]
        if not G:
            continue
        rows = [summarise(r["m"]) for d in G for r in d["rows"]]
        reps = [summarise(p["m"]) for d in G for p in d["pairs"]]
        out.append(f"\n### {prog}: {len(G)} groups, {len(rows)} captures, {len(reps)} capture pairs")
        out.append("- octave decay (-5/-25 dB) signed ratio - 1 %, R: " + " ".join(f"{f}: {100 * med(rows, f'oct_t{f}'):+.1f} [{100 * med(reps, f'oct_t{f}'):+.1f}]" for f in OCT))
        out.append("- octave decay |ratio - 1| %, R: " + " ".join(f"{f}: {100 * med(rows, f'oct_t{f}', abs):.1f} [{100 * med(reps, f'oct_t{f}', abs):.1f}]" for f in OCT))
        out.append("- octave late decay (-25/-40) signed %, R: " + " ".join(f"{f}: {100 * med(rows, f'oct_late{f}'):+.1f} [{100 * med(reps, f'oct_late{f}'):+.1f}]" for f in OCT))
        out.append("- octave EDT signed %, R: " + " ".join(f"{f}: {100 * med(rows, f'oct_edt{f}'):+.1f} [{100 * med(reps, f'oct_edt{f}'):+.1f}]" for f in OCT))
        out.append("- largest EDC difference to -30 dB (dB), R: " + " ".join(f"{f}: {med(rows, f'oct_dev{f}'):.2f} [{med(reps, f'oct_dev{f}'):.2f}]" for f in OCT))
        out.append("- third-octave decay signed %, R: " + " ".join(f"{f}: {100 * med(rows, f'third_t{f}'):+.1f} [{100 * med(reps, f'third_t{f}'):+.1f}]" for f in THIRD))
        out.append(f"- echo density |diff|: {med(rows, 'ned'):.3f} [{med(reps, 'ned'):.3f}]")
        out.append("- G6, reverb-specific (static model vs capture [capture pair]):")
        out += g6_lines(rows, reps, indent="  - ")
        out.append("- L/R coherence model - capture [capture - capture] (capture value): " + " ".join(f"{f}: {med(rows, f'coh_{f}'):+.3f} [{med(reps, f'coh_{f}'):+.3f}] ({med(rows, f'cohv_{f}'):.2f})" for f in THIRD[::2]))
        out.append("- fine structure R (L): " + "; ".join(f"{lo}-{hi}: {med(rows, f'fine_{lo}-{hi}'):.3f} ({med(rows, f'fine_{lo}-{hi}_L'):.3f}) [{med(reps, f'fine_{lo}-{hi}'):.3f}]" for lo, hi in FS_BANDS))
        ws = ["first_pass", "recirc1", "recirc2", "recirc3"] + [f"{a}-{b}" for a, b in ML_TW]
        mn = {w: [r["null"][w] for d in G for r in d["rows"] if r["null"].get(w) is not None] for w in ws}
        pn = {w: [p["null"][w] for d in G for p in d["pairs"] if p["null"].get(w) is not None] for w in ws}
        out.append("- decorrelation: null per window, model (static) [capture pair]: " + "; ".join(
            f"{w}: {np.median(mn[w]):.1f} [{np.median(pn[w]):.1f}]" if mn[w] and pn[w] else f"{w}: --" for w in ws))
        for g, fn in (("RT short", lambda d: d["rt"] <= {"M": 1.0, "L": 3.0}[prog]), ("RT long", lambda d: d["rt"] >= {"M": 3.0, "L": 8.0}[prog]),
                      ("LF 0.1", lambda d: d["lf"] == 0.1), ("HF <= 0.2", lambda d: d["hf"] <= 0.2)):
            rr = [summarise(r["m"]) for d in G if fn(d) for r in d["rows"]]
            pp = [summarise(p["m"]) for d in G if fn(d) for p in d["pairs"]]
            if rr:
                out.append(f"  - {g} ({len(rr)} captures): octave decay signed %: " + " ".join(f"{f}: {100 * med(rr, f'oct_t{f}'):+.1f} [{100 * med(pp, f'oct_t{f}'):+.1f}]" for f in OCT))
    return out


# ---------- the gate: a short pulse check after any model change ----------
GATE_LH = [(0.1, 0.2), (2.0, 1.0), (0.4, 0.6), (1.0, 0.4), (0.1, 1.0), (2.0, 0.2), (0.4, 0.4), (1.0, 0.6)]
GATE_S = [f"{lf}_{rt}_{hf}" for i, rt in enumerate((0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5))
          for lf, hf in [GATE_LH[(2 * i) % 8], GATE_LH[(2 * i + 1) % 8]]]
GATE_ML = {"M": ["1.2_5.0_0.9", "2.0_0.8_0.9", "0.7_3.0_0.2", "0.1_1.8_0.4"], "L": ["1.2_20.0_1.0", "1.2_4.5_0.8", "2.0_2.6_0.5", "0.1_10.0_0.3"]}
GATE_UNITS = ["1.6"]                                            # + CLA
from defects import DEFECTS  # noqa: E402  bad AMPL-0 captures and slipped files (library_notes 2, model_s_results 28)


def gate_stage(workers, tag):
    """16 S settings (two per RT step, every LF and HF step), CLA and the Altiverb RT 1.6 preset, 8 M / L groups, into
    analysis/pulse/gate_<tag>/; then PASS / FAIL against the full run (REF) on the same keys and the absolute targets."""
    global OUT
    OUT = os.path.join(REF, f"gate_{tag}")
    s_stage(workers, GATE_S)
    if not os.path.exists(os.path.join(OUT, "units.json")):
        units_stage(only=GATE_UNITS)
    ml_stage(["M", "L"], workers, keys_=GATE_ML)
    lines, ok = [f"## Pulse gate '{tag}' (model: {os.path.basename(CV)}, EMT_PMAX {os.environ.get('EMT_PMAX')}) against {os.path.relpath(REF, ROOT)}"], True

    def check(name, val, ref, tol, absolute=None, val_abs=None):
        """val against ref + tol (regression); val_abs (default val) against the absolute target."""
        nonlocal ok
        va = val if val_abs is None else val_abs
        good = (val <= ref + tol) and (absolute is None or va <= absolute)
        ok &= good
        lines.append(f"- {'PASS' if good else 'FAIL'} {name}: {val:.2f} (reference {ref:.2f}, allowed +{tol:g}" +
                     (f"; target <= {absolute:g} on {'the 90th percentile ' + format(va, '.2f') if val_abs is not None else 'this'}" if absolute is not None else "") + ")")
    S = lambda root, k: json.load(open(os.path.join(root, "s", k + ".json")))
    for w in ("50-150", "150-300", "300-600"):
        v = [r["null"][w][1] for k in GATE_S for r in S(OUT, k)["rows"] if r["null"].get(w)]
        rf = [r["null"][w][1] for k in GATE_S for r in S(REF, k)["rows"] if r["null"].get(w)]
        check(f"S median R null {w} ms (dB)", np.median(v), np.median(rf), 0.5)
    clean = lambda root, k: [r for r in S(root, k)["rows"] if (k, r["amp"]) not in DEFECTS]   # tails without data defects (29)
    for f in OCT:
        v = [abs(summarise(r["m"])[f"oct_t{f}"] or 0) * 100 for k in GATE_S for r in clean(OUT, k)]
        rf = [abs(summarise(r["m"])[f"oct_t{f}"] or 0) * 100 for k in GATE_S for r in clean(REF, k)]
        check(f"S octave decay {f} Hz, largest |error| (%)", max(v), max(rf), 0.5, 3.0, np.percentile(v, 90))
    v = [summarise(r["m"])[f"oct_dev{f}"] or 0 for k in GATE_S for r in clean(OUT, k) for f in OCT]
    rf = [summarise(r["m"])[f"oct_dev{f}"] or 0 for k in GATE_S for r in clean(REF, k) for f in OCT]
    check("S largest EDC difference to -30 dB, any octave (dB)", max(v), max(rf), 0.1, 0.5, np.percentile(v, 90))
    U, UR = json.load(open(os.path.join(OUT, "units.json"))), json.load(open(os.path.join(REF, "units.json")))
    for name in ["CLA"] + [f"Altiverb RT {k}" for k in GATE_UNITS]:
        for w in ("50-150", "300-600"):
            check(f"{name} R null {w} ms (dB)", U[name]["null"][w][1], UR[name]["null"][w][1], 0.5)
        v = max(abs(summarise(U[name]["m"])[f"oct_t{f}"] or 0) * 100 for f in OCT)
        check(f"{name} octave decay, largest |error| (%)", v, max(abs(summarise(UR[name]["m"])[f"oct_t{f}"] or 0) * 100 for f in OCT), 0.5, 3.0)
    for prog in "ML":
        M_ = lambda root: [json.load(open(os.path.join(root, "ml", f"{prog}_{k}.json"))) for k in GATE_ML[prog]]
        for w in ("first_pass", "recirc1"):
            v = [r["null"][w] for d in M_(OUT) for r in d["rows"] if r["null"].get(w) is not None]
            rf = [r["null"][w] for d in M_(REF) for r in d["rows"] if r["null"].get(w) is not None]
            check(f"{prog} median null {w} (dB)", np.median(v), np.median(rf), 0.5)
        v = [abs(summarise(r["m"])[f"oct_t{f}"] or 0) * 100 for d in M_(OUT) for r in d["rows"] for f in OCT]
        rf = [abs(summarise(r["m"])[f"oct_t{f}"] or 0) * 100 for d in M_(REF) for r in d["rows"] for f in OCT]
        check(f"{prog} early octave decay (-5/-25 dB), median |error| (%)", np.median(v), np.median(rf), 0.3, 3.0)
    lines.append(f"**{'PASS' if ok else 'FAIL'}**")
    open(os.path.join(OUT, "gate.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------- report ----------
def _ratio(x):
    return None if x is None else x[0] / x[1] - 1


def summarise(m):
    """One comparison -> scalar summaries (R channel; L in the *_L fields)."""
    s = {}
    for c, tag in ((1, ""), (0, "_L")):
        for band, fs_ in (("oct", OCT), ("third", THIRD)):
            for f in fs_:
                b = m[band][f if f in m[band] else str(f)][c]
                s[f"{band}_t{f}{tag}"] = _ratio(b["t"])
                s[f"{band}_edt{f}{tag}"] = _ratio(b["edt"])
                s[f"{band}_late{f}{tag}"] = _ratio(b.get("late"))
                s[f"{band}_dev{f}{tag}"] = b["edc_dev"]
        if "fine" in m:
            for w, x in m["fine"].items():
                s[f"fine_{w}{tag}"] = x[c]
    s["ned"] = m["ned"]["dev"]
    if "rv" in m:                                   # G6: differences A - B (model - capture), ratios for times
        a, b = m["rv"]["a"], m["rv"]["b"]
        le = [b[k] for k in b if k.startswith("lev_e_")]
        ll = [b[k] for k in b if k.startswith("lev_l_")]
        empty = bool(le and ll) and np.mean(le) < np.mean(ll) - 40   # no early sound (L at AMPL 0: no cluster before the
        early = ("c50_", "c80_", "d50_", "ts_", "lev_e_", "iacc_e", "lr_e")   # 86 ms onset): early ratios are meaningless
        for k in a:
            if k in ("edr", "ned1") or k not in b or a[k] is None or b[k] is None or (empty and k.startswith(early)):
                continue
            s[f"rv_{k}"] = (a[k] / b[k] - 1) if k == "mix_ms" else a[k] - b[k]
            s[f"rvv_{k}"] = b[k]
        if a.get("edr") and b.get("edr"):
            ta, tb = np.array([np.nan if x is None else x for x in a["edr"]]), np.array([np.nan if x is None else x for x in b["edr"]])
            q = np.abs(ta / tb - 1)
            q = q[np.isfinite(q)]
            if len(q) > 20:
                s["rv_edr_p90"], s["rv_edr_max"] = float(np.percentile(q, 90)), float(np.max(q))

            def ring(t):                            # the ringiest band: its decay time over its 1/3-octave neighbourhood's median
                v = [t[i] / np.nanmedian(t[max(0, i - 4):i + 5]) for i in range(len(t)) if np.isfinite(t[i])]
                return float(np.nanmax(v)) if v else np.nan
            s["rv_ring"], s["rvv_ring"] = ring(ta) - ring(tb), ring(tb)
    if m.get("nbenv"):
        q = np.array([np.nan if x is None else x for x in m["nbenv"]])
        if np.isfinite(q).sum() > 20:
            s["rv_nb_p90"], s["rv_nb_max"] = float(np.nanpercentile(q, 90)), float(np.nanmax(q))
    if "fine_sd" in m:
        for w, (x, y) in m["fine_sd"].items():
            s[f"rv_fsd_{w}"], s[f"rvv_fsd_{w}"] = x - y, y
    if "coh" in m:
        for f in THIRD:
            a, b = m["coh"]["a"][str(f)], m["coh"]["b"][str(f)]
            s[f"coh_{f}"] = a - b
            s[f"cohv_{f}"] = b
    return s


def g6_lines(rows, reps, indent="- "):
    """G6 summary lines: median model - capture [capture - capture] (lists of summarise dicts)."""
    def md(rr, k, f=lambda x: x):
        x = [f(r[k]) for r in rr if r.get(k) is not None and np.isfinite(r[k])]
        return np.median(x) if x else np.nan
    cell = lambda k, f=abs, fmt="{:.2f}", sc=1: f"{fmt.format(sc * md(rows, k, f))} [{fmt.format(sc * md(reps, k, f))}]"
    out = [indent + "C50 |diff| dB per octave: " + " ".join(f"{f}: {cell(f'rv_c50_{f}')}" for f in OCT),
           indent + "C80 |diff| dB per octave: " + " ".join(f"{f}: {cell(f'rv_c80_{f}')}" for f in OCT),
           indent + "C80 signed dB (capture value): " + " ".join(f"{f}: {md(rows, f'rv_c80_{f}'):+.2f} ({md(rows, f'rvv_c80_{f}'):.1f})" for f in OCT),
           indent + "D50 |diff| per octave: " + " ".join(f"{f}: {cell(f'rv_d50_{f}', fmt='{:.3f}')}" for f in OCT),
           indent + "centre time Ts |diff| ms per octave (capture value): " + " ".join(f"{f}: {cell(f'rv_ts_{f}', fmt='{:.1f}')} ({md(rows, f'rvv_ts_{f}'):.0f})" for f in OCT),
           indent + "tonal balance, third-octave level |diff| dB, early 0-80 ms: " + " ".join(f"{f}: {cell(f'rv_lev_e_{f}', fmt='{:.2f}')}" for f in THIRD[::2]),
           indent + "tonal balance, late 80 ms on: " + " ".join(f"{f}: {cell(f'rv_lev_l_{f}', fmt='{:.2f}')}" for f in THIRD[::2]),
           indent + "tonal balance signed dB (model - capture), early / late: " + " ".join(f"{f}: {md(rows, f'rv_lev_e_{f}'):+.2f} / {md(rows, f'rv_lev_l_{f}'):+.2f}" for f in THIRD[::3]),
           indent + f"EDR, decay per 1/24 octave 100 Hz-10 kHz, |ratio - 1| %: 90th percentile over bands {cell('rv_edr_p90', fmt='{:.1f}', sc=100)}, "
                    f"worst band {cell('rv_edr_max', fmt='{:.1f}', sc=100)}; ringiest band's decay / its 1/3-octave median, model - capture "
                    f"{cell('rv_ring', f=lambda x: x, fmt='{:+.3f}')} (capture {md(rows, 'rvv_ring'):.3f})",
           indent + f"narrow-band envelope (1/24 octave, first 30 dB of decay), largest |level difference| dB: 90th percentile over bands "
                    f"{cell('rv_nb_p90', fmt='{:.2f}')}, worst band {cell('rv_nb_max', fmt='{:.2f}')}",
           indent + "tail coloration, fine-structure sd dB, model - capture (capture value): " + " ".join(f"{lo}-{hi}: {cell(f'rv_fsd_{lo}-{hi}', f=lambda x: x, fmt='{:+.2f}')} ({md(rows, f'rvv_fsd_{lo}-{hi}'):.2f})" for lo, hi in FS_BANDS),
           indent + f"mixing time (NED >= 0.9), model / capture - 1 %: {cell('rv_mix_ms', fmt='{:.1f}', sc=100)} (capture {md(rows, 'rvv_mix_ms'):.0f} ms)",
           indent + f"IACC |diff| early {cell('rv_iacc_e', fmt='{:.3f}')} (capture {md(rows, 'rvv_iacc_e'):.3f}), late {cell('rv_iacc_l', fmt='{:.3f}')} (capture {md(rows, 'rvv_iacc_l'):.3f}); "
                    "late per octave: " + " ".join(f"{f}: {cell(f'rv_iacc_l_{f}', fmt='{:.3f}')}" for f in (500, 1000, 2000, 4000)),
           indent + f"L/R level |diff| dB early {cell('rv_lr_e')}, late {cell('rv_lr_l')}"]
    return out


def report_stage():
    D = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(OUT, "s", "*.json")))]
    rows = [(d, summarise(r["m"]), r) for d in D for r in d["rows"]]
    reps = [(d, summarise(p["m"]), p) for d in D for p in d["pairs"] if p["cross"]]
    intk = [(d, summarise(p["m"]), p) for d in D for p in d["pairs"] if not p["cross"]]
    out = [f"## Pulse check, program S: {len(D)} settings, {len(rows)} captures held out by RT step (§32.2 baseline); model vs capture "
           f"[cross-take capture vs capture, {len(reps)} pairs] (in-take pairs in braces where shown)"]

    def med(rr, k, f=lambda x: x):
        x = [f(s[k]) for _, s, _ in rr if s.get(k) is not None and np.isfinite(s[k])]
        return np.median(x) if x else np.nan

    groups = [("All", lambda d: True), ("RT <= 0.6", lambda d: d["rt"] <= 0.6), ("RT 0.8-1.4", lambda d: 0.8 <= d["rt"] <= 1.4),
              ("RT >= 1.8", lambda d: d["rt"] >= 1.8), ("LF 0.1", lambda d: d["lf"] == 0.1), ("LF 2.0", lambda d: d["lf"] == 2.0),
              ("HF 0.2", lambda d: d["hf"] == 0.2), ("HF 1.0", lambda d: d["hf"] == 1.0)]
    sel = lambda rr, g: [x for x in rr if g(x[0])]
    out.append("\n### G2: octave decay time from the EDC (-5 to -25 dB), |ratio - 1| in %, R; median")
    out.append("| Group | " + " | ".join(f"{f}" for f in OCT) + " |\n|---|" + "---|" * len(OCT))
    for g, fn in groups:
        out.append(f"| {g} | " + " | ".join(f"{100 * med(sel(rows, fn), f'oct_t{f}', abs):.1f} [{100 * med(sel(reps, fn), f'oct_t{f}', abs):.1f}]" for f in OCT) + " |")
    out.append("\n### G2: signed octave decay-time ratio - 1 (model / capture), %, R; median (bias)")
    out.append("| Group | " + " | ".join(f"{f}" for f in OCT) + " |\n|---|" + "---|" * len(OCT))
    for g, fn in groups:
        out.append(f"| {g} | " + " | ".join(f"{100 * med(sel(rows, fn), f'oct_t{f}'):+.1f} [{100 * med(sel(reps, fn), f'oct_t{f}'):+.1f}]" for f in OCT) + " |")
    out.append("\n### G2: EDT |ratio - 1| %, and largest EDC difference to -30 dB (dB), R; median")
    out.append("| Group | " + " | ".join(f"{f}" for f in OCT) + " |\n|---|" + "---|" * len(OCT))
    for g, fn in groups:
        out.append(f"| {g} EDT | " + " | ".join(f"{100 * med(sel(rows, fn), f'oct_edt{f}', abs):.1f} [{100 * med(sel(reps, fn), f'oct_edt{f}', abs):.1f}]" for f in OCT) + " |")
        out.append(f"| {g} EDC | " + " | ".join(f"{med(sel(rows, fn), f'oct_dev{f}'):.2f} [{med(sel(reps, fn), f'oct_dev{f}'):.2f}]" for f in OCT) + " |")
    out.append("\n### G3: third-octave decay time, signed ratio - 1 (%), R, all settings: median model [repeat]; and the share of captures "
               "whose |model error| is within the setting's largest cross-take |difference|")
    within = {}
    for f in THIRD:
        ok, n = 0, 0
        for d in D:
            rep = [abs(summarise(p["m"])[f"third_t{f}"] or np.nan) for p in d["pairs"] if p["cross"]]
            tol = np.nanmax(rep) if np.any(np.isfinite(rep)) else np.nan
            for r in d["rows"]:
                x = summarise(r["m"])[f"third_t{f}"]
                if x is None or not np.isfinite(tol):
                    continue
                n += 1; ok += abs(x) <= max(tol, 0.03)
        within[f] = (ok, n)
    out.append(" | ".join(f"{f}: {100 * med(rows, f'third_t{f}'):+.1f} [{100 * med(reps, f'third_t{f}'):+.1f}] {within[f][0]}/{within[f][1]}" for f in THIRD))
    out.append("\n### G4: echo density mean |difference| 20-300 ms (R) and L/R coherence (model - capture) per third octave, 50-300 ms")
    for g, fn in groups:
        out.append(f"- {g}: NED {med(sel(rows, fn), 'ned'):.3f} [{med(sel(reps, fn), 'ned'):.3f}] {{{med(sel(intk, fn), 'ned'):.3f}}}")
    out.append("- coherence, all: " + "; ".join(f"{f}: {med(rows, f'coh_{f}'):+.3f} [{med(reps, f'coh_{f}'):+.3f}] (capture {med(rows, f'cohv_{f}'):.3f})" for f in THIRD))
    out.append("\n### G5: spectral fine-structure correlation 50-300 ms, R (L)")
    for g, fn in groups:
        out.append(f"- {g}: " + "; ".join(f"{w}: {med(sel(rows, fn), f'fine_{w}'):.3f} ({med(sel(rows, fn), f'fine_{w}_L'):.3f}) [{med(sel(reps, fn), f'fine_{w}'):.3f}] {{{med(sel(intk, fn), f'fine_{w}'):.3f}}}"
                                          for w in (f"{lo}-{hi}" for lo, hi in FS_BANDS)))
    out.append("\n### G6: reverb-specific metrics, R unless stated, median model - capture [cross-take capture - capture]")
    for g, fn in groups:
        out.append(f"- {g}:")
        out += g6_lines([x[1] for x in sel(rows, fn)], [x[1] for x in sel(reps, fn)], indent="  - ")
    out.append("\n### G1: nulls 100 Hz-10 kHz, R, median model [cross-take repeat] {in-take}")
    for g, fn in groups:
        cells = []
        for w in ("50-150", "150-300", "300-600", "600-1000"):
            m_ = [r["null"][w][1] for d_, _, r in sel(rows, fn) if r["null"].get(w)]
            p_ = [r["null"][w][1] for d_, _, r in sel(reps, fn) if r["null"].get(w)]
            i_ = [r["null"][w][1] for d_, _, r in sel(intk, fn) if r["null"].get(w)]
            cells.append(f"{np.median(m_):.1f} [{np.median(p_):.1f}] {{{np.median(i_):.1f}}}" if m_ and p_ and i_ else "--")
        out.append(f"- {g}: " + "; ".join(f"{w}: {c}" for w, c in zip(("50-150", "150-300", "300-600", "600-1000"), cells)))
    out.append("\n### G1: octave nulls (R), all settings, median per window 50-150 / 150-300 / 300-600")
    for f in OCT:
        out.append(f"- {f}: " + " / ".join(f"{np.median([r['null_oct'][str(f)][j] for _, _, r in rows if r['null_oct'][str(f)][j] is not None]):.1f}" for j in range(3)))
    out += units_report() + units_report(chain=True)
    out += ml_report()
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["s", "units", "ml", "gate", "report"], required=True)
    ap.add_argument("--only", default=None, help="comma-separated setting keys (lf_rt_hf)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chain", action="store_true", help="units: the model through CLA's capture chain relative to Tomkin's (chain_rel)")
    ap.add_argument("--tag", default="check", help="gate: output folder suffix (analysis/pulse/gate_<tag>)")
    ap.add_argument("--prog", default="ML", help="ml: programs, e.g. M, L or ML")
    ap.add_argument("--n", type=int, default=None, help="ml: smoke test, this many groups per program")
    a = ap.parse_args()
    if a.stage == "s":
        s_stage(a.workers, a.only.split(",") if a.only else None)
    elif a.stage == "units":
        units_stage(a.chain)
    elif a.stage == "gate":
        gate_stage(a.workers, a.tag)
    elif a.stage == "ml":
        ml_stage(list(a.prog), a.workers, a.n)
    else:
        report_stage()


if __name__ == "__main__":
    main()
