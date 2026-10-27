"""NEXT_STEPS §0 'Next steps (S low band, 27 October)': the S / M / L models over the plugin's whole band, 20 Hz to the
top of the captures (12.3 kHz in the S / M / L caches, the deconvolution's band window; above it see --stage wide).
Every null so far is 100 Hz-10 kHz. Results: docs/model_s_results.md §41.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/memguard.py --limit-gb 8 -- \\
        .venv/bin/python tools/sml_fullband.py --stage s [--only keys] [--workers 8]
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds \\
        SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_pallf/cells.json SHELF_TAG=_pallf \\
        EMT_AP230=ap230_g0.7071_chain2.json .venv/bin/python tools/sml_fullband.py --stage takefit [--workers 8]
    EMT_KERNEL=chain3 .venv/bin/python tools/sml_fullband.py --stage ml [--prog SML] [--every 2] [--workers 6]
    EMT_KERNEL=chain3 .venv/bin/python tools/sml_fullband.py --stage warp
    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/sml_fullband.py --stage report

  s        Program S, the bench's model (pulse_check.setting_G: the §38.1 fold set, held out by RT step) and the bench's
           per-capture fit (q2_injection.fit_capture: clock, shift, gain on 20 Hz-11 kHz), every clean capture
           (tools/defects.py). The model gets the pipeline's own response above 11 kHz: the reverb x P(f) (the
           deconvolution's band window and the resampling to the device rate, an impulse through both, as
           outofband.py), the cluster x the resampling's part only (its atoms already carry the window). Then
           (a) strict nulls per band (BANDS, 3rd-order zero-phase band-passes) and window (WINS), with the cross-take
               repeat (AMPL 0-30 and 60-90, as pulse_check) and the take-to-take correlation of the model's residual
               (the residual of an AMPL 0 / 90 capture against that of the AMPL 30 / 60 capture of the same pair:
               near 1 if the model's error is one fixed waveform in both takes, i.e. deterministic and fittable);
           (b) the same as spectral band energies in sixth octaves (FINE, Hann window over FWINS): where in frequency
               the model's excess over the repeat starts and stops.
           One JSON per setting in analysis/fullband/s/ (resumable).
  takefit  Where model_s_results §32.9's per-setting refit (12 numbers per setting, one take fitted, the other scored,
           lowband_shelf --stage takefit on the _pallf baseline, its stored fits) gains, per band and window: the
           same captures and model machinery as that run (lowband_shelf.prep / takefit_base, run in its environment,
           above), base and refit scored in BANDS. -> analysis/fullband/takefit/
  ml       Programs M and L (and S through the same model family): ml_lines.Model as outofband.py (the law stage's
           in-band fit per capture, rendered full band x P(f), the repeat capture fitted the same way), per band in
           cluster / gap, first pass and recirculation 1, and the sixth-octave profile over the first pass.
           Every --every-th group. -> analysis/fullband/ml/
  warp     The clock warps the S / M / L tools apply (q1.warp at 128 / beta 12 and at 32 / beta 8, ml_lines.warp
           48 / beta 10), against a reference warp (512 / beta 14), on a real S capture and a model render at the
           library's clock ratios: the artefact per band and window, as dB below the band's energy (NonLin_program
           §11.6 found nl_late.unwarp's ×32 Kaiser(5) upsampler leaves a growing error at 11-12.8 kHz). -> warp.json
  report   Medians -> analysis/fullband/report.{log,json}.
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import pulse_check as pc  # noqa: E402
from core import _band_window, to_device_grid  # noqa: E402
from defects import DEFECTS  # noqa: E402

q2, q1 = pc.q2, pc.q1
FS, PRE = q1.FS, q1.PRE
OUT = os.path.join(ROOT, "analysis", "fullband")
EDGES = [20, 40, 80, 160, 270, 400, 800, 1600, 3200, 6400, 10000, 11000, 12300]
BANDS = [(a, b) for a, b in zip(EDGES[:-1], EDGES[1:])]
BNAME = [f"{a}-{b}" for a, b in BANDS]
SOS = [butter(3, (a, b), btype="band", fs=FS, output="sos") for a, b in BANDS]
WINS = [(0, 50), (50, 150), (150, 300), (300, 600), (600, 1000)]
FINE = 20.0 * 2 ** (np.arange(0, 6 * np.log2(12300 / 20) + 1) / 6)   # sixth-octave edges, 20 Hz-12.3 kHz
FWINS = [(50, 150), (150, 300), (300, 600)]
CROSS = pc.CROSS                                                       # (0, 30), (60, 90): the cross-take pairs


# ---------- the pipeline's response above 11 kHz ----------
def pipeline(N, window=True):
    """On an N-point device-rate rfft grid: a unit impulse at 96 kHz through the deconvolution's band window (or not)
    and to_device_grid, over the ideal device-rate impulse; unit gain at 1 kHz (outofband.pipeline_response)."""
    fs = 96000.0
    n = 1 << 16
    x = np.zeros(n)
    i0 = 20000
    x[i0] = 1.0
    if window:
        x = np.fft.irfft(np.fft.rfft(x) * _band_window(np.fft.rfftfreq(n, 1 / fs)), n)
    k0 = int(np.floor(i0 / fs * FS))
    h = to_device_grid(x, fs, float(i0), k0, q1.ms.CFG)
    P = np.fft.rfft(h, N) * np.exp(2j * np.pi * np.fft.rfftfreq(N) * k0)
    return P / np.abs(P[int(round(1000 / FS * N))])


def apply_f(x, H):
    """x (n, 2) through the response H on an rfft grid, one (N/2+1,) or per channel (N/2+1, 2)."""
    N = 2 * (len(H) - 1)
    return np.fft.irfft(np.fft.rfft(x, N, axis=0) * (H[:, None] if H.ndim == 1 else H), N, axis=0)[:len(x)]


T_REV = os.environ.get("SML_TREV") == "1"         # tfix: the transfer on the reverb (loop paths) only, not the cluster
T_PATH = os.environ.get("SML_TPATH", "")           # tfix: "first" = on each loop path's first pass only (not the recirculations)


def full_model(J, v, G, n, amp, mdl, T=None):
    """mdl (fit_capture's model) with the pipeline response: cluster x resampling, reverb x (window x resampling).
    T: a per-channel transfer (2, TNB) on the whole model, or with SML_TREV=1 on the reverb only."""
    cl, f, xs = J.split(v)
    F_, X_ = J._spectra(f, xs, G["N"])
    r = J.reverb(G, F_, X_)[:n]
    c = mdl - r
    N = 1 << int(np.ceil(np.log2(n + 8192)))
    Pc, Pr = pipeline(N, False)[:, None] * np.ones((1, 2)), pipeline(N, True)[:, None] * np.ones((1, 2))
    if T is not None and T_PATH == "first":          # the transfer on each loop path's first pass only (X_i PB, before any
        rf = J.reverb(dict(G, R=[G["PB"]] * len(G["R"])), F_, X_)[:n]   # recirculation; S's first-pass FIR is zero)
        return apply_f(c, Pc) + apply_f(r - rf, Pr) + apply_f(rf, Pr * t_response(T, N))
    if T is not None:
        H = t_response(T, N)
        Pr = Pr * H
        if not T_REV:
            Pc = Pc * H
    return apply_f(c, Pc) + apply_f(r, Pr)


# ---------- energies per band and window ----------
def win_idx(a_ms, b_ms, te, n):
    i0, i1 = PRE + int(a_ms / 1000 * FS), min(PRE + int(b_ms / 1000 * FS), PRE + int(te * FS), n)
    return (i0, i1) if i1 - i0 >= int(0.05 * FS) else None


def band_energies(sigs, te):
    """sigs: {name: (n, 2)}; per band and window, per name [E_L, E_R] (None where the window is under 50 ms)."""
    n = min(len(x) for x in sigs.values())
    out = {}
    for b, sos in zip(BNAME, SOS):
        f = {k: sosfiltfilt(sos, x[:n], axis=0) for k, x in sigs.items()}
        row = {}
        for a_ms, b_ms in WINS:
            ix = win_idx(a_ms, b_ms, te, n)
            if ix is None:
                continue
            i0, i1 = ix
            row[f"{a_ms}-{b_ms}"] = {k: np.sum(x[i0:i1] ** 2, 0).tolist() for k, x in f.items()}
            if "ea" in f and "eb" in f:
                row[f"{a_ms}-{b_ms}"]["cab"] = np.sum(f["ea"][i0:i1] * f["eb"][i0:i1], 0).tolist()
        out[b] = row
    return out


def fine_energies(sigs, te, wins=FWINS, k0=PRE):
    """Hann-windowed spectra per window: band energies in the sixth octaves FINE; 'cab' the real cross spectrum of ea, eb."""
    n = min(len(x) for x in sigs.values())
    out = {}
    for a_ms, b_ms in wins:
        i0, i1 = k0 + int(a_ms / 1000 * FS), min(k0 + int(b_ms / 1000 * FS), k0 + int(te * FS), n)
        if i1 - i0 < int(0.05 * FS):
            continue
        h = np.hanning(i1 - i0)[:, None]
        N = 1 << int(np.ceil(np.log2(4 * (i1 - i0))))
        f = np.fft.rfftfreq(N, 1 / FS)
        idx = np.digitize(f, FINE) - 1
        ok = (idx >= 0) & (idx < len(FINE) - 1)
        S = {k: np.fft.rfft(x[i0:i1] * h, N, axis=0) for k, x in sigs.items()}
        row = {}
        for k, X in S.items():
            row[k] = [np.bincount(idx[ok], np.abs(X[ok, c]) ** 2, len(FINE) - 1).tolist() for c in (0, 1)]
        if "ea" in S and "eb" in S:
            row["cab"] = [np.bincount(idx[ok], np.real(S["ea"][ok, c] * np.conj(S["eb"][ok, c])), len(FINE) - 1).tolist() for c in (0, 1)]
        out[f"{a_ms}-{b_ms}"] = row
    return out


# ---------- stage s ----------
TN, TBIN, TTOP = 1 << 16, 2.0, 12300.0      # transfer: rfft length (2.5 s), bin width (Hz), top
TNB = int(TTOP / TBIN)
TFOLDS = {"A": (0.3, 0.8, 1.4, 2.2), "B": (0.6, 1.0, 1.8, 2.5)}   # tfix: T fitted on one fold's RT steps, applied to the other's


def transfer_acc(y, mg, te, t0=0.15):
    """Over the reverb [t0, te]: per 2 Hz bin and channel, sum Y M*, |M|^2, |Y|^2 (M the model x fade ramp)."""
    i0, i1 = PRE + int(t0 * FS), min(PRE + int(te * FS), len(y))
    Y, M = np.fft.rfft(y[i0:i1], TN, axis=0), np.fft.rfft(mg[i0:i1], TN, axis=0)
    idx = (np.fft.rfftfreq(TN, 1 / FS) / TBIN).astype(int)
    ok = idx < TNB
    acc = lambda z: np.stack([np.bincount(idx[ok], z[ok, c], TNB) for c in (0, 1)])
    ym = Y * np.conj(M)
    return acc(ym.real) + 1j * acc(ym.imag), acc(np.abs(M) ** 2), acc(np.abs(Y) ** 2)


def t_response(T, N):
    """A per-channel transfer on the 2 Hz grid (2, TNB) as an rfft-grid response (N, 2): interpolated, 1 outside 15 Hz-12.3 kHz."""
    f = np.fft.rfftfreq(N, 1 / FS)
    fc = (np.arange(TNB) + 0.5) * TBIN
    H = np.ones((len(f), 2), complex)
    m = (f >= 15) & (f <= TTOP - TBIN)
    for c in (0, 1):
        H[m, c] = np.interp(f[m], fc, T[c].real) + 1j * np.interp(f[m], fc, T[c].imag)
    return H


def refit_capture(mdl, c, ends, n):
    """q2_injection.fit_capture's clock / shift / gain fit against a given model (n samples): y, fade ramp g, te."""
    yb = q2._band(c["file"])[:n]
    if c["fade2"] >= 0.15:
        _, (e, sh, gain) = q1._clock_job((yb, mdl, c["fade2"], c["eps"], c["amp"] == 0))
    else:
        g0 = q2.ramp_gain(n, ends[c["file"]])
        _, (e, sh, gain) = q1._clock_job((yb, mdl * g0, q2.comp_end(yb, g0), c["eps"], c["amp"] == 0))
    y = q1.apply_clock(yb, e, sh)[:n] / gain
    g = q2.ramp_gain(n, [PRE + (x + sh - PRE) / (1 + e * 1e-6) for x in ends[c["file"]]])
    return y, g, q2.comp_end(y, g)


def hf_correct(m, c):
    """sum_k c[ch, k] HP(m)(n + k) per channel (the kfix form: KF_HPS zero-phase, c on m(n + k), k = -K..K)."""
    K = (c.shape[1] - 1) // 2
    mh = sosfiltfilt(KF_HPS, m, axis=0)
    out = np.zeros_like(m)
    for ch in (0, 1):
        for j, k in enumerate(range(-K, K + 1)):
            out[:, ch] += c[ch, j] * np.roll(mh[:, ch], -k)
    return out


SHF_WIN = (0, 50)                                    # shf: the window the S correction is fitted on (ms after tap 0)


def _shf_job(a):
    """One S setting: per clean capture the 0-50 ms window (+-K) of the band-limited capture, the full model and its
    high-passed copy, both channels (as _s_job's model and fit)."""
    key, lab, b230, law, fold, fc, v, caps, ends, K = a
    q1.EXT = 1700
    G = pc.setting_G(lab, b230, law, fold, fc)
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
    out = []
    for c in caps:
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        mf = full_model(J, v, G, len(y), c["amp"], mdl) * g
        yf, mff = sosfiltfilt(KF_BAND, y, axis=0), sosfiltfilt(KF_BAND, mf, axis=0)
        mh = sosfiltfilt(KF_HPS, mff, axis=0)
        i0, i1 = PRE + int(SHF_WIN[0] / 1000 * FS), PRE + int(SHF_WIN[1] / 1000 * FS)
        out.append((key, c["amp"]) + tuple(x[i0 - K:i1 + K].astype(np.float32) for x in (yf, mff, mh)))
    return out


def shf_fit(workers, every=1):
    """S's top octave in the cluster window: per channel an 11-tap correction above KF_HP fitted on one take's 0-50 ms
    (AMPL 0 + 90 or 30 + 60, every setting), least squares on the full-band residual -> analysis/fullband/shf_fit.json."""
    from concurrent.futures import ProcessPoolExecutor
    K = 5
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    folds = {rt: json.load(open(os.path.join(pc.CV, f"fold_{rt}.json"))) for rt in sorted({g["rt"] for g in info.values()})}
    jobs = []
    for k in sorted(info)[::every]:
        g = info[k]
        caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(g) if (k, a_) not in DEFECTS]
        if caps:
            jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], law, folds[g["rt"]], cells.get(k), v, caps, ends, K))
    rows = []
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for r in ex.map(_shf_job, jobs):
            rows += r
    seg = lambda m: np.stack([np.roll(m, -k)[K:-K] for k in range(-K, K + 1)], 1)
    res = {}
    for take, amps in (("A", (0, 90)), ("B", (30, 60))):
        c = np.zeros((2, 2 * K + 1))
        for ch in (0, 1):
            A = np.vstack([seg(mh[:, ch].astype(float)) for _, a_, y, m, mh in rows if a_ in amps])
            b = np.concatenate([(y[:, ch].astype(float) - m[:, ch].astype(float))[K:-K] for _, a_, y, m, mh in rows if a_ in amps])
            c[ch] = np.linalg.lstsq(A, b, rcond=None)[0]
        res[take] = c.tolist()
    res["meta"] = dict(win_ms=SHF_WIN, hp=dict(order=KF_ORDER, hz=KF_HP), taps=2 * K + 1, settings=len(jobs), captures=len(rows))
    json.dump(res, open(os.path.join(OUT, "shf_fit.json"), "w"), indent=1)
    print(json.dumps({k: np.round(np.array(x), 3).tolist() if k != "meta" else x for k, x in res.items()}), flush=True)


def _s_job(a):
    key, lab, b230, law, fold, fc, v, caps, ends, path = a[:10]
    T = a[10] if len(a) > 10 else None
    HFC = a[11] if len(a) > 11 else None                 # shf: {take: (2, 2K+1) correction above KF_HP}, applied held out by take
    q1.EXT = 1700                                    # before any Joint: its excitation layout (fit_capture's ext)
    G = pc.setting_G(lab, b230, law, fold, fc)
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
    P = {}
    acc = [0, 0, 0]
    for c in caps:
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        n = len(y)
        mf = full_model(J, v, G, n, c["amp"], mdl, T)
        if HFC is not None:                              # the correction fitted on the other take
            cc_ = HFC["B" if c["amp"] in (0, 90) else "A"]
            mf = mf + hf_correct(mf, cc_)
            if T_REFIT:
                y, g, te = refit_capture(mf, c, ends, n)
        if T is not None and T_REFIT:                # as the bench: clock, shift and gain fitted against the corrected model
            y, g, te = refit_capture(mf, c, ends, n)
        if T is None:
            acc = [x + z for x, z in zip(acc, transfer_acc(y, mf * g, te))]
        P[c["amp"]] = dict(y=y, m=mf, g=g, te=te, e=y - mf * g, check=q1.nulls(y, mdl * g, te, q2.W_COMP, narrow=True))
    if T is None:
        np.savez(path[:-5] + ".npz", sym=acc[0].astype(np.complex64), smm=acc[1].astype(np.float32), syy=acc[2].astype(np.float32))
    out = dict(key=key, **lab, caps={}, pairs={})
    for amp, d in P.items():
        sig = dict(y=d["y"], e=d["e"])
        out["caps"][amp] = dict(te=d["te"], check=d["check"], bands=band_energies(sig, d["te"]), fine=fine_energies(sig, d["te"]))
    for a_, b_ in CROSS:
        if a_ not in P or b_ not in P:
            continue
        A, B = P[a_], P[b_]
        n = min(len(A["y"]), len(B["y"]))
        te = min(A["te"], B["te"])
        sig = dict(ya=A["y"][:n] * B["g"][:n], d=A["y"][:n] * B["g"][:n] - B["y"][:n] * A["g"][:n], ea=A["e"][:n], eb=B["e"][:n])
        out["pairs"][f"{a_}-{b_}"] = dict(te=te, bands=band_energies(sig, te), fine=fine_energies(sig, te))
    json.dump(out, open(path, "w"))
    print(key, flush=True)
    return key


def fold_T(fold, smooth=2):
    """The fixed transfer capture / model fitted on the settings of RT fold `fold` at LF >= 1.0 (no LF shelf, so no decay
    law in it), pooled over settings and captures, smoothed over +-smooth bins: (2, TNB) complex."""
    info = q2.q1_tables()[0]
    sym, smm = 0, 0
    for k, g in info.items():
        p = os.path.join(OUT, "s", k + ".npz")
        if g["rt"] in TFOLDS[fold] and g["lf"] >= 1.0 and os.path.exists(p):
            z = np.load(p)
            sym, smm = sym + z["sym"].astype(complex), smm + z["smm"].astype(float)
    k_ = np.ones(2 * smooth + 1)
    num = np.stack([np.convolve(sym[c], k_, "same") for c in (0, 1)])
    den = np.stack([np.convolve(smm[c], k_, "same") for c in (0, 1)])
    return num / np.maximum(den, 1e-30)


T_FORM = os.environ.get("SML_TFORM", "")          # tfix: a compact form of the fold transfer (T_FORMS), else the free one
T_REFIT = os.environ.get("SML_TREFIT") == "1"     # tfix: each capture's clock, shift and gain refitted with T in the model
T_FORMS = {"short": dict(lags=10), "short_ap": dict(lags=10, echo=(230, 460)), "short_ap_lf": dict(lags=10, echo=(230, 460), lf=100.0)}


def compact_T(T, form):
    """T's impulse response (T - 1, both channels) kept only at |lag| <= lags, at +-echo lags (+-2 each), plus the part of
    T - 1 below lf Hz (a smooth low-band term: lf / TBIN complex bins); back on the 2 Hz grid. Returns (T_c, numbers)."""
    cfg = T_FORMS[form]
    N = 1 << 16
    H = t_response(T, N) - 1.0
    lag = np.r_[np.arange(N // 2 + 1), np.arange(-N // 2 + 1, 0)]
    keep = np.abs(lag) <= cfg["lags"]
    nnum = 2 * cfg["lags"] + 1
    for e in cfg.get("echo", ()):
        keep |= (np.abs(np.abs(lag) - e) <= 2)
        nnum += 10
    fr = np.fft.rfftfreq(N, 1 / FS)
    out = np.zeros_like(H)
    for c in (0, 1):
        h = np.fft.irfft(H[:, c], N)
        out[:, c] = np.fft.rfft(h * keep, N)
        if cfg.get("lf"):
            lo = fr < cfg["lf"]
            out[lo, c] = H[lo, c]                   # below lf: T itself (cos-free switch: the kernel's part there is ~0)
    fc = (np.arange(TNB) + 0.5) * TBIN
    Tc = np.stack([1 + np.interp(fc, fr, out[:, c].real) + 1j * np.interp(fc, fr, out[:, c].imag) for c in (0, 1)])
    if cfg.get("lf"):
        nnum += 2 * int(cfg["lf"] / TBIN)
    return Tc, 2 * nnum


def s_stage(workers, only=None, tfix=False, sub=None):
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, (("s_tfixrev" if T_REV else "s_tfix") + (f"_{T_FORM}" if T_FORM else "") + (f"_{T_PATH}" if T_PATH else "") + ("_refit" if T_REFIT else "")) if tfix
                     else ("s_xcorr" + os.environ.get("SML_XTAG", "")) if os.environ.get("EMT_XCORR") else "s")
    d = os.path.join(OUT, sub) if sub else d
    os.makedirs(d, exist_ok=True)
    Ts = {f: fold_T(f) for f in TFOLDS} if tfix else None
    if os.environ.get("SML_SHF") == "1":            # shf scoring: the correction of the other take on every capture
        z = json.load(open(os.path.join(OUT, "shf_fit.json")))
        HFC = {t: np.array(z[t]) for t in ("A", "B")}
        d = os.path.join(OUT, "s_shf" + ("_refit" if T_REFIT else ""))
        os.makedirs(d, exist_ok=True)
    else:
        HFC = None
    if tfix and T_FORM:
        Ts = {f: compact_T(x, T_FORM)[0] for f, x in Ts.items()}
        print(f"{T_FORM}: {compact_T(fold_T('A'), T_FORM)[1]} numbers", flush=True)
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    folds = {rt: json.load(open(os.path.join(pc.CV, f"fold_{rt}.json"))) for rt in sorted({g["rt"] for g in info.values()})}
    jobs = []
    for k, g in info.items():
        path = os.path.join(d, k + ".json")
        if (only and k not in only) or os.path.exists(path):
            continue
        lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
        caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(g) if (k, a_) not in DEFECTS]
        if not caps:
            continue
        job = (k, lab, ap["steps"][ap["settings"][k]["step"]], law, folds[g["rt"]], cells.get(k), v, caps, ends, path)
        if tfix:                                    # held out: the transfer of the fold that does not contain this RT step
            job = job + (Ts["A"] if g["rt"] in TFOLDS["B"] else Ts["B"],)
        if HFC is not None:
            job = job + ((None,) if not tfix else ()) + (HFC,)
        jobs.append(job)
    print(f"{len(jobs)} settings to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_s_job, jobs):
            pass


# ---------- stage takefit ----------
def _takefit_job(a):
    key, path, src = a
    import lowband_shelf as ls
    f = json.load(open(src))
    d, G, G1, pb, hf1 = ls.takefit_base(key)
    z1 = np.exp(-1j * d["w"])
    caps = [c for c in d["caps"] if (key, c["amp"]) not in DEFECTS]
    takes = {"A": [c for c in caps if c["amp"] in (0, 90)], "B": [c for c in caps if c["amp"] in (30, 60)]}
    nfree = 8 if hf1 else 12
    if not caps:                                     # 1.0_0.3_0.4: no clean capture (tools/defects.py)
        json.dump(dict(key=key, **d["lab"], caps={}), open(path, "w"))
        return key
    r0 = ls.reverb(d, G)
    N = 1 << int(np.ceil(np.log2(max(len(c["y"]) for c in caps) + 8192)))
    PW, PD = pipeline(N, True), pipeline(N, False)
    out = dict(key=key, **d["lab"], caps={})
    for tr, te in (("A", "B"), ("B", "A")):
        if tr not in f["fit"]:
            continue
        u = np.r_[f["fit"][tr], np.zeros(12 - nfree)]
        r1 = ls.reverb(d, ls.takefit_Gs(G1, pb, hf1, z1, d["w"], d["law"]["a"], u))
        for c in takes[te]:
            n = len(c["y"])
            cl = (c["amp"] / 90) * d["clus"][:n]
            m0 = apply_f(cl, PD) + apply_f(r0[:n], PW)
            m1 = apply_f(cl, PD) + apply_f(r1[:n], PW)
            sig = dict(y=c["y"], e0=c["y"] - m0 * c["g"][:n], e1=c["y"] - m1 * c["g"][:n])
            out["caps"][c["amp"]] = dict(te=c["te"], check=ls.nulls_broad(c["y"], ls.model(d, c, r0) * c["g"], c["te"]),
                                         bands=band_energies(sig, c["te"]), fine=fine_energies(sig, c["te"]))
    json.dump(out, open(path, "w"))
    print(key, flush=True)
    return key


def takefit_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    src = os.path.join(ROOT, "analysis", "lowband", "chain3_shelf_takefit_cap0.8475" + os.environ.get("SHELF_TAG", "_pallf"))
    d = os.path.join(OUT, "takefit" + os.environ.get("SHELF_TAG", "_pallf"))
    os.makedirs(d, exist_ok=True)
    jobs = []
    for p in sorted(glob.glob(os.path.join(src, "*.json"))):
        k = os.path.basename(p)[:-5]
        path = os.path.join(d, k + ".json")
        if (only and k not in only) or os.path.exists(path):
            continue
        jobs.append((k, path, p))
    print(f"{len(jobs)} settings to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_takefit_job, jobs):
            pass


# ---------- stage ml ----------
ML_WINS = ("cluster_gap", "first_pass", "recirc1")
ML_TFIX = os.environ.get("SML_TFIX") == "1"


def _ml_job(a):
    prog, key, path = a
    from scipy.optimize import minimize
    import outofband as ob
    cc, mlines, ml = ob.cc, ob.mlines, ob.ml
    M = mlines.Model(prog)
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, amps = ob.captures4(prog, key, n)
    n = len(caps[amps[-1]][0])
    W = mlines.windows(prog)
    Craw, Rraw = M.parts(n, key, band=False)
    C, R = (sosfiltfilt(ml.BAND, x, axis=0) for x in (Craw, Rraw))
    N = 1 << int(np.ceil(np.log2(n + 4096)))
    P = pipeline(N, True)
    if ML_TFIX:                                     # S's fixed transfer (both folds averaged) on the M / L model: the chain?
        P = P[:, None] * t_response(0.5 * (fold_T("A") + fold_T("B")), N)
    rows = []
    for amp in amps:
        if (prog, key, amp) in __import__("defects").ML_SLIPPED:
            continue
        yb, yr = caps[amp]
        end = mlines.usable_end(yb, g[amp], prog)
        sl, nf = mlines.fit_window(prog, amp > 0)
        _, u, gs, (Cw, _) = mlines.fit_warp(yb, [C, R], sl, nf)
        Cwr, Rwr = mlines.warp(Craw, *u), mlines.warp(Rraw, *u)
        cl_raw = apply_f(Cwr * np.array([gs[ch][0] for ch in (0, 1)])[None], P)
        mraw = cl_raw + apply_f(Rwr * np.array([gs[ch][1] for ch in (0, 1)])[None], P)
        other = max((x for x in amps if x != amp), key=lambda x: min(mlines.usable_end(caps[x][0], g[x], prog)))
        cl_b = Cw * np.array([gs[ch][0] for ch in (0, 1)])[None]
        obb = caps[other][0]

        def cost(vv):
            mdl, _ = mlines.gains(yb[:nf], [cl_b[:nf], mlines.warp(obb[:nf], *vv)], sl)
            return sum(np.sum((yb[sl[ch], ch] - mdl[sl[ch], ch]) ** 2) for ch in (0, 1))
        vv = minimize(cost, (0.0, 0.0), method="Nelder-Mead",
                      options=dict(xatol=1e-3, fatol=1e-14, initial_simplex=[(0, 0), (0.2, 0), (0, 10)])).x
        _, gr = mlines.gains(yb[:nf], [cl_b[:nf], mlines.warp(obb[:nf], *vv)[:nf]], sl)
        orw = mlines.warp(caps[other][1], *vv)
        rep = np.stack([gr[ch][0] * cl_raw[:, ch] + gr[ch][1] * orw[:, ch] for ch in (0, 1)], 1)
        oe = mlines.usable_end(obb, g[other], prog)
        e_m, e_r = yr - mraw, yr - rep
        row = dict(key=key, amp=int(amp), bands={}, fine={})
        filt = {b: (sosfiltfilt(s, yr, axis=0), sosfiltfilt(s, e_m, axis=0), sosfiltfilt(s, e_r, axis=0)) for b, s in zip(BNAME, SOS)}
        for w in ML_WINS:
            if w == "cluster_gap" and amp == 0:
                continue
            for ch in (0, 1):
                a_ = cc.K0 + W[ch][w][0]
                b_ = min(cc.K0 + W[ch][w][1], end[ch], oe[ch])
                if b_ - a_ <= 200:
                    continue
                for b, (fy, fm, fr) in filt.items():
                    row["bands"].setdefault(b, {}).setdefault(w, [None, None])[ch] = [float(np.sum(x[a_:b_, ch] ** 2)) for x in (fy, fm, fr)]
                if w != "cluster_gap":
                    h = np.hanning(b_ - a_)
                    NN = 1 << int(np.ceil(np.log2(4 * (b_ - a_))))
                    ff = np.fft.rfftfreq(NN, 1 / FS)
                    idx = np.digitize(ff, FINE) - 1
                    ok = (idx >= 0) & (idx < len(FINE) - 1)
                    row["fine"].setdefault(w, [None, None])[ch] = [
                        np.bincount(idx[ok], np.abs(np.fft.rfft(x[a_:b_, ch] * h, NN)[ok]) ** 2, len(FINE) - 1).tolist() for x in (yr, e_m, e_r)]
        rows.append(row)
    json.dump(rows, open(path, "w"))
    print(prog, key, flush=True)
    return key


def ml_stage(progs, every, workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    import ml_carryover as cc
    import ml_lines as mlines
    d = os.path.join(OUT, ("ml_tfix" if ML_TFIX else "ml") + mlines.SUFFIX + ("_kf" if os.environ.get("ML_KFIX_FILE") else ""))   # ML_HK / ML_FP / ML_KFIX_FILE: own folder
    os.makedirs(d, exist_ok=True)
    jobs = []
    for p in progs:
        ks = sorted(cc.groups(p))[::every]
        if only:
            ks = ks[:only]
        jobs += [(p, k, os.path.join(d, f"{p}_{k}.json")) for k in ks if not os.path.exists(os.path.join(d, f"{p}_{k}.json"))]
    print(f"{len(jobs)} groups to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_ml_job, jobs):
            pass


# ---------- stage kfix: the M / L output kernel over the whole band ----------
def _kfix_job(a):
    """One undamped M / L group: per capture the R channel's full-band capture and full-band model x P(f) (as _ml_job,
    the law stage's in-band clock and gains), over the cluster / gap and first-pass windows plus K samples each side."""
    prog, key, K = a
    import outofband as ob
    cc, mlines, ml = ob.cc, ob.mlines, ob.ml
    M = mlines.Model(prog)
    if M.b.get(key, 0.0) >= 0.05:
        return key, []
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, amps = ob.captures4(prog, key, n)
    n = len(caps[amps[-1]][0])
    Craw, Rraw = M.parts(n, key, band=False)
    C, R = (sosfiltfilt(ml.BAND, x, axis=0) for x in (Craw, Rraw))
    N = 1 << int(np.ceil(np.log2(n + 4096)))
    P = pipeline(N, True)
    W = mlines.windows(prog)[1]
    out = []
    for amp in amps:
        if amp == 0 or (prog, key, amp) in __import__("defects").ML_SLIPPED:
            continue
        yb, yr = caps[amp]
        end = mlines.usable_end(yb, g[amp], prog)
        sl, nf = mlines.fit_window(prog, True)
        _, u, gs, _ = mlines.fit_warp(yb, [C, R], sl, nf)
        mraw = apply_f(mlines.warp(Craw, *u) * np.array([gs[ch][0] for ch in (0, 1)])[None]
                       + mlines.warp(Rraw, *u) * np.array([gs[ch][1] for ch in (0, 1)])[None], P)
        a0, a1 = cc.K0 + W["cluster_gap"][0], min(cc.K0 + W["first_pass"][1], end[1])
        if a1 - a0 > 200:                            # filtered over the whole capture, then cut (no edge clicks)
            yf, mf = (sosfiltfilt(KF_BAND, x[:, 1]) for x in (yr, mraw))
            mh = sosfiltfilt(KF_HPS, mf)
            out.append((amp, a0, a1) + tuple(x[a0 - K:a1 + K].astype(np.float32) for x in (yf, mf, mh)))
    return key, out


KF_HP = float(os.environ.get("SML_KF_HP", 9000))
KF_ORDER = int(os.environ.get("SML_KF_ORDER", 4))
KF_HPS = butter(KF_ORDER, KF_HP, btype="high", fs=FS, output="sos")    # kfix: the correction acts above KF_HP only
KF_COST = os.environ.get("SML_KF_COST", "hp")                   # kfix: "full" scores the fit on the full-band residual (leakage costs)
KF_HPLEN = int(os.environ.get("SML_KF_FIR", 48))                 # its zero-phase high-pass as an FIR of +-KF_HPLEN taps
KF_BAND = butter(3, (20, 12300), btype="band", fs=FS, output="sos")


def kfix_stage(workers, every):
    """Per program: an 11-tap correction to the current output kernel (as tools/ml_kernel_fix.py), fitted on the R channel's
    cluster / gap and first pass of the undamped groups in the full band 20 Hz-12.3 kHz (captures uncut, model x P(f)),
    held out by LF step; written as a candidate kernel file (ML_KFIX_FILE) -> analysis/fullband/kernel_fix_full.{json,log}."""
    from concurrent.futures import ProcessPoolExecutor
    import ml_carryover as cc
    import ml_lines as mlines
    hp = KF_HPS
    res, lines = {}, []
    for prog in ("M", "L"):
        h = np.array(mlines.Model(prog).h)
        K = (len(h) - 1) // 2
        with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:   # the groups the ml stage does not score (odd ones)
            got = dict(ex.map(_kfix_job, [(prog, k, K) for k in sorted(cc.groups(prog))[1::every]]))
        caps = [(k, a0, a1, y.astype(float), m.astype(float), mh.astype(float)) for k, v in got.items() for _, a0, a1, y, m, mh in v]
        seg = lambda m: [np.roll(m, -k)[K:-K] for k in range(-K, K + 1)]      # m(n + k) on the window

        def design(sel):
            A, b = [], []
            for (k, a0, a1, y, m, mh), s_ in zip(caps, sel):
                if s_:
                    yy, mm = y[K:-K], m[K:-K]
                    gg = np.dot(yy, mm) / np.dot(mm, mm)
                    A.append(np.stack(seg(mh), 1))
                    b.append((y / gg - m)[K:-K] if KF_COST == "full" else sosfiltfilt(hp, y / gg - m)[K:-K])
            return np.vstack(A), np.concatenate(b)

        def score(c, sel):
            rows = []
            for (k, a0, a1, y, m, mh), s_ in zip(caps, sel):
                if s_:
                    yy = y[K:-K]
                    mm = m[K:-K] + np.stack(seg(mh), 1) @ c
                    rows.append({b_: _seg_null(yy, mm, sos) for b_, sos in zip(BNAME, SOS)})
            return {b_: float(np.median([r[b_] for r in rows])) for b_ in BNAME}
        hi = np.array([float(k.split("_")[0]) > 1.0 for k, *_ in caps])
        for tr, te in ((~hi, hi), (hi, ~hi)):
            A, b = design(tr)
            c = np.linalg.lstsq(A, b, rcond=None)[0]
            s0, s1 = score(np.zeros(2 * K + 1), te), score(c, te)
            lines.append(f"{prog}: held out by LF ({tr.sum()} -> {te.sum()} captures), R cluster / gap + first pass, per band: "
                         + "; ".join(f"{b_} {s0[b_]:.1f} -> {s1[b_]:.1f}" for b_ in BNAME))
        A, b = design(np.ones(len(caps), bool))
        c = np.linalg.lstsq(A, b, rcond=None)[0]
        d_ = np.zeros(2 * KF_HPLEN + 1)
        d_[KF_HPLEN] = 1.0
        hpf = sosfiltfilt(KF_HPS, np.r_[np.zeros(512), d_, np.zeros(512)])[512:512 + len(d_)]   # the zero-phase high-pass as an FIR
        corr = np.convolve(c[::-1], hpf)                 # the fit's regressors are m(n + k): as a convolution kernel c reversed
        corr[KF_HPLEN + K] += 1.0
        hn = np.convolve(h, corr)                       # centred: (len(h) + len(corr) - 1) taps
        res[prog] = dict(kernel=hn.tolist(), base=h.tolist(), correction=c.tolist(), hp=dict(order=KF_ORDER, hz=KF_HP, fir_half=KF_HPLEN), captures=len(caps),
                         numbers=2 * K + 1)
        lines.append(f"{prog}: {len(caps)} captures; correction above {KF_HP:.0f} Hz {np.round(c, 4).tolist()}")
        print(lines[-3], lines[-2], lines[-1], sep="\n", flush=True)
    tag = ("" if KF_COST == "hp" else "_" + KF_COST) + ("" if (KF_ORDER, KF_HP) == (4, 9000.0) else f"_o{KF_ORDER}_{KF_HP:.0f}")
    json.dump(res, open(os.path.join(OUT, f"kernel_fix_full{tag}.json"), "w"), indent=1)
    open(os.path.join(OUT, f"kernel_fix_full{tag}.log"), "w").write("\n".join(lines) + "\n")


def _seg_null(y, m, sos):
    yf, mf = sosfiltfilt(sos, y), sosfiltfilt(sos, m)
    gg = 1.0
    return float(10 * np.log10(np.sum((yf - gg * mf) ** 2) / np.sum(yf ** 2)))


# ---------- stage warp ----------
def ksinc(x, t, half, beta, block=2048):
    """x (n, c) read at the positions t, Kaiser-windowed sinc of half-length half (q1.warp's kernel)."""
    i = np.floor(t).astype(int)
    fr = t - i
    k = np.arange(-half + 1, half + 1)
    xp = np.concatenate([np.zeros((half,) + x.shape[1:]), x, np.zeros((half + 2,) + x.shape[1:])])
    out = np.zeros((len(t),) + x.shape[1:])
    for s0 in range(0, len(t), block):
        dd = fr[s0:s0 + block, None] - k[None, :]
        wgt = np.sinc(dd) * np.i0(beta * np.sqrt(np.clip(1 - (dd / (half + 0.5)) ** 2, 0, None))) / np.i0(beta)
        idx = np.clip(i[s0:s0 + block, None] + k[None, :] + half, 0, len(xp) - 1)
        out[s0:s0 + block] = np.einsum("nk,nk...->n...", wgt, xp[idx])
    return out


WARP_REF = (512, 14.0)


def warp_stage():
    """Each warp against the same read-out map through a 512 / beta 14 kernel, per band and window: dB of
    (test - reference) below the reference's band energy. Signals: an S capture as the S tools warp it (1.0_1.4_1.0
    AMPL 90, chain3 cache; 12.3 kHz band window) and an M model render as ml_lines warps it (full band, then x P(f)),
    at the 5th / 50th / 95th percentile of the library's clock ratios."""
    import ml_lines as mlines
    info, fade, _, _, eps0 = q2.q1_tables()
    key = "1.0_1.4_1.0"
    nm = [f for a_, f in q2.files_of(info[key]) if a_ == 90][0]
    n = PRE + int(1.5 * FS)
    xs = q1.load(q2.Q1, nm)[:n]
    Mm = mlines.Model("M")
    mkey = "1.0_1.4_0.7" if "1.0_1.4_0.7" in __import__("ml_carryover").groups("M") else sorted(__import__("ml_carryover").groups("M"))[0]
    xm = np.sum(Mm.parts(n, mkey, band=False), 0)
    N = 1 << int(np.ceil(np.log2(n + 8192)))
    P = pipeline(N, True)
    eps = sorted(eps0.values())
    ppms = [float(np.percentile(eps, 5)), float(np.median(eps)), float(np.percentile(eps, 95))]
    nn = np.arange(n, dtype=float)
    k0 = mlines.K0
    tests = [("q1.warp 128/12 (S apply_clock)", xs, lambda x, p: q1.warp(x, p), lambda p: PRE + (nn - PRE) * (1 + p * 1e-6), False),
             ("q1.warp 32/8 (S clock fit)", xs, lambda x, p: q1.warp(x, p, half=32, beta=8.0), lambda p: PRE + (nn - PRE) * (1 + p * 1e-6), False),
             ("ml_lines.warp 48/10 (M / L model)", xm, lambda x, p: mlines.warp(x, 0.0, p), lambda p: k0 + (nn - k0) / (1 + p * 1e-6), True)]
    bands = list(zip(BNAME, SOS)) + [("12300-12800", butter(3, (12300, 12800), btype="band", fs=FS, output="sos"))]
    out = dict(key=key, file=nm, mkey=mkey, ppm=ppms, ref=WARP_REF, res={})
    lines = [f"warp artefact, dB below the band (S capture {nm}; M model {mkey}); windows 50-150 / 300-600 / 1000-1500 ms"]
    for p in ppms:
        for name, x, fn, tmap, post in tests:
            T, R = fn(x, p), ksinc(x, tmap(p), *WARP_REF)
            if post:
                T, R = apply_f(T, P), apply_f(R, P)
            row = {}
            for b, sos in bands:
                fb, fd = sosfiltfilt(sos, R, axis=0), sosfiltfilt(sos, T - R, axis=0)
                row[b] = {}
                for a_ms, b_ms in WINS + [(1000, 1500)]:
                    ix = win_idx(a_ms, b_ms, 1.5, n)
                    if ix is None:
                        continue
                    i0, i1 = ix
                    row[b][f"{a_ms}-{b_ms}"] = float(10 * np.log10(max(np.sum(fd[i0:i1] ** 2), 1e-30) / max(np.sum(fb[i0:i1] ** 2), 1e-30)))
            out["res"][f"{name} @ {p:.1f} ppm"] = row
            lines.append(f"{name} @ {p:+.1f} ppm: " + "; ".join(
                f"{b} " + " / ".join(f"{row[b].get(w, float('nan')):.0f}" for w in ("50-150", "300-600", "1000-1500"))
                for b in ("20-40", "3200-6400", "10000-11000", "11000-12300", "12300-12800")))
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, "warp.json"), "w"), indent=1)
    open(os.path.join(OUT, "warp.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------- stage wide: 12.3-12.8 kHz ----------
WIDE_WIN = (12000.0, 12800.0)                    # N's wide analysis band (NonLin_program 11.3); the unit is -43 dB at 12.8 kHz
WBANDS = [(6400, 10000), (10000, 11000), (11000, 12000), (12000, 12300), (12300, 12800)]
WBNAME = [f"{a}-{b}" for a, b in WBANDS]
WSOS = [butter(3, ab, btype="band", fs=FS, output="sos") for ab in WBANDS]


def _with_window(lohi, fn, *a):
    import core
    old = core.F_ROLL, core.F_HI
    core.F_ROLL, core.F_HI = lohi
    try:
        return fn(*a)
    finally:
        core.F_ROLL, core.F_HI = old


def device_ir(path):
    """ml_structure._make in memory (the same deconvolution, anchor and resampling), with the current band window."""
    import soundfile as sf
    from core import PAD, Kernel, deconvolve
    ms = q1.ms
    x, fs = sf.read(path, always_2d=True)
    a = np.abs(x).max(1)
    main = int(np.argmax(a > a.max() * 0.1))
    lo = max(main - int(0.002 * fs), 0)
    onset = lo + int(np.argmax(a[lo:main + 1] > a.max() * 1e-2))
    d = deconvolve(np.vstack([np.zeros((PAD, 2)), x]), fs, Kernel.load(ms.KERNEL))
    anchor = float(onset + PAD)
    k0 = int(np.floor(anchor / fs * ms.FSD))
    h = np.stack([to_device_grid(d[:, 0], fs, anchor + ms.LR_US * 1e-6 * fs, k0, ms.CFG),
                  to_device_grid(d[:, 1], fs, anchor, k0, ms.CFG)], 1)
    return h[max(k0 - ms.K0, 0):]


def q1_align(h, amp, shift_lr):
    """q1.data_stage's per-channel alignment from its stored shifts (AMPL 90: shifted in place; others: PRE zeros first)."""
    if amp != 90:
        h = np.vstack([np.zeros((PRE, 2)), h])
        shift_lr = [s_ - PRE for s_ in shift_lr]
    return np.stack([q1.fshift(h[:, c], shift_lr[c]) for c in (0, 1)], 1)


def regularised(N):
    """The deconvolution's own response to a flat input on an N-point device-rate grid (2 channels): |G|^2 / (|G|^2 +
    eps max|G|^2) of the chain kernel (nl_structure.processing_response without the window), unit at 1 kHz."""
    from core import Kernel
    k = Kernel.load(q1.ms.KERNEL)
    M = 1 << 18
    ff = np.fft.rfftfreq(M, 1 / k.fs)
    G = k.response(M, k.fs)
    G = G if G.ndim == 2 else np.stack([G, G], 1)
    g2 = np.abs(G) ** 2
    R = g2 / (g2 + 1e-3 * np.max(np.abs(G), 0) ** 2)
    R = R / R[np.argmin(np.abs(ff - 1000))]
    f = np.fft.rfftfreq(N, 1 / FS)
    return np.stack([np.interp(f, ff, R[:, c]) for c in (0, 1)], 1)


def _wide_job(a):
    key, lab, b230, law, fold, fc, v, caps, ends, path = a
    from scipy.signal import sosfiltfilt as sff
    q1.EXT = 1700
    G = pc.setting_G(lab, b230, law, fold, fc)
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
    cl, f_, xs = J.split(v)
    P = {}
    chk = []
    for c in caps:
        yb = q2._band(c["file"])
        n = min(len(yb), G["N"])
        yb = yb[:n]
        mdl = J.render(v, G, n, c["amp"])
        if c["fade2"] >= 0.15:
            _, (e, sh, gain) = q1._clock_job((yb, mdl, c["fade2"], c["eps"], c["amp"] == 0))
        else:
            g0 = q2.ramp_gain(n, ends[c["file"]])
            _, (e, sh, gain) = q1._clock_job((yb, mdl * g0, q2.comp_end(yb, g0), c["eps"], c["amp"] == 0))
        y = q1.apply_clock(yb, e, sh)[:n] / gain
        g = q2.ramp_gain(n, [PRE + (x + sh - PRE) / (1 + e * 1e-6) for x in ends[c["file"]]])
        te = q2.comp_end(y, g)
        hw = q1_align(_with_window(WIDE_WIN, device_ir, c["path"]), c["amp"], c["shift"])
        if not chk:                                  # the same path at the harness's window reproduces q1's cache
            hn = q1_align(device_ir(c["path"]), c["amp"], c["shift"])
            ref = q1.load(q2.Q1, c["file"])
            m_ = min(len(hn), len(ref))
            chk.append(float(10 * np.log10(np.sum((hn[:m_] - ref[:m_]) ** 2) / np.sum(ref[:m_] ** 2))))
        hn_ = q1.load(q2.Q1, c["file"])                # the wide deconvolution's peak normalisation differs (a wider band
        m_ = min(len(hn_), len(hw))                    # window, a taller impulse): one scalar to the harness cache in band
        bw, bn = (sff(q1.ms.BAND, z[:m_], axis=0) for z in (hw, hn_))
        hw = hw * float(np.sum(bw * bn) / np.sum(bw * bw))
        yw = q1.apply_clock(sff(q1.LFBAND, hw[:n] if len(hw) >= n else np.vstack([hw, np.zeros((n - len(hw), 2))]), axis=0), e, sh)[:n] / gain
        r = J.reverb(G, *J._spectra(f_, xs, G["N"]))[:n]
        cd = np.zeros((n, 2))
        for j in range(8):
            cd[PRE + q1.m.CLUSTER[j // 4][j % 4], j // 4] += (c["amp"] / 90) * J.cs * cl[j]
        cd = sff(q1.LFBAND, cd, axis=0)
        N = 1 << int(np.ceil(np.log2(n + 8192)))
        Rg = regularised(N)
        Pw = _with_window(WIDE_WIN, pipeline, N, True)[:, None] * Rg
        Pn = pipeline(N, True)[:, None] * Rg
        mw = apply_f(cd + r, Pw)
        mn = apply_f(cd + r, Pn)
        P[c["amp"]] = dict(y=yw, g=g, te=te, e=yw - mw * g, en=y - mn * g, yn=y)
    out = dict(key=key, **lab, cache_check_db=chk, caps={}, pairs={})
    for amp, d in P.items():
        sig = dict(y=d["y"], e=d["e"], yn=d["yn"], en=d["en"])
        out["caps"][amp] = dict(te=d["te"], bands=_wbands(sig, d["te"]))
    for a_, b_ in CROSS:
        if a_ in P and b_ in P:
            A, B = P[a_], P[b_]
            n = min(len(A["y"]), len(B["y"]))
            sig = dict(ya=A["y"][:n] * B["g"][:n], d=A["y"][:n] * B["g"][:n] - B["y"][:n] * A["g"][:n])
            out["pairs"][f"{a_}-{b_}"] = dict(bands=_wbands(sig, min(A["te"], B["te"])))
    json.dump(out, open(path, "w"))
    print(key, chk, flush=True)
    return key


def _wbands(sigs, te):
    n = min(len(x) for x in sigs.values())
    out = {}
    for b, sos in zip(WBNAME + ["1000-4000"], WSOS + [butter(3, (1000, 4000), btype="band", fs=FS, output="sos")]):
        f = {k: sosfiltfilt(sos, x[:n], axis=0) for k, x in sigs.items()}
        out[b] = {}
        for a_ms, b_ms in WINS[:4]:
            ix = win_idx(a_ms, b_ms, te, n)
            if ix is not None:
                out[b][f"{a_ms}-{b_ms}"] = {k: np.sum(x[ix[0]:ix[1]] ** 2, 0).tolist() for k, x in f.items()}
    return out


def wide_stage(workers, every, only=None):
    from concurrent.futures import ProcessPoolExecutor
    d = os.path.join(OUT, "wide")
    os.makedirs(d, exist_ok=True)
    info, fade, _, _, eps0 = q2.q1_tables()
    paths = {os.path.basename(p): p for p in glob.glob(os.path.join(ROOT, "irs", "EMT 246 captures", "96kHz", "Short RT", "*.wav"))}
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    folds = {rt: json.load(open(os.path.join(pc.CV, f"fold_{rt}.json"))) for rt in sorted({g["rt"] for g in info.values()})}
    jobs = []
    for k in sorted(info)[::every]:
        g = info[k]
        path = os.path.join(d, k + ".json")
        if (only and k not in only) or os.path.exists(path):
            continue
        lab = dict(lf=g["lf"], rt=g["rt"], hf=g["hf"])
        caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm], path=paths[nm], shift=g["files"][str(a_)]["shift_LR"])
                for a_, nm in q2.files_of(g) if (k, a_) not in DEFECTS]
        if caps:
            jobs.append((k, lab, ap["steps"][ap["settings"][k]["step"]], law, folds[g["rt"]], cells.get(k), v, caps, ends, path))
    print(f"{len(jobs)} settings to run", flush=True)
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_wide_job, jobs):
            pass


# ---------- stage xfit: a per-loop correction of S's loop paths ----------
X_LAGS = list(range(-10, 11)) + [s_ * (e + d) for e in (230, 460) for d in range(-2, 3) for s_ in (1, -1)]
X_WINS = [(0, 50), (50, 150), (150, 300)]           # fitted where the first pass and the first recirculation lie, each equally


def _xfit_job(a):
    """One S setting on the baseline (EMT_XCORR unset): per capture the residual and each loop path's contribution r_ci
    (both x the pipeline response), shifted by every lag in X_LAGS, over X_WINS (rows scaled per window by the capture's
    energy there): the normal equations A'A, A'b of the per-loop correction, shared by L and R."""
    key, lab, b230, law, fold, fc, v, caps, ends = a
    q1.EXT = 1700
    G = pc.setting_G(lab, b230, law, fold, fc)
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((4096, 2)), w=np.ones(4096), amp=0)])])
    cl, f_, xs = J.split(v)
    F_, X_ = J._spectra(f_, xs, G["N"])
    nl, ncol = len(G["R"]), len(G["R"]) * len(X_LAGS)
    ATA, ATb = np.zeros((ncol, ncol)), np.zeros(ncol)
    for c in caps:
        y, mdl, g, te = q2.fit_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)
        n = len(y)
        mf = full_model(J, v, G, n, c["amp"], mdl)
        e = y - mf * g
        N = 1 << int(np.ceil(np.log2(n + 8192)))
        Pr = pipeline(N, True)
        for ch in (0, 1):
            ri = [apply_f(np.fft.irfft(X_[ch, i] * G["R"][i], G["N"])[:n, None] * np.ones((1, 2)), Pr)[:, 0] for i in range(nl)]
            for a_ms, b_ms in X_WINS:
                ix = win_idx(a_ms, b_ms, te, n)
                if ix is None:
                    continue
                i0, i1 = ix
                sc = 1.0 / np.sqrt(np.sum(y[i0:i1, ch] ** 2))
                A = np.stack([np.roll(r, lag)[i0:i1] * g[i0:i1, ch] for r in ri for lag in X_LAGS], 1) * sc
                bb = e[i0:i1, ch] * sc
                ATA += A.T @ A
                ATb += A.T @ bb
    print(key, flush=True)
    return key, lab["rt"], ATA, ATb


def xfit_stage(workers, only=None):
    from concurrent.futures import ProcessPoolExecutor
    assert not os.environ.get("EMT_XCORR"), "fit on the baseline: unset EMT_XCORR"
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    cells = json.load(open(os.path.join(pc.CV, "cells.json")))
    folds = {rt: json.load(open(os.path.join(pc.CV, f"fold_{rt}.json"))) for rt in sorted({g["rt"] for g in info.values()})}
    jobs = []
    for k, g in sorted(info.items()):
        caps = [dict(amp=a_, file=nm, eps=eps0[nm], fade2=fade[nm]) for a_, nm in q2.files_of(g) if (k, a_) not in DEFECTS]
        if caps and (not only or k in only):
            jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], law, folds[g["rt"]], cells.get(k), v, caps, ends))
    acc = {}
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for key, rt, ATA, ATb in ex.map(_xfit_job, jobs):
            for name, rts in list(TFOLDS.items()) + [("all", tuple(sorted({g["rt"] for g in info.values()})))]:
                if rt in rts:
                    A0, b0 = acc.get(name, (0, 0))
                    acc[name] = (A0 + ATA, b0 + ATb)
    nl = 4
    out = dict(lags=X_LAGS, windows=X_WINS, numbers=nl * len(X_LAGS), folds={}, all=None)
    for name, (A, bb) in acc.items():
        h = np.linalg.solve(A + 1e-9 * np.trace(A) / len(A) * np.eye(len(A)), bb).reshape(nl, len(X_LAGS))
        entry = dict(train_rt=list(TFOLDS[name]) if name in TFOLDS else "all", h=h.tolist())
        if name == "all":
            out["all"] = entry
        else:
            out["folds"][name] = entry
        print(name, "per loop |h| at lag 0 / sum |h|:", [f"{x[X_LAGS.index(0)]:+.5f} / {np.abs(x).sum():.4f}" for x in h], flush=True)
    hA, hB = np.array(out["folds"]["A"]["h"]), np.array(out["folds"]["B"]["h"])
    out["fold_agreement"] = float(np.sum(hA * hB) / np.sqrt(np.sum(hA ** 2) * np.sum(hB ** 2)))
    print("fold agreement (cosine)", round(out["fold_agreement"], 3), flush=True)
    json.dump(out, open(os.path.join(OUT, "xcorr" + ("_smoke" if only else "") + ".json"), "w"), indent=1)


# ---------- stage gate: the full-band sanity gate ----------
GATE_TOL = 0.5                                       # dB a band-window median may get worse


def gate_stage(workers, tag):
    """pulse_check's gate keys (16 S settings, 8 M / L groups) scored in BANDS: S nulls per band at 50-150 / 150-300 /
    300-600 ms, M / L per band in cluster / gap and first pass, L and R. --tag ref writes the reference (the adopted model);
    any other tag is PASS / FAIL against it: no band-window median worse by more than GATE_TOL. Run beside
    pulse_check --stage gate after any model change (its nulls are 100 Hz-10 kHz only). -> analysis/fullband/gate_<tag>/"""
    global OUT
    from concurrent.futures import ProcessPoolExecutor
    base = OUT
    OUT = os.path.join(base, f"gate_{tag}")
    try:
        s_stage(workers, list(pc.GATE_S), sub="s")
        d = os.path.join(OUT, "ml")
        os.makedirs(d, exist_ok=True)
        jobs = [(p, k, os.path.join(d, f"{p}_{k}.json")) for p in "ML" for k in pc.GATE_ML[p] if not os.path.exists(os.path.join(d, f"{p}_{k}.json"))]
        with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
            for _ in ex.map(_ml_job, jobs):
                pass
    finally:
        OUT = base
    if tag == "ref":
        print("full-band gate reference written", flush=True)
        return
    G, R = os.path.join(base, f"gate_{tag}"), os.path.join(base, "gate_ref")
    lines, ok = [f"## Full-band gate '{tag}' against gate_ref (BANDS 20 Hz-12.3 kHz; median null per band and window, L / R; allowed +{GATE_TOL} dB)"], True
    for w in ("50-150", "150-300", "300-600"):
        for b in BNAME:
            for ch in (0, 1):
                val = lambda root: _med([_db(c["bands"][b][w]["e"][ch], c["bands"][b][w]["y"][ch]) for k in pc.GATE_S
                                         for c in json.load(open(os.path.join(root, "s", k + ".json")))["caps"].values() if w in c["bands"][b]])
                v, r = val(G), val(R)
                if v is not None and r is not None and v > r + GATE_TOL:
                    ok = False
                    lines.append(f"- FAIL S {w} ms {b} {'LR'[ch]}: {v:.2f} (reference {r:.2f})")
    for p in "ML":
        for w in ("cluster_gap", "first_pass"):
            for b in BNAME:
                for ch in (0, 1):
                    def val(root):
                        xs = [r_["bands"][b][w][ch] for k in pc.GATE_ML[p] for r_ in json.load(open(os.path.join(root, "ml", f"{p}_{k}.json")))
                              if b in r_["bands"] and w in r_["bands"][b] and r_["bands"][b][w][ch]]
                        return _med([_db(x[1], x[0]) for x in xs])
                    v, r = val(G), val(R)
                    if v is not None and r is not None and v > r + GATE_TOL:
                        ok = False
                        lines.append(f"- FAIL {p} {w} {b} {'LR'[ch]}: {v:.2f} (reference {r:.2f})")
    lines.append(f"**{'PASS' if ok else 'FAIL'}** ({'no band-window median worse by more than ' + str(GATE_TOL) + ' dB' if ok else 'see above'})")
    open(os.path.join(G, "gate.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------- report ----------
def _db(num, den):
    return float(10 * np.log10(num / den)) if num > 0 and den > 0 else None


def _med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(xs)) if xs else None


def _fmt(x, w=6):
    return f"{x:{w}.1f}" if x is not None else " " * (w - 1) + "-"


def report_stage():
    lines, rep = [], {}
    # S: strict nulls per band and window, cross-take repeat, excess, take-to-take correlation of the residual
    _report_s("s", lines, rep)
    if glob.glob(os.path.join(OUT, "s", "*.npz")):
        fr = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 2000, 4000, 8000, 10000, 11000, 12000]
        TA, TB = fold_T("A"), fold_T("B")
        rep["T"] = {}
        lines.append("\n## The fixed transfer capture / model over the reverb (150 ms on), LF >= 1.0, fold A ; fold B, dB / deg")
        for c in (0, 1):
            cells = []
            for x in fr:
                k = int(x / TBIN)
                a_, b_ = TA[c, k], TB[c, k]
                cells.append(f"{x:g}: {20 * np.log10(abs(a_)):+.2f} / {np.degrees(np.angle(a_)):+.1f} ; {20 * np.log10(abs(b_)):+.2f} / {np.degrees(np.angle(b_)):+.1f}")
                rep["T"].setdefault("LR"[c], {})[str(x)] = [20 * np.log10(abs(a_)), np.degrees(np.angle(a_)), 20 * np.log10(abs(b_)), np.degrees(np.angle(b_))]
            lines.append(f"- {'LR'[c]}: " + "; ".join(cells))
    for vdir in sorted(os.path.basename(x) for x in glob.glob(os.path.join(OUT, "s_tfix*"))):
        _report_variant("s", vdir, lines, rep)
    _report_rest(lines, rep)


def _report_variant(base, vdir, lines, rep):
    """A held-out variant against the baseline, capture by capture: median null per band and window (R / L) and the
    captures more than 0.5 dB worse (the strict rule, per band)."""
    S = {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(os.path.join(OUT, base, "*.json"))}
    V = {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(os.path.join(OUT, vdir, "*.json"))}
    keys = sorted(set(S) & set(V))
    if not keys:
        return
    rep["vs_" + vdir] = {}
    lines.append(f"\n## {vdir} against {base}: {len(keys)} settings; median null R base -> variant [L], captures worse by > 0.5 dB (R+L)")
    lines.append("| band | " + " | ".join(f"{a}-{b}" for a, b in WINS) + " |")
    lines.append("|---" * (len(WINS) + 1) + "|")
    for b in BNAME:
        cells = []
        for a_, b_ in WINS:
            wn = f"{a_}-{b_}"
            x0, x1, worse, n = {0: [], 1: []}, {0: [], 1: []}, 0, 0
            for k in keys:
                for amp, c in S[k]["caps"].items():
                    cv = V[k]["caps"].get(amp)
                    if cv is None or wn not in c["bands"][b] or wn not in cv["bands"][b]:
                        continue
                    for ch in (0, 1):
                        u = _db(c["bands"][b][wn]["e"][ch], c["bands"][b][wn]["y"][ch])
                        w = _db(cv["bands"][b][wn]["e"][ch], cv["bands"][b][wn]["y"][ch])
                        x0[ch].append(u); x1[ch].append(w)
                        worse += w > u + 0.5
                        n += 1
            rep["vs_" + vdir].setdefault(wn, {})[b] = [_med(x0[1]), _med(x1[1]), _med(x0[0]), _med(x1[0]), worse, n]
            cells.append(f"{_fmt(_med(x0[1])).strip()} -> {_fmt(_med(x1[1])).strip()} [{_fmt(_med(x1[0])).strip()}] {worse}/{n}")
        lines.append(f"| {b} | " + " | ".join(cells) + " |")


def _report_s(sdir, lines, rep):
    S = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(OUT, sdir, "*.json")))]
    if S:
        rep[sdir] = {}
        lines.append(f"\n## {sdir}: S, {len(S)} settings, clean captures: median strict null [cross-take repeat] (excess), residual correlation A-B; L / R")
        for wn in [f"{a}-{b}" for a, b in WINS]:
            lines.append(f"\n### {wn} ms")
            lines.append("| band | model L / R | repeat L / R | excess L / R | rho L / R |")
            lines.append("|---|---|---|---|---|")
            for b in BNAME:
                cell = {}
                for ch in (0, 1):
                    mod = [_db(c["bands"][b][wn]["e"][ch], c["bands"][b][wn]["y"][ch]) for s_ in S for c in s_["caps"].values() if wn in c["bands"][b]]
                    rp = [_db(q["bands"][b][wn]["d"][ch], q["bands"][b][wn]["ya"][ch]) for s_ in S for q in s_["pairs"].values() if wn in q["bands"][b]]
                    ex = []
                    for s_ in S:                                   # per setting: its model median minus its repeat median
                        mm = _med([_db(c["bands"][b][wn]["e"][ch], c["bands"][b][wn]["y"][ch]) for c in s_["caps"].values() if wn in c["bands"][b]])
                        rr = _med([_db(q["bands"][b][wn]["d"][ch], q["bands"][b][wn]["ya"][ch]) for q in s_["pairs"].values() if wn in q["bands"][b]])
                        if mm is not None and rr is not None:
                            ex.append(mm - rr)
                    rho = [q["bands"][b][wn]["cab"][ch] / np.sqrt(q["bands"][b][wn]["ea"][ch] * q["bands"][b][wn]["eb"][ch])
                           for s_ in S for q in s_["pairs"].values() if wn in q["bands"][b]]
                    cell[ch] = (_med(mod), _med(rp), _med(ex), _med(rho))
                rep[sdir].setdefault(wn, {})[b] = {k: [cell[0][i], cell[1][i]] for i, k in enumerate(("model", "repeat", "excess", "rho"))}
                lines.append(f"| {b} | {_fmt(cell[0][0])} / {_fmt(cell[1][0])} | {_fmt(cell[0][1])} / {_fmt(cell[1][1])} | "
                             f"{_fmt(cell[0][2])} / {_fmt(cell[1][2])} | " + " / ".join(f"{x:.2f}" if x is not None else "-" for x in (cell[0][3], cell[1][3])) + " |")
        # the sixth-octave profile: excess (model over repeat) and rho, R, per window
        lines.append(f"\n## {sdir}: S, sixth-octave profile (Hann-windowed band energies), R: median excess over the cross-take repeat dB / residual rho")
        fc = np.sqrt(FINE[:-1] * FINE[1:])
        rep[sdir + "_fine"] = {"fc": fc.tolist()}
        hdr = "| f (Hz) | " + " | ".join(f"{a}-{b}" for a, b in FWINS) + " |"
        lines += [hdr, "|---" * (len(FWINS) + 1) + "|"]
        prof = {}
        for wn in [f"{a}-{b}" for a, b in FWINS]:
            for ch in (0, 1):
                ex, rh, mo = [], [], []
                for k in range(len(fc)):
                    e_k, r_k, m_k = [], [], []
                    for s_ in S:
                        mm = _med([_db(c["fine"][wn]["e"][ch][k], c["fine"][wn]["y"][ch][k]) for c in s_["caps"].values() if wn in c["fine"]])
                        rr = _med([_db(q["fine"][wn]["d"][ch][k], q["fine"][wn]["ya"][ch][k]) for q in s_["pairs"].values() if wn in q["fine"]])
                        if mm is not None and rr is not None:
                            e_k.append(mm - rr)
                            m_k.append(mm)
                        r_k += [q["fine"][wn]["cab"][ch][k] / np.sqrt(q["fine"][wn]["ea"][ch][k] * q["fine"][wn]["eb"][ch][k])
                                for q in s_["pairs"].values() if wn in q["fine"] and q["fine"][wn]["ea"][ch][k] > 0 and q["fine"][wn]["eb"][ch][k] > 0]
                    ex.append(_med(e_k)); rh.append(_med(r_k)); mo.append(_med(m_k))
                prof[(wn, ch)] = (ex, rh, mo)
                rep[sdir + "_fine"].setdefault(wn, {})["LR"[ch]] = dict(excess=ex, rho=rh, model=mo)
        for k in range(len(fc)):
            lines.append(f"| {fc[k]:.0f} | " + " | ".join(
                f"{_fmt(prof[(f'{a}-{b}', 1)][0][k])} / {prof[(f'{a}-{b}', 1)][1][k]:.2f}" if prof[(f'{a}-{b}', 1)][1][k] is not None
                else f"{_fmt(prof[(f'{a}-{b}', 1)][0][k])} / -" for a, b in FWINS) + " |")


def _report_rest(lines, rep):
    # takefit: base vs refit per band, held out by take
    for td in sorted(glob.glob(os.path.join(OUT, "takefit*"))):
        T = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(td, "*.json")))]
        if not T:
            continue
        name = os.path.basename(td)
        rep[name] = {}
        chk = [c["check"]["150-300"][1] for t in T for c in t["caps"].values() if c["check"].get("150-300")]
        lines.append(f"\n## {name}: {len(T)} settings, held out by take, R median null base -> refit (lead dB); "
                     f"check: base 100 Hz-10 kHz 150-300 ms median {_med(chk):.2f} dB")
        lines.append("| band | " + " | ".join(f"{a}-{b}" for a, b in WINS[1:4]) + " |")
        lines.append("|---" * 4 + "|")
        for b in BNAME:
            row = []
            for a_, b_ in WINS[1:4]:
                wn = f"{a_}-{b_}"
                m0 = [_db(c["bands"][b][wn]["e0"][1], c["bands"][b][wn]["y"][1]) for t in T for c in t["caps"].values() if wn in c["bands"][b]]
                m1 = [_db(c["bands"][b][wn]["e1"][1], c["bands"][b][wn]["y"][1]) for t in T for c in t["caps"].values() if wn in c["bands"][b]]
                lead = [x - y for x, y in zip(m0, m1) if x is not None and y is not None]
                rep[name].setdefault(wn, {})[b] = [_med(m0), _med(m1), _med(lead)]
                row.append(f"{_fmt(_med(m0))} -> {_fmt(_med(m1))} ({_fmt(_med(lead), 4)})")
            lines.append(f"| {b} | " + " | ".join(row) + " |")
        # energy share of the lead: where the refit removes residual energy (summed over captures, 100 Hz-10 kHz = the fit band)
        lines.append(f"- share of the residual energy the refit removes, by band (R, summed over captures):")
        for a_, b_ in WINS[1:4]:
            wn = f"{a_}-{b_}"
            dE = {b: sum(c["bands"][b][wn]["e0"][1] - c["bands"][b][wn]["e1"][1] for t in T for c in t["caps"].values() if wn in c["bands"][b]) for b in BNAME}
            tot = sum(v for b, v in dE.items() if 100 <= int(b.split("-")[0]) and int(b.split("-")[1]) <= 10000 or b in ("80-160",))
            lines.append(f"  - {wn} ms: " + ", ".join(f"{b} {v / max(tot, 1e-30) * 100:+.0f} %" for b, v in dE.items()))
    # M / L (and S) through ml_lines, first pass etc.
    for mdir in ("ml", "ml_kf", "ml_tfix"):
        _report_ml(mdir, lines, rep)
    for vdir in ("ml_kf", "ml_tfix"):
        _report_ml_variant("ml", vdir, lines, rep)
    _report_wide(lines, rep)
    os.makedirs(OUT, exist_ok=True)
    json.dump(rep, open(os.path.join(OUT, "report.json"), "w"), indent=1)
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def _report_wide(lines, rep):
    Wf = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(OUT, "wide", "*.json")))]
    if not Wf:
        return
    chk = [x for w in Wf for x in w["cache_check_db"]]
    lines.append(f"\n## wide: S, {len(Wf)} settings (every 4th), deconvolved with the 12.0-12.8 kHz window (cache check: the harness "
                 f"window through the same path reproduces q1's cache to {max(chk):.0f} dB). Median R [L]: level re 1-4 kHz; model null with "
                 f"the wide band; the same model in the harness's band; cross-take repeat (wide)")
    rep["wide"] = {}
    lines.append("| band | window | level | model wide | model harness band | repeat |")
    lines.append("|---|---|---|---|---|---|")
    for b in WBNAME:
        for w in ("50-150", "150-300", "300-600"):
            row = []
            for ch in (1, 0):
                lv = [_db(c["bands"][b][w]["y"][ch], c["bands"]["1000-4000"][w]["y"][ch]) for x in Wf for c in x["caps"].values() if w in c["bands"][b]]
                mw = [_db(c["bands"][b][w]["e"][ch], c["bands"][b][w]["y"][ch]) for x in Wf for c in x["caps"].values() if w in c["bands"][b]]
                mn = [_db(c["bands"][b][w]["en"][ch], c["bands"][b][w]["yn"][ch]) for x in Wf for c in x["caps"].values() if w in c["bands"][b]]
                rp = [_db(q["bands"][b][w]["d"][ch], q["bands"][b][w]["ya"][ch]) for x in Wf for q in x["pairs"].values() if w in q["bands"][b]]
                row.append([_med(lv), _med(mw), _med(mn), _med(rp)])
            rep["wide"].setdefault(b, {})[w] = row
            lines.append(f"| {b} | {w} | " + " | ".join(f"{_fmt(row[0][i])} [{_fmt(row[1][i]).strip()}]" for i in range(4)) + " |")


def _report_ml_variant(base, vdir, lines, rep):
    """M / L: a variant's per-capture change against the base run (same groups and captures): median null base -> variant,
    R [L], and the channel-nulls more than 0.5 dB worse."""
    for prog in "ML":
        A = {os.path.basename(p): json.load(open(p)) for p in glob.glob(os.path.join(OUT, base, f"{prog}_*.json"))}
        B = {os.path.basename(p): json.load(open(p)) for p in glob.glob(os.path.join(OUT, vdir, f"{prog}_*.json"))}
        keys = sorted(set(A) & set(B))
        if not keys:
            continue
        lines.append(f"\n## {vdir} against {base}, {prog}: {len(keys)} groups; median null R base -> variant [L base -> variant], channel-nulls worse by > 0.5 dB")
        lines.append("| band | " + " | ".join(ML_WINS) + " |")
        lines.append("|---" * (len(ML_WINS) + 1) + "|")
        for b in BNAME:
            cells = []
            for w in ML_WINS:
                x0, x1, worse, n = {0: [], 1: []}, {0: [], 1: []}, 0, 0
                for k in keys:
                    for ra, rb in zip(A[k], B[k]):
                        for ch in (0, 1):
                            ca = ra["bands"].get(b, {}).get(w, [None, None])[ch]
                            cb = rb["bands"].get(b, {}).get(w, [None, None])[ch]
                            if not ca or not cb:
                                continue
                            u, v_ = _db(ca[1], ca[0]), _db(cb[1], cb[0])
                            x0[ch].append(u); x1[ch].append(v_)
                            worse += v_ > u + 0.5
                            n += 1
                rep.setdefault(f"vs_{vdir}", {}).setdefault(prog, {}).setdefault(w, {})[b] = [_med(x0[1]), _med(x1[1]), _med(x0[0]), _med(x1[0]), worse, n]
                cells.append(f"{_fmt(_med(x0[1])).strip()} -> {_fmt(_med(x1[1])).strip()} [{_fmt(_med(x0[0])).strip()} -> {_fmt(_med(x1[0])).strip()}] {worse}/{n}")
            lines.append(f"| {b} | " + " | ".join(cells) + " |")


def _report_ml(mdir, lines, rep):
    MLf = sorted(glob.glob(os.path.join(OUT, mdir, "*.json")))
    if MLf:
        rows = {}
        for p in MLf:
            prog = os.path.basename(p)[0]
            rows.setdefault(prog, []).extend(json.load(open(p)))
        rep[mdir] = {}
        for prog, R in rows.items():
            lines.append(f"\n## {mdir}: {prog} (ml_lines model, {len(R)} captures): median model null [repeat] (excess per capture), L / R")
            lines.append("| band | " + " | ".join(ML_WINS) + " |")
            lines.append("|---" * (len(ML_WINS) + 1) + "|")
            for b in BNAME:
                cells = []
                for w in ML_WINS:
                    cc_ = []
                    for ch in (0, 1):
                        xs_ = [r_["bands"][b][w][ch] for r_ in R if b in r_["bands"] and w in r_["bands"][b] and r_["bands"][b][w][ch]]
                        mo = [_db(x[1], x[0]) for x in xs_]
                        rp = [_db(x[2], x[0]) for x in xs_]
                        ex = [a_ - b_ for a_, b_ in zip(mo, rp) if a_ is not None and b_ is not None]
                        cc_.append((_med(mo), _med(rp), _med(ex)))
                    rep[mdir].setdefault(prog, {}).setdefault(w, {})[b] = cc_
                    cells.append(f"{_fmt(cc_[0][0])} / {_fmt(cc_[1][0])} [{_fmt(cc_[0][1])} / {_fmt(cc_[1][1])}] ({_fmt(cc_[0][2], 4)} / {_fmt(cc_[1][2], 4)})")
                lines.append(f"| {b} | " + " | ".join(cells) + " |")
            fc = np.sqrt(FINE[:-1] * FINE[1:])
            for w in ("first_pass", "recirc1"):
                prof = []
                for ch in (0, 1):
                    ex = []
                    for k in range(len(fc)):
                        v_ = [_db(r_["fine"][w][ch][1][k], r_["fine"][w][ch][0][k]) - _db(r_["fine"][w][ch][2][k], r_["fine"][w][ch][0][k])
                              for r_ in R if w in r_["fine"] and r_["fine"][w][ch] and min(r_["fine"][w][ch][j][k] for j in range(3)) > 0]
                        ex.append(_med(v_))
                    prof.append(ex)
                rep[mdir].setdefault(prog, {})[w + "_fine"] = dict(fc=fc.tolist(), excess=prof)
                lines.append(f"- {prog} {w}, sixth-octave excess over repeat (dB) L / R: " + ", ".join(
                    f"{fc[k]:.0f} {_fmt(prof[0][k], 1).strip()}/{_fmt(prof[1][k], 1).strip()}" for k in range(len(fc))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["s", "tfix", "takefit", "ml", "warp", "wide", "kfix", "shf", "gate", "xfit", "report"], required=True)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--prog", default="SML")
    ap.add_argument("--every", type=int, default=2)
    ap.add_argument("--tag", default="ref")
    a = ap.parse_args()
    if a.stage in ("s", "tfix"):
        s_stage(a.workers, a.only, a.stage == "tfix")
    elif a.stage == "takefit":
        takefit_stage(a.workers, a.only)
    elif a.stage == "ml":
        ml_stage(list(a.prog), a.every, a.workers, int(a.only[0]) if a.only else None)
    elif a.stage == "xfit":
        xfit_stage(a.workers, a.only)
    elif a.stage == "gate":
        gate_stage(a.workers, a.tag)
    elif a.stage == "shf":
        shf_fit(a.workers, a.every)
    elif a.stage == "kfix":
        kfix_stage(a.workers, a.every)
    elif a.stage == "wide":
        wide_stage(a.workers, a.every, a.only)
    elif a.stage == "warp":
        warp_stage()
    else:
        report_stage()


if __name__ == "__main__":
    main()
