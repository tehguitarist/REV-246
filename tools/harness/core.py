"""Per-file measurements for the EMT 246 library harness (spec stages A1–A4, B1, B4, C2).

Everything here works on one IR file and returns a flat dict (one row of files.csv) plus a
short device-rate signature used later for duplicate / repeat-capture comparison.
"""
import hashlib
import json
import os
import sys
from fractions import Fraction

import numpy as np
import soundfile as sf
from scipy.signal import butter, find_peaks, resample_poly, sosfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.dirname(HERE))
from analyze_ir import edc, rt_from_edc  # noqa: E402

VERSION = "2"
F_ROLL, F_HI = 11000.0, 12300.0     # analysis band edge (device Nyquist is 13.19 kHz)
F_LO = (3.0, 12.0)                  # low band edge, only for kernels with high-pass corners (nothing inverted below 3 Hz)
PAD = 4096                          # zeros prepended so taps at t≈0 don't wrap in the FFT
SIG_S = 0.3                         # signature length (s, device rate)
BANDS = {"lf": (44, 88), "mid": (707, 1414), "hf": (5900, 8300)}
OCTAVES = [63, 125, 250, 500, 1000, 2000, 4000, 8000]


def load_config(path=None):
    with open(path or os.path.join(HERE, "config.json")) as f:
        return json.load(f)


def fs_dev(cfg):
    return cfg["fs_dev_num"] / cfg["fs_dev_den"]


def dev_ratio(fs, cfg):
    """(up, down) with fs * up / down == fs_dev exactly."""
    fr = Fraction(cfg["fs_dev_num"], cfg["fs_dev_den"] * int(round(fs)))
    return fr.numerator, fr.denominator


class Kernel:
    """Impulse response of the analogue output chain (one cluster tap), at a given capture rate.

    g is (n,) for one kernel used on both channels (the library kernel) or (n, 2) per channel (L, R); `lead` samples of
    it come before the tap's time (index 0 otherwise), and `hp_hz` are first-order analogue high-pass corners applied
    on top (analytic, so their long tails are not truncated)."""

    def __init__(self, g, fs, source="", lead=0, hp_hz=()):
        self.g, self.fs, self.source = np.asarray(g, float), float(fs), str(source)
        self.lead, self.hp_hz = int(lead), [float(v) for v in hp_hz]
        self._cache = {}

    @classmethod
    def load(cls, path):
        z = np.load(path if os.path.isabs(path) else os.path.join(ROOT, path))
        return cls(z["g"], float(z["fs"]), str(z["source"]), int(z["lead"]) if "lead" in z else 0,
                   z["hp_hz"] if "hp_hz" in z else ())

    def save(self, path):
        np.savez(path, g=self.g, fs=self.fs, source=self.source, lead=self.lead, hp_hz=np.array(self.hp_hz))

    def at(self, fs):
        fs = int(round(fs))
        if fs not in self._cache:
            fr = Fraction(fs, int(round(self.fs)))
            if fr != 1 and self.lead:
                raise ValueError("a kernel with a lead is only defined at its own rate")
            self._cache[fs] = self.g if fr == 1 else resample_poly(self.g, fr.numerator, fr.denominator, axis=0)
        return self._cache[fs]

    def response(self, N, fs=None):
        """Frequency response on an N-point rfft grid at rate fs: (N/2+1,) or (N/2+1, 2)."""
        fs = self.fs if fs is None else fs
        f = np.fft.rfftfreq(N, 1 / fs)
        G = np.fft.rfft(self.at(fs), N, axis=0)
        H = np.exp(2j * np.pi * f * self.lead / fs)
        for fc in self.hp_hz:
            H = H * (1j * f / fc) / (1 + 1j * f / fc)
        return G * (H[:, None] if G.ndim == 2 else H)


# ---------------------------------------------------------------- signal helpers
def _band_window(f):
    return np.where(f < F_ROLL, 1.0, np.where(f < F_HI, 0.5 * (1 + np.cos(np.pi * (f - F_ROLL) / (F_HI - F_ROLL))), 0.0))


def subsample_peak(d, i, up=32, half=24):
    """Sub-sample position and signed value of the band-limited peak near index i."""
    lo = max(i - half, 0)
    seg = d[lo:i + half + 1]
    u = resample_poly(seg, up, 1)
    c = (i - lo) * up
    j = c - up + int(np.argmax(np.abs(u[c - up:c + up + 1])))
    if 0 < j < len(u) - 1:
        a, b, cc = u[j - 1], u[j], u[j + 1]
        den = a - 2 * b + cc
        delta = 0.5 * (a - cc) / den if den != 0 else 0.0
        val = b - 0.25 * (a - cc) * delta
    else:
        delta, val = 0.0, u[j]
    return lo + (j + delta) / up, float(val)


def deconvolve(x, fs, kernel, eps_rel=1e-3):
    """Remove the analogue chain: each cluster tap becomes a unit band-limited impulse. A per-channel kernel deconvolves
    each channel with its own; a kernel with high-pass corners also gets a low band edge (F_LO, zero phase)."""
    n = len(x)
    N = 1 << int(np.ceil(np.log2(2 * (n + len(kernel.g) + 2 * kernel.lead))))
    f = np.fft.rfftfreq(N, 1 / fs)
    G = kernel.response(N, fs)
    G2 = G if G.ndim == 2 else G[:, None]
    win = _band_window(f)
    if kernel.hp_hz:
        win = win * np.sin(np.pi / 2 * np.clip((f - F_LO[0]) / (F_LO[1] - F_LO[0]), 0, 1)) ** 2
    W = np.conj(G2) / (np.abs(G2) ** 2 + eps_rel * np.max(np.abs(G2), 0) ** 2) * win[:, None]
    d = np.fft.irfft(np.fft.rfft(x, N, axis=0) * (W if x.ndim == 2 else W[:, 0]), N, axis=0)[:n]
    norm = []
    for c in range(G2.shape[1]):
        dg = np.roll(np.fft.irfft(G2[:, c] * W[:, c], N), 64)[:256]
        norm.append(abs(subsample_peak(dg, int(np.argmax(np.abs(dg))))[1]))
    return d / (np.array(norm) if G.ndim == 2 else norm[0])


def to_device_grid(dc, fs, anchor, k0, cfg):
    """Resample one deconvolved channel so device sample k0 lands exactly on capture index `anchor`."""
    up, down = dev_ratio(fs, cfg)
    s = anchor - k0 * fs / fs_dev(cfg)
    N = 1 << int(np.ceil(np.log2(2 * len(dc))))
    fr = np.fft.rfftfreq(N)
    shifted = np.fft.irfft(np.fft.rfft(dc, N) * np.exp(2j * np.pi * fr * s), N)[:len(dc)]
    return resample_poly(shifted, up, down)


def band_rt(s, fs, lo, hi):
    """T30 (T20 fallback) and EDT of one band; band-pass applied backwards in time to avoid filter-ring bias."""
    sos = butter(4, [lo, min(hi, 0.45 * fs)], btype="band", fs=fs, output="sos")
    y = sosfilt(sos, s[::-1])[::-1]
    d, nf = edc(y, fs)
    t30 = rt_from_edc(d, fs, -5, -35)
    return {"T": t30 if t30 is not None else rt_from_edc(d, fs, -5, -25),
            "EDT": rt_from_edc(d, fs, 0, -10), "edc_noise_db": nf}


def echo_density_mix_time(s, fs, start, win=0.02, thresh=0.8):
    """Time (s, from start index) where Abel–Huang normalised echo density first exceeds thresh."""
    w = int(win * fs)
    hw = np.hanning(w); hw /= hw.sum()
    for i in range(start, len(s) - w, w // 4):
        seg = s[i:i + w]
        sd = np.sqrt(np.sum(hw * seg ** 2))
        if sd > 0 and np.sum(hw * (np.abs(seg) > sd)) / 0.3173 > thresh:
            return (i - start) / fs + win / 2
    return None


def find_loops(h, start, end, rt, fsd, n_sigma=8.0, max_lag=20000, tol=1):
    """Recirculation periods (device samples) common to both channels of the decay-flattened tail."""
    per = []
    t = np.arange(len(h)) / fsd
    max_lag = min(max_lag, (end - start) // 3)
    if max_lag < 600:
        return [], {}
    acfs = {}
    for c in range(h.shape[1]):
        a = (h[:, c] * 10 ** (3 * t / rt))[start:end]
        a = a - a.mean()
        N = 1 << int(np.ceil(np.log2(2 * len(a))))
        ac = np.fft.irfft(np.abs(np.fft.rfft(a, N)) ** 2, N)[:max_lag]
        ac /= ac[0]
        sd = np.std(ac[100:])
        pk, _ = find_peaks(ac, height=n_sigma * sd)
        fund = []
        for p in sorted(int(q) for q in pk if q > 200):  # ascending, so multiples are recognised
            if not any(round(p / q) >= 2 and abs(p - round(p / q) * q) <= tol * round(p / q) for q in fund):
                fund.append(p)
        per.append(fund)
        acfs[c] = {p: float(ac[p]) for p in fund}
    common = [p for p in per[0] if all(any(abs(p - q) <= tol for q in other) for other in per[1:])]
    return common, acfs


def match_loopset(periods, known, tol=1):
    for name, ref in known.items():
        if len(periods) == len(ref) and all(abs(a - b) <= tol for a, b in zip(sorted(periods), sorted(ref))):
            return name
    return ("U:" + "-".join(str(p) for p in sorted(periods))) if periods else "?"


def _match_template(offsets, signs, tmpl):
    """+1 / -1 if the tap offsets match the template with same / inverted polarity, 0 if mixed signs, None if no match."""
    if len(offsets) != len(tmpl["offsets"]) or any(abs(a - b) > 1 for a, b in zip(offsets, tmpl["offsets"])):
        return None
    prod = [a * b for a, b in zip(signs, tmpl["signs"])]
    return 1 if all(p > 0 for p in prod) else (-1 if all(p < 0 for p in prod) else 0)


def classify_cluster(taps, templates):
    """A4: channel swap / polarity fingerprint from per-channel (offsets, signs)."""
    (oL, sL), (oR, sR) = taps
    TL, TR = templates["L"], templates["R"]
    mL, mR = _match_template(oL, sL, TL), _match_template(oR, sR, TR)
    xL, xR = _match_template(oL, sL, TR), _match_template(oR, sR, TL)
    if mL in (1, -1) and mR in (1, -1):
        if mL == 1 and mR == 1:
            return "ok"
        if mL == -1 and mR == -1:
            return "inverted_both"
        return "inverted_L" if mL == -1 else "inverted_R"
    if xL in (1, -1) and xR in (1, -1):
        extra = "" if (xL == 1 and xR == 1) else "+" + ("inverted_both" if xL == xR == -1 else ("inverted_L" if xL == -1 else "inverted_R"))
        return "swapped" + extra
    if 0 in (mL, mR, xL, xR):
        return "sign_mismatch"
    return "novel_cluster"


# ---------------------------------------------------------------- per-file analysis
def file_key(path, kernel, lr_offset_us):
    st = os.stat(path)
    raw = f"{os.path.abspath(path)}|{st.st_size}|{st.st_mtime_ns}|{VERSION}|{kernel.source}|{len(kernel.g)}|{lr_offset_us}"
    return hashlib.sha1(raw.encode()).hexdigest()[:16]


def analyze_file(path, cfg, kernel, lr_offset_us=None):
    """Return (row: dict, signature: np.ndarray | None)."""
    fsd = fs_dev(cfg)
    lr_offset_us = cfg["default_lr_offset_us"] if lr_offset_us is None else lr_offset_us
    info = sf.info(path)
    x, fs = sf.read(path, always_2d=True)
    row = {"path": path, "fs": fs, "channels": x.shape[1], "subtype": info.subtype,
           "duration_s": len(x) / fs, "flags": []}
    flags = row["flags"]

    # ---- A1 format
    if x.shape[1] == 1:
        flags.append("mono")
        x = np.repeat(x, 2, axis=1)
    elif x.shape[1] > 2:
        flags.append(f"{x.shape[1]}ch_used_first_2")
        x = x[:, :2]
    if fs < 32000:
        flags.append("low_sample_rate")
    clip = int(np.sum(np.abs(x) >= 0.999))
    row["clipped_samples"] = clip
    if clip:
        flags.append("clipped")
    for c, nm in enumerate("LR"):
        row[f"peak_dbfs_{nm}"] = float(20 * np.log10(np.abs(x[:, c]).max() + 1e-20))
        row[f"dc_dbfs_{nm}"] = float(20 * np.log10(abs(np.mean(x[:, c])) + 1e-20))
    if max(row["dc_dbfs_L"], row["dc_dbfs_R"]) > -60:
        flags.append("dc_offset")
    if np.abs(x).max() == 0:
        flags.append("silent")
        return row, None

    # ---- onset, noise, tail shape (A2/A3)
    # onset: first -40 dB crossing within 2 ms before the first -20 dB crossing (robust to stray pre-onset noise)
    a = np.abs(x).max(axis=1)
    main = int(np.argmax(a > a.max() * 0.1))
    lo = max(main - int(0.002 * fs), 0)
    onset = lo + int(np.argmax(a[lo:main + 1] > a.max() * 10 ** (-40 / 20)))
    row["onset_ms"] = onset / fs * 1000
    nz = np.nonzero(np.any(x != 0, axis=1))[0]
    zero_tail = len(x) - 1 - nz[-1]
    row["zero_tail_ms"] = zero_tail / fs * 1000
    xe = x[:len(x) - zero_tail] if zero_tail else x
    w = int(0.02 * fs)
    for c, nm in enumerate("LR"):
        e = np.convolve(xe[:, c] ** 2, np.ones(w) / w, "same")
        row[f"noise_dbfs_{nm}"] = float(10 * np.log10(np.mean(xe[-int(0.15 * fs):, c] ** 2) + 1e-30))
        seg = e[-int(0.3 * fs):-w] if len(e) > int(0.4 * fs) else e
        tt = np.arange(len(seg)) / fs
        row[f"end_slope_db_s_{nm}"] = float(np.polyfit(tt, 10 * np.log10(seg + 1e-30), 1)[0]) if len(seg) > 10 else None
    pre_end = onset - int(0.0005 * fs)
    if pre_end > int(0.002 * fs):
        pre = 10 * np.log10(np.mean(x[:pre_end] ** 2) + 1e-30)
        row["pre_onset_minus_noise_db"] = float(pre - min(row["noise_dbfs_L"], row["noise_dbfs_R"]))
        if row["pre_onset_minus_noise_db"] > 6:
            flags.append("pre_onset_energy")
    else:
        row["pre_onset_minus_noise_db"] = None
        flags.append("no_pre_roll")

    # ---- decay (C2): measured on the capture itself
    for c, nm in enumerate("LR"):
        d_bb, nf = edc(x[:, c], fs)
        row[f"T30_broadband_{nm}"] = rt_from_edc(d_bb, fs, -5, -35) or rt_from_edc(d_bb, fs, -5, -25)
        for b, (lo, hi) in BANDS.items():
            r = band_rt(x[:, c], fs, lo, hi)
            row[f"T_{b}_{nm}"] = r["T"]
            row[f"EDT_{b}_{nm}"] = r["EDT"]
        for fc in OCTAVES:
            row[f"T_oct{fc}_{nm}"] = band_rt(x[:, c], fs, fc / np.sqrt(2), fc * np.sqrt(2))["T"]
    mean = lambda *v: float(np.mean([q for q in v if q is not None])) if any(q is not None for q in v) else None
    rt_bb = mean(row["T30_broadband_L"], row["T30_broadband_R"])
    row["rt_mid"] = mean(row["T_mid_L"], row["T_mid_R"])
    lf, hf = mean(row["T_lf_L"], row["T_lf_R"]), mean(row["T_hf_L"], row["T_hf_R"])
    row["lf_ratio"] = lf / row["rt_mid"] if (lf and row["rt_mid"]) else None
    row["hf_ratio"] = hf / row["rt_mid"] if (hf and row["rt_mid"]) else None

    # ---- A2 tail classification
    if rt_bb:
        expected = -60 / rt_bb
        slopes = [s for s in (row["end_slope_db_s_L"], row["end_slope_db_s_R"]) if s is not None]
        s_end = min(slopes) if slopes else 0.0
        if zero_tail > int(0.001 * fs):
            flags.append("zero_tail")
        if s_end < 1.5 * expected:
            flags.append("fade_out")
        elif s_end < 0.5 * expected:
            flags.append("truncated_before_noise_floor")
        if row["duration_s"] < onset / fs + 1.0 * rt_bb:
            flags.append("short_file")

    # ---- device grid (B1/B4) and cluster fingerprint (A4)
    nwin = min(len(x), onset + int(cfg["analysis_window_s"] * fs))
    xw = np.vstack([np.zeros((PAD, 2)), x[:nwin]])
    d = deconvolve(xw, fs, kernel)
    on_p = onset + PAD
    lo_i, hi_i = on_p - int(0.001 * fs), on_p + int(0.016 * fs)
    dist = int(20 * fs / fsd)
    lobe = int(1.5 * fs / fsd) + 1
    taps, sparsity = [], []
    for c in range(2):
        seg = d[lo_i:hi_i, c]
        pk, _ = find_peaks(np.abs(seg), height=0.3 * np.abs(seg).max(), distance=dist)
        top = sorted(pk, key=lambda p: -abs(seg[p]))[:4]
        mask = np.zeros(len(seg), bool)
        for p in top:
            mask[max(p - lobe, 0):p + lobe + 1] = True
        sparsity.append(float(np.sum(seg[mask] ** 2) / (np.sum(seg ** 2) + 1e-30)))
        taps.append(sorted(subsample_peak(d[:, c], lo_i + p) for p in top))
    row["cluster_sparsity"] = min(sparsity)
    has_cluster = row["cluster_sparsity"] >= 0.55 and all(len(t) == 4 for t in taps)
    row["has_cluster"] = has_cluster
    if has_cluster:
        fp = []
        for t in taps:
            t0 = t[0][0]
            offs = [(ti - t0) / fs * fsd for ti, _ in t]
            fp.append(([int(round(o)) for o in offs], [int(np.sign(v)) for _, v in t]))
        row["cluster_status"] = classify_cluster(fp, cfg["cluster_templates"])
        row["cluster_offsets_L"], row["cluster_offsets_R"] = fp[0][0], fp[1][0]
        row["cluster_signs_L"], row["cluster_signs_R"] = fp[0][1], fp[1][1]
        gains = [abs(v) for t in taps for _, v in t]
        row["cluster_gain_db"] = float(20 * np.log10(np.median(gains)))
        row["cluster_gain_spread_db"] = float(20 * np.log10(max(gains) / min(gains)))
        # B1: device clock from tap spacings (both channels, LS through origin)
        n_all, dt_all, frac = [], [], []
        for t in taps:
            for ti, _ in t[1:]:
                dt = (ti - t[0][0]) / fs
                n = round(dt * fsd)
                n_all.append(n); dt_all.append(dt); frac.append(abs(dt * fsd - n))
        n_all, dt_all = np.array(n_all), np.array(dt_all)
        fs_meas = float(np.sum(n_all * dt_all) / np.sum(dt_all ** 2))
        row["clock_ppm"] = (fs_meas / fsd - 1) * 1e6
        row["grid_max_frac"] = float(max(frac))
        if row["grid_max_frac"] > 0.25:
            flags.append("off_device_grid")
        # B4: L/R output offset from the shared tap 0
        row["lr_offset_us"] = (taps[0][0][0] - taps[1][0][0]) / fs * 1e6
        row["tap0_ms"] = (taps[1][0][0] - PAD) / fs * 1000
        anchor_R, off = taps[1][0][0], row["lr_offset_us"]
    else:
        row["cluster_status"] = "no_cluster"
        anchor_R, off = float(on_p), lr_offset_us
    k0 = int(np.floor(anchor_R / fs * fsd))
    h = np.stack([to_device_grid(d[:, 0], fs, anchor_R + off * 1e-6 * fs, k0, cfg),
                  to_device_grid(d[:, 1], fs, anchor_R, k0, cfg)], 1)

    # ---- mixing time, reverb level, cluster/reverb ratio
    mix = echo_density_mix_time(x[:, 1], fs, onset)
    row["mix_ms"] = mix * 1000 if mix is not None else None
    if mix is not None:
        km = k0 + int((onset + PAD + mix * fs - anchor_R) / fs * fsd)
        rv = h[km:km + int(0.05 * fsd)]
        row["reverb_level_db"] = float(10 * np.log10(np.mean(rv ** 2) + 1e-30))
        if has_cluster:
            row["cluster_to_reverb_db"] = row["cluster_gain_db"] - row["reverb_level_db"]

    # ---- loops (C2 program identity)
    if mix is not None and rt_bb:
        start = km + int(0.02 * fsd)
        d_bb, nf = edc(x[:, 1], fs)
        stop_db = max(-45.0, nf + 10)
        t_stop = (np.argmax(d_bb < stop_db) if np.any(d_bb < stop_db) else len(d_bb)) / fs
        end = min(len(h), k0 + int((onset + PAD + t_stop * fs - anchor_R) / fs * fsd))
        if end - start >= 6000:
            periods, acfs = find_loops(h, start, end, rt_bb, fsd)
            row["loops"] = periods
            row["program_measured"] = match_loopset(periods, cfg["known_loopsets"])
        else:
            row["loops"], row["program_measured"] = [], "?"
            flags.append("tail_too_short_for_loops")

    # ---- low-pass corner (C2)
    row["lp_corner_khz"], row["lp_method"] = _lp_corner(d, x, fs, taps, has_cluster, onset, mix)

    # ---- signature (A5/A7)
    sig = h[k0 - 20:k0 - 20 + int(SIG_S * fsd)].astype(np.float32) if k0 >= 20 else None
    # noise relative to the signature window, both measured in the device band (< F_HI) like the signature
    sos = butter(8, F_HI, fs=fs, output="sos")
    xb = sosfilt(sos, xe, axis=0)
    nb = max(np.mean(xb[-int(0.15 * fs):, c] ** 2) for c in range(2))
    sig_win = xb[onset:onset + int(SIG_S * fs)]
    row["noise_rel_db"] = float(10 * np.log10(nb + 1e-30) - 10 * np.log10(np.mean(sig_win ** 2) + 1e-30))
    row["anchored"] = has_cluster
    return row, sig


def _lp_corner(d, x, fs, taps, has_cluster, onset, mix):
    """-3 dB corner of the low-pass. Cluster tap 0 (R) is a clean view of LP x analogue chain."""
    if has_cluster:
        t0 = int(taps[1][0][0])
        gap = int((taps[1][1][0] - taps[1][0][0]) * 0.9)
        seg = d[t0 - int(0.001 * fs):t0 + gap, 1] * np.hanning(int(0.001 * fs) + gap)
        method = "cluster_tap"
        thr = 10 ** (-3 / 20)
    else:
        if mix is None:
            return None, None
        s0 = onset + int(mix * fs)
        seg = x[s0:s0 + int(0.1 * fs), 1] * np.hanning(int(0.1 * fs))
        method = "reverb_spectrum"
        thr = 10 ** (-6 / 20)
    N = 1 << 16
    S = np.abs(np.fft.rfft(seg, N))
    f = np.fft.rfftfreq(N, 1 / fs)
    k = max(int(200 / (f[1] - f[0])), 1)
    S = np.convolve(S, np.ones(k) / k, "same")
    ref = np.median(S[(f > 300) & (f < 1500)])
    above = np.where((f > 1000) & (f < F_HI) & (S >= thr * ref))[0]
    if not len(above):
        return None, method
    fc = f[above[-1]] / 1000
    return (float(fc) if fc < F_ROLL / 1000 * 0.98 else float(F_ROLL / 1000)), method + ("" if fc < F_ROLL / 1000 * 0.98 else "_at_band_edge")
