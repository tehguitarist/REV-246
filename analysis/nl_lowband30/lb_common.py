"""Shared loader for the NonLin first-30-ms low-band investigation (read-only on tools/ and analysis/nl/)."""
import os, sys, json, glob
os.environ.setdefault("NL_BK", "0")
import numpy as np
from scipy.signal import butter, sosfiltfilt
ROOT = "/Users/leighpierce/Git/246"
sys.path.insert(0, os.path.join(ROOT, "tools"))
import nl_model as M          # noqa: E402  (import only; nothing is written)
import nl_early as ne         # noqa: E402
import nl_late as nlt         # noqa: E402
import nl_structure as nl     # noqa: E402

FS, K0 = M.FS, M.K0
BENCH = os.path.join(ROOT, "analysis", "nl", os.environ.get("LB_BENCH", "bench"))      # read only
OUTD = os.path.join(ROOT, "analysis", "nl_lowband30")
NMAX = 4096                                                # samples from file start that every analysis needs


def bands(lo, hi, order=4):
    return butter(order, [lo, hi], btype="band", fs=FS, output="sos")


def bench_all():
    out = {}
    for f in sorted(glob.glob(os.path.join(BENCH, "*.json"))):
        k = tuple(float(v) for v in os.path.basename(f)[:-5].split("_"))
        out[k] = json.load(open(f))
    return out


def load(key, r, n=NMAX, proc=True, band=False):
    """(y, m, h): capture (unwarped by the capture's clock) and model at the bench alignment, each (n, 2), full band,
    channel gains from the bench, so e = y - m is the bench residual. h = raw cache (not unwarped)."""
    raw = nl.ir(nl.table()[key])
    rawp = np.vstack([raw[:n + 64], np.zeros((max(0, n + 64 - len(raw)), 2))])
    y = nlt.unwarp(rawp, r["ppm"])[:n]
    m = M.render(key, n, r["o"], r["front_shift"], band=band, proc=proc) * np.array(r["gain"])
    return y, m, raw


def render_parts(key, n, o, fr, proc=True):
    """Model split into the front part (front taps + R direct term) and the branch part (AP(124) branch: main taps + late), each (n, 2)."""
    N = 1 << int(np.ceil(np.log2(n + 2 * ne.TAPS[-1] + 4096)))
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    R = M._proc(N) if proc else np.ones((len(ph), 2))
    f, b = [], []
    for ch in (0, 1):
        F, Br = M.spectra(key, ch, N)
        sh = np.exp(ph * (K0 + fr)) * R[:, ch]
        f.append(np.fft.irfft(F * sh, N)[:n]); b.append(np.fft.irfft(Br * np.exp(ph * o) * sh, N)[:n])
    return np.array(f).T, np.array(b).T


def frame(y, fr, n=3000):
    """y re-sampled so that the front-frame zero (file index K0 + fr) sits at index K0 (fractional shift in the FFT domain)."""
    N = 8192; f = np.fft.rfftfreq(N)
    Y = np.fft.rfft(np.vstack([y[:n], np.zeros((N - n, y.shape[1]))]), axis=0) * np.exp(2j * np.pi * f * fr)[:, None]
    return np.fft.irfft(Y, N, axis=0)[:n]


def stage_columns(key, ch, N, r, with_deriv=True):
    """Model split per stage: list of complex spectra (N//2+1) at the bench alignment (front part separate, R direct term in it), through proc.
    Returns dict name -> spectrum (unit gain scale 'as modelled'); derivative columns (time-shift) as name+'_d'."""
    span, early, late = key
    ns = ne.N_STAGES[span]
    ph = -2j * np.pi * np.fft.rfftfreq(N)
    R = M._proc(N)[:, ch]
    sh = np.exp(ph * (K0 + r["front_shift"])) * R
    fdb, sdb = ne.law_db(span, early, ch)
    d = ne.diffuser(ch, N // 2); Df = np.fft.rfft(d, N)
    out = {}
    front = sum(ne.FRONT_AMP[ch][i] * 10 ** (fdb[i] / 20) * np.exp(ph * f) for i, f in enumerate(ne.FRONT[ch]))
    out["front"] = front * sh
    for t, s, i in ne.taps(ch, ns):
        j = [jj for jj, st in enumerate(ne.STAGES) if i in st][0]
        g = s * 10 ** (sdb[j] / 20)
        if ch == 1 and i == 0:
            out["front"] = out["front"] + g * d[0] * np.exp(ph * ne.R_FIRST_DIRECT) * sh
            c = g * (Df - d[0]) * np.exp(ph * (t + r["o"])) * sh
        else:
            c = g * Df * np.exp(ph * (t + r["o"])) * sh
        k = f"s{j}"
        out[k] = out.get(k, 0) + c
        if with_deriv: out[k + "_d"] = out.get(k + "_d", 0) + c * ph * 1.0       # d/dshift (samples)
    return out
