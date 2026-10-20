"""EMT 246 impulse-response analyser.

Usage:
    .venv/bin/python tools/analyze_ir.py irs/<file>.wav [--rt-guess 1.0] [--out analysis/<name>]

Produces <out>/summary.json, <out>/report.png and <out>/h_dev.npy (device-rate IR).

Pipeline (each stage is a function so it can be reused in batch / fitting scripts):
  1. load + basic stats (onset, peak, noise floor, out-of-band energy)
  2. device-grid recovery: deconvolve the analogue chain using the first cluster tap as the
     kernel and resample on the unit's own sample clock (default 19.2 MHz / 728)
  3. cluster extraction: integer-sample taps + gains of the 8-reflection cluster
  4. decay analysis: Schroeder EDC, T20/T30/EDT per octave band, noise-compensated
  5. echo density profile (Abel & Huang normalised echo density)
  6. loop identification: autocorrelation of the decay-flattened tail -> recirculation periods,
     validated by checking that spectral modes sit on each loop's harmonic series
  7. per-loop L/R output offset (circular cross-correlation of periodic components)
  8. time-invariance test: spectral fine-structure correlation between decay windows
"""
import argparse, json, os, sys
import numpy as np
import soundfile as sf
import scipy.sparse as sp
from scipy.signal import butter, sosfilt, sosfiltfilt, find_peaks
from scipy.sparse.linalg import lsqr

sys.path.insert(0, os.path.dirname(__file__))
from deconv_grid import recover  # noqa: E402

FS_DEV = 19.2e6 / 728
OCTAVES = [63, 125, 250, 500, 1000, 2000, 4000, 8000]


# ---------- 1. basics ----------
def basics(x, fs):
    out = {"fs": fs, "channels": x.shape[1], "duration_s": len(x) / fs, "per_channel": []}
    for c in range(x.shape[1]):
        s = x[:, c]
        thr = np.abs(s).max() * 10 ** (-40 / 20)
        tail = s[-int(0.1 * fs):]
        out["per_channel"].append({
            "onset_ms": float(np.argmax(np.abs(s) > thr) / fs * 1000),
            "peak": float(np.abs(s).max()),
            "tail_rms_dbfs": float(20 * np.log10(np.sqrt(np.mean(tail ** 2)) + 1e-20)),
        })
    X = np.fft.rfft(x[:, 0]); F = np.fft.rfftfreq(len(x), 1 / fs)
    ref = np.mean(np.abs(X[(F > 500) & (F < 2000)]) ** 2)
    out["energy_above_13k2_db"] = float(10 * np.log10(np.mean(np.abs(X[F > 13200]) ** 2) / ref))
    return out


# ---------- 3. cluster ----------
def cluster_taps(h, thresh=0.5, span=1200):
    res = []
    for c in range(h.shape[1]):
        a = np.abs(h[:, c]); k0 = int(np.argmax(a > thresh))
        pk, _ = find_peaks(a[k0 - 2:k0 + span], height=thresh, distance=4)  # start early so tap 0 is a peak
        res.append({"first_tap_index": k0,
                    "taps": [{"n": int(p) - 2, "gain": round(float(h[k0 + p - 2, c]), 3)} for p in pk]})
    return res


# ---------- 4. decay ----------
def edc(s, fs, noise_ms=150):
    e = s ** 2
    nf = np.mean(e[-int(noise_ms / 1000 * fs):])
    w = int(0.01 * fs); sm = np.convolve(e, np.ones(w) / w, "same")
    above = np.where(sm > nf * 2)[0]
    cut = above[-1] if len(above) else len(e)
    ee = np.clip(e[:cut] - nf, 0, None)
    d = np.cumsum(ee[::-1])[::-1]
    return 10 * np.log10(d / d[0] + 1e-30), float(10 * np.log10(nf / sm.max()))


def rt_from_edc(d, fs, hi, lo):
    t = np.arange(len(d)) / fs; m = (d <= hi) & (d >= lo)
    if m.sum() < 10:
        return None
    return float(-60 / np.polyfit(t[m], d[m], 1)[0])


def decay_bands(x, fs, f_cap=12500):
    out = []
    for c in range(x.shape[1]):
        bands = {}
        d, nf = edc(x[:, c], fs)
        bands["broadband"] = {"T20": rt_from_edc(d, fs, -5, -25), "T30": rt_from_edc(d, fs, -5, -35),
                              "EDT": rt_from_edc(d, fs, 0, -10), "noise_db": nf}
        for fc in OCTAVES:
            sos = butter(4, [fc / np.sqrt(2), min(fc * np.sqrt(2), f_cap)], btype="band", fs=fs, output="sos")
            d, nf = edc(sosfilt(sos, x[:, c]), fs)
            bands[str(fc)] = {"T20": rt_from_edc(d, fs, -5, -25), "T30": rt_from_edc(d, fs, -5, -35),
                              "EDT": rt_from_edc(d, fs, 0, -10), "noise_db": nf}
        out.append(bands)
    return out


# ---------- 5. echo density ----------
def echo_density(s, fs, win=0.02):
    w = int(win * fs); hw = np.hanning(w); hw /= hw.sum(); t, v = [], []
    for i in range(0, len(s) - w, w // 4):
        seg = s[i:i + w]; sd = np.sqrt(np.sum(hw * seg ** 2))
        t.append(i / fs + win / 2); v.append(np.sum(hw * (np.abs(seg) > sd)) / 0.3173)
    return np.array(t), np.array(v)


# ---------- 6/7. loops ----------
def flatten(s, fs, rt):
    return s * 10 ** (3 * (np.arange(len(s)) / fs) / rt)


def find_loops(h, k0, rt, fs=FS_DEV, start=600, end_s=0.9, max_lag=8000, n_sigma=8):
    """Return fundamental recirculation periods (samples) shared by both channels."""
    found = []
    for c in range(h.shape[1]):
        a = flatten(h[:, c], fs, rt)[k0 + start:k0 + int(end_s * fs)]; a = a - a.mean()
        N = 1 << int(np.ceil(np.log2(2 * len(a))))
        ac = np.fft.irfft(np.abs(np.fft.rfft(a, N)) ** 2, N)[:max_lag]; ac /= ac[0]
        sd = np.std(ac[100:])
        pk, _ = find_peaks(ac, height=n_sigma * sd)
        pk = [int(p) for p in pk if p > 200]
        fund = []
        for p in sorted(pk):  # ascending, so multiples of a loop are recognised as harmonics
            if not any(abs(p / f - round(p / f)) * f <= 2 and round(p / f) >= 2 for f in fund):
                fund.append(p)
        found.append({"periods": sorted(fund), "acf": {str(p): round(float(ac[p]), 3) for p in sorted(fund)}, "acf_sigma": float(sd)})
    common = sorted(set(found[0]["periods"]) & set(found[-1]["periods"]))
    return common, found


def harmonic_check(h, k0, rt, periods, fs=FS_DEV, fmax=3000):
    """dB by which the tail spectrum at each loop's harmonics exceeds random frequencies."""
    rng = np.random.default_rng(1); res = []
    for c in range(h.shape[1]):
        s = flatten(h[:, c], fs, rt)[k0 + 600:k0 + int(1.1 * fs)]
        N = 1 << 21
        X = 20 * np.log10(np.abs(np.fft.rfft(s * np.hanning(len(s)), N)) + 1e-12); df = fs / N
        lvl = lambda fr: np.mean(X[np.round(fr[(fr > 100) & (fr < fmax)] / df).astype(int)])
        base = lvl(rng.uniform(100, fmax, 20000)); row = {}
        for p in periods:
            best = max((lvl(np.arange(1, 600) * fs / (p + dp)) - base, p + dp) for dp in np.linspace(-1.5, 1.5, 61))
            row[str(p)] = {"excess_db": round(float(best[0]), 2), "effective_len": round(float(best[1]), 2)}
        res.append(row)
    return res


def periodic_decomposition(h, k0, rt, periods, fs=FS_DEV, f_lp=3000, end_s=0.8):
    """Least-squares split of the flattened (<f_lp) tail into one periodic component per loop."""
    sos = butter(6, f_lp, fs=fs, output="sos"); comps = {}; explained = []
    for c in range(h.shape[1]):
        y = flatten(sosfiltfilt(sos, h[:, c]), fs, rt)[k0 + 500:k0 + int(end_s * fs)]
        n = np.arange(len(y)) + 500
        rows, cols, off = [], [], 0
        for p in periods:
            rows.append(np.arange(len(n))); cols.append(off + n % p); off += p
        A = sp.csr_matrix((np.ones(len(n) * len(periods)), (np.concatenate(rows), np.concatenate(cols))), shape=(len(n), off))
        sol = lsqr(A, y, iter_lim=2000)[0]
        explained.append(float(1 - np.var(y - A @ sol) / np.var(y)))
        o = 0
        for p in periods:
            comps[(c, p)] = sol[o:o + p]; o += p
    offsets = {}
    for p in periods:
        a, b = comps[(0, p)], comps[(1, p)]
        X = np.fft.irfft(np.fft.rfft(a) * np.conj(np.fft.rfft(b)), p) / np.sqrt(np.sum(a ** 2) * np.sum(b ** 2))
        i = int(np.argmax(np.abs(X)))
        offsets[str(p)] = {"lag": i if i < p / 2 else i - p, "corr": round(float(X[i]), 3)}
    return explained, offsets


# ---------- 8. time invariance ----------
def fine_structure_corr(s, fs, w1=(0.05, 0.35), w2=(0.35, 0.65), band=(200, 8000)):
    def spec(a, b):
        seg = s[int(a * fs):int(b * fs)]; seg = seg * np.hanning(len(seg))
        N = 1 << 20; X = np.abs(np.fft.rfft(seg, N)); f = np.fft.rfftfreq(N, 1 / fs)
        m = (f > band[0]) & (f < band[1]); L = 20 * np.log10(X[m] + 1e-12)
        k = int(50 / (f[1] - f[0])); return L - np.convolve(L, np.ones(k) / k, "same")
    return float(np.corrcoef(spec(*w1), spec(*w2))[0, 1])


# ---------- plotting ----------
def plot(x, fs, h, k0, rt, decay, ned, loops, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 2, figsize=(14, 12)); names = ["L", "R"]
    t = np.arange(len(x)) / fs
    for c in range(2):
        ax[0, 0].plot(t * 1000, 20 * np.log10(np.abs(x[:, c]) + 1e-9) - 3 * c, lw=0.4, label=names[c])
    ax[0, 0].set(xlim=(0, 120), ylim=(-70, 3), xlabel="ms", ylabel="dB (R offset -3)", title="Onset (96 kHz capture)"); ax[0, 0].legend()
    n = np.arange(1500)
    for c in range(2):
        ax[0, 1].plot(n, h[k0 + n, c] + (2.2 if c == 0 else 0), lw=0.6, color=f"C{c}")
    ax[0, 1].set(title="Device-rate IR after deconvolution (first 1500 samples from first tap; L raised)", xlabel=f"samples @ {FS_DEV:.1f} Hz")
    for c in range(2):
        d, _ = edc(x[:, c], fs); ax[1, 0].plot(np.arange(len(d)) / fs, d, label=f"{names[c]} broadband")
    ax[1, 0].set(ylim=(-80, 0), xlabel="s", ylabel="dB", title="Schroeder EDC"); ax[1, 0].legend()
    for c in range(2):
        ax[1, 1].semilogx(OCTAVES, [decay[c][str(f)]["T30"] or np.nan for f in OCTAVES], "o-", label=names[c])
    ax[1, 1].set(xlabel="Hz", ylabel="T30 (s)", title="Decay time per octave"); ax[1, 1].legend(); ax[1, 1].grid(True, which="both", alpha=.3)
    for c in range(2):
        ax[2, 0].plot(ned[c][0] * 1000, ned[c][1], label=names[c])
    ax[2, 0].set(xlim=(0, 300), ylim=(0, 1.3), xlabel="ms", title="Normalised echo density"); ax[2, 0].legend()
    a = flatten(h[:, 0], FS_DEV, rt)[k0 + 600:k0 + int(0.9 * FS_DEV)]; a -= a.mean()
    N = 1 << 17; ac = np.fft.irfft(np.abs(np.fft.rfft(a, N)) ** 2, N)[:7000]; ac /= ac[0]
    ax[2, 1].plot(ac, lw=0.5)
    for p in loops:
        ax[2, 1].axvline(p, color="r", ls=":", lw=0.8)
    ax[2, 1].set(ylim=(-0.1, 0.35), xlabel="lag (device samples)", title="Tail autocorrelation (L, decay-flattened); red = detected loops")
    fig.tight_layout(); fig.savefig(path, dpi=110)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wav"); ap.add_argument("--out"); ap.add_argument("--fs-dev", type=float, default=FS_DEV)
    a = ap.parse_args()
    out = a.out or os.path.join("analysis", os.path.splitext(os.path.basename(a.wav))[0])
    os.makedirs(out, exist_ok=True)
    x, fs = sf.read(a.wav, always_2d=True)
    summary = {"file": a.wav, "basics": basics(x, fs)}
    h, t0, offs = recover(x, fs, a.fs_dev)
    np.save(os.path.join(out, "h_dev.npy"), h)
    summary["device_grid"] = {"fs_dev": a.fs_dev, "channel_offsets_us": [o * 1e6 for o in offs]}
    cl = cluster_taps(h); summary["cluster"] = cl; k0 = cl[0]["first_tap_index"]
    summary["predelay_to_first_tap_ms_from_file_start"] = k0 / a.fs_dev * 1000
    decay = decay_bands(x, fs); summary["decay"] = decay
    rt = decay[0]["1000"]["T30"] or decay[0]["broadband"]["T30"]
    ned = [echo_density(x[:, c], fs) for c in range(2)]
    summary["echo_density_reaches_0.9_ms"] = [float(t[np.argmax(v > 0.9)] * 1000) for t, v in ned]
    loops, raw = find_loops(h, k0, rt); summary["loops"] = {"periods": loops, "per_channel": raw}
    if loops:
        summary["loops"]["harmonic_check"] = harmonic_check(h, k0, rt, loops)
        expl, offsets = periodic_decomposition(h, k0, rt, loops)
        summary["loops"]["periodic_model_explained_var"] = expl
        summary["loops"]["LR_offset_per_loop"] = offsets
    summary["time_invariance_fine_corr"] = [fine_structure_corr(x[:, c], fs) for c in range(2)]
    plot(x, fs, h, k0, rt, decay, ned, loops, os.path.join(out, "report.png"))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(json.dumps({k: summary[k] for k in ("device_grid", "cluster", "loops", "time_invariance_fine_corr")}, indent=1, default=float))


if __name__ == "__main__":
    main()
