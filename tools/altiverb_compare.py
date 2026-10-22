"""Altiverb 246 IRs against the Tomkin library: cluster timing, decay, spectrum and a reverb null.
Results: docs/library_notes.md §6, docs/EMT246_unit_facts.md §5c.

    .venv/bin/python tools/altiverb_compare.py [--stage all|deconv|cluster|decay|spectrum|null] [--out analysis/altiverb]

Input: sweep recordings of the plugin in irs/Altiverb EMT256 IRs/ (not committed), made with the stimulus of
tools/make_sweep.py (1 s lead-in, 10 s log sweep, any tail; 48 kHz). Writes <out>/summary.json and caches the
deconvolved IRs as <out>/cache/*.npy (not committed).

  deconv    Deconvolve each recording with the sweep's inverse filter (lag 0 at LEAD + sweep length). Reports
            the first arrival and the clipped-sample count; a pre-IR energy figure that mixes distortion with the
            method's leakage floor (about -47...-55 dB re IR RMS for a perfectly linear unit).
  cluster   Tap positions in device samples (19.2 MHz/728) for the NON-LIN SHORT (S), MED (M) and LONG (L)
            presets: peaks of the 16x-upsampled IR in the cluster window, relative to the first tap.
  decay     Octave-band T60 (Schroeder, -5...-25 dB) of each preset against the Tomkin files with the same label RT.
            Tomkin LF and HF of the Altiverb presets are unknown, so Tomkin is shown at the bracketing grid values.
  spectrum  Single-tap response (+-0.5 ms Hann window) against the CLA reference kernel, and octave-band levels of
            the 20-150 ms and 150-600 ms windows against the nearest Tomkin setting.
  null      Program S reverb only (Tomkin AMPL 0): Altiverb against Tomkin, aligned on the first Altiverb tap, with
            gain, shift and clock ratio fitted per channel on 20-150 ms (3 parameters); later windows held out.
            100 Hz-10 kHz. Program M and L reverbs are time-variant (ml_structure.md §4): no null is attempted.

Tomkin 96 kHz files are decimated to 48 kHz (they carry nothing above 14 kHz).
"""
import argparse, glob, json, os, sys
import numpy as np
import soundfile as sf
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from scipy.signal import butter, fftconvolve, find_peaks, resample_poly, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from make_sweep import LEAD, sweep_and_inverse  # noqa: E402

FS = 48000
FSD = 19.2e6 / 728  # device rate
ALT = os.path.join(ROOT, "irs", "Altiverb EMT256 IRs")
TOM = os.path.join(ROOT, "irs", "EMT 246 captures", "96kHz")
KERNEL = os.path.join(ROOT, "ref", "kernel_eugtone_96k.npz")
CF = [125, 250, 500, 1000, 2000, 4000, 8000]
# preset family -> program, and the cluster window (s) after the first tap
FAMILY = {"S": "NON-LIN SHORT RT", "M": "MED RT", "L": "LONG RT"}
SPAN = {"S": 0.020, "M": 0.040, "L": 0.080}
TEMPLATE = {  # device-sample taps (ml_structure.md §1; unit facts §5, §5b)
    "S": {"L": [0, 54, 194, 288], "R": [0, 143, 194, 353]},
    "M": {"L": [0, 108, 388, 576], "R": [0, 285, 388, 704]},
    "L": {"L": [0, 216, 776, 1152], "R": [0, 570, 776, 1407]},
}
BP = butter(4, [100, 10000], btype="band", fs=FS, output="sos")
db = lambda v: 20 * np.log10(v + 1e-30)


# ---------------------------------------------------------------- loading
def altiverb_ir(name, out):
    """Deconvolved Altiverb IR (2 x N, 48 kHz), cached. Lag 0 is at LEAD + sweep length (zero plugin latency assumed)."""
    cache = os.path.join(out, "cache", name + ".npy")
    if os.path.exists(cache):
        return np.load(cache)
    x, fs = sf.read(os.path.join(ALT, f"EMT256_{name}.wav"), always_2d=True)
    assert fs == FS, (name, fs)
    _, inv = sweep_and_inverse(FS)
    k0 = int(LEAD * FS) + len(inv) - 1
    ir = np.array([fftconvolve(x[:, c], inv)[k0:] for c in (0, 1)])
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.save(cache, ir)
    return ir


def tomkin(prog, amp, lf, rt, hf):
    d = {"S": "Short RT", "M": "Med RT", "L": "Long RT"}[prog]
    x, fs = sf.read(os.path.join(TOM, d, f"12_0_{amp}_{lf}_{rt}_{hf}_dc.wav"), always_2d=True)
    assert fs == 96000
    return resample_poly(x, 1, 2, axis=0).T


def first_peak(x):
    a = np.abs(x)
    return int(np.argmax(a > 0.4 * a[: int(0.4 * FS)].max()))


def presets():
    return sorted(os.path.basename(p)[7:-4] for p in glob.glob(os.path.join(ALT, "EMT256_*.wav")))


# ---------------------------------------------------------------- stages
def stage_deconv(out):
    res = {}
    _, inv = sweep_and_inverse(FS)
    for n in presets():
        x, _ = sf.read(os.path.join(ALT, f"EMT256_{n}.wav"), always_2d=True)
        ir = altiverb_ir(n, out)
        k0 = int(LEAD * FS) + len(inv) - 1
        fa = [float(np.argmax(np.abs(c) > np.abs(c).max() * 10 ** (-40 / 20)) / FS * 1e3) for c in ir]
        pre = []
        for c in (0, 1):
            d = fftconvolve(x[:, c], inv)
            w = d[k0 - 9 * FS: k0 - int(0.05 * FS)]
            pre.append(float(db(np.sqrt((w ** 2).mean())) - db(np.sqrt((ir[c][: 14 * FS] ** 2).mean()))))
        res[n] = {"clipped": int((np.abs(x) >= 0.99999).sum()), "peak_dbfs": float(db(np.abs(x).max())),
                  "first_arrival_ms": fa, "pre_ir_db": pre}
        print(f"{n:26s} clipped {res[n]['clipped']:6d}  peak {res[n]['peak_dbfs']:6.1f} dBFS  "
              f"first arrival {fa[0]:6.1f}/{fa[1]:6.1f} ms  pre-IR {pre[0]:6.1f}/{pre[1]:6.1f} dB")
    return res


def stage_cluster(out):
    res = {}
    for prog, pre in FAMILY.items():
        names = [n for n in presets() if n.startswith(pre + " ")]
        taps = {"L": [], "R": []}
        for n in names:
            ir = altiverb_ir(n, out)
            for ci, c in enumerate("LR"):
                x = ir[ci][: int(0.4 * FS)]
                fa = int(np.argmax(np.abs(x) > np.abs(x).max() * 10 ** (-40 / 20)))
                y = resample_poly(x[max(fa - 96, 0): fa + int(SPAN[prog] * FS * 1.4)], 16, 1)
                pk, _ = find_peaks(np.abs(y), height=0.4 * np.abs(y).max(), distance=int(0.0003 * FS * 16))
                if len(pk) == 4:
                    taps[c].append(((pk - pk[0]) / (FS * 16) * FSD, np.sign(y[pk]) * np.sign(y[pk[0]])))
        for c in "LR":
            T = np.array([t for t, _ in taps[c]])
            Sg = np.array([s for _, s in taps[c]])
            res[f"{prog}{c}"] = {"files": len(names), "with_4_peaks": len(T), "taps": T.mean(0).round(3).tolist(),
                                 "std": T.std(0).round(4).tolist(), "signs": Sg.mean(0).tolist(), "template": TEMPLATE[prog][c]}
            print(f"{prog} {c}: {len(T)}/{len(names)} files  taps {np.round(T.mean(0), 2)}  template {TEMPLATE[prog][c]}  "
                  f"max |diff| {np.abs(T.mean(0) - TEMPLATE[prog][c]).max():.3f}  signs {Sg.mean(0).astype(int)}")
    return res


def band_t60(ir):
    def t60(x, lo=-5, hi=-25):
        e = np.cumsum((x ** 2)[::-1])[::-1]
        e = 10 * np.log10(e / e[0] + 1e-30)
        tt = np.nonzero((e < lo) & (e > hi))[0]
        if len(tt) < 200 or tt[-1] / FS > len(x) / FS - 0.02:
            return np.nan  # decay does not reach -25 dB before the end of the file (Tomkin fades)
        return -60 / np.polyfit(tt / FS, e[tt], 1)[0]
    out = []
    for fc in CF:
        sos = butter(4, [fc / 2 ** 0.5, min(fc * 2 ** 0.5, 0.45 * FS)], btype="band", fs=FS, output="sos")
        out.append(float(np.nanmean([t60(sosfiltfilt(sos, ir[c][: 12 * FS])) for c in (0, 1)])))
    return out


# Altiverb preset -> program, Tomkin label RT, Tomkin variants shown (amp, LF, HF); LF/HF of the presets are unknown.
DECAY_CASES = [
    ("S", "NON-LIN SHORT RT 0.6s", "06", [(0, "10", "10"), (0, "10", "06")]),
    ("S", "NON-LIN SHORT RT 1.0s", "10", [(0, "10", "10"), (0, "10", "06")]),
    ("M", "MED RT 1.8s", "18", [(10, lf, hf) for lf in ("07", "12") for hf in ("07", "09")]),
    ("M", "MED RT 3.0s", "30", [(10, lf, hf) for lf in ("07", "12") for hf in ("07", "09")]),
    ("M", "MED RT 5.0s", "50", [(10, lf, hf) for lf in ("07", "12") for hf in ("07", "09")]),
    ("L", "LONG RT 6.0s", "60", [(0, lf, hf) for lf in ("08", "12") for hf in ("03", "05", "08", "10")]),
    ("L", "LONG RT 8.0s", "80", [(0, lf, hf) for lf in ("08", "12") for hf in ("03", "05", "08", "10")]),
]


def stage_decay(out):
    res = {}
    print("octave-band T60 (s) at", CF)
    for prog, alt, rt, variants in DECAY_CASES:
        r = {"Altiverb": band_t60(altiverb_ir(alt, out))}
        print(alt)
        print(f"  {'Altiverb':30s}", " ".join(f"{v:6.2f}" for v in r["Altiverb"]))
        for amp, lf, hf in variants:
            r[f"Tomkin AMPL{amp} LF{lf} HF{hf}"] = band_t60(tomkin(prog, amp, lf, rt, hf))
            print(f"  {'Tomkin AMPL%s LF%s HF%s' % (amp, lf, hf):30s}", " ".join(f"{v:6.2f}" for v in r[f'Tomkin AMPL{amp} LF{lf} HF{hf}']))
        res[alt] = r
    return res


def band_levels(ir, t1, t2):
    i0 = first_peak(ir[0])
    a, b = int(i0 + t1 * FS), int(i0 + t2 * FS)
    v = []
    for fc in CF:
        sos = butter(4, [fc / 2 ** 0.5, min(fc * 2 ** 0.5, 0.45 * FS)], btype="band", fs=FS, output="sos")
        v.append(10 * np.log10(np.mean([np.sum(sosfiltfilt(sos, ir[c][: 12 * FS])[a:b] ** 2) for c in (0, 1)])))
    v = np.array(v)
    return (v - v[3]).tolist()  # dB re the 1 kHz band


def stage_spectrum(out):
    res = {"tap": {}, "bands": {}}
    fr = np.array([1000, 4000, 8000, 10000, 11000, 12000, 13000, 14000, 16000.])

    def tap_response(x, fs):
        half = int(0.5e-3 * fs)
        k = int(np.argmax(np.abs(x)))
        seg = x[k - half: k + half + 1] * np.hanning(2 * half + 1)
        X = np.abs(np.fft.rfft(seg, 1 << 16))
        r = np.interp(fr, np.fft.rfftfreq(1 << 16, 1 / fs), X)
        return (20 * np.log10(r / r[0])).round(2).tolist()

    z = np.load(KERNEL)
    res["tap"]["CLA kernel"] = tap_response(np.pad(z["g"], (100, 100)), int(z["fs"]))
    for n in ("NON-LIN SHORT RT 1.0s", "MED RT 3.0s", "LONG RT 6.0s"):
        ir = altiverb_ir(n, out)[0][: int(0.4 * FS)]
        fa = int(np.argmax(np.abs(ir) > np.abs(ir).max() * 10 ** (-40 / 20)))
        up = resample_poly(ir[fa - 48: fa + 96], 8, 1)
        res["tap"][n] = tap_response(up[: 48 * 8 * 2], FS * 8)
    print("tap response dB re 1 kHz at", (fr / 1000).tolist(), "kHz")
    for k, v in res["tap"].items():
        print(f"  {k:24s}", " ".join(f"{x:7.2f}" for x in v))
    for nm, alt, prog, tk in (("S RT1.0", "NON-LIN SHORT RT 1.0s", "S", (0, "10", "10", "10")),
                              ("M RT3.0", "MED RT 3.0s", "M", (10, "12", "30", "07")),
                              ("L RT6.0", "LONG RT 6.0s", "L", (0, "12", "60", "05"))):
        A, T = altiverb_ir(alt, out), tomkin(prog, *tk)
        for w in ((0.02, 0.15), (0.15, 0.6)):
            res["bands"][f"{nm} {w}"] = {"Altiverb": band_levels(A, *w), f"Tomkin AMPL{tk[0]} LF{tk[1]} HF{tk[3]}": band_levels(T, *w)}
    print("band levels dB re 1 kHz at", CF)
    for k, v in res["bands"].items():
        for who, lv in v.items():
            print(f"  {k:24s} {who[:24]:24s}", " ".join(f"{x:6.2f}" for x in lv))
    return res


def reverb_null(alt, prog, tk, out, start=0.020, fit_end=0.15):
    ir = altiverb_ir(alt, out)
    T = tomkin(prog, *tk)
    wins = {"fit 20-150 ms": (start, fit_end), "HO 150-300": (0.15, 0.30), "HO 300-600": (0.30, 0.60), "HO 600-1000": (0.60, 1.0)}
    res = {}
    for c in (0, 1):
        a, b = ir[c][: 4 * FS], T[c][: 4 * FS]
        t0a = first_peak(a)
        A, B = sosfiltfilt(BP, a), sosfiltfilt(BP, b)
        lo, hi = int(t0a + start * FS), int(t0a + fit_end * FS)
        seg = A[lo:hi]
        cc = fftconvolve(B, seg[::-1], mode="valid")
        k = int(np.argmax(np.abs(cc)))
        t0b = k - (lo - t0a)
        cs = CubicSpline(np.arange(len(B)), B)

        def model(p, l, h):
            s, eps, g = p
            n = np.arange(l, h)
            return g * cs(np.clip(t0b + (n - t0a) * (1 + eps) + s, 0, len(B) - 1))

        cost = lambda p: np.sum((A[lo:hi] - model(p, lo, hi)) ** 2) / np.sum(A[lo:hi] ** 2)
        g0 = np.sign(cc[k]) * np.linalg.norm(seg) / np.linalg.norm(B[k: k + len(seg)])
        best = min((minimize(cost, [s0, 0, g0], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=6000))
                    for s0 in (-0.4, 0, 0.4)), key=lambda r: r.fun)
        o = {}
        for nm, (t1, t2) in wins.items():
            l, h = int(t0a + t1 * FS), int(t0a + t2 * FS)
            o[nm] = float(db(np.linalg.norm(A[l:h] - model(best.x, l, h)) / np.linalg.norm(A[l:h])))
        res["LR"[c]] = {"shift_samples": float(best.x[0]), "clock_ppm": float(best.x[1] * 1e6), "gain": float(best.x[2]), "null_db": o}
        print(f"  {'LR'[c]}: shift {best.x[0]:+.3f} smp, clock {best.x[1] * 1e6:+.1f} ppm, gain {best.x[2]:+.0f} | "
              + " | ".join(f"{k} {v:.1f}" for k, v in o.items()))
    return res


def stage_null(out):
    res = {}
    for alt, rt in (("NON-LIN SHORT RT 1.0s", "10"), ("NON-LIN SHORT RT 0.6s", "06")):
        print(f"{alt} vs Tomkin S AMPL 0 LF1.0 HF1.0 RT{rt} (Altiverb HF is nearer 0.9; 3 fitted parameters per channel)")
        res[alt] = reverb_null(alt, "S", (0, "10", rt, "10"), out)
    return res


STAGES = {"deconv": stage_deconv, "cluster": stage_cluster, "decay": stage_decay, "spectrum": stage_spectrum, "null": stage_null}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=["all", *STAGES])
    ap.add_argument("--out", default=os.path.join(ROOT, "analysis", "altiverb"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "summary.json")
    summary = json.load(open(path)) if os.path.exists(path) else {}
    for name, fn in STAGES.items():
        if a.stage in ("all", name):
            print(f"== {name}")
            summary[name] = fn(a.out)
            json.dump(summary, open(path, "w"), indent=1)


if __name__ == "__main__":
    main()
