"""Absolute AMPL of every Altiverb factory preset from the measured AMPL law (tools/ml_ampl_check.py: linear through zero,
plugin = hardware to 0.01 dB after the M / L fix): the preset IR's cluster-to-reverb ratio against the plugin's render at
the preset's own settings at a reference AMPL (S / M / L 90, E 80), AMPL = A_ref x 10^((ratio_IR - ratio_render) / 20).

    .venv/bin/python tools/ml_ampl_presets.py [--only <substring>]

Measure (48 kHz, 100 Hz-10 kHz, the IR and the render through the same function):
  cluster  R channel rms over its direct taps, from 1 ms before tap 0 to 1 ms after the last R tap (before R's reverb onset:
           S / M / L R taps end at 353 / 704 / 1407, R onsets 572 / 1134 / 2268 device samples; E R taps 0 / 429 / 718,
           its reverb from about 26 ms (fitted below the first read))
  reverb   rms of both channels in the R first pass [onset - 3, onset + Pmin - 10) (R channel) and in a late window from
           max(40 ms, last cluster tap + 2 ms) to 120 ms (both channels)
  bands    the same ratio in 100-1000 Hz and 1-10 kHz: a capture-chain tilt would make the two disagree
Writes analysis/ml_ampl/preset_ampl.json and preset_measure.json.
"""
import argparse, json, os, subprocess, sys
import numpy as np
import soundfile as sf
from scipy.signal import butter, resample_poly, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import altiverb_presets as ap  # noqa: E402

FS, FSD = 48000, 19.2e6 / 728
OUT = os.path.join(ROOT, "analysis", "ml_ampl")
REN = os.path.join(OUT, "preset_renders")
EXE = os.path.join(ROOT, "plugin", "build-render", "REV246Render_artefacts", "Release", "REV246Render")
PRESETS = json.load(open(os.path.join(ROOT, "analysis", "altiverb", "presets.json")))["presets"]
BANDS = {"full": (100, 10000), "lo": (100, 1000), "hi": (1000, 10000)}
SOS = {k: butter(4, v, btype="band", fs=FS, output="sos") for k, v in BANDS.items()}
R_END = {"S": 353, "M": 704, "L": 1407, "E": 0}        # E: tap 0 only (its reverb starts ~8 ms after it, before tap 718)
CL_END = {"S": 353, "M": 704, "L": 1407, "E": 718 + 314}
R_ON = {"S": 572, "M": 1134, "L": 2268, "E": 690}      # E: placeholder, replaced by the measured onset below
PMIN = {"S": 1321, "M": 2641, "L": 5279, "E": 2215}
AREF = {"S": 90, "M": 90, "L": 90, "E": 80}


def tap0(x):
    """Sub-sample time (48 kHz samples) of R's first direct tap (the first peak above half the early maximum)."""
    a = np.abs(x[1, : int(0.6 * FS)])
    i = int(np.argmax(a > 0.5 * a.max()))
    return ap.peak_time(x[1], i + 3)[0]


def measure(x, prog, e_onset=None):
    x = np.asarray(x, float)
    t0 = tap0(x)
    d = lambda n: n / FSD * FS                                   # device samples -> 48 kHz samples
    on = e_onset if prog == "E" else d(R_ON[prog])
    out = {"t0_ms": t0 / FS * 1e3}
    for b, sos in SOS.items():
        y = sosfiltfilt(sos, np.hstack([np.zeros((2, 4000)), x[:, : int(t0) + FS]]), axis=1)
        T = t0 + 4000
        c = np.sqrt(np.mean(y[1, int(T - 0.001 * FS): int(T + d(R_END[prog]) + 0.001 * FS)] ** 2))
        fp = y[1, int(T + on - d(3)): int(T + on + d(PMIN[prog] - 10))]
        w0 = max(0.040 * FS, d(CL_END[prog]) + 0.002 * FS)
        lt = y[1, int(T + w0): int(T + 0.120 * FS)]       # R only: the Altiverb's L channel runs ~2 dB hot (tap 0 L / R)
        out[b] = {"fp": float(20 * np.log10(c / np.sqrt(np.mean(fp ** 2)))),
                  "late": float(20 * np.log10(c / np.sqrt(np.mean(lt ** 2))))}
    return out


def e_reverb_onset(x):
    """E: the first sample after the last R direct tap where R exceeds -30 dB of its 0-0.3 s peak (direct taps masked)."""
    x = np.asarray(x, float)
    t0 = tap0(x)
    a = np.abs(sosfiltfilt(SOS["full"], x[1, : int(t0 + 0.3 * FS)]))
    m = a.copy()
    m[: int(t0 + 718 / FSD * FS + 0.001 * FS)] = 0
    return float(np.argmax(m > m.max() * 10 ** (-30 / 20)) - t0)


def render(name, st, prog, ampl, domain="host", mod=0, pre=0.0):
    os.makedirs(REN, exist_ok=True)
    tag = name.replace(" ", "_") + f"_a{ampl}_m{mod}_p{pre}" + ("" if domain == "host" else "_dev")
    path = os.path.join(REN, tag + ".wav")
    a = [EXE, f"out={path}", f"rate={FS}", f"program={prog}", f"lf={st['lf']}", f"hf={st['hf']}", f"ampl={ampl}",
         f"delay={st['delay']}", f"lp={st['lowpass']}", "sixteen=0", f"mod={mod}", f"len={0.8 + pre}", f"pre={pre}", f"domain={domain}"]
    if st.get("rt") is not None:
        a.append(f"rt={st['rt']}")
    o = subprocess.run(a, capture_output=True, text=True, check=True).stdout
    lat = int(o.split("latency=")[1].split()[0]) if domain == "host" else 0
    x, fs = sf.read(path, always_2d=True)
    x = x[lat:].T.astype(np.float64)                             # the printed latency already includes pre
    if domain != "host":                                          # device rate file -> 48 kHz
        x = _to48(x, fs)
    return x


def _to48(x, fs):
    from scipy.signal import resample
    return resample(np.hstack([np.zeros((2, 256)), x, np.zeros((2, 256))]), int(round((x.shape[1] + 512) * FS / fs)), axis=1)


def snap(a):
    return int(min(100, max(0, 10 * round(a / 10))))


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--only", nargs="*")
    args = ap_.parse_args()
    res, meas = {}, {}
    for name, p in sorted(PRESETS.items()):
        prog = p.get("program")
        if not p.get("included") or prog not in ("S", "M", "L", "E"):
            continue
        if args.only and not any(o in name for o in args.only):
            continue
        st = {k: p["settings"][k]["value"] for k in ("lf", "rt", "hf", "delay", "lowpass")}
        old = p["settings"]["ampl"]["value"]
        ir = np.asarray(ap.load(name), float)
        y = render(name, st, prog, AREF[prog])
        eo = None
        if prog == "E":                                      # E's moving reads carry ~6 dB of its reverb: mod=1, 3 phases
            eo = e_reverb_onset(y)
            ys = [render(name, st, prog, AREF[prog], mod=1, pre=q) for q in (0.0, 0.37, 0.9)]
            ms = [measure(v, prog, eo) for v in ys]
            mr = {b: {w: float(-10 * np.log10(np.mean([10 ** (-q[b][w] / 10) for q in ms]))) for w in ("fp", "late")} for b in BANDS}
        else:
            mr = measure(y, prog)
        mi = measure(ir, prog, eo)
        dlt = {b: {w: mi[b][w] - mr[b][w] for w in ("fp", "late")} for b in BANDS}
        amp = {b: {w: AREF[prog] * 10 ** (v / 20) for w, v in dw.items()} for b, dw in dlt.items()}
        meas[name] = dict(prog=prog, ir=mi, render=mr, delta_db=dlt, ampl=amp)
        a_fp, a_lt = amp["full"]["fp"], amp["full"]["late"]
        band_sp = max(abs(dlt["lo"][w] - dlt["hi"][w]) for w in ("fp", "late"))
        win_sp = abs(dlt["full"]["fp"] - dlt["full"]["late"])
        m = float(np.sqrt(a_fp * a_lt))
        reasons = []
        if win_sp > 1.0:
            reasons.append(f"first pass and late window disagree by {win_sp:.1f} dB")
        if band_sp > 1.0:
            reasons.append(f"100-1k vs 1-10k bands disagree by {band_sp:.1f} dB")
        nearest = abs(m - snap(m))
        if nearest > 3.0:
            reasons.append(f"{m:.1f} % is {nearest:.1f} from the step")
        if prog == "E":
            reasons.insert(0, "E's AMPL law (direct x AMPL / 80) is not measured")
        lab = "H" if reasons else "S"
        why = (f"first pass {a_fp:.1f} %, 40-120 ms {a_lt:.1f} % vs the plugin render at AMPL {AREF[prog]} "
               f"(linear law, plugin = hardware 0.01 dB); bands lo/hi {amp['lo']['fp']:.1f}/{amp['hi']['fp']:.1f} %")
        res[name] = dict(program=prog, measured=round(m, 1), step=snap(m), old_step=old, label=lab,
                         reason=why + ("; " + "; ".join(reasons) if reasons else ""))
        print(f"{name:24s} {prog} old {old:>3} new {snap(m):>3} ({m:6.1f}; fp {a_fp:6.1f} late {a_lt:6.1f}; "
              f"lo {amp['lo']['fp']:6.1f} hi {amp['hi']['fp']:6.1f}) {lab}", flush=True)
    tag = "_smoke" if args.only else ""
    json.dump(res, open(os.path.join(OUT, f"preset_ampl{tag}.json"), "w"), indent=1)
    json.dump(meas, open(os.path.join(OUT, f"preset_measure{tag}.json"), "w"))


if __name__ == "__main__":
    main()
