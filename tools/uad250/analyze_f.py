"""F: parameter changes while ringing (docs/UAD250_test_plan.md §3 F, §4 item 4).

    .venv/bin/python tools/uad250/analyze_f.py [case ...]

Reads the F_auto renders and their inputs. For each case it reports the event, the local decay before / after, the
gain-envelope step and any sample jump at the event, and for the predelay / delay cases the instantaneous-frequency
and delay tracks. The aim is the implementation (instant, ramped or recomputed) and the delay-change behaviour
(a click, a crossfade or a pitch ramp), not a null.
"""
import os, sys, json
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt, hilbert

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "F_auto")
INP = os.path.join(ROOT, "analysis", "uad250", "captures", "inputs")
R = 48000
BANDS = [125, 250, 500, 1000, 2000, 4000, 8000]


def load(name):
    y, r = sf.read(os.path.join(CAP, name + ".wav"), always_2d=True)
    meta = json.load(open(os.path.join(CAP, name + ".json")))
    x, _ = sf.read(os.path.join(INP, "F_auto__" + name + ".wav"), always_2d=True)
    return y.astype(float), x.astype(float), meta, r


def bandpass(x, lo, hi):
    return sosfiltfilt(butter(4, [lo, hi], "band", fs=R, output="sos"), x, axis=0)


def envelope(x, win=0.005):
    """RMS in win-second blocks, at block centres (seconds)."""
    n = max(1, int(win * R))
    k = len(x) // n
    e = np.sqrt(np.mean(x[:k * n].reshape(k, n, -1) ** 2, axis=1))
    if e.ndim == 2 and e.shape[1] == 1:
        e = e[:, 0]
    return (np.arange(k) + 0.5) * n / R, e


def local_t60(x, t0, t1):
    """T60 from the log-RMS slope over [t0, t1], blocks within 20 dB of that window's own peak (dB/s -> s)."""
    t, e = envelope(x, 0.01)
    m = (t >= t0) & (t <= t1)
    if m.sum() < 8:
        return None
    m &= e > e[m].max() * 10 ** (-20 / 20)
    if m.sum() < 8:
        return None
    p = np.polyfit(t[m], 10 * np.log10(e[m] ** 2 + 1e-30), 1)
    return float(-60 / p[0]) if p[0] < 0 else None


def jump_at(y, t, half=0.02):
    """Largest |y[n]-y[n-1]| within t +- half, and the local RMS there (the click ratio in dB)."""
    i = int(round(t * R))
    a, b = max(1, i - int(half * R)), min(len(y) - 1, i + int(half * R))
    d = np.abs(np.diff(y[a:b], axis=0)).max()
    return float(20 * np.log10(d / (np.sqrt(np.mean(y[a:b] ** 2)) + 1e-30) + 1e-30))


def demod_delay(y, x, f):
    """Delay of y relative to x at f (seconds), from the analytic phase difference, and the IF of y."""
    n = min(len(y), len(x))
    y, x = y[:n], x[:n]

    def an(a):
        z = hilbert(bandpass(a, f * 0.9, f * 1.1), axis=0)[:, 0]
        return np.abs(z), np.unwrap(np.angle(z * np.exp(-2j * np.pi * f * np.arange(len(z)) / R)))
    ay, py = an(y)
    ax, px = an(x)
    d = -(py - px) / (2 * np.pi * f)
    ifv = f + np.gradient(py) / (2 * np.pi) * R
    return d, ifv, ay


def smooth(v, n):
    return np.convolve(v, np.ones(n) / n, "same")


def f2(v):
    return "--" if v is None else f"{v:.2f}"


def one(name):
    y, x, meta, _ = load(name)
    ev = meta.get("events") or []
    print(f"\n== {name}  ({meta['params'].get('Program', 'Reverb')}, {len(y) / R:.1f} s)")
    for e in (ev if len(ev) <= 6 else ev[:2] + ev[-1:]):
        print(f"   event t={e['t']:.2f} s  {e['param']} -> {e['value']}")
    if len(ev) > 6:
        print(f"   ({len(ev)} events; the first, second and last shown)")

    if name.startswith(("time_", "LFDecay_", "HFDecay_")):
        t = ev[0]["t"]
        before = (max(0.12, t - 1.2), t - 0.1)
        if before[1] - before[0] < 0.08:
            before = (max(0.12, t - 0.18), t - 0.02)
        after = (t + 0.03, min(len(y) / R, t + 2.0))
        print(f"   T60 before {before[0]:.2f}-{before[1]:.2f} s / after {after[0]:.2f}-{after[1]:.2f} s (broadband, L/R): "
              + " / ".join(f"{f2(local_t60(y[:, c], *before))} -> {f2(local_t60(y[:, c], *after))}" for c in (0, 1)))
        if name.startswith(("LFDecay_", "HFDecay_")):
            for fc in BANDS:
                lo, hi = fc / 1.414, min(fc * 1.414, 0.45 * R)
                lb = bandpass(y, lo, hi)[:, 0]
                print(f"   {fc:5d} Hz T60 before -> after L: {f2(local_t60(lb, *before))} -> {f2(local_t60(lb, *after))} s")
        print(f"   sample jump at the event (dB re local RMS): {jump_at(y, t):.1f}")
        # gain envelope across the event, in 10 ms blocks: the step at the event and the slope after it
        te, ee = envelope(y, 0.01)
        a = (te > t - 0.15) & (te < t - 0.02)
        b = (te > t + 0.02) & (te < t + 0.15)
        print(f"   envelope 20-150 ms before/after: {10 * np.log10(np.mean(ee[a] ** 2) / np.mean(ee[b] ** 2)):+.1f} dB (negative = falls)")

    elif name.startswith("predelay_"):
        t = ev[0]["t"]
        f = 1000.0
        d, ifv, a = demod_delay(y, x, f)
        i0, i1 = int((t - 0.8) * R), int((t + 0.8) * R)
        db = np.median(d[int((t - 0.6) * R):int((t - 0.05) * R)])
        da = np.median(d[int((t + 0.1) * R):int((t + 0.6) * R)])
        print(f"   1 kHz delay before/after: {db * 1e3:+.2f} / {da * 1e3:+.2f} ms (step {1e3 * (da - db):+.2f} ms); "
              f"IF 5-95% during the event {np.percentile(ifv[i0:i1], 5):.1f}-{np.percentile(ifv[i0:i1], 95):.1f} Hz")
        print(f"   sample jump at the event: {jump_at(y, t):.1f} dB")

    elif name.startswith("delay_L"):
        t = ev[0]["t"]
        _, ifv, _ = demod_delay(y, x, 1000.0)
        i0, i1 = int((t - 0.2) * R), int((t + 1.0) * R)
        base = np.median(ifv[int((t - 0.15) * R):int((t - 0.01) * R)])
        print(f"   1 kHz IF before {base:.1f} Hz, during 0.2 s after min/max {ifv[i0:i1].min():.1f} / {ifv[i0:i1].max():.1f} Hz; "
              f"sample jump at the event {jump_at(y, t):.1f} dB")

    elif name.startswith("mix_ramp"):
        t, e = envelope(y, 0.002)
        m = (t > 1.0) & (t < 5.1)
        step = np.abs(np.diff(10 * np.log10(e[m] ** 2 + 1e-30)))
        print(f"   2 ms envelope steps over the ramp: max {step.max():.1f} dB, 99th {np.percentile(step, 99):.1f} dB, median {np.median(step):.2f} dB")

    elif name == "power_toggle_midtail":
        for t in (1.0, 2.0):
            print(f"   level at t={t:.1f} s +-0.02 s: {10 * np.log10(np.mean(y[int((t - .02) * R):int((t + .02) * R)] ** 2)):.1f} dBFS")
        print(f"   sample jump at the off / on events: {jump_at(y, 1.0):.1f} / {jump_at(y, 2.0):.1f} dB")

    elif name == "program_switch_midtail":
        for t in (1.0, 2.0):
            print(f"   level at t={t:.1f} s +-0.02 s: {10 * np.log10(np.mean(y[int((t - .02) * R):int((t + .02) * R)] ** 2)):.1f} dBFS")
        print(f"   sample jump at the switch / back: {jump_at(y, 1.0):.1f} / {jump_at(y, 2.0):.1f} dB")


def main():
    names = sys.argv[1:] or [os.path.basename(p)[:-4] for p in sorted(os.listdir(CAP)) if p.endswith(".wav")]
    for n in names:
        one(n)


if __name__ == "__main__":
    main()
