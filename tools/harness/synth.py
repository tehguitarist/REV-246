"""Build a synthetic 'library' from one real EMT 246 IR, with known settings changes and injected faults.

Used by tests/test_harness.py to check that every stage-A–C test fires. The setting changes are
approximations (envelope reshaping of the real IR), good enough to exercise the measurements.

    .venv/bin/python tools/harness/synth.py <source.wav> <out_dir>
"""
import os
import sys

import numpy as np
import soundfile as sf
from scipy.signal import butter, resample_poly, sosfiltfilt

BASE = dict(rt=1.0, lf=0.9, hf=0.9, d=30, a=10)  # arbitrary test labels (not the CLA IR's real settings)


def name(rt=None, lf=None, hf=None, d=None, a=None):
    s = dict(BASE, **{k: v for k, v in dict(rt=rt, lf=lf, hf=hf, d=d, a=a).items() if v is not None})
    return f"EMT246_S_RT{s['rt']}_LF{s['lf']}_HF{s['hf']}_D{s['d']}_A{s['a']}"


def onset(x):
    a = np.abs(x).max(axis=1)
    return int(np.argmax(a > a.max() * 0.01))


def reshape_decay(x, fs, t_start, factor_db_per_s, band=None):
    """Multiply the signal after t_start by an extra exponential decay (optionally only in a band)."""
    t = np.maximum(np.arange(len(x)) / fs - t_start, 0)
    env = 10 ** (factor_db_per_s * t / 20)[:, None]
    if band is None:
        return x * env
    kind, fc = band
    sos = butter(4, fc, btype="low" if kind == "lf" else "high", fs=fs, output="sos")
    part = sosfiltfilt(sos, x, axis=0)
    return x - part + part * env


def build(src, out):
    x, fs = sf.read(src, always_2d=True)
    on = onset(x)
    t_rev = on / fs + 0.019          # after the cluster (+ ringing), before the reverb (+22 ms)
    rt0 = BASE["rt"]
    files = {}

    def put(rel, y, sr=fs, labels=None):
        p = os.path.join(out, rel + ".wav")
        os.makedirs(os.path.dirname(p), exist_ok=True)
        sf.write(p, np.clip(y, -1, 1) if y.ndim == 1 else y, sr, subtype="PCM_24")
        files[rel] = p

    put(name(), x)
    for rt in (0.5, 0.7):  # shorter RT: extra decay of 60/rt - 60/rt0 dB/s
        put(name(rt=rt), reshape_decay(x, fs, t_rev, -(60 / rt - 60 / rt0)))
    put(name(hf=0.5), reshape_decay(x, fs, t_rev, -(60 / (0.5 * rt0) - 60 / (0.9 * rt0)), band=("hf", 5000)))
    put(name(lf=0.5), reshape_decay(x, fs, t_rev, -(60 / (0.5 * rt0) - 60 / (0.9 * rt0)), band=("lf", 150)))
    cl = slice(on - int(0.001 * fs), on + int(0.019 * fs))
    y = x.copy(); y[cl] = 0
    put(name(a=0), y)
    y = x.copy(); y[cl] *= 0.5
    put(name(a=5), y)
    for extra in (10, 50, 100):
        put(name(d=30 + extra), np.vstack([np.zeros((int(extra / 1000 * fs), 2)), x]))
    put("48k/" + name(), resample_poly(x, 1, 2, axis=0), 48000)
    put("44k/" + name(), resample_poly(x, 147, 320, axis=0), 44100)
    rng = np.random.default_rng(0)
    put("repeat/" + name(), x + rng.standard_normal(x.shape) * 10 ** (-97 / 20))
    put("qc/" + name() + "_swapped", x[:, ::-1])
    put("qc/" + name() + "_invR", x * np.array([1, -1]))
    put("qc/" + name() + "_louder", x * 10 ** (1.2 / 20))
    fade = np.ones(len(x)); n = int(0.8 * fs); fade[-n:] = np.linspace(1, 0, n) ** 2
    put("qc/" + name() + "_faded", x * fade[:, None])
    y = x.copy(); y[on - int(0.006 * fs):on - int(0.005 * fs)] += rng.standard_normal((int(0.001 * fs), 2)) * 10 ** (-50 / 20)
    put("qc/" + name() + "_preecho", y)
    put("qc/" + name() + "_mono", x[:, 1])
    put("qc/" + name() + "_clipped", np.clip(x * 1.6, -1, 1))
    put("qc/" + name() + "_trimmed", x[on - int(0.0001 * fs):])
    put("qc/" + name() + "_skew150ppm", resample_poly(x, 20003, 20000, axis=0))  # capture clock 150 ppm fast
    put("qc/" + name(rt=1.0) + "_mislabelled", reshape_decay(x, fs, t_rev, -(60 / 0.7 - 60 / rt0)))
    put("qc/" + name(rt=2.0) + "_dup_other_label", x)
    put("qc/random_unparseable_name", x)
    return files


if __name__ == "__main__":
    for k, v in build(sys.argv[1], sys.argv[2]).items():
        print(v)
