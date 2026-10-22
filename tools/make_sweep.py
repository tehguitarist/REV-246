"""Stimulus files for sweep-capturing a convolution plugin (the Altiverb 246): exponential sine sweep (Farina),
its inverse filter, and a single-sample impulse. Used by tools/altiverb_compare.py. Notes: docs/library_notes.md §6.

    .venv/bin/python tools/make_sweep.py [--out sweep_request]

Writes mono 32-bit float WAVs and inverse_filter_<rate>k.npy (not committed: regenerate with this script).
  EMT246_logsweep_<48|96>k_20-20k_10s_m<6|12|18>dBFS.wav   1 s silence + 10 s sweep (20 Hz-20 kHz) + 6 s silence
  EMT246_logsweep_48k_20-20k_10s_m<12|18>dBFS_15s-tail.wav  the same with 15 s of silence after (for L's long tails)
  EMT246_impulse_48k_m6dBFS.wav                             1 s silence, one sample, 4 s silence
The inverse is built from the -6 dBFS sweep and is the same for every level (output 12 dB lower at -18 dBFS).
"""
import argparse, os
import numpy as np
import soundfile as sf
from scipy.signal import fftconvolve

F1, F2, T = 20.0, 20000.0, 10.0
LEAD, TAIL, PEAK_DB = 1.0, 6.0, -6.0
FADE = 0.05  # s, raised-cosine at both ends


def sweep_and_inverse(fs):
    """Unit-peak-scaled (-6 dBFS) exponential sweep and its inverse filter at rate fs."""
    t = np.arange(int(T * fs)) / fs
    L = T / np.log(F2 / F1)
    s = np.sin(2 * np.pi * F1 * L * (np.exp(t / L) - 1))
    nf = int(FADE * fs)
    w = 0.5 * (1 - np.cos(np.pi * np.arange(nf) / nf))
    s[:nf] *= w
    s[-nf:] *= w[::-1]
    s *= 10 ** (PEAK_DB / 20)
    inv = s[::-1] * np.exp(-t / L)  # -6 dB/oct amplitude slope: sweep * inverse is a flat band-limited pulse
    return s, inv


def stimulus(fs, level_db, tail_s=TAIL):
    s, _ = sweep_and_inverse(fs)
    return np.concatenate([np.zeros(int(LEAD * fs)), s * 10 ** ((level_db - PEAK_DB) / 20), np.zeros(int(tail_s * fs))])


def write(path, x, fs):
    sf.write(path, x.astype(np.float32), fs, subtype="FLOAT")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="sweep_request")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    for fs in (48000, 96000):
        k = fs // 1000
        s, inv = sweep_and_inverse(fs)
        np.save(f"{out}/inverse_filter_{k}k.npy", inv)
        for db in (-6, -12, -18):
            write(f"{out}/EMT246_logsweep_{k}k_20-20k_10s_m{-db}dBFS.wav", stimulus(fs, db), fs)
        # self-test: sweep * inverse is a flat pulse at the expected sample
        d = fftconvolve(stimulus(fs, PEAK_DB), inv)
        i = int(np.argmax(np.abs(d)))
        assert i == int(LEAD * fs) + len(inv) - 1, i
        P = np.abs(np.fft.rfft(d[i - 2000:i + 2000]))
        f = np.fft.rfftfreq(4000, 1 / fs)
        m = (f > 60) & (f < 16000)
        print(f"{fs} Hz: pulse at the expected sample; flat to {20 * np.log10(P[m].max() / P[m].min()):.2f} dB (60 Hz-16 kHz)")
    for db in (-12, -18):
        write(f"{out}/EMT246_logsweep_48k_20-20k_10s_m{-db}dBFS_15s-tail.wav", stimulus(48000, db, 15.0), 48000)
    imp = np.zeros(5 * 48000)
    imp[48000] = 10 ** (-6 / 20)
    write(f"{out}/EMT246_impulse_48k_m6dBFS.wav", imp, 48000)


if __name__ == "__main__":
    main()
