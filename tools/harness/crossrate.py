"""Cross-rate comparison of same-setting files (e.g. 48 kHz vs 96 kHz) in the capture domain.

The library-level pair test in library.py compares short device-rate signatures after an onset anchor.
For the Tomkin library that failed for ambiguous onsets (reversed Non-Lin, low-AMPL files), so this tool
compares the two files directly: the higher-rate file is decimated to the lower rate, aligned to a
fractional lag, gain-fitted, and the residual is reported per time window relative to the energy in
that window. It also reports the correlation of the out-of-band (14.5-22 kHz) tail content: independent
captures have independent noise there (correlation near 0); one capture rendered at two rates shares it.

  .venv/bin/python tools/harness/crossrate.py "irs/EMT 246 captures" --out analysis/library/tomkin/crossrate.csv
Expects <root>/<N>kHz/<same relative path> for each rate; N is read from the folder name.
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
import soundfile as sf
from concurrent.futures import ProcessPoolExecutor
from scipy.signal import butter, resample_poly, sosfiltfilt

WINDOWS = [(0.0, 0.3), (0.3, 1.0), (1.0, 3.0)]


def _align(a, b):
    n = min(len(a), len(b)); a, b = a[:n], b[:n]
    N = 1 << int(np.ceil(np.log2(2 * n)))
    A, B = np.fft.rfft(a, N, axis=0), np.fft.rfft(b, N, axis=0)
    xc = np.fft.irfft(np.sum(A * np.conj(B), axis=1), N)
    k = int(np.argmax(xc)); k = k if k < N // 2 else k - N
    f = np.fft.rfftfreq(N); m = min(n, 48000)
    best = None
    for tau in np.arange(k - 1, k + 1.001, 0.05):
        bs = np.fft.irfft(B * np.exp(-2j * np.pi * f * tau)[:, None], N, axis=0)[:n]
        r = np.sum((a[:m] - bs[:m]) ** 2)
        if best is None or r < best[0]:
            best = (r, tau, bs)
    _, tau, bs = best
    g = np.sum(a * bs) / np.sum(bs * bs)
    return a, g * bs, tau, g


def compare(job):
    root, rel, lo, hi = job
    try:
        a, fs_a = sf.read(os.path.join(root, f"{lo}kHz", rel), always_2d=True)
        b, fs_b = sf.read(os.path.join(root, f"{hi}kHz", rel), always_2d=True)
        a, b = a[:, :2], b[:, :2]
        q = fs_b // fs_a
        if q > 1:
            b = resample_poly(b, 1, q, axis=0)
        a, bs, tau, g = _align(a, b)
        row = {"relpath": rel, "fs_lo": fs_a, "fs_hi": fs_b, "lag_samples": round(float(tau), 3), "gain": round(float(g), 5)}
        for t0, t1 in WINDOWS:
            s = slice(int(t0 * fs_a), int(t1 * fs_a))
            ea = np.sum(a[s] ** 2)
            row[f"res_db_{t0:g}_{t1:g}s"] = round(float(10 * np.log10(np.sum((a[s] - bs[s]) ** 2) / ea + 1e-30)), 1) if ea > 0 else None
        sos = butter(6, [14500, 22000], btype="band", fs=fs_a, output="sos")
        s = slice(int(0.5 * fs_a), None)
        if len(a[s]) > 100:   # sosfiltfilt needs more than its pad length; files under about 0.5 s have no tail
            ha, hb = sosfiltfilt(sos, a[s], axis=0), sosfiltfilt(sos, bs[s], axis=0)
            row["hf_tail_corr"] = round(float(np.sum(ha * hb) / np.sqrt(np.sum(ha ** 2) * np.sum(hb ** 2) + 1e-300)), 3)
            row["hf_tail_dbfs"] = round(float(10 * np.log10(np.mean(ha ** 2) + 1e-30)), 1)
        return row
    except Exception as e:  # noqa: BLE001
        return {"relpath": rel, "error": repr(e)[:80]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root"); ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    rates = sorted(int(m.group(1)) for d in os.listdir(a.root) if (m := re.fullmatch(r"(\d+)\s*kHz", d, re.I)))
    lo, hi = rates[0], rates[-1]
    rels = sorted(os.path.relpath(os.path.join(d, f), os.path.join(a.root, f"{lo}kHz"))
                  for d, _, fs in os.walk(os.path.join(a.root, f"{lo}kHz")) for f in fs
                  if f.lower().endswith(".wav") and os.path.exists(os.path.join(a.root, f"{hi}kHz", os.path.relpath(os.path.join(d, f), os.path.join(a.root, f"{lo}kHz")))))
    rels = rels[: a.limit]
    with ProcessPoolExecutor(a.workers) as ex:
        rows = list(ex.map(compare, [(a.root, r, lo, hi) for r in rels], chunksize=4))
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(f"{len(rows)} pairs ({lo} vs {hi} kHz) -> {a.out}")


if __name__ == "__main__":
    main()
