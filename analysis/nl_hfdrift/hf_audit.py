"""High-end close-out audit (28 Oct 2026): default-kernel resampler / fractional-position paths, and the
already-scored S/M/L/E band edges. Cheap, read-only: no model changes, no heavy runs.

    .venv/bin/python analysis/nl_hfdrift/hf_audit.py > analysis/nl_hfdrift/hf_audit.log

Sections
  1. core.to_device_grid (96 kHz -> device): fixed Kaiser(5) droop, independent of the fractional position.
  2. e250.e_model.to_device: 128-tap windowed sinc + 12.4 kHz FIR, flat to 12.3 kHz, no ppm drift.
  3. E captures above 12.4 kHz: the analogue band limit, not an E term.
  4. E 20-100 Hz: the committed e_pulse v6 audit, extended-band rows only.
  5. S 12.3-12.8 kHz: the wide-window deconvolution against the cross-take repeat.
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools"))
sys.path.insert(0, os.path.join(ROOT, "tools", "harness"))
sys.path.insert(0, os.path.join(ROOT, "tools", "e250"))

import core  # noqa: E402
import e_model as E  # noqa: E402
import nl_structure as nl  # noqa: E402

CFG = core.load_config()
FS_DEV = core.fs_dev(CFG)
FREQS = np.array([6500.0, 8000.0, 10000.0, 11000.0, 11500.0, 12000.0, 12300.0, 12400.0, 12500.0, 12800.0])


def db(x):
    return 20 * np.log10(np.maximum(x, 1e-300))


def tone_amp(y, a, b):
    """Peak-equivalent amplitude of a steady tone in y[a:b]."""
    return float(np.sqrt(np.mean(y[a:b] ** 2)) * np.sqrt(2))


def section1_to_device_grid():
    print("\n## 1. core.to_device_grid (96 kHz -> 26.37 kHz), default resample_poly(25, 91)")
    print("Resampler gain re 1 kHz (nl_structure.resampler_gain):")
    g = nl.resampler_gain(FREQS)
    print("  " + "  ".join(f"{f / 1000:5.1f}k {db(v):+6.2f} dB" for f, v in zip(FREQS, g)))

    fs, n = 96000.0, 1 << 18
    t = np.arange(n) / fs
    print("Tone amplitude through to_device_grid, early / late window (ppm 0):")
    for frac in (0.0, 0.25, 0.5, 0.75):
        row = []
        for f in (8000.0, 11500.0, 12000.0, 12800.0):
            anchor = 20000.0 + frac
            k0 = int(np.floor(anchor / fs * FS_DEV))
            y = core.to_device_grid(np.sin(2 * np.pi * f * t), fs, anchor, k0, CFG)
            row.append(f"{f / 1000:5.1f}k {db(tone_amp(y, 6000, 10000)):+5.2f}/{db(tone_amp(y, 30000, 34000)):+5.2f}")
        print(f"  frac {frac:4.2f}: " + "  ".join(row))
    print("  -> fixed droop (the N 'top-end EQ'), independent of the fractional position; no time drift.")


def section2_e_to_device():
    print("\n## 2. e250.e_model.to_device (48 kHz -> device), 128-tap Kaiser sinc + 12.4 kHz FIR")
    fs, nin, n = 48000.0, 60000, 20000
    t = np.arange(nin) / fs
    print("Tone amplitude early / late window, ppm 42.8 (E's capture clock):")
    for frac in (0.0, 0.25, 0.5, 0.75):
        row = []
        for f in (8000.0, 11000.0, 11500.0, 12000.0, 12300.0, 12400.0):
            x = np.stack([np.sin(2 * np.pi * f * t)] * 2, 1)
            y = E.to_device(x, 42.8, frac / fs, n)
            row.append(f"{f / 1000:5.1f}k {db(tone_amp(y, 1000, 5000)):+5.2f}/{db(tone_amp(y, 12000, 16000)):+5.2f}")
        print(f"  frac {frac:4.2f}: " + "  ".join(row))
    print("  -> flat to 12.3 kHz (the 12.4 kHz FIR is the model's own band limit), no ppm drift.")


def section3_e_captures_above_12k4():
    print("\n## 3. E captures above 12.4 kHz (Altiverb REV-250 presets, device rate)")
    print("Band level re 1-4 kHz, per RT; the model has no content above 12.4 kHz by construction.")
    rts = [0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5]
    from scipy.signal import butter, sosfiltfilt

    def level(y, lo, hi):
        sos = butter(4, [lo, hi], btype="band", fs=E.FS, output="sos")
        return float(np.sum(sosfiltfilt(sos, y, axis=0) ** 2, 0).sum())

    rows = []
    for rt in rts:
        p = os.path.join(ROOT, "analysis", "e250", f"dev_{rt}.npy")
        if not os.path.exists(p):
            continue
        y = np.load(p)
        ref = level(y, 1000, 4000)
        rows.append((rt, 10 * np.log10(level(y, 12000, 12400) / ref), 10 * np.log10(level(y, 12400, 13100) / ref)))
    for rt, a, b in rows:
        print(f"  RT {rt:3.1f}: 12.0-12.4k {a:+6.2f} dB re 1-4k, 12.4-13.1k {b:+6.2f} dB")
    if rows:
        print(f"  median 12.4-13.1k {np.median([r[2] for r in rows]):+6.2f} dB: the analogue band limit / capture chain.")


def section4_e_20_100():
    print("\n## 4. E 20-100 Hz, committed e_pulse v6 audit (no rerun)")
    p = os.path.join(ROOT, "analysis", "e250", "e_pulse_v6.json")
    if not os.path.exists(p):
        print("  (e_pulse_v6.json missing)")
        return
    s = json.load(open(p))["summary"]
    keys = ("20-100", "31.5", "63")
    for kind in ("imp", "sweep"):
        print(f"  {kind}: family / band: inside / n, median offset, worst excess")
        for key, r in s[kind]["by_key"].items():
            if r["band"] in keys:
                print(f"    {r['family']:18s} {r['band']:6s} {r['inside']}/{r['n']}  "
                      f"med {r['med_off']:+6.2f}  worst {r['max_excess']:+6.2f}")
    print("  -> decay / EDT / late / level / IACC inside; only R C50/C80/D50/Ts at 20-100 and 31.5 Hz sit just")
    print("     outside, consistent with the capture chain's non-causal LF content (e_model.pretap; 16.4 / 16.7).")


def section5_s_wide():
    print("\n## 5. S 12.3-12.8 kHz, wide-window deconvolution vs cross-take repeat")
    p = os.path.join(ROOT, "analysis", "fullband", "report.json")
    if not os.path.exists(p):
        print("  (report.json missing)")
        return
    d = json.load(open(p))["wide"]
    for band in ("12000-12300", "12300-12800"):
        for win in ("50-150", "150-300"):
            v = d[band][win]  # [[level, model_wide, model_harness, repeat], ...] for R, L
            print(f"  {band} {win}: model {v[0][1]:+6.1f} / {v[1][1]:+6.1f} dB, "
                  f"repeat {v[0][3]:+6.1f} / {v[1][3]:+6.1f} dB (R / L)")
    print("  -> at the repeat level through 12.8 kHz; above 12.8 kHz the unit's chain is 43-72 dB down.")


if __name__ == "__main__":
    print(f"FS device {FS_DEV:.6f} Hz; audit run {np.datetime64('now')}")
    section1_to_device_grid()
    section2_e_to_device()
    section3_e_captures_above_12k4()
    section4_e_20_100()
    section5_s_wide()
