"""A1: static structure of the UAD EMT 250 from the Noise-off impulse grid (docs/UAD250_test_plan.md §3 A).

    .venv/bin/python tools/uad250/analyze_a1.py [--only 'imp_T2.5*']

Per setting: first arrival per channel (ms, relative to the input impulse, plugin latency subtracted), early taps (R),
octave-band T60 (harness band_rt), recirculation loops (harness find_loops at the 48 kHz capture rate, so lengths are in
ms and converted to 26,373.6 Hz device samples for comparison with docs/EMT250_program.md §3). Writes
analysis/uad250/a1/a1.csv. Loop detection needs RT >= ~3 s (EMT250_program.md §3), so it is attempted at every setting
but only trusted where the ACF peaks clear 8 sigma in both channels.
"""
import argparse, csv, glob, json, os, re, sys
import numpy as np, soundfile as sf
from scipy.signal import find_peaks

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "tools", "harness"))
from harness.core import band_rt, find_loops  # noqa: E402

CAP = os.path.join(ROOT, "analysis", "uad250", "captures", "A_static")
OUT = os.path.join(ROOT, "analysis", "uad250", "a1")
FDEV = 19.2e6 / 728
OCT = [125, 250, 500, 1000, 2000, 4000, 8000]
LAT, IMP_AT = 523, 100


def analyse(path):
    y, r = sf.read(path); y = y.astype(float)
    m = re.match(r"imp_T([\d.]+)_LF([\d.]+)_HF(\S+)\.wav", os.path.basename(path))
    T, LF, HF = float(m[1]), float(m[2]), m[3]
    row = dict(T=T, LF=LF, HF=HF, peak_db=20 * np.log10(abs(y).max() + 1e-30))
    off = IMP_AT + LAT
    for c, ch in enumerate("LR"):
        a = np.abs(y[:, c]); th = 0.02 * a.max()
        row[f"first_{ch}_ms"] = (np.argmax(a > th) - off) / r * 1e3
    # early taps on R: local maxima above 25 % of the first 80 ms peak, at least 2 ms apart
    e = np.abs(y[:int(0.08 * r), 1]); pk, _ = find_peaks(e, height=0.25 * e.max(), distance=int(0.002 * r))
    row["early_R_ms"] = ";".join(f"{(p - off) / r * 1e3:.2f}" for p in pk[:8])
    row["early_R_gain"] = ";".join(f"{e[p] / e[pk[0]]:.2f}" for p in pk[:8])
    for f in OCT:
        row[f"T60_{f}"] = band_rt(y[:, 0], r, f / 1.414, f * 1.414)["T"]
        row[f"T60R_{f}"] = band_rt(y[:, 1], r, f / 1.414, f * 1.414)["T"]
    # loops: decay-flattened ACF from 0.3 s to where the tail reaches ~ -60 dB
    s0, s1 = int(0.3 * r), int(min(len(y), (0.3 + 1.2 * T) * r))
    loops, acfs = find_loops(y, s0, s1, T, float(r), n_sigma=8.0, max_lag=int(0.4 * r))
    row["loops_ms"] = ";".join(f"{p / r * 1e3:.2f}" for p in sorted(loops))
    row["loops_dev"] = ";".join(f"{p / r * FDEV:.1f}" for p in sorted(loops))
    return row


def main():
    ap = argparse.ArgumentParser(usage=__doc__); ap.add_argument("--only", default="*")
    a = ap.parse_args(); os.makedirs(OUT, exist_ok=True)
    rows = [analyse(p) for p in sorted(glob.glob(os.path.join(CAP, f"{a.only}.wav")))
            if os.path.basename(p).startswith("imp_")]
    rows.sort(key=lambda d: (d["T"], d["LF"], d["HF"]))
    if not rows: sys.exit("no files")
    with open(os.path.join(OUT, "a1.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, rows[0].keys()); w.writeheader(); w.writerows(rows)
    f2 = lambda v: "  -  " if v is None else f"{v:5.2f}"
    print("T60 (L) by octave band, LF x1.0, HF x0.50 [label RT vs 1 kHz]")
    print("  RT  " + " ".join(f"{f:>5d}" for f in OCT))
    for d in rows:
        if d["LF"] == 1.0 and d["HF"] == "0.50":
            print(f"{d['T']:5.1f} " + " ".join(f2(d[f"T60_{f}"]) for f in OCT))
    print("\nfirst arrival ms / early taps R / loops (ms) / loops (device samples), LF1.0 HF0.50")
    for d in rows:
        if d["LF"] == 1.0 and d["HF"] == "0.50":
            print(f"T{d['T']:<4} L {d['first_L_ms']:6.2f} R {d['first_R_ms']:6.2f}  early {d['early_R_ms']}  gains {d['early_R_gain']}"
                  f"  loops {d['loops_ms'] or '-'}  dev {d['loops_dev'] or '-'}")
    print(f"\n{len(rows)} settings -> {OUT}/a1.csv")


if __name__ == "__main__":
    main()
