"""Which plugin character settings sit closest to each IR? FLOAT / 16-BIT x XFMR off / on, per A/B example.

    .venv/bin/python -I tools/ab_variants.py   -> analysis/ab_samples/variants.md

For each example of tools/ab_samples.py the plugin wet is rendered four ways (sixteen = 0/1, xfmr = 0/1) and nulled
(best gain, re the IR wet) against the IR's wet after the 2-number clock fit of tools/ab_analysis.py, in three bands.
A static IR has no quantisation noise or A/D clip, so 16-BIT can only move a null away from it; XFMR's low-band
cascade was fitted to the CLA capture (docs/model_s_results.md 37), so it should help the first IR most.
"""
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import butter, sosfiltfilt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
import ab_analysis as A
import ab_samples as S

BANDS = {"32-100 Hz": (32, 100), "100 Hz-4 kHz": (100, 4000), "4-12 kHz": (4000, 12500)}


def bandpass(x, lo, hi):
    return sosfiltfilt(butter(4, [lo, hi], btype="band", fs=S.RATE, output="sos"), x, axis=0)


def main():
    lines = ["# Plugin character settings against each IR", "", "Null (dB, re the IR wet; lower = closer) per band after the best gain, "
             "against the IR wet, whose 2-number clock fit is redone for each variant. See `tools/ab_variants.py`.", "",
             "| example | 16-bit | XFMR | " + " | ".join(BANDS) + " |", "|---|---|---|---|---|---|"]
    for ex in S.EXAMPLES:
        n = ex["name"]
        z = np.load(S.WORK / f"{n}_wets.npz")
        mono = sf.read(S.WORK / f"{n}_mono.wav")[0]
        pl0 = z["plugin"].astype(float)
        n_out = len(pl0)
        for sixteen in (0, 1):
            for xfmr in (0, 1):
                e = dict(ex, args=dict(ex["args"], sixteen=sixteen, xfmr=xfmr))
                w, _ = S.render_plugin(e, S.WORK / f"{n}_mono.wav", n_out, f"{n}_v{sixteen}{xfmr}")
                ir, _ = A.clock_fit(mono, z["kernel"].astype(float), float(z["gain"]), w)  # refit per variant: each gets its own 2 numbers
                cells = []
                for lo, hi in BANDS.values():
                    cells.append(f"{A.null_db(bandpass(ir, lo, hi), bandpass(w, lo, hi))[0]:.1f}")
                lines.append(f"| {n} | {'on' if sixteen else 'off'} | {'on' if xfmr else 'off'} | " + " | ".join(cells) + " |")
        print(n, flush=True)
    (S.WORK / "variants.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
