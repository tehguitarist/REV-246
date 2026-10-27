"""Setting-level tally of a held-out S variant against the baseline (model_s_results.md §42).

    .venv/bin/python -I tools/xcorr_tally.py [variant_dir=s_xcorr]

Per setting and window: the mean over its captures, both channels and the bands of each group of the change in strict
null (variant minus baseline, dB; negative = better) from analysis/fullband/<variant>/*.json against analysis/fullband/s.
Prints the median per LF step and how many settings are worse / better by more than 0.5 dB (the strict rule is
"none worse"). The per-band, per-capture counts are in sml_fullband.py --stage report; this is the same rule pooled.
"""
import collections
import glob
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis", "fullband")
WINS = ["0-50", "50-150", "150-300", "300-600", "600-1000"]
GROUPS = {
    "<3.2k": ["20-40", "40-80", "80-160", "160-270", "270-400", "400-800", "800-1600", "1600-3200"],
    "3.2-6.4k": ["3200-6400"],
    "6.4-10k": ["6400-10000"],
}


def _db(e, y):
    return 10 * math.log10(max(e, 1e-30) / max(y, 1e-30))


def main():
    vdir = sys.argv[1] if len(sys.argv) > 1 else "s_xcorr"
    res = collections.defaultdict(lambda: collections.defaultdict(list))
    for p in sorted(glob.glob(os.path.join(ROOT, "s", "*.json"))):
        q = os.path.join(ROOT, vdir, os.path.basename(p))
        if not os.path.exists(q):
            continue
        S, V = json.load(open(p)), json.load(open(q))
        for grp, bands in GROUPS.items():
            for w in WINS:
                d = []
                for amp, c in S["caps"].items():
                    cv = V["caps"].get(amp)
                    if cv is None:
                        continue
                    for b in bands:
                        if w not in c["bands"][b] or w not in cv["bands"][b]:
                            continue
                        for ch in (0, 1):
                            d.append(_db(cv["bands"][b][w]["e"][ch], cv["bands"][b][w]["y"][ch])
                                     - _db(c["bands"][b][w]["e"][ch], c["bands"][b][w]["y"][ch]))
                if d:
                    res[(grp, w)][S["lf"]].append(float(np.mean(d)))
    for grp in GROUPS:
        print(f"{grp}: median change per setting (dB, negative = better) by LF; settings worse / better by > 0.5 dB / total")
        for w in WINS:
            allv, row = [], []
            for lf in sorted(res[(grp, w)]):
                v = res[(grp, w)][lf]
                allv += v
                row.append(f"LF{lf}:{np.median(v):+.2f}")
            a = np.array(allv)
            print(f"  {w:9s} {' '.join(row)} | worse {int(np.sum(a > 0.5))} better {int(np.sum(a < -0.5))} / {len(a)}")


if __name__ == "__main__":
    main()
