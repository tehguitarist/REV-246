"""lowband_shelf.py --stage hfbb by (HF, RT) cell: the capped-law candidate (no ramp, no HF table) at eps = 0 against the
baseline, each setting's own best eps and the cost there. model_s_results §25.15.

    .venv/bin/python tools/loop_read/hfbb_report.py"""
import collections, glob, json, os
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = [json.load(open(p)) for p in glob.glob(os.path.join(ROOT, "analysis/lowband/chain3_shelf_hfbb/*.json"))]
C = collections.defaultdict(list)
for d in D:
    C[(d["hf"], d["rt"])].append(d)
print("cell        | candidate eps=0 re baseline (dB, from 150 ms / 50 ms), LF 0.1 0.4 1.0 2.0 | own best eps | at own best (150 ms)")
for k in sorted(C):
    u = sorted(C[k], key=lambda d: d["lf"])
    g = [np.array(d["grid"]) for d in u]
    z = [x[np.argmin(abs(x[:, 0]))] for x in g]
    best = [x[np.argmin(x[:, 2])] for x in g]
    print(f"HF {k[0]} RT {k[1]:<3} | " + " ".join(f"{10 * np.log10(r[2] / d['base'][1]):+.2f}/{10 * np.log10(r[1] / d['base'][0]):+.2f}" for r, d in zip(z, u)) +
          " | " + " ".join(f"{b[0]:+.3f}" for b in best) + " | " + " ".join(f"{10 * np.log10(b[2] / d['base'][1]):+.2f}" for b, d in zip(best, u)))
