"""Held-out test of a per-step loop-gain table at LATE 0.2 (NonLin_program 10.4). Train = every other clean LATE 0.2
capture (passes_0.2.json ratios, median of L and R pooled, captures whose 3 ratios are all within 10 % of the pooled
median), test = the rest, scored with the bench's own per-capture alignment (o, clock, front shift; gains by least squares):
window nulls nominal vs table. Also every LATE 0.2 capture with a trimmed start (strict rule) and LATE 0.4 controls.
    heldout.py  -> analysis/nl/late02/heldout.json"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from late02 import Cap, keys_for, kstr, nl, nlt, nm, defects   # noqa (sets NL_BK, EMT_KERNEL)
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ps = json.load(open(os.path.join(HERE, "passes_0.2.json")))
keys = keys_for(0.2)
train, test = keys[0::2], keys[1::2]
R = np.array([ps["captures"][kstr(k)][str(c)]["ratio"] for k in train for c in (0, 1)])
med = np.median(R, 0)
ok = np.all(np.abs(R / med - 1) < 0.1, 1)
TABLE = {P: float(v) for P, v in zip(nlt.LOOPS, np.median(R[ok], 0))}
print("train", len(train), "captures,", int(ok.sum()), "of", len(R), "channel rows kept; table", {P: round(v, 5) for P, v in TABLE.items()},
      "nominal", {P: round(nlt.loop_gain(P, 0.2), 5) for P in nlt.LOOPS}, flush=True)
_orig = nlt.loop_gain


def use(table):
    nlt.loop_gain = (lambda P, late: TABLE[P] if late == 0.2 else _orig(P, late)) if table else _orig
    nm.spectra.cache_clear()


def score(k):
    c = Cap(k)
    m0 = _nulls(c)
    use(True)
    m1 = _nulls(c)
    use(False)
    return dict(nominal=m0, table=m1, floor={w: v[2] for w, v in c.wins.items()})


def _nulls(c):
    m = nm.render(c.key, c.n, c.b["o"], c.b["front_shift"])
    g = [(c.y[:, ch] @ m[:, ch]) / (m[:, ch] @ m[:, ch]) for ch in (0, 1)]
    m = m * np.array(g)
    return {w: [float(10 * np.log10(np.sum((c.y[a:b, ch] - m[a:b, ch]) ** 2) / np.sum(c.y[a:b, ch] ** 2))) for ch in (0, 1)]
            for w, (a, b, _) in c.wins.items()}


trimmed = []
for k in sorted(nl.table()):
    if k[2] != 0.2 or k in defects.NL_DEFECTS or k in keys:
        continue
    b = json.load(open(os.path.join(nm.BENCH, kstr(k) + ".json")))
    trimmed.append(k)
ctrl = keys_for(0.4, 8)
out = dict(table={str(P): v for P, v in TABLE.items()}, nominal={str(P): _orig(P, 0.2) for P in nlt.LOOPS},
           train=[kstr(k) for k in train], test={}, train_scores={}, other={}, control={})
for name, ks in (("test", test), ("train_scores", train), ("other", trimmed), ("control", ctrl)):
    for k in ks:
        out[name][kstr(k)] = s = score(k)
        print(name, kstr(k), "R nominal", " ".join(f"{v[1]:.1f}" for v in s["nominal"].values()), "| table",
              " ".join(f"{v[1]:.1f}" for v in s["table"].values()), flush=True)
json.dump(out, open(os.path.join(HERE, "heldout.json"), "w"))
