"""Report for tools/ml_fix_diag.py --stage sw3 (item 1 smoke: the audit's M / L metrics with the plugin measured through
ml_variance's sweep at 3 start phases): analysis/ml_fix/sw3_report.log, sw3.json.

    .venv/bin/python tools/ml_fix_sw3_report.py

Columns, medians (R; |ratio - 1| %, EDC dev dB, nulls dB):
  SW   sweep-measured plugin vs capture, 3 phases (0 / 120 / 240 deg)   AUD  the audit (impulse, phase 0; plugin_pulse)
  SC   static plugin (impulse, mod off) vs capture                      MC   pulse_check's static model vs capture
  CC   capture vs capture (every pair; the tolerance)                   SPP  sweep phase vs sweep phase (x ~1.22 / +1.8 dB
       against random pairs: 120 deg spacing)
  in-tol: share of captures (x phases) whose |error| is within the group's largest capture-pair |difference| (floor 3 % for
  the decay ratios, as plugin_pulse._within; no floor for EDC dev), SW / AUD / MC; and within the group's largest SPP.
"""
import glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ.setdefault("EMT_KERNEL", "chain3")
import numpy as np  # noqa: E402
import pulse_check as pc  # noqa: E402
from defects import ML_SLIPPED  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_fix")
REF = os.path.join(ROOT, "analysis", "pulse", "ml")
AUD = os.path.join(ROOT, "analysis", "plugin_pulse", "ml")
MET = [f"oct_t{f}" for f in pc.OCT] + [f"oct_edt{f}" for f in (125, 250, 500, 1000, 2000)] + \
      [f"oct_late{f}" for f in (125, 250, 500, 1000)] + [f"oct_dev{f}" for f in (125, 250, 500, 1000)] + \
      ["rv_mix_ms", "ned_def", "ned"]
NUL = ["first_pass", "recirc1", "recirc2", "recirc3", "0.3-0.6", "0.6-1.0", "1.0-2.0", "2.0-4.0"]
COLS = ["SW", "AUD", "SC", "MC", "CC", "SPP"]


def ned_def(m):
    try:
        a, b = (np.array(m["rv"][x]["ned1"], float) for x in ("a", "b"))
    except (KeyError, TypeError):
        return None
    n = min(len(a), len(b))
    if n < 10:
        return None
    da, db = np.sum(np.maximum(0, 1 - a[:n])), np.sum(np.maximum(0, 1 - b[:n]))
    return float(da / db - 1) if db > 0 else None


def summ(m):
    s = pc.summarise(m)
    s["ned_def"] = ned_def(m)
    return s


def sc(k):
    return 1.0 if k.startswith(("oct_dev", "ned")) and k != "ned_def" else 100.0


def med(x, f=abs):
    x = [f(v) for v in x if v is not None and np.isfinite(v)]
    return float(np.median(x)) if x else None


def load():
    G = []
    for fn in sorted(glob.glob(os.path.join(OUT, "sw3", "*.json"))):
        d = json.load(open(fn))
        prog, key = d["prog"], d["key"]
        bad = lambda a: (prog, key, a) in ML_SLIPPED
        ref = json.load(open(os.path.join(REF, os.path.basename(fn))))
        aud = json.load(open(os.path.join(AUD, os.path.basename(fn))))
        g = dict(prog=prog, key=key, rt=d["rt"], lf=d["lf"], v={c: [] for c in COLS}, n={c: [] for c in COLS}, per={})
        for r in d["rows"]:
            if bad(r["amp"]):
                continue
            for v in ("s0", "s1", "s2"):
                g["v"]["SW"].append(r["vs_cap"][v]); g["n"]["SW"].append(r["null_cap"][v])
            g["v"]["SC"].append(r["vs_cap"]["static"]); g["n"]["SC"].append(r["null_cap"]["static"])
            for k, x in r["pairs"].items():
                g["v"]["SPP"].append(x); g["n"]["SPP"].append(r["null_pairs"][k])
        for r in ref["rows"]:
            if not bad(r["amp"]):
                g["v"]["MC"].append(summ(r["m"])); g["n"]["MC"].append(r["null"])
        for p in ref["pairs"]:
            if not any(bad(a) for a in p["pair"]):
                g["v"]["CC"].append(summ(p["m"])); g["n"]["CC"].append(p["null"])
        for r in aud["rows"]:
            if not bad(r["amp"]):
                g["v"]["AUD"].append(summ(r["m"])); g["n"]["AUD"].append(r["null"])
        G.append(g)
    return G


def within(gs, m, col, ref="CC"):
    ok = n = 0
    floor = 0.0 if m.startswith(("oct_dev", "ned")) and m != "ned_def" else 0.03
    for g in gs:
        rep = [abs(s.get(m)) for s in g["v"][ref] if s.get(m) is not None]
        if not rep:
            continue
        tol = max(max(rep), floor)
        x = [abs(s.get(m)) for s in g["v"][col] if s.get(m) is not None]
        ok += sum(v <= tol for v in x)
        n += len(x)
    return ok, n


def main():
    G = load()
    out, res = [], {}
    for prog in ("M", "L"):
        gs = [g for g in G if g["prog"] == prog]
        if not gs:
            continue
        r = res.setdefault(prog, dict(groups=[g["key"] for g in gs], metrics={}, nulls={}, per_group={}))
        out.append(f"\n## {prog}: {len(gs)} groups: {', '.join(g['key'] for g in gs)}")
        out.append("metric | " + " | ".join(COLS) + " | in CC tol SW / AUD / MC | in SPP tol SW / AUD / MC")
        for m in MET:
            v = {c: med([s.get(m) for g in gs for s in g["v"][c]]) for c in COLS}
            w = {c: within(gs, m, c) for c in ("SW", "AUD", "MC")}
            w2 = {c: within(gs, m, c, "SPP") for c in ("SW", "AUD", "MC")}
            k = sc(m)
            r["metrics"][m] = dict(median={c: None if v[c] is None else v[c] * k for c in COLS}, within_CC=w, within_SPP=w2)
            f = lambda x: "-" if x is None else f"{x * k:.2f}"
            fw = lambda d: " / ".join(f"{a}/{b}" for a, b in d.values())
            out.append(f"{m} | " + " | ".join(f(v[c]) for c in COLS) + f" | {fw(w)} | {fw(w2)}")
        out.append("null | " + " | ".join(COLS))
        for wn in NUL:
            v = {c: med([x.get(wn) for g in gs for x in g["n"][c]], f=lambda z: z) for c in COLS}
            r["nulls"][wn] = v
            out.append(f"{wn} | " + " | ".join("-" if v[c] is None else f"{v[c]:.1f}" for c in COLS))
        fb = 125 if prog == "M" else 250
        out.append(f"\n### {prog} per group: {fb} Hz decay, 500 Hz decay, EDT 500 (|%| SW / AUD / MC / CC / SPP); recirc1 SW / AUD / CC / SPP; ned_def SW / AUD / CC / SPP")
        for g in gs:
            q = lambda c, m: med([s.get(m) for s in g["v"][c]])
            f = lambda x: "-" if x is None else f"{100 * x:.1f}"
            cells = []
            for m in (f"oct_t{fb}", "oct_t500", "oct_edt500"):
                cells.append(" / ".join(f(q(c, m)) for c in ("SW", "AUD", "MC", "CC", "SPP")))
            n1 = [med([x.get("recirc1") for x in g["n"][c]], f=lambda z: z) for c in ("SW", "AUD", "CC", "SPP")]
            nd = " / ".join(f(q(c, "ned_def")) for c in ("SW", "AUD", "CC", "SPP"))
            out.append(f"- {g['key']}: " + "; ".join(cells) + "; rc1 " + " / ".join("-" if x is None else f"{x:.1f}" for x in n1) + f"; ned_def {nd}")
            r["per_group"][g["key"]] = dict(cells=cells, recirc1=n1)
    open(os.path.join(OUT, "sw3_report.log"), "w").write("\n".join(out) + "\n")
    json.dump(res, open(os.path.join(OUT, "sw3.json"), "w"), indent=1)
    print("\n".join(out))


if __name__ == "__main__":
    main()
