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
      [f"third_t{f}" for f in pc.THIRD] + ["rv_mix_ms", "ned_def", "ned"]
LEV = [f"rv_{x}_{f}" for x in ("c50", "c80") for f in (125, 500, 2000)] + [f"rv_lev_e_{f}" for f in (100, 200, 400, 800, 1600, 3200, 6400)] + \
      [f"rv_lev_l_{f}" for f in (100, 200, 400, 800, 1600, 3200, 6400)] + ["rv_lr_e", "rv_lr_l", "rv_iacc_e", "rv_iacc_l"] + \
      [f"coh_{f}" for f in (200, 400, 800, 1600, 3200)] + [f"rv_fsd_{w}" for w in ("100-1000", "1000-4000", "4000-10000")]
DIRN = os.environ.get("SW3_DIR", "sw3")
SFX = "" if DIRN == "sw3" else DIRN[3:]
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
    return 1.0 if (k.startswith(("oct_dev", "ned", "rv_", "coh_")) and k not in ("ned_def", "rv_mix_ms")) else 100.0


def med(x, f=abs):
    x = [f(v) for v in x if v is not None and np.isfinite(v)]
    return float(np.median(x)) if x else None


def load():
    G = []
    for fn in sorted(glob.glob(os.path.join(OUT, DIRN, "*.json"))):
        d = json.load(open(fn))
        prog, key = d["prog"], d["key"]
        bad = lambda a: (prog, key, a) in ML_SLIPPED
        ref = json.load(open(os.path.join(REF, os.path.basename(fn))))
        aud = json.load(open(os.path.join(AUD, os.path.basename(fn))))
        g = dict(prog=prog, key=key, rt=d["rt"], lf=d["lf"], hf=d["hf"], fit=[r.get("fit") for r in d["rows"]], v={c: [] for c in COLS}, n={c: [] for c in COLS}, per={})
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
    floor = 0.03 if sc(m) == 100.0 else 0.0
    for g in gs:
        rep = [abs(s.get(m)) for s in g["v"][ref] if s.get(m) is not None]
        if not rep:
            continue
        tol = max(max(rep), floor)
        x = [abs(s.get(m)) for s in g["v"][col] if s.get(m) is not None]
        ok += sum(v <= tol for v in x)
        n += len(x)
    return ok, n


SUBSETS = {"RT short": lambda g: g["rt"] <= {"M": 1.0, "L": 3.0}[g["prog"]], "RT long": lambda g: g["rt"] >= {"M": 3.0, "L": 8.0}[g["prog"]],
           "LF 0.1": lambda g: g["lf"] == 0.1, "LF 2.0": lambda g: g["lf"] == 2.0, "HF <= 0.2": lambda g: g["hf"] <= 0.2}
KEY_MET = {"M": ["oct_t125", "oct_t250", "oct_t500", "oct_edt125", "oct_edt250", "oct_edt500", "oct_dev125"],
           "L": ["oct_t125", "oct_t250", "oct_t500", "oct_edt125", "oct_edt500", "oct_dev250", "rv_mix_ms"]}


def extra(G, res):
    out = []
    for prog in ("M", "L"):
        gs = [g for g in G if g["prog"] == prog]
        if not gs:
            continue
        r = res[prog]
        out.append(f"\n## {prog}: before (AUD, phase-0 impulse) / after (SW, sweep x 3 phases) [CC], by subset; in CC tol after / before")
        for name, fn in [("All", lambda g: True)] + list(SUBSETS.items()):
            sub = [g for g in gs if fn(g)]
            if not sub:
                continue
            cells = []
            for m in KEY_MET[prog]:
                k = sc(m)
                v = {c: med([s.get(m) for g in sub for s in g["v"][c]]) for c in ("AUD", "SW", "CC")}
                w = {c: within(sub, m, c) for c in ("SW", "AUD")}
                f = lambda x: "-" if x is None else f"{x * k:.2f}"
                pct = lambda t: "-" if not t[1] else f"{100 * t[0] / t[1]:.0f}"
                cells.append(f"{m} {f(v['AUD'])} -> {f(v['SW'])} [{f(v['CC'])}] ({pct(w['SW'])}/{pct(w['AUD'])} %)")
                r.setdefault("subsets", {}).setdefault(name, {})[m] = dict(AUD=v["AUD"], SW=v["SW"], CC=v["CC"], within_SW=w["SW"], within_AUD=w["AUD"])
            out.append(f"- {name} ({len(sub)} groups): " + "; ".join(cells))
        out.append(f"\n### {prog} level / spectrum / L-R (dB or abs, |plugin - capture| median): SW | AUD | SC | MC [CC]  (early metrics' CC is the AMPL change, not a tolerance)")
        for m in LEV:
            v = {c: med([s.get(m) for g in gs for s in g["v"][c]]) for c in ("SW", "AUD", "SC", "MC", "CC")}
            r.setdefault("levels", {})[m] = v
            out.append(f"{m} | " + " | ".join("-" if v[c] is None else f"{v[c]:.3f}" for c in ("SW", "AUD", "SC", "MC")) + f" [{'-' if v['CC'] is None else format(v['CC'], '.3f')}]")
        ratio = [20 * np.log10(f["static"]["gains"][1][0] / f["static"]["gains"][1][1]) for g in gs for f in g["fit"] if f and len(f["static"]["gains"][1]) == 2]
        swr = [20 * np.log10(f["s0"]["gains"][1][-1] / f["static"]["gains"][1][-1]) for g in gs for f in g["fit"] if f]
        if ratio:
            out.append(f"- fitted cluster / reverb gain, static impulse, R: median {np.median(ratio):+.3f} dB [10-90 % {np.percentile(ratio, 10):+.3f} {np.percentile(ratio, 90):+.3f}] (audit build: +11.58)")
            out.append(f"- sweep / impulse reverb gain (the deconvolution's scale, fitted out): {np.median(swr):+.2f} dB [{np.percentile(swr, 10):+.2f} {np.percentile(swr, 90):+.2f}]")
            r["gain_ratio_db"] = [float(np.median(ratio)), float(np.percentile(ratio, 10)), float(np.percentile(ratio, 90))]
    sw = os.path.join(OUT, "static_windows.json")
    if os.path.exists(sw):
        S = json.load(open(sw))
        out.append("\n## Static windows on the static impulse renders (this build), dB, median L / R; audit build (plugin_check ml_fix) in brackets")
        for prog in ("M", "L"):
            ks = [k for k in S if k.startswith(prog + "_")]
            cells = []
            for w in ("cluster_gap", "first_pass"):
                x = np.array([r["null"][w] for k in ks for r in S[k] if w in r["null"] and None not in r["null"][w] and (prog, k[2:], r["amp"]) not in ML_SLIPPED])
                old = []
                for k in ks:
                    fo = os.path.join(ROOT, "analysis", "plugin_check", "ml_fix", f"{k}.json")
                    if os.path.exists(fo):
                        old += [r["plugin"][w] for r in json.load(open(fo))["rows"] if w in r["plugin"] and None not in r["plugin"][w] and (prog, k[2:], r["amp"]) not in ML_SLIPPED]
                old = np.array(old)
                if len(x):
                    cells.append(f"{w} {np.median(x[:, 0]):.1f} / {np.median(x[:, 1]):.1f} (R 90th pct {np.percentile(x[:, 1], 90):.1f}; n {len(x)})"
                                 + (f" [{np.median(old[:, 0]):.1f} / {np.median(old[:, 1]):.1f}]" if len(old) else ""))
                    res[prog].setdefault("static_windows", {})[w] = [float(np.median(x[:, 0])), float(np.median(x[:, 1])), float(np.percentile(x[:, 1], 90))]
            out.append(f"- {prog} ({len(ks)} groups): " + "; ".join(cells))
    return out


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
    out += extra(G, res)
    open(os.path.join(OUT, f"sw3{SFX}_report.log"), "w").write("\n".join(out) + "\n")
    json.dump(res, open(os.path.join(OUT, f"sw3{SFX}.json"), "w"), indent=1)
    print("\n".join(out))


if __name__ == "__main__":
    main()
