"""Report for tools/ml_fix_diag.py: analysis/ml_fix/diag.json and diag.log.

    .venv/bin/python tools/ml_fix_report.py

Columns (medians over captures x phases / pairs, R channel):
  PC   plugin (impulse, modulation on, 3 phases) vs capture        PC0  phase 0 only (= the audit's render, timing apart)
  SC   plugin with modulation off vs capture                       MC   pulse_check's static model vs capture
  CC   capture vs capture (every pair of the group: the tolerance)
  PP   plugin vs plugin, phases 120 deg apart (same capture frame and window)
  PS   modulated vs static plugin (signed mean = what the modulation does to the decay on average: the bias test)
  SW*  the plugin measured through ml_variance's sweep (6 phases, 60 deg apart): SWS vs the static plugin, SWP sweep vs
       sweep, SWC sweep vs capture; IC the impulse phases vs capture in the same (sweep-capped) window
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
MET = ["oct_t125", "oct_t250", "oct_t500", "oct_t1000", "oct_t2000", "oct_edt125", "oct_edt250", "oct_edt500",
       "oct_dev125", "oct_dev250", "oct_dev500"] + [f"third_t{f}" for f in pc.THIRD if f <= 400] + [
       "rv_mix_ms", "ned"]
NUL = ["recirc1", "recirc2", "recirc3", "0.3-0.6", "0.6-1.0", "1.0-2.0"]
SC_ = lambda k: 1.0 if k.startswith(("oct_dev", "ned")) else 100.0


def med(x, f=abs):
    x = [f(v) for v in x if v is not None and np.isfinite(v)]
    return float(np.median(x)) if x else None


def mean(x):
    x = [v for v in x if v is not None and np.isfinite(v)]
    return float(np.mean(x)) if x else None


def sem(x):
    x = [v for v in x if v is not None and np.isfinite(v)]
    return float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 2 else None


def load():
    G = {}
    for fn in sorted(glob.glob(os.path.join(OUT, "ml", "*.json"))):
        d = json.load(open(fn))
        prog, key = d["prog"], d["key"]
        bad = lambda a: (prog, key, a) in ML_SLIPPED
        ref = json.load(open(os.path.join(REF, os.path.basename(fn))))
        aud = json.load(open(os.path.join(AUD, os.path.basename(fn))))
        g = dict(prog=prog, key=key, rt=d["rt"], lf=d["lf"], hf=d["hf"], v={k: [] for k in ("PC", "PC0", "SC", "MC", "CC", "PP", "PS", "AUD")},
                 n={k: [] for k in ("PC", "PC0", "SC", "MC", "CC", "PP", "PS")})
        for r in d["rows"]:
            if bad(r["amp"]):
                continue
            for v in ("p0", "p1", "p2"):
                g["v"]["PC"].append(r["vs_cap"][v]); g["n"]["PC"].append(r["null_cap"][v])
                g["v"]["PS"].append(r["vs_static"][v]); g["n"]["PS"].append(r["null_static"][v])
            g["v"]["PC0"].append(r["vs_cap"]["p0"]); g["n"]["PC0"].append(r["null_cap"]["p0"])
            g["v"]["SC"].append(r["vs_cap"]["static"]); g["n"]["SC"].append(r["null_cap"]["static"])
            for k, x in r["pairs"].items():
                g["v"]["PP"].append(x); g["n"]["PP"].append(r["null_pairs"][k])
        for r in ref["rows"]:
            if not bad(r["amp"]):
                g["v"]["MC"].append(pc.summarise(r["m"])); g["n"]["MC"].append(r["null"])
        for p in ref["pairs"]:
            if not any(bad(a) for a in p["pair"]):
                g["v"]["CC"].append(pc.summarise(p["m"])); g["n"]["CC"].append(p["null"])
        for r in aud["rows"]:
            if not bad(r["amp"]):
                g["v"]["AUD"].append(pc.summarise(r["m"]))
        G[f"{prog}_{key}"] = g
    return G


def sweep_rows():
    p = os.path.join(OUT, "sweep.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def main():
    G, SW = load(), sweep_rows()
    port = json.load(open(os.path.join(OUT, "port.json"))) if os.path.exists(os.path.join(OUT, "port.json")) else {}
    out, res = [], dict(groups={}, pooled={}, port={}, sweep={})
    cols = ["PC0", "AUD", "PC", "SC", "MC", "CC", "PP", "PS"]
    for prog in ("M", "L"):
        gs = [g for g in G.values() if g["prog"] == prog]
        if not gs:
            continue
        out.append(f"\n## {prog}: {len(gs)} groups ({', '.join(g['key'] for g in gs)})")
        out.append("metric | " + " | ".join(cols) + " | PC signed mean (sem) | SC signed | PS signed mean (sem)")
        pool = {}
        for m in MET:
            sc = SC_(m)
            vals = {c: [s.get(m) for g in gs for s in g["v"][c]] for c in cols}
            row = {c: med(vals[c]) for c in cols}
            sg = dict(PC=(mean(vals["PC"]), sem(vals["PC"])), SC=(mean(vals["SC"]), sem(vals["SC"])), PS=(mean(vals["PS"]), sem(vals["PS"])))
            pool[m] = dict(abs_median={c: (None if row[c] is None else row[c] * sc) for c in cols},
                           signed={c: [None if x is None else x * sc for x in v] for c, v in sg.items()})
            f = lambda x: "-" if x is None else f"{x * sc:.2f}"
            fs = lambda t: "-" if t[0] is None else f"{t[0] * sc:+.2f} ({(t[1] or 0) * sc:.2f})"
            out.append(f"{m} | " + " | ".join(f(row[c]) for c in cols) + f" | {fs(sg['PC'])} | {fs(sg['SC'])} | {fs(sg['PS'])}")
        out.append("null (dB, median) | " + " | ".join(c for c in cols if c != "AUD"))
        for w in NUL:
            v = {c: med([n.get(w) for g in gs for n in g["n"][c]], f=lambda x: x) for c in cols if c != "AUD"}
            pool[f"null_{w}"] = v
            out.append(f"{w} | " + " | ".join("-" if v[c] is None else f"{v[c]:.1f}" for c in v))
        res["pooled"][prog] = pool
        # per group: the low band and recirc1
        fb = 125 if prog == "M" else 250
        out.append(f"\n### {prog} per group: {fb} Hz decay |%| PC0 / PC / SC / MC / CC / PP, PS signed mean; EDT 500 PC / MC / CC / PP; recirc1 null PC / CC / PP; mix % PC / MC / CC / PP")
        for g in gs:
            q = lambda c, m, sc=100: med([s.get(m) for s in g["v"][c]])
            r = {c: q(c, f"oct_t{fb}") for c in ("PC0", "PC", "SC", "MC", "CC", "PP")}
            ps = mean([s.get(f"oct_t{fb}") for s in g["v"]["PS"]])
            e5 = {c: q(c, "oct_edt500") for c in ("PC", "MC", "CC", "PP")}
            mx = {c: q(c, "rv_mix_ms") for c in ("PC", "MC", "CC", "PP")}
            n1 = {c: med([n.get("recirc1") for n in g["n"][c]], f=lambda x: x) for c in ("PC", "CC", "PP")}
            f = lambda x: "-" if x is None else f"{100 * x:.1f}"
            out.append(f"- {g['key']}: " + " / ".join(f(r[c]) for c in r) + f", PS {'-' if ps is None else f'{100 * ps:+.1f}'}; EDT500 "
                       + " / ".join(f(e5[c]) for c in e5) + "; recirc1 " + " / ".join("-" if n1[c] is None else f"{n1[c]:.1f}" for c in n1)
                       + "; mix " + " / ".join(f(mx[c]) for c in mx))
            res["groups"][f"{prog}_{g['key']}"] = dict(decay_fb=r, PS_signed=ps, edt500=e5, recirc1=n1, mix=mx, band=fb)
    # sweep vs impulse
    if SW:
        out.append("\n## Measurement: the plugin through the sweep (6 phases) against its impulse (3 phases), the 8 SIM_KEYS groups")
        for prog in ("M", "L"):
            ks = [k for k in SW if k.startswith(prog + "_")]
            if not ks:
                continue
            out.append(f"\n### {prog} ({len(ks)} groups): metric | SWS |.| | SWS signed mean (sem) | imp vs static |.| | imp signed mean (sem) | SWP |.| | IPP |.| | SWC |.| | IC |.| | SWC signed | IC signed | CC |.|")
            sres = {}
            for m in [x for x in MET if x.startswith(("oct_t", "oct_edt", "third_t", "oct_dev"))]:
                sc = SC_(m)
                sws = [SW[k]["vs_static"][f"sw{i}"].get(m) for k in ks for i in range(6)]
                imp = [SW[k]["vs_static"][v].get(m) for k in ks for v in ("p0", "p1", "p2")]
                def pairs(names):
                    o = []
                    for k in ks:
                        x = [SW[k]["vs_static"][v].get(m) for v in names]
                        for i in range(len(x)):
                            for j in range(i + 1, len(x)):
                                if x[i] is not None and x[j] is not None:
                                    o.append((x[i] - x[j]) if m.startswith("oct_dev") else (1 + x[i]) / (1 + x[j]) - 1)
                    return o
                swp, ipp = pairs([f"sw{i}" for i in range(6)]), pairs(["p0", "p1", "p2"])
                swc = [SW[k]["vs_cap"][a][f"sw{i}"].get(m) for k in ks for a in SW[k]["vs_cap"] if (prog, k[2:], int(a)) not in ML_SLIPPED for i in range(6)]
                ic = [SW[k]["vs_cap"][a][v].get(m) for k in ks for a in SW[k]["vs_cap"] if (prog, k[2:], int(a)) not in ML_SLIPPED for v in ("p0", "p1", "p2")]
                cc = [s.get(m) for k in ks if k in G for s in G[k]["v"]["CC"]]
                f = lambda x: "-" if x is None else f"{x * sc:.2f}"
                fs = lambda x: "-" if x is None else f"{x * sc:+.2f}"
                out.append(f"{m} | {f(med(sws))} | {fs(mean(sws))} ({f(sem(sws))}) | {f(med(imp))} | {fs(mean(imp))} ({f(sem(imp))}) | {f(med(swp))} | {f(med(ipp))} | "
                           f"{f(med(swc))} | {f(med(ic))} | {fs(mean(swc))} | {fs(mean(ic))} | {f(med(cc))}")
                sres[m] = dict(SWS_abs=med(sws), SWS_mean=mean(sws), SWS_sem=sem(sws), IMP_abs=med(imp), IMP_mean=mean(imp), IMP_sem=sem(imp),
                               SWP_abs=med(swp), IPP_abs=med(ipp), SWC_abs=med(swc), IC_abs=med(ic), SWC_mean=mean(swc), IC_mean=mean(ic), CC_abs=med(cc), scale=sc)
            nl = {w: (med([SW[k]["null_static"][f"sw{i}"].get(w) for k in ks for i in range(6)], f=lambda x: x),
                      med([SW[k]["null_static"][v].get(w) for k in ks for v in ("p0", "p1", "p2")], f=lambda x: x)) for w in ("first_pass", "recirc1", "recirc2", "recirc3")}
            out.append("null vs static, sweep / impulse: " + "; ".join(f"{w} {a:.1f} / {b:.1f}" for w, (a, b) in nl.items() if a is not None and b is not None))
            sres["null_vs_static"] = nl
            res["sweep"][prog] = sres
    if port:
        out.append("\n## Port check: C++ core (device domain) vs the Python reference (ml_variance.core), dB, worst window")
        for k, r in port.items():
            st = max(r["static"].values())
            md = max(max(r[v]["null"].values()) for v in ("p0", "p1", "p2"))
            vs = {v: r[v]["vs_static"] for v in ("p0", "p1", "p2")}
            out.append(f"- {k}: static {st:.1f}; modulated (theta origin d = {[r[v]['d'] for v in ('p0', 'p1', 'p2')]}) {md:.1f}")
            res["port"][k] = dict(static_worst=st, modulated_worst=md, d=[r[v]["d"] for v in ("p0", "p1", "p2")])
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(out) + "\n")
    json.dump(res, open(os.path.join(OUT, "diag.json"), "w"), indent=1)
    print("\n".join(out))


if __name__ == "__main__":
    main()
