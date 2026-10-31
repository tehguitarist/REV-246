"""Report and law fit for tools/ml_fix_depth.py: analysis/ml_fix/depth/report.log, depth.json.

    .venv/bin/python tools/ml_fix_depth_report.py

Per M RT step and theta variant (Python core through the sweep, 3 phases): SW = model vs capture, SPP = phase vs phase,
d0 = theta 0 through the sweep (the static model, = MC), [CC] = Tomkin's capture pairs (analysis/pulse/ml). Nulls: SPP is
corrected by -1.76 dB (3 evenly spaced phases: mean 1 - cos = 1.5 against random pairs' 1). Law fit: per RT step, the
variant minimising a cost on the EVEN groups (sorted keys), scored on the ODD ones.
Cost = mean over the decorrelation windows of |CC - SPP_corrected| (dB) / 1 dB + mean over the low-band EDC metrics of
(SW - CC) / CC, i.e. depth from the decorrelation curve, with the metrics as a check against over-fitting one of them.
"""
import glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ.setdefault("EMT_KERNEL", "chain3")
import numpy as np  # noqa: E402
import pulse_check as pc  # noqa: E402
from defects import ML_SLIPPED  # noqa: E402

D = os.path.join(ROOT, "analysis", "ml_fix", "depth")
REF = os.path.join(ROOT, "analysis", "pulse", "ml")
SW3 = os.path.join(ROOT, "analysis", "ml_fix", "sw3_all")
MET = ["oct_t125", "oct_t250", "oct_t500", "oct_t1000", "oct_t2000", "oct_edt125", "oct_edt250", "oct_edt500", "oct_edt1000",
       "oct_late125", "oct_late250", "oct_late500", "oct_dev125", "oct_dev250", "oct_dev500"]
LOW = ["oct_t125", "oct_t250", "oct_t500", "oct_edt125", "oct_edt250", "oct_edt500"]
NUL = ["recirc1", "recirc2", "recirc3", "0.3-0.6", "0.6-1.0", "1.0-2.0", "2.0-4.0"]
DEC = ["recirc1", "recirc2", "0.3-0.6", "0.6-1.0"]             # the depth estimator's windows (small-theta regime first)
CORR = 10 * np.log10(1.5)
RUN = os.environ.get("DEPTH_DIR", "runs")
TAG = "" if RUN == "runs" else "_" + RUN
sc = lambda m: 1.0 if m.startswith("oct_dev") else 100.0


def med(x, f=abs):
    x = [f(v) for v in x if v is not None and np.isfinite(v)]
    return float(np.median(x)) if x else None


def load():
    G = []
    for fn in sorted(glob.glob(os.path.join(D, RUN, "M_*.json"))):
        d = json.load(open(fn))
        key = d["key"]
        bad = lambda a: ("M", key, a) in ML_SLIPPED
        ref = json.load(open(os.path.join(REF, os.path.basename(fn))))
        vs = list(d["variants"])
        g = dict(key=key, rt=d["rt"], lf=d["lf"], hf=d["hf"], vars=vs, SW={v: [] for v in vs + ["d0"]}, SWn={v: [] for v in vs + ["d0"]},
                 SPP={v: [] for v in vs}, SPPn={v: [] for v in vs}, CC=[], CCn=[], PL=[], PLn=[])
        for r in d["rows"]:
            if bad(r["amp"]):
                continue
            for v in vs:
                for k in range(3):
                    g["SW"][v].append(r["vs_cap"][f"{v}_{k}"]); g["SWn"][v].append(r["null_cap"][f"{v}_{k}"])
            g["SW"]["d0"].append(r["vs_cap"]["d0"]); g["SWn"]["d0"].append(r["null_cap"]["d0"])
            for pk_, x in r["pairs"].items():
                v = pk_.split("_")[0]
                g["SPP"][v].append(x); g["SPPn"][v].append(r["null_pairs"][pk_])
        for p in ref["pairs"]:
            if not any(bad(a) for a in p["pair"]):
                g["CC"].append(pc.summarise(p["m"])); g["CCn"].append(p["null"])
        sw3 = os.path.join(SW3, os.path.basename(fn))      # the plugin (sweep, 3 phases, this build) for the check
        if os.path.exists(sw3):
            for r in json.load(open(sw3))["rows"]:
                if not bad(r["amp"]):
                    for v in ("s0", "s1", "s2"):
                        g["PL"].append(r["vs_cap"][v]); g["PLn"].append(r["null_cap"][v])
        G.append(g)
    return G


def within(gs, m, rows_of):
    ok = n = 0
    floor = 0.0 if m.startswith("oct_dev") else 0.03
    for g in gs:
        rep = [abs(s.get(m)) for s in g["CC"] if s.get(m) is not None]
        if not rep:
            continue
        tol = max(max(rep), floor)
        x = [abs(s.get(m)) for s in rows_of(g) if s.get(m) is not None]
        ok += sum(v <= tol for v in x); n += len(x)
    return ok / n if n else None


def stats(gs, v):
    """Medians for variant v (or 'd0', 'PL') over groups gs."""
    if v == "PL":
        rows, nrows, spp, sppn = (lambda g: g["PL"]), (lambda g: g["PLn"]), None, None
    else:
        rows, nrows = (lambda g: g["SW"][v]), (lambda g: g["SWn"][v])
        spp = (lambda g: g["SPP"][v]) if v != "d0" else None
        sppn = (lambda g: g["SPPn"][v]) if v != "d0" else None
    o = dict(m={}, tol={}, n={}, spp={}, sppn={})
    for m in MET:
        o["m"][m] = med([s.get(m) for g in gs for s in rows(g)])
        o["tol"][m] = within(gs, m, rows)
        if spp:
            o["spp"][m] = med([s.get(m) for g in gs for s in spp(g)])
    for w in NUL:
        o["n"][w] = med([x.get(w) for g in gs for x in nrows(g)], f=lambda z: z)
        if sppn:
            x = med([x.get(w) for g in gs for x in sppn(g)], f=lambda z: z)
            o["sppn"][w] = None if x is None else x - CORR
    return o


def cc_stats(gs):
    return dict(m={m: med([s.get(m) for g in gs for s in g["CC"]]) for m in MET},
                n={w: med([x.get(w) for g in gs for x in g["CCn"]], f=lambda z: z) for w in NUL})


def cost(gs, v):
    s, c = stats(gs, v), cc_stats(gs)
    dec = [abs(c["n"][w] - s["sppn"][w]) for w in DEC if c["n"][w] is not None and s["sppn"].get(w) is not None]
    met = [(s["m"][m] - c["m"][m]) / c["m"][m] for m in LOW if s["m"][m] is not None and c["m"][m]]
    return float(np.mean(dec) + np.mean(met)) if dec and met else None, s, c


def fc(x):
    return "-" if x is None else f"{x:.2f}"


def main():
    G = load()
    vs = G[0]["vars"]
    rts = sorted({g["rt"] for g in G})
    out, res = [], dict(rt={}, law={}, heldout={})
    f2 = lambda x, k=100.0: "-" if x is None else f"{x * k:.2f}"
    out.append(f"# M theta depth scan: {len(G)} groups, variants {vs} (rms rad), Python core through the sweep, 3 phases")
    # check: the Python d13 against the plugin's own sw3 (the same theta, sweep and phases)
    pl, py = stats(G, "PL"), stats(G, "d13") if "d13" in vs else None
    if py:
        out.append("check, plugin sw3 / python d13 (all M): " + "; ".join(f"{m} {f2(pl['m'][m], sc(m))} / {f2(py['m'][m], sc(m))}" for m in LOW)
                   + "; recirc1 " + f"{pl['n']['recirc1']:.1f} / {py['n']['recirc1']:.1f}")
    for rt in rts:
        gs = [g for g in G if g["rt"] == rt]
        c = cc_stats(gs)
        out.append(f"\n## RT {rt} ({len(gs)} groups): variant | " + " | ".join(LOW) + " | dev125 | tol 500 / edt250 / edt500 | SPP 500 / edt500 | null gap CC - SPPcorr " + "/".join(DEC) + " | cost")
        out.append(f"[CC] | " + " | ".join(f2(c["m"][m], sc(m)) for m in LOW) + f" | {f2(c['m']['oct_dev125'], 1)} | | | " + " / ".join(f"{c['n'][w]:.1f}" for w in DEC if c["n"][w] is not None))
        res["rt"][str(rt)] = {}
        for v in ["d0"] + vs:
            if v == "d0":
                s = stats(gs, v); cst = None
            else:
                cst, s, _ = cost(gs, v)
            gap = [None if s["sppn"].get(w) is None or c["n"][w] is None else c["n"][w] - s["sppn"][w] for w in DEC] if v != "d0" else []
            res["rt"][str(rt)][v] = dict(m=s["m"], tol=s["tol"], spp=s["spp"], nulls=s["n"], sppn=s["sppn"], gap=gap, cost=cst)
            out.append(f"{v} | " + " | ".join(f2(s["m"][m], sc(m)) for m in LOW) + f" | {f2(s['m']['oct_dev125'], 1)} | "
                       + " / ".join("-" if s["tol"][m] is None else f"{100 * s['tol'][m]:.0f}" for m in ("oct_t500", "oct_edt250", "oct_edt500"))
                       + (f" | {f2(s['spp'].get('oct_t500'))} / {f2(s['spp'].get('oct_edt500'))} | " + " / ".join("-" if x is None else f"{x:+.1f}" for x in gap) + f" | {'-' if cst is None else f'{cst:.2f}'}" if v != "d0" else ""))
    # law: per RT step the best variant on the even groups, scored on the odd
    ev = lambda g: sorted(x["key"] for x in G if x["rt"] == g["rt"]).index(g["key"]) % 2 == 0
    law = {}
    out.append("\n## Law fit: per RT step, best variant on even groups (cost), held out on odd")
    for rt in rts:
        gE = [g for g in G if g["rt"] == rt and ev(g)]
        gO = [g for g in G if g["rt"] == rt and not ev(g)]
        cE = {v: cost(gE, v)[0] for v in vs}
        cO = {v: cost(gO, v)[0] for v in vs}
        best = min(cE, key=lambda v: cE[v] if cE[v] is not None else 9e9)
        bestO = min(cO, key=lambda v: cO[v] if cO[v] is not None else 9e9)
        law[rt] = best
        out.append(f"- RT {rt}: even " + " ".join(f"{v} {fc(cE[v])}" for v in vs) + f" -> {best}; odd " + " ".join(f"{v} {fc(cO[v])}" for v in vs) + f" (odd's own best {bestO})")
    res["law"] = {str(k): v for k, v in law.items()}
    # held-out before / after per RT: odd groups, d13 (current) against the law's variant
    out.append("\n## Held out (odd groups): current d13 -> law, per RT; metric medians |%| [CC], in-tol %, nulls vs [CC]")
    for rt in rts:
        gO = [g for g in G if g["rt"] == rt and not ev(g)]
        a, b, c = stats(gO, "d13"), stats(gO, law[rt]), cc_stats(gO)
        pct = lambda x: "-" if x is None else f"{100 * x:.0f}"
        cells = [f"{m} {f2(a['m'][m], sc(m))} -> {f2(b['m'][m], sc(m))} [{f2(c['m'][m], sc(m))}] ({pct(a['tol'][m])}->{pct(b['tol'][m])}%)" for m in MET]
        nl = " ".join(f"{w} {a['n'][w]:.1f}->{b['n'][w]:.1f} [{c['n'][w]:.1f}]" for w in NUL if a["n"][w] is not None and c["n"][w] is not None)
        out.append(f"- RT {rt} ({len(gO)} groups, {law[rt]}): " + "; ".join(cells) + f"; nulls {nl}")
        res["heldout"][str(rt)] = dict(variant=law[rt], before=a, after=b, cc=c)
    open(os.path.join(D, f"report{TAG}.log"), "w").write("\n".join(out) + "\n")
    json.dump(res, open(os.path.join(D, f"depth{TAG}.json"), "w"), indent=1, default=lambda x: None)
    print("\n".join(out))


if __name__ == "__main__":
    main()
