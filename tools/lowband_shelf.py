"""S's LF shelf below 400 Hz (model_s_results §25): what limits the §22.9 baseline at 25-200 Hz at LF 0.4 / 0.1.

    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/lowband_shelf.py --stage free [--only keys]
    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/lowband_shelf.py --stage dc [--keep-lfrt]
    EMT_KERNEL=chain3 .venv/bin/python tools/lowband_shelf.py --stage report [--keep-lfrt]

free    per setting, the baseline's shelf replaced by one with pole a + da and the DC / Nyquist per-pass gains changed by
        ddc / dhi dB x P_i / mean P (length-scaled, as the laws), 3 numbers fitted on the 25-1000 Hz residual from 150 ms
        to each capture's compensated end, held out by capture (fit on the other captures, score the one left out).
        At LF 1.0 there is no shelf and da is undetermined (a control only).
dc      LF 0.4: the law's shelf DC ratio (0.412) replaced by one constant on a grid, the lfrt LF 0.4 term removed
        (--keep-lfrt: kept); per setting the 25-1000 Hz cost and the band nulls of every capture at each grid value.
report  free: the fitted numbers and nulls by LF / RT; dc: the constant chosen leaving one RT step out (median of the
        per-setting best over the other RT steps), each RT step scored with its own constant, against the baseline and
        the repeat capture (lowband.py's score stage on chain3).

One result file per setting (analysis/lowband/<kernel>_shelf_<stage>/), skipped on rerun. Model, captures and fits as
lowband.py (q2.final_G with the fold's loop terms, q2.fit_capture); the loop responses are re-rendered in the frequency
domain (reverb = F PB + sum_i X_i PB / (1 - G_i z^-P_i), the cluster added in time), exact against fit_capture's model.
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lowband as lb  # noqa: E402

q2, q1 = lb.q2, lb.q1
m = q1.m
FS, PRE = q2.FS, q2.PRE
SOSF = butter(3, (25, 1000), btype="band", fs=FS, output="sos")    # the fitting band
PBAR = float(np.mean(q1.P))
GRID = np.round(np.arange(0.405, 0.43001, 0.00025), 5)               # dc stage: LF 0.4's DC ratio
LF_DC = 0.4
WS, BANDS = ("150-300", "300-600", "600-1000"), ("25-50", "50-100", "100-200", "100-10k")


def outdir(stage, keep=False):
    return os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_{stage}" + ("K" if keep else ""))


def base_r0(lab, law):
    """Per loop: the shelf's effective DC ratio exactly as q2.G_capped builds it (None at LF 1.0)."""
    rt, L = lab["rt"], law["lf"][lab["lf"]]
    out = []
    for Pi in q1.P:
        kk = -3 * Pi / q1.FS
        g = 10 ** (kk / (L["hi"] * rt))
        if L["dc"] == L["hi"]:
            out.append(None)
            continue
        r0 = 10 ** (kk * (1 / (L["dc"] * rt) - 1 / (L["hi"] * rt)))
        if g * r0 < q2.DC_MIN:
            g = g * 10 ** (-q2.FINAL_BETA * 20 * np.log10(q2.DC_MIN / (g * r0)) / 20)
        out.append(max(r0, q2.DC_MIN / g))
    return out


def prep(key, terms=None):
    """The baseline for one setting: loop responses, excitation spectra, cluster, and every capture fitted (clock,
    shift, gain, fade) as lowband.py. terms: a replacement for lowband.fold_terms."""
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, lb.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=lb.lf_.COMPACT)
    q1.EXT = 1700
    lab = dict(lf=info[key]["lf"], rt=info[key]["rt"], hf=info[key]["hf"])
    b230 = ap["steps"][ap["settings"][key]["step"]]
    keep = []
    G = q2.final_G(lab, b230, law, **(terms or lb.fold_terms)(lab["rt"]), Gmod=lambda l_, w_, ww, Gs: keep.append(Gs) or Gs)
    caps = []
    for amp, nm in q2.files_of(info[key]):
        y, mdl, g, te = q2.fit_capture(G, v, amp, nm, fade[nm], eps0[nm], ends[nm], 1700)
        n = len(y)
        i0, i1 = PRE + int(0.15 * FS), min(PRE + int(te * FS), n)
        caps.append(dict(amp=amp, y=y, mdl=mdl, g=g, te=te, norm=float(np.sqrt(np.sum(sosfiltfilt(SOSF, y, axis=0)[i0:i1] ** 2)))))
    N = G["N"]
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((N, 2)), w=np.ones(N), amp=0)])])
    cl, f, xs = J.split(v)
    F, X = J._spectra(f, xs, N)
    return dict(key=key, lab=lab, law=law, G=G, N=N, w=2 * np.pi * np.fft.rfftfreq(N), F=F, X=X, Gs0=keep[0],
                clus=J.cluster(cl), caps=caps, r0=base_r0(lab, law))


def loops(d, th):
    """Loop responses with the shelf replaced: pole a + da, DC and Nyquist per-pass gain x 10^(dB x P_i / Pbar / 20)."""
    da, ddc, dhi = th
    a = d["law"]["a"]
    out = []
    for Gi, Pi, r0 in zip(d["Gs0"], q1.P, d["r0"]):
        s = Pi / PBAR
        r0 = 1.0 if r0 is None else r0
        out.append(Gi * m.shelf(d["w"], r0 * 10 ** (ddc * s / 20), 10 ** (dhi * s / 20), a + da) / m.shelf(d["w"], r0, 1.0, a))
    return out


def reverb(d, Gs):
    w = d["w"]
    S = d["F"] * d["G"]["PB"][None]
    for i, (Gi, Pi) in enumerate(zip(Gs, q1.P)):
        S = S + d["X"][:, i] * (d["G"]["PB"] / (1 - Gi * np.exp(-1j * w * Pi)))[None]
    return np.fft.irfft(S, d["N"], axis=1).T


def model(d, c, r):
    n = len(c["y"])
    return (c["amp"] / 90) * d["clus"][:n] + r[:n]


def residual(d, Gs, caps, t0=0.15):
    r = reverb(d, Gs)
    out = []
    for c in caps:
        n = len(c["y"])
        i0, i1 = PRE + int(t0 * FS), min(PRE + int(c["te"] * FS), n)
        e = sosfiltfilt(SOSF, c["y"] - model(d, c, r) * c["g"][:n], axis=0)[i0:i1]      # g is (n, 2)
        out.append(e.ravel() / c["norm"])
    return np.concatenate(out)


def free(key, path):
    d = prep(key)
    res = {}
    for ho, c in enumerate(d["caps"]):
        tr = [x for j, x in enumerate(d["caps"]) if j != ho]
        fit = least_squares(lambda th: residual(d, loops(d, th), tr), np.zeros(3), x_scale=[1e-3, 0.05, 0.05], diff_step=1e-3, max_nfev=60)
        r = reverb(d, loops(d, fit.x))
        res[c["amp"]] = dict(th=fit.x.tolist(), base=lb.band_nulls(c["y"], c["mdl"] * c["g"], c["te"]),
                             fit=lb.band_nulls(c["y"], model(d, c, r) * c["g"], c["te"]))
    json.dump(dict(key=key, **d["lab"], res=res), open(path, "w"))


def dc(key, path, keep):
    def terms(rt):
        t = lb.fold_terms(rt)
        t["lfrt"] = json.loads(json.dumps(t["lfrt"]))
        if not keep:
            t["lfrt"]["coefs"].pop(str(LF_DC), None)
        return t
    base = prep(key)
    d = prep(key, terms)
    law, lab = d["law"], d["lab"]
    L = law["lf"][LF_DC]
    res = dict(key=key, **lab, law_dc=L["dc"], amps=[c["amp"] for c in d["caps"]], grid=[],
               base=[lb.band_nulls(c["y"], c["mdl"] * c["g"], c["te"]) for c in base["caps"]])
    for x in GRID:
        r0n = base_r0(lab, dict(law, lf={**law["lf"], LF_DC: dict(L, dc=float(x))}))
        Gs = [Gi * m.shelf(d["w"], rn, 1.0, law["a"]) / m.shelf(d["w"], ro, 1.0, law["a"]) for Gi, rn, ro in zip(d["Gs0"], r0n, d["r0"])]
        r = reverb(d, Gs)
        res["grid"].append(dict(dc=float(x), cost=float(np.sum(residual(d, Gs, d["caps"]) ** 2)),
                                nulls=[lb.band_nulls(c["y"], model(d, c, r) * c["g"], c["te"]) for c in d["caps"]]))
    json.dump(res, open(path, "w"))


def run(stage, only, keep):
    info = q2.q1_tables()[0]
    od = outdir(stage, keep)
    os.makedirs(od, exist_ok=True)
    keys = sorted(k for k in info if (only is None or k in only) and (stage == "free" or info[k]["lf"] == LF_DC))
    keys = [k for k in keys if not os.path.exists(os.path.join(od, k + ".json"))]
    print(f"{stage}: {len(keys)} settings to do", flush=True)
    for i, k in enumerate(keys):                                         # one process, ~0.6 GB (measured)
        free(k, os.path.join(od, k + ".json")) if stage == "free" else dc(k, os.path.join(od, k + ".json"), keep)
        print(f"  {i + 1}/{len(keys)} {k}", flush=True)


def report_free():
    D = [json.load(open(p)) for p in glob.glob(os.path.join(outdir("free"), "*.json"))]
    med = lambda u, f, b, w: np.median([r[f][b][w][1][0] for d in u for r in d["res"].values() if w in r[f][b]])
    out = ["## Free shelf per setting, held out by capture (3 numbers per setting); median over HF and captures; "
           "R null base -> fit, 150-300 | 300-600 ms"]
    for lf in sorted({d["lf"] for d in D}):
        for rt in sorted({d["rt"] for d in D}):
            u = [d for d in D if d["lf"] == lf and d["rt"] == rt]
            th = np.median([r["th"] for d in u for r in d["res"].values()], 0)
            out.append(f"- LF {lf} RT {rt}: da {th[0]:+.5f} ddc {th[1]:+.4f} dhi {th[2]:+.4f} dB/pass; " + "; ".join(
                f"{b} " + " | ".join(f"{med(u, 'base', b, w):.1f} -> {med(u, 'fit', b, w):.1f}" for w in WS[:2]) for b in BANDS[1:]))
    return out


def report_dc(keep):
    D = {d["key"]: d for d in (json.load(open(p)) for p in glob.glob(os.path.join(outdir("dc", keep), "*.json")))}
    rep = {(k, r["amp"]): r["repeat"] for k in D for r in json.load(open(os.path.join(lb.outdir("base"), k + ".json")))}
    best = {k: d["grid"][int(np.argmin([g["cost"] for g in d["grid"]]))]["dc"] for k, d in D.items()}
    rts = sorted({d["rt"] for d in D.values()})
    out = [f"## LF 0.4 shelf DC ratio, one constant, lfrt LF 0.4 {'kept' if keep else 'removed'}",
           "in-sample best dc by RT (median over HF): " + "; ".join(f"{rt}: {np.median([best[k] for k in D if D[k]['rt'] == rt]):.5f}" for rt in rts)]
    for gate in ((0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5), (0.8, 1.0, 1.4, 1.8, 2.2, 2.5)):
        out.append(f"\n### dc chosen leaving one RT step out (median of the per-setting best over training RT steps in {gate}); R, median over captures")
        rows = {"base": [], "held": []}
        for rt in rts:
            x = float(np.median([best[k] for k in D if D[k]["rt"] != rt and D[k]["rt"] in gate]))
            for k, d in D.items():
                if d["rt"] == rt:
                    gi = int(np.argmin(np.abs(np.array([g["dc"] for g in d["grid"]]) - x)))
                    for j, amp in enumerate(d["amps"]):
                        rows["base"].append((rt, amp, k, d["base"][j]))
                        rows["held"].append((rt, amp, k, d["grid"][gi]["nulls"][j]))
            out.append(f"  RT {rt}: dc {x:.5f}")

        def cell(rs, rt, b, w, excess):
            x = [n[b][w][1][0] - (rep[(k, a)][b][w][1][0] if excess else 0) for r_, a, k, n in rs
                 if (rt is None or r_ == rt) and w in n[b] and (not excess or w in rep[(k, a)][b])]
            return np.median(x) if x else np.nan
        for b in BANDS:
            out.append(f"#### {b} Hz  null base -> held out [excess over repeat base -> held out], " + " | ".join(WS))
            for rt in [None] + rts:
                out.append(f"  {'all' if rt is None else 'RT ' + str(rt):7s} " + " | ".join(
                    f"{cell(rows['base'], rt, b, w, False):.1f} -> {cell(rows['held'], rt, b, w, False):.1f} "
                    f"[{cell(rows['base'], rt, b, w, True):+.1f} -> {cell(rows['held'], rt, b, w, True):+.1f}]" for w in WS))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["free", "dc", "report"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--keep-lfrt", action="store_true")
    a = ap.parse_args()
    if a.stage == "report":
        out = report_free() + [""] + report_dc(False) + [""] + report_dc(True)
        open(os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_report.log"), "w").write("\n".join(out) + "\n")
        print("\n".join(out))
    else:
        run(a.stage, a.only.split(",") if a.only else None, a.keep_lfrt)


if __name__ == "__main__":
    main()
