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
floor   LF 0.1 (and LF 0.4 / RT 0.3, the other setting under the floor), on the §25.4 baseline: the law's LF 0.1 DC
        ratio, DC_MIN, beta and gamma on a joint grid, the 25-1000 Hz and 100 Hz-10 kHz costs per setting;
floorscore  each setting with the constants chosen on LF 0.1 leaving its RT step out (mean of the two relative costs);
cellscore / cellreport  the same with the constants chosen per (LF, RT) cell on its other HF steps (a per-step table);
floorreport  the choices per fold and the held-out nulls -> analysis/lowband/<kernel>_shelf_floor_report.log.
perloop  LF 0.1 and LF 0.4 / RT <= 1.0 on the §25.4 baseline: each loop's DC and Nyquist gain changed freely (8 numbers),
        held out by capture: the unit's DC loss per loop against the law's (perloopreport).
hfcell  HF < 1 on the §25.6 baseline: the HF pole shifted on a grid (100 Hz-10 kHz cost and nulls); hfcellreport chooses
        it per (HF, RT) cell on the cell's other LF steps.
bbcell / bbscore / bbreport  all settings on the §25.7 baseline: the broadband rate and the shelf's DC end per (LF, RT)
        cell, chosen on the cell's other HF steps.
octaves / octreport  all settings on the §25.7 baseline: excess over repeat per octave band, window and control step.
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
SOSB = butter(3, (100, 10000), btype="band", fs=FS, output="sos")  # floor stage: the broadband guard
PBAR = float(np.mean(q1.P))
GRID = np.round(np.arange(0.405, 0.43001, 0.00025), 5)               # dc stage: LF 0.4's DC ratio
LF_DC = 0.4
WS, BANDS = ("150-300", "300-600", "600-1000"), ("25-50", "50-100", "100-200", "100-10k")


def outdir(stage, keep=False):
    nyq = stage.startswith(("floor", "cell")) and not GAMMA_DC
    cap = f"_cap{os.environ['EMT_PMAX']}" if os.environ.get("EMT_PMAX") else ""        # 25.18: a refit around another HF cap
    return os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_{stage}" + ("K" if keep else "") + ("N" if nyq else "") + cap)


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


def prep(key, terms=None, lfdc=None):
    """The baseline for one setting: loop responses, excitation spectra, cluster, and every capture fitted (clock,
    shift, gain, fade) as lowband.py. terms: a replacement for lowband.fold_terms; lfdc: {lf: law DC ratio}."""
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, lb.lf_.AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    v, law, _ = q2.final_injection(compact=lb.lf_.COMPACT)
    for lf_k, x in (lfdc or {}).items():
        law["lf"][float(lf_k)]["dc"] = float(x)
    q1.EXT = 1700
    lab = dict(lf=info[key]["lf"], rt=info[key]["rt"], hf=info[key]["hf"])
    b230 = ap["steps"][ap["settings"][key]["step"]]
    keep = []
    ft = (terms or lb.fold_terms)(lab["rt"])
    G = q2.final_G(lab, b230, law, **ft, Gmod=lambda l_, w_, ww, Gs: keep.append(Gs) or Gs)
    caps = []
    for amp, nm in q2.files_of(info[key]):
        y, mdl, g, te = q2.fit_capture(G, v, amp, nm, fade[nm], eps0[nm], ends[nm], 1700)
        n = len(y)
        i0, i1 = PRE + int(0.15 * FS), min(PRE + int(te * FS), n)
        caps.append(dict(amp=amp, y=y, mdl=mdl, g=g, te=te, norm=float(np.sqrt(np.sum(sosfiltfilt(SOSF, y, axis=0)[i0:i1] ** 2))),
                         normB=float(np.sqrt(np.sum(sosfiltfilt(SOSB, y, axis=0)[i0:i1] ** 2)))))
    N = G["N"]
    J = q1.Joint([dict(G, files=[dict(y=np.zeros((N, 2)), w=np.ones(N), amp=0)])])
    cl, f, xs = J.split(v)
    F, X = J._spectra(f, xs, N)
    return dict(key=key, lab=lab, law=law, G=G, N=N, w=2 * np.pi * np.fft.rfftfreq(N), F=F, X=X, Gs0=keep[0],
                clus=J.cluster(cl), caps=caps, r0=base_r0(lab, law), terms=ft, b230=b230)


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


def residual(d, Gs, caps, t0=0.15, sos=SOSF, norm="norm", r=None):
    r = reverb(d, Gs) if r is None else r
    out = []
    for c in caps:
        n = len(c["y"])
        i0, i1 = PRE + int(t0 * FS), min(PRE + int(c["te"] * FS), n)
        e = sosfiltfilt(sos, c["y"] - model(d, c, r) * c["g"][:n], axis=0)[i0:i1]      # g is (n, 2)
        out.append(e.ravel() / c[norm])
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


WORKERS = 1


def run(stage, only, keep):
    info = q2.q1_tables()[0]
    od = outdir(stage, keep)
    os.makedirs(od, exist_ok=True)
    sel = {"free": lambda g: True, "dc": lambda g: g["lf"] == LF_DC,
           "floor": lambda g: g["lf"] == 0.1 or (g["lf"] == 0.4 and g["rt"] == 0.3),
           "perloop": lambda g: g["lf"] == 0.1 or (g["lf"] == 0.4 and g["rt"] <= 1.0),
           "hfcell": lambda g: g["hf"] < 1.0, "bbcell": lambda g: True, "octaves": lambda g: True, "hfbb": lambda g: g["hf"] < 1.0}
    sel["bbscore"] = sel["bbcell"]
    sel["floorscore"] = sel["cellscore"] = sel["floor"]
    keys = sorted(k for k in info if (only is None or k in only) and sel[stage](info[k]))
    keys = [k for k in keys if not os.path.exists(os.path.join(od, k + ".json"))]
    print(f"{stage}: {len(keys)} settings to do", flush=True)
    if WORKERS > 1 and stage in ("hfbb", "hfcell"):                            # one process, a pool of workers (~0.6 GB each)
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(WORKERS, max_tasks_per_child=4) as ex:
            for i, fu in enumerate(as_completed([ex.submit(_run_one, (stage, k, os.path.join(od, k + ".json"), keep)) for k in keys])):
                print(f"  {i + 1}/{len(keys)} {fu.result()}", flush=True)
        return
    for i, k in enumerate(keys):                                         # one process, ~0.6 GB (measured)
        path = os.path.join(od, k + ".json")
        if stage == "dc":
            dc(k, path, keep)
        else:
            {"free": free, "floor": floor, "floorscore": floorscore, "cellscore": cellscore, "perloop": perloop, "hfcell": hfcell, "bbcell": bbcell, "bbscore": bbscore, "octaves": octaves, "hfbb": hfbb}[stage](k, path)
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


# ---------------------------------------------------------------- the DC floor and LF 0.1's DC ratio (model_s_results 25.5)
CV04 = os.environ.get("EMT_CV_FOLDS") or os.path.join(ROOT, "analysis", "loop_cv", "lf04dc")   # the §25.4 baseline's folds (EMT_CV_FOLDS: another set, 25.18)
FL_DC = np.round(np.arange(0.1300, 0.13601, 0.0002), 5)               # LF 0.1's DC ratio (law 0.1317)
FL_DCMIN = np.round(np.arange(0.170, 0.2001, 0.005), 4)               # DC_MIN (0.18)
FL_BETA = (0.0, 0.01, 0.019, 0.03, 0.04)                              # beta (0.019)
FL_GAMMA = (-0.05, -0.04, -0.03, -0.02, -0.01, 0.0)                   # gamma (fold: -0.02 / -0.03)
FL_GRID = [(a, b, c, e) for a in FL_DC for b in FL_DCMIN for c in FL_BETA for e in FL_GAMMA]
GAMMA_DC = True                                                       # --gamma-nyq: gamma on the Nyquist end only (§25.6)


def cv_fold(rt):
    return json.load(open(os.path.join(CV04, f"fold_{rt}.json")))


def cv_terms(rt):
    """The §25.4 fold's loop terms for RT step rt, as final_G arguments."""
    f = cv_fold(rt)
    return dict(gamma=f["gamma"], lfrt=f["lfrt"], hfdp=f["hfdp"], pramp=tuple(f["pramp"]))


def floor_Gs(d, dc01, dc_min, beta, gamma):
    """The §25.4 loop responses with LF 0.1's DC ratio and the floor's constants replaced."""
    law = dict(d["law"], lf={**d["law"]["lf"], 0.1: dict(d["law"]["lf"][0.1], dc=float(dc01))})
    t = dict(d["terms"], gamma=gamma)
    return q2.final_Gs(d["lab"], law, d["w"], lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=d["b230"],
                       dc_min=dc_min, beta=beta, gamma_dc=GAMMA_DC)


def floored(lab, law, dc01, dc_min):
    """Per loop: is the DC floor acting (G_capped's test)?"""
    L = dict(law["lf"][lab["lf"]], **({"dc": dc01} if lab["lf"] == 0.1 else {}))
    out = []
    for Pi in q1.P:
        kk = -3 * Pi / q1.FS
        g = 10 ** (kk / (L["hi"] * lab["rt"]))
        out.append(L["dc"] != L["hi"] and g * 10 ** (kk * (1 / (L["dc"] * lab["rt"]) - 1 / (L["hi"] * lab["rt"]))) < dc_min)
    return out


def floor(key, path):
    """Per setting: the 25-1000 Hz and 100 Hz-10 kHz costs (150 ms to the compensated end, all captures) over FL_GRID.
    Grid points that leave every loop's response unchanged share one evaluation (key: per loop, floored -> the floor's
    constants and the DC ratio, which sets the floor's lift and with it beta's share; else the DC ratio only)."""
    rt = q2.q1_tables()[0][key]["rt"]
    d = prep(key, cv_terms, cv_fold(rt)["lfdc"])
    lab, law = d["lab"], d["law"]
    G0 = floor_Gs(d, law["lf"][0.1]["dc"], q2.DC_MIN, q2.FINAL_BETA, d["terms"]["gamma"])
    chk = max(float(np.max(np.abs(a - b))) for a, b in zip(G0, d["Gs0"]))           # 0 unless --gamma-nyq (then the change itself)
    r0 = reverb(d, d["Gs0"])
    base = [float(np.sum(residual(d, None, d["caps"], r=r0) ** 2)), float(np.sum(residual(d, None, d["caps"], sos=SOSB, norm="normB", r=r0) ** 2))]
    cache, rows = {}, []
    for dc01, dc_min, beta, gamma in FL_GRID:
        fl = floored(lab, law, dc01, dc_min)
        dk = (dc01,) if lab["lf"] == 0.1 else ()
        k = tuple(dk + (dc_min, beta, gamma) if f else dk for f in fl)
        if k not in cache:
            r = reverb(d, floor_Gs(d, dc01, dc_min, beta, gamma))
            cache[k] = [float(np.sum(residual(d, None, d["caps"], r=r) ** 2)), float(np.sum(residual(d, None, d["caps"], sos=SOSB, norm="normB", r=r) ** 2))]
        rows.append(cache[k])
    json.dump(dict(key=key, **lab, check=chk, base=base, n_eval=len(cache), grid=[list(x) for x in FL_GRID], cost=rows), open(path, "w"))


def floor_choice(D, train):
    """Over the training settings, mean dB change of each cost re the §25.4 baseline; among the grid points that do not
    worsen the broadband (mean <= 0 dB), the one with the lowest low-band change. Returns it and (low, broad) dB."""
    u = [d for d in D if d["rt"] in train]
    low = np.mean([10 * np.log10(np.array(d["cost"])[:, 0] / d["base"][0]) for d in u], 0)
    brd = np.mean([10 * np.log10(np.array(d["cost"])[:, 1] / d["base"][1]) for d in u], 0)
    ok = brd <= 1e-9
    i = int(np.argmin(np.where(ok, low, np.inf))) if ok.any() else int(np.argmin(brd))
    return FL_GRID[i], (float(low[i]), float(brd[i]))


def floor_load():
    return [json.load(open(p)) for p in glob.glob(os.path.join(outdir("floor"), "*.json"))]


def floorscore(key, path):
    """Per setting: band nulls of every capture with the floor constants chosen leaving its RT step out (on LF 0.1)."""
    D = [d for d in floor_load() if d["lf"] == 0.1]
    rt = q2.q1_tables()[0][key]["rt"]
    (dc01, dc_min, beta, gamma), _ = floor_choice(D, [r for r in sorted({d["rt"] for d in D}) if r != rt])
    d = prep(key, cv_terms, cv_fold(rt)["lfdc"])
    r = reverb(d, floor_Gs(d, dc01, dc_min, beta, gamma))
    json.dump(dict(key=key, **d["lab"], chosen=[dc01, dc_min, beta, gamma], amps=[c["amp"] for c in d["caps"]],
                   base=[lb.band_nulls(c["y"], c["mdl"] * c["g"], c["te"]) for c in d["caps"]],
                   held=[lb.band_nulls(c["y"], model(d, c, r) * c["g"], c["te"]) for c in d["caps"]]), open(path, "w"))


def cell_choice(D, key):
    """Per (LF, RT) cell: the floor_choice rule over the cell's other settings (held out by HF step)."""
    me = next(d for d in D if d["key"] == key)
    return floor_choice([d for d in D if d["lf"] == me["lf"] and d["rt"] == me["rt"] and d["key"] != key], [me["rt"]])


def cellscore(key, path):
    """Per setting: band nulls of every capture with the floor constants of its (LF, RT) cell chosen on the cell's other
    HF steps (a per-step table, model_s_results 25.6)."""
    (dc01, dc_min, beta, gamma), _ = cell_choice(floor_load(), key)
    rt = q2.q1_tables()[0][key]["rt"]
    d = prep(key, cv_terms, cv_fold(rt)["lfdc"])
    r = reverb(d, floor_Gs(d, dc01, dc_min, beta, gamma))
    json.dump(dict(key=key, **d["lab"], chosen=[dc01, dc_min, beta, gamma], amps=[c["amp"] for c in d["caps"]],
                   base=[lb.band_nulls(c["y"], c["mdl"] * c["g"], c["te"]) for c in d["caps"]],
                   held=[lb.band_nulls(c["y"], model(d, c, r) * c["g"], c["te"]) for c in d["caps"]]), open(path, "w"))


def report_floor(score="floorscore"):
    D = floor_load()
    D1 = [d for d in D if d["lf"] == 0.1]
    rts = sorted({d["rt"] for d in D1})
    out = [f"## The DC floor and LF 0.1's DC ratio ({len(D)} settings; grid {len(FL_GRID)} points; law dc 0.1317, DC_MIN 0.18, beta 0.019; "
           f"gamma {'on the whole loop' if GAMMA_DC else 'on the Nyquist end only'})",
           f"render check (grid point at the baseline's constants vs the baseline): max |dG| {max(d['check'] for d in D):.1e}"]
    for rt in rts + ["all"]:
        tr = [r for r in rts if r != rt]
        (a, b, c, e), (dl, db) = floor_choice(D1, tr)
        own = [floor_choice([d], [d["rt"]])[0] for d in D1 if d["rt"] == rt]
        out.append(f"- {'all RT steps' if rt == 'all' else f'hold RT {rt}'}: dc {a:.4f}, DC_MIN {b:.3f}, beta {c:g}, gamma {e:+.2f} "
                   f"(training mean cost change: low {dl:+.2f} dB, broad {db:+.2f} dB)" +
                   ("" if rt == "all" else "; own best per setting: " + ", ".join(f"({x[0]:.4f}, {x[1]:.3f}, {x[2]:g}, {x[3]:+.2f})" for x in own)))
    S = [json.load(open(p)) for p in glob.glob(os.path.join(outdir(score), "*.json"))]
    if not S:
        return out
    rep = {(k, r["amp"]): r["repeat"] for k in {s_["key"] for s_ in S} for r in json.load(open(os.path.join(lb.outdir("base"), k + ".json")))}

    def cell(ss, b, w, f, excess):
        x = [s_[f][j][b][w][1][0] - (rep[(s_["key"], a)][b][w][1][0] if excess else 0) for s_ in ss for j, a in enumerate(s_["amps"])
             if w in s_[f][j][b] and (not excess or w in rep[(s_["key"], a)][b])]
        return np.median(x) if x else np.nan
    out.append("\n### Held out by " + ("RT step" if score == "floorscore" else "HF step within each (LF, RT) cell") +
               ": R null, median over captures, §25.4 baseline -> chosen constants [excess over repeat]")
    groups = [("LF 0.1", lambda s_: s_["lf"] == 0.1)] + [(f"LF 0.1 RT {rt}", lambda s_, rt=rt: s_["lf"] == 0.1 and s_["rt"] == rt) for rt in rts] + \
             [("LF 0.4 RT 0.3", lambda s_: s_["lf"] == 0.4)]
    for b in BANDS:
        out.append(f"#### {b} Hz, " + " | ".join(WS))
        for g, sel in groups:
            ss = [s_ for s_ in S if sel(s_)]
            if ss:
                out.append(f"  {g:14s} " + " | ".join(f"{cell(ss, b, w, 'base', False):.1f} -> {cell(ss, b, w, 'held', False):.1f} "
                                                      f"[{cell(ss, b, w, 'base', True):+.1f} -> {cell(ss, b, w, 'held', True):+.1f}]" for w in WS))
    return out


# ---------------------------------------------------------------- per-loop DC and Nyquist gains (model_s_results 25.6)
def law_dc_loss(lab, law):
    """Per loop: the law's DC loss per pass in dB (before the floor) and the broadband (Nyquist-end) loss."""
    L = law["lf"][lab["lf"]]
    return ([20 * 3 * Pi / (q1.FS * L["dc"] * lab["rt"]) for Pi in q1.P], [20 * 3 * Pi / (q1.FS * L["hi"] * lab["rt"]) for Pi in q1.P])


def perloop(key, path):
    """Per setting on the §25.4 baseline: each loop's response x a first-order section (the shelf's pole) changing its DC
    gain by x_i dB and its Nyquist gain by y_i dB (8 numbers), fitted on the 25-1000 Hz residual, held out by capture."""
    rt = q2.q1_tables()[0][key]["rt"]
    d = prep(key, cv_terms, cv_fold(rt)["lfdc"])
    a = d["law"]["a"]
    sh = lambda x, y: [m.shelf(d["w"], 10 ** (x[i] / 20), 10 ** (y[i] / 20), a) for i in range(4)]
    Gs = lambda th: [Gi * c for Gi, c in zip(d["Gs0"], sh(th[:4], th[4:]))]
    res = {}
    for ho, c in enumerate(d["caps"]):
        tr = [x for j, x in enumerate(d["caps"]) if j != ho]
        fit = least_squares(lambda th: residual(d, Gs(th), tr), np.zeros(8), x_scale=0.05, diff_step=1e-3, max_nfev=80)
        r = reverb(d, Gs(fit.x))
        res[c["amp"]] = dict(th=fit.x.tolist(), base=lb.band_nulls(c["y"], c["mdl"] * c["g"], c["te"]),
                             fit=lb.band_nulls(c["y"], model(d, c, r) * c["g"], c["te"]))
    ldc, lhi = law_dc_loss(d["lab"], d["law"])
    json.dump(dict(key=key, **d["lab"], law_dc_db=ldc, law_hi_db=lhi, base_dc_db=[float(-20 * np.log10(abs(G[0]))) for G in d["Gs0"]],
                   base_ny_db=[float(-20 * np.log10(abs(G[-1]))) for G in d["Gs0"]], res=res), open(path, "w"))


def report_perloop():
    D = [json.load(open(p)) for p in glob.glob(os.path.join(outdir("perloop"), "*.json"))]
    out = ["## Per-loop DC / Nyquist changes on the §25.4 baseline, held out by capture (8 numbers per setting); per loop: "
           "law DC loss -> baseline (floor, gamma) -> unit (baseline - x), dB per pass; Nyquist change y; median over HF and folds"]
    rows = []
    for lf in sorted({d["lf"] for d in D}):
        for rt in sorted({d["rt"] for d in D if d["lf"] == lf}):
            u = [d for d in D if d["lf"] == lf and d["rt"] == rt]
            th = np.median([r["th"] for d in u for r in d["res"].values()], 0)
            ldc = np.median([d["law_dc_db"] for d in u], 0)
            bdc = np.median([d["base_dc_db"] for d in u], 0)
            cells = [f"{ldc[i]:.2f} -> {bdc[i]:.2f} -> {bdc[i] - th[i]:.2f} (y {th[4 + i]:+.3f})" for i in range(4)]
            g = lambda f, b, w: np.median([r[f][b][w][1][0] for d in u for r in d["res"].values() if w in r[f][b]])
            out.append(f"- LF {lf} RT {rt}: " + "; ".join(cells) + f" | R 100-200 Hz 150-300 {g('base', '100-200', '150-300'):.1f} -> "
                       f"{g('fit', '100-200', '150-300'):.1f}, 100-10k {g('base', '100-10k', '150-300'):.1f} -> {g('fit', '100-10k', '150-300'):.1f}")
            rows += [(lf, rt, ldc[i], bdc[i] - th[i]) for i in range(4)]
    out.append("\n### Unit DC loss against the law's, all loops (sorted by law): law -> unit [LF, RT]")
    out.append("  " + "; ".join(f"{l:.1f}->{t:.2f} [{lf:g},{rt:g}]" for lf, rt, l, t in sorted(rows, key=lambda r: r[2])))
    return out


# ---------------------------------------------------------------- the HF pole per (HF, RT) cell (model_s_results 25.7)
CELLS_JSON = os.path.join(ROOT, "analysis", "loop_cv", "cells", "cells.json")   # the §25.6 floor table
HF_DP = np.round(np.arange(-0.012, 0.01201, 0.0005), 4)                          # absolute HF pole shift on every loop


HF_CELLS_JSON = os.path.join(ROOT, "analysis", "loop_cv", "hfcells", "cells.json")   # the §25.7 floor and HF tables


def base_kw(key, table=CELLS_JSON):
    """The baseline's final_Gs arguments for one setting (fold terms, LF 0.4 dc, the setting's table entry: §25.6's floor
    cells, or with table=HF_CELLS_JSON §25.7's floor and HF cells)."""
    rt = q2.q1_tables()[0][key]["rt"]
    t = cv_terms(rt)
    fc = json.load(open(table)).get(key, {})
    lfdc = dict(cv_fold(rt)["lfdc"])
    kw = {}
    if "dc_min" in fc:
        lfdc["0.1"] = fc["dc"]
        t["gamma"] = fc["gamma"]
        kw = dict(dc_min=fc["dc_min"], beta=fc["beta"], gamma_dc=fc.get("gamma_dc", True))
    if "hf_dp" in fc:
        t["hfdp"] = {str(q2.q1_tables()[0][key]["hf"]): fc["hf_dp"], "p_above": 0.0, "p_below": 1.0, "exclude_lf": []}
    return t, lfdc, kw


def nulls_broad(y, m, te):
    """100 Hz-10 kHz strict nulls per window and channel (band_nulls' 100-10k part only)."""
    yf, mf = sosfiltfilt(SOSB, y, axis=0), sosfiltfilt(SOSB, m, axis=0)
    d = {}
    for a_ms, b_ms in q2.W_COMP:
        i0, i1 = PRE + int(a_ms / 1000 * FS), min(PRE + int(b_ms / 1000 * FS), PRE + int(te * FS), len(y))
        if i1 - i0 >= int(0.05 * FS):
            d[f"{a_ms}-{b_ms}"] = [float(10 * np.log10(np.sum((yf[i0:i1, c] - mf[i0:i1, c]) ** 2) / np.sum(yf[i0:i1, c] ** 2))) for c in (0, 1)]
    return d


def hfcell(key, path):
    """Per setting (HF < 1): the §25.6 baseline with the HF one-pole's pole shifted by dp on every loop (replacing the
    fold's hfdp), on HF_DP: the 100 Hz-10 kHz cost and nulls of every capture."""
    t, lfdc, kw = base_kw(key)
    d = prep(key, lambda rt: dict(t), lfdc)
    G0 = q2.final_Gs(d["lab"], d["law"], d["w"], lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=d["b230"], **kw)
    d["Gs0"] = G0                                          # the floor cell is not in prep's final_G call
    r0 = reverb(d, G0)
    hf = str(d["lab"]["hf"])
    base = float(np.sum(residual(d, None, d["caps"], sos=SOSB, norm="normB", r=r0) ** 2))
    rows = []
    for dp in HF_DP:
        Gs = q2.final_Gs(d["lab"], d["law"], d["w"], lfrt=t["lfrt"], gamma=t["gamma"], pramp=t["pramp"], b230=d["b230"],
                         hfdp={hf: float(dp), "p_above": 0.0, "p_below": 1.0, "exclude_lf": []}, **kw)
        r = reverb(d, Gs)
        rows.append(dict(dp=float(dp), cost=float(np.sum(residual(d, None, d["caps"], sos=SOSB, norm="normB", r=r) ** 2)),
                         nulls=[nulls_broad(c["y"], model(d, c, r) * c["g"], c["te"]) for c in d["caps"]]))
    json.dump(dict(key=key, **d["lab"], base=base, base_nulls=[nulls_broad(c["y"], model(d, c, r0) * c["g"], c["te"]) for c in d["caps"]],
                   amps=[c["amp"] for c in d["caps"]], grid=rows), open(path, "w"))


def report_hfcell():
    D = [json.load(open(p)) for p in glob.glob(os.path.join(outdir("hfcell"), "*.json"))]
    rep = {(k, r["amp"]): r["repeat"]["100-10k"] for k in {d["key"] for d in D} for r in json.load(open(os.path.join(lb.outdir("base"), k + ".json")))}
    rel = {d["key"]: np.array([10 * np.log10(g["cost"] / d["base"]) for g in d["grid"]]) for d in D}
    out = [f"## HF pole shift per (HF, RT) cell on the §25.6 baseline ({len(D)} settings; dp {HF_DP[0]:+.4f}..{HF_DP[-1]:+.4f}); 100 Hz-10 kHz, R"]
    choice = {}
    for d in D:                                            # held out across LF: the cell's other LF steps
        oth = [rel[e["key"]] for e in D if e["hf"] == d["hf"] and e["rt"] == d["rt"] and e["key"] != d["key"]]
        choice[d["key"]] = int(np.argmin(np.mean(oth, 0)))
    for hf in sorted({d["hf"] for d in D}):
        for rt in sorted({d["rt"] for d in D}):
            u = [d for d in D if d["hf"] == hf and d["rt"] == rt]
            own = [HF_DP[int(np.argmin(rel[d["key"]]))] for d in sorted(u, key=lambda d: d["lf"])]
            held = [HF_DP[choice[d["key"]]] for d in sorted(u, key=lambda d: d["lf"])]
            out.append(f"- HF {hf} RT {rt}: own best dp by LF 0.1/0.4/1.0/2.0 " + "/".join(f"{x:+.4f}" for x in own) +
                       "; held-out choice " + "/".join(f"{x:+.4f}" for x in held) +
                       "; held-out cost change " + "/".join(f"{rel[d['key']][choice[d['key']]]:+.2f}" for d in sorted(u, key=lambda d: d["lf"])) + " dB")

    def cell(ss, w, held, excess):
        x = []
        for d in ss:
            nn = d["grid"][choice[d["key"]]]["nulls"] if held else d["base_nulls"]
            for j, a in enumerate(d["amps"]):
                if w in nn[j] and (not excess or w in rep[(d["key"], a)]):
                    x.append(nn[j][w][1] - (rep[(d["key"], a)][w][1][0] if excess else 0))
        return np.median(x) if x else np.nan
    out.append("\n### Held out across LF steps within each (HF, RT) cell: R null, median over captures, §25.6 -> chosen dp [excess over repeat]")
    groups = [("all HF < 1", lambda d: True)] + [(f"HF {hf}", lambda d, hf=hf: d["hf"] == hf) for hf in sorted({d["hf"] for d in D})] + \
             [(f"HF 0.2 RT {rt}", lambda d, rt=rt: d["hf"] == 0.2 and d["rt"] == rt) for rt in sorted({d["rt"] for d in D})]
    for g, sel in groups:
        ss = [d for d in D if sel(d)]
        out.append(f"  {g:14s} " + " | ".join(f"{w}: {cell(ss, w, False, False):.1f} -> {cell(ss, w, True, False):.1f} "
                                              f"[{cell(ss, w, False, True):+.1f} -> {cell(ss, w, True, True):+.1f}]" for w in ("50-150", "150-300", "300-600", "600-1000")))
    return out


# ---------------------------------------------------------------- broadband rate and shelf DC per (LF, RT) cell (model_s_results 25.8)
BB_EPS = np.round(np.arange(-0.010, 0.01001, 0.001), 4)       # broadband decay rate x (1 + eps)
BB_DDC = np.round(np.arange(-0.03, 0.03001, 0.01), 3)         # DC per-pass gain change, dB x P_i / Pbar (a shelf with the law's pole)


def bb_prep(key):
    t, lfdc, kw = base_kw(key, HF_CELLS_JSON)
    d = prep(key, lambda rt: dict(t), lfdc)
    d["Gs0"] = q2.final_Gs(d["lab"], d["law"], d["w"], lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=d["b230"], **kw)
    return d


def bb_Gs(d, eps, ddc):
    rate = 60 / (d["law"]["lf"][d["lab"]["lf"]]["hi"] * d["lab"]["rt"])
    return [Gi * 10 ** (eps * rate * Pi / (20 * q1.FS)) * m.shelf(d["w"], 10 ** (ddc * Pi / PBAR / 20), 1.0, d["law"]["a"])
            for Gi, Pi in zip(d["Gs0"], q1.P)]


def bbcell(key, path):
    """Per setting on the §25.7 baseline: the broadband rate x (1 + eps) and the DC end's per-pass gain + ddc on a grid,
    the 25-1000 Hz and 100 Hz-10 kHz costs over all captures."""
    d = bb_prep(key)
    r0 = reverb(d, d["Gs0"])
    base = [float(np.sum(residual(d, None, d["caps"], r=r0) ** 2)), float(np.sum(residual(d, None, d["caps"], sos=SOSB, norm="normB", r=r0) ** 2))]
    rows = []
    for e in BB_EPS:
        for x in BB_DDC:
            r = reverb(d, bb_Gs(d, e, x))
            rows.append([float(e), float(x), float(np.sum(residual(d, None, d["caps"], r=r) ** 2)),
                         float(np.sum(residual(d, None, d["caps"], sos=SOSB, norm="normB", r=r) ** 2))])
    json.dump(dict(key=key, **d["lab"], base=base, grid=rows), open(path, "w"))


HFBB_EPS = np.round(np.arange(-0.15, 0.0501, 0.0125), 4)       # hfbb: broadband decay rate x (1 + eps) on the capped-pole loops


def hfbb_Gs(d, t1, kw, eps):
    """Candidate strong-HF loops (model_s_results 25.15): HF pole min(p_law, P_MAX), no ramp and no HF-cell shift (the
    fold's own law terms, floor cells and LF 0.4 dc kept), broadband rate x (1 + eps)."""
    G = q2.final_Gs(d["lab"], d["law"], d["w"], lfrt=t1["lfrt"], gamma=t1["gamma"], hfdp=t1["hfdp"], pramp=None, b230=d["b230"], **kw)
    rate = 60 / (d["law"]["lf"][d["lab"]["lf"]]["hi"] * d["lab"]["rt"])
    return [Gi * 10 ** (eps * rate * Pi / (20 * q1.FS)) for Gi, Pi in zip(G, q1.P)]


def hfbb(key, path):
    """Per setting (HF < 1): 100 Hz-10 kHz costs from 50 ms and from 150 ms, all captures, of the §25.7 baseline and of
    the candidate loops over the eps grid."""
    d = bb_prep(key)
    t, lfdc, kw = base_kw(key, HF_CELLS_JSON)
    t1 = dict(t, hfdp=cv_terms(d["lab"]["rt"])["hfdp"])               # the fold's law hfdp, not the HF cell's
    cost = lambda r: [float(np.sum(residual(d, None, d["caps"], t0=t0, sos=SOSB, norm="normB", r=r) ** 2)) for t0 in (0.05, 0.15)]
    base = cost(reverb(d, d["Gs0"]))
    rows = [[float(e)] + cost(reverb(d, hfbb_Gs(d, t1, kw, e))) for e in HFBB_EPS]
    json.dump(dict(key=key, **d["lab"], b230=d["b230"], base=base, grid=rows), open(path, "w"))


def _run_one(a):
    stage, k, path, keep = a
    {"hfbb": hfbb, "hfcell": hfcell}[stage](k, path)
    return k


def bb_choice(D, key):
    """The cell's other HF steps: lowest mean broadband change (dB) among the points that do not worsen the low band."""
    me = next(d for d in D if d["key"] == key)
    u = [d for d in D if d["lf"] == me["lf"] and d["rt"] == me["rt"] and d["key"] != key]
    lo = np.mean([10 * np.log10(np.array(d["grid"])[:, 2] / d["base"][0]) for d in u], 0)
    br = np.mean([10 * np.log10(np.array(d["grid"])[:, 3] / d["base"][1]) for d in u], 0)
    ok = lo <= 1e-9
    i = int(np.argmin(np.where(ok, br, np.inf))) if ok.any() else int(np.argmin(br))
    return me["grid"][i][:2], i


def bbscore(key, path):
    D = [json.load(open(p)) for p in glob.glob(os.path.join(outdir("bbcell"), "*.json"))]
    (e, x), _ = bb_choice(D, key)
    d = bb_prep(key)
    r0, r = reverb(d, d["Gs0"]), reverb(d, bb_Gs(d, e, x))
    json.dump(dict(key=key, **d["lab"], chosen=[e, x], amps=[c["amp"] for c in d["caps"]],
                   base=[lb.band_nulls(c["y"], model(d, c, r0) * c["g"], c["te"]) for c in d["caps"]],
                   held=[lb.band_nulls(c["y"], model(d, c, r) * c["g"], c["te"]) for c in d["caps"]]), open(path, "w"))


def report_bb():
    D = [json.load(open(p)) for p in glob.glob(os.path.join(outdir("bbcell"), "*.json"))]
    S = {s_["key"]: s_ for s_ in (json.load(open(p)) for p in glob.glob(os.path.join(outdir("bbscore"), "*.json")))}
    out = [f"## Broadband rate and shelf DC per (LF, RT) cell on the §25.7 baseline ({len(D)} settings; eps {BB_EPS[0]:+.3f}..{BB_EPS[-1]:+.3f}, "
           f"ddc {BB_DDC[0]:+.2f}..{BB_DDC[-1]:+.2f} dB/pass); choice on the cell's other HF steps"]
    for lf in sorted({d["lf"] for d in D}):
        for rt in sorted({d["rt"] for d in D}):
            u = sorted([d for d in D if d["lf"] == lf and d["rt"] == rt], key=lambda d: d["hf"])
            own = [d["grid"][int(np.argmin(np.array(d["grid"])[:, 3]))][:2] for d in u]
            held = [bb_choice(D, d["key"]) for d in u]
            ch = [10 * np.log10(d["grid"][i][3] / d["base"][1]) for d, (_, i) in zip(u, held)]
            out.append(f"- LF {lf} RT {rt}: own best (eps, ddc) by HF " + " ".join(f"({a:+.3f},{b:+.2f})" for a, b in own) +
                       "; held-out " + " ".join(f"({a:+.3f},{b:+.2f})" for (a, b), _ in held) + "; broadband change " + "/".join(f"{c:+.2f}" for c in ch))
    if not S:
        return out
    rep = {(k, r["amp"]): r["repeat"] for k in S for r in json.load(open(os.path.join(lb.outdir("base"), k + ".json")))}

    def cell(ss, b, w, f, excess):
        x = [s_[f][j][b][w][1][0] - (rep[(s_["key"], a)][b][w][1][0] if excess else 0) for s_ in ss for j, a in enumerate(s_["amps"])
             if w in s_[f][j][b] and (not excess or w in rep[(s_["key"], a)][b])]
        return np.median(x) if x else np.nan
    out.append("\n### Held out by HF step within each (LF, RT) cell: R null, median over captures, §25.7 -> chosen [excess over repeat]")
    groups = [("all", lambda s_: True)] + [(f"LF {lf}", lambda s_, lf=lf: s_["lf"] == lf) for lf in (0.1, 0.4, 1.0, 2.0)] + \
             [("flat", lambda s_: s_["lf"] == 1.0 and s_["hf"] == 1.0), ("RT <= 1.0", lambda s_: s_["rt"] <= 1.0), ("RT >= 1.4", lambda s_: s_["rt"] >= 1.4)]
    for b in ("100-200", "100-10k"):
        out.append(f"#### {b} Hz, " + " | ".join(WS))
        for g, sel in groups:
            ss = [s_ for s_ in S.values() if sel(s_)]
            out.append(f"  {g:10s} " + " | ".join(f"{cell(ss, b, w, 'base', False):.1f} -> {cell(ss, b, w, 'held', False):.1f} "
                                                  f"[{cell(ss, b, w, 'base', True):+.1f} -> {cell(ss, b, w, 'held', True):+.1f}]" for w in WS))
    return out


# ---------------------------------------------------------------- where the §25.7 baseline's excess lives (model_s_results 25.9)
OCT = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
SOS_OCT = {f: butter(3, (f / 2 ** 0.5, min(f * 2 ** 0.5, 11000)), btype="band", fs=FS, output="sos") for f in OCT}
W_OCT = [(50, 150), (150, 300), (300, 600), (600, 1000)]


def octaves(key, path):
    """Per setting on the §25.7 baseline, per capture (R): null per octave band and window, the same against the next
    AMPL capture (repeat), and the band's level re 1 kHz."""
    d = bb_prep(key)
    r = reverb(d, d["Gs0"])
    caps, rows = d["caps"], []
    for j, c in enumerate(caps):
        cb = caps[(j + 1) % len(caps)]
        n = min(len(c["y"]), len(cb["y"]))
        mm = model(d, c, r)[:n] * c["g"][:n]
        te, res = min(c["te"], cb["te"]), {}
        for f, sos in SOS_OCT.items():
            y, m_, yb = (sosfiltfilt(sos, x[:n], axis=0)[:, 1] for x in (c["y"], mm, cb["y"]))
            for a, b in W_OCT:
                i0, i1 = PRE + int(a / 1000 * FS), min(PRE + int(b / 1000 * FS), PRE + int(te * FS), n)
                if i1 - i0 >= int(0.05 * FS):
                    e = lambda u, v: float(10 * np.log10(np.sum((u - v) ** 2) / np.sum(u ** 2)))
                    res[f"{f}_{a}-{b}"] = [e(y[i0:i1], m_[i0:i1]), e(y[i0:i1], yb[i0:i1]), float(10 * np.log10(np.sum(y[i0:i1] ** 2)))]
        rows.append(res)
    json.dump(dict(key=key, **d["lab"], rows=rows), open(path, "w"))


def report_octaves():
    D = [json.load(open(p)) for p in glob.glob(os.path.join(outdir("octaves"), "*.json"))]
    out = [f"## Where the §25.7 baseline's excess lives ({len(D)} settings): R, per octave and window, median over captures: "
           "excess over repeat [band level re 1 kHz]"]
    groups = [("all", lambda d: True), ("flat", lambda d: d["lf"] == 1.0 and d["hf"] == 1.0)] + \
             [(f"LF {x}", lambda d, x=x: d["lf"] == x) for x in (0.1, 0.4, 1.0, 2.0)] + [(f"HF {x}", lambda d, x=x: d["hf"] == x) for x in (0.2, 0.4, 0.6, 1.0)] + \
             [("RT <= 1.0", lambda d: d["rt"] <= 1.0), ("RT >= 1.4", lambda d: d["rt"] >= 1.4)]
    for a, b in W_OCT:
        w = f"{a}-{b}"
        out.append(f"\n### {w} ms\n| group | " + " | ".join(f"{f}" for f in OCT) + " |\n|---|" + "---|" * len(OCT))
        for g, sel in groups:
            cells = []
            for f in OCT:
                k = f"{f}_{w}"
                v = np.array([r[k] for d in D if sel(d) for r in d["rows"] if k in r and f"1000_{w}" in r])
                lv = np.array([r[k][2] - r[f"1000_{w}"][2] for d in D if sel(d) for r in d["rows"] if k in r and f"1000_{w}" in r])
                cells.append(f"{np.median(v[:, 0] - v[:, 1]):+.1f} [{np.median(lv):+.0f}]" if len(v) else "-")
            out.append(f"| {g} | " + " | ".join(cells) + " |")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["free", "dc", "floor", "floorscore", "cellscore", "perloop", "hfcell", "bbcell", "bbscore", "octaves", "report", "octreport", "floorreport", "cellreport", "perloopreport", "hfcellreport", "bbreport", "hfbb"], required=True)
    ap.add_argument("--workers", type=int, default=1, help="hfbb: a pool of workers in this one process")
    ap.add_argument("--only", default=None)
    ap.add_argument("--keep-lfrt", action="store_true")
    ap.add_argument("--gamma-nyq", action="store_true", help="floor stages: gamma on the Nyquist end only, DC held at DC_MIN")
    a = ap.parse_args()
    global GAMMA_DC, WORKERS
    GAMMA_DC, WORKERS = not a.gamma_nyq, a.workers
    if a.stage == "octreport":
        out = report_octaves()
        open(os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_octaves_report.log"), "w").write("\n".join(out) + "\n")
        print("\n".join(out))
    elif a.stage == "bbreport":
        out = report_bb()
        open(os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_bb_report.log"), "w").write("\n".join(out) + "\n")
        print("\n".join(out))
    elif a.stage == "hfcellreport":
        out = report_hfcell()
        open(os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_hfcell_report.log"), "w").write("\n".join(out) + "\n")
        print("\n".join(out))
    elif a.stage == "perloopreport":
        out = report_perloop()
        open(os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_perloop_report.log"), "w").write("\n".join(out) + "\n")
        print("\n".join(out))
    elif a.stage in ("floorreport", "cellreport"):
        out = report_floor("floorscore" if a.stage == "floorreport" else "cellscore")
        name = "floor" if a.stage == "floorreport" else "cell"
        open(os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_{name}{'' if GAMMA_DC else 'N'}_report.log"), "w").write("\n".join(out) + "\n")
        print("\n".join(out))
    elif a.stage == "report":
        out = report_free() + [""] + report_dc(False) + [""] + report_dc(True)
        open(os.path.join(ROOT, "analysis", "lowband", f"{q1.ms.KTAG or 'library'}_shelf_report.log"), "w").write("\n".join(out) + "\n")
        print("\n".join(out))
    else:
        run(a.stage, a.only.split(",") if a.only else None, a.keep_lfrt)


if __name__ == "__main__":
    main()
