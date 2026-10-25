"""NEXT_STEPS §5 item 5, S's recirculations at strong settings: what shape do the loops need where the laws fail
(LF 0.1, HF 0.2, short RT)? Tomkin 96 kHz library, EMT_KERNEL=chain2. Results: docs/model_s_results.md §22.

    EMT_KERNEL=chain2 .venv/bin/python tools/loop_fir.py --stage fit [--only keys] [--workers 4]
    EMT_KERNEL=chain2 .venv/bin/python tools/loop_fir.py --stage report

The injection is frozen: the compact model refitted on chain2 (q2_injection --compact, model_s_results §21). Per
setting, with the clock, shift and gain of each capture fitted once against the capped law, three loop forms are
fitted on the AMPL 30 + 90 captures over the whole window from 50 ms to each capture's fade (time weight
10^(1.5 t / RT) as q2's loops stage), and scored on the AMPL 0 + 60 captures (held out by capture; fade-compensated
nulls as the final stage):

  law    the capped law of the final model (HF pole <= 0.84, per-pass DC gain >= 0.18, beta), nothing fitted;
  five   q1's 5-number loop form (loop RT, LF shelf DC and Nyquist RT ratios and pole, HF lowpass ratio), shared by
         the four loops and length-scaled, refitted (q2's `filters` / `loops` stages fitted it on the first 4096
         samples only, about 50-150 ms: Joint's dummy 4096-sample capture capped every window);
  fir    the capped law times a free correction per loop, C_i(z) = sum_k c_ik z^-k, k = -KC..KC (C_i = 1 is the law):
         G_i -> G_i C_i. 4 (2 KC + 1) numbers, exact Jacobian.

floor: the DC floor's constants (DC_MIN, beta) and a constant extra loss under the floor (gamma) on full windows, a grid
on the LF 0.1 / 0.4 settings, chosen on TRAIN_RT, reported held out by setting -> <kernel>_floor.json, _floor_best.json.

hfcap: HF 0.2 / 0.4 at RT <= 1.4, the HF pole's cap and an offset on a grid, chosen leaving one RT step out -> <kernel>_hfcap.json.
probe: which part of a setting's own correction carries its gain (common, magnitude, phase, per-loop) -> <kernel>_probe.json.
tilt: per LF step the training settings' median common magnitude shape (dB/s) -> <kernel>_tilt.json (q2 final --tilt).

shelf: on top of gamma, one shift of the LF shelf's pole (all LF steps) and each LF step's top end at fixed DC, a grid on
the LF 0.1 / 0.4 / 2.0 settings, chosen on TRAIN_RT -> <kernel>_shelf.json, _shelf_best.json.

interp: held out by setting, the corrections of RT 0.6 / 1.0 / 1.8 predicted by interpolating each loop's correction
linearly in RT from the neighbouring steps (same LF and HF) -> analysis/loop_fir/<kernel>_interp.json.

One result file per setting in analysis/loop_fir/<kernel>/ (skipped on rerun); report reads them all: nulls and excess
over the repeat-capture null (nulls_free_<kernel>), and the corrections' shape per pass: magnitude (dB) and phase delay
(samples) against frequency, per loop and by control. coefs: each correction read as changes of the law's own
coefficients -> analysis/loop_fir/<kernel>_coefs.json.
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import q2_injection as q2  # noqa: E402

q1 = q2.q1
FS, PRE = q2.FS, q2.PRE
KC = int(os.environ.get("LOOP_FIR_KC", 4))          # --kc (via the environment, so the pool's workers see it)
COMPACT = "diffuser_tied_fix-0.707107_ap230new_chain2"
AP230 = os.environ.get("EMT_AP230", "ap230_g0.7071_chain2.json")   # S's AP230 table in analysis/q2/; EMT_AP230=ap230_mltable.json: the damped
                                                                   # steps on the M / L table (model_s_results 38, tried, not adopted)
OUT = os.path.join(ROOT, "analysis", "loop_fir", q1.ms.KTAG or "library")
FREQS = [63, 125, 250, 500, 1000, 2000, 4000, 6600, 8000, 10000, 11500]


CURRENT = os.environ.get("LOOP_FIR_CURRENT") == "1"   # --current: the correction on top of the current S baseline's loops (model_s_results 25.14)
if CURRENT or KC != 4:
    OUT = OUT + ("_current" if CURRENT else "") + (f"_kc{KC}" if KC != 4 else "")


def _fit_setting(a):
    key, lab, b230, law, caps, ends, v0 = a
    q1.EXT = 1700
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    v, _, _ = q2.final_injection(compact=f"diffuser_tied_fix-0.707107_ap230new_{q1.ms.KTAG}" if CURRENT else COMPACT)
    G0 = q2.branch_spectra(lab, law, [0.0], b230, q2.FINAL_DIFF)
    N = G0["N"]
    w = 2 * np.pi * np.fft.rfftfreq(N)
    Z = [np.exp(-1j * w * Pi) for Pi in q1.P]
    J = q1.Joint([dict(G0, files=[dict(y=np.zeros((N, 2)), w=np.ones(N), amp=0)])])     # full length (not 4096)
    cl, f, xs = J.split(v)
    F, X = J._spectra(f, xs, N)
    Cl = J.cluster(cl)
    PB = G0["PB"]
    if CURRENT:                                     # §25.11: fold constants held out by RT step, LF 0.4 dc, floor and HF cells
        import dense_mid as dm
        G_law = dm.loop_G(key, N)
    else:
        G_law = q2.G_capped(lab, law, w, beta=q2.FINAL_BETA)

    def reverb(Gs):
        R = np.array([PB / (1 - Gi * Zi) for Gi, Zi in zip(Gs, Z)])
        return np.fft.irfft(F * PB[None] + np.einsum("cik,ik->ck", X, R), N, axis=1).T

    r_law = reverb(G_law)
    data = []
    for c in caps:
        yb = q2._band(c["file"])
        n = min(len(yb), N)
        mdl = (c["amp"] / 90) * Cl[:n] + r_law[:n]
        t_end = max(c["fade2"], 0.15)
        _, (e, sh, gain) = q1._clock_job((yb[:n], mdl, t_end, c["eps"], c["amp"] == 0))
        y = q1.apply_clock(yb[:n], e, sh) / gain
        i0, i1 = PRE + int(0.05 * FS), min(PRE + int(t_end * FS), n)
        t = (np.arange(i0, i1) - PRE) / FS
        wt = 10 ** (1.5 * np.minimum(t, q1.T_CAP) / lab["rt"])[:, None]
        data.append(dict(c, y=y, i0=i0, i1=i1, wt=wt / np.sqrt(np.sum((y[i0:i1] * wt) ** 2))))
    train = [d for d in data if d["amp"] in (30, 90)]

    def err_r(r):
        return np.concatenate([(((d["amp"] / 90) * Cl[d["i0"]:d["i1"]] + r[d["i0"]:d["i1"]]) - d["y"][d["i0"]:d["i1"]]) * d["wt"]
                               for d in train]).ravel()

    # five: q1's 5-number form, shared, length-scaled (as q2's loops stage, but over the whole window)
    idx = [0, 4] + ([1, 2, 3] if lab["lf"] != 1.0 else [])
    lo = np.array([0.05, 0.02, 0.2, -0.5, 0.05, -300]); hi = np.array([20, 20, 5, 0.99999, 1.0, 300])
    xsc = np.array([0.001 * lab["rt"], 0.01, 0.01, 0.001, 0.01, 1.0])[idx]
    vv = np.array(v0, float); vv[5] = 0.0; vv[4] = min(vv[4], 1.0)

    def five(u):
        p = vv.copy(); p[idx] = u
        return q1.G_loops(p, w)
    if CURRENT:                                     # no 5-number refit: only the correction on the baseline
        f5 = type("r", (), dict(x=vv[idx], cost=0.5 * float(np.sum(err_r(r_law) ** 2)), nfev=0))
    else:
        f5 = least_squares(lambda u: err_r(reverb(five(u))), vv[idx], bounds=(lo[idx], hi[idx]), x_scale=xsc, diff_step=1e-6, max_nfev=100)
    p5 = vv.copy(); p5[idx] = f5.x

    # fir: G_i C_i with C_i a (2 KC + 1)-tap FIR, exact Jacobian
    lags = np.arange(-KC, KC + 1)
    E = np.exp(-1j * np.outer(lags, w))                                          # (taps, bins)
    nt = len(lags)

    def Cs(u):
        return [u[i * nt:(i + 1) * nt] @ E for i in range(4)]

    def res_fir(u):
        return err_r(reverb([Gi * Ci for Gi, Ci in zip(G_law, Cs(u))]))

    def jac_fir(u):
        C = Cs(u)
        cols = []
        for i in range(4):
            D = 1 - G_law[i] * C[i] * Z[i]
            base = X[:, i] * (PB * G_law[i] * Z[i] / D ** 2)[None]                 # (2, bins)
            for k in range(nt):
                dr = np.fft.irfft(base * E[k][None], N, axis=1).T                # (N, 2)
                cols.append(np.concatenate([(dr[d["i0"]:d["i1"]] * d["wt"]).ravel() for d in train]))
        return np.stack(cols, 1)
    u0 = np.zeros(4 * nt); u0[[i * nt + KC for i in range(4)]] = 1.0
    ff = least_squares(res_fir, u0, jac=jac_fir, x_scale=0.01, max_nfev=60)

    out = dict(key=key, **lab, b230=b230, five=p5[:5].tolist(), fir=ff.x.reshape(4, nt).tolist(), lags=lags.tolist(),
               cost=dict(law=float(np.sum(err_r(r_law) ** 2)), five=float(2 * f5.cost), fir=float(2 * ff.cost)),
               nfev=dict(five=int(f5.nfev), fir=int(ff.nfev)), rows=[])
    for name, Gs in ((("law", G_law),) + ((("five", five(f5.x)),) if not CURRENT else ()) + (("fir", [Gi * Ci for Gi, Ci in zip(G_law, Cs(ff.x))]),)):
        G = q2.make_G(lab, law, Gs, b230, q2.FINAL_DIFF)
        for c in caps:
            out["rows"].append(dict(key=key, **lab, amp=c["amp"], loops=name, held_out=c["amp"] not in (30, 90),
                                    **q2.score_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)))
    return out


def fit_stage(workers, only=None):
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    fits = json.load(open(os.path.join(q2.Q1, "loops.json")))
    law = q1.fit_laws(fits)
    os.makedirs(OUT, exist_ok=True)
    keys = [k for k in info if only is None or k in only]
    todo = [k for k in keys if not os.path.exists(os.path.join(OUT, k + ".json"))]
    print(f"loop_fir: {len(keys)} settings, {len(keys) - len(todo)} done, {len(todo)} to fit", flush=True)
    jobs = []
    for k in todo:
        g = info[k]
        caps = [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in q2.files_of(g)]
        v0 = np.median([r["v"] for r in fits if r["file"] in {c["file"] for c in caps}], 0)
        jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
                     caps, ends, v0))
    with ProcessPoolExecutor(workers, max_tasks_per_child=2) as ex:
        for fu in as_completed([ex.submit(_fit_setting, j) for j in jobs]):
            r = fu.result()
            json.dump(r, open(os.path.join(OUT, r["key"] + ".json"), "w"), indent=0)
            print(f"{r['key']}: cost law {r['cost']['law']:.3g} -> five {r['cost']['five']:.3g} -> fir {r['cost']['fir']:.3g}", flush=True)


def corr_shape(c, lags, f):
    """Magnitude (dB) and phase delay (samples) per pass of a correction FIR at frequencies f (Hz)."""
    w = 2 * np.pi * np.asarray(f) / FS
    H = np.exp(-1j * np.outer(w, lags)) @ np.asarray(c)
    return 20 * np.log10(np.abs(H)), -np.unwrap(np.angle(H)) / w


def report_stage():
    res = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(OUT, "*.json")))]
    rep = {(r["key"], r["amp"]): r for r in json.load(open(next(p_ for p_ in (os.path.join(q2.OUT, q1.ms.tagged("nulls_free") + ".json"), os.path.join(q2.OUT, "nulls_free_chain2.json")) if os.path.exists(p_))))}
    ws = ["50-150", "150-300", "300-600", "600-1000"]
    rows = [r for s in res for r in s["rows"] if r["held_out"]]

    def ex(rr, w):
        e = [r["null_comp"][w][1] - rep[(r["key"], r["amp"])]["repeat_comp"][w][1] for r in rr
             if r["null_comp"].get(w) and rep.get((r["key"], r["amp"])) and rep[(r["key"], r["amp"])]["repeat_comp"].get(w)]
        return np.median(e) if e else np.nan

    def nul(rr, w):
        x = [r["null_comp"][w][1] for r in rr if r["null_comp"].get(w)]
        return np.median(x) if x else np.nan
    groups = [("all", lambda r: True), ("flat", lambda r: r["lf"] == 1.0 and r["hf"] == 1.0), ("LF 0.1", lambda r: r["lf"] == 0.1),
              ("HF 0.2", lambda r: r["hf"] == 0.2), ("RT <= 1.0", lambda r: r["rt"] <= 1.0), ("RT >= 1.4", lambda r: r["rt"] >= 1.4)]
    print(f"## Loop forms, held out by capture (AMPL 0 + 60; fitted on 30 + 90), {len(res)} settings; R, fade-compensated null [excess over repeat]")
    print("| settings | loops | " + " | ".join(ws) + " |\n|---|---|" + "---|" * len(ws))
    for g, sel in groups:
        for lp in ("law", "five", "fir"):
            rr = [r for r in rows if r["loops"] == lp and sel(r)]
            if rr:
                print(f"| {g} | {lp} | " + " | ".join(f"{nul(rr, w):.1f} [{ex(rr, w):+.1f}]" for w in ws) + " |")
    print(f"\n## The fir corrections per pass, median over the settings of a group and the four loops (magnitude dB / phase delay samples)")
    print("| settings | n | " + " | ".join(f"{f} Hz" for f in FREQS) + " |\n|---|---|" + "---|" * len(FREQS))
    for g, sel in groups + [("RT 0.6", lambda r: r["rt"] == 0.6), ("RT 1.0", lambda r: r["rt"] == 1.0), ("RT 2.5", lambda r: r["rt"] == 2.5)]:
        ss = [s for s in res if sel(s)]
        if not ss:
            continue
        sh = [corr_shape(c, s["lags"], FREQS) for s in ss for c in s["fir"]]
        m = np.median([x[0] for x in sh], 0); d = np.median([x[1] for x in sh], 0)
        print(f"| {g} | {len(ss)} | " + " | ".join(f"{a:+.3f} / {b:+.3f}" for a, b in zip(m, d)) + " |")
    print("\n## Per loop (median over all settings): magnitude dB / phase delay samples")
    for i, Pi in enumerate(q1.P):
        sh = [corr_shape(s["fir"][i], s["lags"], FREQS) for s in res]
        print(f"- loop {Pi}: " + "; ".join(f"{f}: {a:+.3f} / {b:+.3f}" for f, a, b in zip(FREQS, np.median([x[0] for x in sh], 0), np.median([x[1] for x in sh], 0))))


def law_coefs(lab, law):
    """Per loop the capped law's coefficients as G_capped builds them: (g, r0, a, p); r0 = 1 without a shelf, p = 0
    without the HF lowpass."""
    L = law["lf"][lab["lf"]]
    out = []
    for Pi in q1.P:
        kk = -3 * Pi / FS
        g = 10 ** (kk / (L["hi"] * lab["rt"]))
        r0 = 1.0
        if L["dc"] != L["hi"]:
            r0 = 10 ** (kk * (1 / (L["dc"] * lab["rt"]) - 1 / (L["hi"] * lab["rt"])))
            if g * r0 < q2.DC_MIN:
                g = g * 10 ** (-q2.FINAL_BETA * 20 * np.log10(q2.DC_MIN / (g * r0)) / 20)
            r0 = max(r0, q2.DC_MIN / g)
        p = 0.0
        if lab["hf"] < 1.0:
            p = min(float(q1.onepole_for(np.array(10 ** (kk / lab["rt"] * L["k"] * (1 / lab["hf"] - 1))), law["f_ref"])), q2.P_MAX)
        out.append((g, r0, law["a"], p))
    return out


def log_G(theta, w):
    """ln G(w) of the law form g . shelf(r0, 1, a) . (1 - p)/(1 - p z^-1) . z^-d, theta = (g, r0, a, p, d)."""
    import model_s as m
    g, r0, a, p, d = theta
    z1 = np.exp(-1j * w)
    return np.log(g) + np.log(m.shelf(w, r0, 1.0, a)) + np.log((1 - p) / (1 - p * z1)) - 1j * w * d


def coefs_stage():
    """Each loop's fir correction read as changes of the law's own coefficients: C - 1 ~ sum_j dtheta_j dlnG/dtheta_j
    (complex, 100 Hz-10 kHz), theta = (g, r0, a, p, d); only the parameters the setting uses (r0, a where LF != 1.0;
    p where HF < 1.0). Reports how much of each correction that explains and the coefficients the unit used."""
    law = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    f = np.linspace(100, 10000, 400)
    w = 2 * np.pi * f / FS
    out = []
    for path in sorted(glob.glob(os.path.join(OUT, "*.json"))):
        s = json.load(open(path))
        lab = dict(lf=s["lf"], rt=s["rt"], hf=s["hf"])
        for i, ((g, r0, a, p), c) in enumerate(zip(law_coefs(lab, law), s["fir"])):
            C = np.exp(-1j * np.outer(w, s["lags"])) @ np.asarray(c)
            y = np.log(C)
            th = np.array([g, r0, a, p, 0.0])
            use = [0, 4] + ([1, 2] if lab["lf"] != 1.0 else []) + ([3] if lab["hf"] < 1.0 else [])
            cols = []
            for j in use:
                h = 1e-6 * max(abs(th[j]), 1e-3)
                tp = th.copy(); tp[j] += h
                cols.append((log_G(tp, w) - log_G(th, w)) / h)
            A = np.stack(cols, 1)
            Ar = np.vstack([A.real, A.imag]); yr = np.r_[y.real, y.imag]
            d = np.linalg.lstsq(Ar, yr, rcond=None)[0]
            expl = 1 - np.sum((yr - Ar @ d) ** 2) / np.sum(yr ** 2)
            dt = np.zeros(5); dt[use] = d
            out.append(dict(key=s["key"], **lab, loop=int(q1.P[i]), law=th[:4].tolist(), delta=dt.tolist(), explained=float(expl),
                            corr_db_8k=float(20 * np.log10(abs(np.exp(-1j * 2 * np.pi * 8000 / FS * np.asarray(s["lags"])) @ np.asarray(c))))))
    json.dump(out, open(OUT + "_coefs.json", "w"), indent=0)
    names = ["g", "r0", "a", "p", "d"]
    print(f"\n## fir corrections as law-coefficient changes ({len(out)} loops x settings): share of the correction explained, median "
          f"{np.median([o['explained'] for o in out]):.2f} (IQR {np.percentile([o['explained'] for o in out], 25):.2f}..{np.percentile([o['explained'] for o in out], 75):.2f})")
    for j, nm in enumerate(names):
        rel = [o["delta"][j] / (o["law"][j] if j < 4 and o["law"][j] else 1.0) for o in out if o["delta"][j] != 0.0]
        if rel:
            print(f"- d{nm}{'' if j == 4 else ' / ' + nm}: median {np.median(rel):+.2e}, IQR {np.percentile(rel, 25):+.2e}..{np.percentile(rel, 75):+.2e}, "
                  f"abs median {np.median(np.abs(rel)):.2e} (n {len(rel)})")
    for c, vals in (("rt", None), ("hf", None), ("lf", None)):
        for v_ in sorted({o[c] for o in out}):
            oo = [o for o in out if o[c] == v_]
            print(f"  {c} {v_}: " + "; ".join(f"d{nm} {np.median([o['delta'][j] / (o['law'][j] if j < 4 and o['law'][j] else 1.0) for o in oo if o['delta'][j] != 0.0]):+.1e}"
                                          for j, nm in enumerate(names) if any(o["delta"][j] != 0.0 for o in oo)))


FLOOR_GRID = dict(dc_min=[0.15, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28],
                  beta=[0.0, 0.019, 0.04, 0.08], gamma=[0.0, -0.01, -0.02, -0.03, -0.04])
TRAIN_RT = (0.3, 0.8, 1.4, 2.2, 2.5)
SHELF_GRID = dict(da=[-0.0006, -0.0003, -0.0001, 0.0, 0.0001, 0.0003, 0.0006, 0.001, 0.002], hi=[-0.02, -0.01, 0.0, 0.01])
SHELF_GAMMA = -0.02
HFCAP_GRID = dict(pmax=[0.84, 0.88, 0.92, 0.999], dp=[-0.02, -0.01, -0.005, 0.0, 0.005, 0.01, 0.02, 0.03])


def G_shelf(lab, law, w, da, hi):
    """The capped law with gamma = SHELF_GAMMA, the shelf pole a + da (shared by every LF step) and the LF step's
    Nyquist-end RT ratio hi x (1 + hi_scale) (the top end moves, the DC end stays: g r0 does not depend on hi)."""
    lw = json.loads(json.dumps(law)); lw["lf"] = {float(k): v for k, v in lw["lf"].items()}
    lw["a"] = law["a"] + da
    if lab["lf"] != 1.0:
        lw["lf"][lab["lf"]]["hi"] = law["lf"][lab["lf"]]["hi"] * (1 + hi)
    return q2.G_capped(lab, lw, w, beta=q2.FINAL_BETA, gamma=SHELF_GAMMA)


def G_floor(lab, law, w, dc_min, beta, gamma):
    """q2.G_capped with the floor's constants free and gamma (a constant extra broadband loss where the floor acts)."""
    return q2.G_capped(lab, law, w, dc_min=dc_min, beta=beta, gamma=gamma)


def _floor_setting(a, mode="floor"):
    """One setting: the data prepared once (clock, shift, gain against the current law, as _fit_setting), then the
    weighted residual over the whole window to the fade for every (dc_min, beta, gamma) in FLOOR_GRID."""
    key, lab, b230, law, caps = a
    q1.EXT = 1700
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    v, _, _ = q2.final_injection(compact=COMPACT)
    G0 = q2.branch_spectra(lab, law, [0.0], b230, q2.FINAL_DIFF)
    N = G0["N"]
    w = 2 * np.pi * np.fft.rfftfreq(N)
    Z = [np.exp(-1j * w * Pi) for Pi in q1.P]
    J = q1.Joint([dict(G0, files=[dict(y=np.zeros((N, 2)), w=np.ones(N), amp=0)])])
    cl, f, xs = J.split(v)
    F, X = J._spectra(f, xs, N)
    Cl = J.cluster(cl)
    PB = G0["PB"]

    def reverb(Gs):
        R = np.array([PB / (1 - Gi * Zi) for Gi, Zi in zip(Gs, Z)])
        return np.fft.irfft(F * PB[None] + np.einsum("cik,ik->ck", X, R), N, axis=1).T
    r_law = reverb(q2.G_capped(lab, law, w, beta=q2.FINAL_BETA))
    data = []
    for c in caps:
        yb = q2._band(c["file"])
        n = min(len(yb), N)
        t_end = max(c["fade2"], 0.15)
        _, (e, sh, gain) = q1._clock_job((yb[:n], (c["amp"] / 90) * Cl[:n] + r_law[:n], t_end, c["eps"], c["amp"] == 0))
        y = q1.apply_clock(yb[:n], e, sh) / gain
        i0, i1 = PRE + int(0.05 * FS), min(PRE + int(t_end * FS), n)
        t = (np.arange(i0, i1) - PRE) / FS
        wt = 10 ** (1.5 * np.minimum(t, q1.T_CAP) / lab["rt"])[:, None]
        data.append(dict(amp=c["amp"], y=y, i0=i0, i1=i1, wt=wt / np.sqrt(np.sum((y[i0:i1] * wt) ** 2))))
    def cost_of(Gs):
        r = reverb(Gs)
        return sum(float(np.sum(((((d["amp"] / 90) * Cl[d["i0"]:d["i1"]] + r[d["i0"]:d["i1"]]) - d["y"][d["i0"]:d["i1"]]) * d["wt"]) ** 2))
                   for d in data)
    grid = []
    if mode == "probe":                              # parts of the setting's own fir correction (model_s_results 22.4)
        fir = json.load(open(os.path.join(OUT, key + ".json")))
        lags = np.array(fir["lags"]); E = np.exp(-1j * np.outer(lags, w))
        C = [np.asarray(c) @ E for c in fir["fir"]]
        G_law = q2.G_capped(lab, law, w, beta=q2.FINAL_BETA)
        Cm = np.exp(np.median(np.log(np.abs(C)), 0)) * np.exp(1j * np.median(np.unwrap(np.angle(C), axis=1), 0))
        variants = {"law": [1.0] * 4, "own": C, "common": [Cm] * 4, "magnitude": [np.abs(c) for c in C],
                    "phase": [np.exp(1j * np.angle(c)) for c in C], "per-loop only": [c / Cm for c in C]}
        return dict(key=key, **lab, grid=[dict(variant=k, cost=cost_of([Gi * ci for Gi, ci in zip(G_law, v)])) for k, v in variants.items()])
    if mode == "hfcap":                              # the HF pole near and above the cap (model_s_results 22.8)
        for pm in HFCAP_GRID["pmax"]:
            for dp in HFCAP_GRID["dp"]:
                Gs = q2.G_capped(lab, law, w, beta=q2.FINAL_BETA, p_map=lambda p, pm=pm, dp=dp: min(max(p + dp, 0.0), pm))
                grid.append(dict(pmax=pm, dp=dp, cost=cost_of(Gs)))
        grid.append(dict(pmax=None, dp=None, cost=cost_of(q2.G_capped(lab, law, w, beta=q2.FINAL_BETA))))
        return dict(key=key, **lab, grid=grid)
    if mode == "shelf":
        for da in SHELF_GRID["da"]:
            for hs in SHELF_GRID["hi"]:
                grid.append(dict(da=da, hi=hs, cost=cost_of(G_shelf(lab, law, w, da, hs))))
        grid.append(dict(da=None, hi=None, cost=cost_of(q2.G_capped(lab, law, w, beta=q2.FINAL_BETA))))    # the law, no gamma
        return dict(key=key, **lab, grid=grid)
    for dm in FLOOR_GRID["dc_min"]:
        for be in FLOOR_GRID["beta"]:
            for ga in FLOOR_GRID["gamma"]:
                grid.append(dict(dc_min=dm, beta=be, gamma=ga, cost=cost_of(G_floor(lab, law, w, dm, be, ga))))
    return dict(key=key, **lab, grid=grid)


def _shelf_setting(a):
    return _floor_setting(a, "shelf")


def _hfcap_setting(a):
    return _floor_setting(a, "hfcap")


def hfcap_stage(workers, only=None):
    """HF 0.2 / 0.4 at RT <= 1.4: a grid over the HF pole's cap and an offset, p = min(p_law + dp, pmax), on full windows;
    per HF step the combination chosen leaving one RT step out, scored on that step (cost re the law)."""
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, AP230)))
    law = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    keys = [k for k, g in info.items() if g["hf"] in (0.2, 0.4) and g["rt"] <= 1.4 and (only is None or k in only)]
    jobs = [(k, dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
             [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in q2.files_of(info[k])]) for k in keys]
    with ProcessPoolExecutor(workers, max_tasks_per_child=4) as ex:
        res = list(ex.map(_hfcap_setting, jobs))
    json.dump(res, open(OUT + "_hfcap" + ("_smoke" if only else "") + ".json", "w"), indent=0)
    combos = [(g["pmax"], g["dp"]) for g in res[0]["grid"] if g["pmax"] is not None]

    def rel(s, c):
        base = next(g["cost"] for g in s["grid"] if g["pmax"] is None)
        return next(g["cost"] for g in s["grid"] if (g["pmax"], g["dp"]) == c) / base
    for hf in (0.2, 0.4):
        rts = sorted({s["rt"] for s in res if s["hf"] == hf})
        print(f"HF {hf}: best per RT step on its own (in-sample) and leave-one-RT-out (chosen on the others, dB re the law):")
        for rt in rts:
            own = [s for s in res if s["hf"] == hf and s["rt"] == rt]
            oth = [s for s in res if s["hf"] == hf and s["rt"] != rt]
            t_own = [np.mean([rel(s, c) for s in own]) for c in combos]
            t_oth = [np.mean([rel(s, c) for s in oth]) for c in combos]
            c_own, c_cv = combos[int(np.argmin(t_own))], combos[int(np.argmin(t_oth))]
            print(f"  RT {rt}: own best pmax {c_own[0]} dp {c_own[1]:+} -> {10 * np.log10(min(t_own)):+.2f} dB; "
                  f"chosen on the other RT steps pmax {c_cv[0]} dp {c_cv[1]:+} -> {10 * np.log10(np.mean([rel(s, c_cv) for s in own])):+.2f} dB")
        allc = [np.mean([rel(s, c) for s in res if s["hf"] == hf]) for c in combos]
        print(f"  all RT steps: best pmax {combos[int(np.argmin(allc))][0]} dp {combos[int(np.argmin(allc))][1]:+} -> {10 * np.log10(min(allc)):+.2f} dB")


def _probe_setting(a):
    return _floor_setting(a, "probe")


TILT_F = np.r_[np.arange(0, 1000, 50), np.arange(1000, 13000, 250), 13186.8]


def tilt_stage():
    """Per LF step (not 1.0) the common magnitude correction of the training-RT settings (HF 0.6 and 1.0), in dB per
    second, relative to its value at <= 250 Hz (the broadband part is gamma's): median over settings of the four loops'
    median -> <kernel>_tilt.json (frequencies in Hz, dB/s), for q2_injection final --tilt."""
    res = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(OUT, "*.json")))]
    w = 2 * np.pi * TILT_F / FS
    out = dict(f_hz=TILT_F.tolist(), train_rt=list(TRAIN_RT), note="dB/s re the value at <= 250 Hz; model_s_results 22.4")
    for lf in (0.1, 0.4, 2.0):
        shapes = []
        for s_ in res:
            if s_["lf"] != lf or s_["hf"] not in (0.6, 1.0) or s_["rt"] not in TRAIN_RT:
                continue
            E = np.exp(-1j * np.outer(w, s_["lags"]))
            dbps = np.median([20 * np.log10(np.abs(E @ np.asarray(c))) * FS / Pi for c, Pi in zip(s_["fir"], q1.P)], 0)
            shapes.append(dbps - np.median(dbps[TILT_F <= 250]))
        med = np.median(shapes, 0)
        med[TILT_F > 6000] = med[np.argmin(abs(TILT_F - 6000))]        # the 8-12 kHz pattern is common to every LF step
        out[str(lf)] = med.tolist()
        print(f"LF {lf} ({len(shapes)} settings): " + " ".join(f"{f:.0f}:{v:+.2f}" for f, v in zip(TILT_F[::8], np.median(shapes, 0)[::8])))
    json.dump(out, open(OUT + "_tilt.json", "w"), indent=0)


def probe_stage(workers, only):
    """Which part of a setting's own fir correction carries its gain: the correction in all four loops alike (median),
    magnitude only, phase only, the per-loop deviations only; cost over the whole window, all captures, dB re the law."""
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, AP230)))
    law = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    jobs = [(k, dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
             [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in q2.files_of(info[k])]) for k in only]
    with ProcessPoolExecutor(workers, max_tasks_per_child=4) as ex:
        res = list(ex.map(_probe_setting, jobs))
    json.dump(res, open(OUT + "_probe.json", "w"), indent=0)
    for s in res:
        base = next(g["cost"] for g in s["grid"] if g["variant"] == "law")
        print(f"{s['key']}: " + "; ".join(f"{g['variant']} {10 * np.log10(g['cost'] / base):+.2f}" for g in s["grid"] if g["variant"] != "law") + " dB re the law")


def shelf_stage(workers, only=None):
    """The shelf's pole (one shift for every LF step) and each LF step's top end at fixed DC, on top of gamma: a grid on
    full windows over the LF 0.1 / 0.4 / 2.0 settings, chosen on TRAIN_RT, reported per LF and RT (cost re the law)."""
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, AP230)))
    law = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    keys = [k for k, g in info.items() if g["lf"] != 1.0 and (only is None or k in only)]
    jobs = [(k, dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
             [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in q2.files_of(info[k])]) for k in keys]
    with ProcessPoolExecutor(workers, max_tasks_per_child=4) as ex:
        res = list(ex.map(_shelf_setting, jobs))
    json.dump(res, open(OUT + "_shelf" + ("_smoke" if only else "") + ".json", "w"), indent=0)
    combos = [(g["da"], g["hi"]) for g in res[0]["grid"] if g["da"] is not None]

    def rel(s, combo):
        base = next(g["cost"] for g in s["grid"] if g["da"] is None)
        return next(g["cost"] for g in s["grid"] if (g["da"], g["hi"]) == combo) / base
    best = {}
    for lf in (0.1, 0.4, 2.0):
        tr = [s for s in res if s["lf"] == lf and s["rt"] in TRAIN_RT]
        if not tr:
            continue
        tot = [np.mean([rel(s, c) for s in tr]) for c in combos]
        best[lf] = combos[int(np.argmin(tot))]
        print(f"LF {lf}: chosen on training RT (alone): da {best[lf][0]}, hi x (1 {best[lf][1]:+}); gamma only (da 0, hi 0): "
              f"{10 * np.log10(np.mean([rel(s, (0.0, 0.0)) for s in tr])):+.2f} dB, chosen {10 * np.log10(min(tot)):+.2f} dB re the law")
    tot = [np.mean([rel(s, c) for s in res if s["rt"] in TRAIN_RT]) for c in combos]
    da_all = combos[int(np.argmin(tot))][0]
    print(f"one da for all LF steps (training RT): {da_all}")
    for lf in (0.1, 0.4, 2.0):
        for rt in sorted({s["rt"] for s in res}):
            ss = [s for s in res if s["lf"] == lf and s["rt"] == rt]
            if ss and lf in best:
                print(f"  LF {lf} RT {rt}{' (held out)' if rt not in TRAIN_RT else ''}: gamma only {10 * np.log10(np.mean([rel(s, (0.0, 0.0)) for s in ss])):+.2f}, "
                      f"+ chosen shelf {10 * np.log10(np.mean([rel(s, best[lf]) for s in ss])):+.2f} dB re the law")
    json.dump({str(k): dict(da=v[0], hi=v[1]) for k, v in best.items()}, open(OUT + "_shelf_best.json", "w"), indent=1)


def floor_stage(workers, only=None):
    """The DC floor's constants on full windows: a grid over (dc_min, beta, gamma) on the LF 0.1 and 0.4 settings;
    chosen on the training RT steps (TRAIN_RT), reported on the others (held out by setting)."""
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, AP230)))
    law = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    keys = [k for k, g in info.items() if g["lf"] in (0.1, 0.4) and (only is None or k in only)]
    jobs = [(k, dict(lf=info[k]["lf"], rt=info[k]["rt"], hf=info[k]["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
             [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in q2.files_of(info[k])]) for k in keys]
    with ProcessPoolExecutor(workers, max_tasks_per_child=4) as ex:
        res = list(ex.map(_floor_setting, jobs))
    json.dump(res, open(OUT + "_floor" + ("_smoke" if only else "") + ".json", "w"), indent=0)
    combos = [(g["dc_min"], g["beta"], g["gamma"]) for g in res[0]["grid"]]

    def total(sel):
        out = np.zeros(len(combos))
        for s in res:
            if sel(s):
                base = next(g["cost"] for g in s["grid"] if (g["dc_min"], g["beta"], g["gamma"]) == (q2.DC_MIN, q2.FINAL_BETA, 0.0))
                out += np.array([g["cost"] for g in s["grid"]]) / base      # each setting relative to the current constants
        return out
    for name, sel in (("training RT", lambda s: s["rt"] in TRAIN_RT), ("held-out RT", lambda s: s["rt"] not in TRAIN_RT)):
        for lf in (0.1, 0.4):
            tot = total(lambda s: sel(s) and s["lf"] == lf)
            n = sum(1 for s in res if sel(s) and s["lf"] == lf)
            order = np.argsort(tot)[:5]
            print(f"{name}, LF {lf} ({n} settings): mean cost re current constants, best: " +
                  "; ".join(f"dc {combos[i][0]} beta {combos[i][1]} gamma {combos[i][2]}: {tot[i] / n:.3f}" for i in order))
    tr = total(lambda s: s["rt"] in TRAIN_RT)
    best = combos[int(np.argmin(tr))]
    print(f"\nchosen on the training RT steps (LF 0.1 + 0.4): dc_min {best[0]}, beta {best[1]}, gamma {best[2]}")
    for lf in (0.1, 0.4):
        for rt in sorted({s["rt"] for s in res}):
            ss = [s for s in res if s["lf"] == lf and s["rt"] == rt]
            rel = np.mean([next(g["cost"] for g in s["grid"] if (g["dc_min"], g["beta"], g["gamma"]) == best) /
                           next(g["cost"] for g in s["grid"] if (g["dc_min"], g["beta"], g["gamma"]) == (q2.DC_MIN, q2.FINAL_BETA, 0.0)) for s in ss])
            print(f"  LF {lf} RT {rt}{' (held out)' if rt not in TRAIN_RT else ''}: cost {10 * np.log10(rel):+.2f} dB re current")
    json.dump(dict(dc_min=best[0], beta=best[1], gamma=best[2], train_rt=TRAIN_RT), open(OUT + "_floor_best.json", "w"), indent=1)


INTERP = {0.6: (0.3, 0.8), 1.0: (0.8, 1.4), 1.8: (1.4, 2.2)}      # held-out RT step: its fitted neighbours


def _interp_setting(a):
    """One held-out setting: each loop's correction interpolated linearly in RT from the same LF / HF at the two
    neighbouring RT steps (their fir fits), applied to the capped law; every capture scored (none of this setting's
    captures entered the neighbours' fits)."""
    key, lab, b230, law, caps, ends, c_lo, c_hi, t = a
    q1.EXT = 1700
    law["lf"] = {float(k): x for k, x in law["lf"].items()}
    v, _, _ = q2.final_injection(compact=COMPACT)
    w = 2 * np.pi * np.fft.rfftfreq(q1.grid_n(lab["rt"]))
    lags = np.arange(-KC, KC + 1)
    E = np.exp(-1j * np.outer(lags, w))
    G_law = q2.G_capped(lab, law, w, beta=q2.FINAL_BETA)
    C = [((1 - t) * np.asarray(lo) + t * np.asarray(hi)) @ E for lo, hi in zip(c_lo, c_hi)]
    G = q2.make_G(lab, law, [Gi * Ci for Gi, Ci in zip(G_law, C)], b230, q2.FINAL_DIFF)
    return [dict(key=key, **lab, amp=c["amp"], loops="interp", held_out=True,
                 **q2.score_capture(G, v, c["amp"], c["file"], c["fade2"], c["eps"], ends[c["file"]], 1700)) for c in caps]


def interp_stage(workers, only=None):
    info, fade, _, _, eps0 = q2.q1_tables()
    ap = json.load(open(os.path.join(q2.OUT, AP230)))
    ends = json.load(open(os.path.join(q2.OUT, "ramp.json")))
    law = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
    fit = {json.load(open(p))["key"]: json.load(open(p)) for p in glob.glob(os.path.join(OUT, "*.json"))}
    jobs = []
    for k, g in info.items():
        if g["rt"] not in INTERP or (only and k not in only):
            continue
        r0, r1 = INTERP[g["rt"]]
        k0, k1 = f"{g['lf']}_{r0}_{g['hf']}", f"{g['lf']}_{r1}_{g['hf']}"
        caps = [dict(amp=a, file=nm, eps=eps0[nm], fade2=fade[nm]) for a, nm in q2.files_of(g)]
        jobs.append((k, dict(lf=g["lf"], rt=g["rt"], hf=g["hf"]), ap["steps"][ap["settings"][k]["step"]], json.loads(json.dumps(law)),
                     caps, ends, fit[k0]["fir"], fit[k1]["fir"], (g["rt"] - r0) / (r1 - r0)))
    with ProcessPoolExecutor(workers, max_tasks_per_child=4) as ex:
        rows = [r for rr in ex.map(_interp_setting, jobs) for r in rr]
    json.dump(rows, open(OUT + "_interp" + ("_smoke" if only else "") + ".json", "w"), indent=0)
    rep = {(r["key"], r["amp"]): r for r in json.load(open(next(p_ for p_ in (os.path.join(q2.OUT, q1.ms.tagged("nulls_free") + ".json"), os.path.join(q2.OUT, "nulls_free_chain2.json")) if os.path.exists(p_))))}
    own = {(r["key"], r["amp"], r["loops"]): r for k in {r["key"] for r in rows} for r in fit[k]["rows"]}
    ws = ["50-150", "150-300", "300-600", "600-1000"]

    def cell(rr, w):
        x = [r["null_comp"][w][1] for r in rr if r["null_comp"].get(w)]
        e = [r["null_comp"][w][1] - rep[(r["key"], r["amp"])]["repeat_comp"][w][1] for r in rr
             if r["null_comp"].get(w) and rep.get((r["key"], r["amp"])) and rep[(r["key"], r["amp"])]["repeat_comp"].get(w)]
        return f"{np.median(x):.1f} [{np.median(e):+.1f}]" if x else "--"
    print(f"\n## Held out by setting (RT {sorted(INTERP)}): corrections interpolated in RT from the neighbouring steps; R, compensated null [excess]")
    print("| settings | loops | " + " | ".join(ws) + " |\n|---|---|" + "---|" * len(ws))
    for g, sel in (("all", lambda r: True), ("LF 0.1", lambda r: r["lf"] == 0.1), ("LF 0.4", lambda r: r["lf"] == 0.4),
                   ("LF 1.0", lambda r: r["lf"] == 1.0), ("LF 2.0", lambda r: r["lf"] == 2.0), ("HF 0.2", lambda r: r["hf"] == 0.2),
                   ("HF 1.0", lambda r: r["hf"] == 1.0)):
        for lp in ("law", "interp", "fir"):
            rr = [r for r in rows if sel(r)] if lp == "interp" else [own[(r["key"], r["amp"], lp)] for r in rows if sel(r)]
            print(f"| {g} | {lp} | " + " | ".join(cell(rr, w) for w in ws) + " |")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="fit", choices=["fit", "report", "coefs", "interp", "floor", "shelf", "probe", "tilt", "hfcap"])
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--current", action="store_true", help="fit: the correction on top of the current S baseline (model_s_results 25.14)")
    ap.add_argument("--kc", type=int, default=4, help="fit: the correction's half-length in lags")
    a = ap.parse_args()
    if (a.current, a.kc) != (CURRENT, KC):          # re-run with the environment set, so the pool's workers see it
        os.environ.update(LOOP_FIR_CURRENT="1" if a.current else "0", LOOP_FIR_KC=str(a.kc))
        os.execv(sys.executable, [sys.executable] + sys.argv)
    if a.stage == "fit":
        fit_stage(a.workers, a.only.split(",") if a.only else None)
    if a.stage == "coefs":
        coefs_stage()
    elif a.stage == "tilt":
        tilt_stage()
    elif a.stage == "hfcap":
        hfcap_stage(a.workers, a.only.split(",") if a.only else None)
    elif a.stage == "probe":
        probe_stage(a.workers, a.only.split(","))
    elif a.stage == "shelf":
        shelf_stage(a.workers, a.only.split(",") if a.only else None)
    elif a.stage == "floor":
        floor_stage(a.workers, a.only.split(",") if a.only else None)
    elif a.stage == "interp":
        interp_stage(a.workers, a.only.split(",") if a.only else None)
    else:
        report_stage()


if __name__ == "__main__":
    main()
