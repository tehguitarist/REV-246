"""HF 0.2's per-setting pole residuals (model_s_results 35.3): takefit's per-loop pole corrections on the §32.6 baseline (mean
of the two takes, SHELF_TAKEFIT, default _pallf) against candidates the panel provides that the pole curve ignores, each
as a leave-one-setting-out linear prediction (R^2 = 1 - SSE / sum r^2: the model predicts 0).

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds \\
        SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_pallf/cells.json \\
        [C_EXCL=2.0_0.8_0.2,0.4_0.6_0.2] .venv/bin/python tools/loop_read/pole_resid.py 0.2 <out.json>   # seconds

Candidates: the law pole before vs after the HF floor (RT x HF >= 0.12 s); the HF control's quantisation (d p_unit / d HF
through the curve, and d p_law / d HF: one HF offset per setting would predict r = slope x offset); the shelf's
response at f_ref (the LF step through the shelf-pole interaction); the law's DC request; log LF; the law pole; the
loop's length; 1 / RT; and combinations. C_EXCL leaves settings out (2.0_0.8_0.2: one take carries the 99.9 ms
defect, glitch_scan.py; 0.4_0.6_0.2: the knee setting §35.1 changes)."""
import json, os, sys, glob
sys.path.insert(0, "tools"); sys.path.insert(0, "tools/loop_read")
import numpy as np
import lowband_shelf as ls
q2, q1 = ls.q2, ls.q1
TF = os.environ.get("SHELF_TAKEFIT") or "analysis/lowband/chain3_shelf_takefit_cap0.8475_pallf"
law0 = q1.fit_laws(json.load(open(os.path.join(q2.Q1, "loops.json"))))
pc = json.load(open("analysis/lowband/chain3_pole_cells_allf.json"))
K, C = np.array(pc["knots"]), np.array(pc["all"])
curve = lambda x: float(np.interp(np.clip(x, K[0], K[-1]), K, C)) if x < K[-1] else float(C[-1])
def plaw(lab, L, hf, Pi, floor=True):
    hf_e = max(hf, q2.HFRT_MIN / lab["rt"]) if floor else hf
    return float(q1.onepole_for(np.array(10 ** (-3 * Pi / q1.FS / lab["rt"] * L["k"] * (1 / hf_e - 1))), law0["f_ref"]))
HFSEL = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
rows = []
for f in sorted(glob.glob(os.path.join(TF, f"*_{HFSEL}.json"))):
    t = json.load(open(f))
    if len(t["fit"]) < 2 or t["key"] in os.environ.get("C_EXCL", "").split(","):
        continue
    lab = dict(lf=t["lf"], rt=t["rt"], hf=t["hf"])
    tt, lfdc, kw = ls.base_kw(t["key"], ls.HF_CELLS_JSON)
    law = json.loads(json.dumps(law0)); law["lf"] = {float(k): x for k, x in law["lf"].items()}
    for k_, x in lfdc.items():
        law["lf"][float(k_)]["dc"] = float(x)
    L = law["lf"][lab["lf"]]
    wref = 2 * np.pi * law["f_ref"] / q1.FS
    w = np.array([0.0, wref, np.pi])
    G1 = q2.final_Gs(dict(lab, hf=1.0), law, w, lfrt=tt["lfrt"], gamma=tt["gamma"], hfdp=None, pramp=tt["pramp"], **kw)
    A, B = np.array(t["fit"]["A"]), np.array(t["fit"]["B"])
    for i, Pi in enumerate(q1.P):
        pl, plr = plaw(lab, L, lab["hf"], Pi), plaw(lab, L, lab["hf"], Pi, floor=False)
        h = 0.01
        dpu = (curve(plaw(lab, L, lab["hf"] + h, Pi)) - curve(plaw(lab, L, lab["hf"] - h, Pi))) / (2 * h) if pl < 0.856 else 0.0
        dpl = (plaw(lab, L, lab["hf"] + h, Pi) - plaw(lab, L, lab["hf"] - h, Pi)) / (2 * h)
        shelf_ref = 20 * np.log10(abs(G1[i][1]) / abs(G1[i][2]))      # the shelf's residual at f_ref re Nyquist (dB)
        kk = -3 * Pi / q1.FS
        x = 20 * np.log10(10 ** (kk / (L["hi"] * lab["rt"])) * 10 ** (kk * (1 / (L["dc"] * lab["rt"]) - 1 / (L["hi"] * lab["rt"]))))
        rows.append(dict(key=t["key"], **lab, loop=i, r=1e3 * (A[8 + i] + B[8 + i]) / 2, nz=1e3 * (A[8 + i] - B[8 + i]) / 2,
                         dg=(A[i] + B[i]) / 2, ddc=(A[4 + i] + B[4 + i]) / 2,
                         pl=pl, plr=plr, floorgap=plr - pl, dpu=dpu, dpl=dpl, shelf=shelf_ref, x=x, lf_=np.log(lab["lf"]),
                         Pc=(Pi - np.mean(q1.P)) / np.mean(q1.P), rt_=lab["rt"], dropped=pl >= 0.856))
r = np.array([d["r"] for d in rows]); nz = np.array([d["nz"] for d in rows])
keys = sorted({d["key"] for d in rows})
print(f"HF {HFSEL}: {len(keys)} settings, {len(rows)} loops; residual rms {np.sqrt(np.mean(r**2)):.3f}e-3, take half-difference rms {np.sqrt(np.mean(nz**2)):.3f}e-3 "
      f"(explainable share of sum r^2: {1 - np.sum(nz**2)/np.sum(r**2):.2f})")
print("drop-back loops (law pole >= 0.856):", sum(d["dropped"] for d in rows))
def corr(f):
    v = np.array([f(d) for d in rows]); return np.corrcoef(v, r)[0, 1] if np.std(v) > 0 else np.nan
cands = {
    "pre-floor law pole gap (plr - pl)": lambda d: d["floorgap"],
    "HF quantisation: dp_unit/dHF (curve)": lambda d: d["dpu"],
    "dp_law/dHF": lambda d: d["dpl"],
    "shelf at f_ref (dB)": lambda d: d["shelf"],
    "law DC request x (dB)": lambda d: d["x"],
    "log LF": lambda d: d["lf_"],
    "law pole p_law": lambda d: d["pl"],
    "loop (P - Pbar)/Pbar": lambda d: d["Pc"],
    "1/RT": lambda d: 1 / d["rt_"],
}
print("\ncorrelation with the residual, and leave-one-setting-out R^2 (1 - SSE/sum r^2), no intercept | with intercept:")
def loso(feats, icpt):
    X = np.array([[f(d) for f in feats] + ([1.0] if icpt else []) for d in rows])
    pred = np.zeros_like(r)
    for k in keys:
        te = np.array([d["key"] == k for d in rows])
        co = np.linalg.lstsq(X[~te], r[~te], rcond=None)[0]
        pred[te] = X[te] @ co
    return 1 - np.sum((r - pred) ** 2) / np.sum(r ** 2), pred
res = {}
for nm, f in cands.items():
    a, _ = loso([f], False); b, _ = loso([f], True)
    res[nm] = (corr(f), a, b)
    print(f"  {nm:40s} rho {corr(f):+.2f}   LOSO R2 {a:+.2f} | {b:+.2f}")
a, _ = loso([], True); print(f"  {'intercept only':40s}            LOSO R2 {'':5s} | {a:+.2f}")
combos = {"dpu + shelf": [cands["HF quantisation: dp_unit/dHF (curve)"], cands["shelf at f_ref (dB)"]],
          "floorgap + dpu + shelf": [cands["pre-floor law pole gap (plr - pl)"], cands["HF quantisation: dp_unit/dHF (curve)"], cands["shelf at f_ref (dB)"]],
          "shelf + x + logLF": [cands["shelf at f_ref (dB)"], cands["law DC request x (dB)"], cands["log LF"]]}
for nm, fs in combos.items():
    a, _ = loso(fs, False); b, pr = loso(fs, True)
    print(f"  {nm:40s}            LOSO R2 {a:+.2f} | {b:+.2f}")
print("\nper setting: residual (loops 1-4, e-3) | best single candidate's held-out prediction")
best = max(cands, key=lambda n: max(res[n][1], res[n][2]))
_, pr = loso([cands[best]], res[best][2] > res[best][1])
for k in keys:
    s = [i for i, d in enumerate(rows) if d["key"] == k]
    print(f"  {k:12s} " + " ".join(f"{r[i]:+.2f}" for i in s) + "  | " + " ".join(f"{pr[i]:+.2f}" for i in s) + "   (noise " + " ".join(f"{abs(nz[i]):.2f}" for i in s) + ")")
print("best single:", best)
json.dump(dict(hf=HFSEL, n=len(rows), rms=float(np.sqrt(np.mean(r**2))), noise=float(np.sqrt(np.mean(nz**2))),
               cands={k: dict(rho=v[0], loso=v[1], loso_icpt=v[2]) for k, v in res.items()}), open(sys.argv[2], "w"), indent=1)
