"""HF filter forms at the capped cells, on the current-baseline corrections (loop_fir --current --kc 8): target H_b C_i
(the baseline's HF part times the free correction), over 177 Hz - CAPFORMS_FMAX (default 10 kHz), against forms with the main pole free or fixed
at the 0.84 cap and an extra section. model_s_results §25.15.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/capforms.py"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.optimize import least_squares
import dense_mid as dm
q2, q1, ls = dm.q2, dm.q1, dm.ls
FS = q1.FS
f = np.geomspace(177, float(os.environ.get("CAPFORMS_FMAX", 10000)), 160); w = 2 * np.pi * f / FS; z1 = np.exp(-1j * w)
op = lambda p: (1 - p) / (1 - p * z1)
CAP = q2.P_MAX
FORMS = {                                                           # (function of theta, x0, lower, upper)
    "pole":            (lambda t: t[0] * op(t[1]), [1, 0.8], [0.5, 0], [1.5, 0.999]),
    "cap+pole":        (lambda t: t[0] * op(CAP) * op(t[1]), [1, 0.1], [0.5, -0.9], [1.5, 0.95]),
    "cap+2tap":        (lambda t: t[0] * op(CAP) * ((1 - t[1]) + t[1] * z1), [1, 0.1], [0.5, -0.5], [1.5, 1.0]),
    "cap+sym3":        (lambda t: t[0] * op(CAP) * ((1 - 2 * t[1]) + t[1] * (z1 + 1 / z1)), [1, 0.05], [0.5, -0.3], [1.5, 0.5]),
    "pole+pole":       (lambda t: t[0] * op(t[1]) * op(t[2]), [1, 0.8, 0.1], [0.5, 0, -0.9], [1.5, 0.999, 0.95]),
    "pole+2tap":       (lambda t: t[0] * op(t[1]) * ((1 - t[2]) + t[2] * z1), [1, 0.8, 0.1], [0.5, 0, -0.5], [1.5, 0.999, 1.0]),
    "pole+adv":        (lambda t: t[0] * op(t[1]) * np.exp(1j * w * t[2]), [1, 0.8, 0.0], [0.5, 0, -2], [1.5, 0.999, 2]),
    "pole+interp":     (lambda t: t[0] * op(t[1]) * ((1 - t[2]) + t[2] / z1), [1, 0.8, 0.1], [0.5, 0, -1], [1.5, 0.999, 1]),
}
SUMMARY = os.environ.get("CAPFORMS_SUMMARY") == "1"
law = q2.final_injection(compact=ls.lb.lf_.COMPACT)[1]
keys = [f"{lf}_{rt}_{hf}" for hf in ("0.2", "0.4") for rt in ("0.3", "0.6", "0.8", "1.0") for lf in ("0.1", "0.4", "1.0", "2.0")] if os.environ.get("CAPFORMS_SUMMARY") == "1" else ["1.0_0.6_0.2", "2.0_0.6_0.2", "1.0_0.3_0.2", "2.0_0.3_0.2", "1.0_0.3_0.4", "0.4_0.3_0.4", "2.0_0.8_0.2", "1.0_0.8_0.2", "1.0_1.0_0.2"]
for key in keys:
    d = json.load(open(os.path.join(ROOT, "analysis/loop_fir/chain3_current_kc8", key + ".json")))
    t, lfdc, kw = ls.base_kw(key, ls.HF_CELLS_JSON)
    L = json.loads(json.dumps(law)); L["lf"] = {float(k): x for k, x in L["lf"].items()}
    for k_, x in lfdc.items():
        L["lf"][float(k_)]["dc"] = float(x)
    lab = dict(lf=d["lf"], rt=d["rt"], hf=d["hf"])
    args = dict(lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=0.0, **kw)
    Gb, G1 = q2.final_Gs(lab, L, w, **args), q2.final_Gs(dict(lab, hf=1.0), L, w, **args)
    lags = np.array(d["lags"])
    if not SUMMARY:
        print(f"{key}:", flush=True)
    summ = []
    for i in range(4):
        Hb = Gb[i] / G1[i]
        T = Hb * (np.exp(-1j * np.outer(w, lags)) @ np.array(d["fir"][i]))
        e0 = np.sqrt(np.mean(np.abs(T - Hb) ** 2))
        res = []
        for nm, (fn, x0, lo, hi) in FORMS.items():
            if SUMMARY and nm not in ("pole", "pole+adv", "pole+pole"):
                continue
            r = least_squares(lambda u: np.r_[(fn(u) - T).real, (fn(u) - T).imag], x0, bounds=(lo, hi))
            db = 20 * np.log10(np.sqrt(np.mean(np.abs(fn(r.x) - T) ** 2)) / e0)
            res.append(f"{nm} {db:+.1f} [{', '.join(f'{x:.4f}' for x in r.x)}]")
            if nm == "pole+adv":
                summ.append((r.x[1], r.x[2], r.x[0], db))
            if nm == "pole":
                summ[-1:] = summ[-1:]
                pole_db = db
        if SUMMARY:
            summ[-1] = summ[-1] + (pole_db,)
            continue
        dl = -np.angle(T / Hb)[[np.argmin(abs(f - x)) for x in (500, 2000)]] / w[[np.argmin(abs(f - x)) for x in (500, 2000)]]
        print(f"   loop {i}: correction delay {dl[0]:+.2f} / {dl[1]:+.2f} smp (500 Hz / 2 kHz); " + "; ".join(res), flush=True)
    if SUMMARY:
        u = np.array(summ)
        print(f"{key} b {json.load(open(os.path.join(q2.OUT, ls.lb.lf_.AP230)))['steps'][json.load(open(os.path.join(q2.OUT, ls.lb.lf_.AP230)))['settings'][key]['step']]:.2f} | "
              f"pole+adv: p " + " ".join(f"{x:.4f}" for x in u[:, 0]) + " | adv " + " ".join(f"{x:+.2f}" for x in u[:, 1]) +
              f" | gain {np.median(u[:, 2]):.4f} | fit {np.median(u[:, 3]):+.1f} dB (pole alone {np.median(u[:, 4]):+.1f})", flush=True)
