"""model_s_results §16.17: AP230's damping form where b > 0. Per setting, the first pass (240-1550 samples after
tap 0) of the AMPL-0 capture against the b = 0 reference setting's, x = ref . A(w) / A0(w), as q2_injection's ap230
stage, for placements of the two-tap filter F(b) = (1 - b) + b z^-1 inside the all-pass

    A = (-g + F_ff z^-230) / (1 - g F_fb z^-230),   g = -1/sqrt(2):

  both   F_ff = F_fb = F(b)  (the current form: a delay of 230 + b, linear interpolation)
  fb     F_ff = 1, F_fb = F(b)  (damping in the feedback only)
  ff     F_ff = F(b), F_fb = 1  (damping in the feed-forward only)
  free   F_ff = F(b1), F_fb = F(b2)  (both free: contains the three above)

    EMT_KERNEL=chain2 .venv/bin/python tools/ap230_form.py [--only keys]   -> analysis/ap230_form/forms_<kernel>.{json,log}
"""
import argparse, json, os, sys
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.signal import sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import q1_generalise as q1  # noqa: E402

G = -1 / np.sqrt(2)
M = 230
N = 1 << 14
W = 2 * np.pi * np.fft.rfftfreq(N)
Z = np.exp(-1j * W * M)
Z1 = np.exp(-1j * W)
OUT = os.path.join(ROOT, "analysis", "ap230_form")
Q1 = os.path.join(ROOT, "analysis", "q1")
I0, I1 = q1.PRE + 240, q1.PRE + 1550


def F(b):
    return (1 - b) + b * Z1


def resp(bff, bfb):
    return (-G + F(bff) * Z) / (1 - G * F(bfb) * Z)


FORMS = {"both": lambda p: resp(p[0], p[0]), "fb": lambda p: resp(0.0, p[0]), "ff": lambda p: resp(p[0], 0.0),
         "free": lambda p: resp(p[0], p[1])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    info = json.load(open(os.path.join(Q1, q1.DATA)))
    steps = json.load(open(os.path.join(ROOT, "analysis", "q2", q1.ms.tagged("ap230_g0.7071") + ".json")))
    ref = sosfiltfilt(q1.LFBAND, q1.load(Q1, info[q1.AP230_REF]["files"]["0"]["file"]), axis=0)
    R = np.fft.rfft(ref[:N], axis=0)
    A0 = resp(0.0, 0.0)
    keys = a.only.split(",") if a.only else [k for k, v in steps["settings"].items() if v["step"] > 0]
    res = {}
    for key in keys:
        x = sosfiltfilt(q1.LFBAND, q1.load(Q1, info[key]["files"]["0"]["file"]), axis=0)

        def cost(H):
            mdl = np.fft.irfft(R * (H / A0)[:, None], N, axis=0)
            return np.sum((x[I0:I1] - mdl[I0:I1]) ** 2) / np.sum(x[I0:I1] ** 2)
        row = dict(step=steps["settings"][key]["step"], b_table=steps["settings"][key]["b"], b0=10 * np.log10(cost(A0)))
        for name, f in FORMS.items():
            if name == "free":
                b = row["both"]["b"]
                r = minimize(lambda p: cost(f(p)), [b, b], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-10))
                row[name] = dict(b=[float(v) for v in r.x], null=float(10 * np.log10(r.fun)))
            else:
                r = minimize_scalar(lambda v: cost(f([v])), bounds=(-0.2, 0.95), method="bounded", options={"xatol": 1e-5})
                row[name] = dict(b=float(r.x), null=float(10 * np.log10(r.fun)))
        res[key] = row
        print(f"{key} step {row['step']}: b=0 {row['b0']:.1f} | " + " | ".join(
            f"{n} {row[n]['null']:.1f} (b {np.round(row[n]['b'], 4).tolist()})" for n in FORMS), flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, q1.ms.tagged("forms") + ("_smoke" if a.only else "") + ".json"), "w"), indent=1)
    print("\nmedian first-pass null (dB) by step:")
    for st in sorted({r["step"] for r in res.values()}):
        rr = [r for r in res.values() if r["step"] == st]
        print(f"  step {st} ({len(rr)}): b=0 {np.median([r['b0'] for r in rr]):.1f}; " + "; ".join(
            f"{n} {np.median([r[n]['null'] for r in rr]):.1f}" for n in FORMS)
              + f"; free b_ff/b_fb {np.median([r['free']['b'][0] for r in rr]):.4f}/{np.median([r['free']['b'][1] for r in rr]):.4f}")


if __name__ == "__main__":
    main()
