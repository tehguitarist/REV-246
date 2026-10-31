"""M at short RT: is the plugin's extra low-band spread theta's depth relative to a short decay? (analysis/ml_fix/depth/)

    .venv/bin/python tools/memguard.py --limit-gb 4 -- .venv/bin/python tools/ml_fix_depth.py --stage run [--only k,..] [--workers 2] [--chunk i/n]
    .venv/bin/python tools/ml_fix_depth.py --stage report

The Python core (ml_variance.core, read-out form; = the plugin's C++ core to -139...-148 dB modulated, analysis/ml_fix/
port.json) measured through ml_variance's sweep (ml_variance.sweep_measure: 10 Hz-13 kHz, L 3.4 s, K0 grid, 100 Hz-10 kHz),
theta(t) = A sin(2 pi f t + 2 pi k / 3), k = 0, 1, 2, at several depths (rms) and rates. Per capture of the group: the static
model's cluster part (ml_lines.Model.parts) + the sweep reverb part, ml_lines.fit_warp (as the audit, plugin_check's ml
stage and ml_fix_diag's sw3), then the audit's EDC metrics (pulse_check.band_decay per octave, both channels) and
_ml_nulls; the 3 phases of each variant against each other (SPP) in the highest AMPL's frame. Variant "d0" = theta 0 (one
render: the static model through the sweep, the MC check). One JSON per group (resumable).
"""
import argparse, glob, json, os, sys, time

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("EMT_KERNEL", "chain3")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402
from scipy.signal import sosfiltfilt  # noqa: E402
import pulse_check as pc  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_fix", "depth")
FS = pc.FS
RATE = 0.83
VARIANTS = {"d13": (0.13, RATE), "d11": (0.11, RATE), "d09": (0.09, RATE), "d07": (0.07, RATE)}
if os.environ.get("DEPTH_VARIANTS"):                 # e.g. "r13x2:0.13:1.66,..." name:rms:rate
    VARIANTS = {v.split(":")[0]: (float(v.split(":")[1]), float(v.split(":")[2])) for v in os.environ["DEPTH_VARIANTS"].split(",")}
NPH = 3
BANDS = pc.OCT


def prep(x):
    return {f: sosfiltfilt(pc.SOS[("o", f)], x, axis=0) for f in BANDS}


def metrics(A, B, i0, i1):
    """The audit's EDC metrics, A against B (B the reference), per octave, R and L: ratio - 1 / dB."""
    out = {}
    for f in BANDS:
        for c, tag in ((1, ""), (0, "_L")):
            b = pc.band_decay(A[f][:, c], B[f][:, c], i0, i1)
            r = lambda x: None if x is None else x[0] / x[1] - 1
            out[f"oct_t{f}{tag}"], out[f"oct_edt{f}{tag}"], out[f"oct_late{f}{tag}"] = r(b["t"]), r(b["edt"]), r(b["late"])
            out[f"oct_dev{f}{tag}"] = b["edc_dev"]
    return out


def _job(a):
    prog, key, path = a
    import ml_carryover as cc
    import ml_lines as mll
    import ml_variance as mv
    t = time.time()
    K0 = cc.K0
    M = mll.Model(prog)
    g = cc.groups(prog)[key]
    n = max(len(cc._IR(f)) for f in g.values()) + 200
    caps, _, amps = cc.captures(prog, key, n)
    n = len(caps[amps[-1]])
    C, R0 = M.parts(n, key)
    tail = max(6.0, min(18.0, 0.8 * float(key.split("_")[1])))
    fit = lambda x: np.vstack([x, np.zeros((max(0, n - len(x)), 2))])[:n]
    R = {"d0": fit(mv.sweep_measure(M, key, lambda tt: 0 * tt, tail))}
    for v, (rms, rate) in VARIANTS.items():
        A = rms * np.sqrt(2)
        for k in range(NPH):
            R[f"{v}_{k}"] = fit(mv.sweep_measure(M, key, lambda tt, A=A, rate=rate, ph=2 * np.pi * k / NPH: A * np.sin(2 * np.pi * rate * tt + ph), tail))
    tr = time.time() - t
    i0 = K0 + int(((cc.ONSET[prog][1] + 20) / FS) * FS)
    a_, b_ = R0[K0:K0 + int(2 * FS), 1], R["d0"][K0:K0 + int(2 * FS), 1]
    check = float(10 * np.log10(np.sum((a_ - np.dot(a_, b_) / np.dot(b_, b_) * b_) ** 2) / np.sum(a_ ** 2)))   # gain-free: the sweep's scale
    ends = {a: min(min(mll.usable_end(caps[a], g[a], prog)), K0 + int((tail - 0.5) * FS)) for a in amps}
    pair_amp = max(amps, key=lambda a: (ends[a], a))  # phase pairs in the frame of the capture with the longest usable end
    rows = []
    for amp in amps:
        y = caps[amp]
        i1 = min(min(mll.usable_end(y, g[amp], prog)), K0 + int((tail - 0.5) * FS))
        sl, nf = mll.fit_window(prog, amp > 0)
        Y = prep(y[:i1])
        row = dict(amp=amp, t_end=(i1 - K0) / FS, vs_cap={}, null_cap={}, pairs={}, null_pairs={})
        mdl = {}
        for v, Rv in R.items():
            mdl[v], _, _, _ = mll.fit_warp(y, [C, Rv] if amp > 0 else [Rv], sl, nf)
            P = prep(mdl[v][:i1])
            row["vs_cap"][v] = metrics(P, Y, i0, i1)
            row["null_cap"][v] = pc._ml_nulls(prog, y, mdl[v], i1)
            mdl[v] = (mdl[v], P)
        if amp == pair_amp:
            for v in VARIANTS:
                for i in range(NPH):
                    for j in range(i + 1, NPH):
                        (ma, Pa), (mb, Pb) = mdl[f"{v}_{i}"], mdl[f"{v}_{j}"]
                        row["pairs"][f"{v}_{i}-{j}"] = metrics(Pa, Pb, i0, i1)
                        row["null_pairs"][f"{v}_{i}-{j}"] = pc._ml_nulls(prog, ma, mb, i1)
        rows.append(row)
        del mdl, Y
    lf, rt, hf = (float(x) for x in key.split("_"))
    json.dump(dict(key=key, prog=prog, lf=lf, rt=rt, hf=hf, variants=VARIANTS, sweep0_vs_parts_db=check, pair_amp=pair_amp, render_s=tr, rows=rows), open(path, "w"))
    r = [x for x in rows if x["amp"] == pair_amp][0]["vs_cap"]
    print(f"{prog} {key}: {time.time() - t:.0f} s (sims {tr:.0f}); d0 vs parts {check:.0f} dB; 500 Hz decay vs cap "
          + " ".join(f"{v} {np.mean([abs(r[f'{v}_{k}']['oct_t500'] or np.nan) for k in range(NPH)]) * 100:.1f}" for v in VARIANTS)
          + f" d0 {abs(r['d0']['oct_t500'] or np.nan) * 100:.1f}", flush=True)
    return key


def run(only, workers, chunk, prog="M"):
    from concurrent.futures import ProcessPoolExecutor
    import ml_carryover as cc
    d = os.path.join(OUT, os.environ.get("DEPTH_DIR", "runs"))
    os.makedirs(d, exist_ok=True)
    keys = sorted(cc.groups(prog))
    if chunk:
        i, m = (int(x) for x in chunk.split("/"))
        keys = keys[i::m]
    jobs = [(prog, k, os.path.join(d, f"{prog}_{k}.json")) for k in keys if (not only or k in only)]
    jobs = [j for j in jobs if not os.path.exists(j[2])]
    print(f"{len(jobs)} groups to run", flush=True)
    if workers == 1:
        for j in jobs:
            _job(j)
        return
    with ProcessPoolExecutor(workers, max_tasks_per_child=1) as ex:
        for _ in ex.map(_job, jobs):
            pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["run"], required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--chunk", default=None)
    ap.add_argument("--prog", default="M")
    a = ap.parse_args()
    run(a.only.split(",") if a.only else None, a.workers, a.chunk, a.prog)


if __name__ == "__main__":
    main()
