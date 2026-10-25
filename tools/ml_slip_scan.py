"""Sample slips in the M / L captures (ml_structure.md section 23; the M / L analogue of model_s_results 28's SLIPPED).

Five M captures sat 25-30 dB above their group on the first pass (R -17...-21 dB against -48...-52), and their law-stage
clock fit had landed at -9...+4 ppm against every other M capture's 12.6 +- 0.3 (the fit traded the clock against a lag step
inside its window). This scans every M / L capture against the static model (clock, shift and gains refitted on the cluster
only): the residual / signal per 2 ms bin, L and R, from the first pass's start to where the model is no longer within
-25 dB (the reverb's time variance). A jump is >= 15 dB from one bin to the next that stays (the next 3 bins all >= 12 dB
above the median of the 5 before); it is a slip if the lag step that refits the 15 ms after it is >= 0.5 file samples. For each slip: its time after tap 0 and in the file, and the lag step, from a
shift refit on the 15 ms after the slip with the clock and gains fixed, in device and in 96 kHz file samples.

    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 9 -- .venv/bin/python tools/ml_slip_scan.py [--workers 6]
        -> analysis/ml_lines/slip_scan.{json,log}
"""
import argparse, json, os, sys

for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("EMT_KERNEL", "chain3")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
OUT = os.path.join(ROOT, "analysis", "ml_lines")
BIN = 0.002


def _job(a):
    prog, key = a
    import soundfile as sf
    from scipy.optimize import minimize_scalar
    import ml_lines as mll
    import ml_carryover as cc
    K0, FS = mll.K0, mll.FS
    W = mll.windows(prog)
    g = cc.groups(prog)[key]
    rows = {r["amp"]: r for r in json.load(open(mll.group_path(prog, key)))}
    out = []
    w = int(BIN * FS)
    M = mll.Model(prog)
    caps, _, amps = cc.captures(prog, key, max(len(cc._IR(f)) for f in g.values()) + 200)
    C, R = M.parts(len(caps[amps[-1]]), key)
    nf = mll.fit_window(prog, True)[1]
    for amp in amps:                                # clock, shift and gains refitted where no slip can sit yet: the cluster
        y = caps[amp]                               # (AMPL 0, no cluster: the first 30 ms of the first pass), so a slip
        ue = mll.usable_end(y, g[amp], prog)        # in the first pass can't bend the clock (the law fit's can)
        sl = [slice(K0 + W[ch]["cluster_gap"][0], K0 + W[ch]["cluster_gap"][1]) if amp > 0 else
              slice(K0 + W[ch]["first_pass"][0], K0 + W[ch]["first_pass"][0] + int(0.03 * FS)) for ch in (0, 1)]
        mdl, u, _, _ = mll.fit_warp(y, [C, R], sl, nf, x0=(0.0, float(np.median([rows[x]["clock"][1] for x in amps]))))
        for ch in (0, 1):
            a0 = K0 + W[ch]["first_pass"][0]
            e = y[:, ch] - mdl[:, ch]
            r = []
            for i in range(a0, ue[ch] - w, w):
                s = np.sum(y[i:i + w, ch] ** 2)
                r.append(10 * np.log10(np.sum(e[i:i + w] ** 2) / s) if s > 0 else 0.0)
            r = np.array(r)
            for j in range(5, len(r) - 3):
                pre = np.median(r[j - 5:j])
                if pre > -25:                       # past where the model holds the capture: no slip test
                    break
                if r[j] - r[j - 1] >= 15 and np.all(r[j:j + 4] - pre >= 12):
                    i0 = a0 + j * w                 # the slip's bin; refine to the sample where the residual departs
                    seg = np.abs(e[i0 - w:i0 + w])
                    base = np.sqrt(np.mean(e[i0 - 6 * w:i0 - w] ** 2))
                    i_s = i0 - w + int(np.nonzero(seg > 6 * base)[0][0]) if np.any(seg > 6 * base) else i0
                    ppm = rows[amp]["clock"][1]     # the lag step: the shift that refits the 15 ms after the slip
                    t0, t1 = i_s + int(0.001 * FS), i_s + int(0.016 * FS)
                    gm = mdl[t0:t1, ch]

                    def cost(ds):
                        m = mll.warp(mdl[:, ch:ch + 1], ds, 0.0)[t0:t1, 0]
                        return float(np.sum((y[t0:t1, ch] - np.dot(y[t0:t1, ch], m) / np.dot(m, m) * m) ** 2))
                    ds = float(minimize_scalar(cost, bounds=(-3, 3), method="bounded", options={"xatol": 1e-4}).x)
                    after = 10 * np.log10(cost(ds) / np.sum(y[t0:t1, ch] ** 2))
                    raw, fs = sf.read(g[amp])
                    pk, pky = int(np.argmax(np.abs(raw[:, 1]))), int(np.argmax(np.abs(y[:, 1])))    # R's peak: one reference for both channels
                    tf = (pk + (i_s - pky) * fs / FS) / fs
                    out.append(dict(prog=prog, key=key, amp=amp, ch=ch, file=os.path.basename(g[amp]), t_ms=1000 * (i_s - K0) / FS,
                                    t_file_ms=1000 * tf, jump_db=float(r[j] - pre), before_db=float(pre), after_db=float(np.median(r[j:j + 4])),
                                    lag_dev=ds, lag_file=ds * fs / FS, after_realign_db=float(after), clock_ppm=ppm,
                                    gm=float(np.sqrt(np.mean(gm ** 2)))))
                    break
    return out


def main():
    from concurrent.futures import ProcessPoolExecutor
    import ml_carryover as cc
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", default=None, help="smoke test: prog_key,prog_key (no files written)")
    a = ap.parse_args()
    jobs = [(p, k) for p in ("M", "L") for k in sorted(cc.groups(p))]
    if a.only:
        jobs = [tuple(x.split("_", 1)) for x in a.only.split(",")]
    res = []
    with ProcessPoolExecutor(a.workers, max_tasks_per_child=4) as ex:
        for i, r in enumerate(ex.map(_job, jobs, chunksize=2)):
            res += r
            if r:
                print(jobs[i], [(x["amp"], "LR"[x["ch"]], round(x["t_ms"], 2), round(x["lag_file"], 2)) for x in r], flush=True)
    if a.only:
        print(json.dumps(res, indent=1))
        return
    json.dump(res, open(os.path.join(OUT, "slip_scan.json"), "w"), indent=1)
    slip = [x for x in res if abs(x["lag_file"]) >= 0.5]  # a lag step; a jump without one is a model-time event (L: the onset
    other = [x for x in res if abs(x["lag_file"]) < 0.5]  # of R's recirculation 1 at 286 ms, ml_structure 23)
    lines = [f"{len(jobs)} groups scanned (M {len(cc.groups('M'))}, L {len(cc.groups('L'))}); {len({(x['prog'], x['key'], x['amp']) for x in slip})} captures "
             f"with a slip (a lag step >= 0.5 file samples); {len(other)} jumps without a lag step (listed last)"]
    for x in sorted(slip, key=lambda x: (x["prog"], x["key"], x["amp"], x["ch"])) + sorted(other, key=lambda x: (x["prog"], x["key"], x["amp"], x["ch"])):
        lines.append(f"- {x['prog']} {x['key']} AMPL {x['amp']} {'LR'[x['ch']]} ({x['file']}): {x['t_ms']:.2f} ms after tap 0 = {x['t_file_ms']:.2f} ms file time; "
                     f"residual {x['before_db']:.1f} -> {x['after_db']:.1f} dB (jump {x['jump_db']:.1f}); lag step {x['lag_dev']:+.3f} device = {x['lag_file']:+.2f} file samples, "
                     f"realigned {x['after_realign_db']:.1f} dB; law clock {x['clock_ppm']:+.2f} ppm")
    open(os.path.join(OUT, "slip_scan.log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
