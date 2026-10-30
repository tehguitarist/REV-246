"""Re-score the plugin (phase 5, PLUGIN_CHECK_TAG=fix renders) on the slipped N captures after undoing the library's
slip: the capture's part after a cut in the front/branch gap moved by -n 96 kHz samples (n from the bench's o, an
integer number of 96 kHz samples). Everything else is plugin_check's _n_job unchanged (shift, clock, 2 gains)."""
import json, os, sys
os.environ["PLUGIN_CHECK_TAG"] = "fix"
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
import numpy as np
import plugin_check as pc
import nl_structure as nl
import ml_lines as mll
SCR = "/Users/leighpierce/Git/246/analysis/nl/branch_slip/deslip"
STEP = 96000 / nl.FSD
_ir = nl.ir


def nhost(key):
    b = json.load(open(os.path.join(pc.N_BENCH, "_".join(str(v) for v in key) + ".json")))
    return int(round((b["o"] - 1) * STEP)), b["front_shift"]


_INV = {v: k for k, v in nl.table().items()}


def ir_deslip(key):
    y = _ir(key)
    key = key if isinstance(key, tuple) else _INV[key]
    n, fr = nhost(key)
    if n == 0 or n < -400:
        return y
    s = -n / STEP                                     # dev samples the post part must move later
    cut = nl.K0 + int(round(fr + (466 + 588 - s) / 2))
    post = np.zeros_like(y); post[cut:] = y[cut:]
    out = y.copy(); out[cut:] = 0
    return out + mll.warp(post, s, 0.0)


def job(key):
    nl.ir = ir_deslip
    pc.OUT = SCR
    return pc._n_job(key)


if __name__ == "__main__":
    os.makedirs(os.path.join(SCR, "n_fix"), exist_ok=True)
    pc.OUT = SCR
    keys = [k for k in pc.n_keys() if nhost(k)[0] != 0 and nhost(k)[0] >= -400]
    if len(sys.argv) > 1 and sys.argv[1] == "rest3":
        done = set(os.listdir(os.path.join(SCR, "n_fix")))
        keys = [k for k in keys if "_".join(str(v) for v in k) + ".json" not in done][::3]
    elif len(sys.argv) > 1 and sys.argv[1] != "all":
        keys = keys[:int(sys.argv[1])]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    print(len(keys), "captures", flush=True)
    from concurrent.futures import ProcessPoolExecutor
    if workers == 1:
        for k in keys: job(k)
    else:
        with ProcessPoolExecutor(workers, max_tasks_per_child=10) as ex:
            list(ex.map(job, keys, chunksize=2))
