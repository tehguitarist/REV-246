"""M / L output kernel refitted on the first pass after the undamped all-pass step went to b = 0 (ml_structure §18).

    EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 8 -- .venv/bin/python tools/ml_kernel_fix.py
                                                     # -> analysis/ml_lines/kernel_fix.json (read by ml_lines.Model), .log

The kernel (11 taps per program, ml_carryover's firstpass.json) was fitted with the last all-pass's undamped step at its stale
per-setting b (M ~0.0014, L ~0.004) and absorbed part of that high-frequency loss. Per program: a correction to each of its
taps, least squares on the R channel's first-pass window of the undamped settings (b < 0.05), each capture's model scaled by
one gain per window as in the pulse check; the model is ml_lines.fitted with ML_KFIX=0 (the stored kernel). Held out by LF
step (fit on LF <= 1.0, score on LF > 1.0 and the reverse), then fitted on all for the file. 11 numbers per program, the
same count as before.
"""
import json, os, sys

os.environ["ML_KFIX"] = "0"
for _v in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, "analysis", "ml_lines", "kernel_fix")


def _job(a):
    prog, key = a
    import ml_lines as mll
    M = mll.Model(prog)
    if M.b.get(key, 0.0) >= 0.05:
        return key, None
    W = mll.windows(prog)[1]["first_pass"]
    a0, a1 = mll.K0 + W[0], mll.K0 + W[1]
    out = []
    for amp, (y, m, ue) in mll.fitted(prog, key).items():
        if min(ue) < a1:
            continue
        out.append((y[:a1 + 400, 1].astype(np.float64), m[:a1 + 400, 1].astype(np.float64)))
    return key, out


def design(caps, a0, a1, K):
    A, b = [], []
    for y, m in caps:
        yy, mm = y[a0:a1], m[a0:a1]
        g = np.dot(yy, mm) / np.dot(mm, mm)
        A.append(np.stack([np.roll(m, k)[a0:a1] for k in range(-K, K + 1)], 1))
        b.append(yy / g - mm)
    return np.vstack(A), np.concatenate(b)


def null(y, m, a0, a1):
    yy, mm = y[a0:a1], m[a0:a1]
    g = np.dot(yy, mm) / np.dot(mm, mm)
    return 10 * np.log10(np.sum((yy - g * mm) ** 2) / np.sum(yy ** 2))


def main():
    import ml_carryover as cc
    import ml_lines as mll
    res, lines = {}, []
    for prog in ("M", "L"):
        h = np.array(mll.Model(prog).h)
        K = (len(h) - 1) // 2
        W = mll.windows(prog)[1]["first_pass"]
        a0, a1 = mll.K0 + W[0], mll.K0 + W[1]
        with ProcessPoolExecutor(6, max_tasks_per_child=4) as ex:
            got = dict(ex.map(_job, [(prog, k) for k in sorted(cc.groups(prog))]))
        caps = [(k, y, m) for k, v in got.items() if v for y, m in v]
        hi = np.array([float(k.split("_")[0]) > 1.0 for k, _, _ in caps])
        conv = lambda m, c: m + sum(c[k + K] * np.roll(m, k) for k in range(-K, K + 1))
        for tr, te in ((~hi, hi), (hi, ~hi)):
            A, b = design([(y, m) for (k, y, m), t in zip(caps, tr) if t], a0, a1, K)
            c = np.linalg.lstsq(A, b, rcond=None)[0]
            before = np.median([null(y, m, a0, a1) for (k, y, m), t in zip(caps, te) if t])
            after = np.median([null(y, conv(m, c), a0, a1) for (k, y, m), t in zip(caps, te) if t])
            lines.append(f"{prog}: held out by LF ({tr.sum()} -> {te.sum()} captures): R first pass {before:.2f} -> {after:.2f} dB")
        A, b = design([(y, m) for _, y, m in caps], a0, a1, K)
        c = np.linalg.lstsq(A, b, rcond=None)[0]
        # the correction acts on the model's output, which already carries h: new kernel = h + c * h (convolution)
        hn = h + np.convolve(c, h)[K:K + len(h)]
        res[prog] = dict(kernel=hn.tolist(), old=h.tolist(), correction=c.tolist(), captures=len(caps))
        lines.append(f"{prog}: {len(caps)} undamped captures; kernel {np.round(h, 4).tolist()} -> {np.round(hn, 4).tolist()}")
    json.dump(res, open(OUT + ".json", "w"), indent=1)
    open(OUT + ".log", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
