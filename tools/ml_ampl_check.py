"""The AMPL law: cluster-per-AMPL against the reverb level, hardware (Tomkin 96 kHz AMPL groups, chain3 device IRs; the CLA
IR for S) against the plugin's device-domain renders (mod=0) and the Python model (ml_lines.Model parts, unit cluster =
AMPL 90 by the plugin's convention). Same measure for all, no free gain. Writes analysis/ml_ampl/.

    EMT_KERNEL=chain3 .venv/bin/python tools/ml_ampl_check.py --stage hw|render|plugin|model|report [--only N]

Measures per AMPL group (band 100 Hz-10 kHz, device rate, cluster tap 0 at K0):
  c      cluster tap amplitude per AMPL unit: rms of the extreme pair difference over the cluster span / sqrt(4 taps) / dA
  rev    rms of the cluster-free reverb (y_lo - A_lo d) in the first pass [onset-3, onset+Pmin-10) and in 40-120 ms
  ratio  20 log10(90 c / rev): the AMPL-90 cluster tap over the reverb rms, dB
  alpha  each capture's cluster projected on d over [-30, onset-3) (cluster only): linearity and the AMPL-0 intercept
  g      the gain between the extreme captures on the post-cluster first pass (level consistency of the library)
  fit    (model parts, M / L / S) per capture least-squares cluster and reverb gains, gc / gr / AMPL x 90 in dB
"""
import argparse, glob, json, os, subprocess, sys
import numpy as np
import soundfile as sf
from scipy.signal import sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import ml_structure as ml  # noqa: E402
import ml_carryover as cc  # noqa: E402
from core import subsample_peak  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_ampl")
REN = os.path.join(OUT, "renders")
BIN = os.path.join(ROOT, "plugin", "build-render", "REV246Render_artefacts", "Release", "REV246Render")
K0, FS = cc.K0, cc.FS
CL_END = {"S": (288, 353), "M": (576, 704), "L": (1152, 1407)}       # last cluster tap per channel (L, R)
N_EVAL = {"S": 6000, "M": 9000, "L": 14000}


def keys(prog, only=None):
    ks = sorted(cc.groups(prog))
    return ks[::max(1, len(ks) // only)][:only] if only else ks


def measure(caps, prog):
    """caps: {amp: (n, 2) band-limited, tap 0 at K0}. Returns the per-group measures."""
    amps = sorted(caps)
    lo, hi = amps[0], amps[-1]
    dA = hi - lo
    d = (caps[hi] - caps[lo]) / dA
    o = cc.ONSET[prog]
    Pm = min(cc.LINES[prog])
    res = dict(amps=amps, c=[], rev_fp=[], rev_40=[], ratio_fp=[], ratio_40=[], alpha={}, g=[])
    r = caps[lo] - lo * d
    for ch in (0, 1):
        sl = slice(K0 - 30, K0 + CL_END[prog][ch] + 30)
        c = float(np.sqrt(np.sum(d[sl, ch] ** 2) / 4))
        fp = r[K0 + o[ch] - 3:K0 + o[ch] + Pm - 10, ch]
        w4 = r[K0 + int(0.04 * FS):K0 + int(0.12 * FS), ch]
        rf, r4 = float(np.sqrt(np.mean(fp ** 2))), float(np.sqrt(np.mean(w4 ** 2)))
        res["c"].append(c); res["rev_fp"].append(rf); res["rev_40"].append(r4)
        res["ratio_fp"].append(float(20 * np.log10(90 * c / rf)))
        res["ratio_40"].append(float(20 * np.log10(90 * c / r4)))
        cw = slice(K0 - 30, K0 + o[ch] - 3)
        res["alpha"][ch] = {int(a): float(np.sum(caps[a][cw, ch] * d[cw, ch]) / np.sum(d[cw, ch] ** 2)) for a in amps}
        gw = slice(K0 + max(CL_END[prog]) + 50, K0 + o[ch] + Pm - 10)
        res["g"].append(float(20 * np.log10(np.sum(caps[hi][gw, ch] * caps[lo][gw, ch]) / np.sum(caps[lo][gw, ch] ** 2))))
    return res


def model_fit(caps, prog, key, M):
    """Per capture, cluster and reverb gains of the model parts (fit window as ml_lines' law stage, no warp)."""
    import ml_lines as mll
    n = len(caps[max(caps)])
    C, R = M.parts(n, key)
    out = {}
    for a, y in caps.items():
        sl, nf = mll.fit_window(prog, a > 0)
        _, gs = mll.gains(y[:nf], [C[:nf], R[:nf]], sl)
        out[int(a)] = gs
    return out


# ---------------------------------------------------------------- hardware
def hw_stage(only, fit=False):
    import ml_lines as mll
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for p in "SML":
        res[p] = {}
        M = mll.Model(p) if fit else None
        for k in keys(p, only):
            caps, _, amps = cc.captures(p, k)
            caps = {a: y[:N_EVAL[p]] for a, y in caps.items()}
            res[p][k] = measure(caps, p)
            if fit:
                res[p][k]["fit"] = model_fit(caps, p, k, M)
        print(p, len(res[p]), flush=True)
    # CLA (S, AMPL 80, one capture): R channel cluster span is clean (R onset 572 > 353)
    h = np.load(os.path.join(ROOT, "analysis", "shortRT", "h_dev.npy"))
    yb = sosfiltfilt(ml.BAND, np.vstack([np.zeros((1000, 2)), h]), axis=0)
    res["CLA"] = single_capture(yb, "S", 80)
    tag = "_smoke" if only else ""
    json.dump(res, open(os.path.join(OUT, f"hw{tag}.json"), "w"))


def single_capture(y, prog, amp):
    """One capture: the R channel's cluster rms over its 4 taps (before the R onset) per AMPL, over the reverb rms."""
    i = int(np.argmax(np.abs(y[:4000, 1]) > 0.5 * np.abs(y[:4000, 1]).max()))
    t0 = subsample_peak(y[:, 1], i)[0]
    y = cc.fshift(y, K0 - t0)
    o = cc.ONSET[prog][1]
    Pm = min(cc.LINES[prog])
    c = float(np.sqrt(np.sum(y[K0 - 30:K0 + CL_END[prog][1] + 30, 1] ** 2) / 4)) / amp
    rf = float(np.sqrt(np.mean(y[K0 + o - 3:K0 + o + Pm - 10, 1] ** 2)))
    r4 = float(np.sqrt(np.mean(y[K0 + int(0.04 * FS):K0 + int(0.12 * FS), 1] ** 2)))
    return dict(amp=amp, c=c, ratio_fp=float(20 * np.log10(90 * c / rf)), ratio_40=float(20 * np.log10(90 * c / r4)))


def hw_single(only):
    """The same single-capture R measure on Tomkin S captures (like-for-like with CLA)."""
    out = []
    for k in keys("S", only):
        g = cc.groups("S")[k]
        for a, f in g.items():
            if a > 0:
                out.append(dict(key=k, **single_capture(cc._IR(f, band=True), "S", a)))
    return out


# ---------------------------------------------------------------- plugin
def rpath(p, k, a):
    return os.path.join(REN, p, f"{k}_a{a}.wav")


def render_stage(only):
    lines = []
    for p in "SML":
        for k in keys(p, only):
            lf, rt, hf = k.split("_")
            for a in sorted(cc.groups(p)[k]):
                f = rpath(p, k, a)
                if not os.path.exists(f):
                    os.makedirs(os.path.dirname(f), exist_ok=True)
                    lines.append(f"out={f} program={p} lf={lf} rt={rt} hf={hf} ampl={a} delay=0 lp=12 sixteen=0 mod=0 "
                                 f"len={N_EVAL[p] / FS + 0.1:.3f} domain=device")
    print(len(lines), "renders", flush=True)
    if lines:
        lst = os.path.join(OUT, "render_list.txt" if not only else "render_list_smoke.txt")
        open(lst, "w").write("\n".join(lines) + "\n")
        r = subprocess.run([BIN, lst], capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(r.stderr[-2000:])
        print(r.stdout.splitlines()[:2])


def load_render(f, t0=None):
    y, _ = sf.read(f, always_2d=True)
    y = sosfiltfilt(ml.BAND, np.vstack([np.zeros((K0 + 1000, 2)), y.astype(np.float64)]), axis=0)   # pad first: tap 0 may be sample 0
    if t0 is None:
        i = int(np.argmax(np.abs(y[:6000, 1]) > 0.5 * np.abs(y[:6000, 1]).max()))
        t0 = subsample_peak(y[:, 1], i)[0]
    return cc.fshift(y, K0 - t0), t0


def plugin_stage(only, fit):
    import ml_lines as mll
    res = {}
    for p in "SML":
        res[p] = {}
        M = mll.Model(p) if fit else None
        for k in keys(p, only):
            amps = sorted(cc.groups(p)[k])
            hi = amps[-1]
            yh, t0 = load_render(rpath(p, k, hi))
            caps = {hi: yh[:N_EVAL[p]]}
            for a in amps[:-1]:                     # same tap-0 time as the highest AMPL (the AMPL-0 render has no cluster)
                caps[a] = load_render(rpath(p, k, a), t0)[0][:N_EVAL[p]]
            res[p][k] = measure(caps, p)
            if fit:
                res[p][k]["fit"] = model_fit(caps, p, k, M)
        print(p, len(res[p]), flush=True)
    tag = "_smoke" if only else ""
    json.dump(res, open(os.path.join(OUT, f"plugin{tag}.json"), "w"))


def model_stage(only):
    """ml_lines' model: unit cluster = AMPL 90 (the plugin's convention): captures {0: R, 90: C + R}."""
    import ml_lines as mll
    res = {}
    for p in "SML":
        M = mll.Model(p)
        res[p] = {}
        for k in keys(p, only):
            if k not in M.b and p != "S":
                continue
            C, R = M.parts(N_EVAL[p], k)
            res[p][k] = measure({0: R, 90: C + R}, p)
        print(p, len(res[p]), flush=True)
    tag = "_smoke" if only else ""
    json.dump(res, open(os.path.join(OUT, f"model{tag}.json"), "w"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--only", type=int, default=None)
    ap.add_argument("--fit", action="store_true")
    a = ap.parse_args()
    if a.stage == "hw":
        hw_stage(a.only, a.fit)
        json.dump(hw_single(a.only), open(os.path.join(OUT, "hw_single%s.json" % ("_smoke" if a.only else "")), "w"))
    elif a.stage == "render":
        render_stage(a.only)
    elif a.stage == "plugin":
        plugin_stage(a.only, a.fit)
    elif a.stage == "model":
        model_stage(a.only)


if __name__ == "__main__":
    main()
