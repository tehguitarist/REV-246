"""Program E in the time domain with moving reads, and simulated sweep captures of it (docs/EMT250_program.md §12).

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_sim.py --stage data      # statistics of the presets
    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_sim.py --stage sim [--json params]

The static network is e_model's (§11.2, the 246 law at LF 1.44 / HF 0.72). The time-variant part (§11.3) is M reads per loop,
each at an integer output delay o(t) (as the 250's addresses: no interpolation) moving between lo and hi, up at v and down at vdn samples per second, a gain +-b (relative to the loop's static read), added to R and, with probability 'share', to L with the opposite sign
(otherwise an independent read for L). Motion laws: 'tri' (bounces between lo and hi), 'sine' (centre in [lo, hi], depth A, rate f), 'walk' (constant speed, direction
re-drawn every 'tdir' seconds). A capture is the response to an exponential sweep (sweep seconds 'tsw', unknown for the Altiverb
captures, a nuisance) deconvolved by its inverse filter, as the library captures were, so the moving reads are smeared the same way.

Statistics (stats()): the random part r = capture - static model, against the static part s:
  e_ratio  random / static energy (dB) in 300-805 / 805-2215 / 2215-4400 / 4400-9000 samples after R's tap 1 (R)
  lobes    r's decay-normalised autocorrelation at P_k + (-12..12), mean over the three loops (R)
  lr       correlation of R's and L's random parts at the static alignment (315) and +-1, +-2
  frac     static fraction per octave, 805-9000 samples: rho^2 of the capture against the static model (R)
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import lfilter, butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "tools"))
import e_model as E  # noqa: E402

FS, K0 = E.FS, E.K0
X_STATIC = None                      # e_model parameters (246 law), set by static_params()
DEF = dict(M=6, b=0.45, lo=600, hi=3600, v=60.0, vsd=0.0, vdn=None, f=0.55, A=17.0, fsd=0.0, law="tri", tdir=1.0, share=0.5, tsw=10.0, seed=0, mdel=314)
OCTS = [250, 500, 1000, 2000, 4000]


def static_params():
    """All-data static fit with the 246 law (cached)."""
    global X_STATIC
    f = os.path.join(E.OUT, "static_all_246.json")
    if not os.path.exists(f):
        E.LAWMODE = "246"; E.X0 = np.r_[E.X0[:9], 1.44, 0.72]
        x, kern, _ = E.fitstatic([r for r in E.RTS if r >= 1.0], 1.0)
        json.dump(dict(x=x.tolist()), open(f, "w"))
    E.LAWMODE = "246"
    X_STATIC = np.array(json.load(open(f))["x"])
    return X_STATIC


def loop_iir(rt, Pk, lf, hf):
    import ml_lines as mll
    g, r0, p = mll.line_coefs([Pk], rt, lf, hf)[0]
    a = mll.LAW["a"]; b0, b1 = (r0 * (1 - a) + (1 + a)) / 2, (r0 * (1 - a) - (1 + a)) / 2
    return g * (1 - p) * np.array([b0, b1]), np.convolve([1, -a], [1, -p])


def trajectories(n, prm, rng):
    """Integer output delays (n_reads, n) for each loop's M moving reads, and their gains / L signs."""
    M, lo, hi = prm["M"], prm["lo"], prm["hi"]
    t = np.arange(n) / FS; out = []
    for _ in range(3 * M):
        o0 = rng.uniform(lo, hi); d = rng.choice([-1, 1]); v = prm["v"] * rng.uniform(1 - prm["vsd"], 1 + prm["vsd"])
        if prm["law"] == "tri":         # up at v, down at vdn (asymmetric triangle; vdn None: symmetric)
            span = hi - lo; vd = v * (prm["vdn"] / prm["v"]) if prm["vdn"] else v
            Tu, Td = span / v, span / vd; ph = (rng.uniform(0, Tu + Td) + t) % (Tu + Td)
            o = np.where(ph < Tu, lo + v * ph, hi - vd * (ph - Tu))
        elif prm["law"] == "sine":       # the UAD 250's law (UAD250_test_plan §2): centre uniform in [lo, hi], depth A, rate f (+- fsd)
            f = prm["f"] * rng.uniform(1 - prm["fsd"], 1 + prm["fsd"])
            o = o0 + prm["A"] * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
        else:                            # walk: constant speed, direction re-drawn every tdir s, reflected at lo / hi
            k = int(np.ceil(t[-1] / prm["tdir"])) + 1
            dirs = rng.choice([-1, 1], k); seg = np.minimum((t / prm["tdir"]).astype(int), k - 1)
            o = o0 + np.cumsum(dirs[seg]) * v / FS
            span = hi - lo; ph = (o - lo) % (2 * span); o = lo + np.where(ph < span, ph, 2 * span - ph)
        out.append(np.round(o).astype(int))
    return np.array(out)


def render(x, rt, prm, kern, rng):
    """(n, 2) device-rate output of E for a mono input x (R's tap 1 at the input's own time + K0 - 16)."""
    q = E.unpack(X_STATIC); n = len(x)
    # direct path
    direct = x.copy(); direct[429:] += q["a"][0] * x[:-429]; direct[718:] += q["a"][1] * x[:-718]
    # all-pass chain
    u = x.copy()
    for D, g in zip(E.AP_D, q["g"]):
        bb = np.zeros(D + 1); bb[0] = -g; bb[D] = 1; aa = np.zeros(D + 1); aa[0] = 1; aa[D] = -g
        u = lfilter(bb, aa, u)
    # loops: line[t] = G{ u[t] - line[t - P] }, block-wise (block < P)
    lines = []
    for Pk, (lf_, hf_) in zip(E.P, [tuple(q["m"])] * 3):
        bG, aG = loop_iir(rt, Pk, lf_, hf_); line = np.zeros(n); zi = np.zeros(max(len(aG), len(bG)) - 1); B = 2048
        for t0 in range(0, n, B):
            t1 = min(n, t0 + B); fb = np.zeros(t1 - t0)
            s0 = t0 - Pk
            if t1 - Pk > 0:
                a0 = max(s0, 0); fb[a0 - s0:] = line[a0:t1 - Pk]
            line[t0:t1], zi = lfilter(bG, aG, u[t0:t1] - fb, zi=zi)
        lines.append(line)
    stat = np.zeros(n)
    for line, dk, ck in zip(lines, E.DK, q["c"]):
        stat[dk:] += ck * line[:-dk]
    # moving reads
    trR = trajectories(n, prm, rng); trL = trajectories(n, prm, rng)
    sg = rng.choice([-1, 1], len(trR)); shared = rng.uniform(size=len(trR)) < prm["share"]
    movR = np.zeros(n); movL = np.zeros(n); idx = np.arange(n)
    for i in range(len(trR)):
        k = i // prm["M"]; ck = abs(q["c"][k]); line = lines[k]
        jR = idx - trR[i]; movR += np.where(jR >= 0, line[np.maximum(jR, 0)], 0) * sg[i] * prm["b"] * ck
        if shared[i]:
            movL -= np.where(jR >= 0, line[np.maximum(jR, 0)], 0) * sg[i] * prm["b"] * ck
        else:
            jL = idx - trL[i]; movL += np.where(jL >= 0, line[np.maximum(jL, 0)], 0) * sg[i] * prm["b"] * ck
    R = direct + stat + movR
    L = np.zeros(n)
    md = prm["mdel"]
    L[E.LDEL:] = -direct[:-E.LDEL]; L[E.LREV:] += stat[:-E.LREV]; L[md:] += movL[:-md]; L *= q["gL"]
    out = np.stack([np.convolve(L, kern[0])[:n], np.convolve(R, kern[1])[:n]], 1)
    return out, (direct, stat)


def sweep(tsw, f1=20.0, f2=12500.0):
    t = np.arange(int(tsw * FS)) / FS; Lr = tsw / np.log(f2 / f1)
    s = np.sin(2 * np.pi * f1 * Lr * (np.exp(t / Lr) - 1))
    fade = int(0.01 * FS); s[:fade] *= np.hanning(2 * fade)[:fade]; s[-fade:] *= np.hanning(2 * fade)[fade:]
    return s


def capture(rt, prm, kern, tail=2.5):
    """Simulated sweep capture: (n, 2) IR with R's tap 1 at K0, and the static model rendered on the same grid."""
    rng = np.random.default_rng(prm["seed"])
    s = sweep(prm["tsw"]); lead = 4096; n = lead + len(s) + int(tail * FS)
    x = np.zeros(n); x[lead:lead + len(s)] = s
    y, _ = render(x, rt, prm, kern, rng)
    N = 1 << int(np.ceil(np.log2(n + len(s))))
    S = np.fft.rfft(s, N); inv = np.conj(S) / (np.abs(S) ** 2 + 1e-6 * np.max(np.abs(S)) ** 2)
    h = np.fft.irfft(np.fft.rfft(y, N, axis=0) * inv[:, None], N, axis=0)
    # the chain kernel puts tap 1's peak at index 16 of the kernel; the input starts at lead
    i0 = lead + 16 - K0
    return h[i0:i0 + K0 + int(tail * FS)]


def stats(Ys, Ss):
    """Statistics of a set of captures Ys (list of (n, 2)) against their static parts Ss (same shapes)."""
    out = {"e_ratio": [], "lobes": [], "lr": [], "frac": []}
    for y, s in zip(Ys, Ss):
        r = y - s; R = r[:, 1]
        out["e_ratio"].append([10 * np.log10(np.sum(R[K0 + a:K0 + b] ** 2) / np.sum(s[K0 + a:K0 + b, 1] ** 2))
                               for a, b in ((300, 805), (805, 2215), (2215, 4400), (4400, 9000))])
        env = np.sqrt(np.convolve(R ** 2, np.ones(801) / 801, "same")) + 1e-15; rn = R / env
        a = rn[K0 + 400:K0 + 6000]
        out["lobes"].append(np.mean([[np.dot(a, rn[K0 + 400 + Pk + d:K0 + 6000 + Pk + d]) / len(a) for d in range(-12, 13)] for Pk in E.P], 0))
        rr = r[K0 + 700:K0 + 6000, 1]
        out["lr"].append([np.dot(rr, r[K0 + 700 + E.LREV + d:K0 + 6000 + E.LREV + d, 0])
                          / np.sqrt(np.dot(rr, rr) * np.dot(r[K0 + 700 + E.LREV + d:K0 + 6000 + E.LREV + d, 0] ** 2, np.ones(5300))) for d in (-2, -1, 0, 1, 2)])
        fr = []
        for f in OCTS:
            sos = butter(4, [f / 1.414, f * 1.414], btype="band", fs=FS, output="sos")
            yb = sosfiltfilt(sos, y[:, 1])[K0 + 805:K0 + 9000]; sb = sosfiltfilt(sos, s[:, 1])[K0 + 805:K0 + 9000]
            fr.append(np.dot(yb, sb) ** 2 / (np.dot(yb, yb) * np.dot(sb, sb)))
        out["frac"].append(fr)
    return {k: np.mean(v, 0).tolist() for k, v in out.items()} | {k + "_sd": np.std(v, 0).tolist() for k, v in out.items()}


def show(st, tag):
    f = lambda v: " ".join(f"{x:+.2f}" for x in v)
    print(f"{tag}: e_ratio {f(st['e_ratio'])} | lr {f(st['lr'])} | frac {f(st['frac'])}")
    lo = st["lobes"]; print(f"   lobes -12..12: {' '.join(f'{x:+.3f}' for x in lo)}", flush=True)


def data_stats():
    static_params(); rts = [r for r in E.RTS if r >= 1.0]; kern = E.kernels(rts)
    n = K0 + int(2.5 * FS)
    Ys = [E.dev(r)[:n] for r in rts]; Ss = [E.render(X_STATIC, r, kern, n) for r in rts]
    st = stats(Ys, Ss); json.dump(st, open(os.path.join(E.OUT, "rand_stats_data.json"), "w"), indent=1)
    return st


def sim_stats(prm, rts=None):
    static_params(); rts = rts or [r for r in E.RTS if r >= 1.0]; kern = E.kernels([r for r in E.RTS if r >= 1.0])
    Ys, Ss = [], []
    for i, rt in enumerate(rts):
        p = dict(prm); p["seed"] = prm["seed"] * 100 + i
        y = capture(rt, p, kern); Ys.append(y); Ss.append(E.render(X_STATIC, rt, kern, len(y)))
    return stats(Ys, Ss)


def cost(st, ref):
    """Mismatch of simulated statistics st against the data's ref, each term in its own unit."""
    d = lambda k: np.array(st[k]) - np.array(ref[k])
    return float(np.sum((d("e_ratio")[1:] / 1.0) ** 2) + (d("e_ratio")[0] / 3.0) ** 2 + np.sum((d("lobes") / 0.01) ** 2)
                 + np.sum((d("lr") / 0.05) ** 2) + np.sum((d("frac") / 0.03) ** 2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--json", default="{}")
    ap.add_argument("--rts", default="")
    a = ap.parse_args()
    if a.stage == "data":
        show(data_stats(), "data")
    elif a.stage == "sim":
        prm = DEF | json.loads(a.json)
        show(sim_stats(prm, [float(v) for v in a.rts.split(",")] if a.rts else None), "sim " + json.dumps(json.loads(a.json)))
