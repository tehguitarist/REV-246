"""NEXT_STEPS §0 item 2: the low-band transfer between the S model and other units (model_s_results.md §26.2), and which
side carries it. Results: docs/model_s_results.md §27.

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 .venv/bin/python tools/unit_transfer.py --stage measure|report

  measure  CLA and the five Altiverb S presets against the S model exactly as pulse_check's units stage (all-RT constants,
           nothing refitted; per IR clock, shift, gain and the 14-number channel FIR on the cluster-only region), in
           variants of what the model passes through before the comparison:
             plain    the digital model (pulse_check --stage units);
             chain3   the model through the library's chain kernel chain3 (FIR x the schematic's high-pass, normalised
                      to 0 dB / 0 delay at 300-3000 Hz), the chain between the digital taps and Tomkin's files, fitted on
                      Tomkin's cluster-only windows. Nothing is fitted on CLA or the Altiverb: if the shared transfer
                      is the unit's own analogue chain (removed from the library by chain3 and so absent from the
                      model), this closes it; if it is the capture chains, it does not.
           Per variant and IR: nulls per window (L, R), per octave 63 Hz-8 kHz and window (L, R), and the transfer
           unit / model per sixth octave 40 Hz-2 kHz over 50-600 ms (complex ratio and coherence, per channel).
           chain3's (and chain2's) own response per sixth octave, normalised the same way -> analysis/unit_transfer/.
  report   tables -> analysis/unit_transfer/report.log
"""
import argparse, glob, json, os, sys

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import pulse_check as pc  # noqa: E402
from core import Kernel  # noqa: E402

q1, q2 = pc.q1, pc.q2
FS, PRE = pc.FS, pc.PRE
OUT = os.path.join(ROOT, "analysis", "unit_transfer")
KDIR = os.path.join(ROOT, "analysis", "library", "tomkin96", "kernels")
SIXTH = 40 * 2 ** (np.arange(0, 33) / 6)                    # 40 Hz - 1.6 kHz
OCT = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
OCT_SOS = {f: butter(3, [f / np.sqrt(2), min(f * np.sqrt(2), 11000)], btype="band", fs=FS, output="sos") for f in OCT}
WIN = [(50, 150), (150, 300), (300, 600), (600, 1000)]
T_WIN = (0.05, 0.6)                                         # transfer window (s after tap 0)


def kernel_resp(name, f):
    """Chain kernel response at frequencies f (Hz) per channel, normalised to unit gain and zero delay at 300-3000 Hz
    (a gain and a linear phase fitted there, per channel). (len(f), 2) complex."""
    K = Kernel.load(os.path.join(KDIR, name + ".npz"))
    N = 1 << 22
    ff = np.fft.rfftfreq(N, 1 / K.fs)
    G = K.response(N)
    if G.ndim == 1:
        G = np.stack([G, G], 1)
    out = np.empty((len(f), 2), complex)
    sel = (ff >= 300) & (ff <= 3000)
    for c in (0, 1):
        ph = np.unwrap(np.angle(G[:, c]))
        a = np.polyfit(ff[sel], ph[sel], 1)                 # linear phase (delay) and offset in the mid band
        mag = np.abs(G[:, c]) / np.exp(np.mean(np.log(np.abs(G[sel, c]))))
        resid = ph - np.polyval(a, ff)
        out[:, c] = np.interp(f, ff, mag) * np.exp(1j * np.interp(f, ff, resid))
    return out


def through(x, H):
    """x (n, 2) on the device grid through H(f) (callable on Hz -> (len, 2) complex)."""
    n = len(x)
    N = 1 << int(np.ceil(np.log2(n + 16384)))
    f = np.fft.rfftfreq(N, 1 / FS)
    return np.fft.irfft(np.fft.rfft(x, N, axis=0) * H(f), N, axis=0)[:n]


def transfer(y, m, i1):
    """unit / model per sixth octave over T_WIN (Tukey window), per channel: (ratio, coherence) lists."""
    a, b = PRE + int(T_WIN[0] * FS), min(PRE + int(T_WIN[1] * FS), i1)
    w = np.hanning(2 * int(0.02 * FS))
    tw = np.ones(b - a)
    tw[:len(w) // 2], tw[-(len(w) // 2):] = w[:len(w) // 2], w[len(w) // 2:]
    N = 1 << int(np.ceil(np.log2(4 * (b - a))))
    f = np.fft.rfftfreq(N, 1 / FS)
    Y, M = (np.fft.rfft(x[a:b] * tw[:, None], N, axis=0) for x in (y, m))
    out = []
    for c in (0, 1):
        r, co = [], []
        for fc in SIXTH:
            s = (f >= fc * 2 ** (-1 / 12)) & (f < fc * 2 ** (1 / 12))
            ym, mm, yy = np.sum(Y[s, c] * np.conj(M[s, c])), np.sum(np.abs(M[s, c]) ** 2), np.sum(np.abs(Y[s, c]) ** 2)
            r.append([float((ym / mm).real), float((ym / mm).imag)])
            co.append(float(np.abs(ym) ** 2 / (mm * yy)))
        out.append(dict(ratio=r, coh=co))
    return out


def oct_nulls(y, m, te):
    out = {}
    for f, sos in OCT_SOS.items():
        yb, mb = sosfiltfilt(sos, y, axis=0), sosfiltfilt(sos, m, axis=0)
        row = []
        for a, b in WIN:
            i0, i1 = PRE + int(a / 1000 * FS), min(PRE + int(b / 1000 * FS), PRE + int(te * FS), len(y))
            row.append(None if i1 - i0 < int(0.05 * FS) else
                       [float(10 * np.log10(np.sum((yb[i0:i1, c] - mb[i0:i1, c]) ** 2) / np.sum(yb[i0:i1, c] ** 2))) for c in (0, 1)])
        out[str(f)] = row
    return out


# ---------------------------------------------------------------- compaction: min-phase biquads at the device rate
def peak_sos(f0, g_db, q):
    """RBJ peaking EQ at the device rate."""
    A, w = 10 ** (g_db / 40), 2 * np.pi * f0 / FS
    al = np.sin(w) / (2 * q)
    b = np.array([1 + al * A, -2 * np.cos(w), 1 - al * A])
    a = np.array([1 + al / A, -2 * np.cos(w), 1 - al / A])
    return np.r_[b / a[0], a / a[0]]


def hp1_sos(fc):
    """First-order high-pass (bilinear) at the device rate."""
    k = np.tan(np.pi * fc / FS)
    return np.r_[np.array([1, -1, 0]) / (1 + k), 1, (k - 1) / (1 + k), 0]


def lowshelf_sos(f0, g_db, q):
    """RBJ low shelf at the device rate."""
    A, w = 10 ** (g_db / 40), 2 * np.pi * f0 / FS
    al, c = np.sin(w) / (2 * q), np.cos(w)
    sA = 2 * np.sqrt(A) * al
    b = np.array([A * ((A + 1) - (A - 1) * c + sA), 2 * A * ((A - 1) - (A + 1) * c), A * ((A + 1) - (A - 1) * c - sA)])
    a = np.array([(A + 1) + (A - 1) * c + sA, -2 * ((A - 1) + (A + 1) * c), (A + 1) + (A - 1) * c - sA])
    return np.r_[b / a[0], a / a[0]]


def bq_resp(x, f):
    """x = [hp fc, shelf (f0, g_db, q), (f0, g_db, q) x 2]: a first-order high-pass, a low shelf and two peaking
    sections (minimum phase, the transfer's form: a lagging low-frequency boost, a lead from a capture-type corner below
    it, ripple to 400 Hz); the response at f (Hz)."""
    z = np.exp(-2j * np.pi * np.asarray(f) / FS)
    H = np.ones(len(z), complex)
    for s_ in (hp1_sos(x[0]), lowshelf_sos(*x[1:4]), peak_sos(*x[4:7]), peak_sos(*x[7:10])):
        H *= (s_[0] + s_[1] * z + s_[2] * z ** 2) / (s_[3] + s_[4] * z + s_[5] * z ** 2)
    return H


X0 = [10.0, 70.0, 4.0, 0.7, 160.0, -0.7, 1.5, 300.0, 0.3, 2.0]   # 10 numbers
LO = [1.0, 30.0, -1.0, 0.3, 100.0, -6.0, 0.3, 230.0, -6.0, 0.5]
HI = [40.0, 200.0, 12.0, 2.0, 220.0, 6.0, 8.0, 350.0, 6.0, 8.0]    # below ~50 Hz the data can't separate shelf and corner


def plain_T(R, name):
    t = R["irs"][name]["plain"]["transfer"]
    return np.array([np.mean([complex(*t[c]["ratio"][i]) for c in (0, 1)]) for i in range(len(R["sixth"]))])


def fit_bq(name, band=(40, 1600)):
    """The cascade fitted to one IR's plain transfer (mean of L and R, complex, sixth octaves in band, normalised to its
    1.6 kHz value), held to 1 at 2-10 kHz (where the transfer is flat and the mid-band nulls are at library level).
    A delay is free on top (the comparison refits the shift). Returns (x, rms complex residual in band)."""
    from scipy.optimize import least_squares
    R = json.load(open(os.path.join(OUT, "measure.json")))
    f, T = np.array(R["sixth"]), plain_T(R, name)
    T = T / T[-1]
    sel = (f >= band[0]) & (f <= band[1])
    fa = np.array([2000, 3000, 4000, 6000, 8000, 10000.0])

    def res(v):
        r = bq_resp(v[1:], f[sel]) * np.exp(-2j * np.pi * f[sel] * v[0] / FS) - T[sel]
        ra = bq_resp(v[1:], fa) * np.exp(-2j * np.pi * fa * v[0] / FS) - 1
        return np.r_[r.real, r.imag, ra.real, ra.imag]
    s_ = least_squares(res, np.r_[0.0, X0], bounds=([-2] + LO, [2] + HI))
    return s_.x[1:].tolist(), float(np.sqrt(np.mean(s_.fun[:2 * sel.sum()] ** 2)))


def variant_H(vn):
    """plain: none; chain3: the full chain kernel; bq_<IR>: the biquads fitted on that IR (bq.json); np_<IR>: that IR's
    measured sixth-octave transfer itself (complex, interpolated in log f; 1 above 1.6 kHz)."""
    if vn == "plain":
        return None
    if vn == "chain3":
        return lambda f: kernel_resp("chain3", f)
    if vn.startswith("np_"):
        R = json.load(open(os.path.join(OUT, "measure.json")))
        lf, T = np.log(R["sixth"]), plain_T(R, vn[3:])
        T = T / T[-1]

        def H(f):
            u = np.log(np.maximum(f, 1.0))
            z = np.interp(u, lf, np.abs(T), left=abs(T[0]), right=1.0) * np.exp(1j * np.interp(u, lf, np.unwrap(np.angle(T)), left=np.angle(T[0]), right=0.0))
            return np.repeat(z[:, None], 2, 1)
        return H
    x = json.load(open(os.path.join(OUT, "bq.json")))[vn[3:]]["x"]
    return lambda f: np.repeat(bq_resp(x, f)[:, None], 2, 1)


def compact_stage():
    out = {}
    for name in ("CLA", "Altiverb RT 1.6"):
        x, rms = fit_bq(name)
        out[name] = dict(x=x, rms=rms)
        print(f"{name}: high-pass {x[0]:.1f} Hz; shelf {x[1]:.1f} Hz {x[2]:+.2f} dB Q {x[3]:.2f}; peaks "
              + "; ".join(f"{x[4 + 3 * i]:.0f} Hz {x[5 + 3 * i]:+.2f} dB Q {x[6 + 3 * i]:.2f}" for i in range(2))
              + f"; rms complex residual {rms:.4f}", flush=True)
    json.dump(out, open(os.path.join(OUT, "bq.json"), "w"), indent=1)


VARIANTS = ["plain", "chain3", "np_CLA", "bq_CLA", "bq_Altiverb RT 1.6"]


def measure_stage(variants):
    import altiverb_predict as av
    os.makedirs(OUT, exist_ok=True)
    q1.EXT = 1700
    v, law, _ = q2.final_injection(compact=pc.COMPACT)
    ap = json.load(open(os.path.join(q2.OUT, pc.lf_.AP230)))
    fold = json.load(open(os.path.join(ROOT, "analysis", "loop_cv", "fold_all.json")))
    fold["pramp"] = [float(np.median([json.load(open(p))["pramp"][0] for p in
                                      glob.glob(os.path.join(ROOT, "analysis", "loop_cv", f"cap{os.environ['EMT_PMAX']}_folds", "fold_*.json"))])), 0.62]
    irs = [("CLA", 1.0, av.cla_h())] + [(f"Altiverb RT {k}", rt, av.altiverb_h(k)[0]) for k, rt in av.PRESETS.items()]
    path = os.path.join(OUT, "measure.json")
    res = json.load(open(path)) if os.path.exists(path) else {
        "sixth": SIXTH.tolist(), "irs": {},
        "kernels": {k: [[[float(z.real), float(z.imag)] for z in kernel_resp(k, SIXTH)[:, c]] for c in (0, 1)] for k in ("chain3", "chain2")}}
    for name, rt, h in irs:
        G = pc.setting_G(dict(lf=1.0, rt=rt, hf=0.9), ap["steps"][0], law, fold, None)
        yb = sosfiltfilt(q1.LFBAND, h * (np.sign(h[PRE].sum()) or 1.0), axis=0)
        n = min(len(yb) - 600, G["N"])
        m0 = q1.Joint([dict(G, files=[dict(y=np.zeros((n, 2)), w=np.ones(n), amp=80)])]).render(v, G, n, 80)
        res["irs"].setdefault(name, {"rt": rt})
        for vn in variants:
            H = variant_H(vn)
            mdl = m0 if H is None else through(m0, H)
            eg, shg, gain = av.fit_clock(yb, mdl, starts=[98.4] if name == "CLA" else [42.8])
            y = q1.apply_clock(yb, eg, shg)[:n] / gain
            mdl2, hs = av.channel_fir(y, mdl)
            eg2, shg2, gain2 = av.fit_clock(yb, mdl2, starts=[eg])
            y2 = q1.apply_clock(yb, eg2, shg2)[:n] / gain2
            te = min(1.5, q2.comp_end(y2, np.ones_like(y2)))
            res["irs"][name][vn] = dict(eps_ppm=eg2, t_end=te, null=q1.nulls(y2, mdl2, te, q2.W_COMP[:6], narrow=True),
                                        null_oct=oct_nulls(y2, mdl2, te), transfer=transfer(y2, mdl2, PRE + int(te * FS)),
                                        fir=[np.round(x, 4).tolist() for x in hs])
            r = res["irs"][name][vn]
            print(f"{name} {vn}: null R " + " ".join(f"{w}:{x[1]:.1f}" for w, x in r["null"].items() if x) + " | oct R 50-150 "
                  + " ".join(f"{f}:{r['null_oct'][str(f)][0][1]:.1f}" for f in OCT[:5]), flush=True)
    json.dump(res, open(path, "w"))


def report_stage():
    R = json.load(open(os.path.join(OUT, "measure.json")))
    f = np.array(R["sixth"])
    pick = [i for i, x in enumerate(f) if round(x) in (40, 50, 63, 80, 101, 127, 160, 202, 254, 320, 403, 508, 640, 806, 1016)]
    db = lambda z: 20 * np.log10(abs(z))
    deg = lambda z: np.degrees(np.angle(z))
    cx = lambda p: complex(*p)
    out = ["# Low-band unit transfer (tools/unit_transfer.py)", "",
           "## Transfer unit / model, 50-600 ms, sixth octaves (dB / deg; coherence), mean of L and R; chain3's own response "
           "(normalised at 300-3000 Hz)", "  Hz      " + " ".join(f"{f[i]:>13.0f}" for i in pick)]
    k3 = [np.mean([cx(R["kernels"]["chain3"][c][i]) for c in (0, 1)]) for i in range(len(f))]
    out.append("  chain3  " + " ".join(f"{db(k3[i]):+6.2f}/{deg(k3[i]):+5.1f}" for i in pick))
    for name, d in R["irs"].items():
        for vn in [v for v in VARIANTS if v in d]:
            t = d[vn]["transfer"]
            z = [np.mean([cx(t[c]["ratio"][i]) for c in (0, 1)]) for i in range(len(f))]
            out.append(f"  {name[:13]:13s} {vn:6s}".ljust(10) + "\n    " + " ".join(f"{db(z[i]):+6.2f}/{deg(z[i]):+5.1f}" for i in pick)
                       + "\n    coh " + " ".join(f"{np.mean([t[c]['coh'][i] for c in (0, 1)]):13.4f}" for i in pick))
    out += ["", "## Nulls, R (L): broadband per window, then per octave at 50-150 / 150-300 / 300-600 ms"]
    for name, d in R["irs"].items():
        for vn in [v for v in VARIANTS if v in d]:
            r = d[vn]
            out.append(f"- {name} {vn}: " + " ".join(f"{w}: {x[1]:.1f} ({x[0]:.1f})" for w, x in r["null"].items() if x))
            out.append("    " + "; ".join(f"{fq}: " + "/".join("--" if x is None else f"{x[1]:.1f}" for x in r["null_oct"][str(fq)][:3]) for fq in OCT))
    print("\n".join(out))
    open(os.path.join(OUT, "report.log"), "w").write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["measure", "compact", "report"], required=True)
    ap.add_argument("--variants", default="plain,chain3", help="measure: comma-separated, from " + ", ".join(VARIANTS))
    o = ap.parse_args()
    if o.stage == "measure":
        measure_stage(o.variants.split(","))
    if o.stage == "compact":
        compact_stage()
    report_stage()


if __name__ == "__main__":
    main()
