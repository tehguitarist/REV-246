"""Best-guess front-panel settings for every Altiverb EMT 246 preset, for the plugin's factory presets.
Results: analysis/altiverb/presets.json (cache/features.json and cache/n_grid.npz are not committed), docs/library_notes.md section 7.

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/altiverb_presets.py [--stage all|features|ngrid|fit]
        [--only <substring> ...] [--workers 3]

Needs the deconvolved IRs of tools/altiverb_compare.py (analysis/altiverb/cache/*.npy, --stage deconv) and its Tomkin
feature cache (cache/tomkin_features.json, --stage settings). The preset files carry only a name; every setting is read
from the IR, relative to the preset families whose RT is in the name (S "NON-LIN SHORT RT", M "MED RT", L "LONG RT",
E "REV-250"; their LF / HF / AMPL from altiverb_compare --stage settings and docs/EMT250_program.md).

  features  per preset: the program from the first 70 ms (cluster / direct / front-tap templates in device samples, the
            L-R reverb lag where there is no cluster), the onset (tap 0, E's R direct tap, N's first front tap, or L's
            reverb onset - 948 samples), octave-band T60 (altiverb_compare.band_t60), the cluster-tap / reverb ratio
            (altiverb_compare.tap_ratio's method at the measured tap 0), E's direct-tap / reverb ratio, and the tap-0
            response against its family's (1-12 kHz) -> features.json.
  ngrid     program N: the model of tools/nl_model.py (device side, no capture terms) rendered over SPAN 0-0.5 (library
            steps, 0.1), EARLY (16 steps) and LATE (0-2.0 by 0.2): 5 ms band-energy envelopes per channel -> n_grid.npz
            (checkpointed; about 1,000 renders).
  fit       the settings:
            DELAY    onset minus the family's onset, plus the family's DELAY (S 10 ms: the CLA IR, DELAY 10 on its display,
                     has its first tap at 10.07 ms, the S family at 10.08; M / E 20; L 39, N 19: the plugin's own DELAY-0 onset
                     convention, against the round 40 / 20; every preset lands on whole ms);
            RT/LF/HF the 246's loop law (ml_lines.line_coefs / rt_of at the program's line lengths; E at its three loops)
                     times a per-band calibration taken from the family (median measured / law T60 at the family's RT
                     label and LF / HF), fitted in log T60 over 125 Hz-8 kHz with the family's per-band scatter as
                     weights; checked by leave-one-out on the families (each family preset fitted with the others'
                     calibration);
            AMPL     S / M / L: the cluster-tap / reverb ratio over Tomkin's ratio-per-AMPL law at the fitted RT / LF / HF
                     (as altiverb_compare --stage settings); E: the direct tap / reverb ratio against the family's at the
                     same RT (log-RT interpolation), the family taken as AMPL 80 (e_plugin_spec controls.AMPL);
            LOW PASS the tap-0 response over the family's, against the plugin's LOW PASS (2nd-order Butterworth at the
                     device rate, rv::dsp::DeviceLowPass): the step with the least rms dB error, 1-12 kHz; the 2-tap
                     device-rate FIR [1 - a, a] that fits the measured excess is reported beside it;
            N        SPAN / EARLY / LATE: least squares on the 5 ms envelopes (dB, 100 Hz-10 kHz, frames within 50 dB of
                     the peak), one gain per channel, against n_grid.
            Each value is snapped to the plugin's panel steps (plugin/src/params/Parameters.cpp) and labelled C / S / H.
"""
import argparse, glob, json, os, sys
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.optimize import minimize
from scipy.signal import butter, find_peaks, resample_poly, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import altiverb_compare as ac  # noqa: E402

FS, FSD = ac.FS, ac.FSD
OUT = os.path.join(ROOT, "analysis", "altiverb")
CACHE = os.path.join(OUT, "cache")
CF = ac.CF

# ---------------------------------------------------------------- the plugin's panel steps (Parameters.cpp)
LP_STEPS = [2, 3, 4, 5, 6, 8, 10, 12]
LF_STEPS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.8, 2.0]
HF_STEPS = [round(0.1 * (i + 1), 1) for i in range(10)]
RT_STEPS = {
    "S": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5],
    "M": [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.3, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0],
    "L": [2.0, 2.2, 2.6, 3.1, 3.6, 4.0, 4.5, 5.2, 6.0, 6.9, 8.0, 10.0, 12.0, 15.0, 17.0, 20.0],
    "E": [0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.3, 3.6, 4.0, 4.5],
}
N_EARLY = [-3.0, -2.0, -1.0, -0.7, -0.5, -0.2, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 10.0]
N_LATE = [round(0.2 * i, 1) for i in range(11)]
N_SPAN_LIB = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]          # library steps (the plugin's 0.05 steps between them are H)
DEFAULT_STEP = {"lowpass": 7, "delay": 10, "ampl": 8, "lf": 9, "rt": 8, "hf": 8}   # joyDefaultStep

# ---------------------------------------------------------------- programs and families
LINES = {"S": [1321, 1389, 1461, 1587], "M": [2641, 2775, 2919, 3171], "L": [5279, 5547, 5837, 6339],
         "E": [2215, 2507, 2776]}
FAMILY = {"S": "NON-LIN SHORT RT ", "M": "MED RT ", "L": "LONG RT ", "E": "REV-250 "}
FAMILY_REF = {"S": "NON-LIN SHORT RT 1.0s", "M": "MED RT 3.0s", "L": "LONG RT 6.0s", "E": "REV-250 2.0s", "N": "NON-LIN"}
# family settings (value, label, reason): altiverb_compare --stage settings, EMT250_program.md 16.5 / 17.3, NonLin 11.7
FAMILY_SET = {
    "S": {"lf": (1.0, "S", "RT-matched Tomkin slices 0.99 / 1.03 / 1.01 (RT 0.6 / 1.0 / 2.2); CLA's panel (same S setting) reads LF 1.0"),
          "hf": (0.9, "S", "RT-matched slices 0.90-0.92, pooled fit 0.88-0.93 at every RT; CLA's panel reads HF 0.9"),
          "ampl": (80, "S", "cluster / reverb ratio on Tomkin's linear AMPL law: 79.7-85.6 %; CLA's panel reads AMPL 80")},
    "M": {"lf": (1.1, "H", "RT-matched slices 1.10-1.20 (RT 1.8-5.0): between the 1.1 and 1.2 steps"),
          "hf": (0.7, "S", "RT-matched slices 0.69-0.71, pooled 0.69-0.75"),
          "ampl": (80, "S", "cluster / reverb ratio: 78.1-85.9 %")},
    "L": {"lf": (1.2, "H", "RT-matched slices 1.19 / 1.20 / 1.30 (RT 6 / 8 / 10); pooled fit unreliable for LF"),
          "hf": (0.4, "S", "RT-matched slices 0.42-0.42, pooled 0.38-0.52"),
          "ampl": (60, "H", "cluster / reverb ratio 57.5-64.0 %, scale resting on Tomkin's single L AMPL (80)")},
    "E": {"lf": (1.49, "S", "fitted LF 1.49 on the 63-500 Hz decay (EMT250_program 16.5), the EMT 250's LF 1.5 step"),
          "hf": (0.69, "S", "fitted HF 0.69 on the 2-12.4 kHz decay, held out 0.6-2.1 % per octave (EMT250_program 17.3)"),
          "ampl": (80, "H", "E's AMPL law is not in the data; the family is taken as AMPL 80 (e_plugin_spec)")},
}
FAMILY_DELAY = {"S": (10, "S"), "M": (20, "S"), "E": (20, "S"), "L": (39, "H"), "N": (19, "H")}
# L / N: the plugin's own DELAY-0 onset (L tap 0 1.13 ms, N front tap 1.57 ms after S's) puts the Altiverb onset at 39 / 19
TIMING_NOTE = {"L": "39 matches the Altiverb's onset in the plugin's timing (L's DELAY-0 onset 1.13 ms after S's), against the round 40",
               "N": "19 matches the Altiverb's onset in the plugin's timing (N's front tap 1.57 ms after S's tap 0), against the round 20"}
FAMILY_DELAY_WHY = ("family onset {on:.2f} ms; S's anchor: the CLA IR (DELAY 10 on its display) has its first tap at 10.07 ms "
                    "and the S family at 10.08; M / L / E / N families sit at x2 / x4 / x2 / x2 of S and every numbered preset "
                    "lands on whole ms from there")
# cluster / direct templates (device samples from the first peak of each channel)
TPL = {"S": ([54, 194, 288], [143, 194, 353]), "M": ([108, 388, 576], [285, 388, 704]),
       "L": ([216, 776, 1152], [570, 776, 1407])}
E_TPL = [429, 718]            # E's R direct taps after the first; L's first 314 later (EMT250_program 11)
N_TPL = ([369, 422], [237, 422])   # N's front taps 41 / 410 / 463 (R 278) from the first
REV_ONSET = {"S": (240, 572), "M": (474, 1134), "L": (948, 2268)}   # reverb onset after cluster tap 0 (unit facts, ml 2)
SCALE = {"S": 1, "M": 2, "L": 4, "E": 2}


def presets():
    return sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(CACHE, "*.npy")))


def load(name):
    return np.load(os.path.join(CACHE, name + ".npy"), mmap_mode="r")


def family_of(name):
    for p, pre in FAMILY.items():
        if name.startswith(pre):
            return p
    return "N" if name == "NON-LIN" else None


# ---------------------------------------------------------------- features
def first_arrival(x, rel=0.01):
    a = np.abs(x)
    return int(np.argmax(a > a.max() * rel))


def peaks_dev(x, t0, span_s=0.07, rel=0.3):
    """Peak positions (device samples, from t0 in 48 kHz samples) of the 8x-upsampled channel in [t0 - 1 ms, + span)."""
    k0 = max(t0 - 48, 0)
    y = resample_poly(np.array(x[k0: k0 + int(span_s * FS)]), 8, 1)
    pk, _ = find_peaks(np.abs(y), height=rel * np.abs(y).max(), distance=int(0.0003 * FS * 8))
    return [((p / 8 + k0) - t0) / FS * FSD for p in pk], [float(y[p]) for p in pk]


def peak_time(x, k):
    """Sub-sample time (48 kHz samples) of the largest peak within +-10 samples of k (8x upsampled)."""
    a = max(k - 24, 0)
    y = resample_poly(np.array(x[a: k + 24]), 8, 1)
    lo, hi = (k - 10 - a) * 8, (k + 10 - a) * 8
    i = lo + int(np.argmax(np.abs(y[lo:hi])))
    return a + i / 8, float(y[i])


def match(pos, tpl, tol=2.0):
    return all(any(abs(p - t) < tol for p in pos) for t in tpl)


def classify(ir):
    """(program, onset in ms, how): from the first 70 ms of both channels."""
    fa = [first_arrival(np.array(ir[c][: int(0.6 * FS)])) for c in (0, 1)]
    pos = []
    for c in (0, 1):
        t, _ = peak_time(ir[c], fa[c] + 6)
        p, _ = peaks_dev(ir[c], int(round(t)), rel=0.2)
        pos.append(([q - (t - round(t)) / FS * FSD for q in p], t))
    (pl, tl), (pr, tr) = pos
    for prog, (tL, tR) in TPL.items():
        if abs(tl - tr) / FS * FSD < 3 and match(pl, tL) and match(pr, tR):
            return prog, min(tl, tr) / FS * 1e3, "cluster taps " + "/".join(map(str, tL)) + " (L), " + "/".join(map(str, tR)) + " (R)"
    dlr = (tl - tr) / FS * FSD
    if abs(dlr - 314) < 3 and match(pr, E_TPL):
        return "E", tr / FS * 1e3, "R direct taps 0/429/718, L 314 samples later"
    if abs(dlr) < 3 and match(pl, N_TPL[0]) and match(pr, N_TPL[1]):
        return "N", min(tl, tr) / FS * 1e3 - 41 / FSD * 1e3, "front taps 41/410/463 (L), 41/278/463 (R)"
    for prog, (rl, rr) in REV_ONSET.items():
        if abs(-dlr - (rr - rl)) < 4 and prog == "L":   # no cluster: the reverb's own L-R lag
            return prog, tl / FS * 1e3 - rl / FSD * 1e3, f"no cluster; reverb onset R - L = {-dlr:.0f} samples (L's 1320)"
    return None, min(tl, tr) / FS * 1e3, f"no 246 S/M/L/E/N template; first peaks L {tl / FS * 1e3:.2f} / R {tr / FS * 1e3:.2f} ms"


def tap_ratio_at(ir, prog, t0):
    """altiverb_compare.tap_ratio with the tap 0 given (48 kHz samples) instead of found."""
    amps = []
    for c in (0, 1):
        x = np.array(ir[c][: int(t0 + 0.6 * FS)])
        y = resample_poly(x[max(t0 - 96, 0): t0 + int(0.12 * FS * ac.SCALE[prog])], 16, 1)
        base = min(t0, 96) * 16
        for off in ac.TEMPLATE[prog]["LR"[c]]:
            k = base + int(round(off / FSD * FS * 16))
            amps.append(np.abs(y[max(k - 8, 0): k + 9]).max())
    w = slice(int(t0 + 0.020 * ac.SCALE[prog] * FS), int(t0 + 0.120 * ac.SCALE[prog] * FS))
    rev = np.sqrt(np.mean([np.mean(np.array(ir[c][w]) ** 2) for c in (0, 1)]))
    return float(np.mean(amps) / rev)


def e_ratio(ir, t0):
    """E: R's first direct tap over the reverb RMS 40-240 ms after it (both channels)."""
    y = resample_poly(np.array(ir[1][t0 - 48: t0 + 48]), 16, 1)
    w = slice(int(t0 + 0.040 * FS), int(t0 + 0.240 * FS))
    rev = np.sqrt(np.mean([np.mean(np.array(ir[c][w]) ** 2) for c in (0, 1)]))
    return float(np.abs(y).max() / rev)


HALF = {"S": 0.9e-3, "M": 1.8e-3, "L": 3.5e-3, "E": 4e-3, "N": 4e-3}


def tap_spectrum(ir, ch, k, half):
    seg = np.array(ir[ch][k - int(half * FS): k + int(half * FS) + 1]) * np.hanning(2 * int(half * FS) + 1) ** 0.5
    return np.abs(np.fft.rfft(seg, 1 << 15))


F_SPEC = np.fft.rfftfreq(1 << 15, 1 / FS)


def features_one(name):
    ir = load(name)
    prog, onset, how = classify(ir)
    r = {"program": prog, "onset_ms": onset, "how": how, "t60": ac.band_t60(np.array(ir[:, : 12 * FS]))}
    if prog in ("S", "M", "L"):
        t0 = int(round(onset * 1e-3 * FS))
        has_cluster = not how.startswith("no cluster")
        r["has_cluster"] = has_cluster
        r["tap_ratio"] = tap_ratio_at(ir, prog, t0) if has_cluster else 0.0
    if prog == "E":
        r["e_ratio"] = e_ratio(ir, int(round(onset * 1e-3 * FS)))
    if prog is not None and (prog != "L" or r.get("has_cluster", True)):
        ch = 1 if prog == "E" else 0
        k = int(round(onset * 1e-3 * FS + (41 / FSD * FS if prog == "N" else 0)))
        k = k + int(np.argmax(np.abs(np.array(ir[ch][k - 3: k + 4])))) - 3
        r["tap_spec"] = tap_spectrum(ir, ch, k, HALF[prog])[: int(12500 / (FS / (1 << 15)))].tolist()
    return name, r


def stage_features(only, workers):
    path = os.path.join(CACHE, "features.json")
    feats = json.load(open(path)) if os.path.exists(path) else {}
    todo = [n for n in presets() if (not only or any(o in n for o in only))]
    with ProcessPoolExecutor(workers) as ex:
        for n, r in ex.map(features_one, todo):
            feats[n] = r
            print(f"{n:24s} {str(r['program']):4s} onset {r['onset_ms']:7.2f} ms  {r['how']}")
    json.dump(feats, open(path, "w"))
    return feats


# ---------------------------------------------------------------- N grid
NFRAME = 0.005
NLEN = 1.6


def n_envelope_model(key):
    import nl_model as nm
    y = nm.render(key, int(NLEN * nm.FS) + 200, o=36.0, proc=False, band=True)   # Altiverb: main taps at +36 (NonLin 11.7)
    k0 = nm.K0
    fl = NFRAME * nm.FS
    out = []
    for c in (0, 1):
        e = [np.sum(y[k0 + int(i * fl): k0 + int((i + 1) * fl), c] ** 2) for i in range(int(NLEN / NFRAME))]
        out.append(10 * np.log10(np.array(e) / fl + 1e-30))
    return key, np.array(out)


NBAND = butter(4, [100, 10000], btype="band", fs=FS, output="sos")


def n_envelope_data(ir, onset_ms):
    """Altiverb N preset: 5 ms band energies from the front frame's zero (first front tap - 41 samples)."""
    k0 = onset_ms * 1e-3 * FS
    out = []
    for c in (0, 1):
        x = sosfiltfilt(NBAND, np.array(ir[c][: int(k0 + (NLEN + 0.2) * FS)]))
        fl = NFRAME * FS
        e = [np.sum(x[int(k0 + i * fl): int(k0 + (i + 1) * fl)] ** 2) for i in range(int(NLEN / NFRAME))]
        out.append(10 * np.log10(np.array(e) / fl + 1e-30))
    return np.array(out)


def stage_ngrid(workers):
    path = os.path.join(CACHE, "n_grid.npz")
    done = {}
    if os.path.exists(path):
        z = np.load(path)
        done = {tuple(k): v for k, v in zip(z["keys"], z["env"])}
    keys = [(s, e, l) for s in N_SPAN_LIB for e in N_EARLY for l in N_LATE if (s, e, l) not in done]
    print(f"N grid: {len(done)} done, {len(keys)} to render")
    with ProcessPoolExecutor(workers) as ex:
        for i, (k, env) in enumerate(ex.map(n_envelope_model, keys, chunksize=4)):
            done[k] = env
            if (i + 1) % 100 == 0 or i + 1 == len(keys):
                np.savez(path, keys=np.array(list(done)), env=np.array(list(done.values())))
                print(f"  {i + 1}/{len(keys)}")
    return done


def n_fit(env_d, grid):
    """Best (SPAN, EARLY, LATE): dB least squares on frames within 50 dB of each channel's peak, one gain per channel."""
    res = []
    for k, env_m in grid.items():
        c = 0.0
        for ch in (0, 1):
            d, m = env_d[ch], env_m[ch]
            sel = d > d.max() - 50
            r = d[sel] - np.maximum(m[sel], m.max() - 90)
            c += np.sum((r - r.mean()) ** 2)
        res.append((c / 2, k))
    res.sort()
    return res


# ---------------------------------------------------------------- the decay law
def law_t60(prog, rt, lf, hf):
    import ml_lines as mlx
    P = LINES[prog]
    return mlx.rt_of(mlx.line_coefs(P, rt, lf, hf), P, CF)


def rt_label(name):
    return float(name.split(" ")[-1][:-1])


def calibration(prog, feats, leave=None):
    """Per-band log(measured / law) over the family (RT >= the second-shortest step; E from RT 1.0): median and scatter."""
    fam = [n for n in feats if family_of(n) == prog and n != leave and n != FAMILY[prog].strip()]
    rts = sorted(rt_label(n) for n in fam)
    fam = [n for n in fam if rt_label(n) >= (1.0 if prog == "E" else rts[1])]
    lf, hf = FAMILY_SET[prog]["lf"][0], FAMILY_SET[prog]["hf"][0]
    R = np.array([np.log(np.array(feats[n]["t60"]) / law_t60(prog, rt_label(n), lf, hf)) for n in fam])
    return np.nanmedian(R, 0), np.maximum(np.nanstd(R, 0), 0.02)


def fit_decay(prog, t60, cal, sd, rt0):
    t = np.log(np.array(t60))
    ok = np.isfinite(t)

    def cost(q):
        pred = cal + np.log(law_t60(prog, np.exp(q[0]), np.exp(q[1]), q[2]))
        return np.sum(((t - pred)[ok] / sd[ok]) ** 2)
    best = min((minimize(cost, [np.log(rt0), np.log(lf0), hf0], method="L-BFGS-B",
                         bounds=[(np.log(RT_STEPS[prog][0] * 0.8), np.log(RT_STEPS[prog][-1] * 1.25)), (np.log(0.1), np.log(2.0)), (0.1, 1.0)])
                for lf0 in (0.5, 1.0, 1.6) for hf0 in (0.3, 0.6, 0.95)), key=lambda r: r.fun)
    return float(np.exp(best.x[0])), float(np.exp(best.x[1])), float(best.x[2]), float(best.fun / max(ok.sum() - 3, 1))


# ---------------------------------------------------------------- AMPL (S / M / L): Tomkin's ratio law
def tomkin_ratio_law():
    feat = json.load(open(os.path.join(CACHE, "tomkin_features.json")))
    law = {}
    for prog in "SML":
        rr = [(ac.parse_tomkin(f), v) for f, v in feat["ratio"].items() if f.split(os.sep)[-2] == ac.TOM_DIR[prog]]
        X = np.array([[1, np.log(rt), np.log(lf), hf] for (amp, lf, rt, hf), v in rr])
        y = np.array([np.log(v / amp) for (amp, lf, rt, hf), v in rr])
        c = np.linalg.lstsq(X, y, rcond=None)[0]
        law[prog] = (c, float(np.std(y - X @ c)))
    return law


# ---------------------------------------------------------------- LOW PASS
def bw2_db(f, fc):
    return -10 * np.log10(1 + (np.tan(np.pi * f / FSD) / np.tan(np.pi * fc * 1e3 / FSD)) ** 4)


def fir2_db(f, a):
    return 10 * np.log10(1 - 2 * a * (1 - a) * (1 - np.cos(2 * np.pi * f / FSD)))


FIR_A_FAMILY = {"S": 0.0, "M": 0.047, "L": 0.23, "E": 0.23, "N": 0.0}   # the program's output kernel (ml_structure 7.5)


def lowpass_fit(spec, ref, prog):
    f = F_SPEC[: len(spec)]
    d = 20 * np.log10(np.array(spec) / np.array(ref))
    d -= np.interp(1000, f, d)
    m = (f > 1000) & (f < 12000)
    errs = {s: float(np.sqrt(np.mean((d[m] - (0 if s == 12 else bw2_db(f[m], s))) ** 2))) for s in LP_STEPS}
    # the step is chosen on the linear amplitude error (the pass band counts, a deep stop band does not)
    lin = {s: float(np.sqrt(np.mean((10 ** (d[m] / 20) - 10 ** ((0 if s == 12 else bw2_db(f[m], s)) / 20)) ** 2))) for s in LP_STEPS}
    best = min(lin, key=lin.get)
    ag = np.arange(0, 0.5, 0.002)
    ce = [np.sqrt(np.mean((d[m] - fir2_db(f[m], a) + fir2_db(f[m], FIR_A_FAMILY[prog])) ** 2)) for a in ag]
    return best, errs, float(ag[int(np.argmin(ce))]), float(min(ce)), {k: float(np.interp(k * 1e3, f, d)) for k in (4, 8, 10, 11)}


# ---------------------------------------------------------------- snapping and naming
def nearest(steps, v, log=False):
    if log:
        return int(np.argmin([abs(np.log(s / v)) for s in steps]))
    return int(np.argmin([abs(s - v) for s in steps]))


def val(value, step, label, reason, inferred=None):
    d = {"value": value, "step": int(step), "label": label, "reason": reason}
    if inferred is not None:
        d["inferred"] = round(float(inferred), 4)
    return d


SHORT = {"PIANO": "Piano", "CATHEDRL": "Cathedral", "BRASS": "Brass", "SM SECT1": "Sm Sect 1", "LG SECT1": "Lg Sect 1",
         "WOODWIND": "Woodwind", "STRINGS": "Strings", "VOCALS": "Vocals", "GATED": "Gated", "PERC": "Percussion",
         "ECHOROOM": "Echo Room", "CHAMBER": "Chamber", "CONCERT": "Concert", "VERY LG": "Very Large"}


def plugin_name(n, prog):
    if n.startswith("PRESET"):
        num, rest = n[6:].split(" ", 1)
        title = SHORT.get(rest) or " ".join(w.capitalize() if w not in ("SM", "LG", "MED", "REF") else
                                            {"SM": "Sm", "LG": "Lg", "MED": "Med", "REF": "Ref"}[w] for w in rest.split())
        return f"Altiverb P{num} {title}"
    if n == "NON-LIN":
        return "Altiverb N Non-Lin"
    if n == "REV-250":
        return "Altiverb REV-250 (bare)"
    rt = rt_label(n)
    return f"Altiverb {prog} RT {rt:g}" if rt >= 10 else f"Altiverb {prog} RT {rt:.1f}"


# ---------------------------------------------------------------- fit stage
def stage_fit(only):
    feats = json.load(open(os.path.join(CACHE, "features.json")))
    grid = None
    gpath = os.path.join(CACHE, "n_grid.npz")
    if os.path.exists(gpath):
        z = np.load(gpath)
        grid = {tuple(float(x) for x in k): v for k, v in zip(z["keys"], z["env"])}
    law = tomkin_ratio_law()
    fam_onset = {p: feats[FAMILY_REF[p]]["onset_ms"] for p in FAMILY_REF}
    cal = {p: calibration(p, feats) for p in "SMLE"}
    out, loo = {}, {}
    # leave-one-out on the families (the check of the decay fit)
    for p in "SMLE":
        lf_f, hf_f = FAMILY_SET[p]["lf"][0], FAMILY_SET[p]["hf"][0]
        for n in [n for n in feats if family_of(n) == p and n != "REV-250"]:
            c, s = calibration(p, feats, leave=n)
            rt, lf, hf, chi = fit_decay(p, feats[n]["t60"], c, s, rt_label(n))
            loo[n] = {"rt_label": rt_label(n), "rt": rt, "lf": lf, "hf": hf, "chi2_dof": chi,
                      "rt_step_ok": nearest(RT_STEPS[p], rt, True) == nearest(RT_STEPS[p], rt_label(n), True),
                      "lf_step": LF_STEPS[nearest(LF_STEPS, lf)], "hf_step": HF_STEPS[nearest(HF_STEPS, hf)]}
            print(f"LOO {n:24s} RT {rt_label(n):5.2f} -> {rt:5.2f}  LF {lf_f:4.2f} -> {lf:4.2f}  HF {hf_f:4.2f} -> {hf:4.2f}  chi2/dof {chi:5.2f}")
    ef = {n: feats[n]["e_ratio"] for n in feats if family_of(n) == "E" and n != "REV-250"}
    e_rt = np.log([rt_label(n) for n in ef])
    e_rr = np.log(list(ef.values()))
    o = np.argsort(e_rt)
    e_rt, e_rr = e_rt[o], e_rr[o]

    for n in sorted(feats):
        if only and not any(q in n for q in only):
            continue
        f = feats[n]
        prog = f["program"]
        fam = family_of(n)
        e = {"altiverb": n, "plugin_name": plugin_name(n, prog), "program": prog, "how": f["how"], "onset_ms": round(f["onset_ms"], 3)}
        if prog is None:
            e["included"] = False
            e["why_not"] = ("a bare two-tap delay (L 20.04 / R 24.99 ms), no reverb: probably the 246's delay program "
                            "(EMT250_program 14.2); not one of the plugin's programs") if n == "REV-250" else \
                           ("a two-channel recirculating echo (L every 99.8 ms x 0.49, R every 65.8 ms x 0.62; T60 about "
                            "1.0 s flat): the 246's delay / echo program, not one of the plugin's programs")
            out[n] = e
            print(f"{n:24s} --  {e['why_not'][:80]}")
            continue
        e["included"] = True
        s = {}
        # DELAY
        fd, fdl = FAMILY_DELAY[prog]
        dv = fd + (f["onset_ms"] - fam_onset[prog])
        st = int(round(dv))
        if fam == prog:
            why = FAMILY_DELAY_WHY.format(on=f["onset_ms"])
            if prog in TIMING_NOTE:
                why = f"family onset {f['onset_ms']:.2f} ms; DELAY " + TIMING_NOTE[prog]
            s["delay"] = val(st, st, fdl, why, dv)
        else:
            s["delay"] = val(st, st, "S" if fdl == "S" else "H",
                             f"onset {f['onset_ms']:.2f} ms = family + {f['onset_ms'] - fam_onset[prog]:.2f} ms "
                             f"(off a whole ms by {dv - st:+.2f}); family DELAY {fd} ({fdl})"
                             + (f"; family {TIMING_NOTE[prog]}, so {st + 1} by the round reading" if prog in TIMING_NOTE else ""), dv)
        # LOW PASS
        if fam == prog:
            s["lowpass"] = val(12, LP_STEPS.index(12), "S", "family preset: its tap 0 is the program's own output kernel "
                               "(the same as on Tomkin, LP 12): open")
        elif "tap_spec" in f:
            best, errs, a, ae, pts = lowpass_fit(f["tap_spec"], feats[FAMILY_REF[prog]]["tap_spec"], prog)
            if best == 12 and abs(pts[10]) < 0.3:
                lab, why = "C", f"tap 0 {pts[10]:+.2f} dB at 10 kHz re the family's: open"
            elif best == 12:
                lab = "H"
                why = (f"tap 0 {pts[8]:+.1f} / {pts[10]:+.1f} / {pts[11]:+.1f} dB at 8 / 10 / 11 kHz re the family, nearer open than "
                       f"any plugin LP step (2nd-order Butterworth); the excess fits a 2-tap device FIR a = {a:.3f} (rms {ae:.2f} dB), "
                       "a roll-off the plugin's LP cannot make")
            else:
                lab = "H"
                why = (f"tap 0 {pts[8]:+.1f} / {pts[10]:+.1f} / {pts[11]:+.1f} dB at 8 / 10 / 11 kHz re the family; "
                       f"nearest plugin LP (2nd-order Butterworth, least amplitude error 1-12 kHz) {best} kHz; "
                       f"the excess fits a 2-tap device FIR a = {a:.3f} (rms {ae:.2f} dB): the unit's LP shape is unknown")
            s["lowpass"] = val(best, LP_STEPS.index(best), lab, why)
            e["lowpass_fit"] = {"rms_db_per_step": errs, "fir2_a": a, "fir2_rms_db": ae, "excess_db": pts}
        else:
            s["lowpass"] = val(12, LP_STEPS.index(12), "H", "no cluster tap to read the LP from: open (the default)")
        if prog in "SMLE":
            # RT / LF / HF
            rt0 = rt_label(n) if fam == prog else float(f["t60"][3])
            rt, lf, hf, chi = fit_decay(prog, f["t60"], *cal[prog], rt0)
            e["decay_fit"] = {"rt": rt, "lf": lf, "hf": hf, "chi2_dof": chi, "t60_measured": f["t60"]}
            if fam == prog:
                rl = rt_label(n)
                s["rt"] = val(rl, nearest(RT_STEPS[prog], rl, True), "C", "from the preset name", rt)
                for k in ("lf", "hf"):
                    v, lab, why = FAMILY_SET[prog][k]
                    steps = LF_STEPS if k == "lf" else HF_STEPS
                    i = nearest(steps, v)
                    s[k] = val(steps[i], i, lab, why + f" (this file alone: {(lf if k == 'lf' else hf):.2f})")
            else:
                i = nearest(RT_STEPS[prog], rt, True)
                lo = RT_STEPS[prog][max(i - 1, 0)]
                hi = RT_STEPS[prog][min(i + 1, 15)]
                gap = min(abs(np.log(rt / lo)), abs(np.log(hi / rt))) / max(abs(np.log(hi / lo)) / 2, 1e-9)
                s["rt"] = val(RT_STEPS[prog][i], i, "S" if chi < 3 else "H",
                              f"decay law fit {rt:.2f} s (1 kHz T60 {f['t60'][3]:.2f} s; chi2/dof {chi:.2f}); nearest step", rt)
                il = nearest(LF_STEPS, lf)
                s["lf"] = val(LF_STEPS[il], il, "H", f"law fit {lf:.2f} from the 125-500 Hz decay against the family's "
                              f"({FAMILY_SET[prog]['lf'][0]}); family LOO scatter below", lf)
                ih = nearest(HF_STEPS, hf)
                s["hf"] = val(HF_STEPS[ih], ih, "S" if chi < 3 else "H", f"law fit {hf:.2f} from the 2-8 kHz decay against the "
                              f"family's ({FAMILY_SET[prog]['hf'][0]})", hf)
            # AMPL
            if prog in "SML":
                if not f.get("has_cluster", True) or f.get("tap_ratio", 0) == 0:
                    s["ampl"] = val(0, 0, "S", "no cluster taps above the reverb (the reverb's own L-R onset lag identifies L)")
                else:
                    c, csd = law[prog]
                    rr = s["rt"]["value"] if fam == prog else rt
                    cp = float(np.exp(np.array([1, np.log(rr), np.log(lf), hf]) @ c))
                    a = f["tap_ratio"] / cp
                    if fam == prog:
                        v, lab, why = FAMILY_SET[prog]["ampl"]
                        s["ampl"] = val(v, v // 10, lab, why + f" (this file: {a:.1f} %)", a)
                    else:
                        i = int(np.clip(round(a / 10), 0, 10))
                        s["ampl"] = val(i * 10, i, "S" if prog != "L" and a < 105 else "H",
                                        f"cluster / reverb ratio on Tomkin's AMPL law: {a:.1f} % (law sd {100 * (np.exp(csd) - 1):.0f} %)"
                                        + ("; L's scale rests on Tomkin's single L AMPL" if prog == "L" else ""), a)
            else:
                r = f["e_ratio"]
                rr = rt_label(n) if fam == prog else rt
                ref = float(np.exp(np.interp(np.log(rr), e_rt, e_rr)))
                a = 80 * r / ref
                if fam == prog:
                    v, lab, why = FAMILY_SET["E"]["ampl"]
                    s["ampl"] = val(v, v // 10, lab, why)
                else:
                    i = int(np.clip(round(a / 10), 0, 10))
                    s["ampl"] = val(i * 10, i, "H", f"direct tap / reverb ratio {r / ref:.2f} x the family's at RT {rr:.2f} "
                                    f"-> {a:.1f} % (family = 80, H)", a)
        else:   # N
            env = n_envelope_data(load(n), f["onset_ms"])
            ranked = n_fit(env, grid)
            (c0, (sp, ea, la)) = ranked[0]
            alt = [k for _, k in ranked[1:6]]
            e["n_fit"] = {"best": [sp, ea, la], "cost": c0, "next": [list(k) for k in alt], "next_cost": [c for c, _ in ranked[1:6]]}
            close = [k for c, k in ranked[1:] if c < 1.15 * c0]
            def lab_of(i):
                return "S" if all(k[i] == (sp, ea, la)[i] for k in close) else "H"
            s["ampl"] = val(sp, int(round(sp / 0.05)), lab_of(0), f"SPAN: envelope fit over the library steps (cost {c0:.0f}; "
                            f"runners-up {alt[:3]})")
            s["lf"] = val(ea, N_EARLY.index(ea), lab_of(1), "EARLY RT: the stages' gain slope, envelope fit")
            s["hf"] = val(la, N_LATE.index(la), lab_of(2), "LATE RT: the flat late decay, envelope fit")
            s["rt"] = val(None, DEFAULT_STEP["rt"], "C", "RT is blank in N; the step is left at the default")
        e["settings"] = s
        e["label_string"] = label_string(prog, s)
        out[n] = e
        print(f"{n:24s} {prog}  " + e["label_string"])
    res = {"method": "tools/altiverb_presets.py (docstring); docs/library_notes.md section 7", "joy_order": ["lowpass", "delay", "ampl", "lf", "rt", "hf"],
           "presets": out, "family_leave_one_out": loo,
           "calibration": {p: {"log_ratio": cal[p][0].tolist(), "sd": cal[p][1].tolist()} for p in cal}}
    path = os.path.join(OUT, "presets.json")
    if only and os.path.exists(path):
        old = json.load(open(path))
        old["presets"].update(out)
        res["presets"] = old["presets"]
    json.dump(res, open(path, "w"), indent=1)
    return res


def label_string(prog, s):
    def g(k, fmt):
        return f"{fmt.format(s[k]['value'])} ({s[k]['label']})"
    if prog == "N":
        parts = [f"SPAN {g('ampl', '{:.2f}')}", f"EARLY {g('lf', '{:g}')}", f"LATE {g('hf', '{:.1f}')}"]
    else:
        rtv = s["rt"]["value"]
        parts = [f"RT {rtv:.1f} ({s['rt']['label']})" if rtv < 10 else f"RT {rtv:g} ({s['rt']['label']})", f"LF {g('lf', '{:.1f}')}", f"HF {g('hf', '{:.1f}')}", f"AMPL {g('ampl', '{:d}')}"]
    parts.append(f"DELAY {g('delay', '{:d}')}")
    if s["lowpass"]["value"] != 12:
        parts.append(f"LP {g('lowpass', '{:d}')}")
    return f"{prog}: " + ", ".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=["all", "features", "ngrid", "fit"])
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    if a.stage in ("all", "features"):
        stage_features(a.only, a.workers)
    if a.stage in ("all", "ngrid"):
        stage_ngrid(a.workers)
    if a.stage in ("all", "fit"):
        stage_fit(a.only)


if __name__ == "__main__":
    main()
