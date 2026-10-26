"""Task 2: what is the unlabelled Altiverb `REV-250` preset? (docs/EMT250_program.md section 6)

    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/e250/rev250_unlabelled.py

Brings the preset to the 246's device grid (e_model.to_device, capture clock 42.8 ppm, L's first tap at K0) and tests it against the
labelled E family (REV-250 <RT>s: R taps 0 / 429 / 718, gains +1 / -0.709 / +0.589, L the same 314 later and inverted, all-passes
46 / 146 / 392 g 0.71, loops P = 2215 / 2507 / 2776):
  1. arrival times, gains, polarity per channel; first-tap output kernel against S / MED / LONG / E (lp_chain.py machinery)
  2. search for further taps (high-passed matched filter), bounds at the E positions
  3. tail level per octave and time window against the E family's, and the tail's L / R coherence
  4. E static network: LS gain of the unlabelled tail on the E family's static part (mean of the nine presets), and the all-pass
     deconvolution sparsity (positive control: the E presets)
Reads analysis/altiverb/cache (not committable); derived numbers go to analysis/e250/rev250_unlabelled.json.
"""
import json, os, sys
import numpy as np
from scipy.signal import butter, sosfiltfilt, fftconvolve, find_peaks
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import e_model as em   # noqa: E402
import lp_chain as lc  # noqa: E402

FS, K0 = em.FS, em.K0
OUT = em.OUT
RES = {}


def tap_time(h, ch):
    a = np.abs(h[:, ch]); i = int(np.argmax(a > 0.5 * a.max())); i = i - 3 + int(np.argmax(a[i - 3:i + 4]))
    y0, y1, y2 = h[i - 1:i + 2, ch]
    return i + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)


def hp(x, f=300.0):
    return sosfiltfilt(butter(4, f, btype="high", fs=FS, output="sos"), x, axis=0)


def stage_arrivals(h, y, tL, tR):
    r = RES["arrivals"] = {}
    K = 48000 / FS * (1 + 42.8e-6)
    r["tap_ms"] = {"L": tL / 48.0, "R": tR / 48.0}
    r["R_minus_L_device_samples"] = float((tR - tL) / K)
    r["R_minus_L_ms"] = float((tR - tL) / 48.0)
    r["gain_LF_L_over_R_dB_unlabelled_vs_S"] = None
    fr = np.fft.rfftfreq(lc.NFFT, 1 / lc.FA); m = (fr > 200) & (fr < 1000)
    g = {}
    for fam in ("S", "M", "L", "E", "U"):
        v = []
        for ch in (0, 1):
            x, _, _, _ = lc.family_tap(fam, ch); v.append(np.mean(np.abs(np.fft.rfft(lc.window(x, 1.7))[m])))
        g[fam] = 20 * np.log10(v[0] / v[1])
    r["LF_gain_L_over_R_dB"] = g
    r["peak_sign"] = {"L": float(np.sign(y[K0, 0])), "R": float(np.sign(y[int(round((tR - tL) / K)) + K0, 1]))}
    r["labelled_E_peak_sign"] = {"L": float(np.sign(em.dev(2.0)[K0 + 314, 0])), "R": float(np.sign(em.dev(2.0)[K0, 1]))}
    print("arrivals", json.dumps(r, indent=1))


def stage_kernel():
    r = RES["kernel"] = {}
    T = 1 / FS
    for ch in (0, 1):
        for b in ("S", "M", "L"):
            R, lev, ph = lc.ratio("U", b, ch, 3.0 if b != "S" or ch == 1 else 1.7)
            Rf = R[lc.IFIT]
            fn = lambda q, f: (1 - q[0]) + q[0] * np.exp(-2j * np.pi * f * T)
            fit = least_squares(lambda q: (lambda e: np.concatenate([e.real, np.angle(np.exp(1j * e.imag))]))(np.log(fn(q, lc.FIT) / Rf)), [0.1])
            e = np.log(fn(fit.x, lc.FIT) / Rf)
            band = (lc.FIT > 100) & (lc.FIT < 10000)
            r[f"{'LR'[ch]} U/{b}"] = dict(a_two_tap=float(fit.x[0]), mag_dB_at_8k=float(20 * np.log10(abs(R[np.argmin(abs(lc.FREQS - 8000))]))),
                                          mag_dB_at_10k=float(20 * np.log10(abs(R[np.argmin(abs(lc.FREQS - 10000))]))),
                                          max_abs_dB_100_10k=float(np.max(np.abs(20 * np.log10(np.abs(Rf[band]))))), LF_phase_deg=float(ph))
            print(f"{'LR'[ch]} U/{b}: a={fit.x[0]:+.3f}  8k {r[f'{chr(76) if ch == 0 else chr(82)} U/{b}']['mag_dB_at_8k']:+.2f} dB  max |dB| 100-10k {r[f'{chr(76) if ch == 0 else chr(82)} U/{b}']['max_abs_dB_100_10k']:.3f}  LF phase {ph:+.1f} deg")


def stage_taps(y):
    """Band-pass 300 Hz-6 kHz (no chain skirt, no 8-13 kHz burst), matched filter with the tap's own band-passed waveform
    (peak = 1 at the tap), bounds at the E positions, rms of the matched-filter output as the noise reference."""
    r = RES["taps"] = {}
    yb = sosfiltfilt(butter(4, [300, 6000], btype="band", fs=FS, output="sos"), y, axis=0)
    kk = (K0, K0 + int(round(RES["arrivals"]["R_minus_L_device_samples"])))
    for ch in (0, 1):
        k = kk[ch]; seg = yb[k - 12:k + 40, ch].copy(); tpl = seg / np.linalg.norm(seg)
        full = fftconvolve(yb[:, ch], tpl[::-1], mode="full")            # full[j] = sum_i x[i] tpl[i - j + len - 1]
        mf = full[len(tpl) - 1 - 12:len(tpl) - 1 - 12 + len(yb)]         # mf[n] = <x[n - 12 ... n + 40), tpl>: peak at n = k
        rel = mf / mf[k]
        lag = np.arange(len(rel)) - k
        d = r["LR"[ch]] = {}
        for name, off in (("E tap2 +429", 429), ("E tap3 +718", 718), ("E tap2 +429+314", 429 + 314), ("E tap3 +718+314", 718 + 314), ("E L-tap +314", 314)):
            d[name] = float(np.max(np.abs(rel[k + off - 3:k + off + 4])))
        sel = (lag > int(0.004 * FS)) & (lag < int(1.5 * FS))
        d["largest_after_4ms"] = [(int(lag[t]), float(rel[t])) for t in (np.nonzero(sel)[0][np.argsort(-np.abs(rel[sel]))[:3]])]
        d["mf_rms_4_30ms"] = float(np.sqrt(np.mean(rel[(lag > 0.004 * FS) & (lag < 0.03 * FS)] ** 2)))
        d["mf_rms_30_300ms"] = float(np.sqrt(np.mean(rel[(lag > 0.03 * FS) & (lag < 0.3 * FS)] ** 2)))
        print("taps", "LR"[ch], {k_: (v if not isinstance(v, list) else v[:3]) for k_, v in d.items()})


def band_levels(sig, k, pk, lo, hi, wins):
    x = sosfiltfilt(butter(4, [lo, hi], btype="band", fs=FS, output="sos"), sig)
    return [20 * np.log10(np.sqrt(np.mean(x[k + int(a * FS):k + int(b * FS)] ** 2)) / pk + 1e-30) for a, b in wins]


def stage_tail(y):
    r = RES["tail"] = {}
    wins = [(0.003, 0.01), (0.01, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.0), (1.0, 2.0)]
    bands = [(125, 250), (250, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 8000), (8000, 12500)]
    r["windows_s"] = wins; r["bands_Hz"] = bands
    kk = (K0, K0 + int(round(RES["arrivals"]["R_minus_L_device_samples"])))
    un = {}
    for ch in (0, 1):
        pk = abs(y[kk[ch], ch])
        un["LR"[ch]] = [band_levels(y[:, ch], kk[ch], pk, lo, hi, wins) for lo, hi in bands]
    r["unlabelled_dB_re_tap_peak"] = un
    # the E presets, R channel, window 0.1-0.2 s onwards (their direct taps end at 27 ms): dB re R tap peak
    eref = {}
    for rt in (1.0, 2.5, 4.5):
        d = em.dev(rt); pk = abs(d[K0, 1])
        eref[str(rt)] = [band_levels(d[:, 1], K0, pk, lo, hi, wins) for lo, hi in bands]
    r["E_R_dB_re_tap_peak"] = eref
    # L/R coherence of the tail (R = L delayed by the arrival gap)
    yh = hp(y, 200.0); gap = kk[1] - kk[0]
    coh = {}
    for a, b in ((0.02, 0.04), (0.04, 0.08), (0.08, 0.16), (0.16, 0.3), (0.3, 0.6)):
        i0, i1 = int(a * FS), int(b * FS)
        u = yh[kk[0] + i0:kk[0] + i1, 0]; v = yh[kk[1] + i0:kk[1] + i1, 1]
        coh[f"{a}-{b}s"] = float(np.dot(u, v) / np.sqrt(np.dot(u, u) * np.dot(v, v)))
    r["LR_tail_corr_at_gap_shift"] = coh
    ed = em.dev(2.5)
    cE = {}
    for a, b in ((0.04, 0.08), (0.08, 0.16), (0.16, 0.3)):   # E: L's reverb is R's 315 later, decorrelated random parts
        i0, i1 = int(a * FS), int(b * FS)
        u = hp(ed, 200.0)[K0 + 314 + i0:K0 + 314 + i1, 0]; v = hp(ed, 200.0)[K0 + i0:K0 + i1, 1]
        cE[f"{a}-{b}s"] = float(np.dot(u, v) / np.sqrt(np.dot(u, u) * np.dot(v, v)))
    r["E2.5_LR_corr_at_314_shift"] = cE
    print("tail unlabelled (rows = octave bands 125..8k; cols = windows)\n", json.dumps({k_: np.round(v, 0).tolist() for k_, v in un.items()}))
    print("E 2.5 s, R:", np.round(eref["2.5"], 0).tolist())
    print("coherence unl:", coh, " E:", cE)


def sparsity(x):
    x = np.asarray(x); return float(np.sum(np.abs(x)) / (np.sqrt(len(x)) * np.sqrt(np.sum(x ** 2)) + 1e-300))


def ap_inverse_spec(g=0.71):
    """conj(H) of the three all-passes on the rfft grid em.W (|H| = 1: the inverse is the anti-causal all-pass)."""
    H = np.ones(len(em.W), complex)
    for D in em.AP_D:
        z = np.exp(-2j * np.pi * em.W * D)
        H *= (-g + z) / (1 - g * z)
    return np.conj(H)


def stage_static(y):
    r = RES["static"] = {}
    kk = (K0, K0 + int(round(RES["arrivals"]["R_minus_L_device_samples"])))
    # E static template: mean of the nine labelled presets (device, R tap = 1)
    E = np.mean([em.dev(rt) for rt in em.RTS], 0)
    yh = hp(y, 300.0); Eh = hp(E, 300.0)
    wins = [(0.030, 0.165), (0.165, 0.5), (0.5, 1.5)]
    combos = [("unl L", 0, "E R", 1, K0), ("unl R", 1, "E R", 1, K0), ("unl L", 0, "E L", 0, K0 + 314), ("unl R", 1, "E L", 0, K0 + 314)]
    r["ls_gain"] = {}
    for name, uc, ename, ec, ek in combos:
        row = []
        for a, b in wins:
            i0, i1 = int(a * FS), int(b * FS)
            u = yh[kk[uc] + i0:kk[uc] + i1, uc] / abs(y[kk[uc], uc]); t = Eh[ek + i0:ek + i1, ec]
            g = np.dot(u, t) / np.dot(t, t); rho = np.dot(u, t) / np.sqrt(np.dot(u, u) * np.dot(t, t))
            sg = np.sqrt(np.mean((u - g * t) ** 2) / np.dot(t, t))
            row.append((float(g), float(sg), float(rho)))
        r["ls_gain"][f"{name} on {ename}"] = row
        print(f"{name} on {ename}: " + "  ".join(f"[{a:.3f}-{b:.3f}s] g={g:+.2e}+-{sg:.1e} rho={rho:+.3f}" for (a, b), (g, sg, rho) in zip(wins, row)))
    # positive control: a labelled preset held out, predicted by the mean of the other eight (same measure, same filter)
    ctrl = []
    for rt in (1.0, 2.5, 4.5):
        others = np.mean([em.dev(q) for q in em.RTS if q != rt], 0); Oh = hp(others, 300.0); dh = hp(em.dev(rt), 300.0)
        row = []
        for a, b in wins:
            i0, i1 = int(a * FS), int(b * FS)
            u = dh[K0 + i0:K0 + i1, 1]; t = Oh[K0 + i0:K0 + i1, 1]
            row.append(float(np.dot(u, t) / np.sqrt(np.dot(u, u) * np.dot(t, t))))
        ctrl.append((rt, row))
    r["control_E_LOO_rho_R"] = ctrl
    print("control (an E preset vs the mean of the other eight, R):", ctrl)
    # all-pass deconvolution sparsity: x / (chain kernel), then x / (kernel AP46 AP146 AP392)
    ci = ap_inverse_spec()
    def sp_pair(sig, kern, k, a, b):
        X = np.fft.rfft(sig, em.NFS); Kf = np.fft.rfft(kern, em.NFS)
        Xd = X * np.conj(Kf) / (np.abs(Kf) ** 2 + 1e-3 * np.max(np.abs(Kf)) ** 2)
        base = np.fft.irfft(Xd, em.NFS); apd = np.fft.irfft(Xd * ci, em.NFS)
        i0, i1 = k + int(a * FS), k + int(b * FS)
        return sparsity(base[i0:i1]), sparsity(apd[i0:i1])
    r["sparsity_before_after_AP_division"] = {}
    for a, b in ((0.03, 0.3), (0.3, 1.0)):
        out = {}
        # E mean (static) and one E preset (positive controls), kernel = E R tap-1 response (em.kernels)
        kern = em.kernels(em.RTS)[1]
        sig = np.zeros(em.NFS); sig[:len(E)] = E[:, 1]
        out["E mean R"] = sp_pair(sig[:len(E)], kern, K0 - 16, a, b)
        d = em.dev(2.5)[:, 1]; out["E 2.5 R"] = sp_pair(d, kern, K0 - 16, a, b)
        # unlabelled: kernel = its own first tap (tap -16 .. +410), from the high-passed? no: raw
        for ch in (0, 1):
            k = kk[ch]; kern_u = y[k - 16:k + em.KL, ch].copy(); kern_u = kern_u / kern_u[16]
            out["unlabelled " + "LR"[ch]] = sp_pair(y[:, ch] / y[k, ch], kern_u, k - 16, a, b)
        r["sparsity_before_after_AP_division"][f"{a}-{b}s"] = out
        print(f"sparsity (L1/(sqrt(N) L2)) before -> after dividing the E all-passes, window {a}-{b}s after tap:", {k_: tuple(round(x, 3) for x in v) for k_, v in out.items()})


def main():
    h = np.load(os.path.join(em.CACHE, "REV-250.npy")).T
    tL, tR = tap_time(h, 0), tap_time(h, 1)
    n = int(6 * FS)
    y = em.to_device(h, 42.8, tL - K0 * 48000 / FS, n)
    lc.FAM["U"] = ["REV-250"]; lc.NEXT[("U", 0)] = 48000; lc.NEXT[("U", 1)] = 48000; lc.NDEV[("U", 0)] = 528; lc.NDEV[("U", 1)] = 659
    stage_arrivals(h, y, tL, tR)
    stage_kernel()
    stage_taps(y)
    stage_tail(y)
    stage_static(y)
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "rev250_unlabelled.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
