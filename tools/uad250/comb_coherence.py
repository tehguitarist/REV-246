"""Comb coherence per octave band at k x loop period: single-capture IRs, plugin (known delay modulation) versus Altiverb REV-250.

    .venv/bin/python tools/uad250/comb_coherence.py

Metric of EMT250_program.md §4 made per band: decay-flatten the band-limited tail (0.3 s to 0.3 + 1.2 RT), take its autocorrelation, and read it at
k x the 167 ms loop (k = 1..8, peak within +-(3+k) samples). A modulation that moves the loop delay by d breaks coherence at frequency f by about
(2 pi f d)^2/2, so delay modulation lowers the high bands first and more with k; a gain/rotation modulation (the 246's M/L, flat across frequency)
lowers all bands alike. Sources at RT 4.5: Altiverb REV-250 (one sweep capture), the plugin's impulse snapshots (24, no sweep smear), and the plugin
through sweeps of 5-30 s deconvolved the way a library capture is (smear included).
"""
import os, sys
import numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "tools", "harness"))
from harness.core import band_rt  # noqa: E402
sys.path.insert(0, os.path.join(ROOT, "tools", "uad250"))
import uad250_suite as su  # noqa: E402
r = 48000; BANDS = [250, 500, 1000, 2000, 4000]; RT = 4.5; KS = range(1, 9)

def first(x, th=0.02): return int(np.argmax(np.abs(x) > th * np.abs(x).max()))

def coherence(h, L0_ms):
    """h: (n,) IR at 48 kHz. Returns {band: [rho_k, k=1..8]} and the loop period (ms)."""
    ia = first(h); out = {}; Ls = {}
    for f in BANDS:
        sos = butter(4, [f / 1.414, f * 1.414], "band", fs=r, output="sos"); x = sosfiltfilt(sos, h)
        T = band_rt(h, r, f / 1.414, f * 1.414)["T"] or RT
        s0 = ia + int(0.3 * r); s1 = ia + int((0.3 + 1.2 * RT) * r); t = (np.arange(s0, s1) - ia) / r
        a = x[s0:s1] * 10 ** (3 * t / T); a -= a.mean()
        N = 1 << int(np.ceil(np.log2(2 * len(a)))); ac = np.fft.irfft(np.abs(np.fft.rfft(a, N)) ** 2, N)[:int(1.5 * r)]; ac /= ac[0]
        c0 = int(L0_ms * r / 1000); w = ac[c0 - 70:c0 + 71]; L = c0 - 70 + int(np.argmax(w)); Ls[f] = L
        out[f] = [float(np.max(ac[int(k * L) - (3 + k):int(k * L) + (3 + k) + 1])) for k in KS]
        out[(f, "sd")] = float(np.std(ac[int(0.02 * r):int(0.15 * r)]))
    return out, float(np.mean(list(Ls.values()))) / r * 1e3

def deconv(y, x, r, eps=1e-3):
    N = 1 << int(np.ceil(np.log2(len(y) + len(x)))); X = np.fft.rfft(x, N); Y = np.fft.rfft(y, N)
    H = Y * np.conj(X) / (np.abs(X) ** 2 + eps * (np.abs(X) ** 2).max()); return np.fft.irfft(H, N)

def show(name, res, ch):
    print(f"\n{name} ({ch})  rho at k = 1..8 per band (sd of ACF noise in brackets)")
    for f in BANDS: print(f"  {f:5d} Hz " + " ".join(f"{v:5.2f}" for v in res[f]) + f"   [{res[(f, 'sd')]:.3f}]")

def avg(rs):
    return {**{f: list(np.mean([q[f] for q in rs], 0)) for f in BANDS}, **{(f, "sd"): float(np.mean([q[(f, "sd")] for q in rs])) for f in BANDS}}

if __name__ == "__main__":
    out = {}
    h = np.load(os.path.join(ROOT, "analysis", "altiverb", "cache", "REV-250 4.5s.npy")); f0, f1 = first(h[0]), first(h[1]); Ri = 0 if f0 < f1 else 1
    for nm, c in (("R", Ri), ("L", 1 - Ri)):
        res, L = coherence(h[c][:int(8 * r)], 168.0); out[("Altiverb", nm)] = res; show(f"Altiverb REV-250 4.5 s  (loop {L:.2f} ms)", res, nm)
    # plugin impulse snapshots (24, absolute times 0.1 + 12k s), no sweep smear
    y, _ = sf.read(os.path.join(su.CAP, "C_mod", "train12s_T4.5_noise1.wav"))
    for nm, c in (("R", 1), ("L", 0)):
        rs = []
        for k in range(24):
            a = int(round((0.1 + 12 * k) * r)); rs.append(coherence(y[a:a + int(8 * r), c].astype(float), 167.4)[0])
        out[("plugin imp", nm)] = avg(rs); show("plugin impulse snapshots, mean of 24 (T4.5)", out[("plugin imp", nm)], nm)
    # plugin through sweeps of 5-30 s, deconvolved
    for sec in (5, 10, 20, 30):
        x = su.ess(r, sec=float(sec)); y, _ = sf.read(os.path.join(su.CAP, "C_sweep", f"ess{sec}s_x4_T4.5.wav"))
        for nm, c in (("R", 1), ("L", 0)):
            rs = []
            for k in range(4):
                a = int(0.5 * r) + k * (len(x) + 15 * r); seg = y[a:a + len(x) + 14 * r, c].astype(float)
                ir = deconv(seg, x, r); ir = ir[:int(9 * r)]; rs.append(coherence(ir, 167.4)[0])
            out[(f"plugin ess{sec}s", nm)] = avg(rs); show(f"plugin through a {sec} s sweep, deconvolved, mean of 4 (T4.5)", out[(f"plugin ess{sec}s", nm)], nm)
    print("\nsummary: rho at k=1 / k=4 / k=8, R channel, bands 250 / 1k / 4k")
    for key in [k for k in out if k[1] == "R"]:
        v = out[key]; print(f"  {key[0]:18s} " + " | ".join(f"{b:4d}: {v[b][0]:.2f} {v[b][3]:.2f} {v[b][7]:.2f}" for b in (250, 1000, 4000)))
