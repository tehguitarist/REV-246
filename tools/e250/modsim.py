"""Synthetic three-comb system with a modulated delay, measured the way a library capture is (ESS + deconvolution), for the E-program
modulation-depth fit (docs/UAD250_test_plan.md §2.8).

The comb structure is deliberately plain (three combs of the measured lengths, a one-pole loop filter set from the preset's T60 at 250 Hz and
4 kHz); the modulation is a common 0.552 Hz sinusoid on every loop plus a 0.74 Hz one on the longest. What it is for is the *ratio* of comb
coherence at k passes to k = 1 per band, which depends on the accumulated delay change and not on the read-tap detail.
"""
import ctypes, os
import numpy as np
from scipy.optimize import fsolve
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
_lib = ctypes.CDLL(os.path.join(HERE, "libmodcomb.dylib"))
_dp = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
_lib.modcomb.argtypes = [_dp, _dp, ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, _dp, _dp, _dp, ctypes.c_double]
LOOPS_S = (0.16798, 0.19012, 0.21052)   # Altiverb REV-250 loop lengths (EMT250_program.md §3)

def loop_filter(L_s, t60_lo, t60_hi, f_lo=250.0, f_hi=4000.0, r=48000):
    """(g, p): one-pole loop filter whose per-pass loss matches 60 L / T60 dB at f_lo and f_hi."""
    tgt = np.array([60 * L_s / t60_lo, 60 * L_s / t60_hi])
    def loss(v):
        g, p = v; w = 2 * np.pi * np.array([f_lo, f_hi]) / r
        return -20 * np.log10(abs(g) * (1 - p) / np.abs(1 - p * np.exp(-1j * w))) - tgt
    g, p = fsolve(loss, [0.8, 0.2]); return float(g), float(np.clip(p, 0, 0.95))

def simulate(x, r, A_ms, rng, t60_lo=5.67, t60_hi=3.67, lines=((0.552, 1.0, None), (0.74, 0.4, 2))):
    """x: input. A_ms: depth of the 0.552 Hz line (ms). Each line is (f, relative depth, loop index or None for all loops)."""
    y = np.zeros_like(x); ph = {i: rng.uniform(0, 2 * np.pi, len(lines)) for i in range(3)}
    for j, L_s in enumerate(LOOPS_S):
        g, p = loop_filter(L_s, t60_lo, t60_hi, r=r)
        ls = [(f, a, k) for (f, a, k) in lines if k is None or k == j]
        f = np.array([l[0] for l in ls]); A = np.array([l[1] * A_ms * 1e-3 * r for l in ls]); pj = np.array([ph[0][lines.index(l)] for l in ls])
        out = np.zeros_like(x); _lib.modcomb(np.ascontiguousarray(x), out, len(x), L_s * r, g, p, len(f), f, A, pj, float(r)); y += out
    return y
