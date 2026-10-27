"""Measured low-band transfer correction of the Tomkin 96 kHz pipeline (analysis/nl_lowband30/REPORT.md parts 3-5).

The captures carry a small fixed transfer error relative to every model that goes through `nl_structure.processing_response`
(the harness's own deconvolution response): 0.9 % at 20-100 Hz falling to 0.1 % at 1 kHz, fine structure with a period of
about 25 Hz (+-0.05 dB, +-0.5 deg). Below about 100 Hz Program N and Program S share it (joint table within 0-1.7 dB of each
program's own, held out on both); above that it is program-specific.

    NL_LOWBAND=0   off
    NL_LOWBAND=1   the shared table, 12-130 Hz (1 below 13 Hz and above 130 Hz)           analysis/lowband/chain3_lowband_joint.json
    NL_LOWBAND=2   the shared table below 100 Hz, Program N's own table above 130 Hz       analysis/nl/lowband_N_fine.npz
                   (linear blend over 100-130 Hz; N only). DEFAULT for Program N since 27 October 2026 (the user's decision).
    SML_LOWBAND=1  S only (sml_fullband.full_model): the shared table; off by default.

`response(f, mode)` -> (len(f), 2) complex data/model transfer (L, R) to multiply onto the model's response. Bench side only:
it is a property of the captures' processing, not of the unit, and does not enter the plugin.
"""
import functools, json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHARED = os.environ.get("LOWBAND_SHARED") or os.path.join(ROOT, "analysis", "lowband", "chain3_lowband_joint.json")   # LOWBAND_SHARED: an alternative table (tests)
NFINE = os.path.join(ROOT, "analysis", "nl", "lowband_N_fine.npz")


def _blend(f):
    """1 up to 100 Hz, linear to 0 at 130 Hz."""
    return np.clip((130.0 - f) / 30.0, 0.0, 1.0)


@functools.lru_cache(maxsize=1)
def _shared():
    z = json.load(open(SHARED))
    f = np.array(z["f_hz"])
    return f, [10 ** (np.array(z[f"db{c}"]) / 20) * np.exp(1j * np.radians(np.array(z[f"deg{c}"]))) for c in (0, 1)]


@functools.lru_cache(maxsize=1)
def _nfine():
    z = np.load(NFINE)
    return z["fr"], [z["T0"], z["T1"]]


def shared(f):
    """The shared transfer on f (Hz): 1 outside 12-130 Hz, tapered at both ends (the table is stored already tapered)."""
    f = np.asarray(f, float)
    fg, T = _shared()
    out = np.ones((len(f), 2), complex)
    m = (f >= fg[0]) & (f <= fg[-1])
    for c in (0, 1):
        out[m, c] = np.interp(f[m], fg, T[c].real) + 1j * np.interp(f[m], fg, T[c].imag)
    return out


def n_fine(f):
    """Program N's own fine table (1 below 100 Hz, tapered to 1 by 3.5 kHz)."""
    f = np.asarray(f, float)
    fg, T = _nfine()
    out = np.ones((len(f), 2), complex)
    m = (f >= fg[0]) & (f <= fg[-1])
    for c in (0, 1):
        out[m, c] = np.interp(f[m], fg, T[c].real) + 1j * np.interp(f[m], fg, T[c].imag)
    return out


def response(f, mode):
    f = np.asarray(f, float)
    if mode == 1:
        return shared(f)
    if mode == 2:
        w = _blend(f)[:, None]
        return w * shared(f) + (1 - w) * n_fine(f)
    return np.ones((len(f), 2), complex)
