"""Empirical low-band correction T(f) for Program N's model (analysis/nl_lowband30): complex data/model transfer ratio per channel on a log grid,
measured on the 96 kHz Tomkin captures.  Use: R_model(f) *= T(f) in nl_structure.processing_response (bench side; the device/plugin is unaffected).

    import lowband_fix as lf;  T = lf.interp(f_hz, ch)       # complex array, 1 outside the table
    lf.build(keys_excluded_span=None)                       # (re)estimate from the bench alignment (needs nl_model + bench dir)
"""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, "lowband_T.json")
_T = None


def grid():
    f1 = 14 * 1.06 ** np.arange(0, 80)
    f1 = f1[f1 <= 400]
    f2 = f1[-1] * 1.12 ** np.arange(1, 20)
    g = np.r_[f1, f2[f2 <= 1600]]
    return g


def interp(f, ch, lo=14.0, hi=1600.0):
    global _T
    if _T is None:
        _T = json.load(open(TABLE))
    F = np.array(_T["f_hz"]); db = np.array(_T[f"db{ch}"]); ph = np.radians(np.array(_T[f"deg{ch}"]))
    ph = np.unwrap(ph)
    x = np.log(np.clip(np.asarray(f, float), 1e-3, None)); xf = np.log(np.r_[F[0] / 1.2, F, F[-1] * 1.2])
    d = np.interp(x, xf, np.r_[0, db, 0], left=0, right=0); p = np.interp(x, xf, np.r_[0, ph, 0], left=0, right=0)
    return 10 ** (d / 20) * np.exp(1j * p)


FINE = os.path.join(HERE, "fineT_all_sm6.npz")


def fine(f, ch):
    """The recommended table: pooled per-bin complex T (Hann-smoothed +-6 bins = +-2.4 Hz), 2^16-point grid at the device rate 26,373.6 Hz; 1 outside 14 Hz-3.5 kHz
    (blended).  f in Hz; returns complex T(f) by linear interpolation of the stored grid."""
    z = np.load(FINE)
    fr = z["fr"]; T = z["T0"] if ch == 0 else z["T1"]
    return np.interp(f, fr, T.real) + 1j * np.interp(f, fr, T.imag)
