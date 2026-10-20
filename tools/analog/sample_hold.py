"""The two sample-and-hold stages of the EMT 246, as transcribed from the 246S service manual.

Input S&H (analog input board 7 246 121, sheet 2/3, printed p. 41; text p. 37; timing p. 39):
    Z102 (buffer after the input filter) -> R114 787 R -> node A.
    Node A -> TL191 switches 9-8 and 16-1 (closed while integrating) -> Z104 (LF356) inverting input,
    C105 1 nF in the feedback. Feedback T: Z104 out -> R116 1 k -> node B (C104 1 nF to ground)
    -> R115 1 k -> node A. While holding, switches 3-4 and 6-5 short node A to ground and the
    integrator input is open, so its output holds. Gain -2 R115/R114 = -2.54.
    Window: ADSHA high = integrate, 9.1 us of the 37.92 us period (about 44 of 182 timing steps),
    read off the logic-analyser plot on p. 39; the SAR clocks run while ADSHA is low (hold).
    Then C106 10 uF (tantalum) -> R117 4.75 k -> Z112 buffer -> summing point: a 3.35 Hz high-pass.

Output S&H (analog output board 7 246 122, sheet 2/3, printed p. 48), for reference:
    DAC (I/U Z207) -> R202 1.21 k -> node X; R203 1.21 k from Z209's output back to X;
    X -> P202 (charge current, 5 k) -> LF13333 switch -> Z209 (LF356) with C205 1 nF.
    While integrating: tau = (R203 + 2 P202) C205 (R202 = R203), gain -1. The DAC value is
    constant inside the window, so the stage is the discrete one-pole y[n] = p y[n-1] + (1-p) x[n].

All functions return complex responses at the frequencies f (Hz), normalised to unity at DC.
"""
import numpy as np
from scipy.linalg import expm

FS_DEV = 19.2e6 / 728
T_DEV = 1 / FS_DEV
STEP = 4 / 19.2e6                      # one timing-PROM step (CLKOUT A), 208.3 ns

k, n = 1e3, 1e-9
IN_SH = {"R114": 787.0, "R115": 1 * k, "R116": 1 * k, "C104": 1 * n, "C105": 1 * n,
         "Ti": 9.1e-6, "Ron": 0.0}
IN_COUPLING_HZ = 1 / (2 * np.pi * 4.75e3 * 10e-6)        # C106 / R117, 3.35 Hz
IN_AMP_COUPLING_HZ = 1 / (2 * np.pi * (1e6 + 100e3 * 47.5e3 / 147.5e3) * 680e-9)  # C101 / R105 (+R101||R102), 0.23 Hz


def _input_state(c):
    """x = [Vo, VB]; dx/dt = A x + b Vin while integrating (node A through the switches, R_on, to virtual ground)."""
    ron = max(c["Ron"], 1e-6)
    g = 1 / c["R114"] + 1 / c["R115"] + 1 / ron
    a = 1 / (c["C105"] * ron * g)                    # dVo/dt = -a (Vin/R114 + VB/R115)
    A = np.array([[0.0, -a / c["R115"]],
                  [1 / (c["C104"] * c["R116"]),
                   (-(1 / c["R116"] + 1 / c["R115"]) + 1 / (g * c["R115"] ** 2)) / c["C104"]]])
    b = np.array([-a / c["R114"], 1 / (g * c["R115"] * c["R114"] * c["C104"])])
    return A, b


def input_sh(f, Ti=None, Ron=None, comp=IN_SH):
    """Sampled response (continuous input -> held sample sequence) of the integrate/hold input stage.

    Each period: VB starts settled at Vo/2 (the hold lasts about 29 us against R115||R116 * C104 = 0.5 us),
    then the window of length Ti is solved exactly. Vo[n] = c Vo[n-1] + G(f) e^{j w t_n}, so
    H(f) = G(f) / (1 - c e^{-j w T}). Aliasing is ignored (the input filter removes it).
    """
    c = dict(comp)
    if Ti is not None: c["Ti"] = Ti
    if Ron is not None: c["Ron"] = Ron
    A, b = _input_state(c)
    Ti = c["Ti"]
    x0 = np.array([1.0, 0.5])
    cpole = float((expm(A * Ti) @ x0)[0])
    f = np.atleast_1d(np.asarray(f, float))
    H = np.empty(len(f), complex)
    I = np.eye(2)
    for i, fr in enumerate(f):
        s = 2j * np.pi * fr
        M = A - s * I                                # G = [(A - sI)^-1 (e^{(A - sI) Ti} - I) b]_0
        G = np.linalg.solve(M, (expm(M * Ti) - I) @ b)[0]
        H[i] = G / (1 - cpole * np.exp(-s * T_DEV))
    s0 = np.linalg.solve(A, (expm(A * Ti) - I) @ b)[0] / (1 - cpole)
    return H / s0, cpole


def output_sh(f, pole):
    """Output integrate/hold as the discrete one-pole at the device rate (DAC held through the window)."""
    w = 2 * np.pi * np.asarray(f, float) / FS_DEV
    return (1 - pole) / (1 - pole * np.exp(-1j * w))


def zoh_db(f, ref=1e3):
    """Zero-order hold of the output staircase (value held between S&H windows), dB re ref."""
    f = np.asarray(f, float)
    return 20 * np.log10(np.abs(np.sinc(f / FS_DEV)) / abs(np.sinc(ref / FS_DEV)))


def highpass(f, corners):
    """Cascade of first-order analogue high-passes (corners in Hz)."""
    f = np.asarray(f, float)
    H = np.ones(f.shape, complex)
    for fc in corners:
        H = H * (1j * f / fc) / (1 + 1j * f / fc)
    return H
