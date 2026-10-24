"""Output-domain view of the --base increment: per test setting, the 250 Hz octave output of fold A's and fold B's
increments, their agreement per window, and the split into the loops' common part vs the loop-specific remainder.
model_s_results §25.12 (dense_mid.py --base gain_band200-1500 fits first).

    EMT_KERNEL=chain3 .venv/bin/python tools/resid_common/out250.py <key,key,...>"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt
import dense_mid as dm
dm.SUBSET.update(b0=1.0, rtmin=0.6); dm.FIT_BAND = (200, 1500); dm.BASE = "gain_band200-1500"
dm.SOS_FIT = butter(3, dm.FIT_BAND, btype="band", fs=dm.FS, output="sos")
dm.set_win()
q1 = dm.q1
th = {f: np.load(os.path.join(dm.OUT, f"theta_{dm.tag(f)}.npy")) for f in "AB"}
Lx = [int(p) + q1.EXT for p in q1.P]
off = 8 + 2 * q1.LFIR
n = min(Lx)

def common(x):
    """Each loop's correction replaced by the mean over loops on the first n lags (the part the loops share)."""
    y = np.zeros_like(x); T = sum(Lx)
    for c in (0, 1):
        segs = [x[T * c + sum(Lx[:i]):T * c + sum(Lx[:i]) + n] for i in range(4)]
        m = np.mean(segs, axis=0)
        for i in range(4):
            y[T * c + sum(Lx[:i]):T * c + sum(Lx[:i]) + n] = m
    return y

sos = butter(4, (177, 354), btype="band", fs=dm.FS, output="sos")
W = [(50, 150), (150, 300), (300, 600)]
keys = sys.argv[1].split(",")
for k in keys:
    g = dm.load(k, full=True)
    J = q1.Joint([dict(g, R=g["R_full"])])
    R = np.array(g["R_full"])
    def outp(x):
        _, _, xs = J.split(np.r_[np.zeros(off), x])
        F, X = J._spectra(np.zeros((2, q1.LFIR)), xs, g["N"])
        return sosfiltfilt(sos, np.fft.irfft(np.einsum("cik,ik->ck", X, R), g["N"], axis=1).T[:, 1])
    o = {f: outp(th[f]) for f in "AB"}
    oc = {f: outp(common(th[f])) for f in "AB"}
    c0 = g["caps"][0]
    yres = sosfiltfilt(sos, c0["y"][:, 1] - c0["m"][:, 1])
    line = []
    for a, b in W:
        i0, i1 = dm.PRE + int(a / 1000 * dm.FS), dm.PRE + int(b / 1000 * dm.FS)
        if i1 > min(len(yres), len(o["A"])):
            continue
        cc = lambda u, v: float(u[i0:i1] @ v[i0:i1] / np.sqrt(u[i0:i1] @ u[i0:i1] * v[i0:i1] @ v[i0:i1]))
        en = lambda u: float(u[i0:i1] @ u[i0:i1])
        spec = {f: o[f] - oc[f] for f in "AB"}
        line.append(f"{a}-{b}: A~B {cc(o['A'], o['B']):+.2f}, common A~B {cc(oc['A'], oc['B']):+.2f}, specific A~B {cc(spec['A'], spec['B']):+.2f}, "
                    f"common share {np.mean([en(oc[f]) / en(o[f]) for f in 'AB']):.2f}; vs residual: A {cc(o['A'], yres):+.2f} B {cc(o['B'], yres):+.2f}, "
                    f"common {np.mean([cc(oc[f], yres) for f in 'AB']):+.2f}, specific {np.mean([cc(spec[f], yres) for f in 'AB']):+.2f}")
    print(k, g["lab"]); print("   " + "\n   ".join(line), flush=True)
