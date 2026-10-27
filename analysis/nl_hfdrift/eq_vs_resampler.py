import sys, json, numpy as np
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
from common import *
from scipy.signal import firwin
H = json.load(open(OUT + "/heldout.json")); knots = np.array(H["knots"])
def gain(up, down, fs_in, f_hz):
    """Gain of resample_poly(up, down)'s default FIR (Kaiser 5.0, half_len 10*max) at input-rate frequency f_hz, DC-normalised."""
    mr = max(up, down); h = firwin(2 * 10 * mr + 1, 1.0 / mr, window=("kaiser", 5.0))
    N = 1 << 22; Hf = np.abs(np.fft.rfft(h, N)); fr = np.fft.rfftfreq(N) * fs_in * up
    return np.interp(f_hz, fr, Hf) / Hf[0]
f = knots[:-1]
gdg = 20 * np.log10(gain(25, 91, 96000, f))          # to_device_grid's polyphase filter (96 kHz -> device rate)
guw = 20 * np.log10(gain(32, 1, FS, f))              # default unwarp's upsampling kernel (static part; the real unwarp effect is time dependent)
eq_bench = np.interp(f, *(lambda e: (e["f_hz"], e["db"]))(json.load(open("/Users/leighpierce/Git/246/analysis/nl/out_eq.json"))))
g_new = 20 * np.log10(np.abs(np.array(H["g_new_tr"])[:-1])); g_new_te = 20 * np.log10(np.abs(np.array(H["g_new_te"])[:-1]))
cn = np.array([complex(a, b) for a, b in H["c_new_tr"]])[:-1]
print("f kHz            ", np.round(f / 1000, 2))
print("to_device_grid FIR", np.round(gdg, 2))
print("fitted EQ new-unwarp train (real)", np.round(g_new, 2))
print("fitted EQ new-unwarp test  (real)", np.round(g_new_te, 2))
print("fitted EQ complex |.| train        ", np.round(20 * np.log10(np.abs(cn)), 2))
print("fitted minus to_device_grid (dB), train real:", np.round(g_new - gdg, 2))
print("                              test real :", np.round(g_new_te - gdg, 2))
print("bench out_eq.json interp          ", np.round(eq_bench, 2))
print("default unwarp kernel (static)   ", np.round(guw, 2))
