"""Build the deliverable tables from the measured transfers:
  analysis/lowband/chain3_lowband_joint.json   shared joint N+S transfer 12-130 Hz (2 Hz bins), tapered 1 -> T over 12-16 Hz and T -> 1 over 100-130 Hz
  analysis/nl/lowband_N_fine.npz               N's own fine table (fineT_all_sm6.npz of the current model), <= 3.5 kHz (below 100 Hz the shared table takes over in NL_LOWBAND=2)"""
import json, numpy as np, os
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.dirname(os.path.dirname(H))
z = json.load(open(os.environ.get("LB_JOINT_IN") or os.path.join(H, "lowband_joint_T.json")))
f = np.array(z["f_hz"]); m = (f >= 12) & (f <= 130)
w = np.clip(np.minimum((f - 12) / 4.0, (130 - f) / 30.0), 0, 1)[m]
out = {"f_hz": f[m].tolist(), "note": "shared (N+S) data/model complex transfer, tapered; see analysis/nl_lowband30/REPORT.md part 5 (current N model warp2 + branch kernel 0.7; S _apml)"}
for c in (0, 1):
    T = 10 ** (np.array(z[f"joint_db{c}"])[m] / 20) * np.exp(1j * np.radians(np.array(z[f"joint_deg{c}"])[m]))
    T = 1 + (T - 1) * w
    out[f"db{c}"] = (20 * np.log10(abs(T))).tolist(); out[f"deg{c}"] = np.degrees(np.angle(T)).tolist()
os.makedirs(os.path.join(R, "analysis", "lowband"), exist_ok=True)
json.dump(out, open(os.environ.get("LB_JOINT_OUT") or os.path.join(R, "analysis", "lowband", "chain3_lowband_joint.json"), "w"))
n = np.load(os.path.join(H, "fineT_all_sm6.npz")); fr = n["fr"]; keep = fr <= 3500.0
np.savez_compressed(os.path.join(R, "analysis", "nl", "lowband_N_fine.npz"), fr=fr[keep].astype(np.float32), T0=n["T0"][keep].astype(np.complex64), T1=n["T1"][keep].astype(np.complex64))
print("tables written")
