"""HF control scan on the Altiverb LONG RT presets (ml_structure 8.6): the L carry-over model (cluster FIR, b from
model_s_results 18.2, shift and gains fitted on the cluster-to-first-pass window) with the effective HF scanned, and the
recirculation 1 / 2 nulls (R) per HF. The presets give no HF (inferred 0.4); the scan finds the HF each preset behaves
as. Run with the floors off so the scanned HF is the one used:

    HS=0.1,0.2,0.3,0.4,0.5,0.6 ML_HFRT=0 ML_L22=0 EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 4 -- \\
        .venv/bin/python tools/altiverb_floor.py 2.0 3.0 4.0     # ALT_LF=1.2 by default
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.signal import sosfiltfilt
import altiverb_predict as ap
import ml_lines as ml
import ml_structure as mls
import ml_carryover as mc
prog = "L"; P = ml.cc.LINES[prog]; M = ml.Model(prog)
LF = float(os.environ.get("ALT_LF", "1.2")); HS = [float(x) for x in os.environ["HS"].split(",")]
B = {"2.0": 0.483, "3.0": 0.155}                      # model_s_results 18.2 (b fitted per preset; 0 at longer RT)
KH = mc.KH; T = 4.0
def device_long(nm, k0_target, n):
    ir = np.load(os.path.join(ap.ROOT, "analysis", "altiverb", "cache", nm + ".npy")).T
    fa = int(np.argmax(np.abs(ir[:, 0]) > np.abs(ir[:, 0]).max() * 10 ** (-40 / 20)))
    h, t0, offs = ap.ai.recover(ir[fa - 240: fa + int(T * ap.FS48)], ap.FS48, ap.ai.FS_DEV)
    k0 = ap.ai.cluster_taps(h)[0]["first_tap_index"]
    h = h * (np.sign(h[k0].sum()) or 1.0)
    h = ap.q1.apply_clock(h, ap.E0.get("eps", 42.8), 0.0)
    return sosfiltfilt(mls.BAND, h[k0 - k0_target: k0 - k0_target + n], axis=0)
W = ml.windows(prog)
for rt in sys.argv[1:]:
    nm = f"LONG RT {rt}s"; key = f"{LF}_{rt}_0.4"; M.b[key] = B.get(rt, 0.0)
    n = ml.K0 + W[1]["recirc3"][1] + 400
    y = device_long(nm, ml.K0, n)[:n]
    c0 = mc.render(prog, np.zeros((2, 0)), np.zeros((2, 0)), mc.cluster_taps(prog), n + 50)[0]
    hs = []
    for ch in (0, 1):
        sl = slice(ml.K0 - 30, ml.K0 + mc.ONSET[prog][ch] - 3 - KH)
        X = np.stack([np.roll(c0[:, ch], k)[sl] for k in range(-KH, KH + 1)], 1)
        hs.append(np.linalg.lstsq(X, y[sl, ch], rcond=None)[0])
    fir = lambda A: np.stack([sum(hs[ch][k + KH] * np.roll(A[:, ch], k) for k in range(-KH, KH + 1)) for ch in (0, 1)], 1)
    C0, R0 = [fir(a) for a in M.parts(n, key, coefs=ml.line_coefs(P, float(rt), LF, 0.4, caps=True))]
    sl, nf = ml.fit_window(prog, True)
    _, u, gs, (Cw, _) = ml.fit_warp(y, [C0, R0], sl, nf, x0=(0.0, 0.0))
    out = []
    for h in HS:
        R = fir(M.parts(n, key, coefs=ml.line_coefs(P, float(rt), LF, h, caps=True))[1]); Rw = ml.warp(R, *u)
        mdl = np.stack([gs[ch][0] * Cw[:, ch] + gs[ch][1] * Rw[:, ch] for ch in (0, 1)], 1)
        out.append((h, [ml.null_db(y[:, 1], mdl[:, 1], ml.K0 + W[1][w][0], min(ml.K0 + W[1][w][1], n - 1)) for w in ("first_pass", "recirc1", "recirc2")]))
    best = min(out, key=lambda t: t[1][1] if t[1][1] is not None else 0)
    print(f"{nm} (b {M.b[key]}, LF {LF}) first_pass/recirc1/recirc2 by HF:", " ".join(f"{h}:{s[0]:.1f}/{s[1]:.1f}/{s[2]:.1f}" for h, s in out), f"| best recirc1 HF {best[0]}", flush=True)
