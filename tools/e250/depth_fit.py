"""Fit the E program's delay-modulation depth: comb coherence versus passes, per band (docs/UAD250_test_plan.md §2.8).

    .venv/bin/python tools/e250/depth_fit.py

Simulates the three-comb system with common-line depth A through ESS sweeps of 10 / 20 / 30 s (several repeats, random modulation phase),
deconvolves, measures comb coherence with tools/uad250/comb_coherence.py, and compares the ratio rho_k / rho_1 per band with (a) the plugin's own sweeps
(known depth about 0.4 ms) as a check of the method and (b) the Altiverb REV-250 4.5 s capture.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "tools", "uad250")); sys.path.insert(0, HERE)
import comb_coherence as cc, uad250_suite as su, modsim  # noqa: E402
r = 48000; B = cc.BANDS

def measure(A_ms, sec, nrep, seed=0, ir_len=9):
    rng = np.random.default_rng(seed); x = su.ess(r, sec=float(sec)); res = []
    for k in range(nrep):
        sig = np.zeros(len(x) + 14 * r); sig[:len(x)] = x
        y = modsim.simulate(sig, r, A_ms, rng)
        ir = cc.deconv(y, x, r)[:int(ir_len * r)]; res.append(cc.coherence(ir, 167.98)[0])
    return cc.avg(res)

def ratios(res):  # rho_k / rho_1 for k = 2, 4, 8 per band
    return {f: [res[f][k - 1] / max(res[f][0], 1e-3) for k in (2, 4, 8)] for f in B}

if __name__ == "__main__":
    h = np.load(os.path.join(ROOT, "analysis", "altiverb", "cache", "REV-250 4.5s.npy")); Ri = 0 if cc.first(h[0]) < cc.first(h[1]) else 1
    alt = {nm: cc.coherence(h[c][:int(8 * r)], 168.0)[0] for nm, c in (("R", Ri), ("L", 1 - Ri))}
    print("data  rho_1 (R/L) and rho_k/rho_1 at k=2,4,8, Altiverb:")
    for f in B: print(f"  {f:5d} Hz  rho_1 {alt['R'][f][0]:.2f}/{alt['L'][f][0]:.2f}  R: " + " ".join(f"{alt['R'][f][k - 1] / alt['R'][f][0]:.2f}" for k in (2, 4, 8)) + "   L: " + " ".join(f"{alt['L'][f][k - 1] / alt['L'][f][0]:.2f}" for k in (2, 4, 8)))
    for sec in (10, 20, 30):
        print(f"\n== simulated sweep {sec} s: rho_1 and ratios at k=2,4,8 (mean of 6 random phases) for common-line depth A")
        print("   A ms | " + " | ".join(f"{f:5d} Hz: rho1  k2   k4   k8 " for f in (500, 1000, 2000, 4000)))
        for A in (0.0, 0.02, 0.05, 0.1, 0.2, 0.4):
            res = measure(A, sec, 6); rt = ratios(res)
            print(f"  {A:5.2f} | " + " | ".join(f"        {res[f][0]:.2f} " + " ".join(f"{v:4.2f}" for v in rt[f]) for f in (500, 1000, 2000, 4000)), flush=True)
