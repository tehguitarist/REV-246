"""Joint (f_m, A) fit of the E program's delay modulation from the Altiverb REV-250 4.5 s comb coherence (docs/UAD250_test_plan.md §2.8).

    .venv/bin/python tools/e250/depth_grid.py [--rebuild]

Table: for a single common sinusoidal delay modulation (f_m, A) on all three combs, the per-band ratio rho_k / rho_1, k = 2..8, through a 20 s ESS
and deconvolution (mean of 3 random phases). Fit: least squares of the data ratios (R and L) to the table, bands 1 / 2 / 4 kHz (500 Hz and below
are noise-limited), with a self-test that recovers known (f_m, A) from simulated data of a different seed.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "tools", "uad250")); sys.path.insert(0, HERE)
import comb_coherence as cc, uad250_suite as su, modsim  # noqa: E402
r = 48000; BANDS = [1000, 2000, 4000]; KS = list(range(2, 9))
FS = [0.2, 0.3, 0.4, 0.552, 0.74, 1.0, 1.5, 2.0, 3.0, 5.0]; AS = [0.01, 0.02, 0.03, 0.045, 0.07, 0.1, 0.15, 0.25]
TAB = os.path.join(ROOT, "analysis", "e250", "depth_table.npz")

def vec(res): return np.array([[res[f][k - 1] / max(res[f][0], 1e-3) for k in KS] for f in BANDS])

def sim_vec(fm, A, seed, nrep=3, sec=20):
    rng = np.random.default_rng(seed); x = su.ess(r, sec=float(sec)); out = []
    for _ in range(nrep):
        sig = np.zeros(len(x) + 14 * r); sig[:len(x)] = x
        y = modsim.simulate(sig, r, A, rng, lines=((fm, 1.0, None),))
        out.append(cc.coherence(cc.deconv(y, x, r)[:int(9 * r)], 167.98)[0])
    return vec(cc.avg(out))

def build():
    T = np.zeros((len(FS), len(AS), len(BANDS), len(KS)))
    for i, fm in enumerate(FS):
        for j, A in enumerate(AS):
            T[i, j] = sim_vec(fm, A, seed=100 + i * 10 + j); print(f"table f_m {fm} A {A}", flush=True)
    np.savez(TAB, T=T, FS=FS, AS=AS, BANDS=BANDS, KS=KS)

def fit(v, T, w=None):
    err = ((T - v) ** 2).sum((2, 3)); i, j = np.unravel_index(np.argmin(err), err.shape)
    return i, j, err

if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(TAB): build()
    z = np.load(TAB); T = z["T"]
    # self-test: truths simulated with a different seed
    print("self-test (recover truth from a fresh simulation, joint f_m and A):")
    for fm, A in [(0.552, 0.045), (1.0, 0.07), (0.3, 0.03), (3.0, 0.02)]:
        v = sim_vec(fm, A, seed=999); i, j, err = fit(v, T)
        order = np.argsort(err.ravel())[:3]
        print(f"  truth f_m {fm} A {A}: best f_m {FS[i]} A {AS[j]} (rms err {np.sqrt(err[i, j] / v.size):.3f}); next best " +
              ", ".join(f"({FS[o // len(AS)]},{AS[o % len(AS)]})" for o in order[1:]))
    h = np.load(os.path.join(ROOT, "analysis", "altiverb", "cache", "REV-250 4.5s.npy")); Ri = 0 if cc.first(h[0]) < cc.first(h[1]) else 1
    print("\nAltiverb REV-250 4.5 s (data ratios rho_k/rho_1, k=2..8 at 1 / 2 / 4 kHz):")
    dat = {}                                     # the two measured channel vectors, then their mean
    for nm, c in (("R", Ri), ("L", 1 - Ri), ("mean", None)):
        if c is None:
            v = (dat["R"] + dat["L"]) / 2
        else:
            v = vec(cc.coherence(h[c][:int(8 * r)], 168.0)[0])
        dat[nm] = v
        i, j, err = fit(v, T); rms = np.sqrt(err / v.size)
        good = [(FS[a], AS[b], round(float(rms[a, b]), 3)) for a, b in zip(*np.where(rms <= rms[i, j] * 1.25))]
        print(f"  {nm:4s}: best f_m {FS[i]} Hz, A {AS[j]} ms ({AS[j] * 1e-3 * 26373.6:.1f} device samples), rms err {rms[i, j]:.3f}; within 25 % of best: {good[:8]}")
    # best A per f_m for the mean data
    err = ((T - dat["mean"]) ** 2).sum((2, 3)); print("\nbest A per f_m (mean of R and L), rms error:")
    for a, fm in enumerate(FS): b = int(np.argmin(err[a])); print(f"  f_m {fm:5.3f} Hz: A {AS[b]:.3f} ms ({AS[b] * 1e-3 * 26373.6:4.1f} dev smp), rms {np.sqrt(err[a, b] / dat['mean'].size):.3f}")
