"""Per capture: measured step c (cache - LTI model, pre-roll 8:40), simulated finite-record effect of the model (sim_trunc - sim_long),
step produced by the file's own pre-roll DC offset through the pipeline (c_dc), in 1e-4 of the front-tap peak."""
import sys, json, time; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
from sim import _SOS
np.set_printoptions(linewidth=220, precision=2, suppress=True)
keys = [tuple(float(v) for v in a.split(",")) for a in sys.argv[2:]]
out = {}
for key in keys:
    t0 = time.time()
    r = align(key)
    x, fs = sf.read(file_path(key), always_2d=True); nf = len(x)
    cache = nl.ir(key); nlti = len(cache)
    xs0 = make_x(key, r, nf, 1.0); A = fit_scale(x, xs0)
    x_tr = make_x(key, r, nf, A); x_long = make_x(key, r, nf, A, 6.0)
    mp = M.render(key, nlti, r["o"], r["front_shift"], band=False, proc=True) * np.array(r["gain"])
    lti = forward_warp(mp, r["ppm"])
    z = int(round(K0 + r["front_shift"])); u = 1e4 / np.abs(mp[z:z + 800, 0]).max()
    sl = slice(8, 40)
    c_meas = ((cache - lti)[sl].mean(0)) * u
    h_tr, _ = pipeline(x_tr); h_lg, _ = pipeline(x_long)
    c_fin = ((h_tr[:nlti] - h_lg[:nlti])[sl].mean(0)) * u
    # first arrival in the sim file
    xb = sosfiltfilt(_SOS, x_tr, axis=0); a0 = int(np.argmax(np.abs(xb).max(1) > 1e-3 * np.abs(xb).max())) - 20
    n0 = max(30, a0)
    dc = x[:n0].mean(0)
    h_dc, _ = pipeline(np.ones((nf, 2)) * dc)
    c_dc = (h_dc[:nlti][sl].mean(0)) * u
    # step by 1 unit DC for the unit response
    # LF residual (< 100 Hz) of the whole file through the pipeline
    xr = x - x_tr; X = np.fft.rfft(xr, axis=0); f = np.fft.rfftfreq(nf, 1 / FS96)
    lo = np.fft.irfft(X * (f < 100)[:, None], nf, axis=0)
    h_lo, _ = pipeline(lo); c_lo = (h_lo[:nlti][sl].mean(0)) * u
    out["_".join(map(str, key))] = dict(c_meas=c_meas.tolist(), c_fin=c_fin.tolist(), c_dc=c_dc.tolist(), c_lo=c_lo.tolist(), dc=dc.tolist(), n0=n0, nf=nf, fr=r["front_shift"])
    print(key, f"nf {nf/FS96:.2f}s fr {r['front_shift']:.1f} n0 {n0} dc {dc*1e6} uV-ish(1e-6) | c_meas {c_meas} | c_fin(sim finite record) {c_fin} | c_dc {c_dc} | c_lf<100 {c_lo} | {time.time()-t0:.0f}s", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
