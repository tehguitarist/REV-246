"""Does the 100 Hz-1 kHz model error live in per-stage gain / timing?  Per capture and channel regress the spectral residual (Hann, z+0.05..0.45 s, 100-1000 Hz)
on the model's stage columns (gain) and their time-shift derivatives; report in-sample null, and the consistency of the fitted relative stage-gain deviations
across EARLY; leave-one-EARLY-out with SHARED relative deviations (14 gain + 14 shift numbers per channel)."""
import sys, collections
from lb_common import *
B = bench_all()
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS); lo, hi = 100, 1000; mk = (fr >= lo) & (fr < hi)
keys = sorted(k for k in B if k[0] == 0.5 and k[1] >= 0.7 and k[2] <= 1.0 and not B[k]["trimmed"] and abs(B[k]["o"] - 1.0) < 0.05 and abs(B[k]["front_shift"] - 22.94) < .3)
print(len(keys), "captures", sorted({k[1] for k in keys}))
data = {}
for k in keys:
    r = B[k]; z = K0 + r["front_shift"]; a, b = int(z + 0.05 * FS), int(z + 0.45 * FS)
    n = b + 64
    y, m, _ = load(k, r, n=n, band=False)
    N = 1 << int(np.ceil(np.log2(n + 2 * ne.TAPS[-1] + 4096)))
    win = np.hanning(b - a)
    Y = np.fft.rfft(y[a:b] * win[:, None], NF, axis=0)
    for c in (0, 1):
        cols = stage_columns(k, c, N, r)
        # model column time series, windowed spectra
        sp = {}
        for nm, S in cols.items():
            x = np.fft.irfft(S, N)[a:b] * r["gain"][c]
            sp[nm] = np.fft.rfft(x * win, NF)[mk]
        data[(k, c)] = (Y[mk, c], sp)
print("built")
names = [f"s{j}" for j in range(14)]
def design(sp, kind):
    cols = []
    for nm in (["front"] + names):
        if nm in sp: cols.append(sp[nm])
    if kind == "gd":
        for nm in names:
            if nm + "_d" in sp: cols.append(sp[nm + "_d"])
    return cols
for c in (0, 1):
    # in-sample per capture: free stage gains only
    out = []; dev = collections.defaultdict(list)
    for (k, cc), (Y, sp) in data.items():
        if cc != c: continue
        base = sum(sp[nm] for nm in sp if not nm.endswith("_d"))
        e0 = 10 * np.log10(np.sum(abs(Y - base) ** 2) / np.sum(abs(Y) ** 2))
        cl = [sp[nm] for nm in (["front"] + names) if nm in sp]
        A = np.stack(cl, 1); Ar = np.vstack([A.real, A.imag]); br = np.r_[Y.real, Y.imag]
        x = np.linalg.lstsq(Ar, br, rcond=None)[0]
        e1 = 10 * np.log10(np.sum(abs(Y - A @ x) ** 2) / np.sum(abs(Y) ** 2))
        out.append((e0, e1))
        for i, nm in enumerate(["front"] + [n_ for n_ in names if n_ in sp]): dev[(k[1], nm)].append(x[i] - 1)
    out = np.array(out)
    print(f"ch{c}: 100-1k null  model {np.median(out[:,0]):6.1f}  ->  free per-stage gains per capture (in-sample, {len(cl)} numbers) {np.median(out[:,1]):6.1f}")
    print("   fitted relative stage-gain deviation (median over captures), by EARLY (rows) and stage:")
    for e in sorted({k[1] for k in keys}):
        row = [np.median(dev[(e, nm)]) * 1e3 if (e, nm) in dev else np.nan for nm in ["front"] + names]
        print(f"    EARLY {e:4.1f}: " + " ".join(f"{v:6.1f}" for v in row), " (x1e-3)")
