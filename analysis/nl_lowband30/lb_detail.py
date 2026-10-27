"""Spectrum / phase / time course of the front-section residual e = y - m per capture group (session = front-shift cluster).
usage: lb_detail.py [maxperGroup]"""
import sys, collections, pickle
from lb_common import *
B = bench_all()
maxper = int(sys.argv[1]) if len(sys.argv) > 1 else 24
# groups by front shift (rounded cluster), untrimmed, SPAN-agnostic
def grp(fs):
    for name, c in (("fs0", -0.13), ("fs16", 16.35), ("fs19.6", 19.65), ("fs22.9", 22.94), ("fs36", 36.13)):
        if abs(fs - c) < 0.1: return name
    return None
G = collections.defaultdict(list)
for k, r in B.items():
    if r["trimmed"] or k[1] < 0.3: continue
    g = grp(r["front_shift"])
    if g: G[g].append(k)
rng = np.random.default_rng(1)
sel = {}
for g, ks in G.items():
    ks = sorted(ks)
    sel[g] = [ks[i] for i in sorted(rng.choice(len(ks), min(maxper, len(ks)), replace=False))]
    print(g, len(ks), "selected", len(sel[g]), flush=True)
NF = 8192
freqs = np.fft.rfftfreq(NF, 1 / FS)
data = {}
for g, ks in sel.items():
    for k in ks:
        r = B[k]
        y, m, _ = load(k, r, n=2200)
        z = K0 + r["front_shift"]
        b = int(z + 0.025 * FS)           # 25 ms: just before the first main tap (685)
        a = 0
        tap = np.ones(b - a)
        t = 55                            # 2 ms taper at the end
        tap[-t:] = 0.5 * (1 + np.cos(np.pi * np.arange(t) / t))
        out = {}
        for c in (0, 1):
            Y = np.fft.rfft(y[a:b, c] * tap, NF); Mo = np.fft.rfft(m[a:b, c] * tap, NF)
            out[c] = (Y, Mo)
        data[(g, k)] = dict(out=out, z=z, e=(y - m)[:, :], m=m, y=y, fs=r["front_shift"])
pickle.dump({"data": data, "freqs": freqs}, open("detail.pkl", "wb"))
print("saved", len(data))
