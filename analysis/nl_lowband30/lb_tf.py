import pickle, collections, numpy as np
FS = 26373.626373626373
S = pickle.load(open("low_sel.pkl", "rb"))
NF = 1 << 16; fr = np.fft.rfftfreq(NF, 1 / FS)
def spec(d, t0=0.05, t1=0.45):
    a, b = int(d["z"] + t0 * FS), int(d["z"] + t1 * FS); w = np.hanning(b - a)[:, None]
    return np.fft.rfft(d["y"][a:b] * w, NF, axis=0), np.fft.rfft(d["m"][a:b] * w, NF, axis=0)
sp = {k: spec(d) for k, d in S.items()}
grp = collections.defaultdict(list)
for k, d in S.items(): grp[("fs0" if abs(d["fs"]) < 0.5 else "fs23", k[1])].append(k)
F = [4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 200, 300, 500]
def T(keys, c, f):
    mk = (fr >= f / 1.12) & (fr < f * 1.12)
    num = sum(np.sum(sp[k][0][mk, c] * sp[k][1][mk, c].conj()) for k in keys); den = sum(np.sum(abs(sp[k][1][mk, c]) ** 2) for k in keys)
    return num / den
print("Pooled transfer ratio T(f) = <Y M*>/<|M|^2> (data over model; 1 = perfect), mag dB / phase deg, +-12% band, per channel and group")
print("group            ch | " + " ".join(f"{f:>9d}" for f in F))
for g in sorted(grp):
    for c in (0, 1):
        row = [T(grp[g], c, f) for f in F]
        print(f"{g[0]} E{g[1]:<5} n{len(grp[g]):2d} ch{c} | " + " ".join(f"{20*np.log10(abs(t)):5.1f}/{np.degrees(np.angle(t)):4.0f}" for t in row))
allk = list(S)
for c in (0, 1):
    row = [T(allk, c, f) for f in F]
    print(f"ALL n{len(allk)} ch{c} | " + " ".join(f"{20*np.log10(abs(t)):5.1f}/{np.degrees(np.angle(t)):4.0f}" for t in row))
