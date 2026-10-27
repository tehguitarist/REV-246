import pickle, numpy as np, json
exec(open("lb_tfcorr.py").read().split("def interp")[0])
fs0 = [k for k, d in S.items() if abs(d["fs"]) < .5]; fs23 = [k for k in S if k not in fs0]
Fs = np.geomspace(6, 400, 28)
def est2(keys, c, f):
    mk = (fr >= f / 1.12) & (fr < f * 1.12)
    return sum(np.sum(sp[k][0][mk, c] * sp[k][1][mk, c].conj()) for k in keys) / sum(np.sum(abs(sp[k][1][mk, c]) ** 2) for k in keys)
print("f(Hz) | fs23 L, R | fs0 L, R   : T = data/model, dB / deg")
res = {}
for f in Fs:
    row = []
    for g in (fs23, fs0):
        for c in (0, 1):
            t = est2(g, c, f); row.append(t)
    res[float(f)] = [[float(t.real), float(t.imag)] for t in row]
    print(f"{f:6.1f} | " + " | ".join(f"{20*np.log10(abs(t)):6.2f}/{np.degrees(np.angle(t)):6.1f}" for t in row))
json.dump(res, open("tf_pooled.json", "w"))
