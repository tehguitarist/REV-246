exec(open("lb_tfcorr.py").read().split("def interp")[0])
F2 = [6, 8, 10, 12, 15, 18, 21, 25, 30, 40, 60]
def T2(keys, c, f):
    mk = (fr >= f / 1.12) & (fr < f * 1.12)
    return sum(np.sum(sp[k][0][mk, c] * sp[k][1][mk, c].conj()) for k in keys) / sum(np.sum(abs(sp[k][1][mk, c]) ** 2) for k in keys)
print("T(f) dB/deg by LATE (pooled over the 33 captures, L)  f:", F2)
for late in sorted({k[2] for k in S}):
    ks = [k for k in S if k[2] == late]
    print(f"LATE {late:3.1f} n{len(ks):2d} " + " ".join(f"{20*np.log10(abs(t)):5.2f}/{np.degrees(np.angle(t)):4.0f}" for t in (T2(ks, 0, f) for f in F2)))
