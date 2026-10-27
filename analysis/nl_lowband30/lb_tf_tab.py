import pickle, numpy as np
D = pickle.load(open("tf_all2.pkl", "rb")); F = D["F"]; rows = [r for r in D["rows"] if abs(r["fs"] - 26.2) > 1]
print(" f(Hz) |  L dB / deg  |  R dB / deg | n=%d" % len(rows))
T = [sum(r["num"][:, c] for r in rows) / sum(r["den"][:, c] for r in rows) for c in (0, 1)]
# bootstrap spread across captures
rng = np.random.default_rng(0); bs = []
for _ in range(200):
    idx = rng.integers(0, len(rows), len(rows)); rr = [rows[i] for i in idx]
    bs.append([sum(r["num"][:, c] for r in rr) / sum(r["den"][:, c] for r in rr) for c in (0, 1)])
bs = np.array(bs)
for i, f in enumerate(F):
    if f > 260: break
    sd = [np.std(np.abs(bs[:, c, i])) for c in (0, 1)]
    print(f"{f:6.1f} | {20*np.log10(abs(T[0][i])):6.2f} {np.degrees(np.angle(T[0][i])):6.1f} (sd {sd[0]*8.7:4.2f} dB) | {20*np.log10(abs(T[1][i])):6.2f} {np.degrees(np.angle(T[1][i])):6.1f} (sd {sd[1]*8.7:4.2f})")
