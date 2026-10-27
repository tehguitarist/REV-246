import pickle, collections, numpy as np
D = pickle.load(open("tf_all2.pkl", "rb")); F = D["F"]; rows = D["rows"]
sel = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33]
print("F:", np.round(F[sel], 1))
grp = collections.defaultdict(list)
for r in rows: grp[("fs0" if abs(r["fs"]) < .5 else "fs%d" % round(r["fs"]), r["key"][0])].append(r)
def pooled(rs, c):
    return sum(r["num"][:, c] for r in rs) / sum(r["den"][:, c] for r in rs)
for g in sorted(grp, key=lambda x: (x[1], x[0])):
    for c in (0, 1):
        T = pooled(grp[g], c)[sel]
        print(f"{g[0]:5s} SPAN {g[1]} n{len(grp[g]):2d} ch{c} " + " ".join(f"{20*np.log10(abs(t)):5.2f}/{np.degrees(np.angle(t)):4.0f}" for t in T))
for c in (0, 1):
    T = pooled(rows, c)[sel]; print(f"ALL n{len(rows)} ch{c} " + " ".join(f"{20*np.log10(abs(t)):5.2f}/{np.degrees(np.angle(t)):4.0f}" for t in T))
