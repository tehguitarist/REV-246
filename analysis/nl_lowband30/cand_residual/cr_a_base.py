"""Step 0: model null (bench alignment fixed, and gain+shift refitted in band) vs the repeat-capture level, identical captures / window / filter.
usage: cr_a_base.py [limit]"""
from cr_common import *
import itertools, collections
B = bench_all()
lim = int(sys.argv[1]) if len(sys.argv) > 1 else 10**9
EAR = (0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 10.0)
sel = [k for k in sorted(B) if k[1] in EAR and not B[k]["trimmed"] and (k[0] == 0.5 or (k[0] in (0.2, 0.3, 0.4) and k[2] == 0.0))][:lim]
BD = [(20, 100), (100, 1000), (1000, 5000), (5000, 10000)]
SOS = {b: band_sos(*b) for b in BD}
f = np.fft.rfftfreq(NCACHE) * 2 * np.pi
def deriv(x):
    return np.fft.irfft(np.fft.rfft(x, axis=0) * (-1j * f)[:, None], NCACHE, axis=0)      # d/dn of x(n) (a delay of s samples is -s d/dn)
res = {}
bandsig = {}
t = time.time()
for k in sel:
    r = B[k]; y, m = cached(k, r); a, b = window(k, r)
    d = {}
    for bd in BD:
        yb = sosfiltfilt(SOS[bd], y, axis=0)[a:b]; mb = sosfiltfilt(SOS[bd], m, axis=0)[a:b]
        dmb = sosfiltfilt(SOS[bd], deriv(m), axis=0)[a:b]
        dyb = sosfiltfilt(SOS[bd], deriv(y), axis=0)[a:b]
        rows = []
        for c in (0, 1):
            ey = fpow(yb[:, c])
            fixed = fpow(yb[:, c] - mb[:, c]) / ey
            X = np.c_[mb[:, c], dmb[:, c]]; co = np.linalg.lstsq(X, yb[:, c], rcond=None)[0]
            refit = fpow(yb[:, c] - X @ co) / ey
            rows.append((fixed, refit))
        d[bd] = rows
        bandsig[(k, bd)] = (yb, dyb)
    res[k] = d
print("modelled", len(sel), time.time() - t, flush=True)
# repeat pairs
rng = np.random.default_rng(0)
rep = {}
for sp in (0.5,):
    for e in EAR:
        ks = [k for k in sel if k[0] == sp and k[1] == e]
        pairs = [(p, q) for p, q in itertools.combinations(ks, 2) if abs(B[p]["front_shift"] - B[q]["front_shift"]) < 0.3 and abs(B[p]["o"] - B[q]["o"]) < 0.02]
        if len(pairs) > 20: pairs = [pairs[i] for i in sorted(rng.choice(len(pairs), 20, replace=False))]
        out = {bd: [] for bd in BD}
        for p, q in pairs:
            for bd in BD:
                ya, _ = bandsig[(p, bd)]; yb_, dyb_ = bandsig[(q, bd)]
                r_ = []
                for c in (0, 1):
                    X = np.c_[yb_[:, c], dyb_[:, c]]; co = np.linalg.lstsq(X, ya[:, c], rcond=None)[0]
                    r_.append(fpow(ya[:, c] - X @ co) / fpow(ya[:, c]))
                out[bd].append(r_)
        rep[e] = (len(pairs), out)
def db(x): return 10 * np.log10(x)
print("SPAN 0.5, window 0.05-0.45 s. median over captures of L/R-mean dB. model fixed | model refit(g,shift) | repeat pair raw (g,shift fitted)")
tot = collections.defaultdict(list)
for e in EAR:
    ks = [k for k in sel if k[0] == 0.5 and k[1] == e]
    line = f"EARLY {e:5.1f} n={len(ks):2d} pairs={rep[e][0]:2d}"
    for bd in BD:
        fx = np.median([0.5 * (db(res[k][bd][0][0]) + db(res[k][bd][1][0])) for k in ks])
        rf = np.median([0.5 * (db(res[k][bd][0][1]) + db(res[k][bd][1][1])) for k in ks])
        rp = np.median([0.5 * (db(v[0]) + db(v[1])) for v in rep[e][1][bd]]) if rep[e][1][bd] else float("nan")
        line += f" | {bd[0]}-{bd[1]}: {fx:6.1f} {rf:6.1f} {rp:6.1f}"
        for k in ks: tot[(bd, "fx")].append(0.5 * (db(res[k][bd][0][0]) + db(res[k][bd][1][0])))
        for k in ks: tot[(bd, "rf")].append(0.5 * (db(res[k][bd][0][1]) + db(res[k][bd][1][1])))
        for v in rep[e][1][bd]: tot[(bd, "rp")].append(0.5 * (db(v[0]) + db(v[1])))
    print(line)
print("pooled:", " | ".join(f"{bd[0]}-{bd[1]}: " + " ".join(f"{np.median(tot[(bd, w)]):6.1f}" for w in ("fx", "rf", "rp")) for bd in BD))
for sp in (0.2, 0.3, 0.4):
    ks = [k for k in sel if k[0] == sp]
    print("SPAN", sp, len(ks), "model fixed (window to last tap+20ms):", " | ".join(f"{bd[0]}-{bd[1]}: {np.median([0.5*(db(res[k][bd][0][0])+db(res[k][bd][1][0])) for k in ks]):6.1f}" for bd in BD))
json.dump({tag(k): {f"{bd[0]}-{bd[1]}": v for bd, v in d.items()} for k, d in res.items()}, open(os.path.join(HERE, "a_base.json"), "w"))
