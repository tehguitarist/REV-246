"""(d) model error (bench alignment) split by channel, SPAN and session; window 0.05-0.45 s or to last tap + 20 ms."""
from cr_common import *
import collections
B = bench_all()
EAR = (1.0, 1.5, 2.0, 3.0, 10.0)
BD = [(20, 100), (100, 400), (400, 1000), (1000, 5000)]
SOS = {b: band_sos(*b) for b in BD}
def db(x): return 10 * np.log10(x)
acc = collections.defaultdict(list)
for k in sorted(B):
    if k[1] not in EAR or B[k]["trimmed"] or not (k[0] == 0.5 or (k[0] in (0.2, 0.3, 0.4) and k[2] == 0.0)): continue
    y, m = cached(k, B[k]); a, b = window(k, B[k]); ses = round(B[k]["front_shift"] / 10)
    for bd in BD:
        for c in (0, 1):
            yb = sosfiltfilt(SOS[bd], y[:, c])[a:b]; mb = sosfiltfilt(SOS[bd], m[:, c])[a:b]
            v = db(fpow(yb - mb) / fpow(yb))
            acc[("ch", 'LR'[c], bd)].append(v); acc[("span", k[0], bd)].append(v); acc[("ses", ses, bd)].append(v)
            if k[1] == 2.0: acc[("span@E2", k[0], bd)].append(v)
for grp in ("ch", "span", "span@E2", "ses"):
    print(grp)
    for key in sorted({kk[1] for kk in acc if kk[0] == grp}, key=str):
        print(f"  {key}: " + " | ".join(f"{bd[0]}-{bd[1]}: {np.median(acc[(grp, key, bd)]):6.1f} (n={len(acc[(grp, key, bd)])})" for bd in BD))
