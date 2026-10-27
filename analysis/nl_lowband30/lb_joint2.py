"""Joint low-band transfer below ~100 Hz: table, agreement of T_N and T_S per frequency, and compact (rational) fits of the joint T on 15-100 Hz, scored on both programs
from the stored sums (S: other RT fold; N: leave-one-SPAN-out) as in lb_joint.py, own transfers used above the cut.   usage: lb_joint2.py"""
import sys
sys.argv = [sys.argv[0], "2", "1"]
exec(open("lb_joint.py").read().split("one = np.ones")[0])
one = np.ones((2, TNB), complex)
Sall = tot([S[q] for q in S_keys()]); Nall = tot(list(N.values()))
TJ = T_from([(*Sall, 1.0), (*Nall, 1.0)]); TNa = T_from([(*Nall, 1.0)]); TSa = T_from([(*Sall, 1.0)])
print("f(Hz): T_N, T_S, T_joint  (dB / deg), L then R")
for f in (16, 20, 24, 28, 32, 38, 44, 52, 60, 70, 80, 90, 100, 120):
    i = int(f / 2)
    print(f"{fc[i]:6.1f} | " + " | ".join(f"{20*np.log10(abs(T[c][i])):6.3f}/{np.degrees(np.angle(T[c][i])):6.2f}" for c in (0, 1) for T in (TNa, TSa, TJ)))
# ---- rational fits of T_J on 15-100 Hz
W0 = 2 * np.pi * 30.0
def fit_rat(f, T, n, iters=60):
    s = 1j * 2 * np.pi * f / W0; a = np.zeros(n)
    for it in range(iters):
        Dp = 1 + sum(a[k - 1] * s ** k for k in range(1, n + 1))
        cols = [s ** k / Dp for k in range(n)] + [-T * s ** k / Dp for k in range(1, n)] + [(s ** n - T * s ** n) / Dp]
        A = np.stack(cols, 1); rhs = T / Dp
        x = np.linalg.lstsq(np.vstack([A.real, A.imag]), np.r_[rhs.real, rhs.imag], rcond=None)[0]
        a = x[n:]; b = x[:n]
    return np.r_[b, a[-1]], np.r_[1.0, a]
def H(num, den, f):
    s = 1j * 2 * np.pi * np.asarray(f) / W0
    return sum(num[k] * s ** k for k in range(len(num))) / sum(den[k] * s ** k for k in range(len(den)))
LO, HI = 15.0, 100.0
mfit = (fc >= LO) & (fc <= HI)
taper = np.clip(np.minimum((fc - 12) / 3.0, (130 - fc) / 30.0), 0, 1)          # blend to 1 at 12 Hz and 130 Hz
def combine(Tlow, Town):
    """joint/compact transfer below ~100 Hz, own transfer above 130 Hz, linear blend between"""
    w = np.clip((130 - fc) / 30.0, 0, 1)
    return Tlow * w + Town * (1 - w)
def score(T_of_prog, prog_items, label):
    out = {}
    for bi, (lo, hi) in enumerate(BD[:5]):
        m = (fc >= lo) & (fc < hi); v = []
        for key, d, T in prog_items:
            for c in (0, 1): v.append(10 * np.log10(max(err(d, T, m, c), 1e-30) / d["syy"][c][m].sum()))
        out[bi] = np.median(v)
    return out
variants = {}
for n in (0, 6, 8, 10):                                 # 0 = the joint table itself (2 Hz bins)
    Tlow = {}
    for c in (0, 1):
        if n == 0: Tlow[c] = TJ[c]
        else:
            num, den = fit_rat(fc[mfit], TJ[c][mfit], n); Tlow[c] = np.where(fc < 12, 1.0, H(num, den, fc))
    variants[n] = Tlow
    # score: S (other fold own transfer above 130 Hz), N (LOSO own above)
    sI = []
    for k, v in S.items():
        other = "B" if v["rt"] in FOLDS["A"] else "A"
        Town = T_from([(*tot([S[q] for q in S_keys(other)]), 1.0)])
        Tj = np.stack([combine(Tlow[c] if n else Tlow[c], Town[c]) for c in (0, 1)]) if n == 0 else np.stack([combine(Tlow[c], Town[c]) for c in (0, 1)])
        sI.append((k, v, Tj))
    nI = []
    for sp in (0.2, 0.3, 0.4, 0.5):
        Town = T_from([(*tot([d for d in N.values() if d["span"] != sp]), 1.0)])
        for k, d in N.items():
            if d["span"] == sp: nI.append((k, d, np.stack([combine(Tlow[c], Town[c]) for c in (0, 1)])))
    rs, rn = score(None, sI, "S"), score(None, nI, "N")
    nm = "joint table" if n == 0 else f"rational order {n} ({4*n} numbers)"
    print(f"\n{nm}: below 100 Hz joint/compact, own transfer above 130 Hz")
    print("   band:      " + "  ".join(f"{lo}-{hi}" for lo, hi in BD[:5]))
    print("   S scored:  " + "  ".join(f"{rs[i]:6.1f}" for i in range(5)))
    print("   N scored:  " + "  ".join(f"{rn[i]:6.1f}" for i in range(5)))
