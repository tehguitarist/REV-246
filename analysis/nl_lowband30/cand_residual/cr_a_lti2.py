"""(a2) Does an LTI T(f) fitted on the SAME EARLY (other captures) remove the coherent error, where one fitted on other EARLY does not?
Hann-tapered segment, series pre-band-limited 12-3000 Hz. Groups: 1/6 octave (41 groups up to 1.6 kHz; 2 reals each, per channel)."""
from cr_common import *
B = bench_all()
EAR = (1.0, 1.5, 2.0, 3.0, 10.0)
sel = [k for k in sorted(B) if k[0] == 0.5 and k[1] in EAR and not B[k]["trimmed"]]
NF = 16384
fr = np.fft.rfftfreq(NF, 1 / FS)
PRE = butter(6, [12, 3000], btype='band', fs=FS, output="sos")
BD = {"20-100": (20, 100), "100-400": (100, 400), "400-1000": (400, 1000), "100-1000": (100, 1000)}
def db(x): return 10 * np.log10(x)
S = {}
for k in sel:
    r = B[k]; y, m = cached(k, r); a, b = window(k, r)
    w = np.hanning(b - a)[:, None]
    S[k] = (np.fft.rfft(sosfiltfilt(PRE, y, axis=0)[a:b] * w, NF, axis=0), np.fft.rfft(sosfiltfilt(PRE, m, axis=0)[a:b] * w, NF, axis=0))
edges = 16 * 2 ** (np.arange(0, 60) / 6.0); edges = edges[edges <= 1600]
idx = np.nonzero((fr >= 16) & (fr < 1600))[0]; g = np.digitize(fr[idx], edges) - 1
ng = g.max() + 1
def estimate(train):
    T = np.zeros((2, ng), complex)
    for c in (0, 1):
        num = np.zeros(ng, complex); den = np.zeros(ng)
        for k in train:
            Y, M = S[k]; w = 1.0 / np.sum(np.abs(Y[idx, c]) ** 2)
            np.add.at(num, g, w * Y[idx, c] * np.conj(M[idx, c])); np.add.at(den, g, w * np.abs(M[idx, c]) ** 2)
        T[c] = num / np.maximum(den, 1e-30)
    return T
def score(te, T):
    o = {bd: [[], []] for bd in BD}
    for k in te:
        Y, M = S[k]
        for c in (0, 1):
            Mt = M[:, c].copy(); Mt[idx] *= T[c][g]
            for bd, (lo, hi) in BD.items():
                s_ = (fr >= lo) & (fr < hi); den = np.sum(np.abs(Y[s_, c]) ** 2)
                o[bd][0].append(np.sum(np.abs(Y[s_, c] - M[s_, c]) ** 2) / den); o[bd][1].append(np.sum(np.abs(Y[s_, c] - Mt[s_]) ** 2) / den)
    return {bd: (db(np.median(v[0])), db(np.median(v[1]))) for bd, v in o.items()}
def show(n, r): print(f"  {n:46s}" + " | ".join(f"{bd}: {a:6.1f}->{b:6.1f}" for bd, (a, b) in r.items()))
print("EARLY | model error dB -> after a 1/6-octave complex T estimated as stated (median over test captures; 2 ch)")
for e_ in EAR:
    ks = [k for k in sel if k[1] == e_]
    s0 = [k for k in ks if B[k]["front_shift"] < 10]; s1 = [k for k in ks if B[k]["front_shift"] >= 10]
    print(f"EARLY {e_}: {len(s0)} captures fs<10, {len(s1)} fs>=10")
    others = [k for k in sel if k[1] != e_]
    show("T from other EARLY (all)", score(ks, estimate(others)))
    if s0 and s1:
        show("T from same EARLY, other session", score(s0, estimate(s1))); show("   (reverse)", score(s1, estimate(s0)))
    # leave-one-capture-out inside same EARLY same session
    r_ = {bd: [[], []] for bd in BD}
    for grp in (s0, s1):
        for k in grp:
            tr = [q for q in grp if q != k]
            if len(tr) < 2: continue
            rr = score([k], estimate(tr))
            for bd in BD: r_[bd][0].append(rr[bd][0]); r_[bd][1].append(rr[bd][1])
    if r_["20-100"][0]: show("T from same EARLY+session, other LATEs", {bd: (np.median(v[0]), np.median(v[1])) for bd, v in r_.items()})
    show("per-capture oracle (T fitted on the capture)", {bd: tuple(np.median([score([k], estimate([k]))[bd][i] for k in ks]) for i in (0, 1)) for bd in BD})
