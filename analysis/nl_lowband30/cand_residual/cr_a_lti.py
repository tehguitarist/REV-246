"""(a) Is the residual still deterministic, and how much of it is a fixed LTI transfer T(f) (per channel, pooled)?
 1. coherent part: cross-capture covariance of e = y - m between captures of one (SPAN, EARLY) that differ only in LATE.
 2. pooled per-group T on the corrected model, held out by EARLY / by SPAN / by session (front shift)."""
from cr_common import *
import itertools, collections
B = bench_all()
EAR = (0.7, 1.0, 1.5, 2.0, 3.0, 10.0)
sel = [k for k in sorted(B) if k[1] in EAR and not B[k]["trimmed"] and (k[0] == 0.5 or (k[0] in (0.2, 0.3, 0.4) and k[2] == 0.0))]
NF = 16384
fr = np.fft.rfftfreq(NF, 1 / FS)
BD = {"20-100": (20, 100), "100-1000": (100, 1000), "1000-5000": (1000, 5000)}
def db(x): return 10 * np.log10(x)
S = {}
E = {}
PRE = butter(6, [12, 1600], btype='band', fs=FS, output='sos')
for k in sel:
    r = B[k]; y, m = cached(k, r); a, b = window(k, r)
    yf_, mf_ = sosfiltfilt(PRE, y, axis=0), sosfiltfilt(PRE, m, axis=0)      # band-limit in time first: rect-window leakage of out-of-band error would otherwise dominate the LF bins
    S[k] = (np.fft.rfft(yf_[a:b], NF, axis=0), np.fft.rfft(mf_[a:b], NF, axis=0), b - a)
    sos = {bd: band_sos(*v) for bd, v in BD.items()}
    E[k] = {bd: (sosfiltfilt(s, y, axis=0)[a:b], sosfiltfilt(s, y - m, axis=0)[a:b]) for bd, s in sos.items()}
print("loaded", len(sel), flush=True)
# ---- 1. coherent part of e across captures of one setting (differing only in LATE)
print("\n[1] cross-capture covariance of the residual (captures of one (SPAN, EARLY) differing only in LATE, |dfs|<0.3, |do|<0.02)")
print("   coherent = <e_A, e_B> / |y|^2 (mean over pairs, dB; negative cov -> n/a), rho = <e_A,e_B>/sqrt(|e_A|^2|e_B|^2); total = |e|^2/|y|^2")
rng = np.random.default_rng(1)
for bd in BD:
    print(" band", bd)
    for sp in (0.5,):
        for e_ in EAR:
            ks = [k for k in sel if k[0] == sp and k[1] == e_]
            pairs = [(p, q) for p, q in itertools.combinations(ks, 2) if abs(B[p]["front_shift"] - B[q]["front_shift"]) < 0.3 and abs(B[p]["o"] - B[q]["o"]) < 0.02]
            if not pairs: continue
            coh, rho, tot = [], [], []
            for p, q in pairs:
                for c in (0, 1):
                    yp, ep = E[p][bd]; yq, eq = E[q][bd]
                    coh.append(ep[:, c] @ eq[:, c] / (0.5 * (fpow(yp[:, c]) + fpow(yq[:, c]))))
                    rho.append(ep[:, c] @ eq[:, c] / np.sqrt(fpow(ep[:, c]) * fpow(eq[:, c])))
                    tot.append(0.5 * (fpow(ep[:, c]) / fpow(yp[:, c]) + fpow(eq[:, c]) / fpow(yq[:, c])))
            print(f"   EARLY {e_:5.1f} {len(pairs):3d} pairs: total {db(np.median(tot)):6.1f}  coherent {db(max(np.mean(coh),1e-12)):6.1f}  incoherent {db(max(np.mean(tot)-np.mean(coh),1e-12)):6.1f}  median rho {np.median(rho):5.2f}")
# ---- 2. pooled T by group, held out
def groups(width_hz=None, sixth=False, lo=15.0, hi=1500.0):
    idx = np.nonzero((fr >= lo) & (fr < hi))[0]
    if sixth:
        edges = lo * 2 ** (np.arange(0, 200) / 6.0); edges = edges[edges < hi]
        g = np.digitize(fr[idx], edges)
    else:
        g = ((fr[idx] - lo) // width_hz).astype(int)
    return idx, g
def estimate(train, idx, g):
    ng = g.max() + 1
    T = np.zeros((2, ng), complex)
    for c in (0, 1):
        num = np.zeros(ng, complex); den = np.zeros(ng)
        for k in train:
            Y, M, _ = S[k]
            w = 1.0 / np.sum(np.abs(Y[idx, c]) ** 2)          # each capture counts equally (unit signal energy)
            np.add.at(num, g, w * Y[idx, c] * np.conj(M[idx, c])); np.add.at(den, g, w * np.abs(M[idx, c]) ** 2)
        T[c] = num / np.maximum(den, 1e-30)
    return T
def score(test, T, idx, g, bandmask):
    out = {bd: [[], []] for bd in BD}
    for k in test:
        Y, M, _ = S[k]
        for c in (0, 1):
            for bd, (lo, hi) in BD.items():
                sel_ = (fr >= lo) & (fr < hi)
                den = np.sum(np.abs(Y[sel_, c]) ** 2)
                e0 = np.sum(np.abs(Y[sel_, c] - M[sel_, c]) ** 2) / den
                Mt = M[:, c].copy(); Mt[idx] = Mt[idx] * T[c][g]
                e1 = np.sum(np.abs(Y[sel_, c] - Mt[sel_]) ** 2) / den
                out[bd][0].append(e0); out[bd][1].append(e1)
    return {bd: (db(np.median(v[0])), db(np.median(v[1]))) for bd, v in out.items()}
def show(name, res):
    print(f"   {name:34s}" + " | ".join(f"{bd}: {a:6.1f} -> {b:6.1f} ({b-a:+5.1f})" for bd, (a, b) in res.items()))
for label, kw in (("5 Hz groups (3 bins)", dict(width_hz=4.83)), ("20 Hz groups", dict(width_hz=19.3)), ("1/6 octave", dict(sixth=True))):
    idx, g = groups(**kw); print(f"\n[2] T per group, {label}: {g.max()+1} groups x 2 (complex) x 2 ch = {4*(g.max()+1)} reals up to 1.5 kHz; spectral null (rect window), median over test captures of channel-nulls")
    # leave-one-EARLY-out (SPAN 0.5), train on all other EARLY of SPAN 0.5
    all5 = [k for k in sel if k[0] == 0.5]
    resf = {bd: [[], []] for bd in BD}
    for e_ in EAR:
        tr = [k for k in all5 if k[1] != e_]; te = [k for k in all5 if k[1] == e_]
        T = estimate(tr, idx, g); r_ = score(te, T, idx, g, None)
        show(f"leave EARLY {e_} out (n={len(te)})", r_)
    # leave-SPAN-out: train SPAN 0.5 -> test 0.2/0.3/0.4 LATE 0; train 0.2-0.4 -> test 0.5
    T = estimate(all5, idx, g)
    for sp in (0.2, 0.3, 0.4):
        te = [k for k in sel if k[0] == sp]; show(f"train SPAN 0.5 -> test SPAN {sp} (n={len(te)})", score(te, T, idx, g, None))
    T = estimate([k for k in sel if k[0] < 0.5], idx, g)
    show("train SPAN 0.2-0.4 -> test SPAN 0.5", score(all5, T, idx, g, None))
    # session split (front shift) at SPAN 0.5
    s0 = [k for k in all5 if B[k]["front_shift"] < 10]; s1 = [k for k in all5 if B[k]["front_shift"] >= 10]
    show(f"train fs>=10 ({len(s1)}) -> test fs<10 ({len(s0)})", score(s0, estimate(s1, idx, g), idx, g, None))
    show(f"train fs<10 -> test fs>=10", score(s1, estimate(s0, idx, g), idx, g, None))
    # in-sample for reference
    show("in-sample (all SPAN 0.5)", score(all5, estimate(all5, idx, g), idx, g, None))
    if kw.get("width_hz") == 4.83:
        np.savez(os.path.join(HERE, "a_T5hz.npz"), f=fr[idx][np.r_[True, np.diff(g) > 0]], T=estimate(all5, idx, g))
