"""Can S's per-loop correction be made compact? Cheap tests on the normal equations of `sml_fullband.py --stage xprobe`
(model_s_results.md §42.1). No model runs: every form is a linear least-squares solve on stored blocks.

    .venv/bin/python -I tools/xprobe_analyse.py [analysis/fullband/xprobe/groups.npz]

Blocks: per (RT fold A / B, LF) group, window (0-50 / 50-150 / 150-300 / 300-600 ms) and band (20-400 / 400-3200 /
3200-12300 Hz), the normal equations of the correction h (4 loops x lags, shared by L and R) against the baseline residual,
rows scaled per capture by its own energy in that band and window. The change a form h makes to a block is
10 log10((bTb - 2 h'ATb + h'ATA h) / bTb) dB (negative = better), a time-domain proxy of the bench's strict null with
the captures' clocks and gains NOT refitted. Training minimises the sum over the chosen (window, band) blocks of that
ratio (each block weighted 1 / bTb: a fractional error, so no band dominates by energy). Held out = trained on one RT
fold, evaluated on the other (both directions, averaged).
"""
import sys

import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "analysis/fullband/xprobe/groups.npz"
Z = np.load(path)
GR, LAGS, WINS, BANDS = list(Z["groups"]), list(Z["lags"]), [tuple(w) for w in Z["wins"]], [tuple(b) for b in Z["bands"]]
ATA, ATb, BTB, CNT = Z["ATA"], Z["ATb"], Z["BTB"], Z["CNT"]
NL, NLAG = 4, len(LAGS)
NW, NB = len(WINS), len(BANDS)
NCOL = NL * NLAG
col = lambda loop, lag: loop * NLAG + LAGS.index(lag)
FOLD = [g.split("|")[0] for g in GR]
LF = [float(g.split("|")[1]) for g in GR]
CH = [int(g.split("|")[2]) if g.count("|") == 2 else None for g in GR]      # present only in a XPROBE_CH=1 file
ECHO = [s * (e + d) for e in (230, 460) for d in range(-2, 3) for s in (1, -1)]
SHORT = [l for l in LAGS if abs(l) <= 10 or l in ECHO]


def sel(folds=None, lfs=None, chs=None):
    return [i for i in range(len(GR)) if (folds is None or FOLD[i] in folds) and (lfs is None or LF[i] in lfs) and (chs is None or CH[i] in chs)]


def blk(gs, w, b):
    return ATA[gs, w, b].sum(0), ATb[gs, w, b].sum(0), BTB[gs, w, b].sum()


def fit(gs, ws, bs, basis=None, lam=1e-6, whiten=True):
    M, v = 0, 0
    for w in ws:
        for b in bs:
            A, c, t = blk(gs, w, b)
            if t > 0:
                u = t if whiten else 1.0           # whiten: each block counts as its own relative change
                M, v = M + A / u, v + c / u
    if basis is not None:
        M, v = basis.T @ M @ basis, basis.T @ v
    th = np.linalg.solve(M + lam * np.trace(M) / len(M) * np.eye(len(M)), v)
    return th if basis is None else basis @ th


def change(h, gs, w, b):
    A, c, t = blk(gs, w, b)
    return 10 * np.log10(max((t - 2 * h @ c + h @ A @ h) / t, 1e-9)) if t > 0 else np.nan


def table(h_of, train, test, ws_eval=range(NW), label=""):
    """h_of(train_groups) -> h; evaluated on `test` groups, both directions if train/test are fold names."""
    out = np.full((NW, NB), np.nan)
    for tr, te in (("A", "B"), ("B", "A")):
        h = h_of(sel(folds=tr) if train == "fold" else train)
        gs = sel(folds=te) if test == "fold" else test
        for w in ws_eval:
            for b in range(NB):
                out[w, b] = np.nanmean([out[w, b], change(h, gs, w, b)])
    return out


def show(title, M):
    print(f"\n{title}")
    print("  window ms   " + "  ".join(f"{a}-{b} Hz".rjust(11) for a, b in BANDS) + "   worst cell")
    for w, (a, b) in enumerate(WINS):
        print(f"  {a:>3}-{b:<4}    " + "  ".join(f"{M[w, k]:+11.2f}" for k in range(NB)) + f"   {np.nanmax(M[w]):+.2f}")
    print(f"  mean {np.nanmean(M):+.2f} dB, cells worse than +0.1 dB: {int(np.sum(M > 0.1))} of {M.size}")


def basis_cols(cols):
    B = np.zeros((NCOL, len(cols)))
    for j, c in enumerate(cols):
        B[c, j] = 1
    return B


def basis_common(lags):
    B = np.zeros((NCOL, len(lags)))
    for j, l in enumerate(lags):
        for i in range(NL):
            B[col(i, l), j] = 1
    return B


TRAIN_W = [0, 1, 2]           # fitted where the first pass and first recirculation lie (as xfit); 300-600 is held out in time
ALLB = list(range(NB))
cols_short = [col(i, l) for i in range(NL) for l in SHORT]

print(f"{path}: {len(GR)} groups, {NCOL} columns ({NLAG} lags x {NL} loops)")

# 1. calibration against the run that was scored (164 numbers, windows 0-300, held out by RT fold)
show("1. calibration: the scored form (4 loops x short lags = 164 numbers), trained 0-300 ms on one fold, tested on the other",
     table(lambda g: fit(g, TRAIN_W, ALLB, basis_cols(cols_short)), "fold", "fold"))
print("   (the scored run's measured change, 50-150 ms: 80-400 Hz about -1 dB, 0.4-3.2 kHz -1.5, 3.2-6.4 kHz +1.8, 6.4-10 kHz -2.7)")

# 2. support: add the long lags
show("2. the same, all lags (95 per loop = 380 numbers)", table(lambda g: fit(g, TRAIN_W, ALLB), "fold", "fold"))
cols_long = [col(i, l) for i in range(NL) for l in LAGS if abs(l) > 10 and l not in ECHO]
show("2b. long lags only (12..90, step 3: no short lags, no echoes)", table(lambda g: fit(g, TRAIN_W, ALLB, basis_cols(cols_long)), "fold", "fold"))

# 3. in-sample ceiling: the best any static per-loop FIR of this support can do for ONE window and band, tested on that block
ceil = np.full((NW, NB), np.nan)
for w in range(NW):
    for b in range(NB):
        ceil[w, b] = np.mean([change(fit(sel(folds=tr), [w], [b]), sel(folds=te), w, b) for tr, te in (("A", "B"), ("B", "A"))])
show("3. ceiling: trained on ONE window and band only (all lags), tested on that block, other RT fold", ceil)

# 4. cross-window: train on one window (all bands, whitened), test everywhere
print("\n4. cross-window: trained on one window (all bands, all lags) -> mean change over bands in each test window (other fold)")
print("   train \\ test " + "  ".join(f"{a}-{b}".rjust(8) for a, b in WINS))
for w in range(NW):
    M = table(lambda g: fit(g, [w], ALLB), "fold", "fold")
    print(f"   {WINS[w][0]:>3}-{WINS[w][1]:<4}     " + "  ".join(f"{np.nanmean(M[v]):+8.2f}" for v in range(NW)))

# 5. cross-LF: train on LF >= 1.0 only, test on each LF
print("\n5. cross-LF: trained on LF >= 1.0 groups of one fold (windows 0-300, all lags) -> mean change over bands in each window, per test LF (other fold)")
for lf in sorted(set(LF)):
    out = np.full(NW, np.nan)
    for tr, te in (("A", "B"), ("B", "A")):
        h = fit(sel(folds=tr, lfs=[x for x in set(LF) if x >= 1.0]), TRAIN_W, ALLB)
        gs = sel(folds=te, lfs=[lf])
        if gs:
            out = np.nanmean([out, [np.nanmean([change(h, gs, w, b) for b in range(NB)]) for w in range(NW)]], axis=0)
    print(f"   LF {lf}: " + "  ".join(f"{WINS[w][0]}-{WINS[w][1]} ms {out[w]:+.2f}" for w in range(NW)))

# 6. compactness
print("\n6. compactness: trained 0-300 ms, other fold, mean dB over bands per window (0-50 / 50-150 / 150-300 / 300-600), worst cell")


def row(name, M, n):
    print(f"   {name:46s} {n:4d}  " + " / ".join(f"{np.nanmean(M[w]):+.2f}" for w in range(NW)) + f"   worst {np.nanmax(M):+.2f}  worse cells {int(np.sum(M > 0.1))}/{M.size}")


row("full 380", table(lambda g: fit(g, TRAIN_W, ALLB), "fold", "fold"), NCOL)
row("per-loop short lags (the scored form)", table(lambda g: fit(g, TRAIN_W, ALLB, basis_cols(cols_short)), "fold", "fold"), len(cols_short))
row("common to all loops, short lags", table(lambda g: fit(g, TRAIN_W, ALLB, basis_common(SHORT)), "fold", "fold"), len(SHORT))
row("common to all loops, all lags", table(lambda g: fit(g, TRAIN_W, ALLB, basis_common(LAGS)), "fold", "fold"), len(LAGS))
row("common, |lag| <= 3", table(lambda g: fit(g, TRAIN_W, ALLB, basis_common([l for l in LAGS if abs(l) <= 3])), "fold", "fold"), 7)
row("common, lag 0 only", table(lambda g: fit(g, TRAIN_W, ALLB, basis_common([0])), "fold", "fold"), 1)
row("per-loop lag 0 + 4 echo centres", table(lambda g: fit(g, TRAIN_W, ALLB, basis_cols([col(i, l) for i in range(NL) for l in (0, 230, -230, 460, -460)])), "fold", "fold"), 20)


def omp(gs, k_max, ws=TRAIN_W):
    """Forward selection of columns, exact refit, on the whitened training objective."""
    M, v = 0, 0
    for w in ws:
        for b in ALLB:
            A, c, t = blk(gs, w, b)
            M, v = M + A / t, v + c / t
    S, path_ = [], []
    for _ in range(k_max):
        best, bj = -1, None
        for j in range(NCOL):
            if j in S:
                continue
            T = S + [j]
            sol = np.linalg.solve(M[np.ix_(T, T)] + 1e-9 * np.trace(M) / NCOL * np.eye(len(T)), v[T])
            g_ = v[T] @ sol
            if g_ > best:
                best, bj = g_, j
        S.append(bj)
        path_.append(list(S))
    return path_, M, v


for k in (2, 4, 8, 16, 32):
    M_ = np.full((NW, NB), np.nan)
    for tr, te in (("A", "B"), ("B", "A")):
        paths, M, v = omp(sel(folds=tr), k)
        S = paths[-1]
        h = np.zeros(NCOL)
        h[S] = np.linalg.solve(M[np.ix_(S, S)] + 1e-9 * np.trace(M) / NCOL * np.eye(len(S)), v[S])
        for w in range(NW):
            for b in range(NB):
                M_[w, b] = np.nanmean([M_[w, b], change(h, sel(folds=te), w, b)])
    row(f"greedy sparse selection, {k} columns", M_, k)
    if k == 8:
        paths, M, v = omp(sel(folds="A"), 8)
        print("      fold-A picks (loop, lag): " + ", ".join(f"({c // NLAG},{LAGS[c % NLAG]})" for c in paths[-1]))


# 7. (a XPROBE_CH=1 file only) one correction shared by L and R against one per channel
if CH[0] is not None:
    print("\n7. shared by L and R vs one per channel (trained 0-300 ms on one fold, tested on the other, mean dB over bands per window)")
    forms = [("common to all loops, short lags", basis_common(SHORT)), ("common, |lag| <= 3", basis_common([l for l in LAGS if abs(l) <= 3])),
             ("common, lag 0 only", basis_common([0])), ("per-loop short lags", basis_cols(cols_short))]
    for name, Bm in forms:
        res = {}
        for mode in ("shared", "per channel"):
            M_ = np.full((NW, NB), np.nan)
            for tr, te in (("A", "B"), ("B", "A")):
                for c in (0, 1):
                    h = fit(sel(folds=tr) if mode == "shared" else sel(folds=tr, chs=[c]), TRAIN_W, ALLB, Bm)
                    for w in range(NW):
                        for b in range(NB):
                            M_[w, b] = np.nanmean([M_[w, b], change(h, sel(folds=te, chs=[c]), w, b)])
            res[mode] = M_
        for mode, M_ in res.items():
            print(f"   {name:34s} {mode:12s} " + " / ".join(f"{np.nanmean(M_[w]):+.2f}" for w in range(NW)) + f"   worst {np.nanmax(M_):+.2f}")
    hL = fit(sel(chs=[0]), TRAIN_W, ALLB, basis_common(SHORT))
    hR = fit(sel(chs=[1]), TRAIN_W, ALLB, basis_common(SHORT))
    print(f"   L and R common-short kernels (all data): cosine {hL @ hR / np.sqrt((hL @ hL) * (hR @ hR)):.3f}; lag-0 L {hL[SHORT.index(0)]:+.5f} R {hR[SHORT.index(0)]:+.5f}")


# 8. why the scored run fell short of section 1: xfit did not whiten (raw summed residual energy, dominated by the
#    bands with the largest fractional error). The same forms trained that way:
print("\n8. unwhitened training (blocks summed as they are: the low band's large fractional error dominates), as xfit did; other fold")
for name, Bm in (("per-loop short lags (164)", basis_cols(cols_short)), ("common to all loops, short lags (41)", basis_common(SHORT))):
    for wh in (False, True):
        M_ = table(lambda g: fit(g, TRAIN_W, ALLB, Bm, whiten=wh), "fold", "fold")
        show(f"8. {name}, {'whitened' if wh else 'UNWHITENED'}", M_)
