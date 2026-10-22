"""NEXT_STEPS §5 item 5: the S carry-over test for programs M and L (Tomkin 96 kHz library, device-rate IRs of
tools/ml_structure.py). Results: docs/ml_structure.md §7.

    .venv/bin/python tools/ml_carryover.py [--stage diffuser|taps|routing|firstpass|all] [--workers 10] [--only N]

Every stage works on cluster-free reverbs: per AMPL group, y_lo - A_lo / (A_hi - A_lo) (y_hi - y_lo) with the low-AMPL
capture aligned onto the high one (integer, then sub-sample on the reverb after the cluster), shifted so that cluster
tap 0 of the high-AMPL capture sits exactly at index K0. For L (A_lo = 0) that is the AMPL-0 capture itself.

The diffuser is all-pass, so correlating the first pass with the diffuser's impulse response (a matched filter) leaves
the read taps as isolated spikes when the diffuser is right, and smears them when it is not. That is the test.

Stages:
  diffuser   per program, the all-pass stage lengths that make the matched-filter output sparsest (share of energy in
             the 20 largest isolated peaks of the first pass), on a sample of settings: a local grid of +-3 samples per
             stage around the best set, then the common all-pass gain and a two-tap damping b on the last stage.
  taps       per program, every AMPL group: the 20 taps per channel (sub-sample position re cluster tap 0, sign,
             magnitude) by greedy pursuit over the first pass plus room for its last taps (band-limited pulses shaped by
             the program's output kernel, from the reference setting's cluster), matched to S's tap table (Q2,
             k1_sparsity.json) scaled by 2 / 4; S's taps scaled against M's and L's, and the R - L shift of each tap.
  routing    per program, each tap's recirculation through each line (gain and sub-sample delay re the tap at tap + P_j),
             flattest-HF settings: S's taps recirculate through their own loop only; M's and L's through two lines.
             The routing matrix (line <- the tap's S loop) is the coupling.
  firstpass  held out: cluster + diffuser + 20 equal sign-only taps per channel + the output kernel; positions, signs and
             kernel from 6 training settings (flattest HF, longest RT); per capture one shift, one cluster gain and one
             tap gain; per setting the last all-pass's two-tap damping b (as S's AP230 step), fitted on the high-AMPL
             capture. Nulls (100 Hz-10 kHz) in the cluster-and-gap window [0, onset) and the first pass [onset, onset +
             shortest line), for b = 0, b fitted, and b fitted without the kernel, against the repeat-capture null of the
             group with the same freedoms (another capture as the model).
"""
import argparse, json, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.signal import lfilter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
import ml_structure as ml  # noqa: E402
from core import subsample_peak  # noqa: E402

OUT = os.path.join(ROOT, "analysis", "ml_carryover")
CACHE = ml.CACHE
K0, FS = ml.K0, ml.FSD
SC = {"S": 1, "M": 2, "L": 4}
G_AP = -1 / np.sqrt(2)
DIFF = {"S": (20, 72, 230), "M": (42, 144, 460), "L": (82, 286, 918)}       # from the diffuser stage
LINES = {"S": [1321, 1389, 1461, 1587], "M": [2641, 2775, 2919, 3171], "L": [5279, 5547, 5837, 6339]}
REF = {"S": (1.0, 2.5, 1.0), "M": (1.2, 5.0, 0.9), "L": (1.2, 20.0, 1.0)}   # template setting (lf, rt, hf): long RT, flattest
CLUSTER_END = {"S": 353, "M": 704, "L": 1407}                                # last cluster tap (R), ml_structure.md §1
ONSET = {"S": (240, 572), "M": (474, 1133), "L": (948, 2268)}                 # reverb onset after tap 0 (L, R), §2
NTAP = 20


def fshift(x, s):
    """x delayed by s samples (band-limited FFT shift), along axis 0."""
    N = 1 << int(np.ceil(np.log2(len(x) + abs(s) + 256)))
    ph = np.exp(-2j * np.pi * np.fft.rfftfreq(N) * s)
    X = np.fft.rfft(x, N, axis=0)
    return np.fft.irfft(X * (ph[:, None] if x.ndim == 2 else ph), N, axis=0)[:len(x)]


def ap(x, M, g, b=0.0):
    """All-pass (-g + z^-M F) / (1 - g z^-M F), F = (1 - b) + b z^-1 (b = 0: the plain all-pass)."""
    F = np.r_[np.zeros(M), 1 - b, b]
    num = -g * np.r_[1.0, np.zeros(M + 1)] + F
    den = np.r_[1.0, np.zeros(M + 1)] - g * F
    return lfilter(num, den, x)


def diff_ir(Ms, n, g=G_AP, b=0.0):
    d = np.zeros(n)
    d[0] = 1.0
    for k, M in enumerate(Ms):
        d = ap(d, M, g, b if k == len(Ms) - 1 else 0.0)
    return d


# ---------------------------------------------------------------- data
_IR = ml.IRs(CACHE)


def groups(prog):
    return {f"{k[0]}_{k[1]}_{k[2]}": v for k, v in ml.ampl_groups(prog).items()}


def first_recirc(prog, ch=0):
    """End of the first pass in channel ch: the reverb onset plus the shortest line (R starts later than L)."""
    return ONSET[prog][ch] + min(LINES[prog])


def tap0(h):
    """Sub-sample position of cluster tap 0 (shared by L and R) near K0, from the R channel."""
    i = K0 - 4 + int(np.argmax(np.abs(h[K0 - 4:K0 + 5, 1])))
    return subsample_peak(h[:, 1], i)[0]


def sub_align(a, b, i0, i1):
    """Sub-sample shift s minimising |a - g fshift(b, s)| over [i0, i1) (both channels)."""
    from scipy.optimize import minimize_scalar

    def cost(s):
        x = fshift(b[:i1 + 64], s)[i0:i1]
        return np.sum((a[i0:i1] - np.sum(a[i0:i1] * x) / np.sum(x * x) * x) ** 2)
    return minimize_scalar(cost, bounds=(-1.5, 1.5), method="bounded", options={"xatol": 1e-3}).x


def captures(prog, key, n=None):
    """All captures of one AMPL group on a common grid with cluster tap 0 at K0 (band 100 Hz-10 kHz):
    {amp: y}, plus the cluster-free reverb r and the AMPL values."""
    g = groups(prog)[key]
    amps = sorted(g)
    hi = amps[-1]
    n = n or int((first_recirc(prog) + 3 * min(LINES[prog]) + 8000 * SC[prog]))
    pad = lambda y: np.vstack([y, np.zeros((max(0, n + 4000 * SC[prog] - len(y)), 2))])   # short files end early (faded)
    yh = pad(_IR(g[hi], band=True))
    t0 = tap0(yh)
    out = {hi: fshift(yh, K0 - t0)[:n]}
    i0, i1 = K0 + CLUSTER_END[prog] + 50, K0 + first_recirc(prog)
    for a in amps[:-1]:
        y = pad(_IR(g[a], band=True))
        _, y = ml.align(out[hi], y, i0, i1 + 2000 * SC[prog], maxlag=ONSET[prog][1] + 200)   # AMPL 0: onset = the reverb
        y = pad(y)
        s = sub_align(out[hi], y, i0, i1)
        out[a] = fshift(y, s)[:n]
    lo = amps[0]
    r = out[lo] - lo / (hi - lo) * (out[hi] - out[lo])
    return out, r, amps


def trace(r, Ms, n_out, g=G_AP, b=0.0):
    """Matched filter: r correlated with the diffuser's impulse response, per channel (index = tap position + K0)."""
    ext = 8 * max(Ms) * 10
    d = diff_ir(Ms, ext, g, b)
    N = 1 << int(np.ceil(np.log2(len(r) + ext)))
    return np.fft.irfft(np.fft.rfft(r, N, axis=0) * np.conj(np.fft.rfft(d, N))[:, None], N, axis=0)[K0:K0 + n_out]


def top_peaks(c, k, sep=3):
    e = c ** 2
    taken = []
    for i in np.argsort(e)[::-1]:
        if all(abs(i - j) > sep for j in taken):
            taken.append(int(i))
        if len(taken) == k:
            break
    return sorted(taken), float(sum(e[max(i - 1, 0):i + 2].sum() for i in taken) / e.sum())


def sample_keys(prog, n):
    """n settings spread over the library (every k-th key in label order), always including the reference."""
    ks = sorted(groups(prog))
    ref = "_".join(str(x) for x in REF[prog])
    pick = ks[::max(1, len(ks) // n)][:n]
    return [ref] + [k for k in pick if k != ref]


# ---------------------------------------------------------------- diffuser
def sparsity(r, Ms, prog, g=G_AP, b=0.0):
    c = trace(r, Ms, first_recirc(prog, 1) - 10, g, b)
    return [top_peaks(c[:first_recirc(prog, ch) - 10, ch], NTAP)[1] for ch in (0, 1)]


def _diffuser_job(a):
    prog, key, grid = a
    _, r, _ = captures(prog, key)
    import itertools
    best = DIFF[prog]
    res = dict(key=key, plain_x=sparsity(r, tuple(SC[prog] * m for m in DIFF["S"]), prog))
    if grid:
        sc = {}
        for Ms in itertools.product(*[range(m - 3, m + 4) for m in DIFF[prog]]):
            sc[Ms] = float(np.mean(sparsity(r, Ms, prog)))
        best = max(sc, key=sc.get)
        res["grid_best"] = list(best)
        res["grid_second"] = sorted(sc.values())[-2]
    res["best"] = sparsity(r, best, prog)
    res["gain"] = {f"{g:.3f}": float(np.mean(sparsity(r, best, prog, g))) for g in np.r_[np.arange(-0.80, -0.615, 0.02), G_AP]}
    res["b_last"] = {f"{b:.2f}": float(np.mean(sparsity(r, best, prog, G_AP, b))) for b in (0.0, 0.1, 0.17, 0.25, 0.33, 0.5)}
    return prog, res


def diffuser_stage(workers, n):
    jobs = [(p, k, i < 4) for p in "SML" for i, k in enumerate(sample_keys(p, n))]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_diffuser_job, jobs))
    out = {p: [r for q, r in res if q == p] for p in "SML"}
    json.dump(out, open(os.path.join(OUT, "diffuser.json"), "w"), indent=1)
    for p in "SML":
        rr = out[p]
        print(f"\n## {p}: diffuser {DIFF[p]} (S x {SC[p]} would be {tuple(SC[p] * m for m in DIFF['S'])}), {len(rr)} settings")
        grid = [r["grid_best"] for r in rr if "grid_best" in r]
        print(f"- grid (+-3 per stage) best on {len(grid)} settings: {grid}; runner-up mean sparsity "
              f"{np.mean([r['grid_second'] for r in rr if 'grid_second' in r]):.3f}")
        sp = np.array([r["best"] for r in rr])
        px = np.array([r["plain_x"] for r in rr])
        print(f"- top-{NTAP} energy share L / R: found {np.median(sp[:, 0]):.3f} / {np.median(sp[:, 1]):.3f} "
              f"(min {sp.min(0).round(3).tolist()}); plain x{SC[p]}: {np.median(px[:, 0]):.3f} / {np.median(px[:, 1]):.3f}")
        gk = list(rr[0]["gain"])
        print("- all-pass gain (median share): " + ", ".join(f"{k}: {np.median([r['gain'][k] for r in rr]):.3f}" for k in gk))
        bk = list(rr[0]["b_last"])
        print("- two-tap damping b on the last stage (median share): " + ", ".join(f"{k}: {np.median([r['b_last'][k] for r in rr]):.3f}" for k in bk))
        for r in rr:
            print(f"  {r['key']}: {np.round(r['best'], 3).tolist()} best gain {max(r['gain'], key=r['gain'].get)} best b {max(r['b_last'], key=r['b_last'].get)}")


# ---------------------------------------------------------------- taps
HP = 600
_PULSE = sosfiltfilt(ml.BAND, np.eye(1, 2 * HP + 1, HP)[0])      # one tap as the trace sees it (zero-phase band-limit)


def pursuit(c, k, refine=3, h=None):
    """Greedy pursuit of k taps in one channel of the trace: pick the largest residual sample, refit all amplitudes
    (least squares on band-limited pulses), repeat; then refine each position to sub-sample with the others removed.
    Overlapping taps (two within a few samples) and pulse sidelobes are handled. Returns positions, values and the
    share of the trace energy explained."""
    n = len(c)
    pulse = _PULSE if h is None else np.convolve(_PULSE, h, "same")     # the program's output kernel widens every tap

    def atom(p):
        a = np.zeros(n)
        i = int(np.floor(p))
        pl = fshift(np.r_[pulse, 0.0], p - i)[:2 * HP + 1]
        lo, hi = i - HP, i + HP + 1
        a[max(lo, 0):min(hi, n)] = pl[max(lo, 0) - lo:2 * HP + 1 - max(hi - n, 0)]
        return a
    pos = []
    res = c.copy()
    for _ in range(k):
        cand = np.argsort(np.abs(res))[::-1]
        pos.append(float(next(i for i in cand if all(abs(i - p) >= 2 for p in pos))))   # a fractional tap is one atom, not two
        A = np.stack([atom(p) for p in pos], 1)
        v = np.linalg.lstsq(A, c, rcond=None)[0]
        res = c - A @ v
    for _ in range(refine):
        for j in range(len(pos)):
            part = res + v[j] * atom(pos[j])
            i = int(round(pos[j]))
            i = i - 1 + int(np.argmax(np.abs(part[i - 1:i + 2])))
            pos[j] = subsample_peak(part, i)[0]
            A = np.stack([atom(p) for p in pos], 1)
            v = np.linalg.lstsq(A, c, rcond=None)[0]
            res = c - A @ v
    o = np.argsort(pos)
    return np.array(pos)[o], v[o], float(1 - np.sum(res ** 2) / np.sum(c ** 2))


EXTRA = 8                       # pursuit atoms beyond the 20 taps: early recirculations inside the extended window
OFFS = {"S": (0, 0), "M": (-9, -11), "L": (-27, -33)}   # rough offset of M / L taps from S x 2 / x 4 (L, R), to seed the search


def s_table():
    """S's 20 taps per channel (lag re cluster tap 0, sign, loop index), from model_s_results §16.2 (k1_sparsity.json)."""
    t = json.load(open(os.path.join(ROOT, "analysis", "q2", "k1_sparsity.json")))["tied"]
    L, R = [], []
    for i, Pi in enumerate(LINES["S"]):
        e = t[str(Pi)]
        for (lag, sg), sh in zip(e["L_taps"], e["R_shift"]):
            L.append((lag, sg, i))
            R.append((lag + sh, sg, i))
    return sorted(L), sorted(R)


def tap_template(prog):
    """Expected tap positions: S's table, scaled by 2 / 4 and offset (the search seed; measured positions replace it)."""
    return [[SC[prog] * lag + OFFS[prog][ch] for lag, _, _ in tab] for ch, tab in enumerate(s_table())]


def tap_window(prog, ch):
    """The first pass of channel ch plus room for its last tap (S's last taps come after the first recirculation)."""
    return first_recirc(prog, ch) + 260 * SC[prog]


def measure_taps(c, tmpl, prog, h=None):
    """Per channel: the NTAP taps of the pursuit over that channel's first pass, matched one-to-one to the template
    (within 1.5 samples; NaN when missing), their values, and the share of the first-pass trace energy explained."""
    pos = np.full((2, NTAP), np.nan)
    val = np.full((2, NTAP), np.nan)
    share = []
    for ch in (0, 1):
        p, v, sh = pursuit(c[:tap_window(prog, ch), ch], NTAP + EXTRA, h=h)
        used = set()
        for k, t in enumerate(tmpl[ch]):
            near = [j for j in range(len(p)) if j not in used and abs(p[j] - t) <= 2 + SC[prog]]
            if near:
                j = max(near, key=lambda j: abs(v[j]))             # the tap, not a neighbouring remainder atom
                pos[ch, k], val[ch, k] = p[j], v[j]
                used.add(j)
        share.append(sh)
    return pos, val, share


def _taps_job(a):
    prog, key, tmpl, h = a
    _, r, _ = captures(prog, key)
    c = trace(r, DIFF[prog], tap_window(prog, 1))
    pos, val, share = measure_taps(c, tmpl, prog, h)
    return prog, key, pos, val, share


def taps_stage(workers, only=None):
    tm = {p: tap_template(p) for p in "SML"}
    with ProcessPoolExecutor(workers) as ex:
        hk = {p: output_kernel(p, ["_".join(str(x) for x in REF[p])], ex)[0] for p in "SML"}   # atom shape: the reference cluster's kernel
        jobs = [(p, k, tm[p], hk[p]) for p in "SML" for k in (sorted(groups(p))[:only] if only else sorted(groups(p)))]
        res = list(ex.map(_taps_job, jobs))
    out = {}
    for p in "SML":
        rr = [x for x in res if x[0] == p]
        pos = np.array([x[2] for x in rr])           # (groups, ch, tap); NaN where the pursuit missed a template tap
        val = np.array([x[3] for x in rr])
        share = np.array([x[4] for x in rr])
        found = np.isfinite(pos).all(2).mean(0)
        mag = np.abs(val) / np.nanmedian(np.abs(val), axis=2, keepdims=True)
        sg = np.sign(np.nanmedian(val, 0))
        med = np.nanmedian(pos, 0)
        out[p] = dict(template=[list(map(int, t)) for t in tm[p]], groups=len(rr), all_found=found.tolist(),
                      pos_median=med.tolist(), pos_iqr=(np.nanpercentile(pos, 75, 0) - np.nanpercentile(pos, 25, 0)).tolist(),
                      pos_maxdev=np.nanmax(np.abs(pos - med), 0).tolist(), frac_median=(med - np.round(med)).tolist(),
                      sign=sg.astype(int).tolist(), sign_agree=np.nanmean(np.sign(val) == sg, 0).tolist(),
                      mag_db_spread=np.nanmedian(np.nanstd(20 * np.log10(mag), 2), 0).tolist(),
                      mag_db_median=np.nanmedian(20 * np.log10(mag), 0).tolist(),
                      share=np.median(share, 0).tolist(), share_p10=np.percentile(share, 10, 0).tolist(), keys=[x[1] for x in rr],
                      pos_all=np.round(pos, 3).tolist(), val_all=np.round(val, 5).tolist())
    # S scaled against M and L, and the R - L shifts (R tap paired with its L tap through S's 26 / 332 routing)
    S = np.array(out["S"]["pos_median"])
    for p in "ML":
        P = np.array(out[p]["pos_median"])
        out[p]["minus_scaled_S"] = [(P[ch] - SC[p] * S[ch]).tolist() for ch in (0, 1)]
    shifts = {}
    Ls, Rs = s_table()
    for p in "SML":
        P = np.array(out[p]["pos_median"])
        # each R tap is S's L tap of the same loop, routed through 26 or 332 (model_s_results §16.13): pair by S's table
        pairs = []
        for k, (lag, sg, i) in enumerate(Rs):
            j, s0 = next((j, lag - l) for j, (l, sl, il) in enumerate(Ls) if il == i and lag - l in (26, 332))
            pairs.append(dict(L=float(P[0, j]), R=float(P[1, k]), shift=float(P[1, k] - P[0, j]), s_route=int(s0), L_index=j, loop=i))
        shifts[p] = pairs
    out["shifts"] = shifts
    json.dump(out, open(os.path.join(OUT, "taps" + ("_smoke" if only else "") + ".json"), "w"), indent=1)
    for p in "SML":
        o = out[p]
        P = np.array(o["pos_median"])
        print(f"\n## {p}: {o['groups']} AMPL groups; taps re cluster tap 0 (median sub-sample position, sign)")
        for ch in (0, 1):
            print(f"- {'LR'[ch]}: " + " ".join(f"{x:.2f}{'+' if s > 0 else '-'}" for x, s in zip(P[ch], o["sign"][ch])))
            print(f"  max |dev| over groups {np.max(o['pos_maxdev'][ch]):.2f}, IQR max {np.max(o['pos_iqr'][ch]):.3f}; sign agreement min "
                  f"{np.min(o['sign_agree'][ch]):.3f}; magnitude spread over the 20 taps (std, median over groups) "
                  f"{o['mag_db_spread'][ch]:.2f} dB; all 20 found in {o['all_found'][ch] * 100:.0f} % of groups; "
                  f"energy explained median {o['share'][ch]:.3f} (p10 {o['share_p10'][ch]:.3f})")
        if p != "S":
            d = np.array(o["minus_scaled_S"])
            for ch in (0, 1):
                print(f"- {'LR'[ch]} minus {SC[p]} x S: " + " ".join(f"{x:+.2f}" for x in d[ch]))
        sh = shifts[p]
        print(f"- R - L shifts (route): " + " ".join(f"{x['shift']:.2f}({x['s_route']})" for x in sh))


# ---------------------------------------------------------------- routing
def _routing_job(a):
    prog, key, P = a
    _, r, _ = captures(prog, key)
    n = int(max(P[1]) + max(LINES[prog]) + 20)
    c = trace(r, DIFF[prog], n)
    g = np.zeros((2, NTAP, 4))
    d = np.zeros((2, NTAP, 4))
    v0 = np.zeros((2, NTAP))
    for ch in (0, 1):
        for k, t in enumerate(P[ch]):
            i = int(round(t))
            i = i - 2 + int(np.argmax(np.abs(c[i - 2:i + 3, ch])))
            t_, v0[ch, k] = subsample_peak(c[:, ch], i)
            for j, L in enumerate(LINES[prog]):
                i = int(round(t_ + L))
                i = i - 3 + int(np.argmax(np.abs(c[i - 3:i + 4, ch])))
                tj, vj = subsample_peak(c[:, ch], i)
                g[ch, k, j] = vj / v0[ch, k]
                d[ch, k, j] = tj - t_ - L
    return prog, key, g, d


def routing_stage(workers, only=None):
    tp = json.load(open(os.path.join(OUT, "taps.json")))
    jobs = []
    for p in "SML":
        ks = sorted(groups(p))
        ks = [k for k in ks if float(k.split("_")[2]) >= {"S": 1.0, "M": 0.9, "L": 1.0}[p]]   # flattest HF: clean delays
        jobs += [(p, k, tp[p]["pos_median"]) for k in (ks[:only] if only else ks)]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_routing_job, jobs))
    out = {}
    for p in "SML":
        rr = [x for x in res if x[0] == p]
        G = np.array([x[2] for x in rr])                         # (groups, ch, tap, line)
        D = np.array([x[3] for x in rr])
        Gm = np.median(G, 0)
        # per tap: the lines it feeds (|median gain| > 0.25 x the largest), their signs, and the delays
        feeds = []
        for ch in (0, 1):
            row = []
            for k in range(NTAP):
                a = np.abs(Gm[ch, k])
                js = [j for j in range(4) if a[j] > 0.25 * a.max()]
                row.append(dict(lines=js, sign=[int(np.sign(Gm[ch, k, j])) for j in js], gain=[float(Gm[ch, k, j]) for j in js],
                                sign_agree=[float(np.mean(np.sign(G[:, ch, k, j]) == np.sign(Gm[ch, k, j]))) for j in js],
                                delay=[float(np.median(D[:, ch, k, j])) for j in js]))
            feeds.append(row)
        # line delays: median over clean entries (one line fed with the expected magnitude)
        dl = []
        for j in range(4):
            v = [D[gi, ch, k, j] for ch in (0, 1) for k in range(NTAP) for gi in range(len(rr)) if j in feeds[ch][k]["lines"]]
            dl.append([float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75))])
        # routing matrix: line j <- S-loop i (each tap's loop from S's table), median gain re the tap over taps, channels, groups;
        # taps whose recirculation collides with another tap (|gain| > 1.5 x the median of the clean entries) are left out
        Ls, Rs = s_table()
        loop = [[i for _, _, i in Ls], [i for _, _, i in Rs]]
        F = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                v = np.concatenate([G[:, ch, k, j] for ch in (0, 1) for k in range(NTAP) if loop[ch][k] == i])
                F[j, i] = np.median(v)
        out[p] = dict(groups=len(rr), keys=[x[1] for x in rr], feeds=feeds, line_delay_offset=dl, gain_median=Gm.tolist(),
                      matrix=F.tolist(), rt=[float(k.split("_")[1]) for k in (x[1] for x in rr)])
    json.dump(out, open(os.path.join(OUT, "routing" + ("_smoke" if only else "") + ".json"), "w"), indent=1)
    for p in "SML":
        o = out[p]
        print(f"\n## {p}: routing over {o['groups']} groups; lines {LINES[p]}; delay of the recirculated tap re tap + line (median, IQR): "
              + ", ".join(f"{L}: {x[0]:+.2f} [{x[1]:+.2f}, {x[2]:+.2f}]" for L, x in zip(LINES[p], o["line_delay_offset"])))
        for ch in (0, 1):
            print(f"- {'LR'[ch]}: " + "; ".join(
                "{" + ",".join(f"{LINES[p][j]}{'+' if s > 0 else '-'}{abs(g):.2f}" for j, s, g in zip(f["lines"], f["sign"], f["gain"])) + "}"
                for f in o["feeds"][ch]))
            print(f"  sign agreement min {min(min(f['sign_agree']) for f in o['feeds'][ch]):.3f}")
        print("- routing matrix (rows: line; columns: the tap's S loop 1321 / 1389 / 1461 / 1587), median gain re the tap:")
        for j, row in enumerate(o["matrix"]):
            print(f"  {LINES[p][j]}: " + " ".join(f"{x:+.3f}" for x in row))


# ---------------------------------------------------------------- first pass, held out
KH = 5                          # output kernel half-length (taps -KH..KH)


def render(prog, P, sgn, cl, n, b=0.0, h=None):
    """Unit-gain model parts on the K0 grid: (cluster, first pass) per channel, positions possibly fractional;
    h: an output FIR (taps -KH..KH) applied to both."""
    N = 1 << int(np.ceil(np.log2(n + 8 * max(DIFF[prog]) * 10)))
    w = 2 * np.pi * np.fft.rfftfreq(N)
    D = np.fft.rfft(diff_ir(DIFF[prog], N, G_AP, b))
    C = np.zeros((len(w), 2), complex)
    T = np.zeros((len(w), 2), complex)
    for ch in (0, 1):
        for p_, s in cl[ch]:
            C[:, ch] += s * np.exp(-1j * w * (p_ + K0))
        for p_, s in zip(P[ch], sgn[ch]):
            T[:, ch] += s * np.exp(-1j * w * (p_ + K0))
    if h is not None:
        Hk = sum(hk * np.exp(-1j * w * (k - KH)) for k, hk in enumerate(h))
        C, T = C * Hk[:, None], T * Hk[:, None]
    band = lambda X: sosfiltfilt(ml.BAND, np.fft.irfft(X, N, axis=0)[:n], axis=0)
    return band(C), band(T * D[:, None])


def _kernel_job(a):
    """Least-squares output FIR (taps -KH..KH) on the cluster window of the high-AMPL capture, per channel."""
    prog, key = a
    caps, _, amps = captures(prog, key)
    y = caps[amps[-1]]
    n = K0 + ONSET[prog][1] + 50
    cl = render(prog, np.zeros((2, 0)), np.zeros((2, 0)), cluster_taps(prog), n + 50)[0]
    hs = []
    for ch in (0, 1):
        sl = slice(K0 - 30, K0 + ONSET[prog][ch] - 3 - KH)
        X = np.stack([np.roll(cl[:, ch], k)[sl] for k in range(-KH, KH + 1)], 1)
        hs.append(np.linalg.lstsq(X, y[sl, ch], rcond=None)[0])
    return prog, key, hs


def output_kernel(prog, keys, ex):
    """The program's output FIR from the training groups' clusters (median over captures and channels, unit DC)."""
    H = [h / h.sum() for _, _, hs in ex.map(_kernel_job, [(prog, k) for k in keys]) for h in hs]
    return np.median(H, 0), np.std(H, 0)


def cluster_taps(prog):
    s = json.load(open(os.path.join(ROOT, "analysis", "ml_structure", "summary.json")))["clusters"][prog]
    return [[(round(o), sg) for o, sg in zip(s["offsets"][c], s["signs"][c])] for c in "LR"]


def win(prog):
    """Per channel: (cluster and gap) [0, onset - 3), first pass [onset - 3, onset + shortest line - 10)."""
    return [[(0, ONSET[prog][c] - 3), (ONSET[prog][c] - 3, ONSET[prog][c] + min(LINES[prog]) - 10)] for c in (0, 1)]


def fit_null(y, parts, prog, amp_has_cluster):
    """One shift, cluster gain and tap gain per capture (the gains solved per channel jointly over both windows; from the
    first pass on only where the capture has no cluster); nulls per window and channel."""
    from scipy.optimize import minimize_scalar
    W = win(prog)
    i_end = K0 + W[1][1][1]

    def solve(s):
        ys = fshift(y[:i_end + 200], -s)
        res, mdl = 0.0, np.zeros_like(ys)
        for ch in (0, 1):
            sl = slice(K0 + W[ch][0 if amp_has_cluster else 1][0], K0 + W[ch][1][1])   # no cluster: fit the first pass only
            X = np.stack([p[sl, ch] for p in parts], 1)   # cluster column always: a partner capture's cluster can reach past the onset
            coef = np.linalg.lstsq(X, ys[sl, ch], rcond=None)[0]
            mdl[sl, ch] = X @ coef
            res += np.sum((ys[sl, ch] - mdl[sl, ch]) ** 2)
        return res, ys, mdl
    s = minimize_scalar(lambda s: solve(s)[0], bounds=(-1.0, 1.0), method="bounded", options={"xatol": 1e-3}).x
    _, ys, mdl = solve(s)
    out = {}
    for ch in (0, 1):
        for name, (a, b) in zip(("cluster_gap", "first_pass"), W[ch]):
            if name == "cluster_gap" and not amp_has_cluster:
                continue
            sl = slice(K0 + a, K0 + b)
            out.setdefault(name, [None, None])[ch] = float(10 * np.log10(np.sum((ys[sl, ch] - mdl[sl, ch]) ** 2) / np.sum(ys[sl, ch] ** 2)))
    return out, float(s)


def _firstpass_job(a):
    """One group: b (the last all-pass's two-tap damping, as S's AP230 step) fitted on the high-AMPL capture's first pass;
    then every capture nulled with b = 0 and with the fitted b, and against another capture of the group (repeat)."""
    from scipy.optimize import minimize_scalar
    prog, key, P, sgn, train, h = a
    caps, _, amps = captures(prog, key)
    n = K0 + win(prog)[1][1][1] + 400
    cl = cluster_taps(prog)
    hi = amps[-1]

    def cost(b):
        null, _ = fit_null(caps[hi][:n], render(prog, P, sgn, cl, n, b=b, h=h), prog, True)
        return sum(10 ** (x / 10) for x in null["first_pass"])
    b = float(minimize_scalar(cost, bounds=(0.0, 0.6), method="bounded", options={"xatol": 0.002}).x)
    parts = {"b0": render(prog, P, sgn, cl, n, b=0.0, h=h), "bfit": render(prog, P, sgn, cl, n, b=b, h=h),
             "bfit_nokernel": render(prog, P, sgn, cl, n, b=b)}
    rows = []
    for amp in amps:
        row = dict(key=key, amp=amp, train=train, b=b, b_fitted_on=hi)
        for v, pt in parts.items():
            row[v], row["shift_" + v] = fit_null(caps[amp][:n], pt, prog, amp > 0)
        # repeat capture as the model: another capture of the group (a different AMPL), same freedoms (shift, cluster gain, scale)
        other = [x for x in amps if x != amp][0]
        row["repeat"], _ = fit_null(caps[amp][:n], (parts["b0"][0], caps[other][:n]), prog, amp > 0)
        row["repeat_with"] = other
        rows.append(row)
    return prog, rows


def firstpass_stage(workers, only=None):
    tp = json.load(open(os.path.join(OUT, "taps.json")))
    jobs = []
    for p in "SML":
        ks = sorted(groups(p))
        if only:
            ks = ks[::max(1, len(ks) // only)][:only]
        # training settings: positions and signs from the taps stage restricted to these (the flattest HF, longest RTs)
        hf_max = max(float(k.split("_")[2]) for k in groups(p))
        tr = sorted([k for k in groups(p) if float(k.split("_")[2]) == hf_max], key=lambda k: -float(k.split("_")[1]))[:6]
        idx = [tp[p]["keys"].index(k) for k in tr if k in tp[p]["keys"]]
        jobs += [(p, k, None, None, k in tr) for k in ks]
        tp[p]["train"] = tr
        tp[p]["idx"] = idx
    # positions from the training groups only: re-measure there (taps.json holds medians over all groups)
    with ProcessPoolExecutor(workers) as ex:
        tm = {p: tap_template(p) for p in "SML"}
        Hk = {p: output_kernel(p, tp[p]["train"], ex) for p in "SML"}
        tr_jobs = [(p, k, tm[p], Hk[p][0]) for p in "SML" for k in tp[p]["train"]]
        tr_res = list(ex.map(_taps_job, tr_jobs))
        Pp = {p: np.median([x[2] for x in tr_res if x[0] == p], 0) for p in "SML"}
        Sg = {p: np.sign(np.median([x[3] for x in tr_res if x[0] == p], 0)) for p in "SML"}
        jobs = [(p, k, Pp[p], Sg[p], t, Hk[p][0]) for p, k, _, _, t in jobs]
        res = list(ex.map(_firstpass_job, jobs))
    out = {p: dict(train=tp[p]["train"], positions=Pp[p].tolist(), signs=Sg[p].astype(int).tolist(),
                   kernel=Hk[p][0].tolist(), kernel_std=Hk[p][1].tolist(),
                   rows=[r for q, rr in res if q == p for r in rr]) for p in "SML"}
    json.dump(out, open(os.path.join(OUT, "firstpass" + ("_smoke" if only else "") + ".json"), "w"), indent=1)
    for p in "SML":
        rr = out[p]["rows"]
        print(f"\n## {p}: first pass, held out ({len({r['key'] for r in rr})} groups, {len(rr)} captures; trained on {out[p]['train']}); "
              f"100 Hz-10 kHz, null L / R [repeat-capture null L / R]")
        print(f"- output kernel (taps -{KH}..{KH}, unit DC): " + " ".join(f"{x:+.4f}" for x in out[p]["kernel"])
              + f"; spread over training captures max {max(out[p]['kernel_std']):.4f}")
        for v in ("b0", "bfit", "bfit_nokernel"):
            for sel_name, sel in (("held out", lambda r: not r["train"]), ("training", lambda r: r["train"])):
                for w in ("cluster_gap", "first_pass"):
                    x = np.array([r[v][w] for r in rr if sel(r) and w in r[v]])
                    q = np.array([r["repeat"][w] for r in rr if sel(r) and w in r["repeat"]])
                    if len(x):
                        print(f"- {v}, {sel_name}, {w}: {np.median(x[:, 0]):.1f} / {np.median(x[:, 1]):.1f} "
                              f"[{np.median(q[:, 0]):.1f} / {np.median(q[:, 1]):.1f}]; p90 {np.percentile(x[:, 0], 90):.1f} / "
                              f"{np.percentile(x[:, 1], 90):.1f}; excess R {np.median(x[:, 1] - q[:, 1]):+.1f}; n {len(x)}")
        hfs = sorted({float(r["key"].split("_")[2]) for r in rr})
        rts = sorted({float(r["key"].split("_")[1]) for r in rr})
        bt = {}
        for r in rr:
            bt.setdefault((float(r["key"].split("_")[1]), float(r["key"].split("_")[2])), []).append(r["b"])
        print("- fitted b (median over LF), rows RT, columns HF " + " ".join(f"{h:>5}" for h in hfs))
        for rt in rts:
            print(f"  {rt:>5}: " + " ".join(f"{np.median(bt[(rt, h)]):5.3f}" if (rt, h) in bt else "    -" for h in hfs))
        for hf in hfs:
            x = np.array([r["bfit"]["first_pass"] for r in rr if not r["train"] and float(r["key"].split("_")[2]) == hf])
            q = np.array([r["repeat"]["first_pass"] for r in rr if not r["train"] and float(r["key"].split("_")[2]) == hf])
            if len(x):
                print(f"  HF {hf}: first pass (b fitted) {np.median(x[:, 0]):.1f} / {np.median(x[:, 1]):.1f} [{np.median(q[:, 0]):.1f} / {np.median(q[:, 1]):.1f}]")


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--stage", default="all")
    ap_.add_argument("--workers", type=int, default=10)
    ap_.add_argument("--only", type=int, default=None, help="smoke test: this many groups per program")
    ap_.add_argument("--n", type=int, default=12, help="diffuser stage: settings per program")
    a = ap_.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.stage in ("diffuser", "all"):
        diffuser_stage(a.workers, a.only or a.n)
    if a.stage in ("taps", "all"):
        taps_stage(a.workers, a.only)
    if a.stage in ("routing", "all"):
        routing_stage(a.workers, a.only)
    if a.stage in ("firstpass", "all"):
        firstpass_stage(a.workers, a.only)


if __name__ == "__main__":
    main()
