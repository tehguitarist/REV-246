"""Grey-box ("Path B", docs/EMT246_discovery.md §5) model of EMT 246 program S, fitted to and
null-tested against a device-rate IR.

    .venv/bin/python tools/model_s.py [--h analysis/shortRT/h_dev.npy] [--out analysis/model_s]

Writes <out>/summary.json, <out>/model_s.png, model_B_dev.npy / model_C_dev.npy / measured_dev.npy
(arrays (n, 2) at fs_dev, index PRE = first cluster tap) and model_B_params.npz. Results: docs/model_s_results.md.

Model, per channel c (time origin n = 0 at the first cluster tap):

    y_c = B * cluster_c  +  f_c  +  sum_i  x_ic * r_i          (per-channel read-out), or
    y_c = B * cluster_c  +  f_c  +  d_c * sum_i x_i * r_i      (common loop sum + output distributor d_c)

  B          band-limit kernel of the analysis chain (tools/deconv_grid.py window, 11-12.3 kHz roll-off)
  cluster_c  4 fixed integer taps with fitted gains
  r_i        loop i impulse response 1 / (1 - G_i(w) e^{-jwP_i}),  P = 1321, 1389, 1461, 1587
  G_i(w)     per-pass response: gain * first-order shelf * fractional delay (fitted; see fit_loops_*)
  x_ic       per-loop excitation, P_i samples starting at S0 (linear least squares)
  f_c        optional first-pass FIR on [S0, S0 + LF) (non-recirculating part; linear LS)

Loop filters are identified without the excitations from the annihilator residual
e = prod_i (1 - G_i z^-P_i) y, which is zero (to noise) once every loop's excitation has ended.
"""
import argparse, json, os, sys
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.sparse.linalg import LinearOperator, lsqr
from scipy.signal import butter, sosfilt

sys.path.insert(0, os.path.dirname(__file__))
from analyze_ir import edc, rt_from_edc, echo_density, OCTAVES  # noqa: E402

FS = 19.2e6 / 728
P = np.array([1321, 1389, 1461, 1587])
K0 = 263                       # first cluster tap index in h_dev
PRE = 96                       # samples kept before the first tap (kernel pre-ripple)
CLUSTER = {0: [0, 54, 194, 288], 1: [0, 143, 194, 353]}
CLUSTER_SIGN = {0: [1, 1, 1, -1], 1: [1, 1, -1, 1]}
S0 = 230                       # excitation start (samples after first tap); L reverb onset ~240, R ~570
NFFT = 1 << 17                 # 4.97 s: comb tails at ~ -300 dB before wrap-around
W = 2 * np.pi * np.fft.rfftfreq(NFFT)          # rad/sample
F = np.fft.rfftfreq(NFFT, 1 / FS)


# ---------- data ----------
def load(path, n_end_s=1.9):
    h = np.load(path)
    return h[K0 - PRE:K0 + int(n_end_s * FS)].copy()       # index i  <->  n = i - PRE


def bl_kernel(half=96):
    """Device-rate image of a single-sample tap after deconv_grid's band-limit window."""
    f = np.linspace(0, FS / 2, 40001); df = f[1] - f[0]
    win = np.where(f < 11000, 1.0, np.where(f < 12300, 0.5 * (1 + np.cos(np.pi * (f - 11000) / 1300)), 0.0))
    n = np.arange(-half, half + 1)
    return n, (2 * np.cos(2 * np.pi * np.outer(n, f) / FS) @ win - win[0]) * df / FS


def cluster_kernels(eps=0.0, kernel=None, half=PRE):
    """Per-tap kernels: the band-limit kernel (analytic, or a measured one on n=-half..half)
    shifted to the tap's sub-sample position t*(1+eps) (capture/device clock ratio)."""
    n = np.arange(-half, half + 1)
    if kernel is None:
        kernel = bl_kernel(half)[1]
    M = 1 << 12; K = np.fft.rfft(np.roll(np.pad(kernel, (0, M - len(kernel))), -half)); w = 2 * np.pi * np.fft.rfftfreq(M)
    kb = {c: [np.roll(np.fft.irfft(K * np.exp(-1j * w * t * eps), M), half)[:len(n)] for t in CLUSTER[c]] for c in (0, 1)}
    return n, kb


def fit_kernel(y, half=PRE, eps=0.0, rounds=4):
    """Measured tap kernel: LS fit of one common kernel to the reverb-free cluster regions
    (L: n < 230, R: n < 560), alternating with the 8 tap gains."""
    regions = {0: 230, 1: 560}; n = np.arange(-half, half + 1)
    gains = {c: np.array(CLUSTER_SIGN[c], float) for c in (0, 1)}
    kern = bl_kernel(half)[1]
    for _ in range(rounds):
        rows, rhs = [], []
        for c in (0, 1):
            seg = np.arange(-half, regions[c]); A = np.zeros((len(seg), len(n)))
            for g, t in zip(gains[c], CLUSTER[c]):
                for j, k in enumerate(n):
                    m = seg == t + k
                    A[m, j] += g
            rows.append(A); rhs.append(y[PRE + seg, c])
        kern = np.linalg.lstsq(np.vstack(rows), np.concatenate(rhs), rcond=None)[0]
        kern /= kern[half]
        for c in (0, 1):
            _, kb = cluster_kernels(eps, kern, half)
            seg = np.arange(-half, regions[c]); A = np.zeros((len(seg), 4))
            for j, t in enumerate(CLUSTER[c]):
                i = t + n; m = (i < regions[c]); A[i[m] + half, j] = kb[c][j][m]
            gains[c] = np.linalg.lstsq(A, y[PRE + seg, c], rcond=None)[0]
    return kern, gains


def cluster_subsample(y):
    """Sub-sample position of every cluster tap (fit of the shifted band-limit kernel)."""
    f = np.linspace(0, FS / 2, 20001); df = f[1] - f[0]
    win = np.where(f < 11000, 1.0, np.where(f < 12300, 0.5 * (1 + np.cos(np.pi * (f - 11000) / 1300)), 0.0))
    kern = lambda t: (2 * np.cos(2 * np.pi * np.outer(t, f) / FS) @ win) * df / FS
    out = {}
    for c in (0, 1):
        res = []
        for t in CLUSTER[c]:
            n = np.arange(-30, 31); seg = y[PRE + t + n, c]
            def cost(d):
                b = kern(n - d); return np.sum((seg - b * (b @ seg) / (b @ b)) ** 2)
            res.append(float(minimize_scalar(cost, bounds=(-0.3, 0.3), method="bounded", options={"xatol": 1e-6}).x))
        out[c] = res
    lags = np.concatenate([CLUSTER[0][1:], CLUSTER[1][1:]]); offs = np.concatenate([out[0][1:], out[1][1:]])
    eps = float(lags @ offs / (lags @ lags))
    return out, eps


# ---------- loop per-pass responses ----------
def shelf(w, r0, rpi, p):
    """First-order section (b0 + b1 z^-1)/(1 - p z^-1) with DC gain r0, Nyquist gain rpi."""
    b0 = (r0 * (1 - p) + rpi * (1 + p)) / 2; b1 = (r0 * (1 - p) - rpi * (1 + p)) / 2
    z1 = np.exp(-1j * w)
    return (b0 + b1 * z1) / (1 - p * z1)


def jot_gain(Pi, rt):
    return 10 ** (-3 * Pi / (FS * rt))


def G_param(theta, kind, w=W):
    """Per-pass responses for the 4 loops from a small parameter vector.
    kind: 'jot'      RT0, RTpi, pole p (shared), Jot-scaled per loop;  tau_i = tau0 + eps*P_i
          'onepole_jot'     one-pole lowpass per loop with Jot-scaled DC/Nyquist gains (p unused)
          'onepole_shared'  one fixed one-pole lowpass, same p in all loops, Jot DC gain (rtpi unused)
          'shared'   same shelf coefficients for every loop (gain still Jot-scaled at DC)
    theta = [rt0, rtpi, p, tau0, eps_ppm]"""
    rt0, rtpi, p, tau0, eps = theta[:5]; eps = eps * 1e-6
    Gs = []
    for Pi in P:
        if kind == "jot":
            H = shelf(w, jot_gain(Pi, rt0), jot_gain(Pi, rtpi), p)
        elif kind == "onepole_jot":   # pure one-pole lowpass per loop, pole set by Jot DC/Nyquist gains (theta p unused)
            r0, rpi = jot_gain(Pi, rt0), jot_gain(Pi, rtpi); pi_ = (r0 - rpi) / (r0 + rpi)
            H = r0 * (1 - pi_) / (1 - pi_ * np.exp(-1j * w))
        elif kind == "onepole_shared":  # fixed one-pole (same p in every loop) after a Jot DC gain (theta rtpi unused)
            H = jot_gain(Pi, rt0) * (1 - p) / (1 - p * np.exp(-1j * w))
        else:  # shared: per-loop DC gain from RT, HF/DC ratio fixed by the mean loop
            Pm = P.mean(); H = jot_gain(Pi, rt0) * shelf(w, 1.0, jot_gain(Pm, rtpi) / jot_gain(Pm, rt0), p)
        Gs.append(H * np.exp(-1j * w * (tau0 + eps * Pi)))
    return Gs


def G_fir(c_taps, K, w=W):
    k = np.arange(-K, K + 1)
    return [np.exp(-1j * np.outer(w, k)) @ ci for ci in c_taps]


def loop_ir(G, Pi):
    return np.fft.irfft(1.0 / (1.0 - G * np.exp(-1j * W * Pi)), NFFT)


# ---------- annihilator identification ----------
def _fwd(y):
    return np.fft.rfft(y, NFFT, axis=0)


def annihilate(Y, Gs, skip=None):
    """prod_{i != skip} (1 - G_i e^{-jwP_i}) applied to spectra Y (NFFT/2+1, ch)."""
    A = np.ones(len(W), complex)
    for i, (G, Pi) in enumerate(zip(Gs, P)):
        if i != skip:
            A = A * (1 - G * np.exp(-1j * W * Pi))
    return np.fft.irfft(Y * A[:, None], NFFT, axis=0)


def tail_weight(n, rt=0.98, t_cap=0.75):
    t = np.minimum(n / FS, t_cap)
    return 10 ** (3 * t / rt)


def fit_loops_fir(y, K=6, win=(0.25, 0.75), sweeps=12, rt=0.98, ridge=1e-6):
    """Nonparametric per-pass FIR G_i (taps k=-K..K around the P_i delay) by alternating LS
    on the annihilator residual over win (s after the first tap); both channels share G_i."""
    Y = _fwd(y); n = np.arange(NFFT) - PRE
    sel = (n >= win[0] * FS) & (n < win[1] * FS); wt = tail_weight(n[sel], rt)[:, None]
    k = np.arange(-K, K + 1)
    taps = [np.where(k == 0, jot_gain(Pi, rt), 0.0) for Pi in P]
    hist = []
    for s in range(sweeps):
        for i in range(4):
            Gs = G_fir(taps, K)
            u = annihilate(Y, Gs, skip=i)                                   # (NFFT, ch)
            idx = np.nonzero(sel)[0]
            X = np.stack([u[idx - P[i] - kk, :] for kk in k], -1)  # (len, ch, 2K+1)
            A = (X * wt[:, :, None]).reshape(-1, len(k)); b = (u[idx] * wt).reshape(-1)
            taps[i] = np.linalg.solve(A.T @ A + ridge * np.trace(A.T @ A) * np.eye(len(k)), A.T @ b)
        e = annihilate(Y, G_fir(taps, K))[sel] * wt
        hist.append(float(np.sum(e ** 2)))
    return taps, hist


def annihilator_cost(y, Gs, win=(0.25, 0.75), rt=0.98):
    Y = _fwd(y); n = np.arange(NFFT) - PRE
    sel = (n >= win[0] * FS) & (n < win[1] * FS)
    e = annihilate(Y, Gs)[sel] * tail_weight(n[sel], rt)[:, None]
    return e


def fit_loops_param(y, kind="jot", theta0=(0.98, 0.85, 0.0, 0.05, 100.0), fix=(), win=(0.25, 0.75)):
    """Parametric per-pass responses by nonlinear LS on the annihilator residual.
    fix: names among ('tau0','eps','p','rtpi') held at theta0."""
    names = ["rt0", "rtpi", "p", "tau0", "eps"]
    free = [j for j, nm in enumerate(names) if nm not in fix]
    th = np.array(theta0, float)
    Y = _fwd(y); n = np.arange(NFFT) - PRE
    sel = (n >= win[0] * FS) & (n < win[1] * FS); wt = tail_weight(n[sel])[:, None]
    def res(v):
        t = th.copy(); t[free] = v
        return (annihilate(Y, G_param(t, kind))[sel] * wt).ravel()
    lo = np.array([0.5, 0.2, -0.95, -2.0, -2000.0])[free]; hi = np.array([2.0, 2.0, 0.95, 2.0, 2000.0])[free]
    r = least_squares(res, th[free], bounds=(lo, hi), x_scale=np.array([0.01, 0.01, 0.01, 0.01, 10.0])[free], diff_step=1e-4)
    th[free] = r.x
    return th, float(np.sum(r.fun ** 2))


# ---------- linear stage: excitations, first-pass FIR, cluster gains ----------
class Linear:
    """Operator for one or two channels. mode 'per_channel' (x_ic) or 'distributor'
    (common x_i, channel FIR d_c over lags 0..KD applied to the loop sum)."""

    def __init__(self, Gs, N, LF=0, mode="per_channel", d=None, weight=None, s0=S0, ext=0, eps=0.0, kernel=None,
                 Ps=None, post=None, post_f=False):
        self.N, self.LF, self.mode, self.s0 = N, LF, mode, s0
        Ps = list(P) if Ps is None else list(Ps)
        self.Lx = [int(Pi) + ext for Pi in Ps]
        self.R = [1.0 / (1.0 - G * np.exp(-1j * W * Pi)) for G, Pi in zip(Gs, Ps)]   # comb spectra
        self.post = np.ones(len(W)) if post is None else post   # series element after the loop sum
        self.R = [R * self.post for R in self.R]
        self.Fpost = self.post if post_f else None
        self.kn, self.kb = cluster_kernels(eps, kernel)
        self.d = d                                   # (KD+1, 2) or None
        self.w = np.ones(N) if weight is None else weight
        self.nx = int(sum(self.Lx))
        self.nch_x = 2 if mode == "per_channel" else 1
        self.cs = 100.0      # cluster columns scaled so min-norm LSQR assigns the taps to the cluster, not to f/x
        if d is not None:
            self.D = np.fft.rfft(d, NFFT, axis=0)

    def size(self):
        return 2 * (4 + self.LF) + self.nch_x * self.nx

    def gains(self, v):
        return v[:8].reshape(2, 4) * self.cs

    def split(self, v):
        o = 0; cl = v[o:o + 8].reshape(2, 4); o += 8
        f = v[o:o + 2 * self.LF].reshape(2, self.LF); o += 2 * self.LF
        xs = []
        for ch in range(self.nch_x):
            row = []
            for L in self.Lx:
                row.append(v[o:o + L]); o += L
            xs.append(row)
        return cl, f, xs

    def loop_sum(self, xs_ch):
        S = np.zeros(NFFT // 2 + 1, complex)
        for x, R in zip(xs_ch, self.R):
            buf = np.zeros(NFFT); buf[PRE + self.s0:PRE + self.s0 + len(x)] = x
            S += np.fft.rfft(buf) * R
        return S

    def forward(self, v):
        cl, f, xs = self.split(v); out = np.zeros((self.N, 2))
        for c in (0, 1):
            for g, t, kb in zip(cl[c], CLUSTER[c], self.kb[c]):
                i = PRE + t + self.kn; m = (i >= 0) & (i < self.N); out[i[m], c] += self.cs * g * kb[m]
            if self.Fpost is None:
                out[PRE + S0:PRE + S0 + self.LF, c] += f[c]
            else:
                buf = np.zeros(NFFT); buf[PRE + S0:PRE + S0 + self.LF] = f[c]
                out[:, c] += np.fft.irfft(np.fft.rfft(buf) * self.Fpost, NFFT)[:self.N]
        if self.mode == "per_channel":
            for c in (0, 1):
                out[:, c] += np.fft.irfft(self.loop_sum(xs[c]), NFFT)[:self.N]
        else:
            S = self.loop_sum(xs[0])
            for c in (0, 1):
                out[:, c] += np.fft.irfft(S * self.D[:, c], NFFT)[:self.N]
        return (out * self.w[:, None]).T.ravel()

    def adjoint(self, r):
        r = r.reshape(2, self.N).T * self.w[:, None]
        Rf = np.fft.rfft(r, NFFT, axis=0)
        cl = np.zeros((2, 4)); f = np.zeros((2, self.LF)); xs = []
        for c in (0, 1):
            for j, t in enumerate(CLUSTER[c]):
                i = PRE + t + self.kn; m = (i >= 0) & (i < self.N); cl[c, j] = self.cs * (r[i[m], c] @ self.kb[c][j][m])
            if self.Fpost is None:
                f[c] = r[PRE + S0:PRE + S0 + self.LF, c]
            else:
                f[c] = np.fft.irfft(Rf[:, c] * np.conj(self.Fpost), NFFT)[PRE + S0:PRE + S0 + self.LF]
        if self.mode == "per_channel":
            for c in (0, 1):
                xs.append([np.fft.irfft(Rf[:, c] * np.conj(R), NFFT)[PRE + self.s0:PRE + self.s0 + L] for R, L in zip(self.R, self.Lx)])
        else:
            Z = sum(Rf[:, c] * np.conj(self.D[:, c]) for c in (0, 1))
            xs.append([np.fft.irfft(Z * np.conj(R), NFFT)[PRE + self.s0:PRE + self.s0 + L] for R, L in zip(self.R, self.Lx)])
        return np.concatenate([cl.ravel(), f.ravel()] + [x for row in xs for x in row])

    def op(self):
        n = self.size()
        return LinearOperator((2 * self.N, n), matvec=self.forward, rmatvec=self.adjoint, dtype=float)

    def solve(self, y, iters=400, x0=None):
        b = (y[:self.N] * self.w[:, None]).T.ravel()
        r = lsqr(self.op(), b, iter_lim=iters, atol=1e-12, btol=1e-12, x0=x0)
        return r[0]

    def render(self, v, w_one=True):
        w = self.w; self.w = np.ones(self.N)
        out = self.forward(v).reshape(2, self.N).T; self.w = w
        return out


def fit_distributor(y, Gs, N, weight, KD=700, rounds=6, iters=300, ridge=1e-4, init=None, **kw):
    """Common loop sum + per-channel output FIR d_c (lags 0..KD), by alternating LS:
    excitations x_i (shared by L and R) given d_c, then d_c given the loop sum s."""
    if init is None:
        d = np.zeros((KD + 1, 2)); d[0, 0] = 1.0; d[26, 1] = 0.5; d[332, 1] = 0.5   # T1 lags (assumed start)
    else:
        d = init.copy()
    hist = []
    for it in range(rounds):
        lin = Linear(Gs, N, mode="distributor", d=d, weight=weight, **kw)
        v = lin.solve(y, iters=iters)
        cl, f, xs = lin.split(v)
        s = np.fft.irfft(lin.loop_sum(xs[0]), NFFT)[:N]
        base = lin.render(np.concatenate([cl.ravel() / lin.cs, f.ravel(), np.zeros(lin.nx)]))
        Xs = np.stack([np.concatenate([np.zeros(j), s[:N - j]]) for j in range(KD + 1)], 1) * weight[:, None]
        AtA = Xs.T @ Xs; AtA += ridge * np.trace(AtA) / len(AtA) * np.eye(len(AtA))
        for c in (0, 1):
            d[:, c] = np.linalg.solve(AtA, Xs.T @ ((y[:N, c] - base[:, c]) * weight))
        lin = Linear(Gs, N, mode="distributor", d=d, weight=weight, **kw)
        v = lin.solve(y, iters=iters)
        r = lin.forward(v) - (y[:N] * weight[:, None]).T.ravel()
        hist.append(float(r @ r))
    return lin, v, d, hist


def fit_extra_loop(y, theta, P5, win=(0.25, 0.75), kind="jot"):
    """Add one more recursion 1 - G5 z^-P5 (gain g5, first-order HF ratio) to the annihilator and
    refit everything. Returns theta', g5, hf5, cost (4 loops), cost (5)."""
    Y = _fwd(y); n = np.arange(NFFT) - PRE
    sel = (n >= win[0] * FS) & (n < win[1] * FS); wt = tail_weight(n[sel])[:, None]
    def res(v):
        A5 = 1 - G5_of(v[:5], v[5], v[6], P5) * np.exp(-1j * W * P5)
        E = Y * A5[:, None]
        return (annihilate(E, G_param(v[:5], kind))[sel] * wt).ravel()
    v0 = np.r_[theta, 0.5, 1.0]; c4 = float(np.sum(res(np.r_[theta, 0.0, 1.0]) ** 2))
    r = least_squares(res, v0, x_scale=np.array([.01, .01, .01, .01, 10, .05, .05]), diff_step=1e-4,
                      bounds=([.5, .2, -.95, -2, -2000, -.99, 0.0], [2, 2, .95, 2, 2000, .99, 1.5]))
    return r.x[:5], float(r.x[5]), float(r.x[6]), c4, float(np.sum(r.fun ** 2))


def G5_of(theta, g5, hf5, P5):
    return g5 * shelf(W, 1.0, hf5, 0.0) * np.exp(-1j * W * theta[4] * 1e-6 * P5)


def pd_excitations(y, Gs, theta, t_fit=(0.1, 0.8)):
    """Spec-literal baseline: one period per loop from least-squares periodic decomposition of the
    decay-flattened (RT0) tail, mapped back to an excitation starting at S0."""
    import scipy.sparse as sps
    eps = theta[4] * 1e-6; n = np.arange(len(y)) - PRE
    m = (n >= t_fit[0] * FS) & (n < t_fit[1] * FS); nn = n[m]
    xs = []
    for c in (0, 1):
        yf = y[m, c] * 10 ** (3 * nn / (FS * theta[0]))
        rows, cols, off = [], [], 0
        for Pi in P:
            ph = np.round(nn / (1 + eps)).astype(int) % Pi       # position within the period on the device clock
            rows.append(np.arange(len(nn))); cols.append(off + ph); off += Pi
        A = sps.csr_matrix((np.ones(len(nn) * 4), (np.concatenate(rows), np.concatenate(cols))), shape=(len(nn), off))
        sol = lsqr(A, yf, iter_lim=3000, atol=1e-10, btol=1e-10)[0]
        o = 0; row = []
        for Pi in P:
            w = sol[o:o + Pi]; o += Pi
            j = np.arange(Pi); g = jot_gain(Pi, theta[0])
            row.append(w[(S0 + j) % Pi] * g ** ((S0 + j) / Pi))
        xs.append(row)
    return xs


# ---------- metrics ----------
WINDOWS_MS = [(0, 50), (50, 300), (300, 1000)]


def null_depth(meas, model, windows=WINDOWS_MS):
    """10 log10(sum (meas - model)^2 / sum meas^2) per window, time from the first cluster tap."""
    out = {}
    for a, b in windows:
        i0, i1 = PRE + int(a / 1000 * FS), PRE + int(b / 1000 * FS)
        out[f"{a}-{b}"] = [float(10 * np.log10(np.sum((meas[i0:i1, c] - model[i0:i1, c]) ** 2) / np.sum(meas[i0:i1, c] ** 2))) for c in (0, 1)]
    return out


def t30_octaves(y, noise_tail=None):
    """T30 per octave at the device rate (8 kHz band capped at 12.3 kHz)."""
    res = []
    for c in range(y.shape[1]):
        row = {}
        for fc in OCTAVES:
            sos = butter(4, [fc / np.sqrt(2), min(fc * np.sqrt(2), 12300)], btype="band", fs=FS, output="sos")
            s = sosfilt(sos, y[:, c])
            d, _ = edc(s, FS)
            row[str(fc)] = rt_from_edc(d, FS, -5, -35)
        res.append(row)
    return res


def ned(y):
    return [echo_density(y[:, c], FS) for c in range(y.shape[1])]


def lr_ratio_test(xs, theta):
    """Common-distributor check: if L = d_L * s and R = d_R * s for every loop, the harmonic ratio
    X_iL/X_iR equals D_L/D_R(f) for all loops, so close harmonics of different loops share it.
    Phase score: mean |e^{j dphi} - 1|^2 / 2 (0 = equal phase, ~1 = unrelated)."""
    eps = theta[4] * 1e-6; fr, rho, lp = [], [], []
    for i, Pi in enumerate(P):
        g = jot_gain(Pi, theta[0]); X = []
        for c in (0, 1):
            x = xs[c][i]; w = x[:Pi].copy(); w[:len(x) - Pi] += x[Pi:] / g
            X.append(np.fft.rfft(np.roll(w, S0)))
        fm = np.arange(len(X[0])) * FS / (Pi * (1 + eps)); m = (fm > 100) & (fm < 10000)
        fr.append(fm[m]); rho.append(X[0][m] / X[1][m]); lp.append(np.full(m.sum(), i))
    fr, rho, lp = (np.concatenate(v) for v in (fr, rho, lp)); o = np.argsort(fr); fr, rho, lp = fr[o], rho[o], lp[o]
    def score(a, b):
        ph = np.abs(np.exp(1j * (np.angle(rho[a]) - np.angle(rho[b]))) - 1) ** 2 / 2
        mag = np.abs(20 * np.log10(np.abs(rho[a]) / np.abs(rho[b])))
        return {"phase_score": float(np.mean(ph)), "median_abs_db_diff": float(np.median(mag)), "n": int(len(a))}
    out = {}
    for tol in (1.0, 5.0):
        a, b = [], []
        for j in range(len(fr)):
            k = j + 1
            while k < len(fr) and fr[k] - fr[j] < tol:
                if lp[k] != lp[j]:
                    a.append(j); b.append(k)
                k += 1
        out[f"cross_loop_within_{tol:g}Hz"] = score(np.array(a), np.array(b))
    rng = np.random.default_rng(0); out["random_pairs"] = score(rng.integers(0, len(fr), 5000), rng.integers(0, len(fr), 5000))
    return out


def residual_mode_test(meas, model, theta, win=(0.1, 0.6), noise_win=(1.4, 1.9)):
    """Is the residual at loop-mode frequencies or broadband?  On the decay-flattened window, per channel:
      on_off_db        spectral density within 0.3 Hz of a loop harmonic over density > 1.2 Hz from every
                       harmonic (0 dB = broadband; the measured IR itself gives the 'fully modal' value)
      mode_shaped_frac share of residual energy that follows the measured IR's modal fine structure:
                       in each 250 Hz band, S_res = a*S_meas + b (a, b >= 0), share = sum a*S_meas / sum S_res
    Also reported for the measurement's own noise tail (noise_win, unflattened) as a broadband control."""
    from scipy.optimize import nnls
    eps = theta[4] * 1e-6; n = np.arange(len(meas)) - PRE
    m = (n >= win[0] * FS) & (n < win[1] * FS); nn = n[m]; fl = 10 ** (3 * nn / (FS * theta[0]))
    Nf = 1 << 21; f = np.fft.rfftfreq(Nf, 1 / FS); band = (f > 100) & (f < 11000)
    harm = np.sort(np.concatenate([np.arange(1, int(12000 * Pi / FS)) * FS / (Pi * (1 + eps)) for Pi in P]))
    j = np.clip(np.searchsorted(harm, f), 1, len(harm) - 1); dist = np.minimum(np.abs(f - harm[j - 1]), np.abs(f - harm[j]))
    on, off = band & (dist < 0.3), band & (dist > 1.2)
    def taper(L):
        tw = np.ones(L); k = L // 20; tw[:k] = np.hanning(2 * k)[:k]; tw[-k:] = np.hanning(2 * k)[k:]; return tw
    tw = taper(len(nn))
    spec = lambda v, w=tw: np.abs(np.fft.rfft(v * w, Nf)) ** 2
    edges = np.arange(100, 11001, 250)
    def shaped(Sr, Sm):
        num = den = 0.0
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (f >= lo) & (f < hi); A = np.stack([Sm[k], np.ones(k.sum())], 1)
            (a_, b_), _ = nnls(A / Sr[k].mean(), Sr[k] / Sr[k].mean())
            num += a_ * Sm[k].sum(); den += Sr[k].sum()
        return float(num / den)
    oo = lambda S: float(10 * np.log10(S[on].mean() / S[off].mean()))
    out = {"window_s": list(win), "bins_on_frac": float(on.sum() / band.sum()), "bins_off_frac": float(off.sum() / band.sum()),
           "measured": [], "residual": [], "noise_tail": []}
    mn = (n >= noise_win[0] * FS) & (n < noise_win[1] * FS)
    for c in (0, 1):
        Sm = spec(meas[m, c] * fl); Sr = spec((meas - model)[m, c] * fl)
        Sn = spec(meas[mn, c], taper(mn.sum()))
        out["measured"].append({"on_off_db": oo(Sm)})
        out["residual"].append({"on_off_db": oo(Sr), "mode_shaped_frac": shaped(Sr, Sm)})
        out["noise_tail"].append({"on_off_db": oo(Sn), "mode_shaped_frac": shaped(Sn, Sm)})
    return out


def band_null(meas, model, a_ms, b_ms, bands=(63, 125, 250, 500, 1000, 2000, 4000, 8000, 11000)):
    from scipy.signal import sosfiltfilt
    i0, i1 = PRE + int(a_ms / 1000 * FS), PRE + int(b_ms / 1000 * FS); out = {}
    for fc in bands:
        lo, hi = (fc / np.sqrt(2), min(fc * np.sqrt(2), 12300)) if fc < 11000 else (10000, 12300)
        sos = butter(4, [lo, hi], btype="band", fs=FS, output="sos")
        yb = sosfiltfilt(sos, meas, axis=0)[i0:i1]; rb = sosfiltfilt(sos, meas - model, axis=0)[i0:i1]
        out[str(fc)] = [float(10 * np.log10(np.sum(rb[:, c] ** 2) / np.sum(yb[:, c] ** 2))) for c in (0, 1)]
    return out


def noise_floor_null(meas, windows):
    n = np.arange(len(meas)) - PRE
    nz = np.mean(meas[(n >= 1.6 * FS) & (n < 1.9 * FS)] ** 2, 0)
    out = {}
    for a, b in windows:
        m = (n >= a / 1000 * FS) & (n < b / 1000 * FS)
        out[f"{a}-{b}"] = [float(10 * np.log10(nz[c] * m.sum() / np.sum(meas[m, c] ** 2))) for c in (0, 1)]
    return out


def plot(path, y, models, fir_tab, metrics):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    main = models["B"]; fig, ax = plt.subplots(3, 2, figsize=(15, 13)); names = ["L", "R"]
    n = np.arange(len(y)) - PRE; t = n / FS * 1000
    blk = int(0.005 * FS)
    def env(s):
        k = len(s) // blk; return np.arange(k) * blk / FS * 1000 + 2.5 - PRE / FS * 1000, 10 * np.log10(np.mean(s[:k * blk].reshape(k, blk) ** 2, 1) + 1e-30)
    for c in range(2):
        a = ax[0, c]
        for s, lab, col in [(y[:, c], "measured", "k"), (main[:, c], "model B", "C0"), (y[:, c] - main[:, c], "residual (B)", "C3"), (y[:, c] - models["A"][:, c], "residual (A: loops only)", "C1")]:
            tt, e = env(s); a.plot(tt, e, color=col, lw=0.9, label=lab)
        a.set(xlim=(-4, 1200), ylim=(-110, 0), xlabel="ms from first cluster tap", ylabel="dB (5 ms RMS)", title=f"{names[c]}: energy envelope, measured vs model"); a.legend(fontsize=8); a.grid(alpha=.3)
        for x in (50, 300, 1000):
            a.axvline(x, color="grey", ls=":", lw=0.8)
    a = ax[1, 0]
    for c, mk in ((0, "o"), (1, "s")):
        a.semilogx(OCTAVES, [metrics["T30"]["measured"][c][str(f)] or np.nan for f in OCTAVES], "k" + mk + "-", label=f"measured {names[c]}", ms=5)
        a.semilogx(OCTAVES, [metrics["T30"]["B"][c][str(f)] or np.nan for f in OCTAVES], "C0" + mk + "--", label=f"model B {names[c]}", ms=5)
    a.set(xlabel="octave centre (Hz)", ylabel="T30 (s)", title="T30 per octave (device rate)"); a.legend(fontsize=8); a.grid(True, which="both", alpha=.3)
    a = ax[1, 1]
    for c, ls in ((0, "-"), (1, "--")):
        tm, vm = metrics["_ned"]["measured"][c]; tb, vb = metrics["_ned"]["B"][c]
        a.plot((tm - PRE / FS) * 1000, vm, "k" + ls, lw=1, label=f"measured {names[c]}"); a.plot((tb - PRE / FS) * 1000, vb, "C0" + ls, lw=1, label=f"model B {names[c]}")
    a.set(xlim=(0, 300), ylim=(0, 1.3), xlabel="ms from first cluster tap", title="Normalised echo density (20 ms window)"); a.legend(fontsize=8); a.grid(alpha=.3)
    a = ax[2, 0]; f = np.array(fir_tab["freq_hz"])
    for i, Pi in enumerate(P):
        a.plot(f, fir_tab["rt_s"][i], color=f"C{i}", label=f"P={Pi} RT from |G|")
    fm = np.linspace(20, 11000, 400); thB = np.array(metrics["loops"]["theta_jot"])
    Gp = G_param(thB, "jot", 2 * np.pi * fm / FS)
    a.plot(fm, -60 * P[0] / (FS * 20 * np.log10(np.abs(Gp[0]))), "k--", lw=1, label="parametric (Jot shelf)")
    a2 = a.twinx()
    for i, Pi in enumerate(P):
        a2.plot(f, fir_tab["tau_samples"][i], color=f"C{i}", ls=":", lw=1)
    a2.set_ylabel("effective length - P (samples, dotted)")
    a.set(xlabel="Hz", ylabel="per-loop RT (s)", title="Loop per-pass response (13-tap FIR estimate)"); a.legend(fontsize=7, loc="lower left"); a.grid(alpha=.3)
    a = ax[2, 1]; bn = metrics["band_null"]
    fb = [int(k) for k in bn["50-300"]]
    for key, col in (("50-300", "C0"), ("300-1000", "C2")):
        for c, mk in ((0, "o-"), (1, "s--")):
            a.semilogx(fb, [bn[key][str(k)][c] for k in fb], col + mk[0] + mk[1:], label=f"{key} ms {names[c]}", ms=4)
    a.set(xlabel="octave centre (Hz; 11k = 10-12.3 kHz)", ylabel="null depth (dB)", title="Model B residual per band"); a.legend(fontsize=8); a.grid(True, which="both", alpha=.3)
    fig.tight_layout(); fig.savefig(path, dpi=110)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", default="analysis/shortRT/h_dev.npy"); ap.add_argument("--out", default="analysis/model_s")
    ap.add_argument("--t-fit", type=float, default=1.0, help="fit data up to this time (s)")
    ap.add_argument("--t-cap", type=float, default=0.3, help="decay-flattening weight capped at this time (s)")
    ap.add_argument("--skip-distributor", action="store_true")
    ap.add_argument("--tighten", action="store_true", help="run only the tightening study (section 8 of the results doc)")
    ap.add_argument("--wav", default="irs/eugtone_EMT246_CLA_shortRT.wav")
    ap.add_argument("--tail", action="store_true", help="run only the long-tail study (section 9 of the results doc)")
    ap.add_argument("--stochastic", action="store_true", help="run only the time-variance / sweep simulation study (section 13)")
    ap.add_argument("--groups", action="store_true", help="run only the group-sparse tap study (section 12; needs tighten.json)")
    ap.add_argument("--allpass", action="store_true", help="run only the all-pass diffuser search (section 11; needs tighten.json)")
    ap.add_argument("--diffuser", action="store_true", help="run only the injection/diffuser study (section 10; needs tighten.json and tail.json)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    if a.tighten:
        return tighten(a)
    if a.tail:
        return tail_study(a)
    if a.diffuser:
        return diffuser_study(a)
    if a.allpass:
        return allpass_study(a)
    if a.groups:
        return group_study(a)
    if a.stochastic:
        return stochastic_study(a)
    y = load(a.h); out = {"h": a.h, "fs_dev": FS, "S0": S0, "t_fit_s": a.t_fit, "t_cap_s": a.t_cap}

    # 1. clock: cluster taps' sub-sample positions
    offs, eps_cl = cluster_subsample(y)
    out["cluster_subsample_offsets"] = {"L": offs[0], "R": offs[1]}; out["eps_cluster_ppm"] = eps_cl * 1e6
    print(f"cluster sub-sample drift: {eps_cl * 1e6:.1f} ppm")

    # 2. nonparametric per-loop per-pass response
    K = 6; taps, hist = fit_loops_fir(y, K=K, sweeps=8)
    fr = np.array([100, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000])
    Gf = G_fir(taps, K, 2 * np.pi * fr / FS)
    tab = {"freq_hz": fr.tolist(), "loss_db_per_pass": [], "tau_samples": [], "rt_s": [], "cost": hist[-1]}
    for i, G in enumerate(Gf):
        L = 20 * np.log10(np.abs(G)); tau = -np.unwrap(np.angle(G)) / (2 * np.pi * fr / FS)
        tab["loss_db_per_pass"].append(L.tolist()); tab["tau_samples"].append(tau.tolist())
        tab["rt_s"].append((-60 * P[i] / (FS * L)).tolist())
    tau_lf = np.array([np.mean(t[:3]) for t in tab["tau_samples"]])
    slope, icpt = np.polyfit(P, tau_lf, 1); resid = tau_lf - (slope * P + icpt)
    se = np.sqrt(np.sum(resid ** 2) / 2 / np.sum((P - P.mean()) ** 2))
    tab["tau_lf_fit"] = {"tau_lf": tau_lf.tolist(), "slope_ppm": slope * 1e6, "slope_se_ppm": se * 1e6, "intercept": icpt}
    out["loops_fir"] = tab
    print(f"tau(LF) per loop {np.round(tau_lf, 3)}; slope {slope * 1e6:.0f} +- {se * 1e6:.0f} ppm, intercept {icpt:.3f}")

    # 3. parametric loop filters: hypotheses for the offset and for Jot vs shared coefficients
    hyp = {}
    for name, kind, fix, th0 in [("jot_free", "jot", (), (0.99, 0.85, 0.3, 0.0, 100)),
                                 ("jot_clock_only", "jot", ("tau0",), (0.99, 0.85, 0.3, 0.0, 100)),
                                 ("jot_constant_offset_only", "jot", ("eps",), (0.99, 0.85, 0.3, 0.15, 0)),
                                 ("shared_free", "shared", (), (0.99, 0.85, 0.3, 0.0, 100)),
                                 ("shared_clock_at_cluster_eps", "shared", ("eps", "tau0"), (0.99, 0.85, 0.3, 0.0, eps_cl * 1e6)),
                                 ("jot_clock_at_cluster_eps", "jot", ("eps", "tau0"), (0.99, 0.85, 0.3, 0.0, eps_cl * 1e6)),
                                 ("onepole_jot", "onepole_jot", ("p", "tau0"), (0.99, 0.85, 0.0, 0.0, 100)),
                                 ("onepole_shared", "onepole_shared", ("rtpi", "tau0"), (0.99, 0.85, 0.04, 0.0, 100)),
                                 ("no_hf_damping", "jot", ("rtpi", "p"), (0.99, 0.99, 0.0, 0.0, 100))]:
        th, c = fit_loops_param(y, kind, th0, fix)
        hyp[name] = {"theta[rt0,rtpi,p,tau0,eps_ppm]": th.tolist(), "cost": c}
        print(f"{name:32s} cost {c:.4f} theta {np.round(th, 4)}")
    th = np.array(hyp["jot_free"]["theta[rt0,rtpi,p,tau0,eps_ppm]"]); eps = th[4] * 1e-6
    th5, g5, hf5, c4, c5 = fit_extra_loop(y, th, 460)
    out["loops"] = {"hypotheses": hyp, "theta_jot": th.tolist(),
                    "extra_460": {"theta": th5.tolist(), "g5": g5, "hf5": hf5, "cost_4": c4, "cost_5": c5,
                                  "rt5_s": float(-60 * 460 / (FS * 20 * np.log10(abs(g5))))}}
    print(f"extra recursion P=460: g={g5:.3f} hf={hf5:.3f} cost {c4:.4f} -> {c5:.4f}")

    # 4. linear stage, model variants
    Nfit = PRE + int(a.t_fit * FS); nn = np.arange(Nfit) - PRE
    wt = tail_weight(np.maximum(nn, 0), t_cap=a.t_cap); Nall = len(y)
    Gs = G_param(th, "jot"); models = {}; desc = {}; gains = {}
    def finish(key, lin_fit, v, d):
        lin_all = Linear(lin_fit._Gs, Nall, LF=lin_fit.LF, mode=lin_fit.mode, d=d, s0=lin_fit.s0, ext=lin_fit.Lx[0] - P[0], eps=eps,
                         Ps=lin_fit._Ps, post=lin_fit._post, post_f=lin_fit._post_f)
        models[key] = lin_all.render(v); gains[key] = lin_fit.gains(v).tolist()
    def run(key, text, Gs_, mode="per_channel", Ps=None, LF=1400, ext=300, d=None, post=None, post_f=False):
        lin = Linear(Gs_, Nfit, LF=LF, mode=mode, d=d, weight=wt, s0=S0, ext=ext, eps=eps, Ps=Ps, post=post, post_f=post_f)
        lin._Gs, lin._Ps, lin._post, lin._post_f = Gs_, Ps, post, post_f
        v = lin.solve(y, iters=300); finish(key, lin, v, d); desc[key] = text; return lin, v
    # A: loops only
    run("A", "joint LS, 4 loops, excitations P_i+300 from S0, no first-pass FIR", Gs, LF=0)
    # B: main model
    linB, vB = run("B", "joint LS, 4 loops + first-pass FIR (1400 samples from S0) [main model]", Gs)
    clB, fB, xsB = linB.split(vB)
    # C: B + parallel 460 loop
    G5 = G5_of(th5, g5, hf5, 460)
    run("C", "as B, refitted loops + parallel 5th recursion P=460 (exploratory)", G_param(th5, "jot") + [G5], Ps=list(P) + [460])
    C5 = 1 / (1 - G5 * np.exp(-1j * W * 460))
    run("Cs", "as C but the 460 recursion in series after the loop sum and the first-pass FIR (exploratory)", G_param(th5, "jot"), post=C5, post_f=True)
    # P: spec-literal periodic-decomposition excitations + measured first-pass residual
    xs_pd = pd_excitations(y, Gs, th)
    linP = Linear(Gs, Nall, LF=0, s0=S0, ext=0, eps=eps)
    vP = np.concatenate([(np.array([CLUSTER_SIGN[c] for c in (0, 1)], float) * np.array(gains["B"]) / np.array([CLUSTER_SIGN[c] for c in (0, 1)])).ravel() / linP.cs,
                         np.concatenate([x for row in xs_pd for x in row])])
    yP = linP.render(vP)
    # least-squares scale per loop and channel for the PD components (amplitudes only)
    comps = []
    for c in (0, 1):
        for i in range(4):
            vv = np.zeros_like(vP); o = 8 + int(np.sum(P)) * c + int(np.sum(P[:i])); vv[o:o + P[i]] = vP[o:o + P[i]]
            comps.append(linP.render(vv)[:, c])
    clus = linP.render(np.concatenate([vP[:8], np.zeros(len(vP) - 8)]))
    yPs = clus.copy(); m = (nn >= 0.1 * FS)
    for c in (0, 1):
        A = np.stack(comps[4 * c:4 * c + 4], 1)[:Nfit][m] * wt[m, None]
        sc = np.linalg.lstsq(A, (y[:Nfit, c] - clus[:Nfit, c])[m] * wt[m], rcond=None)[0]
        yPs[:, c] += np.stack(comps[4 * c:4 * c + 4], 1) @ sc
    fp = np.zeros_like(yPs); i0, i1 = PRE + S0, PRE + S0 + P[0]
    fp[i0:i1] = y[i0:i1] - yPs[i0:i1]           # first-pass FIR = measured minus the loops' first pass
    models["P"] = yPs + fp; desc["P"] = "periodic-decomposition excitations (amplitude-scaled) + measured first-pass residual until the first recirculation"
    # D: common loop sum + output distributor
    if not a.skip_distributor:
        linD, vD, dD, histD = fit_distributor(y, Gs, Nfit, wt, KD=700, rounds=5, s0=S0, ext=300, eps=eps, LF=1400)
        linD._Gs, linD._Ps, linD._post, linD._post_f = Gs, None, None, False
        finish("D", linD, vD, dD); desc["D"] = "common excitations x_i (L=R) + per-channel 701-tap output distributor, first-pass FIR"
        top = {c: [(int(i), float(dD[i, c])) for i in np.argsort(-np.abs(dD[:, c]))[:6]] for c in (0, 1)}
        out["distributor"] = {"als_cost_rel": (np.array(histD) / histD[0]).tolist(), "top_taps": {"L": top[0], "R": top[1]}}

    # 5. per-loop L->R relation on B's folded excitations
    rel = {}
    for i, Pi in enumerate(P):
        g = jot_gain(Pi, th[0]); ws = []
        for c in (0, 1):
            x = xsB[c][i]; w = x[:Pi].copy(); w[:len(x) - Pi] += x[Pi:] / g; ws.append(w)
        X = np.stack([np.roll(ws[0], k) for k in (26, 332)], 1); q = np.linalg.lstsq(X, ws[1], rcond=None)[0]; e = ws[1] - X @ q
        rel[str(Pi)] = {"gains_26_332": q.tolist(), "explained": float(1 - e @ e / (ws[1] @ ws[1]))}
    out["LR_two_lag_relation"] = rel
    out["LR_ratio_test"] = lr_ratio_test(xsB, th)
    fpe = [float(10 * np.log10(np.sum(fB[c] ** 2) / np.sum(y[PRE + S0:PRE + S0 + 1400, c] ** 2))) for c in (0, 1)]
    out["first_pass_fir_energy_rel_db"] = fpe

    # 6. metrics
    windows = WINDOWS_MS + [(0, 9), (9, 50), (300, 600), (600, 1000)]
    met = {"null_db": {}, "T30": {"measured": t30_octaves(y)}, "_ned": {"measured": ned(y)}, "desc": desc, "cluster_gains": gains}
    met["noise_limited_null_db"] = noise_floor_null(y, windows)
    nprof = np.arange(len(y)) - PRE
    met["late_floor_db"] = {f"{t:.1f}": [float(10 * np.log10(np.mean(y[(nprof >= t * FS) & (nprof < (t + 0.1) * FS), c] ** 2))) for c in (0, 1)]
                            for t in np.arange(1.0, 1.85, 0.1)}
    for k, m_ in models.items():
        met["null_db"][k] = null_depth(y, m_, windows)
        met["T30"][k] = t30_octaves(m_); met["_ned"][k] = ned(m_)
    met["residual_modes"] = {k: [residual_mode_test(y, models[k], th, win) for win in ((0.1, 0.6), (0.3, 0.8))] for k in ("B", "C")}
    met["band_null"] = {"50-300": band_null(y, models["B"], 50, 300), "300-1000": band_null(y, models["B"], 300, 1000)}
    met["ned_0.9_ms"] = {k: [float((t[np.argmax(v > 0.9)] - PRE / FS) * 1000) for t, v in met["_ned"][k]] for k in met["_ned"]}
    met["loops"] = out["loops"]
    out["metrics"] = {k: v for k, v in met.items() if not k.startswith("_") and k != "loops"}
    np.save(os.path.join(a.out, "model_B_dev.npy"), models["B"])
    np.savez(os.path.join(a.out, "model_B_params.npz"), theta_jot=th, cluster_gains=np.array(gains["B"]),
             first_pass_fir=fB, S0=S0, PRE=PRE, K0=K0, P=P,
             **{f"x_{'LR'[c]}_{Pi}": xsB[c][i] for c in (0, 1) for i, Pi in enumerate(P)})
    np.save(os.path.join(a.out, "model_C_dev.npy"), models["C"])
    np.save(os.path.join(a.out, "measured_dev.npy"), y)
    plot(os.path.join(a.out, "model_s.png"), y, models, tab, met)
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    for k in models:
        print(k, desc[k]); print("   null", {w: np.round(v, 1).tolist() for w, v in met["null_db"][k].items()})
        print("   T30 L", {f: (round(v, 3) if v else None) for f, v in met["T30"][k][0].items()})
    print("noise-limited", {w: np.round(v, 1).tolist() for w, v in met["noise_limited_null_db"].items()})
    print("T30 measured L", {f: (round(v, 3) if v else None) for f, v in met["T30"]["measured"][0].items()})
    print("T30 measured R", {f: (round(v, 3) if v else None) for f, v in met["T30"]["measured"][1].items()})
    print("T30 B R", {f: (round(v, 3) if v else None) for f, v in met["T30"]["B"][1].items()})
    print("NED 0.9 ms", met["ned_0.9_ms"]); print("band null", met["band_null"]); print("LR rel", rel, "first-pass FIR energy", fpe)



# ---------- capture-domain (96 kHz) comparison ----------
class Capture:
    """Device-rate model -> capture: band-limited interpolation onto the 96 kHz grid with the
    device clock stretched by (1+eps), then a free per-channel analogue-chain FIR K_c at 96 kHz.
    Compared with the raw WAV low-passed at f_cmp (zero-phase), so the deconvolution's
    regularisation and 11-12.3 kHz window play no part."""

    def __init__(self, wav, eps, n_dev, t_fit=0.4, f_cmp=11500.0, k_pre=96, k_post=704, half=24, l_offset=0.0):
        import soundfile as sf
        from deconv_grid import locate_first_tap
        x, fs = sf.read(wav, always_2d=True); self.fs = fs
        self.i0 = locate_first_tap(x[:, 1], fs) - int(0.1e-3 * fs)   # deconv_grid's kernel start = first tap
        Nf = 1 << int(np.ceil(np.log2(len(x) + 1))); f = np.fft.rfftfreq(Nf, 1 / fs)
        lp = np.where(f < f_cmp - 300, 1.0, np.where(f < f_cmp + 300, 0.5 * (1 + np.cos(np.pi * (f - f_cmp + 300) / 600)), 0.0))
        self.lp = lambda s: np.fft.irfft(np.fft.rfft(s, Nf, axis=0) * lp[:, None], Nf, axis=0)[:len(s)]
        self.x = self.lp(x)
        self.n_dev, self.eps, self.k_pre, self.k_post = n_dev, eps, k_pre, k_post
        # device sample i (index in model arrays, n = i - PRE) sits at 96k sample i0 + n*(1+eps)*fs/FS
        self.m_end = min(len(x), self.i0 + int((n_dev - PRE) * (1 + eps) * fs / FS))
        self.M = self.m_end - self.i0 + k_pre                        # 96k samples kept, starting at i0 - k_pre
        self.half = half
        self.S = [self._interp(l_offset), self._interp(0.0)]          # L may sit a fraction of a sample late
        self.ST = [S.T.tocsr() for S in self.S]
        self.Nk = 1 << int(np.ceil(np.log2(self.M + k_pre + k_post + 1)))
        self.K = np.zeros((k_pre + k_post + 1, 2)); self.K[k_pre, :] = 1.0
        self.target = self.x[self.i0 - k_pre:self.i0 - k_pre + self.M]
        t = (np.arange(self.M) - k_pre) / fs
        self.w = tail_weight(np.maximum(t * FS, 0), t_cap=0.3) * (t < t_fit)
        self.t = t

    def _interp(self, offset):
        from scipy.sparse import csr_matrix
        fs = self.fs; pos = (np.arange(self.n_dev) - PRE) * (1 + self.eps) * fs / FS + self.k_pre + offset
        rows, cols, vals = [], [], []
        r = self.half * fs / FS; fc = FS / 2 / fs                     # Kaiser-windowed sinc at fs_dev/2
        for i, p in enumerate(pos):
            m = np.arange(int(np.ceil(p - r)), int(np.floor(p + r)) + 1); m = m[(m >= 0) & (m < self.M)]
            d = m - p
            win = np.i0(8 * np.sqrt(np.clip(1 - (d / r) ** 2, 0, None))) / np.i0(8)
            rows.append(m); cols.append(np.full(len(m), i)); vals.append(2 * fc * np.sinc(2 * fc * d) * win)
        return csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(self.M, self.n_dev))

    def set_l_offset(self, offset):
        self.S[0] = self._interp(offset); self.ST[0] = self.S[0].T.tocsr()

    hp = ()      # first-order high-pass corners (Hz) of the analogue chain, applied before K

    def _hp(self):
        f = np.fft.rfftfreq(self.Nk, 1 / self.fs); H = np.ones(len(f), complex)
        for fc in self.hp:
            p = np.exp(-2 * np.pi * fc / self.fs); z1 = np.exp(-2j * np.pi * f / self.fs)
            H *= (1 + p) / 2 * (1 - z1) / (1 - p * z1)
        return H

    def interp(self, d):
        u = np.stack([self.S[c] @ d[:, c] for c in (0, 1)], 1)
        if self.hp:
            u = np.fft.irfft(np.fft.rfft(u, self.Nk, axis=0) * self._hp()[:, None], self.Nk, axis=0)[:self.M]
        return u

    def to_capture(self, d):
        """d: (n_dev, 2) device-grid signals -> (M, 2) capture-grid model."""
        u = self.interp(d)
        U = np.fft.rfft(u, self.Nk, axis=0); Kf = np.fft.rfft(self.K, self.Nk, axis=0)
        out = np.fft.irfft(U * Kf, self.Nk, axis=0)          # K index 0 = lag -k_pre
        return out[self.k_pre:self.k_pre + self.M]

    def from_capture(self, r):
        """Adjoint of to_capture."""
        buf = np.zeros((self.Nk, 2)); buf[self.k_pre:self.k_pre + self.M] = r
        Kf = np.fft.rfft(self.K, self.Nk, axis=0)
        V = np.fft.rfft(buf, axis=0) * np.conj(Kf)
        if self.hp:
            V = V * np.conj(self._hp())[:, None]
        v = np.fft.irfft(V, self.Nk, axis=0)[:self.M]
        return np.stack([self.ST[c] @ v[:, c] for c in (0, 1)], 1)

    def fit_K(self, d, ridge=1e-6, shared=True, w=None):
        """LS fit of the analogue kernel(s); shared=True fits one kernel for both channels."""
        u = self.interp(d); L = self.K.shape[0]; w = self.w if w is None else w
        w = np.broadcast_to(w if w.ndim == 2 else w[:, None], (self.M, 2))
        AtA, Atb = [], []
        for c in (0, 1):
            rows = np.nonzero(w[:, c] > 0)[0]
            upad = np.concatenate([np.zeros(self.k_post + 1), u[:, c], np.zeros(self.k_pre + 1)])
            X = np.stack([upad[rows + self.k_pre - j + self.k_post + 1] for j in range(L)], 1) * w[rows, c, None]
            AtA.append(X.T @ X); Atb.append(X.T @ (self.target[rows, c] * w[rows, c]))
        if shared:
            A = AtA[0] + AtA[1]; A = A + ridge * np.trace(A) / L * np.eye(L)
            self.K[:, 0] = self.K[:, 1] = np.linalg.solve(A, Atb[0] + Atb[1])
        else:
            for c in (0, 1):
                A = AtA[c] + ridge * np.trace(AtA[c]) / L * np.eye(L); self.K[:, c] = np.linalg.solve(A, Atb[c])

    def _tgt(self, band):
        return band_limit(self.target, band, self.fs) if band else self.target

    def null(self, model, windows, band=None):
        tg = self._tgt(band); md = band_limit(model, band, self.fs) if band else model; out = {}
        for a, b in windows:
            m = (self.t >= a / 1000) & (self.t < b / 1000)
            out[f"{a}-{b}"] = [float(10 * np.log10(np.sum((tg[m, c] - md[m, c]) ** 2) / np.sum(tg[m, c] ** 2))) for c in (0, 1)]
        return out

    def noise_null(self, windows, noise=(1.6, 1.9), band=None):
        tg = self._tgt(band); mn = (self.t >= noise[0]) & (self.t < noise[1])
        nz = np.mean(tg[mn] ** 2, 0); out = {}
        for a, b in windows:
            m = (self.t >= a / 1000) & (self.t < b / 1000)
            out[f"{a}-{b}"] = [float(10 * np.log10(nz[c] * m.sum() / np.sum(tg[m, c] ** 2))) for c in (0, 1)]
        return out


def capture_fit(cap, lin, rounds=3, iters=300, v0=None, shared=True, fit_k=True):
    """Alternate: device linear parameters (LSQR through the capture operator) and the kernels K_c."""
    from scipy.sparse.linalg import LinearOperator
    n = lin.size(); M = cap.M
    def mv(v):
        d = lin.render(v); return (cap.to_capture(d) * cap.w[:, None]).T.ravel()
    def rmv(r):
        r = r.reshape(2, M).T * cap.w[:, None]; dv = cap.from_capture(r)
        return lin.adjoint(dv.T.ravel())
    op = LinearOperator((2 * M, n), matvec=mv, rmatvec=rmv, dtype=float)
    b = (cap.target * cap.w[:, None]).T.ravel()
    v = v0
    hist = []
    for it in range(rounds):
        if v is not None and fit_k:
            cap.fit_K(lin.render(v), shared=shared)
        v = lsqr(op, b, iter_lim=iters, atol=1e-12, btol=1e-12, x0=v)[0]
        r = op.matvec(v) - b; hist.append(float(r @ r))
    if fit_k:
        cap.fit_K(lin.render(v), shared=shared)
    return v, hist


# ---------- tightening study (docs/model_s_results.md §8) ----------
def _excess(nd, floor):
    return {k: [nd[k][c] - floor[k][c] for c in (0, 1)] for k in nd}


def annihilator_structure_tests(y, th5, g5, hf5, win=(0.25, 0.75)):
    """What, beyond 4 loops + the 460 recursion, lowers the annihilator cost?"""
    Y = _fwd(y); n = np.arange(NFFT) - PRE
    sel = (n >= win[0] * FS) & (n < win[1] * FS); wt = tail_weight(n[sel])[:, None]
    A5 = 1 - G5_of(th5, g5, hf5, 460) * np.exp(-1j * W * 460); YY = Y * A5[:, None]
    cost = lambda e: [float(v) for v in np.sum(e ** 2, 0)]
    out = {"loops+460": cost(annihilate(YY, G_param(th5, "jot"))[sel] * wt)}
    # extra recursions
    def rec(g, hf, L, th): return 1 - g * shelf(W, 1.0, hf, 0.0) * np.exp(-1j * W * L * (1 + th[4] * 1e-6))
    for lags in ([230], [920], [690], [230, 690, 920]):
        def res(v):
            A = A5.copy()
            for k, L in enumerate(lags):
                A = A * rec(v[2 * k], v[2 * k + 1], L, th5)
            return (annihilate(Y * A[:, None], G_param(th5, "jot"))[sel] * wt).ravel()
        r = least_squares(res, np.tile([0.1, 1.0], len(lags)), bounds=([-.99, 0] * len(lags), [.99, 1.5] * len(lags)), diff_step=1e-4)
        out[f"+recursions {lags}"] = float(np.sum(r.fun ** 2))
    # free per-loop FIR (K=40, 81 taps)
    y5 = np.fft.irfft(YY, NFFT, axis=0)
    taps, _ = fit_loops_fir(y5, K=40, sweeps=6, win=win); out["per-loop 81-tap FIR"] = cost(annihilate(YY, G_fir(taps, 40))[sel] * wt)
    # per-loop free shelves
    def Gper(v):
        return [shelf(W, jot_gain(Pi, v[1 + 3 * i]), jot_gain(Pi, v[2 + 3 * i]), v[3 + 3 * i]) * np.exp(-1j * W * v[0] * 1e-6 * Pi) for i, Pi in enumerate(P)]
    r = least_squares(lambda v: (annihilate(YY, Gper(v))[sel] * wt).ravel(), np.r_[th5[4], np.tile(th5[:3], 4)],
                      x_scale=np.r_[10, np.tile([.01, .01, .01], 4)], diff_step=1e-4)
    out["per-loop shelves (12 params)"] = float(np.sum(r.fun ** 2))
    # LF shelf in every loop
    def Glf(v):
        return [G * shelf(W, jot_gain(Pi, v[0]) / jot_gain(Pi, th5[0]), 1.0, v[1]) for G, Pi in zip(G_param(th5, "jot"), P)]
    r = least_squares(lambda v: (annihilate(YY, Glf(v))[sel] * wt).ravel(), [th5[0], 0.99], bounds=([.2, 0], [5, .99995]), diff_step=1e-5)
    out["LF shelf"] = {"cost": float(np.sum(r.fun ** 2)), "rt_lf": float(r.x[0]), "corner_hz": float(-np.log(r.x[1]) * FS / (2 * np.pi))}
    # cross-coupling between loops (12 off-diagonal gains): annihilator = det(I - M D)
    G0 = G_param(th5, "jot"); E = [np.exp(-1j * W * Pi) for Pi in P]
    def detA(cv):
        Mx = np.zeros((len(W), 4, 4), complex); k = 0
        for i in range(4):
            for j in range(4):
                if i == j:
                    Mx[:, i, j] = 1 - G0[i] * E[j]
                else:
                    Mx[:, i, j] = -cv[k] * E[j]; k += 1
        return np.linalg.det(Mx)
    r = least_squares(lambda cv: (np.fft.irfft(YY * detA(cv)[:, None], NFFT, axis=0)[sel] * wt).ravel(),
                      np.random.default_rng(1).normal(0, 0.02, 12), x_scale=0.01, diff_step=1e-3)
    out["coupling (12 gains)"] = {"cost": float(np.sum(r.fun ** 2)), "max_abs_c": float(np.max(np.abs(r.x)))}
    # noise-only level of the same objective, from the measured 1.25-1.9 s tail placed at the window
    seg = y[PRE + int(1.25 * FS):PRE + int(1.9 * FS)]
    z = np.zeros((NFFT, 2)); i1 = PRE + int(win[1] * FS); z[max(0, i1 - len(seg)):i1] = seg[-(i1 - max(0, i1 - len(seg))):]
    out["noise_only"] = cost(annihilate(_fwd(z) * A5[:, None], G_param(th5, "jot"))[sel] * wt)
    # residual over noise vs time (50 ms blocks)
    e = annihilate(YY, G_param(th5, "jot")); z2 = np.zeros((NFFT, 2)); z2[PRE:PRE + len(seg)] = seg
    en = annihilate(_fwd(z2) * A5[:, None], G_param(th5, "jot"))
    nz = np.mean(en[(n > 0.3 * FS) & (n < 0.55 * FS)] ** 2, 0)
    out["resid_over_noise_db"] = {f"{t:.2f}": [float(v) for v in 10 * np.log10(np.mean(e[(n >= t * FS) & (n < (t + .05) * FS)] ** 2, 0) / nz)]
                                  for t in np.arange(0.25, 0.9, 0.05)}
    return out


def omp_curve(y, Gs, post, eps, Nfit, wt, checkpoints=(100, 400, 1600), ext=300, LF=1400):
    """Greedy (batch OMP) integer-tap excitations in the Cs basis; nulls at the checkpoints."""
    Nall = len(y)
    R = [np.fft.irfft(post / (1 - G * np.exp(-1j * W * Pi)), NFFT)[:Nall] for G, Pi in zip(Gs, P)]
    Rf = np.fft.irfft(post, NFFT)[:Nall]; kn, kb = cluster_kernels(eps)
    atoms = [(0, i, j) for i in range(4) for j in range(P[i] + ext)] + [(1, 0, j) for j in range(LF)]
    def col(c, a):
        v = np.zeros(Nall); st = PRE + S0 + a[2]; src = R[a[1]] if a[0] == 0 else Rf; v[st:] = src[:Nall - st]; return v
    def ccol(c, k):
        v = np.zeros(Nall); idx = PRE + CLUSTER[c][k] + kn; m = (idx >= 0) & (idx < Nall); v[idx[m]] = kb[c][k][m]; return v
    w2 = np.zeros(NFFT); w2[:Nfit] = wt[:Nfit] ** 2; W2 = np.fft.rfft(w2)
    spec = [np.fft.rfft(r, NFFT) for r in R + [Rf]]; lens = [P[i] + ext for i in range(4)] + [LF]
    nrm = np.sqrt(np.maximum(np.concatenate([np.fft.irfft(W2 * np.conj(np.fft.rfft(r ** 2, NFFT)), NFFT)[PRE + S0:PRE + S0 + L]
                                             for r, L in zip(R + [Rf], lens)]), 1e-30))
    out = {}
    for c in (0, 1):
        cols = [ccol(c, k) for k in range(4)]; sel = []; tgt = y[:, c]; r = tgt.copy(); res = []
        while len(sel) < checkpoints[-1]:
            rw = np.zeros(NFFT); rw[:Nfit] = r[:Nfit] * wt[:Nfit] ** 2; Rr = np.fft.rfft(rw)
            a = np.abs(np.concatenate([np.fft.irfft(Rr * np.conj(S), NFFT)[PRE + S0:PRE + S0 + L] for S, L in zip(spec, lens)])) / nrm
            a[sel] = 0
            for b in np.argsort(-a)[:max(5, len(sel) // 10)]:
                sel.append(int(b)); cols.append(col(c, atoms[b]))
            A = np.stack(cols, 1); coef = np.linalg.lstsq(A[:Nfit] * wt[:Nfit, None], tgt[:Nfit] * wt[:Nfit], rcond=None)[0]
            model = A @ coef; r = tgt - model
            for cp in checkpoints:
                if len(sel) >= cp and all(x[0] != cp for x in res):
                    res.append((cp, len(sel), null_depth(np.stack([tgt, tgt], 1), np.stack([model, model], 1), HOLD_WIN)))
        out["LR"[c]] = [{"target_taps": cp, "taps": k, "null": {w: v[0] for w, v in nd.items()}} for cp, k, nd in res]
    return out


HOLD_WIN = [(0, 50), (50, 300), (300, 400), (400, 600), (600, 1000)]
CMP_BAND = (20.0, 11300.0)   # comparison band: the sweep/deconvolution carries no usable information below ~20 Hz


def band_limit(x, band=CMP_BAND, fs=FS):
    """Zero-phase FFT band-pass with raised-cosine edges (lower edge 10-20 Hz, upper edge +-300 Hz)."""
    Nf = 1 << int(np.ceil(np.log2(2 * len(x)))); f = np.fft.rfftfreq(Nf, 1 / fs); lo, hi = band
    H = np.clip((f - lo / 2) / (lo / 2), 0, 1); H = 0.5 - 0.5 * np.cos(np.pi * H)
    H *= np.where(f < hi - 300, 1.0, np.where(f < hi + 300, 0.5 * (1 + np.cos(np.pi * (f - hi + 300) / 600)), 0.0))
    return np.fft.irfft(np.fft.rfft(x, Nf, axis=0) * (H[:, None] if x.ndim == 2 else H), Nf, axis=0)[:len(x)]


def tighten(a):
    """Everything in docs/model_s_results.md §8. Writes <out>/tighten.json."""
    import time
    y = load(a.h); Nall = len(y); res = {}; T = time.time()
    floor = noise_floor_null(y, HOLD_WIN); res["noise_floor_dev"] = floor
    yb = band_limit(y); floor_b = noise_floor_null(yb, HOLD_WIN); res["noise_floor_dev_band"] = floor_b
    # loops identified from data before 0.4 s only
    th, c4 = fit_loops_param(y, "jot", (0.99, 0.85, 0.03, 0.0, 100), (), win=(0.25, 0.40))
    th5, g5, hf5, _, c5 = fit_extra_loop(y, th, 460, win=(0.25, 0.40))
    res["loops_early"] = {"theta4": th.tolist(), "cost4": c4, "theta5": th5.tolist(), "g5": g5, "hf5": hf5, "cost5": c5}
    print("loops (0.25-0.40 s):", np.round(th5, 4), round(g5, 3), round(hf5, 3), f"{time.time() - T:.0f}s", flush=True)
    # full-window identification for the structure tests
    thF, _ = fit_loops_param(y, "jot", (0.99, 0.85, 0.03, 0.0, 100), ())
    thF5, gF5, hfF5, _, _ = fit_extra_loop(y, thF, 460)
    res["structure"] = annihilator_structure_tests(y, thF5, gF5, hfF5)
    print("structure tests done", f"{time.time() - T:.0f}s", flush=True)
    # cross-channel: loop + 460 parameters from one channel's annihilator, cost on the other
    Y = _fwd(y); n = np.arange(NFFT) - PRE; sel = (n >= .25 * FS) & (n < .75 * FS); wt_a = tail_weight(n[sel])
    def ann_cost(v, ch):
        A5 = 1 - G5_of(v[:5], v[5], v[6], 460) * np.exp(-1j * W * 460)
        return annihilate(Y * A5[:, None], G_param(v[:5], "jot"))[sel, ch] * wt_a
    cc = {}
    for ch in (0, 1):
        r = least_squares(lambda v: ann_cost(v, ch), np.r_[thF5, gF5, hfF5], x_scale=np.array([.01, .01, .01, .01, 10, .05, .05]), diff_step=1e-4)
        cc["LR"[ch]] = r.x
    res["cross_channel"] = {f"params_from_{s}": {"theta+g5,hf5": cc[s].tolist(), "cost_L": float(np.sum(ann_cost(cc[s], 0) ** 2)),
                                                 "cost_R": float(np.sum(ann_cost(cc[s], 1) ** 2))} for s in "LR"}
    # hold-out table (device domain): fit 0-0.4 s, predict 0.4-1.0 s
    Nfit = PRE + int(0.4 * FS); nn = np.arange(Nfit) - PRE; wt = tail_weight(np.maximum(nn, 0), t_cap=0.3); eps = th5[4] * 1e-6
    Gs, Gs5, G5 = G_param(th, "jot"), G_param(th5, "jot"), G5_of(th5, g5, hf5, 460); C5 = 1 / (1 - G5 * np.exp(-1j * W * 460))
    rows = []
    def run(name, Gs_, nlp, Ps=None, post=None, post_f=False, ext=300, LF=1400, eps_=eps, gs=None):
        lin = Linear(Gs_, Nfit, LF=LF, weight=wt, s0=S0, ext=ext, eps=eps_, Ps=Ps, post=post, post_f=post_f); v = lin.solve(y, iters=400)
        la = Linear(Gs_, Nall, LF=LF, s0=S0, ext=ext, eps=eps_, Ps=Ps, post=post, post_f=post_f); out = la.render(v)
        nd = null_depth(y, out, HOLD_WIN); k = 8 + 2 * LF + 2 * sum(la.Lx) + nlp
        ndb = null_depth(yb, band_limit(out), HOLD_WIN)
        rows.append({"model": name, "params": k, "null": nd, "excess": _excess(nd, floor), "null_band": ndb, "excess_band": _excess(ndb, floor_b)})
        print(f"{name:34s} {k:6d}", {w: np.round(v_, 1).tolist() for w, v_ in nd.items()}, flush=True)
        return out
    run("B: 4 loops", Gs, 5)
    run("C: + parallel 460", Gs5 + [G5], 8, Ps=list(P) + [460])
    run("Cs: + series 460", Gs5, 8, post=C5, post_f=True)
    run("Cs, no first-pass FIR", Gs5, 8, post=C5, post_f=True, LF=0)
    run("Cs, R loops from L-only fit", G_param(cc["L"][:5], "jot"), 8,
        post=1 / (1 - G5_of(cc["L"][:5], cc["L"][5], cc["L"][6], 460) * np.exp(-1j * W * 460)), post_f=True, eps_=cc["L"][4] * 1e-6)
    res["holdout_dev"] = rows
    # sparse integer-tap excitations (Cs basis)
    res["omp"] = omp_curve(y, Gs5, C5, eps, Nfit, np.r_[wt, np.zeros(Nall - Nfit)])
    print("omp", {c: [(d["taps"], round(d["null"]["50-300"], 1)) for d in v] for c, v in res["omp"].items()}, f"{time.time() - T:.0f}s", flush=True)
    # capture (96 kHz) domain
    delta = np.zeros(2 * PRE + 1); delta[PRE] = 1.0
    cap = Capture(a.wav, eps, Nall, t_fit=0.4, k_pre=192, k_post=1000)
    reg = np.stack([(cap.t > -0.0015) & (cap.t < 0.009), (cap.t > -0.0015) & (cap.t < 0.021)], 1).astype(float)
    nzc = [np.mean(cap.target[(cap.t > 1.6) & (cap.t < 1.9), c] ** 2) for c in (0, 1)]
    def cl_dev(g):
        d = np.zeros((Nall, 2))
        for c in (0, 1):
            for t, gg in zip(CLUSTER[c], g[c]): d[PRE + t, c] = gg
        return d
    def cl_gains(g, c):
        cols = []
        for t in CLUSTER[c]:
            d = np.zeros((Nall, 2)); d[PRE + t, c] = 1; cols.append(cap.to_capture(d)[:, c])
        m = reg[:, c] > 0; g[c] = np.linalg.lstsq(np.stack(cols, 1)[m], cap.target[m, c], rcond=None)[0]
    def cl_null(g):
        out = cap.to_capture(cl_dev(g)); return [float(10 * np.log10(np.sum((cap.target[reg[:, c] > 0, c] - out[reg[:, c] > 0, c]) ** 2) / np.sum(cap.target[reg[:, c] > 0, c] ** 2))) for c in (0, 1)]
    clus = {"noise_floor": [float(10 * np.log10(nzc[c] * reg[:, c].sum() / np.sum(cap.target[reg[:, c] > 0, c] ** 2))) for c in (0, 1)], "trials": []}
    wR = reg.copy(); wR[:, 0] = 0
    for hp in ((), (10.0,), (10.0, 10.0), (15.0, 15.0)):
        cap.hp = hp; g = {c: np.array(CLUSTER_SIGN[c], float) for c in (0, 1)}
        for _ in range(3):
            cap.fit_K(cl_dev(g), ridge=1e-8, w=wR); cl_gains(g, 1)
        off = minimize_scalar(lambda o: (cap.set_l_offset(o), cl_gains(g, 0), cl_null(g)[0])[2], bounds=(-4, 4), method="bounded", options={"xatol": 1e-3}).x
        cap.set_l_offset(off); cl_gains(g, 0); hold = cl_null(g)
        clus["trials"].append({"hp_hz": list(hp), "R_in_sample": hold[1], "L_holdout": hold[0], "L_offset_dev_samples": float(off * FS / cap.fs),
                               "gains": {s: g[c].tolist() for c, s in ((0, "L"), (1, "R"))}})
        print("cluster hp", hp, np.round(hold, 1), round(off * FS / cap.fs, 3), np.round(g[1], 4), flush=True)
    res["capture_cluster"] = clus
    best = min(clus["trials"], key=lambda d: d["R_in_sample"]); cap.hp = tuple(best["hp_hz"]); cap.set_l_offset(best["L_offset_dev_samples"] * cap.fs / FS)
    g = {0: np.array(best["gains"]["L"]), 1: np.array(best["gains"]["R"])}
    for _ in range(3):
        cap.fit_K(cl_dev(g), ridge=1e-8, w=reg); cl_gains(g, 0); cl_gains(g, 1)
    # full model Cs in the capture domain: integer loops, delta cluster, clock in the interpolation
    th0 = th5.copy(); th0[4] = 0.0; C50 = 1 / (1 - G5_of(th0, g5, hf5, 460) * np.exp(-1j * W * 460))
    lin = Linear(G_param(th0, "jot"), Nall, LF=1400, s0=S0, ext=300, eps=0.0, kernel=delta, post=C50, post_f=True)
    lf = Linear(Gs5, Nfit, LF=1400, weight=wt, s0=S0, ext=300, eps=eps, post=C5, post_f=True); v0 = lf.solve(y, iters=300)
    v0 = v0.copy(); v0[:8] = np.concatenate([g[0], g[1]]) / lin.cs
    win_c = [(0, 9)] + HOLD_WIN; floor_c = cap.noise_null(win_c)
    res["capture_full"] = {"noise_floor": floor_c, "hp_hz": list(cap.hp), "L_offset_dev_samples": best["L_offset_dev_samples"]}
    # (a) kernel fixed from the cluster-only fit; only the device parameters are solved (linear, fit 0-0.4 s)
    v, h = capture_fit(cap, lin, rounds=2, iters=300, v0=v0, fit_k=False)
    mod_c = cap.to_capture(lin.render(v)); nd = cap.null(mod_c, win_c)
    floor_cb = cap.noise_null(win_c, band=CMP_BAND); ndb = cap.null(mod_c, win_c, band=CMP_BAND); res["capture_full"]["noise_floor_band"] = floor_cb
    res["capture_full"]["kernel_fixed"] = {"null": nd, "excess": _excess(nd, floor_c), "null_band": ndb, "excess_band": _excess(ndb, floor_cb),
                                           "cost_hist": h, "cluster_gains": lin.gains(v).tolist()}
    print("capture full, kernel fixed, 20 Hz-11.3 kHz", {w: np.round(x, 1).tolist() for w, x in ndb.items()}, flush=True)
    print("capture full, kernel fixed", {w: np.round(x, 1).tolist() for w, x in nd.items()}, f"{time.time() - T:.0f}s", flush=True)
    np.save(os.path.join(a.out, "capture_model_Cs.npy"), cap.to_capture(lin.render(v)))
    # (b) kernel refitted jointly with the reverb (alternating), for comparison
    for r_ in range(3):
        v, h = capture_fit(cap, lin, rounds=1, iters=150, v0=v)
    mod_c = cap.to_capture(lin.render(v)); nd = cap.null(mod_c, win_c); ndb = cap.null(mod_c, win_c, band=CMP_BAND)
    res["capture_full"]["kernel_joint"] = {"null": nd, "excess": _excess(nd, floor_c), "null_band": ndb, "excess_band": _excess(ndb, floor_cb)}
    print("capture full, kernel joint", {w: np.round(x, 1).tolist() for w, x in nd.items()}, f"{time.time() - T:.0f}s", flush=True)
    # line width: measured vs LTI model (+ measured noise) spectral concentration at loop modes
    mC = Linear(Gs5 + [G5], Nall, LF=1400, s0=S0, ext=300, eps=eps, Ps=list(P) + [460])
    lC = Linear(Gs5 + [G5], Nfit, LF=1400, weight=wt, s0=S0, ext=300, eps=eps, Ps=list(P) + [460])
    model = mC.render(lC.solve(y, iters=400))
    seg = y[PRE + int(1.3 * FS):PRE + int(1.9 * FS)]
    lw = {}
    for w0 in (0.1, 0.4):
        nz = np.zeros_like(model); i0 = PRE + int(w0 * FS); nz[i0:i0 + len(seg)] = seg[:len(nz) - i0]
        rm = residual_mode_test(y, np.zeros_like(y), th5, (w0, w0 + 0.3))["residual"]
        rn = residual_mode_test(model + nz, np.zeros_like(y), th5, (w0, w0 + 0.3))["residual"]
        lw[f"{w0:.1f}-{w0 + .3:.1f}"] = {"measured_on_off_db": [d["on_off_db"] for d in rm], "model_plus_noise_on_off_db": [d["on_off_db"] for d in rn]}
    res["linewidth"] = lw
    with open(os.path.join(a.out, "tighten.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print("linewidth", lw, f"total {time.time() - T:.0f}s")



# ---------- long-tail study (docs/model_s_results.md §9) ----------
Q15 = 1 / 32768


def G_perloop(v, w=W):
    """v = [g1..g4 (per-pass DC gain per loop), kappa (Jot HF excess, dB/s at Nyquist), p, eps_ppm].
    Nyquist/DC ratio of loop i = 10^(-kappa*P_i/(20 fs)), i.e. Jot-scaled HF excess."""
    g, kappa, p, eps = v[:4], v[4], v[5], v[6] * 1e-6
    return [g[i] * shelf(w, 1.0, 10 ** (-kappa * Pi / (20 * FS)), p) * np.exp(-1j * w * eps * Pi) for i, Pi in enumerate(P)]


def fit_tail_loops(y, chans=(0, 1), win=(0.25, 1.8), irls=3, v0=None, block=0.05):
    """Loops + 460 recursion from the annihilator over the whole decay. The residual of a correct model
    is filtered measurement noise (stationary), so it is fitted unweighted, then reweighted by the
    inverse local residual RMS (50 ms blocks) so stretches where the model is off or SNR is low count less.
    Returns params [g1..g4, kappa, p, eps_ppm, g5, hf5], 1-sigma (Jacobian, inflated for residual
    correlation), cost and the weights."""
    Y = _fwd(y); n = np.arange(NFFT) - PRE; sel = (n >= win[0] * FS) & (n < win[1] * FS); ns = n[sel]
    blk = ((ns - ns[0]) // int(block * FS)).astype(int)
    if v0 is None:
        v0 = np.r_[[jot_gain(Pi, 1.0) for Pi in P], 60 * (1 / 0.82 - 1 / 1.0), 0.04, 98.0, 0.44, 0.9]
    wgt = np.ones((sel.sum(), len(chans)))
    def res(v):
        A5 = 1 - G5_of(np.r_[0, 0, 0, 0, v[6]], v[7], v[8], 460) * np.exp(-1j * W * 460)
        return (annihilate(Y * A5[:, None], G_perloop(v))[sel][:, chans] * wgt).ravel()
    v = np.array(v0, float)
    xs = np.r_[[1e-4] * 4, 0.1, 0.01, 1.0, 0.01, 0.01]
    for it in range(irls):
        r = least_squares(res, v, x_scale=xs, diff_step=1e-6); v = r.x
        e = r.fun.reshape(-1, len(chans)) / wgt
        for c in range(len(chans)):
            rms = np.sqrt(np.bincount(blk, e[:, c] ** 2) / np.bincount(blk))
            wgt[:, c] = 1 / np.maximum(rms[blk], 1e-12)
    r = least_squares(res, v, x_scale=xs, diff_step=1e-6)
    J = r.jac; s2 = np.sum(r.fun ** 2) / (len(r.fun) - len(v))
    # integrated autocorrelation time of the whitened residual (per channel, lags up to 400)
    taus = []
    for c in range(len(chans)):
        e = r.fun.reshape(-1, len(chans))[:, c]; Nf = 1 << int(np.ceil(np.log2(2 * len(e))))
        ac = np.fft.irfft(np.abs(np.fft.rfft(e, Nf)) ** 2, Nf)[:401]; ac /= ac[0]
        taus.append(1 + 2 * np.sum(ac[1:]))
    tau = float(max(1.0, np.mean(taus)))
    cov = s2 * tau * np.linalg.inv(J.T @ J)
    return r.x, np.sqrt(np.diag(cov)), float(np.sum(r.fun ** 2)), tau


def noise_model(y, t_ref=(1.4, 1.9), t_prof=(1.0, 1.9)):
    """Additive measurement-noise model per channel: octave-band spectrum from t_ref and a level
    profile in dB vs time (linear fit over t_prof, 100 ms blocks). R is flat; L drifts down."""
    n = np.arange(len(y)) - PRE; out = {"band_db": [], "slope_db_per_s": [], "level_db_at_1.4s": [], "blocks": {}}
    ts = np.arange(t_prof[0], t_prof[1] - 1e-9, 0.1)
    for c in (0, 1):
        lv = [10 * np.log10(np.mean(y[(n >= t * FS) & (n < (t + .1) * FS), c] ** 2)) for t in ts]
        # fit only where the reverb has decayed below the floor (>= 1.3 s for L, >= 1.4 s for R)
        m = ts >= (1.3 if c == 0 else 1.4); sl, ic = np.polyfit(ts[m] + 0.05, np.array(lv)[m], 1)
        out["slope_db_per_s"].append(float(sl)); out["level_db_at_1.4s"].append(float(ic + sl * 1.4))
        out["blocks"]["LR"[c]] = [float(x) for x in lv]
        seg = y[(n >= t_ref[0] * FS) & (n < t_ref[1] * FS), c]; S = np.abs(np.fft.rfft(seg)) ** 2; f = np.fft.rfftfreq(len(seg), 1 / FS)
        out["band_db"].append({str(fc): float(10 * np.log10(S[(f >= fc / np.sqrt(2)) & (f < min(fc * np.sqrt(2), 12300))].sum() / S.sum()))
                               for fc in OCTAVES})
    return out


def noise_energy(y, nm, a_s, b_s, c, drift=True, band=None):
    """Expected noise energy in [a_s, b_s) for channel c under the noise model (optionally in one band:
    fraction of the reference-segment energy in that band)."""
    n = np.arange(len(y)) - PRE; ref = np.mean(y[(n >= 1.4 * FS) & (n < 1.9 * FS), c] ** 2)
    t = np.arange(int(a_s * FS), int(b_s * FS)) / FS
    lvl = 10 ** ((nm["slope_db_per_s"][c] * (t - 1.65)) / 10) if (drift and c == 0) else np.ones(len(t))
    e = ref * lvl.sum()
    if band is not None:
        e *= 10 ** (nm["band_db"][c][str(band)] / 10)
    return float(e)


def tail_study(a):
    """Everything in docs/model_s_results.md §9. Writes <out>/tail.json."""
    from scipy.signal import sosfiltfilt
    import time
    T = time.time(); y = load(a.h); Nall = len(y); res = {}
    # 1. loop parameters from the whole decay, SNR-weighted, with uncertainties
    v, sd, cost, tau = fit_tail_loops(y)
    names = ["g1321", "g1389", "g1461", "g1587", "kappa_db_per_s", "p", "eps_ppm", "g460", "hf460"]
    res["tail_fit_both"] = {"params": dict(zip(names, v.tolist())), "sigma": dict(zip(names, sd.tolist())), "cost": cost, "tau_int": tau,
                            "sigma_in_q15_lsb": {k: float(s / Q15) for k, s in zip(names, sd) if k.startswith("g") or k == "p"},
                            "rt0_per_loop": [float(-3 * Pi / (FS * np.log10(g))) for g, Pi in zip(v[:4], P)]}
    print("tail fit", np.round(v, 6), "sd", np.round(sd, 6), "tau", round(tau, 1), f"{time.time() - T:.0f}s", flush=True)
    splits = {"L_only": dict(chans=(0,)), "R_only": dict(chans=(1,)), "early_0.25-0.7": dict(win=(0.25, 0.7)), "late_0.7-1.8": dict(win=(0.7, 1.8))}
    res["splits"] = {}
    for k, kw in splits.items():
        vs, sds, _, _ = fit_tail_loops(y, v0=v, **kw); res["splits"][k] = {"params": dict(zip(names, vs.tolist())), "sigma": dict(zip(names, sds.tolist()))}
        print(k, np.round(vs, 6), flush=True)
    # 2. model C fitted on 0-0.4 s, compared over the tail (hold-out)
    eps = v[6] * 1e-6; Gs = G_perloop(v); G5 = G5_of(np.r_[0, 0, 0, 0, v[6]], v[7], v[8], 460)
    Nfit = PRE + int(0.4 * FS); nn = np.arange(Nfit) - PRE; wt = tail_weight(np.maximum(nn, 0), t_cap=0.3)
    kw = dict(LF=1400, s0=S0, ext=300, eps=eps, Ps=list(P) + [460])
    model = Linear(Gs + [G5], Nall, **kw).render(Linear(Gs + [G5], Nfit, weight=wt, **kw).solve(y, iters=400))
    np.save(os.path.join(a.out, "model_C_tailfit_dev.npy"), model)
    # 2a. per-band decay-rate and frequency deviation: free 81-tap per-loop FIR (whole decay) vs parametric
    A5 = 1 - G5 * np.exp(-1j * W * 460); y5 = np.fft.irfft(_fwd(y) * A5[:, None], NFFT, axis=0)
    taps, _ = fit_loops_fir(y5, K=40, sweeps=6, win=(0.25, 1.2))
    fr = np.array([125, 250, 500, 1000, 2000, 4000, 8000, 10000]); wf = 2 * np.pi * fr / FS
    Gf, Gp = G_fir(taps, 40, wf), G_perloop(v, wf)
    res["per_band_dev"] = {"freq_hz": fr.tolist(), "decay_db_per_s": [], "freq_mhz": []}
    for i, Pi in enumerate(P):
        res["per_band_dev"]["decay_db_per_s"].append((20 * np.log10(np.abs(Gf[i] / Gp[i])) * FS / Pi).tolist())
        dtau = -np.angle(Gf[i] / Gp[i]) / wf                     # extra effective length (samples)
        res["per_band_dev"]["freq_mhz"].append((-1000 * fr * dtau / Pi).tolist())
    # 2b. noise-compensated EDC per octave, down to -50 dB
    nm = noise_model(y); res["noise_model"] = nm
    ed = {}
    for fc in OCTAVES:
        sos = butter(4, [fc / np.sqrt(2), min(fc * np.sqrt(2), 12300)], btype="band", fs=FS, output="sos")
        row = []
        for c in (0, 1):
            dm, _ = edc(sosfilt(sos, y[:, c]), FS); dmod, _ = edc(sosfilt(sos, model[:, c]), FS)
            L = min(len(dm), len(dmod)); m = dm[:L] >= -50
            row.append({"max_abs_db": float(np.max(np.abs(dm[:L][m] - dmod[:L][m]))), "reach_db": float(dm[:L].min())})
        ed[str(fc)] = row
    res["edc_dev"] = ed
    # 2c. fine-structure correlation, model vs measured, successive windows (+ ceiling with a noise realisation)
    n = np.arange(Nall) - PRE; seg = y[(n >= 1.3 * FS) & (n < 1.9 * FS)][::-1]
    def fs_spec(s):
        X = np.abs(np.fft.rfft(s * np.hanning(len(s)), 1 << 19)); f = np.fft.rfftfreq(1 << 19, 1 / FS)
        Lg = 20 * np.log10(X[(f > 200) & (f < 8000)] + 1e-15); k = int(50 / (f[1] - f[0])); return Lg - np.convolve(Lg, np.ones(k) / k, "same")
    fsc = {}
    for t0 in (0.3, 0.6, 0.9, 1.2):
        m = (n >= t0 * FS) & (n < (t0 + .3) * FS); nz = seg[:m.sum()]
        fsc[f"{t0:.1f}-{t0 + .3:.1f}"] = {"model_vs_meas": [float(np.corrcoef(fs_spec(model[m, c]), fs_spec(y[m, c]))[0, 1]) for c in (0, 1)],
                                          "ceiling_model_vs_model+noise": [float(np.corrcoef(fs_spec(model[m, c]), fs_spec(model[m, c] + nz[:, c]))[0, 1]) for c in (0, 1)]}
    res["fine_structure"] = fsc
    # 3. model error = residual - noise model, per window and band
    wins = [(0.05, 0.3), (0.3, 0.4), (0.4, 0.6), (0.6, 1.0), (1.0, 1.5)]
    me = {}
    for a_, b_ in wins:
        m = (n >= a_ * FS) & (n < b_ * FS); row = {}
        for c in (0, 1):
            er = np.sum((y[m, c] - model[m, c]) ** 2); es = np.sum(y[m, c] ** 2)
            for drift in ((True, False) if c == 0 else (False,)):
                en = noise_energy(y, nm, a_, b_, c, drift)
                key = "LR"[c] + ("" if c else ("_drift" if drift else "_flat"))
                row[key] = {"null_db": float(10 * np.log10(er / es)), "noise_db": float(10 * np.log10(en / es)),
                            "model_error_db": float(10 * np.log10(max(er - en, 1e-30) / es)), "residual_over_noise_db": float(10 * np.log10(er / en))}
        me[f"{int(a_ * 1000)}-{int(b_ * 1000)}"] = row
    res["model_error"] = me
    yb, mb_ = band_limit(y), band_limit(model); nmb = noise_model(yb); meb2 = {}
    for a_, b_ in wins:
        m = (n >= a_ * FS) & (n < b_ * FS); row = {}
        for c in (0, 1):
            er = np.sum((yb[m, c] - mb_[m, c]) ** 2); es = np.sum(yb[m, c] ** 2)
            for drift in ((True, False) if c == 0 else (False,)):
                en = noise_energy(yb, nmb, a_, b_, c, drift)
                key = "LR"[c] + ("" if c else ("_drift" if drift else "_flat"))
                row[key] = {"null_db": float(10 * np.log10(er / es)), "noise_db": float(10 * np.log10(en / es)),
                            "model_error_db": float(10 * np.log10(max(er - en, 1e-30) / es)), "residual_over_noise_db": float(10 * np.log10(er / en))}
        meb2[f"{int(a_ * 1000)}-{int(b_ * 1000)}"] = row
    res["model_error_20Hz_11k3"] = meb2; res["noise_model_20Hz_11k3"] = nmb
    # where the residual sits in frequency (share of total signal energy), 50-300 and 300-600 ms
    split = {}
    for a_, b_ in ((0.05, 0.3), (0.3, 0.6)):
        m = (n >= a_ * FS) & (n < b_ * FS); X = np.fft.rfft((y - model)[m], axis=0); Yf = np.fft.rfft(y[m], axis=0)
        f = np.fft.rfftfreq(m.sum(), 1 / FS); tot = np.sum(np.abs(Yf) ** 2, 0)
        split[f"{int(a_ * 1000)}-{int(b_ * 1000)}"] = {f"{lo}-{hi}": {"residual_db": (10 * np.log10(np.sum(np.abs(X[(f >= lo) & (f < hi)]) ** 2, 0) / tot)).tolist(),
                                                                        "signal_db": (10 * np.log10(np.sum(np.abs(Yf[(f >= lo) & (f < hi)]) ** 2, 0) / tot)).tolist()}
                                                       for lo, hi in ((0, 20), (20, 44), (44, 11300), (11300, 13187))}
    res["residual_by_freq"] = split
    mb = {}
    for a_, b_ in ((0.05, 0.3), (0.3, 0.6), (0.6, 1.0)):
        m = (n >= a_ * FS) & (n < b_ * FS); row = {}
        for fc in OCTAVES:
            sos = butter(4, [fc / np.sqrt(2), min(fc * np.sqrt(2), 12300)], btype="band", fs=FS, output="sos")
            yb = sosfiltfilt(sos, y, axis=0)[m]; rb = sosfiltfilt(sos, y - model, axis=0)[m]
            row[str(fc)] = [float(10 * np.log10(max(np.sum(rb[:, c] ** 2) - noise_energy(y, nm, a_, b_, c, True, fc), 1e-30) / np.sum(yb[:, c] ** 2))) for c in (0, 1)]
        mb[f"{int(a_ * 1000)}-{int(b_ * 1000)}"] = row
    res["model_error_band"] = mb
    with open(os.path.join(a.out, "tail.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("model_error", "edc_dev", "fine_structure")}, default=lambda x: round(x, 2))[:3000])
    print(f"total {time.time() - T:.0f}s")



# ---------- injection / diffuser study (docs/model_s_results.md §10) ----------
def predict_from_state(y, D, T):
    """If the input has ended by T, D(z) y = 0 after T, so y after T follows from y before T alone:
    y_fut = -(1/D) * [(D * y_<T) restricted to n >= T]. D = product of all recursion factors."""
    yt = np.zeros((NFFT, 2)); iT = PRE + int(T * FS); yt[:iT] = y[:iT]
    e = np.fft.irfft(np.fft.rfft(yt, axis=0) * D[:, None], NFFT, axis=0); e[:iT] = 0
    p = -np.fft.irfft(np.fft.rfft(e, axis=0) / D[:, None], NFFT, axis=0)
    out = yt.copy(); out[iT:] = p[iT:]; return out[:len(y)]


def recursion_poly(Gs, extra=()):
    D = np.ones(len(W), complex)
    for G, Pi in zip(Gs, P):
        D *= 1 - G * np.exp(-1j * W * Pi)
    for G5, L in extra:
        D *= 1 - G5 * np.exp(-1j * W * L)
    return D


def bilinear_shared_input(y, R, cl0, eps, Nfit, wt, Lt=300, rounds=8, d0=0):
    """x_{i,c} = d * t_{i,c}: one mono sequence d (from d0 to the end of the fit range) feeding every loop
    through short per-loop, per-channel tap sets t (Lt taps), alternating least squares from a naive start."""
    Nall = len(y); Ld = Nfit - PRE - d0; kn, kb = cluster_kernels(eps)
    CL = np.zeros((Nall, 2, 4))
    for c in (0, 1):
        for k, t_ in enumerate(CLUSTER[c]):
            idx = PRE + t_ + kn; m = (idx >= 0) & (idx < Nall); CL[idx[m], c, k] = kb[c][k][m]
    t = np.random.default_rng(0).normal(0, 0.01, (4, 2, Lt)); t[:, :, 0] = 1.0; d = np.zeros(Ld); cl = cl0.copy()
    F = lambda v: np.fft.rfft(v, NFFT); I = lambda V: np.fft.irfft(V, NFFT)
    def Hc(c): return sum(F(np.r_[t[i, c], np.zeros(NFFT - Lt)]) * R[i] for i in range(4))
    for it in range(rounds):
        H = [Hc(0), Hc(1)]
        tgt = np.concatenate([(y[:Nfit, c] - CL[:Nfit, c] @ cl[c]) * wt for c in (0, 1)])
        def mv(v):
            Dv = np.zeros(NFFT); Dv[PRE + d0:PRE + d0 + Ld] = v; Df = F(Dv); return np.concatenate([I(Df * H[c])[:Nfit] * wt for c in (0, 1)])
        def rmv(r):
            acc = np.zeros(NFFT)
            for c in (0, 1):
                z = np.zeros(NFFT); z[:Nfit] = r[c * Nfit:(c + 1) * Nfit] * wt; acc += I(F(z) * np.conj(H[c]))
            return acc[PRE + d0:PRE + d0 + Ld]
        d = lsqr(LinearOperator((2 * Nfit, Ld), matvec=mv, rmatvec=rmv, dtype=float), tgt, iter_lim=200, x0=d if it else None)[0]
        Dv = np.zeros(NFFT); Dv[PRE + d0:PRE + d0 + Ld] = d; Df = F(Dv); q = [I(Df * R[i])[:Nfit] for i in range(4)]
        for c in (0, 1):
            A = np.stack([np.r_[np.zeros(j), q[i][:Nfit - j]] for i in range(4) for j in range(Lt)] + [CL[:Nfit, c, k] for k in range(4)], 1) * wt[:, None]
            sol = np.linalg.lstsq(A, y[:Nfit, c] * wt, rcond=None)[0]; t[:, c, :] = sol[:4 * Lt].reshape(4, Lt); cl[c] = sol[4 * Lt:]
    Dv = np.zeros(NFFT); Dv[PRE + d0:PRE + d0 + Ld] = d; Df = F(Dv); out = np.zeros((Nall, 2))
    for c in (0, 1):
        out[:, c] = I(Df * Hc(c))[:Nall] + CL[:, c] @ cl[c]
    return out, Ld + 8 * Lt + 8


def diffuser_study(a):
    """Everything in docs/model_s_results.md §10. Writes <out>/diffuser.json."""
    import time
    T0 = time.time(); y = load(a.h); Nall = len(y); res = {}; n = np.arange(NFFT) - PRE; Y = _fwd(y)
    tj = json.load(open(os.path.join(a.out, "tighten.json"))); tl = json.load(open(os.path.join(a.out, "tail.json")))
    le = tj["loops_early"]; th5 = np.array(le["theta5"]); g5, hf5 = le["g5"], le["hf5"]
    tp = tl["tail_fit_both"]["params"]; vt = np.array([tp[k] for k in ["g1321", "g1389", "g1461", "g1587", "kappa_db_per_s", "p", "eps_ppm", "g460", "hf460"]])
    yb = band_limit(y)
    D_early = recursion_poly(G_param(th5, "jot"), [(G5_of(th5, g5, hf5, 460), 460)])
    D_whole = recursion_poly(G_perloop(vt), [(G5_of(np.r_[0, 0, 0, 0, vt[6]], vt[7], vt[8], 460), 460)])
    # 1. where does the input end? predictor from the measured state at T
    pred = {}
    for name, D in (("early", D_early), ("whole", D_whole)):
        for T in (0.30, 0.35, 0.40, 0.45, 0.50, 0.60):
            w = [(int(T * 1000), int(T * 1000) + 200), (400, 600), (600, 1000)] if T <= 0.40 else [(int(T * 1000), int(T * 1000) + 200)]
            nd = null_depth(yb, band_limit(predict_from_state(y, D, T)), w); fl = noise_floor_null(yb, w)
            pred[f"{name}_T{T:.2f}"] = {"null": nd, "excess": _excess(nd, fl)}
    res["state_predictor"] = pred
    print("predictor", {k: {w: np.round(x, 1).tolist() for w, x in v["null"].items()} for k, v in pred.items() if "0.40" in k or "0.45" in k}, flush=True)
    # 2. lag scan of the de-recirculated input, 150-17000 samples
    u = np.fft.irfft(Y * D_whole[:, None], NFFT, axis=0); v = u * ((n >= 0) & (n < 0.75 * FS))[:, None]
    F2 = lambda x: np.fft.rfft(x, 2 * NFFT, axis=0); I2 = lambda X: np.fft.irfft(X, 2 * NFFT, axis=0); scan = {}
    for win in ((0.22, 0.40), (0.30, 0.50)):
        m = ((n >= win[0] * FS) & (n < win[1] * FS)).astype(float)[:, None]; aa = v * m
        num = I2(np.conj(F2(v)) * F2(aa))[:17000]; den = I2(np.conj(F2(v ** 2)) * F2(m))[:17000]; red = (num ** 2 / np.maximum(den, 1e-30)) / np.sum(aa ** 2, 0)
        for lo, hi in ((150, 3500), (3500, 17000)):
            s = red[lo:hi].mean(1); idx = np.argsort(-s)[:8] + lo
            scan[f"{win[0]}-{win[1]}s_lags{lo}-{hi}"] = [(int(L), red[L].tolist(), (num[L] / den[L]).tolist()) for L in idx]
    res["lag_scan"] = scan
    # input envelope (dB over its noise-only level), 10 ms blocks
    seg = y[PRE + int(1.25 * FS):PRE + int(1.9 * FS)]; z = np.zeros((NFFT, 2)); z[PRE:PRE + len(seg)] = seg
    un = np.fft.irfft(np.fft.rfft(z, axis=0) * D_whole[:, None], NFFT, axis=0); sd = np.std(un[PRE + 3000:PRE + 14000], 0)
    res["input_envelope_db_over_noise"] = {"LR"[c]: [float(10 * np.log10(np.mean(u[(n >= t) & (n < t + 264), c] ** 2) / sd[c] ** 2)) for t in range(0, 13200, 264)] for c in (0, 1)}
    # 3. more 460-type stages / a free FIR 460 stage / an AR factor, each judged by the T=0.40 predictor
    Dl = recursion_poly(G_perloop(vt)); eps = vt[6] * 1e-6
    def stage(g, hf, L): return 1 - g * shelf(W, 1.0, hf, 0.0) * np.exp(-1j * W * L * (1 + eps))
    sel = (n >= 0.25 * FS) & (n < 0.75 * FS); wtt = tail_weight(n[sel])[:, None]; var = {}
    def p40(D): return null_depth(yb, band_limit(predict_from_state(y, D, 0.40)), [(400, 600), (600, 1000)])
    for lags in ([460], [460, 460], [460, 920], [460, 230], [460, 460, 230]):
        def rs(vv):
            A = Dl.copy()
            for k, L in enumerate(lags):
                A = A * stage(vv[2 * k], vv[2 * k + 1], L)
            return (np.fft.irfft(Y * A[:, None], NFFT, axis=0)[sel] * wtt).ravel()
        x0 = np.r_[[0.5, 0.9] * len(lags)]
        if len(lags) > 1:
            x0[2] = 0.3
        r = least_squares(rs, x0, bounds=([-.99, 0] * len(lags), [.99, 1.5] * len(lags)), diff_step=1e-4)
        A = Dl.copy()
        for k, L in enumerate(lags):
            A = A * stage(r.x[2 * k], r.x[2 * k + 1], L)
        var["stages " + "+".join(map(str, lags))] = {"cost": float(np.sum(r.fun ** 2)), "gains_hf": r.x.tolist(), "pred_T0.40": p40(A)}
    u_l = np.fft.irfft(Y * Dl[:, None], NFFT, axis=0); k = np.arange(-24, 25); L5 = int(round(460 * (1 + eps)))
    s2 = np.nonzero((n >= 0.25 * FS) & (n < 0.45 * FS))[0]
    X = np.stack([u_l[s2 - L5 - kk, :] for kk in k], -1).reshape(-1, len(k)); b = u_l[s2].reshape(-1); cfir = np.linalg.lstsq(X, b, rcond=None)[0]
    fr = np.array([100, 500, 1000, 2000, 4000, 8000]); G5f = np.exp(-1j * np.outer(W, k)) @ cfir
    var["460 as free 49-tap FIR"] = {"abs_G_at": dict(zip(map(str, fr), np.abs(np.exp(-1j * np.outer(2 * np.pi * fr / FS, k)) @ cfir).tolist())),
                                     "pred_T0.40": p40(Dl * (1 - G5f * np.exp(-1j * W * L5)))}
    u0 = np.fft.irfft(Y * D_whole[:, None], NFFT, axis=0); s3 = np.nonzero((n >= 0.25 * FS) & (n < 0.40 * FS))[0]
    for Lmax in (600, 1200):
        X = np.concatenate([np.stack([u0[s3 - kk, c] for kk in range(1, Lmax + 1)], 1) for c in (0, 1)]); b = np.concatenate([u0[s3, c] for c in (0, 1)])
        AtA = X.T @ X; ar = np.linalg.solve(AtA + 1e-3 * np.trace(AtA) / Lmax * np.eye(Lmax), X.T @ b); e = b - X @ ar
        top = np.argsort(-np.abs(ar))[:12] + 1
        var[f"AR({Lmax}) on the input, 0.25-0.40 s"] = {"explained": float(1 - e @ e / (b @ b)), "top_lags": [(int(kk), float(ar[kk - 1])) for kk in top],
                                                         "pred_T0.40": p40(D_whole * (1 - np.exp(-1j * np.outer(W, np.arange(1, Lmax + 1))) @ ar))}
    var["reference loops+460 (whole)"] = {"pred_T0.40": p40(D_whole)}
    res["recursion_variants"] = var
    print("variants done", f"{time.time() - T0:.0f}s", flush=True)
    # 4. shared long input? (a) spectral fine-structure sharing between the Cs excitations, (b) bilinear model
    Nfit = PRE + int(0.4 * FS); nn = np.arange(Nfit) - PRE; wt = tail_weight(np.maximum(nn, 0), t_cap=0.3); eps5 = th5[4] * 1e-6
    Gs5 = G_param(th5, "jot"); C5 = 1 / (1 - G5_of(th5, g5, hf5, 460) * np.exp(-1j * W * 460))
    lin = Linear(Gs5, Nfit, LF=1400, weight=wt, s0=S0, ext=300, eps=eps5, post=C5, post_f=True); _, _, xs = lin.split(lin.solve(y, iters=400))
    exc = [x for row in xs for x in row]; Nf = 1 << 14; f = np.fft.rfftfreq(Nf, 1 / FS); bd = (f > 200) & (f < 10000); share = {}
    for smh in (20, 100):
        S = []
        for x in exc:
            Lg = 20 * np.log10(np.abs(np.fft.rfft(x, Nf)) + 1e-12)[bd]; kk = int(smh / (f[1] - f[0])); S.append(Lg - np.convolve(Lg, np.ones(kk) / kk, "same"))
        C = np.corrcoef(S); same_loop = [C[i, i + 4] for i in range(4)]; cross = [C[i, j] for i in range(8) for j in range(i + 1, 8) if j != i + 4]
        rng = np.random.default_rng(0); Sr = []
        for x in exc:
            Lg = 20 * np.log10(np.abs(np.fft.rfft(rng.standard_normal(len(x)), Nf)) + 1e-12)[bd]; Sr.append(Lg - np.convolve(Lg, np.ones(kk) / kk, "same"))
        Cr = np.corrcoef(Sr)
        share[f"detrend_{smh}Hz"] = {"same_loop_LR": same_loop, "cross_loop_mean_abs": float(np.mean(np.abs(cross))),
                                     "random_control_mean_abs": float(np.mean(np.abs(Cr[np.triu_indices(8, 1)])))}
    res["fine_structure_sharing"] = share
    R5 = [C5 / (1 - G * np.exp(-1j * W * Pi)) for G, Pi in zip(Gs5, P)]
    outb, kpar = bilinear_shared_input(y, R5, np.array([CLUSTER_SIGN[c] for c in (0, 1)], float) * 1.12, eps5, Nfit, wt)
    ndb = null_depth(yb, band_limit(outb), HOLD_WIN); flb = noise_floor_null(yb, HOLD_WIN)
    res["bilinear"] = {"params": kpar, "null": ndb, "excess": _excess(ndb, flb)}
    print("bilinear", {w: np.round(x, 1).tolist() for w, x in ndb.items()}, f"{time.time() - T0:.0f}s", flush=True)
    with open(os.path.join(a.out, "diffuser.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print(f"total {time.time() - T0:.0f}s")



# ---------- all-pass diffuser search (docs/model_s_results.md §11) ----------
def allpass(Md, g, eps=0.0, w=W):
    z = np.exp(-1j * w * Md * (1 + eps)); return (-g + z) / (1 - g * z)


def allpass_study(a):
    """Everything in docs/model_s_results.md §11. Writes <out>/allpass.json.
    Score of a candidate diffuser: batch-OMP integer taps (in front of the diffuser, per loop and channel)
    fitted on 0-0.4 s; null in-sample (50-300 ms) and held out (400-600 ms), broadband."""
    import time
    T0 = time.time(); y = load(a.h); Nall = len(y); res = {}
    le = json.load(open(os.path.join(a.out, "tighten.json")))["loops_early"]; th5 = np.array(le["theta5"]); g5, hf5 = le["g5"], le["hf5"]; eps = th5[4] * 1e-6
    Gs = G_param(th5, "jot"); Nfit = PRE + int(0.4 * FS); wt = np.r_[tail_weight(np.maximum(np.arange(Nfit) - PRE, 0), t_cap=0.3), np.zeros(Nall - Nfit)]
    comb460 = 1 / (1 - G5_of(th5, g5, hf5, 460) * np.exp(-1j * W * 460)); ap = lambda Md, g: allpass(Md, g, eps)
    def score(post, ch=1, taps=(200,)):
        r = omp_curve(np.stack([y[:, ch]] * 2, 1), Gs, post, eps, Nfit, wt, checkpoints=taps)
        return [(d["taps"], d["null"]["50-300"], d["null"]["400-600"]) for d in r["L"]]
    log = lambda *s: print(*s, f"{time.time() - T0:.0f}s", flush=True)
    # stage 2 on top of the 460 stage (R)
    st = {"comb460": score(comb460), "ap460": score(ap(460, g5))}
    for Md in (20, 30, 52, 72, 100, 144, 230, 690, 920):
        for g in (-0.6, 0.6):
            st[f"ap460*ap{Md}({g})"] = score(ap(460, g5) * ap(Md, g))
    res["stage2_R"] = st; log("stage 2")
    res["grid_20_72_R"] = {f"{g1},{g2}": score(ap(460, g5) * ap(20, g1) * ap(72, g2))
                           for g1, g2 in ((-0.5, -0.5), (-0.6, -0.6), (-0.6, -0.7), (-0.7, -0.6), (-0.7, -0.7), (-0.75, -0.75), (-0.8, -0.8), (-0.9, -0.9))}
    log("grid")
    sharp = {}
    for Md in (18, 19, 21, 22):
        sharp[f"ap{Md} instead of 20"] = score(ap(460, g5) * ap(Md, -0.7) * ap(72, -0.7))
    for Md in (70, 71, 73, 74):
        sharp[f"ap{Md} instead of 72"] = score(ap(460, g5) * ap(20, -0.7) * ap(Md, -0.7))
    for Md in (458, 459, 461, 462):
        sharp[f"ap{Md} instead of 460"] = score(ap(Md, g5) * ap(20, -0.7) * ap(72, -0.7))
    sharp["comb460 instead of ap460"] = score(comb460 * ap(20, -0.7) * ap(72, -0.7))
    sharp["no 460 stage"] = score(ap(20, -0.7) * ap(72, -0.7))
    res["delay_sharpness_R"] = sharp; log("sharpness")
    B = comb460 * ap(20, -0.7) * ap(72, -0.7)
    res["stage3_R"] = {f"ap{Md}({g})": score(B * ap(Md, g)) for Md in (3, 115, 153, 229, 230, 231, 383) for g in (-0.6,)}
    res["stage3_R"].update({f"ap230({g})": score(B * ap(230, g)) for g in (-0.45, -0.7, -0.8, 0.6)})
    B3 = B * ap(230, -0.7)
    res["stage4_R"] = {f"ap{Md}({g})": score(B3 * ap(Md, g)) for Md in (3, 7, 30, 52, 115, 144, 153) for g in (-0.6, 0.6)}
    res["common_gain_R"] = {f"{g}": score(comb460 * ap(20, g) * ap(72, g) * ap(230, g)) for g in (-0.65, -0.6875, -0.7, -0.72, -0.75)}
    log("stages 3-4")
    # L confirmation (chain chosen on R)
    res["L_check"] = {"comb460 only": score(comb460, ch=0), "chain20,72": score(B, ch=0), "chain20,72,230": score(B3, ch=0)}
    log("L check")
    # tap-count curves, both channels, chain vs none
    yb = band_limit(y); flb = noise_floor_null(yb, HOLD_WIN); curves = {}
    for name, post in (("comb460 only", comb460), ("comb460*ap20*ap72*ap230 (-0.7)", B3)):
        r = omp_curve(y, Gs, post, eps, Nfit, wt, checkpoints=(100, 400, 1600, 3200))
        curves[name] = r
    res["curves"] = curves; res["noise_floor_dev"] = noise_floor_null(y, HOLD_WIN)
    with open(os.path.join(a.out, "allpass.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    log("done")



# ---------- group-sparse injection taps (docs/model_s_results.md §12) ----------
def group_omp(y, Gs, post, eps, Nfit, wt, d_max, checkpoints=(4, 8, 16, 32, 64, 128), direct=True, d_min=0):
    """Tap delays d (device samples after cluster tap 0) shared by every slot: 4 loops (+ the direct path)
    x 2 channels, each slot with its own gain. A slot's atom is the response of post*loop (or post alone)
    delayed by d*(1+eps). Greedy group OMP on 0-Nfit (weighted); exact fractional-delay columns in the refit."""
    Nall = len(y); kn, kb = cluster_kernels(eps)
    Hs = [post / (1 - G * np.exp(-1j * W * Pi)) for G, Pi in zip(Gs, P)] + ([post] if direct else [])
    ns = len(Hs)
    hs_int = [np.fft.irfft(H, NFFT) for H in Hs]
    def col(s, d):
        return np.fft.irfft(Hs[s] * np.exp(-1j * W * d * (1 + eps)), NFFT)[np.arange(Nall) - PRE] if False else \
            np.roll(np.fft.irfft(Hs[s] * np.exp(-1j * W * d * (1 + eps)), NFFT), PRE)[:Nall]
    def ccol(c, k):
        v = np.zeros(Nall); idx = PRE + CLUSTER[c][k] + kn; m = (idx >= 0) & (idx < Nall); v[idx[m]] = kb[c][k][m]; return v
    w2 = np.zeros(NFFT); w2[:Nfit] = wt[:Nfit] ** 2; W2 = np.fft.rfft(w2)
    Hspec = [np.fft.rfft(np.roll(h, PRE)) for h in hs_int]
    nrm = [np.sqrt(np.maximum(np.fft.irfft(W2 * np.conj(np.fft.rfft(np.roll(h, PRE) ** 2)), NFFT)[d_min:d_max], 1e-30)) for h in hs_int]
    sel = []; cols = {c: [ccol(c, k) for k in range(4)] for c in (0, 1)}; out = []; models = {}
    r = y.copy()
    while len(sel) < checkpoints[-1]:
        score = np.zeros(d_max - d_min)
        for c in (0, 1):
            rw = np.zeros(NFFT); rw[:Nfit] = r[:Nfit, c] * wt[:Nfit] ** 2; Rr = np.fft.rfft(rw)
            for s in range(ns):
                score += (np.fft.irfft(Rr * np.conj(Hspec[s]), NFFT)[d_min:d_max] / nrm[s]) ** 2
        for d in sel:
            score[max(0, d - d_min - 1):d - d_min + 2] = 0
        d = int(np.argmax(score)) + d_min; sel.append(d)
        model = np.zeros_like(y); gains = {}
        for c in (0, 1):
            cols[c] += [col(s, d) for s in range(ns)]
            A = np.stack(cols[c], 1); coef = np.linalg.lstsq(A[:Nfit] * wt[:Nfit, None], y[:Nfit, c] * wt[:Nfit], rcond=None)[0]
            model[:, c] = A @ coef; gains[c] = coef[4:].reshape(-1, ns)
        r = y - model
        if len(sel) in checkpoints:
            models[len(sel)] = (list(sel), {c: gains[c].copy() for c in (0, 1)}, model.copy())
    for k, (ds, g, model) in models.items():
        out.append({"groups": k, "params": k * ns * 2 + 8, "delays": ds, "gains": {"LR"[c]: g[c].tolist() for c in (0, 1)},
                    "null": null_depth(y, model, HOLD_WIN), "null_band": null_depth(band_limit(y), band_limit(model), HOLD_WIN)})
    return out


def loop_shared_omp(y, Gs, post, eps, Nfit, wt, d_max=2000, checkpoints=(50, 100, 200, 400)):
    """Groups = (loop or direct path, delay) shared by L and R (one gain per channel)."""
    Nall = len(y); kn, kb = cluster_kernels(eps)
    Hs = [post / (1 - G * np.exp(-1j * W * Pi)) for G, Pi in zip(Gs, P)] + [post]; ns = len(Hs)
    hs = [np.roll(np.fft.irfft(H, NFFT), PRE) for H in Hs]; Hspec = [np.fft.rfft(h) for h in hs]
    w2 = np.zeros(NFFT); w2[:Nfit] = wt[:Nfit] ** 2; W2 = np.fft.rfft(w2)
    nrm = [np.sqrt(np.maximum(np.fft.irfft(W2 * np.conj(np.fft.rfft(h ** 2)), NFFT)[:d_max], 1e-30)) for h in hs]
    col = lambda s, d: np.roll(np.fft.irfft(Hs[s] * np.exp(-1j * W * d * (1 + eps)), NFFT), PRE)[:Nall]
    def ccol(c, k):
        v = np.zeros(Nall); idx = PRE + CLUSTER[c][k] + kn; m = (idx >= 0) & (idx < Nall); v[idx[m]] = kb[c][k][m]; return v
    cols = {c: [ccol(c, k) for k in range(4)] for c in (0, 1)}; sel = []; r = y.copy(); out = []
    while len(sel) < checkpoints[-1]:
        sc = np.zeros((ns, d_max))
        for c in (0, 1):
            rw = np.zeros(NFFT); rw[:Nfit] = r[:Nfit, c] * wt[:Nfit] ** 2; Rr = np.fft.rfft(rw)
            for s in range(ns):
                sc[s] += (np.fft.irfft(Rr * np.conj(Hspec[s]), NFFT)[:d_max] / nrm[s]) ** 2
        for s, d in sel:
            sc[s, d] = 0
        s, d = np.unravel_index(np.argmax(sc), sc.shape); sel.append((int(s), int(d))); model = np.zeros_like(y)
        for c in (0, 1):
            cols[c].append(col(s, d)); A = np.stack(cols[c], 1)
            model[:, c] = A @ np.linalg.lstsq(A[:Nfit] * wt[:Nfit, None], y[:Nfit, c] * wt[:Nfit], rcond=None)[0]
        r = y - model
        if len(sel) in checkpoints:
            out.append({"groups": len(sel), "params": 2 * len(sel) + 8, "null": null_depth(y, model, HOLD_WIN),
                        "null_band": null_depth(band_limit(y), band_limit(model), HOLD_WIN)})
    return out


def group_study(a):
    """Everything in docs/model_s_results.md §12. Writes <out>/groups.json."""
    import time
    T0 = time.time(); y = load(a.h); Nall = len(y); res = {}
    le = json.load(open(os.path.join(a.out, "tighten.json")))["loops_early"]; th5 = np.array(le["theta5"]); g5, hf5 = le["g5"], le["hf5"]; eps = th5[4] * 1e-6
    Gs = G_param(th5, "jot"); Nfit = PRE + int(0.4 * FS); wt = np.r_[tail_weight(np.maximum(np.arange(Nfit) - PRE, 0), t_cap=0.3), np.zeros(Nall - Nfit)]
    post = (1 / (1 - G5_of(th5, g5, hf5, 460) * np.exp(-1j * W * 460))) * allpass(20, -.7, eps) * allpass(72, -.7, eps) * allpass(230, -.7, eps)
    res["noise_floor"] = noise_floor_null(y, HOLD_WIN); res["noise_floor_band"] = noise_floor_null(band_limit(y), HOLD_WIN)
    g = group_omp(y, Gs, post, eps, Nfit, wt, 2000, checkpoints=(4, 8, 16, 32, 64))
    for d in g[-1:]:
        ds = np.array(d["delays"]); stats = {}
        for c in "LR":
            gg = np.array(d["gains"][c]); e = 20 * np.log10(np.sqrt(np.mean(gg ** 2, 1)) + 1e-12)
            stats[c] = {"gain_vs_delay_db_per_s": float(np.polyfit(ds / FS, e, 1)[0]), "corr": float(np.corrcoef(ds, e)[0, 1]),
                        "sign_agreement_across_loops": float(np.mean(np.abs(np.sign(gg[:, :4]).mean(1))))}
        res["group_64_gain_stats"] = stats
    res["shared_all_slots"] = [{k: v for k, v in d.items() if k != "gains"} for d in g]
    print("shared all", [(d["groups"], np.round(d["null"]["400-600"], 1).tolist()) for d in g], f"{time.time() - T0:.0f}s", flush=True)
    res["shared_LR_per_loop"] = loop_shared_omp(y, Gs, post, eps, Nfit, wt)
    print("shared LR", [(d["groups"], np.round(d["null"]["400-600"], 1).tolist()) for d in res["shared_LR_per_loop"]], f"{time.time() - T0:.0f}s", flush=True)
    res["independent"] = {}
    for ext in (300, 8500):
        r = omp_curve(y, Gs, post, eps, Nfit, wt, checkpoints=(414, 1600), ext=ext)
        res["independent"][f"ext{ext}"] = r
    print("independent", {k: {c: [(d["taps"], round(d["null"]["400-600"], 1)) for d in v[c]] for c in "LR"} for k, v in res["independent"].items()}, f"{time.time() - T0:.0f}s", flush=True)
    with open(os.path.join(a.out, "groups.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)



# ---------- "stochastic control" measurement simulation (docs/model_s_results.md §13) ----------
def ess(T, f1=20.0, f2=12000.0, fs=FS):
    t = np.arange(int(T * fs)) / fs; L = T / np.log(f2 / f1)
    x = np.sin(2 * np.pi * f1 * L * (np.exp(t / L) - 1)); fade = int(0.01 * fs)
    x[:fade] *= np.hanning(2 * fade)[:fade]; x[-fade:] *= np.hanning(2 * fade)[fade:]
    return x


class Synth246:
    """Synthetic program S from the fitted structure, with sparse injection taps (common to L and R, per loop)
    and sparse per-channel read-out taps (per loop), either fixed or re-drawn / jittered every `period` s.
    Signal path: x -> [pre delay PRE] -> diffuser (comb460*AP20*AP72*AP230) -> injection taps_i(t) -> loop i
    -> read-out taps_{i,c}(t) -> sum + cluster (fixed taps on x)."""

    def __init__(self, Gs, post, n_inj=8, readout_lags=(0, 26, 332), seed=0):
        self.Gs, self.post, self.n_inj, self.rl = Gs, post, n_inj, readout_lags
        self.cluster_amp = 8.0
        self.rng = np.random.default_rng(seed)

    def draw_inj(self):
        return [(S0 + self.rng.integers(0, Pi, self.n_inj), self.rng.normal(0, 1, self.n_inj)) for Pi in P]

    def draw_ro(self):
        return [[self.rng.normal(0, 1, len(self.rl)) for _ in P] for _ in (0, 1)]

    def run(self, x, period=0.2, inj_J=0, inj_gd=0.0, inj_redraw=False, roL_gd=0.0, tail_s=2.5):
        """Injection taps: every `period` s positions jitter by +-inj_J samples and gains by (1 + inj_gd*N(0,1))
        about the base set (or are re-drawn completely if inj_redraw). L read-out gains: (1 + roL_gd*N(0,1)) per period."""
        n = len(x) + int(tail_s * FS); Nf = 1 << int(np.ceil(np.log2(n + 1)))
        w = 2 * np.pi * np.fft.rfftfreq(Nf); X = np.fft.rfft(x, Nf) * np.exp(-1j * w * PRE)
        d = np.fft.irfft(X * self.post(w), Nf)[:n]
        tv = inj_redraw or inj_J or inj_gd or roL_gd
        blocks = [(s, min(n, s + int(period * FS))) for s in range(0, n, int(period * FS))] if tv else [(0, n)]
        inj0 = self.draw_inj(); ro0 = self.draw_ro(); y = np.zeros((n, 2)); cfg_inj, cfg_ro = [], []
        for k, b in enumerate(blocks):
            if inj_redraw and k:
                cfg_inj.append(self.draw_inj())
            else:
                cfg_inj.append([(np.clip(dl + (self.rng.integers(-inj_J, inj_J + 1, len(dl)) if inj_J else 0), S0, None),
                                 g * (1 + inj_gd * self.rng.normal(0, 1, len(g)))) for dl, g in inj0])
            cfg_ro.append([[g * (1 + roL_gd * self.rng.normal(0, 1, len(g))) for g in ro0[0]], ro0[1]])
        for i, (G, Pi) in enumerate(zip(self.Gs, P)):
            z = np.zeros(n)
            for (s, e), cfg in zip(blocks, cfg_inj):
                dl, g = cfg[i]
                for dd, gg in zip(dl, g):
                    lo = max(s, dd)
                    if lo < e:
                        z[lo:e] += gg * d[lo - dd:e - dd]
            yi = np.fft.irfft(np.fft.rfft(z, Nf) / (1 - G(w) * np.exp(-1j * w * Pi)), Nf)[:n]
            for c in (0, 1):
                for (s, e), ro in zip(blocks, cfg_ro):
                    for lag, gg in zip(self.rl, ro[c][i]):
                        lo = max(s, lag)
                        if lo < e:
                            y[lo:e, c] += gg * yi[lo - lag:e - lag]
        xd = np.fft.irfft(X, Nf)[:n]
        for c in (0, 1):
            for t, sg in zip(CLUSTER[c], CLUSTER_SIGN[c]):
                y[t:, c] += sg * self.cluster_amp * xd[:n - t]
        return y


def deconvolve(y, x, f1=20.0, f2=12000.0, n_out=None):
    """Regularised spectral division (what an ESS inverse filter does in-band). Returns the circular IR."""
    Nf = 1 << int(np.ceil(np.log2(len(y) + 1))); f = np.fft.rfftfreq(Nf, 1 / FS)
    X = np.fft.rfft(x, Nf); Y = np.fft.rfft(y, Nf, axis=0)
    win = np.clip((f - f1 / 2) / (f1 / 2), 0, 1) * np.where(f < f2 * 0.95, 1.0, np.where(f < f2, 0.5 * (1 + np.cos(np.pi * (f - 0.95 * f2) / (0.05 * f2))), 0.0))
    eps = 1e-6 * np.max(np.abs(X)) ** 2
    H = Y * (np.conj(X) * win / (np.abs(X) ** 2 + eps))[:, None]
    return np.fft.irfft(H, Nf, axis=0)


def gentle_band(h, lo=(30, 80), hi=(9000, 11000)):
    """Band-limit with wide raised-cosine edges (short ringing), so the state predictor is not fed long band-edge tails."""
    Nf = 1 << int(np.ceil(np.log2(2 * len(h)))); f = np.fft.rfftfreq(Nf, 1 / FS)
    H = 0.5 - 0.5 * np.cos(np.pi * np.clip((f - lo[0]) / (lo[1] - lo[0]), 0, 1)); H *= 0.5 - 0.5 * np.cos(np.pi * np.clip((hi[1] - f) / (hi[1] - hi[0]), 0, 1))
    return np.fft.irfft(np.fft.rfft(h, Nf, axis=0) * (H[:, None] if h.ndim == 2 else H), Nf, axis=0)[:len(h)]


def sim_diagnostics(h, Gs_W, post_W, D, B, eps=0.0, ns_db=None, seed=1, omp_taps=(103, 414), normalise=True):
    """Diagnostics used for both the CLA IR and simulated IRs (tap 0 at index PRE):
    OMP with band-limited atoms (post*B), the state predictor on a gently band-limited IR, cross-band envelopes."""
    h = h[:PRE + int(1.9 * FS)].copy()
    if normalise:
        h /= np.abs(h[:PRE + 400]).max()
    if ns_db is not None:
        h += np.random.default_rng(seed).normal(0, 10 ** (ns_db / 20), h.shape)
    Nall = len(h); Nfit = PRE + int(0.4 * FS)
    wt = np.r_[tail_weight(np.maximum(np.arange(Nfit) - PRE, 0), t_cap=0.3), np.zeros(Nall - Nfit)]
    out = {}
    r = omp_curve(h, Gs_W, post_W * B, eps, Nfit, wt, checkpoints=omp_taps)
    out["omp"] = {c: [(d["taps"], d["null"]["50-300"], d["null"]["400-600"]) for d in r[c]] for c in "LR"}
    g = gentle_band(h); out["predictor"] = {}
    for T in (0.35, 0.40, 0.45):
        w = (int(T * 1000), int(T * 1000) + 200); k = f"{w[0]}-{w[1]}"
        nd = null_depth(g, predict_from_state(g, D, T), [w]); fl = noise_floor_null(g, [w])
        out["predictor"][f"T{T:.2f}"] = {"null": nd[k], "excess": [nd[k][c] - fl[k][c] for c in (0, 1)]}
    out["xband"] = cross_band(g, D)
    return out


def cross_band(h, D, win=(0.05, 0.40), centres=(250, 500, 1000, 2000, 4000, 8000)):
    """Envelope correlation of the de-recirculated input between 1/3-octave bands (per channel),
    as a function of band separation in octaves."""
    n = np.arange(NFFT) - PRE; hh = np.zeros((NFFT, 2)); hh[:len(h)] = h
    u = np.fft.irfft(np.fft.rfft(hh, axis=0) * D[:, None], NFFT, axis=0)
    m = (n >= win[0] * FS) & (n < win[1] * FS); env = {}
    from scipy.signal import hilbert
    for fc in centres:
        sos = butter(4, [fc / 2 ** (1 / 6), fc * 2 ** (1 / 6)], btype="band", fs=FS, output="sos")
        b = sosfilt(sos, u, axis=0); e = np.abs(hilbert(b[m], axis=0))
        k = int(0.005 * FS); env[fc] = np.stack([np.convolve(e[:, c], np.ones(k) / k, "same") for c in (0, 1)], 1)
    res = {}
    for i, f1 in enumerate(centres):
        for f2 in centres[i + 1:]:
            sep = float(np.log2(f2 / f1)); res.setdefault(sep, []).append([float(np.corrcoef(np.log(env[f1][:, c] + 1e-15), np.log(env[f2][:, c] + 1e-15))[0, 1]) for c in (0, 1)])
    return {f"{k:.0f}oct": np.mean(v, 0).tolist() for k, v in sorted(res.items())}


def stochastic_study(a):
    """Everything in docs/model_s_results.md §13. Writes <out>/stochastic.json."""
    import time, soundfile as sf
    from deconv_grid import locate_first_tap
    T0 = time.time(); res = {}
    le = json.load(open(os.path.join(a.out, "tighten.json")))["loops_early"]; th5 = np.array(le["theta5"]); g5, hf5 = le["g5"], le["hf5"]; eps5 = th5[4] * 1e-6
    tl = json.load(open(os.path.join(a.out, "tail.json")))["tail_fit_both"]["params"]
    vt = np.array([tl[k] for k in ["g1321", "g1389", "g1461", "g1587", "kappa_db_per_s", "p", "eps_ppm", "g460", "hf460"]])
    # 1. pre-onset energy and distortion images in the WAV
    x96, fs96 = sf.read(a.wav, always_2d=True); i0 = locate_first_tap(x96[:, 1], fs96); pk = np.abs(x96[i0 - 5:i0 + 5]).max(0)
    tmpl = x96[i0 - 50:i0 + int(0.015 * fs96)]; pre = {}
    for c in (0, 1):
        t_ = tmpl[:, c]; N = 1 << int(np.ceil(np.log2(len(x96) + len(t_))))
        cc = np.fft.irfft(np.fft.rfft(x96[:, c], N) * np.conj(np.fft.rfft(t_, N)), N)[:len(x96)] / np.sqrt(np.sum(t_ ** 2))
        e = np.sqrt(np.convolve(x96[:, c] ** 2, np.ones(len(t_)), "full")[len(t_) - 1:len(t_) - 1 + len(x96)] + 1e-30); ncc = np.abs(cc / e)
        pre["LR"[c]] = {"pre_onset_rms_db_re_tap": float(20 * np.log10(np.sqrt(np.mean(x96[:i0 - 20, c] ** 2)) / pk[c])),
                        "last_0.3s_rms_db_re_tap": float(20 * np.log10(np.sqrt(np.mean(x96[-int(0.3 * fs96):, c] ** 2)) / pk[c])),
                        "best_cluster_copy_pre": float(ncc[:i0 - 200].max()), "best_cluster_copy_end": float(ncc[int(1.0 * fs96):len(x96) - len(t_)].max()),
                        "median_match": float(np.median(ncc[:i0 - 200]))}
    res["preonset"] = {"first_tap_ms": i0 / fs96 * 1000, "channels": pre}
    # 2. synthetic 246 (integer loops) and its measurement through an ESS
    th0 = th5.copy(); th0[4] = 0.0
    Gw = [(lambda w, k=k: G_param(th0, "jot", w)[k]) for k in range(4)]
    def post(w):
        z = lambda d_: np.exp(-1j * w * d_); ap_ = lambda d_, g: (-g + z(d_)) / (1 - g * z(d_))
        return 1 / (1 - g5 * shelf(w, 1.0, hf5, 0.0) * z(460)) * ap_(20, -.7) * ap_(72, -.7) * ap_(230, -.7)
    Gs_W = G_param(th0, "jot"); post_W = post(W); D = recursion_poly(Gs_W, [(g5 * shelf(W, 1.0, hf5, 0.0), 460)])
    f = np.fft.rfftfreq(NFFT, 1 / FS)
    B_sim = np.clip((f - 10) / 10, 0, 1) * np.where(f < 11400, 1.0, np.where(f < 12000, 0.5 * (1 + np.cos(np.pi * (f - 11400) / 600)), 0.0))
    B_cla = np.where(f < 11000, 1.0, np.where(f < 12300, 0.5 * (1 + np.cos(np.pi * (f - 11000) / 1300)), 0.0))
    y = load(a.h)
    postC = (1 / (1 - G5_of(th5, g5, hf5, 460) * np.exp(-1j * W * 460))) * allpass(20, -.7, eps5) * allpass(72, -.7, eps5) * allpass(230, -.7, eps5)
    D_w = recursion_poly(G_perloop(vt), [(G5_of(np.r_[0, 0, 0, 0, vt[6]], vt[7], vt[8], 460), 460)])
    res["CLA"] = sim_diagnostics(y, G_param(th5, "jot"), postC, D_w, B_cla, eps=eps5, normalise=False)
    res["CLA_atoms_without_bandlimit"] = sim_diagnostics(y, G_param(th5, "jot"), postC, D_w, 1.0, eps=eps5, normalise=False)["omp"]
    print("CLA", res["CLA"]["omp"], {k: np.round(v["excess"], 1).tolist() for k, v in res["CLA"]["predictor"].items()}, f"{time.time() - T0:.0f}s", flush=True)
    variants = [("LTI", {}, 10.0), ("inj gain 3% / 0.2 s", dict(inj_gd=0.03), 10.0), ("inj gain 10% / 0.2 s", dict(inj_gd=0.1), 10.0),
                ("inj gain 30% / 0.2 s", dict(inj_gd=0.3), 10.0), ("inj gain 10% / 0.05 s", dict(inj_gd=0.1, period=0.05), 10.0),
                ("inj gain 10% / 1 s", dict(inj_gd=0.1, period=1.0), 10.0), ("inj positions +-1 / 0.2 s", dict(inj_J=1), 10.0),
                ("inj full re-draw / 0.2 s", dict(inj_redraw=True), 10.0), ("inj full re-draw / 1 s", dict(inj_redraw=True, period=1.0), 10.0),
                ("L read-out 1% / 0.2 s", dict(roL_gd=0.01), 10.0), ("L read-out 2% / 0.2 s", dict(roL_gd=0.02), 10.0),
                ("L read-out 3% / 0.2 s", dict(roL_gd=0.03), 10.0), ("L read-out 10% / 0.2 s", dict(roL_gd=0.1), 10.0),
                ("L read-out 2% / 1 s", dict(roL_gd=0.02, period=1.0), 10.0), ("L read-out 2% / 0.2 s, 30 s sweep", dict(roL_gd=0.02), 30.0),
                ("inj 3% + L read-out 2% / 0.2 s", dict(inj_gd=0.03, roL_gd=0.02), 10.0), ("LTI, 30 s sweep", {}, 30.0)]
    sims = {}; xs = {}
    for name, kw, Tsw in variants:
        if Tsw not in xs:
            xs[Tsw] = ess(Tsw)
        S = Synth246(Gw, post, seed=0); h = deconvolve(S.run(xs[Tsw], **kw), xs[Tsw])
        sims[name] = sim_diagnostics(h, Gs_W, post_W, D, B_sim, ns_db=-105)
        if name == "LTI":
            sims[name]["omp_atoms_without_bandlimit"] = sim_diagnostics(h, Gs_W, post_W, D, 1.0, ns_db=-105)["omp"]
        print(name, sims[name]["omp"], {k: np.round(v["excess"], 1).tolist() for k, v in sims[name]["predictor"].items()}, f"{time.time() - T0:.0f}s", flush=True)
    res["sims"] = sims
    with open(os.path.join(a.out, "stochastic.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
