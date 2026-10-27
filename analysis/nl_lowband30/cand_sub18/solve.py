from fitlib import *
import pickle
def pool(D, keys, v, ch, scale_w=None):
    G = np.zeros((NC, NC)); b = np.zeros(NC); yy = 0.0; n = 0
    for k in keys:
        gg = D[k]["g"]; g = (gg["ls"] if "ls" in gg else gg).get(v)
        if g is None: continue
        Gk, bk, y2, nk, sc = g[ch]
        G += Gk; b += bk; yy += y2; n += nk
    return G, b, yy, n

def solve(G, b, lam=1e-3, free=None):
    """ridge toward 0 (= the base model); free: boolean mask of columns allowed (others 0)."""
    idx = np.arange(NC) if free is None else np.nonzero(free)[0]
    Gs = G[np.ix_(idx, idx)]; reg = lam * np.mean(np.diag(Gs)) * np.eye(len(idx))
    c = np.zeros(NC); c[idx] = np.linalg.solve(Gs + reg, b[idx]); return c

def Tof(c):
    return 1 + c[0::2] + 1j * c[1::2]

def gain_db(G, b, yy, c):
    """energy reduction of the (scaled) residual from correction c, in dB re yy (small numbers: whole-spectrum residual)"""
    red = 2 * c @ b - c @ G @ c
    return 10 * np.log10(max(yy - red, 1e-30) / yy)

def sel(D, **kw):
    out = []
    for k, d in D.items():
        r = d["r"]
        if "sess" in kw and d["sess"] not in kw["sess"]: continue
        if "span" in kw and k[0] not in kw["span"]: continue
        if "late0" in kw and (k[2] == 0) != kw["late0"]: continue
        if "late_max" in kw and k[2] > kw["late_max"]: continue
        out.append(k)
    return out

def fmt(T, idxs=None):
    idxs = range(len(KNOTS)) if idxs is None else idxs
    return "  ".join(f"{KNOTS[i]:.0f}:{20*np.log10(abs(T[i])):+5.1f}dB/{np.degrees(np.angle(T[i])):+4.0f}" for i in idxs)
