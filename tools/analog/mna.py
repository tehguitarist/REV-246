"""Small-signal AC analysis of the netlists in netlists.py by modified nodal analysis.

Op-amps are ideal nullors (v+ = v-) by default, or with a single-pole open loop gain
A(s) = 2*pi*GBW/s (TL07x typ. 3 MHz; DC gain 1e5) when gbw is given.
"""
import numpy as np


def response(net, f, gbw=None):
    elems, amps = net["elements"], net["opamps"]
    nodes = sorted({n for e in elems for n in e[1:3]} | {n for a in amps for n in a[:3]} - {"gnd"})
    nodes = [n for n in nodes if n != "gnd"]
    idx = {n: i for i, n in enumerate(nodes)}
    n_n, n_a = len(nodes), len(amps)
    out = net["out"]
    H = np.zeros(len(f), complex)
    for fi, fr in enumerate(f):
        s = 2j * np.pi * fr
        N = n_n + n_a + 1                       # nodes, op-amp output currents, source current
        A = np.zeros((N, N), complex)
        b = np.zeros(N, complex)
        for kind, a_, b_, val, _ in elems:
            y = 1 / val if kind == "R" else s * val
            for p_, q_ in ((a_, b_), (b_, a_)):
                if p_ != "gnd":
                    A[idx[p_], idx[p_]] += y
                    if q_ != "gnd":
                        A[idx[p_], idx[q_]] -= y
        for j, (pl, mi, ou, _) in enumerate(amps):
            r = n_n + j
            A[idx[ou], r] += 1                  # op-amp output injects current r into its node
            if gbw is None:
                if pl != "gnd": A[r, idx[pl]] += 1
                if mi != "gnd": A[r, idx[mi]] -= 1
            else:
                gain = 2 * np.pi * gbw / (s + 2 * np.pi * gbw / 1e5)
                if pl != "gnd": A[r, idx[pl]] += gain
                if mi != "gnd": A[r, idx[mi]] -= gain
                A[r, idx[ou]] -= 1
        A[idx["in"], n_n + n_a] += 1            # source current into node 'in'
        A[n_n + n_a, idx["in"]] = 1
        b[n_n + n_a] = 1.0
        x = np.linalg.solve(A, b)
        H[fi] = x[idx[out]]
    return H
