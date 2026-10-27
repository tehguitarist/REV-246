"""Shared loaders for the chain3 tail-error test (analysis/nl_lowband30/cand_kernel). Read-only on everything else."""
import numpy as np, json, sys
sys.path.insert(0, '/Users/leighpierce/Git/246/tools')
ROOT = '/Users/leighpierce/Git/246'
FS = 96000.0
z = np.load(ROOT + '/analysis/library/tomkin96/kernels/chain3.npz')
g96, LEAD, HP = z['g'], int(z['lead']), [float(v) for v in z['hp_hz']]
def G_of(f, ch):
    """chain3 complex response at f (Hz), direct DFT (FIR incl. lead, times analytic HPs), exactly Kernel.response."""
    f = np.asarray(f, float); n = np.arange(len(g96))
    E = np.exp(-2j * np.pi * np.outer(f, n) / FS)
    G = E @ g96[:, ch]
    H = np.exp(2j * np.pi * f * LEAD / FS)
    for fc in HP: H = H * (1j * f / fc) / (1 + 1j * f / fc)
    return G * H
def J():
    j = json.load(open(ROOT + '/analysis/nl_lowband30/lowband_joint_T.json')); f = np.array(j['f_hz'])
    out = {'f': f}
    for p in ('joint', 'N', 'S'):
        for c in (0, 1):
            out[(p, c)] = 10 ** (np.array(j[f'{p}_db{c}']) / 20) * np.exp(1j * np.radians(np.array(j[f'{p}_deg{c}'])))
    return out
