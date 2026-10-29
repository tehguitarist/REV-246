"""Program E's moving-read reference IRs for the plugin's exact test (plugin/tests/dsp/ETests.cpp).

    .venv/bin/python tools/memguard.py --limit-gb 2.5 -- .venv/bin/python -m tools.plugin_codegen.e_moving_ref

Writes analysis/e250/ref/e_moving_rt<RT>.wav: tools/e250/e_sim.render with e_sim.ADOPTED (the triangle), the moving reads
on, seed 0 (the realisation plugin/src/dsp/generated/e_spec.h carries), L's output kernel only (chain 'device'), AMPL 80,
LF / HF from static_all_246.json. Same format and length as the static references (e_spec.py): float32, device rate
(WAV header 26374), channels (L, R), int(secs * FS) samples with secs = min(3 RT + 0.5, 4.0), R's direct tap at k0 = 64.
The impulse enters render at network time k0 (as ECore delays its input by k0), so the reads' time t = n / FS is ECore's
sample count since reset(); render_ir's own alignment (impulse at 1000) would start the reads 936 samples later.
"""
import os
import sys

import numpy as np
import soundfile as sf

from .common import ROOT

sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "e250"))
import e_sim as S  # noqa: E402

RTS = (1.0, 4.5)
KOFF = 16                                   # the output kernel's offset (spec output_kernel.taps_offset)


def render(rt):
    S.static_params()
    prm = dict(S.DEF) | S.ADOPTED | dict(f=S.LFO_RATE["time"], ampl_scale=1.0)
    k = np.zeros(KOFF + 2); k[KOFF:] = S.L_KERNEL
    n = int(min(3 * rt + 0.5, 4.0) * S.FS)
    x = np.zeros(n + KOFF); x[S.K0] = 1.0
    y, _ = S.render(x, rt, prm, [k, k.copy()], np.random.default_rng(0))
    return y[KOFF:KOFF + n]                 # out[n] = 0.77 y[n - k0] + 0.23 y[n - k0 - 1]


def main():
    out = ROOT / "analysis" / "e250" / "ref"
    os.makedirs(out, exist_ok=True)
    for rt in RTS:
        y = render(rt).astype(np.float32)
        sf.write(out / f"e_moving_rt{rt}.wav", y, int(round(S.FS)), subtype="FLOAT")
        print(f"RT {rt}: {y.shape}, tap {y[S.K0, 1]:+.3f}")


if __name__ == "__main__":
    main()
