"""DeviceBridge constants: the unit's analogue chain as poles and residues, plus the test reference.

    .venv/bin/python -m tools.plugin_codegen.bridge          (from the repository root)
Writes plugin/src/dsp/generated/bridge_spec.h (plugin/docs/ENGINE.md, Signal path; docs/NEXT_STEPS.md section 4).

The chain (docs/EMT246_unit_facts.md section 2b):
  input  = input filter 7 246 131 (band part) x the unit's input coupling high-passes 3.35, 0.23, 0.72 Hz
           (C106/R117, C101/R105, C3109/R3134), normalised to unity magnitude at 1 kHz;
  output = output filter 7 246 130 (band part) x the output coupling high-pass 0.16 Hz (C3009/R3031), unity at 1 kHz.
The S&H lag (a device-rate one-pole) is applied by the C++ itself; only the test reference includes it here.
The filters' own sub-Hz coupling poles from the netlists (0.62 and 0.157 Hz, with their zeros at DC) are replaced by
the schematic corners above, as tools/analog/lf_highpass.py does. Capture-chain terms are not part of the unit.

Each side is written as a partial-fraction expansion H(s) = sum_j R_j / (s - p_j) (strictly proper, distinct poles),
one entry per real pole or per conjugate pair (imag > 0, weight 2), which DeviceBridge discretises exactly for
band-limited signals. The reference response is evaluated from the MNA solver (tools/analog/mna.py), independently of
the pole-zero extraction, at the test frequencies.
"""
import sys
from pathlib import Path

import numpy as np

from .common import ROOT, Header

sys.path.insert(0, str(ROOT / "tools" / "analog"))
from adc_clip import h_analog, pole_zero  # noqa: E402
from mna import response  # noqa: E402
from netlists import INPUT_FILTER, OUTPUT_FILTER  # noqa: E402
from sample_hold import IN_AMP_COUPLING_HZ, IN_COUPLING_HZ  # noqa: E402

TWO_PI = 2 * np.pi
FS_DEV = 19.2e6 / 728
SH_POLE = 0.301  # Tomkin's unit (CLA's: 0.325); the reference uses it, DeviceBridge.h holds the same constant
IN_HP_HZ = [IN_COUPLING_HZ, IN_AMP_COUPLING_HZ, 1 / (TWO_PI * 221e3 * 1e-6)]   # 3.35, 0.23, 0.72 Hz
OUT_HP_HZ = [1 / (TWO_PI * 1e6 * 1e-6)]                                        # 0.16 Hz
REF_HZ = np.r_[20 * 640.0 ** (np.arange(37) / 36), 1e3, 10e3]                 # 20 Hz .. 12.8 kHz, log
ALIAS_HZ = np.array([13.6e3, 14.0e3, 16.0e3])                                  # above device Nyquist (13.19 kHz)


def side(net, hp_hz):
    """Band poles/zeros of a netlist (its own sub-Hz coupling dropped) with the given first-order high-passes."""
    p, z = pole_zero(net)
    own_pole = p[np.abs(p) < TWO_PI * 5]
    p, z = p[np.abs(p) >= TWO_PI * 5], z[np.abs(z) >= TWO_PI * 1e-3]
    p = np.r_[p, [-TWO_PI * f for f in hp_hz]]
    z = np.r_[z, np.zeros(len(hp_hz))]
    s1k = TWO_PI * 1j * np.array([1e3])
    k = 1 / abs(h_analog(s1k, p, z)[0])
    return p, z, k, own_pole


def residues(p, z, k):
    return np.array([k * np.prod(pi - z) / np.prod([pi - pj for j, pj in enumerate(p) if j != i]) for i, pi in enumerate(p)])


def representatives(p, R):
    """One entry per real pole or conjugate pair (imag > 0), weight 1 or 2, sorted by decay rate."""
    out = []
    for pi, ri in zip(p, R):
        if abs(pi.imag) < 1e-9 * abs(pi):
            out.append((complex(pi.real, 0.0), complex(ri.real, 0.0), 1))
        elif pi.imag > 0:
            out.append((pi, ri, 2))
    return sorted(out, key=lambda e: -e[0].real)


def mna_side(net, own_pole, hp_hz, f):
    """MNA response with the netlist's own coupling replaced by hp_hz, normalised to unity magnitude at 1 kHz.

    The input filter inverts (its MNA gain constant is negative); both sides are taken non-inverting here, as the
    pole-zero form is, so the chain is +1 at 1 kHz (the unit's absolute polarity is not modelled)."""
    p_all, z_all = pole_zero(net)
    one = np.array([1e3])
    polarity = np.sign((response(net, one) / h_analog(TWO_PI * 1j * one, p_all, z_all))[0].real)

    def h(fr):
        s = TWO_PI * 1j * fr
        out = response(net, fr) * (s - own_pole[0]) / s
        for fc in hp_hz:
            out = out * s / (s + TWO_PI * fc)
        return out
    return polarity * h(np.asarray(f, float)) / abs(h(one)[0])


def lag(f):
    w = TWO_PI * np.asarray(f, float) / FS_DEV
    h = (1 - SH_POLE) / (1 - SH_POLE * np.exp(-1j * w))
    return h / abs((1 - SH_POLE) / (1 - SH_POLE * np.exp(-1j * TWO_PI * 1e3 / FS_DEV)))


def main():
    hdr = Header("bridge", "tools/analog (netlists.py, the unit's coupling corners)")
    hdr.comment("The analogue chain as partial fractions H(s) = sum_j w_j Re-part of R_j / (s - p_j) (rad/s), one entry per\n"
                "real pole (w = 1) or conjugate pair (imag > 0, w = 2), slowest first. Each side is unity at 1 kHz.")
    for name, net, hp in (("in", INPUT_FILTER, IN_HP_HZ), ("out", OUTPUT_FILTER, OUT_HP_HZ)):
        p, z, k, own = side(net, hp)
        R = residues(p, z, k)
        f = np.array([20.0, 1e3, 10e3, 12.5e3, 20e3])
        s = TWO_PI * 1j * f
        pf = np.array([np.sum(R / (si - p)) for si in s])
        ref = mna_side(net, own, hp, f)
        err = np.max(np.abs(pf - ref))
        assert err < 1e-9, (name, err)
        reps = representatives(p, R)
        hdr.comment(f"{name}: {len(p)} poles; partial fractions match the MNA solver to {err:.1e}")
        hdr.array(f"{name}PoleRe", [e[0].real for e in reps])
        hdr.array(f"{name}PoleIm", [e[0].imag for e in reps])
        hdr.array(f"{name}ResRe", [e[1].real for e in reps])
        hdr.array(f"{name}ResIm", [e[1].imag for e in reps])
        hdr.array(f"{name}Weight", [e[2] for e in reps], ctype="int")
    hdr.comment("Coupling corners (Hz) folded into the residues above: input C106/R117, C101/R105, C3109/R3134; output C3009/R3031.")
    hdr.array("inHighPassHz", IN_HP_HZ)
    hdr.array("outHighPassHz", OUT_HP_HZ)

    # Test reference: the whole chain (R channel) from the MNA solver, S&H lag p = refShPole, unity at 1 kHz.
    p_in, z_in, k_in, own_in = side(INPUT_FILTER, IN_HP_HZ)
    p_out, z_out, k_out, own_out = side(OUTPUT_FILTER, OUT_HP_HZ)
    chain = mna_side(INPUT_FILTER, own_in, IN_HP_HZ, REF_HZ) * mna_side(OUTPUT_FILTER, own_out, OUT_HP_HZ, REF_HZ) * lag(REF_HZ)
    hdr.comment("Test reference (tests/dsp/BridgeTests.cpp): chain response at refHz, R channel, after removing the latency.")
    hdr.scalar("refShPole", SH_POLE)
    hdr.array("refHz", REF_HZ)
    hdr.array("refRe", chain.real)
    hdr.array("refIm", chain.imag)
    # Aliasing: a tone at f above device Nyquist appears at f' = fs_dev - f with amplitude conj(H_in(f)) H_lag(f') H_out(f')
    # (device instant 0 at host time 0).
    fold = FS_DEV - ALIAS_HZ
    alias = np.conj(mna_side(INPUT_FILTER, own_in, IN_HP_HZ, ALIAS_HZ)) * mna_side(OUTPUT_FILTER, own_out, OUT_HP_HZ, fold) * lag(fold)
    hdr.comment("Aliasing reference: input tone aliasHz, output at deviceRate - aliasHz, complex amplitude.")
    hdr.array("aliasHz", ALIAS_HZ)
    hdr.array("aliasRe", alias.real)
    hdr.array("aliasIm", alias.imag)
    hdr.array("aliasInputGainDb", 20 * np.log10(np.abs(mna_side(INPUT_FILTER, own_in, IN_HP_HZ, ALIAS_HZ))))
    path = hdr.write("bridge_spec.h")
    print("wrote", path.relative_to(ROOT))
    print("chain dB at refHz:", np.round(20 * np.log10(np.abs(chain)), 3))
    print("alias dB:", np.round(20 * np.log10(np.abs(alias)), 2))


if __name__ == "__main__":
    main()
