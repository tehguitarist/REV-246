"""The 16-bit storage reference: tools/loop_read/c/s_loops_fx.c run on S's four loops -> plugin/src/dsp/generated/sixteen_ref.h.

    .venv/bin/python tools/memguard.py --limit-gb 2.5 -- .venv/bin/python -m tools.plugin_codegen.sixteen_ref

Compiles s_loops_fx.c (model_s_results §36.5, Q10) and runs it in truncation mode (1) on the three §36.5 settings with
the same per-loop coefficients (analysis/lowband/chain3_fixedpoint_coefs.json, as fixedpoint.q10). The header keeps the
parameters, how the input is made (an impulse, or an integer xorshift32 noise the test regenerates) and a 64-bit FNV-1a
hash of the 16-bit output words, so plugin/tests/dsp/SixteenBitTests.cpp can check ProgramCore's 16-bit storage bit for
bit without committing the renders. The loops' input gain b is §36.5's scaling A or B rounded to a power of two: the
plugin realises a core's internal gain as the storage grid fs / (32768 k), which equals s_loops_fx's Q15 multiplier B
exactly only when b is a power of two (both are exact then).
"""
import os
import subprocess
import tempfile

import numpy as np

from .common import ROOT, Header, load

COEFS = "analysis/lowband/chain3_fixedpoint_coefs.json"
SRC = ROOT / "tools" / "loop_read" / "c" / "s_loops_fx.c"
P = [1321, 1389, 1461, 1587]
SG = [1, -1, 1, -1]
FS_DEV = 19.2e6 / 728


def xorshift_noise(seed, n, amp):
    """Integer noise in [-amp, amp] from xorshift32 (13, 17, 5); the C++ test makes the same sequence."""
    s = seed & 0xFFFFFFFF
    out = np.empty(n, np.int64)
    m = 2 * amp + 1
    for i in range(n):
        s ^= (s << 13) & 0xFFFFFFFF
        s ^= s >> 17
        s ^= (s << 5) & 0xFFFFFFFF
        out[i] = (s % m) - amp
    return out


def fnv1a(words):
    h = 0xCBF29CE484222325
    for w in words:
        for byte in int(w & 0xFFFFFFFF).to_bytes(4, "little"):
            h ^= byte
            h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def main():
    d = load(COEFS)
    # (setting, scaling, input, seconds): impulse at -6 dBFS, low-level noise (the idle case), full-scale noise at a hot
    # scaling (saturates: checks the clamp)
    base_cases = [("1.0_2.5_1.0", "A", ("impulse", 0, 16384), 3.0),
                  ("1.0_1.0_1.0", "A", ("noise", 12345, 8), 2.0),
                  ("2.0_0.8_0.2", "hot", ("noise", 777, 32767), 1.0)]
    cases = [c + (mode,) for mode in (1, 2) for c in base_cases]   # truncation, then rounding
    h = Header("sixteen", f"{COEFS} through tools/loop_read/c/s_loops_fx.c (modes 1 and 2)")
    h.comment("s_loops_fx.c's four-loop network in 16-bit truncation: per case the loop delays, c0, c1, a, p, b, sign;\n"
              "the input (impulse: one sample of `inputAmp` LSB at 0; noise: xorshift32 from `inputSeed`, uniform in\n"
              "[-inputAmp, inputAmp] LSB); the output's 16-bit words hashed with FNV-1a 64 (each word as 4 bytes, little\n"
              "endian), their sum of squares and the number of saturated stores; `mode` 1 truncation, 2 rounding.")
    h.scalar("numCases", len(cases), "int")
    h.array("loopDelay", P, "int")
    h.array("loopSign", SG, "int")
    params, kinds, seeds, amps, lens, hashes, sumsq, clips, bs, modes = [], [], [], [], [], [], [], [], [], []
    with tempfile.TemporaryDirectory() as tmp:
        exe = os.path.join(tmp, "s_loops_fx")
        subprocess.run(["cc", "-O2", "-o", exe, str(SRC), "-lm"], check=True)
        for key, scl, (kind, seed, amp), secs, mode in cases:
            rows = sorted((r for r in d["rows"] if r["key"] == key), key=lambda r: r["loop"])
            base = [r["base"] for r in rows]
            gmax = np.array([max(b["D"], b["N1"]) for b in base])
            n = int(secs * FS_DEV)
            if scl == "A":
                b = 1 / np.sum(1 / (1 - gmax))                       # §36.5 scaling A: no write or output can overflow
            else:
                b = 0.5                                              # hot: full-scale noise saturates
            b = 2.0 ** round(np.log2(b))
            par = [(P[i], base[i]["c0"], base[i]["c1"], d["a"], base[i]["p"] or 0.0, b, SG[i]) for i in range(4)]
            x = np.zeros(n)
            if kind == "impulse":
                x[0] = amp / 32768.0
            else:
                x = xorshift_noise(seed, n, amp) / 32768.0
            pf, xf, of = (os.path.join(tmp, f) for f in ("par.txt", "x.f64", "y.f64"))
            with open(pf, "w") as f:
                f.write(f"{n}\n" + "".join(f"{r[0]} " + " ".join(repr(float(v)) for v in r[1:6]) + f" {r[6]}\n" for r in par))  # P, sign: %d / integers
            x.astype(np.float64).tofile(xf)
            r = subprocess.run([exe, pf, str(mode), xf, of], check=True, capture_output=True, text=True)
            nclip = int(r.stdout.strip() or 0)
            y = np.fromfile(of, np.float64)
            words = np.round(y * 32768.0).astype(np.int64)
            assert np.all(words * 1.0 == y * 32768.0)
            print(f"mode {mode} {key} b {b} {kind}: {n} samples, rms {np.sqrt(np.mean(words ** 2.0)):.2f} LSB, {nclip} saturated")
            params += [v for row in par for v in row[1:5]]
            bs.append(b)
            modes.append(mode)
            kinds.append(1 if kind == "impulse" else 0)
            seeds.append(seed)
            amps.append(amp)
            lens.append(n)
            hashes.append(fnv1a(words))
            sumsq.append(int(np.sum(words ** 2)))
            clips.append(nclip)
    h.array("coef", params)                        # per case and loop: c0, c1, a, p
    h.array("inputGain", bs)                       # b per case (a power of two)
    h.array("mode", modes, "int")
    h.array("impulse", kinds, "int")
    h.array("inputSeed", seeds, "std::uint32_t")
    h.array("inputAmp", amps, "int")
    h.array("length", lens, "int")
    h.raw("inline constexpr std::array<std::uint64_t, %d> outputHash{{%s}};" % (len(hashes), ", ".join(f"0x{v:016x}ull" for v in hashes)))
    h.raw("inline constexpr std::array<std::int64_t, %d> outputSumSq{{%s}};" % (len(sumsq), ", ".join(f"{v}ll" for v in sumsq)))
    h.array("saturated", clips, "long")
    print(h.write("sixteen_ref.h"))


if __name__ == "__main__":
    main()
