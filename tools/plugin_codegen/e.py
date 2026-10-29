"""Program E's constants for the plugin: writes plugin/src/dsp/generated/e_spec.h (plugin/docs/ENGINE.md).

    .venv/bin/python -m tools.plugin_codegen.e

Reads analysis/e250/e_plugin_spec.json for the network, and the 246's loop law that the spec names
(loop_filter.law: tools/ml_lines.line_coefs / line_resp) from its own sources: analysis/q1/laws.json (the LF law and the
shelf pole), analysis/lowband/chain3_pole_cells_allf.json (the HF pole curve) and ml_lines' coefficient limits. The
law's constants are written, not tabulated; the C++ evaluates it. ml_lines.line_coefs itself supplies the test vectors
(per loop g, r0, p at the measured RT steps and at off-step panel values), and e_sim.e_reads' own draw order with
seed 0 the moving reads' rates and phases (one fixed realisation, the start phase at reset()).
"""
import os
import sys

import numpy as np

from .common import ROOT, Header, load

sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "e250"))
for _k in ("ML_POLE_CURVE", "ML_HFRT", "ML_L22"):
    if os.environ.get(_k):
        raise SystemExit(f"{_k} is set: the generator writes the adopted law only")
import ml_lines as mll  # noqa: E402

SPEC = "analysis/e250/e_plugin_spec.json"
LINE_UNIT = 1321          # ml_lines.line_coefs: the HF floor's program-time factor round(min(P) / 1321)
LAW_VECTORS = [           # (RT, LF, HF): the measured steps at the family values, then off-step panel values
    *[(rt, 1.49, 0.69) for rt in (0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5)],
    (0.55, 1.49, 0.69), (5.0, 1.49, 0.69), (1.0, 1.0, 0.9), (1.0, 0.1, 0.2), (0.4, 2.0, 0.1), (2.2, 0.4, 1.0),
    (3.3, 0.7, 0.5), (5.0, 2.0, 1.0), (0.4, 0.1, 0.05), (1.7, 1.2, 1.2), (4.8, 0.25, 0.3), (0.8, 3.0, 0.8),
]


def main():
    sp = load(SPEC)
    law = mll.LAW
    pc = load("analysis/lowband/chain3_pole_cells_allf.json")
    h = Header("e", SPEC)
    num = lambda k: sp[k]["value"]

    h.comment("The network (spec: direct_path, allpasses, loops; EMT250_program §11.2, §17). Delays in device samples after\n"
              "R's direct tap, which the reference IRs put at k0 (16 of those samples are the output kernel's offset).")
    h.scalar("k0", sp["k0"], "int")
    d = num("direct_path")
    h.array("directTaps", d["taps"], "int")
    h.array("directGains", d["gains"])
    h.scalar("lDirectDelay", d["L"]["delay"], "int")
    h.scalar("lDirectSign", d["L"]["sign"])
    h.scalar("lGain", d["L"]["gain"])
    ap = num("allpasses")
    h.array("apDelay", ap["D"], "int")
    h.array("apGain", ap["g"])
    lp = num("loops")
    h.array("loopLength", lp["P"], "int")
    h.array("staticReadDelay", lp["static_read_delay"], "int")
    h.array("staticReadGain", lp["static_read_gain"])
    h.scalar("lReverbDelay", lp["L"]["delay"], "int")
    h.scalar("lReverbSign", lp["L"]["sign"])
    fam = num("loop_filter")["altiverb_family"]
    h.comment("The Altiverb family's LF / HF through the loop law (spec: loop_filter.altiverb_family, §16.5, §17.3)")
    h.scalar("familyLf", fam["lf"])
    h.scalar("familyHf", fam["hf"])
    h.scalar("amplRef", 80.0)

    h.comment("The 246's loop law (ml_lines.line_coefs / line_resp; analysis/q1/laws.json, model_s_results §14.2, §16.15-16.16, §29,\n"
              "§32.6): per pass G(z) = g (b0 + b1 z^-1) / (1 - a z^-1) (1 - p) / (1 - p z^-1).")
    h.scalar("lawShelfPole", law["a"])
    h.scalar("lawFRef", law["f_ref"])
    ks = sorted(law["lf"], key=float)
    h.array("lawLf", [float(k) for k in ks])
    for n in ("hi", "dc", "k"):
        h.array(f"lawLf_{n}", [law["lf"][k][n] for k in ks])
    h.scalar("pMax", mll.P_MAX)
    h.scalar("dcMin", mll.DC_MIN)
    h.scalar("beta", mll.BETA)
    h.comment("HF floor per loop: HF >= hfFloor[k] / RT (HFRT_MIN x round(P / 1321), model_s_results §29)")
    h.array("hfFloor", [mll.HFRT_MIN * round(P / LINE_UNIT) for P in lp["P"]])
    h.comment("The HF pole curve (chain3_pole_cells_allf.json 'all', model_s_results §32.6): piecewise linear below the fold,\n"
              "the drop-back value at or past it")
    h.array("poleKnots", pc["knots"])
    h.array("poleCurve", pc["all"])
    h.scalar("poleFold", pc["fold"])
    h.scalar("poleDropBack", mll.DROP_BACK)
    h.scalar("onepoleT2Max", 0.999999)

    mv = num("moving_reads")
    ad = mv["adopted"]
    reads = mv["reads"]
    h.comment("Moving reads (spec: moving_reads, §16-17): nine reads on the loops' v (before the loop filter), each at\n"
              "centre + A wave(f_j t + phase_j) rounded to an integer (no interpolation); L at the same address with its own final\n"
              "amplitude, mdel later (+ dlSame for the reads whose L and R signs agree); a one-pole mlp on both sums.")
    h.scalar("numReads", len(reads), "int")
    h.array("readCentre", [r["centre"] for r in reads], "int")
    h.array("readLoop", ["ABC".index(r["loop"]) for r in reads], "int")
    h.array("readR", [r["R"] for r in reads])
    h.array("readL", [r["L"] for r in reads])
    h.array("readLExtra", [ad["dl_same"] if r["R"] * r["L"] > 0 else 0 for r in reads], "int")
    h.scalar("readDepth", ad["A"])
    h.scalar("readRate", mv["lfo_rate_hz"]["time"])
    h.scalar("readRateSpread", ad["fsd"])
    h.scalar("movingLDelay", ad["mdel"], "int")
    h.scalar("movingOnePole", ad["mlp"])
    # e_sim.e_reads' draws with seed 0: one common phase (unused, common = false), then per read its rate factor and phase
    rng = np.random.default_rng(0)
    rng.uniform()
    rates, phases = [], []
    for _ in reads:
        rates.append(ad["f"] * rng.uniform(1 - ad["fsd"], 1 + ad["fsd"]))
        phases.append(rng.uniform())
    h.comment("One fixed realisation (H): e_sim.e_reads' draws at seed 0, rate f x U(1 - fsd, 1 + fsd) (Hz) and start phase\n"
              "(cycles) per read; the wave is the spec's triangle 4 |frac(phase) - 0.5| - 1 (e_sim.ADOPTED).")
    h.array("readRateHz", rates)
    h.array("readPhase", phases)

    ok = num("output_kernel")
    h.comment("Output kernel on both channels (spec: output_kernel, §14.1), at k0 and k0 + 1")
    h.array("outKernel", ok["fir"])

    h.comment("Test vectors: ml_lines.line_coefs([P], RT, LF, HF) per loop -> (g, r0, p), for (RT, LF, HF) in lawVectorIn")
    h.array("lawVectorIn", [v for t in LAW_VECTORS for v in t], per_line=3)
    out = []
    for rt, lf, hf in LAW_VECTORS:
        for P in lp["P"]:
            out += list(mll.line_coefs([P], rt, lf, hf)[0])
    h.array("lawVectorOut", out, per_line=3)
    print("wrote", h.write("e_spec.h"))


if __name__ == "__main__":
    main()
