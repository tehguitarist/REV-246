"""Program E's constants for the plugin: writes plugin/src/dsp/generated/e_spec.h (plugin/docs/ENGINE.md).

    .venv/bin/python -m tools.plugin_codegen.e

Reads analysis/e250/e_plugin_spec.json alone. The 246's loop law that the spec names (loop_filter.law:
tools/ml_lines.line_coefs / line_resp) is embedded in it as loop_law (the LF law, the shelf pole, the coefficient limits
and the HF pole curve); the law's constants are written, not tabulated, and the C++ evaluates it. The test vectors (per
loop g, r0, p at the measured RT steps and at off-step panel values) are the spec's law_vectors, and the moving reads'
rates and phases (one fixed realisation, e_sim.e_reads' draws at seed 0; the start phase at reset()) its
moving_reads_draw. tools/e250/e_spec.py (embed_codegen_data) computes all three.
"""
from .common import Header, load

SPEC = "analysis/e250/e_plugin_spec.json"


def main():
    sp = load(SPEC)
    law = sp["loop_law"]["value"]
    pc = law["pole_curve"]
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
    h.scalar("pMax", law["p_max"])
    h.scalar("dcMin", law["dc_min"])
    h.scalar("beta", law["beta"])
    h.comment("HF floor per loop: HF >= hfFloor[k] / RT (HFRT_MIN x round(P / 1321), model_s_results §29)")
    h.array("hfFloor", [law["hfrt_min"] * round(P / law["line_unit"]) for P in lp["P"]])
    h.comment("The HF pole curve (chain3_pole_cells_allf.json 'all', spec loop_law.pole_curve, model_s_results §32.6): piecewise linear below the fold,\n"
              "the drop-back value at or past it")
    h.array("poleKnots", pc["knots"])
    h.array("poleCurve", pc["values"])
    h.scalar("poleFold", pc["fold"])
    h.scalar("poleDropBack", pc["drop_back"])
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
    dr = sp["moving_reads_draw"]["value"]
    rates, phases = dr["rate_hz"], dr["phase"]
    h.comment("One fixed realisation (H): e_sim.e_reads' draws at seed 0, rate f x U(1 - fsd, 1 + fsd) (Hz) and start phase\n"
              "(cycles) per read; the wave is the spec's triangle 4 |frac(phase) - 0.5| - 1 (e_sim.ADOPTED).")
    h.array("readRateHz", rates)
    h.array("readPhase", phases)

    ok = num("output_kernel")
    h.comment("Output kernel on both channels (spec: output_kernel, §14.1), at k0 and k0 + 1")
    h.array("outKernel", ok["fir"])

    lv = sp["law_vectors"]["value"]
    assert lv["loops_P"] == lp["P"]
    h.comment("Test vectors: ml_lines.line_coefs([P], RT, LF, HF) per loop -> (g, r0, p), for (RT, LF, HF) in lawVectorIn (spec: law_vectors)")
    h.array("lawVectorIn", [v for t in lv["inputs"] for v in t], per_line=3)
    h.array("lawVectorOut", [v for t in lv["outputs"] for v in t], per_line=3)
    print("wrote", h.write("e_spec.h"))


if __name__ == "__main__":
    main()
