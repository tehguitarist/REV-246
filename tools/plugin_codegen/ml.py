"""Programs M and L: analysis/ml/ml_plugin_spec.json -> plugin/src/dsp/generated/ml_spec.h (plugin/docs/ENGINE.md).

    .venv/bin/python -m tools.plugin_codegen.ml

Writes the constants of both programs (lines, coupling, injection and direct taps, diffuser, output kernel, the line
law's tables and pole curve, the (RT, HF) cells of the damping table, theta's rate and depth). The laws themselves are
evaluated in C++ (plugin/src/dsp/programs/MLCore.h); the spec's resolved settings are their test vectors
(plugin/tests/dsp/MLTests.cpp reads them from the JSON).
"""
import math

from .common import Header, load

SPEC = "analysis/ml/ml_plugin_spec.json"
S_SHORTEST = 1321           # S's shortest line: the program scale is round(min(P) / 1321) (ml_lines.line_coefs)


def b_table(spec):
    """The shared damping table (b per step), as the distinct b of both programs' resolved settings."""
    vals = sorted({v["b"] for p in "ML" for v in spec[p]["settings"].values()})
    assert len(vals) == 4 and vals[0] == 0.0, vals
    return vals


def cells(sp, table):
    """Per (RT step, HF step) the damping step index; every cell must be resolved and the same at every LF."""
    rts = sp["controls"]["RT"]["value"]["values_s"]
    hfs = sp["controls"]["HF"]["value"]["values_s"]
    grid = {}
    for v in sp["settings"].values():
        step = table.index(v["b"])
        assert grid.setdefault((v["rt"], v["hf"]), step) == step, ("b differs across LF", v)
    assert len(grid) == len(rts) * len(hfs)
    return rts, hfs, [grid[(r, h)] for r in rts for h in hfs]


def program(h, sp, table, name):
    lines, diff, kern = sp["lines"]["value"], sp["diffuser"]["value"], sp["output_kernel"]["value"]
    P = lines["P"]
    scale = round(min(P) / S_SHORTEST)
    h.raw(f"\nnamespace {name}\n{{")
    h.comment(f"Program {name.upper()} [{sp['lines']['source']}]")
    h.array("lineDelay", P, "int")
    h.scalar("programScale", scale, "int")
    h.comment("Diffuser all-pass delays, in series; the last one carries the damping b [" + sp["diffuser"]["source"] + "]")
    h.array("diffuserDelay", diff["delays"], "int")
    for ch in "LR":
        taps = sp["direct_path"]["value"]["taps"][ch]
        h.comment(f"Direct path (cluster), channel {ch}: delay, sign [{sp['direct_path']['source']}]")
        h.array(f"cluster{ch}Delay", [p for p, _ in taps], "int")
        h.array(f"cluster{ch}Sign", [s for _, s in taps], "int")
    for ch in "LR":
        taps = [t for t in sp["injection"]["value"]["taps"] if t["channel"] == ch]
        h.comment(f"Read-out / injection taps, channel {ch}: delay, sign, loop [{sp['injection']['source']}]")
        h.array(f"tap{ch}Delay", [t["delay"] for t in taps], "int", per_line=10)
        h.array(f"tap{ch}Sign", [t["sign"] for t in taps], "int", per_line=10)
        h.array(f"tap{ch}Loop", [t["line"] for t in taps], "int", per_line=10)
    h.comment(f"Output kernel, causal: centre tap at the cluster's tap 0 [{sp['output_kernel']['source']}]")
    h.array("kernel", kern["taps"])
    h.scalar("kernelCentre", kern["centre"], "int")
    rts, hfs, grid = cells(sp, table)
    h.comment("The library's RT and HF steps and each (RT, HF) cell's damping step (row: RT, column: HF)")
    h.array("rtSteps", rts)
    h.array("hfSteps", hfs)
    h.array("dampingStep", grid, "int", per_line=len(hfs))
    h.raw(f"}} // namespace {name}")


def main():
    spec = load(SPEC)
    M, L = spec["M"], spec["L"]
    for key in ("lines", "coupling", "modulation", "baseline"):        # the shared parts must be the same in both
        if key == "lines":
            a, b = dict(M[key]["value"]), dict(L[key]["value"])
            a.pop("P"), b.pop("P")
        else:
            a, b = M[key], L[key]
        assert a == b, key
    lv, base, mod = M["lines"]["value"], M["baseline"], M["modulation"]["value"]
    table = b_table(spec)
    h = Header("ml", SPEC)
    h.comment("Shared by M and L. The line law: S's control laws at each line's own length [" + M["lines"]["source"] + "]")
    keys = sorted(lv["law_tables"], key=float)
    h.array("lawLf", [float(k) for k in keys])
    for n in ("hi", "dc", "k"):
        h.array(f"law{n.capitalize()}", [lv["law_tables"][k][n] for k in keys])
    h.scalar("shelfPole", lv["a"])
    h.scalar("hfRefHz", lv["f_ref"])
    h.comment("The HF pole curve (piecewise linear in the law's pole; past the fold point the drop-back value)")
    h.array("poleKnots", lv["pole_curve"]["knots"])
    h.array("poleValues", lv["pole_curve"]["values"])
    h.scalar("poleFold", lv["pole_curve"]["fold"])
    h.scalar("poleDropBack", lv["pole_curve"]["drop_back"])
    h.comment("The HF control's floor (RT x HF >= hfrtMin x program scale), L's RT 2.2 floor, the DC floor and its lift")
    h.scalar("hfrtMin", base["hrt_min"])
    h.scalar("l22Hf", base["l22_hf"])
    h.scalar("dcMin", base["dc_min"])
    h.scalar("beta", base["beta"])
    h.comment("Coupling F (rows: lines, columns: loops), orthogonal [" + M["coupling"]["source"] + "]")
    h.array("coupling", [x for row in M["coupling"]["value"]["F"] for x in row])
    h.comment("Diffuser all-pass coefficient")
    h.scalar("allpassG", M["diffuser"]["value"]["g"])
    h.comment("The last all-pass's damping b per step, shared by M and L [" + M["diffuser"]["source"] + "]")
    h.array("dampingTable", table)
    h.comment("theta(t) = amplitude sin(2 pi rate t + phase); amplitude = rms x sqrt(2), the adopted depth "
              "[" + M["modulation"]["source"] + "]")
    h.scalar("thetaRateHz", mod["rate_hz"])
    h.scalar("thetaRms", mod["rms_rad"])
    h.scalar("thetaAmplitude", mod["rms_rad"] * math.sqrt(2.0))
    program(h, M, table, "m")
    program(h, L, table, "l")
    print(f"wrote {h.write('ml_spec.h')}")


if __name__ == "__main__":
    main()
