"""Program S: analysis/s/s_plugin_spec.json -> plugin/src/dsp/generated/s_spec.h (plugin/docs/ENGINE.md).

    .venv/bin/python -m tools.plugin_codegen.s

Writes S's constants: the loop delays, the injection (read-out) and direct taps, the diffuser and AP230's damping
steps with their (RT, HF) cell table, the loop-filter law's tables, the RT folds' constants and the per-setting cells.
The law itself is evaluated in C++ (plugin/src/dsp/programs/SCore.h); the spec's 128 resolved settings are its test
vectors (plugin/tests/dsp/STests.cpp reads them from the JSON).

The spec names, but does not copy, the fold and cell files its law reads (loop_filter_law.value.folds / cells, under
baseline.fold_set); they are read from there. The HF pole off the library's steps follows model_s_results §32.6 /
§32.8 (the all-data pole curve); that curve is S's own, embedded in the M / L spec (lines.value.pole_curve), and is
read from there.
"""
import json

from .common import ROOT, Header, load

SPEC = "analysis/s/s_plugin_spec.json"
ML_SPEC = "analysis/ml/ml_plugin_spec.json"
# s_spec.py's baseline switches the fold files assume (tools/q2_injection.py with EMT_PMAX=0.8475, EMT_HFRT_MIN=0.12)
EXPECTED_BASELINE = dict(kernel="chain3", p_max=0.8475, hrt_min=0.12, dc_min=0.18, beta=0.019)


def key(lf, rt, hf):
    return f"{lf}_{rt}_{hf}"


def main():
    sp = load(SPEC)
    base = sp["baseline"]
    for k, v in EXPECTED_BASELINE.items():
        assert base[k] == v, (k, base[k])
    fold_dir = ROOT / base["fold_set"]
    lfs = sp["controls"]["LF"]["value"]["values_s"]
    rts = sp["controls"]["RT"]["value"]["values_s"]
    hfs = sp["controls"]["HF"]["value"]["values_s"]
    assert len(lfs) == 4 and len(rts) == 8 and len(hfs) == 4
    law = sp["loop_filter_law"]["value"]
    loops = sp["loops"]["value"]
    diff = sp["diffuser"]["value"]
    inj = sp["injection"]["value"]
    direct = sp["direct_path"]["value"]
    assert abs(sp["device_rate_hz"] - 19.2e6 / 728) < 1e-9

    h = Header("s", SPEC)
    h.comment("Program S at the device rate [docs/model_s_results.md §17, §22, §32, §38.1]. Delays in device samples.")
    h.scalar("deviceRate", sp["device_rate_hz"])
    h.comment(f"Loop delays P_i [{sp['loops']['source']}]")
    h.array("loopDelay", loops["P"], "int")

    h.comment(f"Injection as read-out taps: delay from tap 0, sign, loop, channel (0 L, 1 R); magnitude [{sp['injection']['source']}]")
    taps = sorted(inj["taps"], key=lambda t: ("LR".index(t["channel"]), t["loop"], t["delay"]))
    h.array("tapDelay", [t["delay"] for t in taps], "int", per_line=10)
    h.array("tapSign", [t["sign"] for t in taps], "int", per_line=10)
    h.array("tapLoop", [t["loop"] for t in taps], "int", per_line=10)
    h.array("tapChannel", ["LR".index(t["channel"]) for t in taps], "int", per_line=10)
    h.scalar("tapMagnitude", inj["magnitude"])

    h.comment(f"Direct path (cluster) per channel: delay and gain at AMPL 90; x AMPL / 90 [{sp['direct_path']['source']}]")
    for ch in "LR":
        h.array(f"direct{ch}Delay", direct["taps"][ch], "int")
        h.array(f"direct{ch}Gain", direct["gains"][ch])
    h.scalar("amplReference", 90.0)

    h.comment(f"Diffuser: AP20 -> AP72 -> AP230 (damping b), one gain g [{sp['diffuser']['source']}]")
    aps = diff["allpasses"]
    assert all(a["g"] == diff["ap230"]["g"] for a in aps)
    h.array("allpassDelay", [a["M"] for a in aps] + [diff["ap230"]["M"]], "int")
    h.scalar("allpassGain", diff["ap230"]["g"])
    steps = diff["ap230"]["b_steps"]
    h.array("ap230BSteps", steps)
    cell_step = {}
    for v in sp["settings"].values():
        s = steps.index(v["ap230_b"])
        assert cell_step.setdefault((v["rt"], v["hf"]), s) == s, ("AP230 step differs across LF", v)
    h.comment("AP230 damping step per (RT step, HF step) cell, index rt * 4 + hf; the same at every LF")
    h.array("ap230Step", [cell_step[(r, f)] for r in rts for f in hfs], "int", per_line=4)

    h.comment(f"The library's steps (s) [{sp['controls']['RT']['source']}]")
    h.array("lfSteps", lfs)
    h.array("rtSteps", rts)
    h.array("hfSteps", hfs)

    h.comment(f"Loop-filter law per LF step: hi, dc, k; shelf pole a; HF reference f_ref [{sp['loop_filter_law']['source']}]")
    nodes = [law["law"][str(lf)] for lf in lfs]
    h.array("lawHi", [n["hi"] for n in nodes])
    h.array("lawDc", [n["dc"] for n in nodes])
    h.array("lawK", [n["k"] for n in nodes])
    h.scalar("shelfPole", law["a"])
    h.scalar("hfRefHz", law["f_ref"])
    h.comment("Baseline limits: DC floor, its lift factor beta, the HF control's floor RT x HF (s) [§22.3, §29]")
    h.scalar("dcMin", base["dc_min"])
    h.scalar("dcBeta", base["beta"])
    h.scalar("hfRtMin", base["hrt_min"])

    # The RT folds (fold_<RT>.json): gamma, LF 0.4's DC ratio, the LF-in-RT correction. hfdp / pramp are dead on the
    # baseline (§32.10: every HF < 1 setting has its cell's hf_poles), so they are not carried.
    gamma, lfdc, lfdc_on, lfrt_a, lfrt_b, lfrt_on, lfrt_clip = [], [], [], [], [], [], []
    for rt in rts:
        f = json.loads((fold_dir / f"fold_{rt}.json").read_text())
        assert f["lfrt"]["form"] == "a + b RT"
        gamma.append(f["gamma"])
        for lf in lfs:
            d = f.get("lfdc", {}).get(str(lf))
            lfdc.append(d if d is not None else 0.0)
            lfdc_on.append(int(d is not None))
            c = f["lfrt"]["coefs"].get(str(lf))
            lfrt_a.append(c[0] if c else 0.0)
            lfrt_b.append(c[1] if c else 0.0)
            lfrt_on.append(int(c is not None))
            lfrt_clip.append(int(str(lf) in f["lfrt"].get("clip", [])))
    h.comment(f"RT folds ({base['fold_set']}/fold_<RT>.json), index rt [* 4 + lf]: gamma (dB per 60 dB of decay, "
              "per pass x P / fs), LF's DC ratio override, the LF correction rel = a + b RT (clipped at 0 where flagged)")
    h.array("foldGamma", gamma)
    h.array("foldLfDc", lfdc)
    h.array("foldLfDcOn", lfdc_on, "int", per_line=16)
    h.array("foldLfRtA", lfrt_a)
    h.array("foldLfRtB", lfrt_b)
    h.array("foldLfRtOn", lfrt_on, "int", per_line=16)
    h.array("foldLfRtClip", lfrt_clip, "int", per_line=16)

    cells = json.loads((fold_dir / "cells.json").read_text())
    kind, dc, dcmin, beta, cgamma, gdc, gain, cdc, poles_on, poles = [], [], [], [], [], [], [], [], [], []
    ff = {n: [] for n in ("F", "c", "tau", "F0", "beta", "gamma")}
    for lf in lfs:
        for rt in rts:
            for hf in hfs:
                c = cells[key(lf, rt, hf)]
                assert not ("dc_min" in c and "floor_form" in c)
                kind.append(1 if "dc_min" in c else 2 if "floor_form" in c else 0)
                dc.append(c.get("dc", 0.0))
                dcmin.append(c.get("dc_min", 0.0))
                beta.append(c.get("beta", 0.0))
                cgamma.append(c.get("gamma", 0.0))
                gdc.append(int(c.get("gamma_dc", True)))
                gain += c["lf_gain"]
                cdc += c["lf_dc"]
                poles_on.append(int("hf_poles" in c))
                poles += c.get("hf_poles", [0.0] * 4)
                fl = c.get("floor_form", {})
                for n in ff:
                    ff[n].append(fl.get(n, 0.0))
    h.comment(f"Per-setting cells ({base['fold_set']}/cells.json), index (lf * 8 + rt) * 4 + hf [x 4 + loop]: kind 0 the "
              "fold's floor, 1 own DC floor (dc for LF 0.1's law, dc_min, beta, gamma, gamma_dc), 2 floor_form; every "
              "cell's per-loop gain and DC shelf (dB); the HF poles at the cell's own step [§32.2-32.6, §35]")
    h.array("cellKind", kind, "int", per_line=16)
    h.array("cellDc", dc)
    h.array("cellDcMin", dcmin)
    h.array("cellBeta", beta)
    h.array("cellGamma", cgamma)
    h.array("cellGammaDc", gdc, "int", per_line=16)
    h.array("cellLoopGainDb", gain)
    h.array("cellLoopDcDb", cdc)
    h.array("cellPolesOn", poles_on, "int", per_line=16)
    h.array("cellPoles", poles)
    for n, v in ff.items():
        h.array("cellFloor" + n[0].upper() + n[1:], v)

    pc = load(ML_SPEC)["M"]["lines"]["value"]["pole_curve"]
    assert pc == load(ML_SPEC)["L"]["lines"]["value"]["pole_curve"]
    h.comment("S's HF pole curve, all-data fit (piecewise linear in the law's pole), and its drop-back past the fold "
              f"point [§32.4, §32.6; from {ML_SPEC} lines.value.pole_curve]")
    h.array("poleCurveKnots", pc["knots"])
    h.array("poleCurveValues", pc["values"])
    h.scalar("poleCurveFold", pc["fold"])
    h.scalar("poleCurveDropBack", pc["drop_back"])
    out = h.write("s_spec.h")
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
