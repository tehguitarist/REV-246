"""Generator for plugin/src/dsp/generated/n_spec.h: program N's constant tables from analysis/nl/nl_plugin_spec.json
(plugin/docs/ENGINE.md, Generated data).

    .venv/bin/python -m tools.plugin_codegen.n

The spec's laws stay laws (the early law, the late gain, the loop gain, the SPAN rule are evaluated in
plugin/src/dsp/programs/NCore.h); this writes their constants, the tap tables, the branch read kernel, and a block of
test vectors: the resolved gains of every library (SPAN, EARLY) pair as tools/nl_early.law_db and tools/nl_model.late_db
give them, which the C++ law is checked against (plugin/tests/dsp/NTests.cpp).
"""
import os
import re
import sys

from .common import ROOT, Header, load

SPEC = "analysis/nl/nl_plugin_spec.json"


def law_vectors(spans, earlies, n_stage_max):
    """Per library (SPAN, EARLY): ns, front dB L (3), front dB R (3), stage dB (padded to n_stage_max), late dB."""
    # the same environment tools/nl_spec.py sets, so importing the Python model does not read the kernel files
    os.environ.update(NL_KERNEL="0", NL_BK="0", NL_WARP2="1", NL_BK_ALPHA="0.7", NL_BK_LP="11000")
    os.environ.pop("NL_BK_FILE", None)
    sys.path.insert(0, str(ROOT / "tools"))
    import nl_early as ne
    import nl_model as nm
    out = []
    for span in spans:
        for early in earlies:
            ns = ne.N_STAGES[span]
            fl, sdb = ne.law_db(span, early, 0)
            fr, _ = ne.law_db(span, early, 1)
            out += [float(ns), *map(float, fl), *map(float, fr), *(list(map(float, sdb)) + [0.0] * (n_stage_max - ns)),
                    float(nm.late_db(span, early))]
    return out


def main():
    s = load(SPEC)
    h = Header("n", SPEC)
    ft, bi, mt = s["front_taps"]["value"], s["branch_input"]["value"], s["main_taps"]["value"]
    late, bk = s["late"]["value"], s["branch_read_kernel"]["value"]
    ctl, el = s["controls"], s["early_law"]["value"]
    span = s["span"]["value"]["stages_running"]

    h.comment("Device rate of the spec (19.2 MHz / 728); NCore uses rv::dsp::deviceRate and the test checks they agree.")
    h.scalar("deviceRateHz", s["device_rate_hz"])

    h.comment("Front taps: plain delays of the input (device samples), amplitude in main-tap units (signs included), and\n"
              "the tick of the early law each gain follows (law time = tick x tickS) [spec front_taps].")
    for ch in "LR":
        v = ft[ch]
        h.array(f"frontDelay{ch}", v["delay"], "int")
        h.array(f"frontAmp{ch}", v["amp"], "double")
        h.array(f"frontTick{ch}", v["law_tick"], "int")

    h.comment("Early law [spec early_law]: the time quantum, and the ticks between two stages (stage j at 4 + 9 j ticks).")
    h.scalar("tickS", el["tick_s"])
    h.scalar("stageTicks", round(el["stage_s"] / el["tick_s"]), "int")
    h.scalar("stageTick0", 4, "int")
    h.scalar("negativeEndTicks", 3, "int")

    h.comment("The branch input all-pass (AP124), and every all-pass's g [spec branch_input, main_taps.out_allpass].")
    h.scalar("apInDelay", bi["allpass"]["D"], "int")
    h.scalar("apG", bi["allpass"]["g"])
    h.scalar("apOutDelayL", mt["out_allpass"]["L"], "int")
    h.scalar("apOutDelayR", mt["out_allpass"]["R"], "int")
    h.scalar("branchDelay", mt["branch_delay"]["value"], "int")

    h.comment("Main taps: delays from the table's own origin (the branch delay is added by the core), signs per channel,\n"
              "and the last tap of each stage (a stage is the taps sharing one gain) [spec main_taps].")
    h.array("mainDelayL", mt["L_delay"], "int", 8)
    h.array("mainDelayR", mt["R_delay"], "int", 8)
    h.array("mainSignL", mt["L_sign"], "int", 27)
    h.array("mainSignR", mt["R_sign"], "int", 27)
    ends = [st[-1] for st in mt["stages"]]
    assert [st for st in mt["stages"]] == [[0]] + [[i, i + 1] for i in range(1, 26, 2)], "stage layout changed"
    h.array("stageEnd", ends, "int", 14)

    h.comment("SPAN: stages running at each library step, and the late feed's position on the branch line there\n"
              "(device samples before the branch's read, one before a main tap) [spec span, late.feed].")
    spans = sorted(float(k) for k in span)
    h.array("libSpan", spans, "double")
    h.array("libStages", [span[str(k)] for k in spans], "int")
    h.array("libFeed", [late["feed"]["delay_on_branch_line"][str(k)] for k in spans], "int")
    r_late = int(re.search(r"\d+", late["feed"]["R"]).group())
    h.comment("R's late feed is this many samples earlier than L's [spec late.feed.R]; NCore reads R's loops that much nearer.")
    h.scalar("lateRShift", r_late, "int")

    h.comment("Late section: the all-pass after the branch line's feed, the three loops, the 4 reads per loop (positions\n"
              "from the feed, signs per channel), and the reads' amplitude [spec late].")
    h.scalar("apLateDelay", late["allpass"]["D"], "int")
    h.array("loopDelay", late["loops"]["P"], "int")
    keys = [str(p) for p in late["loops"]["P"]]
    h.array("readDelay", [d for k in keys for d in late["reads"][k]["delay"]], "int", 4)
    h.array("readSignL", [d for k in keys for d in late["reads"][k]["L_sign"]], "int", 4)
    h.array("readSignR", [d for k in keys for d in late["reads"][k]["R_sign"]], "int", 4)
    h.scalar("readAmp", late["read_amp"])

    h.comment("Panel values the spec lists (tests; the laws accept any value in range) [spec controls].")
    h.array("spanValues", ctl["SPAN (AMPL slot)"]["value"]["values_s"])
    h.array("earlyValues", ctl["EARLY RT (LF slot)"]["value"]["values_s"])
    h.array("lateValues", ctl["LATE RT (HF slot)"]["value"]["values_s"])

    h.comment("Branch read kernel [spec branch_read_kernel]: the stored taps are the fitted, strength-1 kernel; the core\n"
              "uses delta + strength x (taps - delta). Centre = look-ahead = latency in samples.")
    h.array("bkTapsL", bk["taps"]["L"], "double", 4)
    h.array("bkTapsR", bk["taps"]["R"], "double", 4)
    h.scalar("bkCentre", bk["centre"], "int")
    h.scalar("bkStrength", bk["strength"])
    h.scalar("bkBandLimitHz", bk["band_limit_hz"])
    assert len(bk["taps"]["L"]) == 2 * bk["centre"] + 1 == bk["latency_samples"] * 2 + 1

    n_stage_max = len(mt["stages"])
    vals = law_vectors(spans, ctl["EARLY RT (LF slot)"]["value"]["values_s"], n_stage_max)
    h.comment("Test vectors: for each library SPAN (outer) and EARLY (inner) the row\n"
              "  ns, front dB L (3), front dB R (3), stage dB (14, zero beyond ns), late dB\n"
              "from tools/nl_early.law_db and tools/nl_model.late_db (the Python bench's own law).")
    h.scalar("lawVectorRow", 1 + 3 + 3 + n_stage_max + 1, "int")
    h.array("lawVectors", vals, "double", 8)
    print("wrote", h.write("n_spec.h"))


if __name__ == "__main__":
    main()
