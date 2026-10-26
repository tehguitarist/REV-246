"""Program E for the plugin: the specification as data, generated from the reference model's own constants, and the
references to check an implementation against (docs/EMT250_program.md §15).

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_spec.py

Writes analysis/e250/e_plugin_spec.json (every constant with its confidence label and source section) and
analysis/e250/ref/e_static_rt<RT>.wav: Program E's static network (moving reads off), device rate (26,373.6 Hz), float32,
channels (L, R), R's direct tap at index 64 (K0), L's output kernel only (chain 'device': the plugin adds the 246's analogue chain; the unit tap reads 0.77 / 0.23); an
implementation of the static network must null these sample for sample. The moving part is time variant: check it
statistically against the spec's 'targets' with e_sim.stats on renders at several seeds.
"""
import json, os, sys
import numpy as np
import soundfile as sf

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import e_model as E  # noqa: E402
import e_sim as S  # noqa: E402

REF_RTS = (0.4, 1.0, 2.5, 4.5)


def spec():
    x = S.static_params(); q = E.unpack(x)
    lab = lambda v, c, src: dict(value=v, label=c, source=src)
    return {
        "program": "E (the EMT 250 program on the 246), docs/EMT250_program.md §11-15",
        "device_rate_hz": E.FS, "delays": "integer device samples; positions after R's first direct tap",
        "direct_path": lab(dict(taps=list(E.TAPS), gains=[1.0, float(q["a"][0]), float(q["a"][1])],
                                L=dict(delay=E.LDEL, sign=-1, gain=float(q["gL"]))), "S", "§11.2, held out by RT"),
        "allpasses": lab(dict(form="H(z) = (-g + z^-D) / (1 - g z^-D), in series, before the loops", D=list(E.AP_D),
                              g=[float(v) for v in q["g"]]), "S", "§11.2"),
        "loops": lab(dict(P=list(E.P), feedback="negative: line[t] = G{u[t] - line[t - P]}", static_read_delay=list(E.DK),
                          static_read_gain=[float(v) for v in q["c"]], read="after the loop filter (reads the line)",
                          L=dict(delay=E.LREV, sign=+1, gain="as the direct path's L gain")), "S", "§11.2"),
        "loop_filter": lab(dict(law="the 246's own: tools/ml_lines.line_coefs(P, rt, lf, hf) / line_resp (shelf + HF pole, HF floor)",
                                altiverb_family=dict(lf=float(q["m"][0]), hf=float(q["m"][1])),
                                high_band_uncertainty="+-10 % (original capture method unknown, §13.6)"), "S", "§11.2, §13.5-13.6"),
        "moving_reads": lab(dict(adopted={k: v for k, v in S.ADOPTED.items()},
                                 meaning="per loop M reads on the loop line, gain +-b x |static_read_gain|, centre uniform in [lo, hi] "
                                         "samples, offset A x triangle(f t + phase), phase independent per read, rate x (1 +- fsd) per read; "
                                         "a fraction 'share' also added to L with the opposite sign (L delay mdel), the rest independent in L; "
                                         "integer addresses, no interpolation",
                                 lfo_rate_hz=S.LFO_RATE,
                                 plugin_faithful=dict(reads=[dict(centre=c, loop="ABC"[k], R_sign=sr, L_sign=sl, set=st) for c, k, sr, sl, st in S.PLUGIN_READS],
                                                      set_probability=S.SETS, common_triangle=dict(f=0.5535, A=35.2), gain=0.17,
                                                      weak=dict(per_loop=6, gain=0.28, lo=300, hi=3600), dwell_s=0.1,
                                                      note="§14.4-14.5: the UAD plugin's structure; fits E's statistics worse (250 vs 182), a listening alternative")),
                             "S (energy, range, speed, L/R share: E's statistics) / H (triangle and rate: from the UAD plugin)", "§12-14"),
        "output_kernel": lab(dict(fir=list(S.L_KERNEL), rate="device", note="L's output kernel; E has it, the LP lever does not cause it"), "S", "§14.1"),
        "controls": {
            "RT": lab(dict(measured_steps=[0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5], note="the Altiverb family's steps; the 246's E table itself is unknown"), "C (these values exist) / H (full table)", "§0"),
            "LF_HF": lab("the 246's LF x / HF x controls through the loop law", "S", "§11.2"),
            "AMPL": lab(dict(law="direct taps x AMPL / 80", reference="the Altiverb family taken as AMPL 80"), "H", "§15 (S's linear cluster law, model_s_results)"),
            "DELAY": lab("a pure delay of the whole output, 0-999 ms, device samples; the Altiverb family at 0", "H", "§15, operation card"),
            "LOW_PASS": lab(dict(steps_khz=[2, 3, 4, 5, 6, 8, 10, 12], open=12, guess="2nd-order Butterworth at the device rate below 12 kHz"), "H", "§14.1, §15"),
        },
        "not_resolved": ["the L-only detail of E's moving reads (§14.3)", "E's own LFO rate and phase grouping (needs E repeat captures)",
                         "the LP lever's response", "AMPL / DELAY in E"],
        "references": dict(static_ir=[f"analysis/e250/ref/e_static_rt{r}.wav" for r in REF_RTS],
                           static_ir_note="render_ir(rt, moving=False, chain='device'), lf / hf the family's, AMPL 80, DELAY 0, LP open; device rate, (L, R)",
                           targets=json.load(open(os.path.join(E.OUT, "rand_stats_data.json"))),
                           targets_note="e_sim.stats of the eight RT >= 1.0 Altiverb presets (leave-one-out residuals); a plugin's E, rendered at "
                                        "those RTs through e_sim.capture's sweep method and the same function, should land within ~1.3-2 x the "
                                        "data's own half-vs-half cost (100)"),
    }


def main():
    sp = spec(); os.makedirs(os.path.join(E.OUT, "ref"), exist_ok=True)
    json.dump(sp, open(os.path.join(E.OUT, "e_plugin_spec.json"), "w"), indent=1)
    for rt in REF_RTS:
        y = S.render_ir(rt, moving=False, chain="device", secs=min(3 * rt + 0.5, 4.0)).astype(np.float32)
        sf.write(os.path.join(E.OUT, "ref", f"e_static_rt{rt}.wav"), y, int(round(E.FS)), subtype="FLOAT")
        print(f"RT {rt}: {y.shape}, tap {y[S.K0, 1]:+.3f}")
    print("wrote analysis/e250/e_plugin_spec.json")


if __name__ == "__main__":
    main()
