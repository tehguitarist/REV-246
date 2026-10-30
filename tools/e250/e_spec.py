"""Program E for the plugin: the specification as data, generated from the reference model's own constants, and the
references to check an implementation against (docs/EMT250_program.md §15).

    .venv/bin/python tools/memguard.py --limit-gb 3 -- .venv/bin/python tools/e250/e_spec.py
    .venv/bin/python tools/e250/e_spec.py --embed     # only add loop_law / law_vectors / moving_reads_draw to the committed spec

Writes analysis/e250/e_plugin_spec.json (every constant with its confidence label and source section) and
analysis/e250/ref/e_static_rt<RT>.wav: Program E's static network (moving reads off), device rate (26,373.6 Hz), float32,
channels (L, R), R's direct tap at index 64 (K0), L's output kernel [0.77, 0.23] applied to both channels in the
'device' chain (the plugin adds the 246's analogue chain); an implementation of the static network must null these
sample for sample. The moving part is time variant: check it statistically against the spec's 'targets' with
e_sim.stats on renders at several seeds.
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
    static = json.load(open(os.path.join(E.OUT, "static_all_246.json")))
    if "rev_gain" not in static:
        raise RuntimeError(f"{E.OUT}/static_all_246.json has no rev_gain; regenerate it with the adopted E fit")
    lab = lambda v, c, src: dict(value=v, label=c, source=src)
    return {
        "program": "E (the EMT 250 program on the 246), docs/EMT250_program.md §11-17",
        "device_rate_hz": E.FS, "k0": E.K0,
        "delays": "integer device samples; positions after R's first direct tap at index k0",
        "rev_gain_db": lab(float(static["rev_gain"]),
                           "the reverb's level against the direct tap; already included in loops.static_read_gain and "
                           "moving_reads.reads amplitudes - do not apply it again to either",
                           "§17.6-17.7"),
        "direct_path": lab(dict(taps=list(E.TAPS), gains=[1.0, float(q["a"][0]), float(q["a"][1])],
                                L=dict(delay=E.LDEL, sign=-1, gain=float(q["gL"]))), "S", "§11.2, held out by RT"),
        "allpasses": lab(dict(form="H(z) = (-g + z^-D) / (1 - g z^-D), in series, before the loops", D=list(E.AP_D),
                              g=[float(v) for v in q["g"]]), "S", "§11.2"),
        "loops": lab(dict(P=list(E.P), feedback="negative: v[t] = u[t] - line[t - P], line = G{v}", static_read_delay=list(E.DK),
                          static_read_gain=[float(v) for v in q["c"]], read="after the loop filter (reads the line)",
                          L=dict(delay=E.LREV, sign=+1, gain="as the direct path's L gain")), "S", "§11.2"),
        "loop_filter": lab(dict(law="the 246's own: tools/ml_lines.line_coefs(P, rt, lf, hf) then line_resp (shelf + HF pole, "
                                    "HF floor); the panel's LF / HF multiply the loop law",
                                altiverb_family=dict(lf=float(q["m"][0]), hf=float(q["m"][1]),
                                                     note="adopted LF 1.49 / HF 0.69; LF from the 63-500 Hz decay (§16.5), "
                                                          "HF from the 2-12.4 kHz decay (§17.3), held out. Earlier §11.2 / §12.4 / "
                                                          "§13.6 / §15 values (1.44 / 0.72 and 1.38 / 0.74) are superseded"),
                                decay_error_held_out="0.6-2.1 % per octave, 63 Hz-12.4 kHz (§17)"), "S", "§11.2, §16.5, §17.3"),
        "moving_reads": lab(dict(adopted={k: v for k, v in S.ADOPTED.items()}, reads=S.fixed_reads(),
                                 meaning="nine fixed reads on the loop lines, taken BEFORE the loop filter (v[t] = u[t] - line[t - P], "
                                         "line = G{v}; the static reads take line, after G; §17) (centre = output delay after R's direct tap, loop, signed R and L "
                                         "amplitudes, direct tap = 1), each read at centre + A x triangle(f_j t + phase_j), triangle in [-1, 1], "
                                         "phase independent per read and realisation, rate f_j = f x (1 +- fsd) per read; L reads the same "
                                         "address with its own amplitude, mdel samples later, or mdel + dl_same for the reads whose L and R signs agree (§16.7); "
                                         "integer addresses, no interpolation",
                                 triangle="wav = 4 abs(frac(f_j t + phase_j) - 0.5) - 1, in [-1, 1]",
                                 speed_spread="f_j = f x U(1 - fsd, 1 + fsd) per read and realisation (fsd 0.5); "
                                              "phases are independent (common false)",
                                 one_pole="after summing the reads and before L's offset: y[n] = (1 - mlp) x[n] + mlp y[n-1], "
                                          "mlp = 0.07, applied to both the R and L moving sums",
                                 l_gain="reads[].L is the final L output amplitude: do NOT multiply L's moving path by gL "
                                        "again. The 1.114 gL applies to the direct path and the static reverb only; the L "
                                        "moving offset is mdel (314), or 315 for reads whose R and L signs agree",
                                 provenance="the static reads and fixed_reads.json were fitted with render_mean, which omits "
                                            "dl_same and the mlp one-pole; validate the moving part with e_sim.render, not "
                                            "render_mean",
                                 superseded_random_reads={k: v for k, v in S.ADOPTED_13.items()},
                                 lfo_rate_hz=S.LFO_RATE,
                                 plugin_faithful=dict(reads=[dict(centre=c, loop="ABC"[k], R_sign=sr, L_sign=sl, set=st) for c, k, sr, sl, st in S.PLUGIN_READS],
                                                      set_probability=S.SETS, common_triangle=dict(f=0.5535, A=35.2), gain=0.17,
                                                      weak=dict(per_loop=6, gain=0.28, lo=300, hi=3600, share=0.5),
                                                      filter_position="the plugin-faithful reads take the loop line AFTER the loop filter; "
                                                                      "the adopted fixed reads take v before it",
                                                      lfo="the common triangle's 0.5535 Hz holds with lfo='time'; lfo='samples' adds 0.0005",
                                                      dwell_s=0.1,
                                                      note="§14.4-14.5: the UAD plugin's structure; a listening alternative. Its old "
                                                           "comparison (250 vs 182) was against the superseded random-read law; the adopted "
                                                           "fixed-read law is about 86-91")),
                             "S (reads, amplitudes, depth: E's low band, held out; speed spread: E's statistics) / H (triangle shape and rate: the UAD plugin)", "§16-17"),
        "output_kernel": lab(dict(fir=list(S.L_KERNEL), taps_offset=16, rate="device",
                                  note="L's output kernel, applied to BOTH output channels in the 'device' chain at sample "
                                       "offset 16 (taps 16 and 17); the static reference IRs compensate this 16-sample delay "
                                       "so R's direct tap is at k0. E has it, the LP lever does not cause it"),
                             "S", "§14.1"),
        "capture_chain_tail": lab("the Altiverb captures' chain has a slow tail (analysis/e250/kernel_tail.json, B-spline 540-4000 samples); a property of the capture chain, not of E: the plugin's own analogue chain replaces it", "S", "§16"),
        "controls": {
            "RT": lab(dict(measured_steps=[0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5], note="the Altiverb family's steps; the 246's E table itself is unknown"), "C (these values exist) / H (full table)", "§0"),
            "LF_HF": lab("the 246's LF x / HF x controls through the loop law; adopted family values LF 1.49 / HF 0.69", "S (shape) / H (same code)", "§11.2, §16.5, §17.3"),
            "AMPL": lab(dict(law="direct taps x AMPL / 80", reference="the Altiverb family taken as AMPL 80"), "H", "§15 (S's linear cluster law, model_s_results)"),
            "DELAY": lab("a pure delay of the whole output, 0-999 ms, device samples; the Altiverb family at 0", "H", "§15, operation card"),
            "onset": lab("at DELAY 0 E's R direct tap arrives with M's first cluster tap (Altiverb: S 10.17 / M 20.33 / L 40.42 / E 20.33 ms, one capture chain)", "S", "§17.8"),
            "level": lab("E's direct tap 1 at the level of an S / M / L cluster tap at the same AMPL (the Altiverb normalised every preset to equal total energy, 64.3-64.6 dB, so the absolute program levels are not in the data)", "H (calibration guess)", "§17.8"),
            "LOW_PASS": lab(dict(steps_khz=[2, 3, 4, 5, 6, 8, 10, 12], open=12, guess="2nd-order Butterworth at the device rate below 12 kHz"), "H", "§14.1, §15"),
        },
        "capture_pretap": lab("the Altiverb captures carry a deterministic, non-causal low-band component before the direct tap (deconvolution): e_model.pretap, for comparisons with captures only", "S", "§16.7"),
        "not_resolved": ["E's LFO rate within 0.5-0.65 Hz (E's own statistics; 0.553 = the 250's in seconds, adopted, or 0.608 = in samples) and the triangle shape (the UAD plugin; E cannot tell triangle from sine)",
                         "the LP lever's response", "AMPL / DELAY in E"],
        "references": dict(static_ir=[f"analysis/e250/ref/e_static_rt{r}.wav" for r in REF_RTS],
                           wav_sample_rate_hz=int(round(E.FS)),
                           static_ir_note="render_ir(rt, moving=False, chain='device'), adopted LF 1.49 / HF 0.69, AMPL 80, DELAY 0, LP open; "
                                          "device rate (the exact device_rate_hz; the WAV header is the integer wav_sample_rate_hz), "
                                          "(L, R), R's tap at index k0; lengths int(secs*FS) with secs = min(3 RT + 0.5, 4.0), so "
                                          "the RT 4.5 reference is truncated at 0.89 RT",
                           targets=json.load(open(os.path.join(E.OUT, "rand_stats_data.json"))),
                           targets_note="e_sim.stats of the eight RT >= 1.0 Altiverb presets; the adopted fixed-read generator "
                                       "lands at about 86-91 (lower is better) against the data's own half-vs-half cost of 100. "
                                       "lobes and lr use leave-one-out residuals; e_ratio and frac use the static model residual"),
    }


LINE_UNIT = 1321          # ml_lines.line_coefs: the HF floor's program-time factor round(min(P) / 1321)
LAW_VECTORS = [           # (RT, LF, HF): the measured steps at the family values, then off-step panel values
    *[(rt, 1.49, 0.69) for rt in (0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5)],
    (0.55, 1.49, 0.69), (5.0, 1.49, 0.69), (1.0, 1.0, 0.9), (1.0, 0.1, 0.2), (0.4, 2.0, 0.1), (2.2, 0.4, 1.0),
    (3.3, 0.7, 0.5), (5.0, 2.0, 1.0), (0.4, 0.1, 0.05), (1.7, 1.2, 1.2), (4.8, 0.25, 0.3), (0.8, 3.0, 0.8),
]


def embed_codegen_data(sp):
    """Copies what tools/plugin_codegen/e.py needs from outside this spec into it, so the generator reads the spec alone.
    Adds loop_law (the 246's loop law that loop_filter.law names: analysis/q1/laws.json's LF table, shelf pole and f_ref,
    ml_lines' coefficient limits and drop-back, the HF pole curve of analysis/lowband/chain3_pole_cells_allf.json),
    law_vectors (ml_lines.line_coefs per loop at the (RT, LF, HF) test points) and moving_reads_draw (e_sim.e_reads' draws
    at seed 0: rate and start phase per read). Nothing else in the spec changes. Also run alone against the committed
    spec: tools/e250/e_spec.py --embed."""
    for k in ("ML_POLE_CURVE", "ML_HFRT", "ML_L22"):
        if os.environ.get(k):
            raise SystemExit(f"{k} is set: the embedded data is the adopted law only")
    sys.path.insert(0, os.path.join(os.path.dirname(HERE)))
    import ml_lines as mll
    lab = lambda v, c, src: dict(value=v, label=c, source=src)
    law = mll.LAW
    pc = json.load(open(os.path.join(os.path.dirname(os.path.dirname(HERE)), "analysis", "lowband", "chain3_pole_cells_allf.json")))
    sp["loop_law"] = lab(dict(a=law["a"], f_ref=law["f_ref"], lf={k: {n: law["lf"][k][n] for n in ("hi", "dc", "k")} for k in law["lf"]},
                              p_max=mll.P_MAX, dc_min=mll.DC_MIN, beta=mll.BETA, hfrt_min=mll.HFRT_MIN, line_unit=LINE_UNIT,
                              pole_curve=dict(knots=pc["knots"], values=pc["all"], fold=pc["fold"], drop_back=mll.DROP_BACK)),
                         "S", "analysis/q1/laws.json (a, f_ref, lf), tools/ml_lines.py (P_MAX, DC_MIN, BETA, HFRT_MIN, DROP_BACK), "
                              "analysis/lowband/chain3_pole_cells_allf.json (knots, 'all', fold); §14.2, §16.15-16.16, §29, §32.6")
    out = [[float(v) for v in mll.line_coefs([P], rt, lf, hf)[0]] for rt, lf, hf in LAW_VECTORS for P in sp["loops"]["value"]["P"]]
    sp["law_vectors"] = lab(dict(inputs=[list(v) for v in LAW_VECTORS], loops_P=sp["loops"]["value"]["P"], outputs=out,
                                 note="outputs[v * len(loops_P) + k] = (g, r0, p) of loop k at inputs[v] = (RT, LF, HF)"),
                            "S", "ml_lines.line_coefs([P], RT, LF, HF)[0] per loop; the plugin's test vectors")
    ad, reads = sp["moving_reads"]["value"]["adopted"], sp["moving_reads"]["value"]["reads"]
    rng = np.random.default_rng(0)
    rng.uniform()                                     # the common phase (unused: common = false)
    rates, phases = [], []
    for _ in reads:
        rates.append(float(ad["f"] * rng.uniform(1 - ad["fsd"], 1 + ad["fsd"])))
        phases.append(float(rng.uniform()))
    sp["moving_reads_draw"] = lab(dict(seed=0, rate_hz=rates, phase=phases,
                                       note="e_sim.e_reads' draws in its own order: one common phase, then per read rate f x U(1 - fsd, 1 + fsd) and a start phase"),
                                  "H (one fixed realisation)", "e_sim.e_reads, numpy default_rng(0); §16-17")
    return sp


def main():
    if "--embed" in sys.argv:       # add the generator's inputs to the committed spec without regenerating it
        path = os.path.join(E.OUT, "e_plugin_spec.json")
        sp = embed_codegen_data(json.load(open(path)))
        json.dump(sp, open(path, "w"), indent=1)
        print("embedded loop_law, law_vectors, moving_reads_draw in analysis/e250/e_plugin_spec.json")
        return
    sp = spec(); os.makedirs(os.path.join(E.OUT, "ref"), exist_ok=True)
    sp = embed_codegen_data(sp)
    json.dump(sp, open(os.path.join(E.OUT, "e_plugin_spec.json"), "w"), indent=1)
    for rt in REF_RTS:
        y = S.render_ir(rt, moving=False, chain="device", secs=min(3 * rt + 0.5, 4.0)).astype(np.float32)
        sf.write(os.path.join(E.OUT, "ref", f"e_static_rt{rt}.wav"), y, int(round(E.FS)), subtype="FLOAT")
        print(f"RT {rt}: {y.shape}, tap {y[S.K0, 1]:+.3f}")
    print("wrote analysis/e250/e_plugin_spec.json")


if __name__ == "__main__":
    main()
