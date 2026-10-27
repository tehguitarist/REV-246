# High end above 12.3 kHz, all programs: hand-over for a dedicated session (28 Oct 2026)

Read first: `analysis/nl_hfdrift/REPORT.md` (the whole N investigation), `docs/NonLin_program.md` §11.3, §11.6, `analysis/nl_hfdrift/h3_harm.log`. Rules as always: no hardware, 96 kHz Tomkin set only, label claims C / S / H, held-out numbers with excess over the noise floor and the parameter count, every run under `tools/memguard.py` (cap about 4-6 GB, one heavy process), smoke test small first, score **20 Hz - 13 kHz per band**, never only the 100 Hz - 10 kHz bench band, strict rule on every change.

## What this session established (carry forward)

1. **Two harness artefacts, both N-proven, probably everywhere a default-kernel resampler meets a fractional position:**
   - `nl_late.unwarp` (clock correction): scipy's default Kaiser(5) interpolator (3.9 kHz transition band around the device Nyquist 13.19 kHz) gave a delay and amplitude error at 11-12.8 kHz that grew with file time. Fixed (long kernel, default).
   - `core.to_device_grid`'s `resample_poly(25, 91)` carries the **same filter response in every cached IR of every program**: -0.22 / -0.90 / -1.50 / -2.33 / -3.65 dB at 11.5 / 12 / 12.25 / 12.5 / 12.8 kHz (`nl_structure.resampler_gain`). N now models it exactly. S / M / L / Altiverb chain kernels were fitted on the same cached IRs, so they **absorb** it (H, never checked): their high-end shape is partly the harness's, and it must not be double counted if a stage is ever taken from them for the plugin.
2. **What is left in N above 11 kHz is a term of the capture chain, per capture run** (S): a rank-1 band-edge filter shift (9 complex knots, 10.5-12.8 kHz) with one scalar per capture run (12 runs from the WAV `bext` timestamps, 9-12 April 2016), estimated from the run's front taps only; applied by default (`NL_HF`, `analysis/nl/hf_corr_g5.json`, built by `analysis/nl_hfdrift/f_table2.py`). 10-12 kHz -40.8 / -39.9 -> -48.3 / -49.0 dB; 4.8 % of channel-bands worse than the uncorrected model by > 0.5 dB (typical +1.9 dB, still near -46 dB), nothing below 10 kHz moves.
3. **Ceilings (N, 576 LATE-only pairs, alignment fitted on 8-9 kHz so each sub-band is held out):** repeat null 9-11 / 11-11.6 / 11.6-12.3 / 12.3-12.8 kHz = -56.5 / -56.0 / -52.4 / -40.4 dB. 12.3-12.8 kHz is repeat-limited near -40 dB (model best case about -43 dB). Tail noise is 62-76 dB down, so noise is not the limit below 12.3 kHz.
4. **Alignment is not the lever:** refitting o / front shift / ppm on a high band costs 5-10 kHz 4-9 dB for <= 2 dB; the bench alignment is optimal to 1e-4 samples.
5. The **LOW PASS** of N shows no roll-off of its own at its 12 kHz step; the Tomkin library has **no LP or DELAY steps** (`docs/library_notes.md`: the planned sets show as empty), so the filter's shape below 12 kHz is untested in the library.

## Per-program status of the high end (what the code does today)

| Program | Analysis band window (`core.F_ROLL` / `F_HI`) | Clock correction | High-end state |
|---|---|---|---|
| N (Non-Lin) | widened to 12.0-12.8 kHz (`NL_WIDE`) | long-kernel `unwarp` + exact resampler gain + per-run chain term | done to the -48 dB level up to 12 kHz; 12.3-12.8 kHz around -30 dB, ceiling about -43 |
| S, M, L | 11.0-12.3 kHz raised cosine (harness default) | model-side clock (eps in the model; not `unwarp`) | not examined for these artefacts; nothing scored above 12.3 kHz |
| E (250) | own path `e_model.to_device`: 12.4 kHz low-pass + 128-tap windowed sinc | own, not affected by the default-kernel artefact | HF bands are a known open item in `EMT250_program.md` |
| Altiverb / CLA captures | chain kernels of their own | various | same resampler filter in the IRs |

## Proposed tasks, in order (each starts with a small run)

1. **Audit every default-kernel resampler / fractional-position use.** `grep` hits for `resample_poly` / `np.interp(t`: `altiverb_compare.py`, `altiverb_predict.py`, `chain_kernel.py`, `harness/core.py` (`to_device_grid`, `deconvolve`), `harness/crossrate.py`, `harness/synth.py`, `e250/e_sim.py`, `e250/static_fit.py`. Template test (minutes, no data): `analysis/nl_hfdrift/unwarp_test.py` - push a model through the function at a fractional position that advances, and measure the delay and amplitude of the 11-13 kHz bands against 6-8 kHz. Any function that fails gets the long kernel (`common.unwarp2`) behind a switch, as N did.
2. **S / M / L: same chain-term test.** Their cluster taps (four per channel, `config.json` cluster templates) play the role of N's front taps. Estimate the capture chain's complex residual above 11 kHz from the cluster alone per capture, group by capture run (WAV `bext` timestamps; check that S / M / L files carry them and whether the runs are the same days), and test whether the same rank-1 shape fits. Held out on the reverb windows, with the S / M / L chain3 kernel and the harness band window as they are. Expectation (H): the same effect, because the capture chain is common.
3. **Widen the band.** The analysis window ends at 12.3 kHz (S / M / L) and 12.8 kHz (N) while the device Nyquist is 13.19 kHz. Rebuild the IR caches (own cache directories, never over the shared ones) with a resampler that is flat to about 13.0 kHz and a window out to about 13.1 kHz, check the noise floor and the unit's chain level there (the chain is 26 dB down at 12.5 kHz, so little energy), and see whether anything deterministic is above 12.8 kHz. Also look at the raw 96 kHz files above 13.2 kHz (-72...-85 dB per the notes): the DAC's image band.
4. **Per-capture freedom beyond one scalar (N, 2-3.5 dB, fork G):** a free complex gain on 3 or 9 knots per capture gains that much beyond the scalar. Find what it is: a second shape component per run, a late-section term (most failures have LATE > 0), the within-run drift of runs 7 and 11. Test whether a rank-2 shape with a second front-tap scalar helps held out.
5. **Soft shrinkage of the chain term** (c = 0.005 caps the worst case at +1.9 dB for about 1.3 dB less median gain; `h3_harm.py`), if a full bench says the 4.8 % tail matters.
6. **LOW PASS** (needs the unit's description, no data): with no LP steps in the library, derive the plugin's LP from the manual / schematic (12 kHz = flat as far as the library shows), not from the captures. Record this in `EMT246_unit_facts.md` if not already.
7. **Plugin check (read-only, one afternoon):** make sure nothing in `plugin/` takes a high-end shape from the old fitted EQ (`out_eq.json`) or from S / M / L kernels; the N spec's `output_eq` is now "none".

## Acceptance for the session

- Each artefact fix: the synthetic test is clean to +-0.001 samples and 0.01 dB in 11-13 kHz, and the program's bench does not get worse by more than 0.5 dB in any channel-window (strict rule), reported against the committed bench.
- Each model term: held out by capture (or by run), parameter count stated, excess over the noise floor stated, `docs/EMT246_unit_facts.md` updated.
- Keep the **Left open** list of `analysis/nl_hfdrift/REPORT.md` current.

## Tools left behind (all in `analysis/nl_hfdrift/`)

`common.py` (`pair`, `local_delay`, `unwarp2`), `unwarp_test.py`, `predict_drift.py`, `heldout.py`, `resid_eq.py`, `compare_bench.py` (two bench folders, strict rule per window and band), `d_times.py` / `d_front_all.py` (WAV timestamps, front-tap residuals), `f_table.py`, `f_table2.py` (the chain-term tables), `h_est.py`, `h2_est.py`, `h3_harm.py` (estimators and harm distribution), `a_repeat.py`, `a_noise.py` (ceilings), `b_groups.py`, `c_front.py`, `g_oracle.py`, `e_align.py`.
