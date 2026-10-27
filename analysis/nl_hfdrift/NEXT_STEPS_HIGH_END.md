# High end above 12.3 kHz, all programs: hand-over and close-out (28 Oct 2026)

Read first: `analysis/nl_hfdrift/REPORT.md` (the whole N investigation), `docs/NonLin_program.md` §11.3, §11.6, `analysis/nl_hfdrift/h3_harm.log`. Rules as always: no hardware, 96 kHz Tomkin set only, label claims C / S / H, held-out numbers with excess over the noise floor and the parameter count, every run under `tools/memguard.py` (cap about 4-6 GB, one heavy process), smoke test small first, score **20 Hz - 13 kHz per band**, never only the 100 Hz - 10 kHz bench band, strict rule on every change.

## Close-out (28 Oct 2026)

**Decision (user): accept the coverage the data supports; no new model term for the high end.** S is at the cross-take repeat level through 12.8 kHz; M/L are scored to 12.3 kHz and 12.3–12.8 kHz stays declined (plugin roll-off 12.6 kHz); E above 12.4 kHz is the Altiverb capture chain's analogue band limit, and E's 20–100 Hz residual is the same capture chain. The audit below closes the remaining hand-over items without a model change.

- **The only default-kernel resampler in the model paths is `core.to_device_grid`.** Its `resample_poly(25, 91)` Kaiser(5) droop is real (−0.22 / −0.90 / −1.65 / −2.33 / −3.65 dB at 11.5 / 12 / 12.3 / 12.5 / 12.8 kHz) but **fixed**: the tone test through the function gives the same gain in early and late windows and at fractional positions 0 / 0.25 / 0.5 / 0.75 (`hf_audit.py`). It has no `nl_late.unwarp`-style time drift. N compensates it with `nl_structure.resampler_gain`; S/M/L's full-band model carries it explicitly through `sml_fullband.pipeline()`; E does not use this path (`e250.e_model.to_device` is a 128-tap Kaiser sinc after a 12.4 kHz FIR, flat to 12.3 kHz and drift-free at 42.8 ppm). No model changes.
- **S**: the wide-window deconvolution is at the cross-take repeat level in 12.0–12.3 and 12.3–12.8 kHz (e.g. 12.3–12.8 kHz 50–150 ms model −35.7 / −33.1 dB against repeat −34.5 / −34.6 dB, R / L). Above 12.8 kHz the unit's chain is 43–72 dB down. Closed.
- **M / L**: 20 Hz–12.3 kHz is scored and the top-octave kernel fix is adopted (M −57 / −54, L −54 / −46 dB at 10–11 / 11–12.3 kHz). 12.3–12.8 kHz is not scored by decision; the plugin's LP rolls off at 12.6 kHz and the unit is 43–72 dB down above 12.8 kHz. Closed by decision.
- **E**: above 12.4 kHz the captures carry only the analogue band limit (the audit's 12.4–13.1 kHz band is about −43 dB re 1–4 kHz; `EMT250_program` §16.1 gives −46…−48 dB of the energy), left to the plugin's chain. 10–12.4 kHz decay / EDT / late decay are inside the realisation spread; the known C50/Ts band-edge residual and the +1 dB late 10 kHz level are the tested optimum (`EMT250_program` §17.7). Low band: 20–100 Hz reverb is at the data's repeatability; the direct path is capture-chain limited (`e_model.pretap`; §16.4 / §16.7), with R's C50/C80/D50/Ts at 20–100 and 31.5 Hz the only outliers. No model change.
- **Not run, and no longer needed:** the S/M/L chain-term test from the cluster taps (S has no headroom, M/L 12.3–12.8 is declined, and the term is measurement-side, so it cannot improve the plugin); widening beyond 12.8 kHz (chain 43–72 dB down); the N-only per-capture freedom / soft shrinkage (out of scope).

Evidence: `analysis/nl_hfdrift/hf_audit.py` / `hf_audit.log` (cheap, read-only; no heavy run), `analysis/fullband/report.json` (`wide`), `analysis/e250/e_pulse_v6.json`, `docs/EMT250_program.md` §16–18.

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
| S, M, L | 11.0-12.3 kHz raised cosine (harness default) | model-side clock (eps in the model; not `unwarp`) | **closed:** S at the cross-take repeat level through 12.8 kHz (wide window); M/L scored to 12.3 kHz, 12.3-12.8 kHz declined (plugin roll-off 12.6); the `to_device_grid` droop is carried by `pipeline()` |
| E (250) | own path `e_model.to_device`: 12.4 kHz low-pass + 128-tap windowed sinc | own, not affected by the default-kernel artefact | **closed:** flat to 12.3 kHz, 12.4 kHz band edge; above it the capture chain's analogue limit; 20-100 Hz reverb at the data's repeatability |
| Altiverb / CLA captures | chain kernels of their own | various | same resampler filter in the IRs |

## Proposed tasks, in order (closed 28 Oct 2026; kept for the record)

1. **Audit every default-kernel resampler / fractional-position use.** **Done:** `hf_audit.py` / `hf_audit.log`. Only `core.to_device_grid` is in a model path; its droop is fixed and independent of the fractional position. N compensates it (`resampler_gain`), S/M/L carry it through `pipeline()`, E uses its own sinc. Analysis-only uses (`altiverb_*`, `chain_kernel`, `crossrate`, `synth`, `static_fit`) do not feed the plugin.
2. **S / M / L: same chain-term test.** **Closed by decision, not run.** S is at the repeat level through 12.8 kHz, so there is no headroom; M/L 12.3-12.8 kHz is declined; the term is measurement-side and would not improve the plugin.
3. **Widen the band.** **Closed:** S is scored to 12.8 kHz; above 12.8 kHz the chain is 43-72 dB down (the raw files' >13.2 kHz content is the DAC image band at -72...-85 dB). M/L 12.3-12.8 kHz declined; E's 12.4 kHz edge is its own filter and above it is the capture chain.
4. **Per-capture freedom beyond one scalar (N, 2-3.5 dB, fork G).** N-only; out of this scope. Left in `NonLin_program.md` §11.6.
5. **Soft shrinkage of the chain term.** N-only; out of this scope.
6. **LOW PASS.** **Closed by the N work:** the 12 kHz "LOW PASS" was the harness resampler, and the unit is flat as far as the library shows. The plugin's LP is a panel control, not a fitted 12 kHz unit roll-off.
7. **Plugin check (read-only, one afternoon).** **Done** (`REPORT.md` follow-up): `plugin/src/dsp/CaptureEq.h` is the low-band CLA transfer, nothing in `plugin/` reads `out_eq.json`, and the N spec's `output_eq` is none. The M/L output kernel is a port caveat only (it was fitted on caches that carry the harness response; do not double-count it in the plugin chain).

## Acceptance for the session

- Each artefact fix: the synthetic test is clean to +-0.001 samples and 0.01 dB in 11-13 kHz, and the program's bench does not get worse by more than 0.5 dB in any channel-window (strict rule), reported against the committed bench.
- Each model term: held out by capture (or by run), parameter count stated, excess over the noise floor stated, `docs/EMT246_unit_facts.md` updated.
- Keep the **Left open** list of `analysis/nl_hfdrift/REPORT.md` current.

## Tools left behind (all in `analysis/nl_hfdrift/`)

`common.py` (`pair`, `local_delay`, `unwarp2`), `unwarp_test.py`, `predict_drift.py`, `heldout.py`, `resid_eq.py`, `compare_bench.py` (two bench folders, strict rule per window and band), `d_times.py` / `d_front_all.py` (WAV timestamps, front-tap residuals), `f_table.py`, `f_table2.py` (the chain-term tables), `h_est.py`, `h2_est.py`, `h3_harm.py` (estimators and harm distribution), `a_repeat.py`, `a_noise.py` (ceilings), `b_groups.py`, `c_front.py`, `g_oracle.py`, `e_align.py`, `hf_audit.py` / `hf_audit.log` (the close-out audit).
