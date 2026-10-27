# NonLin_program §11.6, 11–12.3 kHz timing drift: report (27 Oct 2026)

All scripts and outputs are in `analysis/nl_hfdrift/` (new files only; nothing else touched, nothing committed). Run under `tools/memguard.py` (peak 0.35 GB).

## Result

**The drift is a software artefact of the harness's own clock correction, `nl_late.unwarp` (S).** It is not the unit, not Tomkin's capture chain, not aliasing, not the clock ratio's precision.

`unwarp(x, ppm)` upsamples ×32 with scipy's default `resample_poly` Kaiser(5) FIR (10 input samples per side), then takes samples at `n(1+ppm·1e-6)·32` with `np.interp`. That FIR has a transition band of about 3.9 kHz around the 13.19 kHz Nyquist (gain −0.22 / −0.90 / −1.65 / −2.72 / −3.65 dB at 11.5 / 12 / 12.3 / 12.6 / 12.8 kHz), so the image of the 11–13 kHz content is not suppressed. Resampling at a fractional position δ then gives `G(f) + G(fs−f)·e^{i2πδ}` instead of a pure delay. δ = ε·n grows with file time, so the band edge gets a phase error ∝ t and an amplitude error ∝ t², then turns round when δ passes about 0.2 (about 0.9 s at 12 ppm). Below 11 kHz G(f)+G(fs−f)=1 and nothing happens.

## Evidence

1. **Reproduced** (capture vs model, bench's own o / ppm / gain, full band, 4096-sample windows). Key (0, 10, 2.0), ppm 12.75, delay of y after m in samples, t = 0.08 / 0.24 / 0.39 / 0.55 / 0.70 / 0.86 / 1.17 s:
   - 11.0–11.6 kHz: +0.001 .. +0.005
   - 11.6–12.3 kHz: +0.010 / 0.020 / 0.025 / 0.033 / 0.037 / 0.043 / 0.026
   - 12.3–12.8 kHz: +0.024 / 0.044 / 0.062 / 0.077 / 0.093 / 0.091 / 0.078
   - ≤ 10.5 kHz: |τ| ≤ 0.001.
   - Same in L and R (±0.003) and in (0.5, 10, 2.0) and (0.3, 1.0, 1.0): independent of SPAN, EARLY, LATE and channel. Also with LATE 0, so not loop passes.
   - |Y/M| at 12.6 kHz falls 1.02 → 0.94 → 0.75 → 0.52 at 0.1 / 0.4 / 0.8 / 1.5 s.
   - Phase across 11.1–12.65 kHz at 0.4 s: −0.002 / −0.020 / −0.054 / −0.113 / −0.186 / −0.289 rad.
2. **File time, not response time.** In trimmed captures (SPAN 0.5, EARLY −0.2, file starts 0.34 s into the response) the drift starts near 0 at the file start (+0.008 at 0.06 s) and grows like the untrimmed ones.
3. **The model alone shows it.** `unwarp(model, ppm)` vs `model` (no capture at all): the 12.3–12.8 kHz band lags the 6–8 kHz band by +0.079 samples at 0.86 s, the 11.6–12.3 kHz band by +0.038 (`unwarp_test.py`). A kernel with 150 input samples per side and Kaiser β = 9 (`common.unwarp2`) gives a common delay in every band to ±0.001 and |Z/M| = 1.000 (`unwarp_test2.py`).
4. **Predicted vs measured** (`predict_drift.py`): measured minus predicted artefact is flat in time (key (0, 10, 2.0), 11.6–12.3 kHz: +0.004 / 0.005 / 0.005 / 0.005 / 0.005; 12.3–12.8 kHz: +0.013 / 0.012 / 0.012 / 0.010 / 0.010). The whole time-growing part is the artefact. What is left is a static offset (point 7).
5. **Fix, held out** (120 random captures, defects excluded, seed 7; the EQ is fitted on 60 and scored on the other 60; o / ppm / gain are the bench's own, 0.05–0.6 s, bench's band metric; `heldout.py`, `heldout120.log`). Median null L / R at 10–12 kHz:
   - bench as published (its EQ, old unwarp), same 60 captures: **−34.6 / −34.0**
   - old unwarp + pooled real EQ refit (12 numbers): −35.6 / −35.0
   - new unwarp, no EQ: −30.9 / −31.2
   - new unwarp + pooled real EQ (12 numbers): **−40.5 / −40.1**
   - new unwarp + pooled complex EQ (24 numbers): **−41.6 / −40.4**
   - Median 11.6–12.3 kHz delay of y after m, all 120 captures: 0.0219 → 0.0045 samples.
   - Sub-bands with new unwarp + complex EQ (test half, L / R): 8–10 kHz −53.4 / −52.9 (old unwarp + real EQ −52.9 / −51.7), 10–11 kHz −49.1 / −46.9 (unchanged), 11–11.6 −41.7 / −40.2, 11.6–12.3 −36.1 / −35.9, 12.3–12.8 −26.9 / −27.6.
   - After the fix the phase and amplitude are flat in time (key (0, 10, 2.0), 12.65 kHz: −0.062 / −0.054 / −0.052 rad and |Y/M| 0.948 / 0.943 / 0.954 at t = 0.1 / 0.4 / 0.8 s, against −0.153 / −0.289 / −0.389 and 0.935 / 0.857 / 0.689 before).
6. **The adopted "top-end EQ" (§11.3) is the harness's own resampler.** `to_device_grid`'s `resample_poly(25, 91)` uses the same kind of Kaiser(5) FIR (3.9 kHz transition). Its gain at 11.5 / 12 / 12.25 / 12.5 / 12.8 kHz is −0.22 / −0.90 / −1.50 / −2.33 / −3.65 dB. The EQ refitted after fixing unwarp is −0.18 / −0.86 / −1.36 / −1.92 / −2.56 dB: within 0.05 dB up to 12 kHz, +0.15 / +0.41 / +1.09 dB above (`eq_vs_resampler.py`). The old EQ (−0.24 / −1.03 / −1.85 / −3.02 at 11.5 / 12 / 12.3 / 12.6) is this plus the old unwarp's average droop. So **no LOW PASS at 12 kHz is needed (S for the numbers; the old H "LOW PASS" is not supported)**: the unit is flat there as far as this library shows.
7. **Static part (fixed filter, not drift).** The pooled complex EQ's phase on the fixed-unwarp data is 0 / −4 / −11 / −17 / −28 / −40 mrad at 11.5 / 11.75 / 12 / 12.25 / 12.5 / 12.8 kHz (delay-equivalent 0.003–0.013 samples). That is a fixed mismatch between the real chain and the chain kernel's regularised inverse (H). It is what limits 11.6–12.8 kHz now.

## Ruled out

- **Clock ratio precision** (`ppmfit.py`, 60 test captures): refitting ppm per capture on 8–11 kHz or 9–11 kHz moves it by sd 0.68 ppm, mean −0.05 / −0.12, and gains ≤ 0.5 dB at 9–12 kHz. Fitting it in 11–12.8 kHz itself (in-sample) gains 1.7 dB there but costs 3–8 dB at 1–10 kHz, correlation between shifts −0.08.
- **Aliasing, and a sweep or capture-clock effect**: no need; the artefact reproduces on the model with no capture in it, and the capture-minus-prediction residual is time-flat.
- **Repeat captures agree trivially**: both carry the same artefact. Repeat null (SPAN 0.5, EARLY 10 / 0.5 / 1 / −1, adjacent LATE pairs, 0.05–0.45 s, delay + gain fitted on 8–11 kHz, 20 pairs) is the same with old and new unwarp (10–12 kHz −53.2 / −48.9 dB, 11.8–12.8 kHz −46.6 / −46.3). My quick repeat ceiling is **not** the −58.9 dB quoted in §11 (different alignment / pairs); against it, the fixed model at −41.6 / −40.4 is about 8–12 dB short, not 25 dB.
- Not done: Altiverb N presets (not needed, the cause is on the harness side); S / M / L tools' clock handling (they apply ε in the model or on other paths; `e250/e_model.to_device` uses a 128-tap windowed sinc, so it does not have this); a full re-bench of 1,044 captures (the bench is running; I used its o / ppm / gain read-only).

## Caveats

- Held out by capture (the structure laws, o / ppm / gain are the bench's, fitted per capture). The EQ's 12 (real) or 24 (complex) numbers are fitted on 60 captures and scored on 60 others.
- The bench's ppm and o were fitted with the old unwarp in 100 Hz–10 kHz; the fix changes that band by < 1 dB, so the ppm fit is unaffected.
- I did not compute a noise floor for 11.6–12.8 kHz; the residual there is a mix of the static phase and noise, not separated.

## Recommended next step

1. In `tools/nl_late.py` replace `unwarp`'s upsampler by a long kernel (`analysis/nl_hfdrift/common.py:unwarp2`, 32×, half 150 input samples, Kaiser β = 9) or apply the clock on the model side; use it in `nl_model.py` (lines 203, 249, 279). Then replace the fitted `out_eq.json` by the exact `resample_poly(25, 91)` FIR gain (0 free numbers) plus, if wanted, the small complex residual above 11.5 kHz (12–24 numbers). Re-run the bench; expect 10–12 kHz about −40…−42 dB (was −34.5), other bands within 1 dB.
2. Docs: close §11.6 first bullet; §11.3 "top-end EQ" and "H: LOW PASS at its 12 kHz step" are a harness artefact. The plugin spec's top-end EQ (`CaptureEq`) should not be presented as a unit property.
3. What is left at 11.6–12.8 kHz (−27…−36 dB held out) is the static phase, 0.01 samples at the edge: next candidate is the regularised inverse (ε·max|G|² = 1e-3 dominates above about 12.1 kHz) against the real chain.

## Follow-up (27 Oct, second half): landed behind a switch, benched in full, then made the default

- `tools/nl_late.py`: `upsample` / `unwarp_from` / `unwarp` with the long kernel (`NL_WARP2`, now default on; `NL_WARP2=0` restores the old warp bit for bit: the legacy bench file for (0.2, 1.0, 0.4) reproduces to 0.0). The clock-independent upsampling is cached per capture in the bench.
- `tools/nl_structure.py`: `resampler_gain` (the exact `resample_poly(25, 91)` gain). `tools/nl_model.py`: with the new warp `_proc` multiplies it in instead of `out_eq.json` (0 fitted numbers); `NL_BK_FILE`, `BK_ALPHA_DEFAULT` 0.7 and bench folder names follow the warp (`bench_warp2*`).
- Full 1,044-capture benches (memguard 6 GB, 4 workers, peak 1.4 GB, about 34 min each): `bench_warp2` (no kernel), `bench_warp2_bk_a0.65` (old kernel, old strength), `bench_warp2_bk_kf` (refit kernel, α 1), `bench_warp2_bk_kf_a0.7` (adopted); subsets `..._a0.5`, `..._a0.6` (151 captures). Comparisons: `compare_bench.py`, outputs `compare_kf_*.txt`, reports `report_default_*.txt`.
- Old kernel on the new warp fails the strict rule (worst +1.56 dB vs `bench_bk_a0.65`): about a quarter of it was the old warp's interpolation error. Refit (`branch_kernel_warp2.json`, 74 / 79 % of the old amplitude): α 1 fails (worst +0.98), **α 0.7 holds vs `bench_bk_a0.65` (worst +0.07 dB at 0–0.6 s)** but loses up to +0.81 dB (22 channel-windows over 0.5 at 0.05–0.15 s) against the same bench without the kernel; α 0.5 holds that on the screened subset (worst +0.49).
- Complex residual above 11.5 kHz (22 numbers, `resid_eq.py`): +0.5…+2.3 dB at 10–12 kHz median, 30 % of captures worse by > 0.5 dB: not kept.
- Plugin (read-only check): `plugin/src/dsp/CaptureEq.h` is the low-band CLA transfer of `model_s_results.md` §27 (10 numbers, biquads), not a 12 kHz EQ; nothing in `plugin/` reads `output_eq` or `out_eq.json`. The spec's `output_eq` is now `null` / "none".

## Second follow-up (27 Oct): the 11-12.8 kHz residual that is left is the capture chain, per capture run (opt-in `NL_HF=1`)

Four parallel read-only checks (scripts and logs `a_*`, `b_*`, `c_*`, `d_*`, `e_*` here), then a full bench:

- **A, ceiling (576 LATE-only pairs, alignment fitted on 8-9 kHz only, so every sub-band held out):** repeat null 9-11 / 11-11.6 / 11.6-12.3 / 12.3-12.8 kHz = -56.5 / -56.0 / -52.4 / -40.4 dB against the model's -52.6 / -42.2 / -36.0 / -28.4 dB (L, before this step): 15-19 dB of headroom above 11 kHz, so it is **model-side, not noise or repeat-limited**. Noise is 62-76 dB down. The 10-12 kHz ceiling -55.8 dB (the -58.9 of §11 is about right; my earlier -53 / -49 dB quick estimate was low).
- **B, structure:** the static residual (phase 0 / -6 / -13 / -20 mrad at 11.5-12.25 kHz, +0.03-0.07 dB) is the same for every time window after 0.05 s, SPAN, LATE, channel and trimming: the chain / regularised inverse, not the model. Group-specific fits never beat the pooled one.
- **C, front taps:** the same residual is measured from the front taps alone (R: -1.5 / -5.7 / -12.1 / -20.8 / -30.3 / -40.0 mrad at 11.5 / 11.75 / 12 / 12.25 / 12.5 / 12.8 kHz against -2.7 / -7.2 / -14.9 / -22.4 / -31.6 / -32.1 mrad on the reverb); one pooled 20-number correction gains 2-4 dB at 10-12 kHz held out, 22 % of captures worse.
- **D, sessions:** the residual is a band-edge filter shift that follows the **capture run** (the WAV `bext` timestamps, 9-12 April 2016; 12 runs), not SPAN, not a delay, not the clock: one scalar per run explains 98 % of the run-to-run variance (rank-1 shape of 9 complex knots), from the run's own front taps. Held out, 259 captures: 10-12 kHz -40.1 / -38.8 -> -47.5 / -46.2 per run, -48.8 / -48.3 from the 15 nearest captures in time.
- **E, alignment:** refitting o / front shift / ppm on a high band costs the low band 4-9 dB and gains <= 2 dB: the bench's alignment is already optimal to 1e-4 samples at 5-10 kHz. The high-band delay (+0.004 samples) is common to the front taps and the branch, so it is a chain phase, not misalignment.
- **Bench:** `tools/nl_model.py` `NL_HF=1` (`hf_corr.json`, built by `f_table.py`: shape from the other runs, scalar from the R front taps of the 15 nearest-in-time captures, the capture's own front taps excluded; unity below 10 kHz and above 12.8 kHz; 36 shared numbers + one scalar per capture run; a measurement-side term of the capture chain like the clock and the gain, not the unit). Full 1,044-capture bench `bench_warp2_bk_kf_a0.7_lb2_hf` against `..._lb2` (the current default, low band folded in by the other session): **10-12 kHz -40.8 / -39.9 -> -48.1 / -48.4 dB (median)**; every window 0-2 s and every band 20 Hz-10 kHz unchanged (worst channel-window +0.49 dB, worst band +0.20 dB; medians 0.00 dB). **But 320 of 2,088 channel-bands at 10-12 kHz get worse by more than 0.5 dB (worst +13.9 dB)**: mostly captures that already nulled well (-47 -> -42 dB), clustered in SPAN 0.3 (115 / 332), SPAN 0.1 (86 / 350) and the trimmed captures (123 / 404; gain only 1.2 dB), plus the whole (0.2, -0.2, *) family. Strength 0.6 (`NL_HF_ALPHA`) on a worst-case-plus-sampled subset does not settle the trade-off. Left opt-in, default off.
- **Open:** why those clusters fail (time-neighbour scalar wrong in a run, little HF energy in trimmed late content, or a late-section term; D saw 37 of 40 failures with LATE > 0), and the per-run drift inside runs 7 and 11.
