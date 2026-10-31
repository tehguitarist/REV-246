# Tomkin library: layout, naming, harness results

The Tomkin Mastering EMT 246 IR library (README credit: Thomas Kienast), in `irs/EMT 246 captures/`. Everything about the files themselves lives here; what to do with them is in `NEXT_STEPS.md`.

*Status: October 2026. Harness run done; service manual and images read.*

## 1. What it is, and the licence

- 4,874 stereo 24-bit WAVs (6 GB), gitignored. **Licence (user): commercial use is allowed; the only restriction is no redistributing the files as they are.** Derived data (parameters, fitted models, reports without audio) may be committed.
- Also in the folder: `README.txt` (the naming rule), `emt246_front.jpg` (front panel in Non-Lin) and `emt246_operation_and_function.jpg` (the operation card). Both images are parsed into `EMT246_unit_facts.md` §3.

## 2. Layout and naming

```
irs/EMT 246 captures/
  48kHz/ and 96kHz/            identical file sets
    Short RT/ (512)  Med RT/ (519)  Long RT/ (360)  Non-Lin/ (1,046)
```

**Working set: 96 kHz only** (`irs/EMT 246 captures/96kHz/`). The 48 kHz files are the same capture as the 96 kHz ones (§5), so they add nothing and cannot serve as a held-out set.

Fields run left to right in front-panel order, and the last three numbers have an implied decimal point:

| Folder | Pattern | Example → setting |
|---|---|---|
| Short / Med / Long RT | `LP_DELAY_AMPL_LF_RT_HF_dc.wav` | `12_0_99_12_08_09` → LP 12 kHz, DELAY 0 ms, AMPL 99 %, LF×1.2, RT 0.8 s, HF×0.9 |
| Non-Lin | `LP_DELAY_SPAN_EARLY_LATE_dc.wav` | `12_0_05_-30_20` → LP 12, DELAY 0, SPAN 0.5 s, EARLY RT −3.0 s (negative = reversed), LATE RT 2.0 s |

- The Non-Lin order is confirmed by the front-panel photo (display order DELAY, SPAN, EARLY, LATE; EARLY −3.0 shown on the unit) (**S**).
- `dc` is probably "DC removed" (**H**).
- The harness parses all 4,874 names (`tools/harness/config.json` → `naming`: `kienast-sml`, `kienast-nonlin`). Program S/M/L comes from the folder name, so run it on the library root or a `<rate>/` folder, never from inside a `* RT` folder.

**Value sets (both rates):** LP is always 12 and DELAY always 0, so **there are no LP or DELAY steps** and no way to see the LP or pre-delay law here.

| | S | M | L |
|---|---|---|---|
| AMPL | 0 / 30 / 60 / 90 | 10 / 60 / 99 | 0 / 80 |
| LF× | 0.1 / 0.4 / 1.0 / 2.0 | 0.1 / 0.4 / 0.7 / 1.2 / 2.0 | 0.1 / 0.8 / 1.2 / 2.0 |
| HF× | 0.2 / 0.4 / 0.6 / 1.0 | 0.1 / 0.2 / 0.4 / 0.7 / 0.9 | 0.1 / 0.3 / 0.5 / 0.8 / 1.0 |
| RT (s) | 0.3–2.5 (8 values) | 0.6–5.0 (8) | 2.2–20 (9) |

Non-Lin: SPAN 0–0.5 s, EARLY RT −3.0…10 s, LATE RT 0–2.0 s. The grid is sparse (not the full 16-step tables).

**File properties** (96 kHz set, 2,437 files): variable length (about 0.7–11 s), trimmed and faded per setting. 1,622 files end in digital silence, 1,751 are flagged `fade_out`, 1,323 have under 2 ms pre-roll. So **RT must come from the label, not the file** (`label_mismatch_rt` on 1,190 files is mostly truncation), and absolute noise floors are not available.

**Fades and bad captures (pulse check, `model_s_results.md` §26, `ml_structure.md` §14):**
- **The fade differs per program (measured per file, `ml_structure.md` §15, `tools/ml_fade.py`).** All are linear amplitude ramps. S: 0.70 s (q2's `ramp.json`, compensated in every S score). **M: 1.0 s** (start 0.98–1.09 s before the end, IQR; some files under 1.5 s 0.3 s; many ramps cut short by the trim). **L: per file, about 0.31 × the file's length** (IQR 0.25–0.40; 0.3–12 s), hand-set near a level ~25 dB below the early reverb (H). `ml_lines.usable_end(y, f, prog)` reads the per-file starts from `analysis/ml_fade/fade.json` (+ 0.05 s guard); before that it assumed 0.72 s for every program, which biased late decays 4–12 % fast and half of L's 125 Hz RT "shortfall".
- **13 S files with a sample slip at 100 ms file time** (`model_s_results.md` §28): from 99.9 ms after the file's start the rest is delayed by ~3–6 samples at 96 kHz (1.2–1.6 device samples). All AMPL 30 / 60, at RT 0.3 / 0.6 with HF 0.2 / 0.4: 0.1_0.3_0.2 60; 0.1_0.3_0.4 60; 0.4_0.3_0.2 60; 0.4_0.3_0.4 30, 60; 1.0_0.3_0.2 30, 60; 1.0_0.3_0.4 30, 60; 1.0_0.6_0.2 30, 60; 2.0_0.3_0.2 30; 2.0_0.6_0.2 30. Leave them out of fits and repeat references (`tools/loop_read/stepband.py` finds them).
- **4 more S files with a departure at the same 99.9 ms file time but no lag step** (`model_s_results.md` §35.2, `tools/loop_read/glitch_scan.py`): `12_0_0_10_03_04` and `12_0_90_10_03_04` (1.0_0.3_0.4 AMPL 0 / 90), `12_0_30_20_03_04` (2.0_0.3_0.4 AMPL 30), `12_0_30_20_08_02` (2.0_0.8_0.2 AMPL 30). From 97 ms after tap 0 they differ from the same filter's other captures at −12…−16 dB (−46…−53 before), with no lag or gain change, so `stepband.py` misses them. Listed as `defects.GLITCH` and, since `model_s_results.md` §35.4, part of `DEFECTS` (= BAD | SLIPPED | GLITCH, 21 files; the clean references, bench and gate summaries use them). 1.0_0.3_0.4 therefore has no clean capture (AMPL 30 / 60 SLIPPED, 0 / 90 GLITCH) and is not scored; 2.0_0.3_0.4 and 2.0_0.8_0.2 keep AMPL 0 / 60 / 90.
- **5 M files with a one-sample slip at 126–139 ms file time** (`ml_structure.md` §23, `tools/ml_slip_scan.py`): `12_0_10_07_30_04` (0.7_3.0_0.4 AMPL 10), `12_0_60_12_18_04` (1.2_1.8_0.4 AMPL 60), `12_0_60_20_06_07`, `12_0_60_20_08_04`, `12_0_60_20_12_04` (2.0_0.6_0.7 / 2.0_0.8_0.4 / 2.0_1.2_0.4 AMPL 60). From the slip on, the rest is 1.00 file samples early (0.27 device samples), in both channels; M's R first pass nulls at −17…−21 dB instead of −48…−52, and the law stage's clock fit traded its clock (−9…+4 ppm against M's 12.6 ± 0.3) for it. None in L (a scan of all 839 M / L captures). Listed as `defects.ML_SLIPPED`, left out of the M / L pulse report.
- **Four bad S AMPL-0 captures:** `12_0_0_10_06_06` (1.0_0.6_0.6), `12_0_0_10_06_10` (1.0_0.6_1.0), `12_0_0_04_10_04` (0.4_1.0_0.4, from ~300 ms) and `12_0_0_01_03_02` (0.1_0.3_0.2) null against their own setting's other captures at only −10…−26 dB (normal −50…−57 dB). Drop them from per-setting fits.

## 3. Harness results (96 kHz set: `analysis/library/tomkin96/report.md`)

- **S clusters** found in 384 files (101 with exactly the S loop set `1321-1389-1461-1587`, so the library unit runs the same firmware as the CLA unit). 15 files have a `novel_cluster` (all Non-Lin). M and L files have no cluster match only because `cluster_templates` holds the S pattern alone (and the 16 ms search window is too short for L's 55 ms cluster); their clusters are S's ×2 and ×4 (`ml_structure.md` §1).
- **116 unknown loop sets**, the M/L candidates. The commonest: `1108-1254-1386` (310 files), `2641` (37), `231-628-803-1108-1292-2032` (37), `628` (35), `628-1108-1386` (28), `5279` (26), `1321-1389-1461` (26, a partial S set), `5279-10826-12176…` (15+15). **Assigned by folder** (`ml_structure.md` §3.1): `2641…` sets are **M** (Med RT), `5279…` sets are **L** (Long RT), and `1108-1254-1386`, `628…` and `231-628-…` are **Non-Lin**. The M/L entries are ACF views (sums of lags) of one five-term recursion per program, not loop sets like S's.
- **AMPL pairs are repeat captures of the reverb.** Settings that differ only in AMPL (S 128 groups of 4, M 160 of 3, L 180 pairs) were captured separately, and AMPL changes only the cluster. They are the library's only repeat captures; M's and L's tails differ between them (time variance), S's don't (`ml_structure.md` §4).
- **Kernel (B2):** the 96k library kernel (8 files, `analysis/library/tomkin/kernels/`; the 48k kernel there is the same capture rendered at 48k) sits within 0.64 dB of the CLA kernel over 100 Hz–11 kHz. The analogue chain around them is modelled in `tools/analog/` (`EMT246_unit_facts.md` §2b).
- **Chain kernel (October 2026, replaces the B2 kernel for the device-domain caches):** `analysis/library/tomkin96/kernels/chain2.npz`, built by `tools/chain_kernel.py` (`ml_structure.md` §12). The B2 kernel is 5.2 ms long, the mean of 8 R tap-0 waveforms, applied to both channels, and its truncation leaves a low-frequency ripple on every impulse (§10 there). `chain2` is per channel: a 19 ms FIR times the schematic's analytic sub-20 Hz high-pass. It is fitted jointly on the cluster-only windows of 884 S/M/L captures, with a clock and shift per capture (`analysis/chain_kernel/shifts.json`) and the M/L output kernels refitted alongside. R is pinned to the B2 kernel's gain and delay at 300–3000 Hz, so levels and timing keep their convention. Use it with `EMT_KERNEL=chain2` (caches `analysis/ml_structure/cache_chain2/`, `analysis/q1/cache_chain2/`, about 1.1 GB, rebuilt in a few minutes with 4 workers, 1.5 GB per worker on the longest file). The harness itself (`run.py`, `files.csv`) still uses the B2 kernel; its per-file measurements don't need the correction.
- **L/R offset (B4):** 19.19 µs median, 19.19–19.21 µs in every one of the 384 S files (about 92 timing steps of 208.3 ns); the CLA unit gives 18.6 µs (about 89 steps).
- **Capture clock:** median −13.2 ppm at 96 kHz (IQR 4.8); 15 outliers. (The 48k files read −15.2 ppm, a constant 2.0 ppm lower in every pair, a resampler artefact.)
- **Pair classes, within the 96k set:** `repeat` / `near` between *different* labels mean identical or near-identical early content (mostly Non-Lin settings that differ only in LATE RT or EARLY RT, plus a few HF-only S pairs); they are not repeat captures. The cross-rate classes from an earlier two-rate run were unreliable (zeroed tails make `noise_ref` meaningless, and the ±40-sample alignment fails on ambiguous onsets) and were retired; use `crossrate.csv`.
- Planned sets `lowpass` and `predelay` show as empty in the coverage tables (no LP or DELAY steps).

## 4. Harness usage

```bash
.venv/bin/python tools/harness/run.py "irs/EMT 246 captures/96kHz" --name tomkin96     # 96 kHz only (the working set); about 15 min cold
.venv/bin/python tools/harness/run.py "irs/EMT 246 captures" --name trial --sample 60  # quick spread-out subset
.venv/bin/python tools/harness/crossrate.py "irs/EMT 246 captures" --out analysis/library/tomkin/crossrate.csv
```

- Per-file results are cached per run name, so an interrupted run resumes. About 0.1 s per file per pass (10 cores), two passes; a cold 96 kHz-only run took about 15 minutes (the first, two-rate run took about 11 minutes warm).
- `--sample N` takes N settings spread evenly with every rate of each; `--limit` takes the first files alphabetically (one folder only).
- `--separate-rates` (config `separate_rate_captures`) exists for libraries whose rates really are separate captures. **Leave it off here.**
- Outputs: `analysis/library/tomkin96/` (report, library.json, coverage, plots; `files.csv`, `pairs.csv` and the cache are not committed, regenerate with the command above) and `analysis/library/tomkin/` (`crossrate.csv` and the library kernels only).

## 5. The 48 kHz and 96 kHz sets are the same capture

**Closed.** Tomkin was asked and said the library was recorded twice at the two rates. The data say otherwise: the two rate sets are renderings of one recording (**S**). We take the data as the truth, do not pursue the discrepancy further, and treat the 96 kHz set as the working set. The direction (which rate was derived from which) is unknown and does not matter for our purposes. Direct comparison in the capture domain (`tools/harness/crossrate.py`, 2,437 pairs):

- Every pair is sample-aligned (lag 0.0) with a constant gain of 0.9994 and agrees to **−67.5 dB median in the first 0.3 s**.
- The out-of-band tail (14.5–22 kHz after 0.5 s, about −150 dBFS) is **correlated between the rates** (median 0.79, minimum 0.37, above 0.9 in 974 pairs). Independent captures would have independent noise there.
- The harness's cluster-tap time for a setting differs between the rates by a constant 0.0466 µs with a scatter of **0.09 ns**, while the same quantity scatters by **2–27 ns between different settings** recorded separately: the rates share their estimation error.
- For 300 random pairs, the 48k file over the decimated 96k file has one fixed response in every band: |H| std 0.001–0.002 dB, phase std 0.01–0.03 mrad.
- The 96k file is a flat 24-bit floor from 14 kHz to 47 kHz, as is the 48k file: neither carries interface noise above the rounding floor.

**What this cannot show:** which file came from which, because the 96k file has no content above 14 kHz. Either way the 96k set contains all the information. The earlier hope of using 48k/96k pairs as repeat captures for time-variance tests is gone, and the question is closed.

## 6. The Altiverb 246 IRs (`irs/Altiverb EMT256 IRs/`, not committed)

Captured for this project by a third party running the Altiverb 246 plugin. No licence is known, so **the audio is not committed** (`.gitignore`). Derived data (`analysis/altiverb/summary.json`) may be. Findings: `EMT246_unit_facts.md` §5c; code: `tools/altiverb_compare.py`.

- **What the files are:** 66 stereo 48 kHz 24-bit WAVs, 26.00 s each, named `EMT256_<preset>.wav`. They are **sweep recordings of the plugin, not IRs**: the stimulus is `tools/make_sweep.py` (1 s lead-in, 10 s log sweep 20 Hz–20 kHz, 15 s tail, played at −18 dBFS). Deconvolve with that script's inverse filter; lag 0 is at LEAD + sweep length (zero plugin latency assumed).
- **Presets:** families `NON-LIN SHORT RT 0.2…2.2 s` (S, six files), `MED RT 0.5…5.0 s` (M, ten), `LONG RT 2…10 s` (L, seven), `REV-250 0.4…4.5 s` (nine, probably the EMT 250 and not used), `NON-LIN`, `PRESET1…10` (PIANO … CATHEDRL) and `PRESET72…90`, `PRESET23`, `PRESET27`, `PRESET54 GATED`. The six `NON-LIN SHORT RT` files are the same IRs the first batch had under `SHORT RT`; the naming is the plugin's and has nothing to do with Tomkin's Non-Lin.
- **Settings are only the preset name.** LF, HF and AMPL are not given; they are not neutral. Inferred from the decay and the cluster level: HF 0.9 / 0.7 / 0.4, LF 1.0 / about 1.1 / about 1.2, AMPL about 80 / 80 / 60 % for S / M / L, the same at every RT (unit facts §5c; the pooled LF fit is unreliable, use the RT-matched slice).
- **It is not 96 kHz and not Tomkin's capture,** so the "96 kHz only" rule does not apply to it, and it is not a repeat capture of Tomkin: separate capture, relative clock about −45 ppm, polarity inverted.
- **Capture lessons (the first batch was spoilt):** at −6 dBFS with a 6 s tail, the dense presets clipped by 0…+3.4 dBFS (many thousands of samples in CATHEDRL and CONCERT, pre-IR distortion 6–15 dB above the method floor) and tails past 6 s were lost. At −18 dBFS with a 15 s tail, no file clips and the longest tails (LONG RT 10, PRESET90) are down to −139 dB by 13–14 s. Where neither run clipped, the two levels agree to −104…−112 dB: the method and the plugin repeat to about −100 dB.
- **Use:** structure and decay checks only (S presets: `pulse_check --stage units`; MED / LONG against the M / L model: `tools/altiverb_ml.py`, `ml_structure.md` §17). One frozen realisation per preset, so no M/L time-variance test and no held-out nulls of the M/L reverb. Not run through the harness (its naming parsers do not know these files).
- **The `REV-250` presets** are the 246's Program E (EMT 250 program); the first measurement pass on them (decay law, onset/stereo, three loops at ~168/190/211 ms, chain) is in `EMT250_program.md`.

## 7. Front-panel settings of every Altiverb preset (factory presets, 31 October 2026)

Best-guess panel settings for all 66 Altiverb files, as the plugin's factory presets. Code: `tools/altiverb_presets.py` (stages `features`, `ngrid`, `fit`; intermediate `cache/features.json`, `cache/n_grid.npz`, not committed); results `analysis/altiverb/presets.json` (value, step index, C / S / H label and a one-line reason per control and preset; the families' leave-one-out check and the decay calibration); generator `tools/plugin_codegen/presets.py` → `plugin/src/presets/FactoryPresets.h` (64 presets).

**Program** from the first 70 ms, in device samples: S / M / L by the cluster taps (§6, unit facts §5c), E by R's direct taps 0 / 429 / 718 with L 314 later, N by the front taps 41 / 410 / 463 (R 278); one preset without a cluster (PRESET90) by the reverb's own L–R onset lag (1320 = L's). Two files are none of the plugin's programs and are left out: the unlabelled `REV-250` (a bare two-tap delay, `EMT250_program.md` §14.2) and `PRESET7 ECHOROOM` (a recirculating echo, L every 99.8 ms × 0.49, R every 65.8 ms × 0.62, flat T60 ≈ 1.0 s): the 246's delay / echo program.

| Family | Files | Program |
|---|---|---|
| `NON-LIN SHORT RT` + PRESET72–76 | 6 + 5 | S |
| `MED RT` + PRESET1, 3, 5, 27, 78, 79, 81, 82 | 10 + 8 | M |
| `LONG RT` + PRESET4, 9, 10, 84, 85, 87–90 | 7 + 9 | L |
| `REV-250 <RT>` + PRESET2, 8, 77, 80, 83, 86 | 9 + 6 | E |
| `NON-LIN` + PRESET6, 23, 54 | 1 + 3 | N |
| `REV-250`, PRESET7 | 2 | delay / echo (left out) |

**Method** (each value relative to its family, whose RT is in the name):
- **DELAY**: onset (tap 0; E's R direct tap; N's front-frame zero) minus the family's, plus the family's DELAY. S family **10 ms (S)**: the CLA IR, DELAY 10 on its display, has its first tap at 10.07 ms, the S family at 10.08, and CLA's other settings are the S family's (LF 1.0, HF 0.9, AMPL 80). The M / L / E / N families sit at ×2 / ×4 / ×2 / ×2 of S (20.33 / 40.43 / 20.33 / 20.5 ms): since the CLA unit shows S has no built-in onset, M and E read as DELAY 20 (S). L and N are set to **39 / 19 (H)**: that matches the Altiverb's onset in the plugin's timing (the plugin's DELAY-0 onset is 1.13 ms (L) / 1.57 ms (N) after S's), against the round 40 / 20 (decision 31 October); their presets follow at −1 ms from the round reading. PRESET6 sits 2 ms before NON-LIN, so N's family DELAY is at least 2 ms. Every numbered preset lands within 0.2 ms of a whole ms from there (S / M / E: 12 / 14 / 16 / 24 / 26 / 30 / 34 / 36 / 40 / 44 / 46 / 50 / 64 / 80 ms; L / N one less than the round 70 / 80 / 100 / 115 / 150 / 180 and 18 / 20 / 40).
- **RT / LF / HF**: the 246's loop law (`ml_lines.line_coefs` / `rt_of` at the program's line lengths; E's three loops) times a per-band calibration from the family (median of measured / law octave T60 at the family's RT label and LF / HF, S 1.0 / 0.9, M 1.1 / 0.7, L 1.2 / 0.4, E 1.49 / 0.69), fitted in log T60 at 125 Hz–8 kHz with the family's scatter as weights (3 numbers per preset). **Leave-one-out on the 32 family presets:** RT within 1–3 % of the label from RT 0.6 (S), 1.0 (M, E), 3.0 (L) (short ends off: S 0.2 → 0.31, M 0.5 → 0.66, E 0.4 → 0.47); HF within ±0.02 (except M 1.0 0.77, L 2.0 0.56, E 0.4 0.88); LF scatters ±0.15 (S 0.84–1.17, M 0.99–1.13, L 1.17–1.32, E 1.32–1.56). So RT and HF are S where the fit is good (χ²/dof < 3), LF is H throughout. PRESET10, 89 and 90 hit LF 2.0 (the range end) with χ²/dof 3–40: everything H there.
- **AMPL**: S / M / L the cluster-tap / reverb ratio over Tomkin's linear ratio-per-AMPL law at the fitted RT / LF / HF (the `--stage settings` method: families 80 / 80 / 60); L's scale rests on Tomkin's one L AMPL, so L is H. E: R's direct tap over the 40–240 ms reverb RMS against the family's at the same RT, the family taken as 80 (H, as `e_plugin_spec`). PRESET90 has no cluster: AMPL 0.
- **LOW PASS**: tap 0's response over the family's tap 0 (1–12 kHz). S and N presets and half the M presets are flat (open, C). The others show a gentle roll-off that is **exactly a 2-tap device-rate FIR [1 − a, a]** (rms 0.01–0.04 dB), the same form as M / L's output kernel (§14.1 of `EMT250_program.md`): a = 0.06 / 0.09 (M; family 0.047), 0.41 / 0.49 (L; family 0.23), 0.33 / 0.41 / 0.49 (E; family 0.23), in a few discrete values that follow the preset (LIVE brighter than REF in the M pairs and L's MED pair), not the DELAY. What sets a is unknown (H: the LP lever, or an output-stage setting the panel does not show). The plugin's LP is a 2nd-order Butterworth, so the step with the least amplitude error is taken: 10 kHz for a ≥ 0.33, open below; H. PRESET23 is different: a steep low-pass, −1.3 / −5.3 / −13 / −39 dB at 8 / 9 / 10 / 11 kHz → LP 8 (H for the step; the cut-off is clear).
- **N** (SPAN / EARLY / LATE): `nl_model.render` (device side, main taps at +36 as on the Altiverb, NonLin §11.7) over the library steps (SPAN 0–0.5 by 0.1, the 16 EARLY steps, LATE 0–2.0 by 0.2: 1,056 renders, `n_grid.npz`), against 5 ms band envelopes (100 Hz–10 kHz, frames within 50 dB of the peak, one gain per channel). It recovers §11.7's NON-LIN (0.3 / 10 / 0.4) and GATED (0.2 / −3 / 0.8) exactly; S where no runner-up is within 15 % of the cost.
- Snapped to the plugin's steps (`Parameters.cpp`: RT tables per program, LF 16 steps, HF 0.1–1.0, AMPL 10 %, DELAY 1 ms, LP table; N: EARLY in the LF slot, LATE in the HF slot, SPAN in the AMPL slot, RT at its default step).

**Families** (RT from the name, C): S LF 1.0 (S) / HF 0.9 (S) / AMPL 80 (S) / DELAY 10 (S); M LF 1.1 (H) / HF 0.7 (S) / AMPL 80 (S) / DELAY 20 (S); L LF 1.2 (H) / HF 0.4 (S) / AMPL 60 (H) / DELAY 39 (H); E LF 1.5 (S) / HF 0.7 (S) / AMPL 80 (H) / DELAY 20 (S); NON-LIN SPAN 0.3 / EARLY 10 / LATE 0.4 (S) / DELAY 19 (H). LOW PASS open in all (S).

**Numbered presets** (label after each value; LP 12 = open):

| Preset | Prog | RT s, LF ×, HF ×, AMPL % (N: SPAN s, EARLY s, LATE s) | DELAY ms | LP kHz |
|---|---|---|---|---|
| 1 PIANO | M | RT 1.5 S, LF 1.2 H, HF 0.7 S, AMPL 70 S | 40 S | 12 C |
| 2 BRASS | E | RT 2.0 S, LF 1.0 H, HF 0.9 S, AMPL 100 H | 30 S | 10 H |
| 3 WOODWIND | M | RT 2.2 S, LF 1.1 H, HF 0.5 S, AMPL 30 S | 50 S | 12 H |
| 4 STRINGS | L | RT 3.6 S, LF 1.5 H, HF 0.4 S, AMPL 90 H | 69 H | 10 H |
| 5 VOCALS | M | RT 2.6 S, LF 1.3 H, HF 0.9 S, AMPL 50 S | 50 S | 12 C |
| 6 PERC | N | SPAN 0.3 S, EARLY 3 S, LATE 1.0 S | 17 H | 12 C |
| 8 CHAMBER | E | RT 1.8 S, LF 1.2 H, HF 0.6 S, AMPL 90 H | 26 S | 10 H |
| 9 CONCERT | L | RT 4.5 S, LF 1.8 H, HF 0.4 S, AMPL 50 H | 99 H | 10 H |
| 10 CATHEDRL | L | RT 5.2 H, LF 2.0 H, HF 0.2 H, AMPL 20 H | 179 H | 10 H |
| 23 SM SECT1 | N | SPAN 0.2 S, EARLY 10 S, LATE 1.4 S | 19 H | 8 H |
| 27 LG SECT1 | M | RT 2.6 S, LF 1.5 H, HF 0.5 S, AMPL 40 S | 64 S | 12 H |
| 54 GATED | N | SPAN 0.2 S, EARLY −3 S, LATE 0.8 S | 39 H | 12 C |
| 72 SM LIVE | S | RT 0.8 S, LF 0.7 H, HF 1.0 S, AMPL 100 H | 12 S | 12 C |
| 73 SM REF | S | RT 0.9 S, LF 0.8 H, HF 0.8 S, AMPL 90 S | 14 S | 12 C |
| 74 SM WARM | S | RT 1.0 S, LF 1.0 H, HF 0.5 S, AMPL 100 S | 16 S | 12 C |
| 75 MED LIVE | S | RT 1.4 S, LF 0.8 H, HF 0.9 S, AMPL 100 S | 24 S | 12 C |
| 76 MED REF | S | RT 1.6 S, LF 1.0 H, HF 0.6 S, AMPL 100 S | 26 S | 12 C |
| 77 MED WARM | E | RT 1.5 S, LF 1.2 H, HF 0.4 S, AMPL 90 H | 26 S | 10 H |
| 78 LG LIVE | M | RT 2.2 S, LF 1.0 H, HF 0.8 S, AMPL 90 S | 34 S | 12 C |
| 79 LG REF | M | RT 1.8 S, LF 1.2 H, HF 0.5 S, AMPL 70 S | 36 S | 12 H |
| 80 LG WARM | E | RT 2.2 S, LF 1.5 H, HF 0.3 S, AMPL 70 H | 36 S | 10 H |
| 81 SM LIVE | M | RT 2.6 S, LF 1.2 H, HF 0.7 S, AMPL 80 S | 44 S | 12 C |
| 82 SM REF | M | RT 2.6 S, LF 1.5 H, HF 0.4 S, AMPL 40 S | 46 S | 12 H |
| 83 SM WARM | E | RT 2.8 S, LF 1.8 H, HF 0.3 S, AMPL 30 H | 46 S | 10 H |
| 84 MED LIVE | L | RT 3.1 S, LF 1.3 H, HF 0.6 S, AMPL 60 H | 79 H | 10 H |
| 85 MED REF | L | RT 3.1 S, LF 1.8 H, HF 0.3 S, AMPL 30 H | 79 H | 10 H |
| 86 MED WARM | E | RT 3.3 S, LF 1.8 H, HF 0.2 S, AMPL 30 H | 80 S | 10 H |
| 87 LG LIVE | L | RT 4.5 S, LF 1.5 H, HF 0.5 S, AMPL 50 H | 114 H | 10 H |
| 88 LG REF | L | RT 4.5 S, LF 1.8 H, HF 0.3 S, AMPL 20 H | 114 H | 10 H |
| 89 LG WARM | L | RT 4.5 H, LF 2.0 H, HF 0.2 H, AMPL 20 H | 114 H | 10 H |
| 90 VERY LG | L | RT 6.0 H, LF 2.0 H, HF 0.1 H, AMPL 0 S | 149 H | 12 H |

**Caveats.** The presets are LIVE / REF / WARM variations: LIVE = higher HF and AMPL, REF = lower HF, WARM = lowest HF. AMPL 100 in four S presets is the range end (the fit gives 97–110 % there; above 105 % marked H). L and N DELAYs follow the plugin's timing (39 / 19 family), not the round panel reading (40 / 20); if the plugin's L / N onset convention changes, add 1 ms back.

**7.x LOW PASS from the presets (31 October 2026, `tools/lowpass_fit.py`, `analysis/lowpass/report.txt`).** Family-relative first taps: the darkened M / L / E presets are not a filter added after the program but a different value of the program's own 2-tap output-kernel coefficient [1 - a, a] (0.002-0.009 dB, tap spacing = the device period, C): M family 0.047, numbered 0.059-0.062 or 0.090-0.092 (no numbered M preset is flat: correction to the table above); L / E family 0.23, numbered 0.333, 0.411 or 0.498 (L and E share values to +-0.001). It doesn't vary with RT / HF / AMPL / DELAY, and a 2-tap can't put -3 dB below 6.6 kHz, so it isn't a 2-12 kHz lever (S): most likely a per-preset kernel variant (H). The "LP 10" guesses are three different kernel values. PRESET23 (N) is a linear-phase cut (-3 dB 8.5 kHz, -39 dB 11 kHz, notches 11.09 / 12.84 kHz, stop band -27...-30 dB): a symmetric 13-tap device-rate FIR fits it to 0.019 dB, no causal IIR does; either the lever at 8 as a delay-compensated FIR or the Altiverb's own processing (H). The LP lever's shape stays H.
