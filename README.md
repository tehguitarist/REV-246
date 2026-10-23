# REV-246: EMT 246 reverb reverse-engineering workspace

Goal: recreate the EMT 246 digital reverb's S/M/L programs (Non-Lin as a bonus) as an efficient, ideally bit-accurate model, working **only from impulse responses and documents**. There is no hardware access.

**Start here:**
1. `docs/NEXT_STEPS.md`: status, working rules, the open questions Q1–Q12 and a suggested order of work.
2. `docs/library_notes.md`: the Tomkin IR library (layout, naming, harness results, the 48k = 96k finding).
3. `docs/EMT246_unit_facts.md`: everything known about how the unit works, with confidence and sources.

## File map

### Docs (`docs/`)

| File | What it is |
|---|---|
| `NEXT_STEPS.md` | Hand-over: status, working rules (96 kHz only), open questions Q1–Q12, analogue stage, order of work, Non-Lin, gaps, modulation, Tomkin questions, UAD 250 test bench, plugin notes |
| `library_notes.md` | The Tomkin library: licence, layout and naming, value sets, harness results, cross-rate (48k = 96k) analysis, harness usage, the Altiverb 246 captures (§6) |
| `EMT246_unit_facts.md` | Consolidated fact sheet: hardware, converters, controls, programs, measured program-S algorithm, budgets, open questions, source index |
| `model_s_results.md` | Program S grey-box model and its null tests. §14–16 the library work (laws, AP230 steps, the 5-tap injection), §17 the consolidated model, §19–21 the deconvolution-kernel fix and refit, §22–22.9 the loop corrections at strong settings, §23–25 the low band. **§25.7 is the current S baseline** (§22.9 plus LF 0.4's shelf DC 0.414, the DC floor per (LF, RT) cell and the HF pole per (HF, RT) cell; all 128 settings held out by setting) |
| `ml_structure.md` | Programs M and L from the library: clusters and onsets (S's ×2 / ×4), the five-term late recursion (§7: S's four loops scaled ×2 / ×4 and coupled; the S carry-over test), lag precision, the time variance found with AMPL pairs as repeat captures, the noise floors and repeat-capture ceiling (S about −57 dB) that bound the null target, the line laws (§8), the time variance: a 0.83 Hz rotation of the coupling (§9), the first-pass gap (§10, the library's deconvolution kernel), the out-of-band check (§11) and the per-channel chain kernel `chain2` (§12, with the first-pass refit §12.7) |
| `EMT246_discovery.md` | The original discovery and viability study. Hardware photos analysis in §2.3a; period sources in §10a; history in §10 |
| `algorithm_research.md` | Reverb-algorithm lineage (EMT 250 patent, Blesser), classic and modern methods, topology hypotheses and tests (T1–T10), the stochastic-control analysis (§7.2b) |
| `period_dsp_practice.md` | 1978–1990 reverb and TMS32010 practice (EMT, Lexicon, patents, German sources) turned into ranked, testable candidates (K1–K7) for the dense per-loop injection/read-out found in Q2, with test recipes, and the offline results for every candidate (§7: K1 holds, five equal signed taps per loop; K2–K7 not supported; in-sample) |
| `analysis_harness_spec.md` | IR-library harness test specification (stages A–G), plus implementation status of A–C |

### Code (`tools/`, `tests/`)

| Path | What it does |
|---|---|
| `tools/analyze_ir.py` | Single-IR deep dive → `analysis/<name>/{summary.json, report.png, h_dev.npy}` |
| `tools/deconv_grid.py` | Recovers the device-rate IR (19.2 MHz/728 ≈ 26,373.6 Hz) from a capture |
| `tools/topology_tests.py` | Topology tests T1/T2 against synthetic controls |
| `tools/model_s.py` | Program S model and studies. Stages: default, `--tighten`, `--tail`, `--diffuser`, `--allpass`, `--groups`, `--stochastic` |
| `tools/ml_structure.py` | M/L structure from the 96 kHz library, with S as control: clusters and onsets from AMPL pairs, held-out tail-predictor (characteristic polynomial) tests, lag precision, AMPL-pair time-variance test, noise floors vs CLA, S repeat-capture null with a per-capture clock fit → `analysis/ml_structure/` |
| `tools/ml_lines.py` | The M/L line gains and filters (`ml_structure.md` §8): the §7 model closed through the coupled lines with S's laws at the M/L line lengths; held-out recirculation nulls against repeat captures, limit scans, decay per octave, and the law fitted per setting. One checkpoint file per AMPL group (resumable; `--chunk i/n`) → `analysis/ml_lines/`. `fitted(prog, key)` rebuilds a setting's captures and models from those files (for residual work) |
| `tools/ml_variance.py` | The M/L time variance (`ml_structure.md` §9): per capture, the residual after `ml_lines` fitted as kernels on candidate derivatives, held out across channels (the coupling rotation and its controls), θ's rate from the passes, ghost lags per octave; `sim` renders the rotating model through a simulated exponential sweep against Tomkin's AMPL pairs. One checkpoint file per AMPL group (resumable; `--chunk i/n`) → `analysis/ml_variance/` |
| `tools/ml_carryover.py` | The S carry-over test for M and L (`NEXT_STEPS` §5 item 5): diffuser search by matched-filter sparsity, the 20 read taps per channel in every AMPL group, each tap's recirculation through the four lines (the coupling matrix), and held-out first-pass nulls against repeat captures (with the AP damping step and the M/L output kernel) → `analysis/ml_carryover/` |
| `tools/chain_kernel.py` | The per-channel chain kernel (`ml_structure.md` §12): cluster-only windows of S/M/L on a 4× device grid, per-capture clock and shift, FIR × the schematic high-pass, held-out variants → `analysis/chain_kernel/` and `analysis/library/tomkin96/kernels/chain2.npz`; `--stage tail` / `write --name chain3` the long-tail kernel chain3 (`ml_structure.md` §13, the current one). With `EMT_KERNEL=chain3` (or `chain2`) the device-domain tools (`ml_structure`, `q1`, `q2`, `ml_*`, `firstpass_gap`, `outofband`) read caches rebuilt with it and write `_chain2` results |
| `tools/q1_generalise.py` | Q1: Path B across the library (loop fits, control laws, AP230 step, shared-injection joint fit, held-out splits, repeat-capture ceiling, CLA check) → `analysis/q1/` |
| `tools/q2_injection.py` | Q2: fade-free windows, Tomkin fade compensation (fixed 0.700 s ramp), per-setting free injections, q1's joint bench with corrected windows, fade-compensated nulls, residue regression, the `branches` bench (shared injection with branch exponents and window length), the `filters` / `loops` stages (loops refitted per setting with an injection frozen), the `tied` / `diffuser` / `ap230` / `capped` stages (the compact 5-tap model, −1/√2 diffuser, AP230 steps, saturating laws), `cla` (cross-unit), and `final` (the consolidated S model, `model_s_results.md` §17; with `--compact` / `--ap230` the chain2 refit of §21 and with `--gamma` / `--lfrt` / `--hfdp` / `--pramp` the loop terms of §22). `free --refs` is checkpointed per setting. With `EMT_KERNEL=<name>` every output gets `_<name>` and the repeat-capture references come from `nulls_free_<name>.json` → `analysis/q2/` |
| `tools/firstpass_gap.py` | The first-pass gap (`ml_structure.md` §10, §12.3): per AMPL group the model parts and fits (`collect`), transfer error and residual spectra (`report`), the deconvolution-kernel correction (`chain`) → `analysis/firstpass_gap/` |
| `tools/outofband.py` | Out-of-band check (`ml_structure.md` §11, §12.5): nulls and decay at 25–100 Hz and 10–12 kHz, the sweep simulation, kernel response → `analysis/outofband/` |
| `tools/ap230_form.py` | AP230's damping placement at b > 0 (`model_s_results.md` §21.2): the two-tap filter in the feed-forward / feedback paths, both free → `analysis/ap230_form/` |
| `tools/loop_fir.py` | S's loops at strong settings (`model_s_results.md` §22–22.8), injection frozen: free per-loop corrections held out by capture (`fit`, `report`, `coefs`), RT interpolation held out by setting (`interp`), grids on cached setups for the DC floor (`floor`), the shelf (`shelf`) and the HF cap / pole offset (`hfcap`), the parts of a correction (`probe`), tilt shapes (`tilt`). One checkpoint per setting in `analysis/loop_fir/<kernel>/` (not committed) |
| `tools/loop_cv.py` | All loop terms with every constant chosen leaving one RT step out (8 folds), each setting scored with its fold's constants: the §22.9 baseline → `analysis/loop_cv/`; with `--lf04dc` §25.4 → `analysis/loop_cv/lf04dc/`; with `--cells` §25.6 → `analysis/loop_cv/cells/`; with `--hfcells` §25.7 (chain2) → `analysis/loop_cv/hfcells/`; with `EMT_KERNEL=chain3 --retake` the **current S baseline** (§25.11: §25.7's constants, compact numbers refitted on chain3) → `analysis/loop_cv/hfcells_chain3/`, and `--xadd <dense_mid tag>` adds a dense_mid correction held out by fold → `hfcells_chain3_<tag>/` |
| `tools/lowband_shelf.py` | S's LF shelf below 400 Hz (`model_s_results.md` §25): a free shelf per setting held out by capture, LF 0.4's DC constant on a grid and chosen leaving one RT step out, the DC floor and LF 0.1's DC ratio on a joint grid (§25.5), per-loop DC / Nyquist gains and the per-cell floor (§25.6), the HF pole per cell (§25.7) → `analysis/lowband/chain3_shelf_*` (reports `chain3_shelf_report.log`, `chain3_shelf_floor_report.log`) |
| `tools/dense_mid.py` | Diagnostic (`model_s_results.md` §25.9–25.10): a correction to the loop excitations fitted on the residual of training settings and scored held out by RT step. It can be dense, shared by the loops, per-tap, a tap level, a read kernel or the diffuser's gains and delays. Options: `--band` (fit band), `--project` (band-limit before scoring), `--test-all` (score every setting), `--stage decode` (read the common kernel) → `analysis/dense_mid/<kernel>/` |
| `tools/memguard.py` | Runs a command under a hard memory cap: polls the physical footprint of its whole process group (workers included), kills it above `--limit-gb` and prints the peak. Use it for every run (NEXT_STEPS §2) |
| `tools/lowband.py` | The S baseline at 25–200 Hz (`model_s_results.md` §23): per-band nulls and repeat references, the capture / model transfer, a fixed-transfer variant held out by RT fold, the cross-program check on S / M / L, third-octave decay errors → `analysis/lowband/` |
| `tools/q2_candidates.py` | Offline tests of the period-practice candidates on the frozen flat injection: in-band sparsity with controls (K1), the comb factor and tied tap model, a recursive-stage scan (K2/K5), tap-position structure (K3), a common distributor (K6), the first-pass FIR, the residual (K5), and per-tap checks across the 128 per-setting injections (K7, K4) → `analysis/q2/` |
| `tools/make_sweep.py` | Stimulus for sweep-capturing a convolution plugin: log sweep (48 / 96 kHz, −6/−12/−18 dBFS, optional 15 s tail), its inverse filter, a single-sample impulse → `sweep_request/` (not committed) |
| `tools/altiverb_compare.py` | Altiverb 246 IRs against Tomkin: sweep deconvolution and clipping check, cluster taps, octave-band T60, tap spectrum vs the CLA kernel, S reverb null (`library_notes.md` §6, unit facts §5c) → `analysis/altiverb/` |
| `tools/altiverb_predict.py` | The §17 S model against the Altiverb S presets, held out (clock, shift, gain + a per-channel FIR on the cluster region), the Altiverb vs CLA comparison, polarity / channel-offset diagnostics (`model_s_results.md` §18) → `analysis/altiverb/predict.json` |
| `tools/analog/` | Analogue chain: `netlists.py` (input 7 246 131 and output 7 246 130 filters transcribed from the service manual), `mna.py` (AC solver, ideal or finite-GBW op-amps), `compare_kernel.py` (vs the reference kernel → `analysis/analog/`) |
| `tools/harness/` | IR-library harness, stages A–C: `run.py` (CLI), `crossrate.py` (direct cross-rate comparison), `core.py` (per file), `library.py` (library level), `naming.py`, `report.py`, `synth.py` (synthetic test library), `config.json` (constants, templates, naming patterns) |
| `tests/test_harness.py` | End-to-end harness tests: `.venv/bin/python -m pytest tests -q` |

### Data

| Path | Contents |
|---|---|
| `irs/EMT 246 captures/` (gitignored; Tomkin Mastering, no redistribution of the files) | The Tomkin/Kienast library: 96 kHz (the working set) and 48 kHz (the same capture), `Short RT`, `Med RT`, `Long RT`, `Non-Lin`, plus its README, the operation card scan and a front-panel photo. Naming scheme in `NEXT_STEPS.md` §2.1 |
| `irs/eugtone_EMT246_CLA_shortRT.wav` | The one hardware IR so far: CLA's setting, captured by eugtone (S, LP 12, DELAY 10, AMPL 80, LF 1.0, RT 1.0, HF 0.9). Third-party |
| `ref/kernel_eugtone_96k.npz` | Reference analogue-chain kernel (one cluster tap) |
| `analysis/shortRT/` | `analyze_ir.py` output for the CLA IR. `h_dev.npy` is regenerated, not committed |
| `analysis/model_s/`, `analysis/model_s_*.log` | Model study outputs (JSON, params, plots, logs) |
| `analysis/library/synthetic/` | Harness run on the synthetic test library |
| `analysis/library/tomkin96/` | Harness run on the Tomkin 96 kHz set: report, library.json, coverage, plots. Large per-file CSVs and the cache are not committed |
| `analysis/library/tomkin/` | Cross-rate result (`crossrate.csv`) and the library kernels |
| `analysis/chain_kernel/` | `chain_kernel.py` output: `fit.{json,log}` (held-out kernel variants), `chain2.json` (the written kernel's response, L/R, refitted M/L output kernels), `shifts.json` (per-capture tap-0 shift and clock); the windows `data_*.npz` are not committed |
| `analysis/ml_lines/` | `ml_lines.py` output: `law/` (one `.json` per AMPL group: nulls per window and variant, clocks, gains, band decay), `report.log`, `fit.json` / `fit.log` |
| `analysis/ml_variance/` | `ml_variance.py` output: `rot/` (one `.json` per AMPL group: held-out nulls per window and candidate, θ rms, kernels, per-octave ghost lags, per-event correlations), `report.log`, `sim/` + `sim.log` (rotating model through a simulated sweep; `sim_075*` the first run at 0.75 Hz / 3.9 s) |
| `analysis/ml_carryover/` | `ml_carryover.py` output: `diffuser`, `taps` (per-group tap positions and values), `routing` (recirculation gains, the coupling matrix) and `firstpass` (held-out nulls, fitted damping steps, output kernels), each `.json` + `.log` |
| `analysis/ml_structure/` | `ml_structure.py` output: `summary.json`, `run.log` (the device-rate IR cache is not committed) |
| `analysis/q2/` | `q2_injection.py` output: corrected windows (`fade2.json`), raw-end map for the fade ramp (`ramp.json`), per-setting free-injection results (`free.json`; the injections in `v/` are not committed), joint splits (`split_*.json`, `v_*.npy`), fade-compensated nulls (`nulls_*.json`), branch runs (`branches_*.json`; the frozen flat injection `branches_flat_0_ext1700.npy`), loop refits (`filters.json`), residue regression (`regress.json`; the residue arrays are regenerated), the candidate tests of `q2_candidates.py` (`k1_sparsity.*`, `k2_stages.*`, `structure.*`, `k5_resid.*`, `settings.*`), logs |
| `analysis/q2/*_chain2*` | The chain2 work (`model_s_results.md` §20–22): AP230 steps and the compact refit (`ap230_g0.7071_chain2`, `diffuser_tied_fix-0.707107_ap230new_chain2`), the free models and repeat-capture references (`free_chain2.json`, `nulls_free_chain2.json`; `v_chain2/` and the `*_groups/` checkpoints are not committed), the final-stage scores (`final_s_chain2_refit*`) |
| `analysis/firstpass_gap/`, `analysis/outofband/` | `firstpass_gap.py` and `outofband.py` outputs (reports, chain corrections; the per-group caches `groups_*` are not committed) |
| `analysis/ap230_form/` | `ap230_form.py` output (`forms_chain2.{json,log}`) |
| `analysis/loop_fir/` | `loop_fir.py` output: logs and summaries (`fit_`, `report_`, `coefs_`, `interp_`, `floor_`, `hfcap_chain2.log`; `chain2_coefs.json`, `chain2_interp.json`, `chain2_floor*.json`, `chain2_hfcap.json`, `chain2_probe*.json`, `chain2_tilt.json`); the per-setting fits in `chain2/` are not committed |
| `analysis/lf_rt/` | Held-out bench runs of the single loop terms (`model_s_results.md` §22.2–22.8): trained constants (`lfrt_train*.json`, `hfdp_train.json`, `split2_*.json`) and their scored settings |
| `analysis/loop_cv/` | `loop_cv.py` output: per fold the constants (`fold_<rt>.json`, `fold_all.json` = the baseline's constants) and scored rows, `constants.log`, `report.log` |
| `analysis/lowband/` | `lowband.py` output: per setting `chain2/`, `chain2_tfix/`, `chain2_transfer/`, `chain2_slopes/`; per group `chain2_hk_mltransfer/`, `chain2_hk_mlcheck/`; the fixed transfer per fold `chain2_tfix_{A,B}.json`; `chain2_{base,tfix}_report.log` |
| `analysis/analog/` | Analogue chain vs the CLA and Tomkin kernels (`summary*.json`, `chain*.png`) |
| `sources/` (local only, gitignored) | Third-party primary sources: data sheet, 252 brochure, catalogues, 245 service manual, **246S service manual (scans)**, EMT Courier, EMT 250 patent, Studio Magazin pages, magazine extracts, studiorepair internal photos. See `sources/README.md` |

## Setup and common commands

```bash
python3 -m venv .venv && .venv/bin/pip install numpy scipy soundfile matplotlib pandas pytest
.venv/bin/python tools/analyze_ir.py irs/<file>.wav --out analysis/<name>      # one IR, deep dive
.venv/bin/python tools/harness/run.py <library_dir> --templates              # inspect naming
.venv/bin/python tools/harness/run.py <library_dir> --name <lib>             # whole library, stages A–C
.venv/bin/python tools/harness/run.py "irs/EMT 246 captures/96kHz" --name tomkin96      # the Tomkin working set (96 kHz only)
.venv/bin/python tools/harness/run.py "irs/EMT 246 captures" --name trial --sample 60   # quick spread-out subset (both rates)
.venv/bin/python tools/harness/crossrate.py "irs/EMT 246 captures" --out analysis/library/tomkin/crossrate.csv   # direct 48k vs 96k comparison
.venv/bin/python tools/ml_structure.py [--stage clusters|loops|lags|pairs|floors|clock]                # programs M/L (and S as control), 96 kHz set
.venv/bin/python tools/ml_carryover.py [--stage diffuser|taps|routing|firstpass] [--only N]           # S carry-over to M/L (needs ml_structure's cache)
.venv/bin/python tools/ml_lines.py [--stage law|report|fit|all] [--only N] [--chunk i/n] [--workers 4]   # M/L line laws (needs ml_carryover's firstpass.json)
.venv/bin/python tools/ml_variance.py [--stage rot|report|sim|all] [--only N] [--chunk i/n] [--workers 4]  # M/L time variance (needs ml_lines' law files)
.venv/bin/python tools/make_sweep.py                                                         # sweep, inverse and impulse stimuli → sweep_request/
.venv/bin/python tools/altiverb_compare.py [--stage deconv|cluster|decay|spectrum|null]      # Altiverb 246 vs Tomkin (needs irs/Altiverb EMT256 IRs/)
.venv/bin/python tools/analog/compare_kernel.py [--kernel <npz> --tag <name>]            # analogue chain vs a kernel
.venv/bin/python tools/model_s.py [--tighten|--tail|--diffuser|--allpass|--groups|--stochastic]
# the S baseline (model_s_results.md §22.9; current: §25.7 with --hfcells), on the chain-kernel IRs (caches: EMT_KERNEL=chain2 tools/q1_generalise.py --stage data)
EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --stage constants                 # fold constants (needs loop_fir's fits and grids)
for rt in 0.3 0.6 0.8 1.0 1.4 1.8 2.2 2.5; do EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --stage run --rt $rt; done
EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --stage report                    # → analysis/loop_cv/report.log
EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --lf04dc --stage constants|run|report   # + LF 0.4's DC 0.414 (model_s_results.md §25.4) → analysis/loop_cv/lf04dc/
EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --cells --stage constants|run|report    # + the DC floor per cell (§25.6; needs lowband_shelf's floor grid, --gamma-nyq) → analysis/loop_cv/cells/
EMT_KERNEL=chain2 .venv/bin/python tools/loop_cv.py --hfcells --stage constants|run|report  # + the HF pole per cell (§25.7; needs lowband_shelf's hfcell grid) → analysis/loop_cv/hfcells/
EMT_KERNEL=chain2 .venv/bin/python tools/lowband.py --stage score [--variant tfix]      # low band, one file per setting (tfix needs transfer + tfit)
EMT_KERNEL=chain2 .venv/bin/python tools/lowband.py --stage report [--variant tfix]     # → analysis/lowband/chain2_<variant>_report.log
.venv/bin/python -m pytest tests -q
```
