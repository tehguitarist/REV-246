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
| `model_s_results.md` | Program S grey-box model: 13 rounds of null tests. §8–13 hold the current state: hold-out, capture domain, long tail, diffuser, group-sparse, time-variance |
| `ml_structure.md` | Programs M and L from the library: clusters and onsets (S's ×2 / ×4), the five-term late recursion (§7: S's four loops scaled ×2 / ×4 and coupled; the S carry-over test), lag precision, the time variance found with AMPL pairs as repeat captures, and the noise floors and repeat-capture ceiling (S about −57 dB) that bound the null target |
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
| `tools/ml_lines.py` | The M/L line gains and filters (`ml_structure.md` §8): the §7 model closed through the coupled lines with S's laws at the M/L line lengths; held-out recirculation nulls against repeat captures, limit scans, decay per octave, and the law fitted per setting. One checkpoint file per AMPL group (resumable; `--chunk i/n`) → `analysis/ml_lines/` |
| `tools/ml_carryover.py` | The S carry-over test for M and L (`NEXT_STEPS` §5 item 5): diffuser search by matched-filter sparsity, the 20 read taps per channel in every AMPL group, each tap's recirculation through the four lines (the coupling matrix), and held-out first-pass nulls against repeat captures (with the AP damping step and the M/L output kernel) → `analysis/ml_carryover/` |
| `tools/q1_generalise.py` | Q1: Path B across the library (loop fits, control laws, AP230 step, shared-injection joint fit, held-out splits, repeat-capture ceiling, CLA check) → `analysis/q1/` |
| `tools/q2_injection.py` | Q2: fade-free windows, Tomkin fade compensation (fixed 0.700 s ramp), per-setting free injections, q1's joint bench with corrected windows, fade-compensated nulls, residue regression, the `branches` bench (shared injection with branch exponents and window length), the `filters` / `loops` stages (loops refitted per setting with an injection frozen), the `tied` / `diffuser` / `ap230` / `capped` stages (the compact 5-tap model, −1/√2 diffuser, AP230 steps, saturating laws), `cla` (cross-unit), and `final` (the consolidated S model and its baseline, `model_s_results.md` §17) → `analysis/q2/` |
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
| `analysis/ml_lines/` | `ml_lines.py` output: `law/` (one `.json` per AMPL group: nulls per window and variant, clocks, gains, band decay), `report.log`, `fit.json` / `fit.log` |
| `analysis/ml_carryover/` | `ml_carryover.py` output: `diffuser`, `taps` (per-group tap positions and values), `routing` (recirculation gains, the coupling matrix) and `firstpass` (held-out nulls, fitted damping steps, output kernels), each `.json` + `.log` |
| `analysis/ml_structure/` | `ml_structure.py` output: `summary.json`, `run.log` (the device-rate IR cache is not committed) |
| `analysis/q2/` | `q2_injection.py` output: corrected windows (`fade2.json`), raw-end map for the fade ramp (`ramp.json`), per-setting free-injection results (`free.json`; the injections in `v/` are not committed), joint splits (`split_*.json`, `v_*.npy`), fade-compensated nulls (`nulls_*.json`), branch runs (`branches_*.json`; the frozen flat injection `branches_flat_0_ext1700.npy`), loop refits (`filters.json`), residue regression (`regress.json`; the residue arrays are regenerated), the candidate tests of `q2_candidates.py` (`k1_sparsity.*`, `k2_stages.*`, `structure.*`, `k5_resid.*`, `settings.*`), logs |
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
.venv/bin/python tools/make_sweep.py                                                         # sweep, inverse and impulse stimuli → sweep_request/
.venv/bin/python tools/altiverb_compare.py [--stage deconv|cluster|decay|spectrum|null]      # Altiverb 246 vs Tomkin (needs irs/Altiverb EMT256 IRs/)
.venv/bin/python tools/analog/compare_kernel.py [--kernel <npz> --tag <name>]            # analogue chain vs a kernel
.venv/bin/python tools/model_s.py [--tighten|--tail|--diffuser|--allpass|--groups|--stochastic]
.venv/bin/python -m pytest tests -q
```
