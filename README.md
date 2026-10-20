# REV-246: EMT 246 reverb reverse-engineering workspace

Goal: recreate the EMT 246 digital reverb's S/M/L programs (Non-Lin as a bonus) as an efficient, ideally bit-accurate model, working **only from impulse responses and documents**. There is no hardware access.

**Start here:**
1. `docs/NEXT_STEPS.md`: status, the Tomkin/Kienast IR library intake, and every open check item.
2. `docs/EMT246_unit_facts.md`: everything known about how the unit works, with confidence and sources.

## File map

### Docs (`docs/`)

| File | What it is |
|---|---|
| `NEXT_STEPS.md` | Hand-over: status, library intake checklist, open questions Q1–Q11, modulation tests, Non-Lin plan, Tomkin questions, UAD 250 test bench, plugin notes |
| `EMT246_unit_facts.md` | Consolidated fact sheet: hardware, converters, controls, programs, measured program-S algorithm, budgets, open questions, source index |
| `model_s_results.md` | Program S grey-box model: 13 rounds of null tests. §8–13 hold the current state: hold-out, capture domain, long tail, diffuser, group-sparse, time-variance |
| `EMT246_discovery.md` | The original discovery and viability study. Hardware photos analysis in §2.3a; period sources in §10a; history in §10 |
| `algorithm_research.md` | Reverb-algorithm lineage (EMT 250 patent, Blesser), classic and modern methods, topology hypotheses and tests (T1–T10), the stochastic-control analysis (§7.2b) |
| `analysis_harness_spec.md` | IR-library harness test specification (stages A–G), plus implementation status of A–C |

### Code (`tools/`, `tests/`)

| Path | What it does |
|---|---|
| `tools/analyze_ir.py` | Single-IR deep dive → `analysis/<name>/{summary.json, report.png, h_dev.npy}` |
| `tools/deconv_grid.py` | Recovers the device-rate IR (19.2 MHz/728 ≈ 26,373.6 Hz) from a capture |
| `tools/topology_tests.py` | Topology tests T1/T2 against synthetic controls |
| `tools/model_s.py` | Program S model and studies. Stages: default, `--tighten`, `--tail`, `--diffuser`, `--allpass`, `--groups`, `--stochastic` |
| `tools/analog/` | Analogue chain: `netlists.py` (input 7 246 131 and output 7 246 130 filters transcribed from the service manual), `mna.py` (AC solver, ideal or finite-GBW op-amps), `compare_kernel.py` (vs the reference kernel → `analysis/analog/`) |
| `tools/harness/` | IR-library harness, stages A–C: `run.py` (CLI), `core.py` (per file), `library.py` (library level), `naming.py`, `report.py`, `synth.py` (synthetic test library), `config.json` (constants, templates, naming patterns) |
| `tests/test_harness.py` | End-to-end harness tests: `.venv/bin/python -m pytest tests -q` |

### Data

| Path | Contents |
|---|---|
| `irs/EMT 246 captures/` (gitignored; Tomkin Mastering, no redistribution of the files) | The Tomkin/Kienast library: 48 kHz and 96 kHz, `Short RT`, `Med RT`, `Long RT`, `Non-Lin`, plus its README, the operation card scan and a front-panel photo. Naming scheme in `NEXT_STEPS.md` §2.1 |
| `irs/eugtone_EMT246_CLA_shortRT.wav` | The one hardware IR so far: CLA's setting, captured by eugtone (S, LP 12, DELAY 10, AMPL 80, LF 1.0, RT 1.0, HF 0.9). Third-party |
| `ref/kernel_eugtone_96k.npz` | Reference analogue-chain kernel (one cluster tap) |
| `analysis/shortRT/` | `analyze_ir.py` output for the CLA IR. `h_dev.npy` is regenerated, not committed |
| `analysis/model_s/`, `analysis/model_s_*.log` | Model study outputs (JSON, params, plots, logs) |
| `analysis/library/synthetic/` | Harness run on the synthetic test library |
| `sources/` (local only, gitignored) | Third-party primary sources: data sheet, 252 brochure, catalogues, 245 service manual, **246S service manual (scans)**, EMT Courier, EMT 250 patent, Studio Magazin pages, magazine extracts, studiorepair internal photos. See `sources/README.md` |

## Setup and common commands

```bash
python3 -m venv .venv && .venv/bin/pip install numpy scipy soundfile matplotlib pandas pytest
.venv/bin/python tools/analyze_ir.py irs/<file>.wav --out analysis/<name>      # one IR, deep dive
.venv/bin/python tools/harness/run.py <library_dir> --templates              # inspect naming
.venv/bin/python tools/harness/run.py <library_dir> --name <lib>             # whole library, stages A–C
.venv/bin/python tools/harness/run.py "irs/EMT 246 captures" --name tomkin   # both rates (they are the same capture; do not pass --separate-rates)
.venv/bin/python tools/harness/crossrate.py "irs/EMT 246 captures" --out analysis/library/tomkin/crossrate.csv   # direct 48k vs 96k comparison
.venv/bin/python tools/harness/run.py "irs/EMT 246 captures" --name trial --separate-rates --sample 60   # quick spread-out subset
.venv/bin/python tools/model_s.py [--tighten|--tail|--diffuser|--allpass|--groups|--stochastic]
.venv/bin/python -m pytest tests -q
```
