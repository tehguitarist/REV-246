# EMT 246 reverb: reverse-engineering workspace

- `docs/EMT246_unit_facts.md`: **consolidated facts about how the unit works** (hardware, controls, programs, measured algorithm, open questions)
- `docs/EMT246_discovery.md`: discovery and viability study
- `docs/model_s_results.md`: program S grey-box model and null tests (`tools/model_s.py`)
- `docs/analysis_harness_spec.md`: test specification for analysing a full IR library
- `docs/algorithm_research.md`: reverb-algorithm lineage (EMT 250 patent, Blesser), classic/modern methods, ranked S/M/L and N topologies, recreation plan
- `tools/topology_tests.py`: topology tests T1/T2 against synthetic controls
- `tools/analyze_ir.py`: IR analysis pipeline → `analysis/<name>/{summary.json, report.png, h_dev.npy}`
- `tools/harness/`: IR-library harness, stages A–C (ingest/QC, calibration, settings inference, coverage); see `docs/analysis_harness_spec.md`
- `tests/`: end-to-end harness tests (`.venv/bin/python -m pytest tests -q`)
- `tools/deconv_grid.py`: recovers the unit's device-rate (19.2 MHz/728 ≈ 26,373.6 Hz) IR from a 96 kHz capture
- `irs/`: source IRs (third-party; check licences before sharing)

```bash
python3 -m venv .venv && .venv/bin/pip install numpy scipy soundfile matplotlib pandas pytest
.venv/bin/python tools/analyze_ir.py irs/<file>.wav --out analysis/<name>          # one IR, deep dive
.venv/bin/python tools/harness/run.py <library_dir> --templates                  # inspect naming
.venv/bin/python tools/harness/run.py <library_dir> --name <lib>                 # whole library
```
