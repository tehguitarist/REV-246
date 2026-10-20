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

## 3. Harness results (96 kHz set: `analysis/library/tomkin96/report.md`)

- **S clusters** found in 384 files (101 with exactly the S loop set `1321-1389-1461-1587`, so the library unit runs the same firmware as the CLA unit). 15 files have a `novel_cluster` (likely Non-Lin or L). M and L files have no cluster match only because `cluster_templates` holds the S pattern alone.
- **116 unknown loop sets**, the M/L candidates. The commonest: `1108-1254-1386` (310 files), `2641` (37), `231-628-803-1108-1292-2032` (37), `628` (35), `628-1108-1386` (28), `5279` (26), `1321-1389-1461` (26, a partial S set), `5279-10826-12176…` (15+15). Not yet assigned to programs.
- **Kernel (B2):** the 96k library kernel (8 files, `analysis/library/tomkin/kernels/`; the 48k kernel there is the same capture rendered at 48k) sits within 0.64 dB of the CLA kernel over 100 Hz–11 kHz. The analogue chain around them is modelled in `tools/analog/` (`EMT246_unit_facts.md` §2b).
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

Tomkin said the library was recorded twice. The data say the two rate sets are renderings of one recording (**S**). Direct comparison in the capture domain (`tools/harness/crossrate.py`, 2,437 pairs):

- Every pair is sample-aligned (lag 0.0) with a constant gain of 0.9994 and agrees to **−67.5 dB median in the first 0.3 s**.
- The out-of-band tail (14.5–22 kHz after 0.5 s, about −150 dBFS) is **correlated between the rates** (median 0.79, minimum 0.37, above 0.9 in 974 pairs). Independent captures would have independent noise there.
- The harness's cluster-tap time for a setting differs between the rates by a constant 0.0466 µs with a scatter of **0.09 ns**, while the same quantity scatters by **2–27 ns between different settings** recorded separately: the rates share their estimation error.
- For 300 random pairs, the 48k file over the decimated 96k file has one fixed response in every band: |H| std 0.001–0.002 dB, phase std 0.01–0.03 mrad.
- The 96k file is a flat 24-bit floor from 14 kHz to 47 kHz, as is the 48k file: neither carries interface noise above the rounding floor.

**What this cannot show:** which file came from which, because the 96k file has no content above 14 kHz. Either way the 96k set contains all the information. The earlier hope of using 48k/96k pairs as repeat captures for time-variance tests is gone.
