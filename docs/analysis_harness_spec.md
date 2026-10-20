# EMT 246 analysis harness: test specification

**Goal:** turn a large third-party IR library (primarily Tomkin Mastering, 2,400+ settings, 96 kHz) into a verified structural model of programs S/M/L. There is no ROM (the service manual's board schematics are now available but do not contain the algorithm) and no mainframe schematics, so the IRs are the only ground truth. The harness has two jobs:

1. **Distrust the data first:** check labels, normalisation, alignment, clocks and capture artefacts.
2. **Then extract the structure,** and prove it holds across *all* settings, not just the ones it was fitted on.

Each test lists the question it answers and its pass/flag criterion. Tests marked ★ are new since the single-IR study; the rest extend what `tools/analyze_ir.py` already does.

Output: one row per IR in `analysis/library/<name>/files.csv` (all measurements and flags), per-program fit reports, and a model-validation report that can be re-run as a regression test.

---

## A. Ingest and data QC (per file)

| # | Test | Question | Flag if |
|---|---|---|---|
| A1 | Format, rate, length, clipping, DC | Usable file? | Clipped samples; DC > −60 dBFS; length < 1.5 × T30 |
| A2 ★ | **Fade / truncation detection.** The EDC tail against the measured noise floor, and the noise floor's own flatness over its final 200 ms | Did the vendor fade or trim the tail? A fade biases RT short | Noise-floor energy falls steadily toward the end of the file |
| A3 ★ | **Pre-onset energy.** Energy before the first tap, compared with the noise floor | Sweep-deconvolution artefacts: harmonic-distortion products appear *before* t = 0 in ESS captures | Pre-onset energy > noise + 6 dB |
| A4 ★ | **Channel-swap and polarity fingerprint.** Match the device-rate cluster against the known L {0,54,194,288 / ++ + −} and R {0,143,194,353 / + + − +} patterns | Are L and R swapped, or polarity inverted? | Pattern matches the other channel or is sign-flipped. Needs amp > 0 |
| A5 ★ | **Duplicate / near-duplicate detection.** Hash, plus normalised cross-correlation of device-rate IRs | Mislabelled copies, which are common in libraries | Correlation > 0.999 between files with different labels, or < 0.99 between files that claim to be the same setting |
| A6 ★ | **Normalisation detection.** Compare absolute cluster-tap level across IRs that share the same amp step | Did the vendor normalise each file? If so, absolute gain laws must be measured as *ratios within an IR* (cluster against reverb, LF against MF), not as levels | Cluster level for an identical amp setting varies by > 0.5 dB across files |
| A7 ★ | **96 / 48 / 44.1 kHz versions of the same setting.** Resample the 96 k file down and null it against the 48/44.1 files | Are they separate captures or resampled copies? Separate captures are **free repeat measurements** for time-variance testing (E1) | — informational |

## B. Calibration (per file, then per library)

| # | Test | Question | Method |
|---|---|---|---|
| B1 | **Clock ratio.** Estimate the capture-to-device clock from cluster spacings (when amp > 0) and from loop-harmonic fits (always) | Converts everything to exact device samples | Report ppm per file. Different interfaces and units give different offsets. Flag any file whose two estimates disagree by > 50 ppm |
| B2 | **Converter kernel per library.** Use the first tap of a high-amp IR | Removes the capture chain plus the 246's analogue output | Compare the Tomkin kernel with the eugtone kernel. What they share is the 246's own output filter (to emulate); the difference is the capture interface (to discard) |
| B3 ★ | **Latency / alignment recovery.** Regress the first-tap time against the pre-delay setting across IRs | Fixes the unknown file alignment. The slope must be 1.000 (else the delay steps aren't 1 ms); the intercept is the capture latency plus any trim | Needs ≥ 3 different pre-delay settings. If every file is trimmed to its onset, report that, and fall back to relative timing only (fine for the model) |
| B4 | **L/R half-sample offset** | Is the 18.6 µs offset universal (a design property) or specific to one unit? | Constant across the library → emulate it |
| B5 ★ | **Unit fingerprint.** Compare noise floor, clock ppm and L/R offset between the eugtone unit and Tomkin's unit | Unit-to-unit variation, and the possibility of different firmware revisions | Loop lengths or cluster positions that differ between units mean firmware differs, so treat them as separate models |

## C. Settings inference and metadata audit

| # | Test | Purpose |
|---|---|---|
| C1 | Parse filenames into (program, RT, LF, HF, delay, amp, LP) | Initial labels |
| C2 | **Measure the same quantities independently.** Program from the loop set; RT from T30 at 1 kHz; LF factor from T30(50–63 Hz) ÷ T30(1 k); HF factor from T30(7 k) ÷ T30(1 k); amp from cluster-to-reverb energy ratio; delay from the first tap | Independent labels |
| C3 | **Snap the measurements to step tables.** RT 16 steps per program, LF 16, HF 10, amp 11; show the clusters in each measured distribution | Recovers the actual display step tables (e.g. is RT log-spaced?). Measured values that fall *between* clusters expose mislabels or bad captures |
| C4 | **Label-against-measurement disagreement report** | Every disagreement is either a vendor labelling error or a model assumption that is wrong; review them by hand |
| C5 ★ | **Coverage map.** Which (program × parameter) cells exist | Shows what Tomkin lacks (low-pass is already known to be missing) → a shopping list for one capture session |

## D. Structural identification (per program, across all its IRs)

| # | Test | Question | Expected if the simple model is right |
|---|---|---|---|
| D1 | **Loop set per IR**, from autocorrelation plus harmonic check | How many loops, what lengths, per program? Do M and L have more or longer loops ("adapted density")? | The same set for every IR of a program. The sum of lengths fits the 64 k-word memory budget once max pre-delay (26,347 words) is subtracted |
| D2 ★ | **Loop-set invariance across RT, LF, HF, amp and delay** | Do the lengths change with settings? | Invariant. If lengths change with RT, the model needs per-step tables |
| D3 | **Cluster pattern per program** | Fixed pattern ("feste Zuordnung") across S/M/L? | Identical positions and signs; only the gain changes with amp |
| D4 ★ | **Cluster → reverb coupling.** Difference IRs: amp 100 % minus amp 0 %, all else equal | Does the cluster feed the loops, or only the output? | Difference is cluster taps only → parallel. A difference that contains a tail → cluster feeds the reverb |
| D5 ★ | **Pre-delay scope.** Difference between delay settings | Does pre-delay move the cluster *and* the reverb together? | A pure time shift of the whole IR |
| D6 | **Per-loop excitation patterns** (periodic decomposition), then correlate each loop's pattern across RT steps | Is the injection fixed, with RT only scaling loop gains? | Pattern correlation > 0.95 across RT steps → the injection network is RT-independent |
| D7 | **Per-loop output tap offsets for L and R**, across IRs | Fixed output tap positions? | Constant lags per loop |
| D8 ★ | **First-pass analysis.** The window from reverb onset until the first recirculation (about 1,300 samples). Run a sparse fit with greedy search over integer all-pass delays and gains | What is the diffusion network: an all-pass chain, multi-tap injection, or something else? | A sparse solution with few stages → structural model. A dense, irreducible solution → measured FIR injection (Path B) |
| D9 ★ | **Per-harmonic decay rates.** Every harmonic of every loop is known exactly; fit each one's decay by tracking at fixed frequency or with ESPRIT | The loop filter magnitude \|Dᵢ(f)\| at hundreds of frequencies, per LF/HF step | Smooth curves that one filter family fits (e.g. one-pole LP plus first-order low shelf) |
| D10 ★ | **Per-harmonic phase / effective length.** Effective loop length against frequency | Loop filter phase → filter order and type | e.g. +0.2 samples at LF tapering toward HF = one-pole |
| D11 ★ | **Gain law.** Fitted gᵢ(1 kHz) against RT step | Is it the textbook gᵢ = 10^(−3Pᵢ / (fs·RT))? | Matched decay across loops; residuals < 1 % |
| D12 ★ | **Brochure reference points.** Check that the LF factor really applies at 50 Hz and the HF factor at 7 kHz | Validates the filter mapping against the published spec | Measured T(50 Hz)/T(1 k) equals the labelled LF factor, and so on |
| D13 ★ | **Coefficient quantisation.** Look at the fitted gains and filter coefficients for clustering on Q15 / 2⁻ⁿ grids or simple fractions | Implementation hints (TMS32010 fixed point; possible shift-and-add) that make exact matching easier | Values that snap to a grid → use the exact values |

## E. Time variance, nonlinearity, noise

| # | Test | Question |
|---|---|---|
| E1 ★ | **Repeat-capture null**, using only what libraries already contain: separate 96/48/44.1 captures (A7) or same-label duplicates. No new captures are possible | Is the system LTI? Reproduces eugtone's L-only-modulation finding: compare the L and R residuals |
| E2 ★ | **Modulation search.** Sidebands around strong isolated modes, phase wander of tracked harmonics, and mode-width excess beyond what decay explains, done separately for L and R | If modulation exists: its rate, depth, and which loop or output it is on |
| E3 ★ | **Noise floor against RT step.** In a recirculating fixed-point system, truncation noise recirculates, so the floor should rise with RT | Confirms 16-bit storage truncation and calibrates the noise model |
| E4 ★ | **Noise spectrum shape and L/R difference** | Emulates the converter and analogue noise; explains L being 6 dB noisier |
| E5 | **Level dependence**, if any library has the same setting at different levels | Nonlinearity: input clipping, 16-bit headroom |

## F. Program N (Non-Lin), only if in scope

Separate structure. Tomkin covers it fully. It needs gate/envelope-shape analysis (early RT, late RT, span) rather than loop identification. Defer until S/M/L pass G.

## G. Model validation (closed loop; run on every model change)

| # | Metric | Target |
|---|---|---|
| G1 | **Device-domain null**: residual energy against time, model against hardware IR | Cluster: exact. Path A: deep null (better than −30 dB over the first 300 ms). Path B: better than −20 dB |
| G2 | **Band EDCs and T30 per octave** | RT within ±3 %; EDC within ±0.5 dB to −40 dB |
| G3 | **Modal match**: loop harmonic frequencies and per-harmonic decays | Frequencies exact; decay within ±3 % |
| G4 | **Echo-density profile, and L/R coherence against frequency** | Curves overlap within the spread between captures |
| G5 | **Spectral fine-structure correlation**, model against hardware | Approaches the hardware's own repeat-capture value (from E1), which is the ceiling |
| G6 ★ | **Hold-out generalisation.** Fit the parameter laws on even steps, test on odd steps; fit on S and M, predict L wherever the structure is shared | Same scores as on the training set. Guards against overfitting a lookup table |
| G7 ★ | **Cross-unit check.** A model fitted on Tomkin's IRs, tested against eugtone's CLA IR | Passes G1–G4. Proves the model is of the *design*, not of one unit |
| G8 | **Listening.** Loudness-matched ABX on drums, vocals and snare, using convolution of the hardware IR as the reference | Not reliably distinguishable |

---

## What the harness can't get from Tomkin alone

- **Low-pass steps:** none in the library. To get them, either do a small capture session (8 IRs plus 1 control), or find other libraries that have them. The fallback is a brochure-based guess (2–12 kHz in 8 steps, filter type unknown), labelled as an approximation.
- **Level-dependence and repeat captures:** only if A7 finds separate captures.
- **Absolute alignment:** only if the files vary pre-delay (B3). The plugin doesn't need it.

## Build order

1. A + B + C (ingest, QC, calibration, settings inference, coverage map). After this, you know what the library really contains.
2. D1–D7 on program S, then on M and L.
3. D9–D13 (loop filters and laws), then the Path-B prototype and the G harness.
4. D8 (diffusion identification), aiming for Path A.
5. E, then F if N is in scope.

---

## Implementation status (stages A–C)

Implemented in `tools/harness/`; tested end to end by `tests/test_harness.py` (14 tests).

| Module | Contents |
|---|---|
| `core.py` | Per-file work: A1–A4, B1, B4 and the C2 measurements |
| `naming.py` | C1 parser and the `--templates` helper |
| `library.py` | A5–A7, B2–B5 and C3–C5 |
| `report.py` | Writes `report.md` and the plots |
| `run.py` | CLI |
| `synth.py` | Builds the synthetic test library |
| `config.json` | Device constants, cluster templates, known loop sets, parameter step tables, planned sets, naming patterns |

```bash
.venv/bin/python tools/harness/run.py <library_dir> --templates        # see the naming scheme first
# add a pattern for it under "naming" in tools/harness/config.json, then:
.venv/bin/python tools/harness/run.py <library_dir> --name tomkin
.venv/bin/python -m pytest tests -q
```

**Outputs** go to `analysis/library/<name>/`:
- `report.md`
- `files.csv`
- `pairs.csv`
- `label_disagreements.csv`
- `library.json`
- `kernels/kernel_<fs>.npz` (the B2 library kernels)
- `coverage_<program>_rt_x_<param>.csv`
- `labels_vs_measured.png`, `calibration.png`

Per-file results are cached, so re-running after a config or naming change is fast. Speed is about 0.5 s per file per pass with all cores in use. The run makes two passes: one with the reference kernel, then one with the library's own kernel.

**Synthetic test library:** 25 files derived from the reference IR. Faults detected:
- swapped and inverted channels;
- fade, pre-echo, mono, clipping, trimmed start;
- a 150 ppm clock skew, measured as −150 ppm;
- a mislabelled RT;
- a duplicate under a different label, and an unlabelled duplicate;
- a repeat capture;
- a cross-rate version.

The delay regression recovers slope 1.000 and the −20 ms intercept; the cluster gain law is measured to 0.02 dB.

**Deviations from the spec and known limits:**
- **A7:** across sample rates, a resampled copy and an independent capture can't be told apart by residual, because SRC ripple sits near the noise floor. Such pairs are tagged `other_rate_version`. Use same-rate repeats for E1.
- **Pair comparisons** run below about 10 kHz in the device domain. Above that, deconvolving the steep converter roll-off amplifies noise differently for each kernel.
- **LF and HF ratios** are measured in a 44–88 Hz band and a 5.9–8.3 kHz band, not exactly at 50 Hz and 7 kHz. Compare against labels through the step tables rather than expecting equality. The per-harmonic method in D9 will be exact.
- **Amplitude under normalisation:** if files look level-normalised (A6), amplitude falls back to `cluster_to_reverb_db`. That measure is confounded by the reverb settings, so C4 only compares it between files with the same RT/LF/HF labels.
- **Program identity** currently knows only the S loop set. The M and L sets will show up as `U:<lengths>`; add them to `known_loopsets` once seen.
- **`pre_onset_energy`** is present on the reference IR itself, probably capture/deconvolution residue. Expect it on every file that derives from it.
- **Presets** are out of scope for now.
