# EMT 246: status and next steps

This is the hand-over for the next session. Read it with `EMT246_unit_facts.md` (what we know) and `model_s_results.md` §8–13 (where the model stands).

*Written October 2026, at the end of the single-IR phase.*

---

## 1. Where we are

- **Unit topology:** **mono in, stereo out**. One input (block diagram, analogue input board), two outputs from one multiplexed DAC. Each stereo IR is the complete 1×2 response.
- **Viability:** confirmed. Program S is a small, deterministic, integer-delay design, and almost all of its structure has been identified from one IR (the CLA "Short RT" capture by eugtone):
  - device sample rate 19.2 MHz/728 = 26,373.6 Hz;
  - the 8-tap reflection cluster, starting at the DELAY time;
  - 4 integer loops (1321/1389/1461/1587) with gains length-scaled to RT and a first-order HF filter at 7 kHz;
  - an all-pass diffuser AP20 · AP72 · AP230 (g ≈ −0.70) plus a 460-sample comb;
  - left/right read-out lags of 26/332 samples;
  - a downstream low-pass, one multiplexed DAC (L half a sample late), and a sub-20 Hz analogue high-pass.
- **Nulls** (details in `model_s_results.md` and the summary below):

  | Window | Null | Note |
  |---|---|---|
  | Cluster / first 50 ms | −50 / −52 dB | In-sample |
  | 400–600 ms | **about −28 dB held out** | |
  | R tail after 0.45 s | within 4–6 dB of the noise floor | Loops + 460 only, no fitted input |
  | L tail | 15–22 dB above noise | Fits about 1 % slow gain jitter on L's read-out (unit or capture chain) |

- **Main gap:** the per-loop injection/read-out structure. Compact fits plateau at about −23 dB where a matching synthetic reaches −46 dB.
- **Ruled out:**
  - coupled FDNs;
  - all-passes inside the loops;
  - delay-length modulation (EMT 250 style);
  - randomised injection taps;
  - a separate late-input generator (the apparent late input was the diffuser tail plus deconvolution edge ringing);
  - a long EMT 250-style loop;
  - a shared input or a shared small tap set across loops.
- **Nothing left to answer from current data.** Every open question below needs library IRs (or optional plugin renders). The single IR has been taken as far as is useful: 13 model rounds, documented in `model_s_results.md`.

---

## 2. The Tomkin / Kienast library (arrived)

The library has 2,437 settings per rate (so 4,874 files; this matches Tomkin Mastering's "2,400+" and the README credits Thomas Kienast: **H** that these are the same product): S/M/L RT plus Non-Lin. It comes at **96 kHz and 48 kHz, captured separately** (confirmed by Tomkin), and has **no low-pass settings**.

### 2.1 Intake (library has arrived; steps 1–3 done)

The library is in `irs/EMT 246 captures/` (6 GB, 4,874 stereo 24-bit WAVs, gitignored). It is **third-party, the Tomkin Mastering library (confirmed by the user). Licence: commercial use is fine; the only restriction is no redistributing the files as they are. README: "(c) Thomas Kienast. Please do not reduplicate"**. It stays gitignored; derived data (parameters, fitted models) is fine to commit. Layout:

```
irs/EMT 246 captures/
  README.txt, emt246_front.jpg, emt246_operation_and_function.jpg
  48kHz/ and 96kHz/            identical file sets (same capture rendered at two rates; see below)
    Short RT/ (512)  Med RT/ (519)  Long RT/ (360)  Non-Lin/ (1,046)
```

**Naming** (README; fields run left to right in front-panel order, and the last three have an implied decimal point):

| Folder | Pattern | Example → setting |
|---|---|---|
| Short / Med / Long RT | `LP_DELAY_AMPL_LF_RT_HF_dc.wav` | `12_0_99_12_08_09` → LP 12 kHz, DELAY 0 ms, AMPL 99 %, LF×1.2, RT 0.8 s, HF×0.9 |
| Non-Lin | `LP_DELAY_SPAN_EARLY_LATE_dc.wav` | `12_0_05_-30_20` → LP 12, DELAY 0, SPAN 0.5 s, EARLY RT −3.0 s (negative = reversed, **H**), LATE RT 2.0 s |

- Non-Lin has five numbers, S/M/L six, so the two patterns don't need the folder. Program S/M/L comes from the folder name only, so run the harness from the library root or a `<rate>/` folder, not from inside a `* RT` folder (the program would then come from the measured loop set instead).
- The Non-Lin field order is inferred (README says "left to right"; the panel photo shows `N`, 12 kHz, delay 0, span 0.4, early −3.0, late 1.2). **H** until a loop-set or EDC check confirms it.
- `dc` is probably "DC removed" (**H**).
- Value sets (both rates): LP 12 and DELAY 0 always. **So there are no LP or DELAY steps.** AMPL: S 0/30/60/90, M 10/60/99, L 0/80. LF: S 0.1/0.4/1.0/2.0, M 0.1/0.4/0.7/1.2/2.0, L 0.1/0.8/1.2/2.0. HF: S 0.2/0.4/0.6/1.0, M 0.1/0.2/0.4/0.7/0.9, L 0.1/0.3/0.5/0.8/1.0. RT: S 0.3–2.5 (8 values), M 0.6–5.0 (8), L 2.2–20 (9). Non-Lin: SPAN 0–0.5, EARLY −3.0…10, LATE 0–2.0. This is a sparse grid, not the full 16-step tables.
- Files are variable length (about 0.7–11 s), so they were trimmed per setting. Expect `zero_tail` / `truncated_before_noise_floor` flags and check the pre-roll (B3).

**Done in the harness** (`tools/harness/config.json` → `naming`: `kienast-sml`, `kienast-nonlin`; all 4,874 files parse):
- Naming patterns support `scale` (implied decimals) and `set` (constants). New label fields `span`, `early`, `late` take part in same-label pairing.
- `--separate-rates` (or `"separate_rate_captures": true` in config) classifies same-label cross-rate pairs as `repeat`, so they feed the modulation tests (spec A7 → E1). The default stays off for the synthetic test library, which contains a resampled 48k copy.

**Trial run (60 settings × 2 rates = 120 files, `--sample 60 --separate-rates`, 19 s, `analysis/library/trial60/`).** The pipeline works end to end. Observations to carry into the full run (the first two are data properties, the last three are harness gaps):
- **Cross-rate pairs are not all independent captures.** Of 60 same-setting 48k/96k pairs: 39 `repeat` (median residual L −64.7 / R −63.2 dB), 11 `copy` (−62 to −77 dB: the 96k file looks like a resampled 48k one, e.g. six Med RT settings and five Non-Lin ones), 3 `near` and 7 `different`. Only genuine repeats feed the modulation tests (E1), and `repeat` vs `copy` is separated only by the residual level, so treat the split as **H** until the full run shows a clear gap.
- **Files are trimmed and faded per setting.** 97/120 end in digital silence, 85/120 are flagged `fade_out`, 90/120 `short_file`, 63/120 have under 2 ms pre-roll. The long-RT labels therefore disagree with the measured RT (L RT 10 reads 4.2 s), so `label_mismatch_rt` is mostly truncation, not mislabelling. RT comes from the label for fits.
- **S clusters are found in 20/20 files with AMPL > 0; M and L never match** because `cluster_templates` has only the S pattern. Expected, and it is the M/L cluster discovery task (Q3). Also 22 unknown loop sets turned up, to review in the full run.
- **Non-Lin with negative EARLY RT does not align across rates.** All 7 `different`/`near` Non-Lin pairs are reversed settings, with the best lag pinned at the ±40-sample search limit. The onset anchor is unsuited to a rising envelope; needs a wider search or another anchor before E1 can use them.
- **2 Long RT pairs are `different` with a −26 dB gain** (RT 10 and 6, AMPL 0). Cause not yet found.

**Full run done** (`analysis/library/tomkin/`, 4,874 files, about 11 min with a warm cache; `report.md` has the detail). Headlines:
- S clusters found in 768 files (190 with the exact S loop set), 31 `novel_cluster` files, 141 unknown loop sets (M/L candidates: e.g. `1108-1254-1386`, `2641`, `5279-10826-12176`, `231-628-803-1108-1292-2032`).
- Kernel B2: both rates agree (max deviation 0.64 dB, spread under 0.02 dB); L/R offset 19.20 µs in every S file; capture clock −15.2 / −13.2 ppm (48k/96k; a constant 2.0 ppm apart in every pair); 30 clock outliers.
- `label_mismatch_rt` on 2,381 files: mostly truncation and fades (the files are trimmed and faded per setting).
- **The pair classes in `report.md` (`copy` / `repeat` / `near` / `different`) are not reliable for this library.** `noise_ref` is meaningless where tails are digitally zeroed, and the ±40-sample alignment window on the short device-rate signature fails for ambiguous onsets (reversed Non-Lin, low-AMPL files). They are superseded by `tools/harness/crossrate.py`, which compares the two files directly in the capture domain (`analysis/library/tomkin/crossrate.csv`).

**Cross-rate result (all 2,437 settings, 48 kHz against 96 kHz decimated, `crossrate.csv`):**
- Every pair is sample-aligned (lag 0.0), has gain 0.9994 (a constant −0.005 dB), and agrees to **−67.5 dB median in the first 0.3 s** (range −67.9 to −66.8 for 75 % of pairs; 5 % worse than −60 dB). Later windows read worse only where the signal has fallen to the 24-bit floor.
- **The out-of-band tail (14.5–22 kHz after 0.5 s, level about −150 dBFS) is correlated between the rates: median 0.79, minimum 0.37, and above 0.9 in 974 pairs, never near 0.** Independent captures would have independent noise there. The correlation is highest for long files (Long RT 0.99) and lowest where the files end near the 24-bit floor.
- **Conclusion (S): the 48 kHz and 96 kHz files are renderings of the same capture data, not independent captures,** and they carry no usable repeat-capture information. This disagrees with the earlier note that the rates were captured separately ("confirmed by Tomkin"); check what Tomkin actually said (they may have been captured separately and then re-derived, or the 48 kHz set is a conversion of the 96 kHz one). The constant 2.0 ppm clock offset is the only sign of a difference and is more likely a resampler artefact than a second capture.
- **Consequence:** the modulation tests in §2.3 cannot use 48k/96k pairs. The −67.5 dB agreement is only an upper bound on any time variance in the *data* (it says the files match, not that the unit is LTI). The `--separate-rates` switch should stay off for this library.
- 7 earlier "different" pairs (e.g. reversed Non-Lin) are fine: they were alignment failures, and all agree to −66 dB once aligned directly.

**Speed:** per-file cost halved (0.207 → 0.101 s) by replacing an O(n·w) `np.convolve` boxcar in `analyze_ir.edc` with a running sum. The harness already runs one process per core (10 here); threads would not help because the work is CPU-bound Python/NumPy. The full run is about 10 min of per-file work (two passes over 4,874 files) plus the pair stage.

**Still to do:**
2. Run `.venv/bin/python tools/harness/run.py "irs/EMT 246 captures" --name tomkin` (run done, without --separate-rates now; results are cached per file, so an interrupted run resumes). Use `--sample N` for a quick spread-out subset with both rates; `--limit` takes the first files alphabetically, which are all one folder.
3. Review `analysis/library/tomkin/report.md`: flags, normalisation (A6), latency regression (B3, but DELAY never varies), clock and L/R offset (B1/B4), unit fingerprint against the CLA unit (B5), step tables (C3), coverage (C5; the `planned_sets` for `lowpass` and `predelay` will show as empty).
4. Add the M and L loop sets to `config.json` → `known_loopsets` once seen. They appear as `U:<lengths>`.
5. **Images re-parsed (done):** folded into `EMT246_unit_facts.md` §3. Results: LF reference is **50 Hz** (Q5, now **S**), HF at 7 kHz, optimum RT ranges, Non-Lin pivot/span diagram, and the panel's Non-Lin layout confirming the file field order.

### 2.6 The EMT 246S service manual (text pages read; schematics not yet)

Two scanned PDFs in `sources/emt-documents/` (`EMT 246S service manual 1.pdf`, 39 pp; `... 2.pdf`, 30 pp), bilingual DE/EN. OCR text (English column) is in `sources/emt-documents/ocr/txt_m<vol>-<page>.txt` (gitignored, made with `tesseract --psm 1`; OCR of numbers can be wrong, so check figures against the page image). The findings so far are in `EMT246_unit_facts.md` §2b. What they change:

| Area | Finding | Use |
|---|---|---|
| Algorithm structure | Every delay is `base counter + TOM[tap]` into one 64 K-word memory; two synchronous TMS320s | Keep the model as integer taps of one circular buffer; tap offsets are per-setting data, so the *tap pattern* (cluster, injection, read-out) is data loaded per setting. Reinforces Q2/Q3: look for a fixed tap table scaled/shifted per setting |
| Fixed-point | 16-bit storage, 32-bit accumulate; **REVERB noise is 15 dB above DELAY** (73 vs 88 dB re nominal, unweighted rms) | Q10: simulate truncation against rounding in the loop model and compare the added noise with 73 dB. Also an absolute check on any model noise: nominal = about −6 dBFS |
| Analogue chain | Both filters transcribed and simulated (`tools/analog/`); the input filter has +5 dB HF pre-emphasis and an elliptic skirt at 12.2 kHz, the output filter is flat then elliptic. The residual is **the output S&H's exponential charging: a discrete one-pole at the device rate, pole 0.325 (CLA) / 0.301 (Tomkin), 0.11 dB rms** (manual p. 48, trimmer P202) | The analogue front end is a fixed structure with one trim parameter (the S&H pole, per unit) plus a small skirt shift. Cheap to implement (one IIR per output channel at 26.37 kHz) |
| LP setting | Through path and kernel agree (12 kHz: −2.5 vs −2.0 dB), so the earlier idea of 5 extra dB from the LP at its 12 kHz step is withdrawn | LP steps still absent from the library; the LP remains an unknown digital stage |
| Calibration | 0 dB = 0.775 V; nominal +6 dBu = −6 dBFS; overload reserve 6 dB at f < 2 kHz | Absolute levels in the IRs are lost in deconvolution and normalisation, so calibration mostly sets headroom and the noise targets |

**Caveats (agreed challenge to the plan):**
- The manual cannot reveal the DSP program (it lives in the PROMs), only the machine it ran on. Treat it as a constraint set, not as an answer.
- The analogue values are **nominal**. Trimmers (P201–P204) and ageing make each unit different, so the fit must keep a few free parameters, and "separating algorithm from analogue" works best when the analogue model is low-order.
- The IR noise floor is not the unit's noise: a deconvolved sweep has processing gain and the capture interface adds its own noise, so compare the spec (73/88 dB) with IRs only through a known sweep level, or via L-versus-R and settings-versus-settings differences.

**Processor and filter schematics (read):** details are in `EMT246_unit_facts.md` §2b. What follows for modelling:
- **Tap budget:** at most **256 TOM entries**, every memory access (read or write) is one entry, and the DSP time (182 instr per DSP per sample) allows far fewer. The count is below.
- **Run-time tap changes are possible** (the TOM can be rewritten by processor A and the tap counter reloaded). A modulated or jittered read-out tap (the L-only "1 % gain jitter" from the CLA IR) is thus a hardware-possible mechanism; the library repeats (48k/96k) are the test.
- **Fixed timing:** the 182-state sequencer places the A/D and D/A strobes; the L/R offset is quantised to 208.3 ns. The library's B4 values are already one value (19.19–19.21 µs in all 768 S files, about 92.1 steps) against the CLA unit's 18.6 µs (about 89 steps), so the two units differ by about 3 timing steps; a per-unit difference of that size is the kind of thing a different timing PROM or analogue S&H trim would produce.
- **Analogue chain:** both filters are nominal-value-checkable; see next.

**Memory-access count for the identified S structure** (each delay element = one TOM entry per read or write; counted from `EMT246_unit_facts.md` §5, **H** for the access pattern):

| Part | Accesses per sample |
|---|---|
| Input write, pre-delay read | 2 |
| Reflection cluster (4 taps per channel, 2 channels) | 8 |
| 4 loops (1321/1389/1461/1587): one read and one write each | 8 |
| 460 comb: read and write | 2 |
| Diffuser AP20 · AP72 · AP230: read and write each | 6 |
| Read-outs (lags 26 and 332 in each of the 4 loops, shared by L and R) | 8 |
| **Identified total** | **about 34** |

So the identified structure uses about 34 of the 256 TOM entries and about half of the 50–70-tap instruction budget (a read, a multiply-accumulate and a write cost several of the 182 cycles per DSP). The **TOM is not the binding limit; the instruction budget is.** What the 256 limit does settle: the sparse injection fits (1,719 integer taps per channel) cannot be the real structure, and the unidentified injection and per-loop read-out must fit in about 20–35 accesses, which is the constraint for Q2.

**Next reads and jobs (not done):**
1. **Filters: done** (see the table above and facts §2b). The lag is the output S&H (see table). Open: the input-side integrate/hold window (boxcar, board 121 p. 41–43) and the finite op-amp GBW question (ideal op-amps fit better than 3–5 MHz GBW models).
2. Sheets 2/5 and 3/5 (p. 55–56; processor A and B wiring, reset, clock sync), the S&H and DAC wiring on boards 121/122 (pp. 40–45, 47–50) for the exact aperture and hold timing, and the interface board 7 246 125 (control CPU link).
3. Test-program pages 26–31 (Tests 1–3, 25, 26) and the front pages (pdf pp. 1–7, 13) not yet read.

### 2.2 Questions to settle (each needs library IRs)

| # | Question | Test |
|---|---|---|
| Q1 | **Does Path B generalise?** Keep the S structure (loops, diffuser, cluster, injection), change only the loop gains and filters, and predict other RT/HF/LF steps | Held-out null against library IRs at other settings. This is the decisive test |
| Q2 | **Per-loop injection/read-out structure** (the −23 dB plateau) | Joint fit across IRs that differ only in RT/HF/LF, with injection, diffuser and cluster shared. The extra constraint should expose the structure |
| Q3 | **Loop sets and cluster patterns for M and L** | Harness D1/D3, loop scan, cluster fingerprint. The 1986 review says S/M/L differ "mainly in the temporal spacing of the early reflections", so expect different cluster spacing per program. Check the memory budget (~39 k words) |
| Q4 | **Control laws and step tables:** RT (16 steps per program), LF× (16), HF× (10), AMPL (11), DELAY | Harness C3 plus per-harmonic loop-filter fits. Seen so far: RT 1.0/3.6/~5/15; LP 2.0/8.0/10/12; AMPL 20/~50/80; DELAY 10/120 |
| Q5 | **LF reference: the card, data-sheet table and magazine say 50 Hz** (only a data-sheet diagram says 250) | LF sweeps: confirm where the factor applies exactly |
| Q6 | **HF ×1.0: is there residual fixed damping?** | HF sweep including the max step |
| Q7 | **Are the loop filters shared or length-scaled?** First order, or higher? | HF ×0.1 and LF ×2 extremes (`model_s_results` §3.2) |
| Q8 | **Is the 460 stage parallel or series? Does its gain vary with any control?** | Across RT and HF steps |
| Q9 | **Exact coefficients (Q15 snapping)** for bit-exactness | Joint fits across many IRs. One IR resolves loop gains to about ±8 LSB; the gain law across RT steps should pin the table |
| Q10 | **Fixed-point behaviour** (16-bit storage, truncation or rounding) | Only once model error nears −60 dB |
| Q11 | **Absolute timing:** pre-delay = DELAY exactly? | Harness B3 regression, if delays vary |

### 2.3 Modulation tests (hypothesis only; there may be none)

> **Update:** the 48k/96k pairs are the same capture rendered twice (§2.1), so the first bullet group below cannot be run on this library. What remains: similar presets one step apart, the per-program comparison, and a repeat capture from another source if one exists.

- **48k/96k pairs (separate captures of the same setting):**
  - Bring both to the device grid; remove level, clock and per-rate kernel differences; compare the residual with the noise floor.
  - How to read the residual:

    | Residual | Meaning |
    |---|---|
    | At noise | LTI |
    | At mode frequencies, flat with frequency, mode pitches unchanged | Gain/tap-weight modulation |
    | Growing with frequency and time, modes shifting | Delay modulation |
    | L only | The ~1 % L read-out jitter: constant across pairs → the unit; varying → the capture chain |
    | First ~300 ms only | Injection |
    | Persisting through the tail | Loops / read-out |

  - **Rate** comes from decorrelation across frequency bands, mapped to time via the sweep. Ask Tomkin for the sweep parameters (§3).
- **Similar presets** (one RT/HF step apart): the shared early structure should match within noise. Residuals across hundreds of pairs give the modulation's depth and location statistics (not its waveform).
- **Per program:** "Grad der Stochastik" differs between S/M/L (1986 review). Compare the L excess and any residual between programs.
- Detail: `algorithm_research.md` §7.2b and `model_s_results.md` §13.5.

### 2.4 Non-Lin (N), the bonus

- **Model:** pre-delay → early section of length SPAN (0.1–0.5 s; can be reversed) → late reverb.
- **Parameters:** probably DELAY, SPAN, EARLY RT and LATE RT; no AMPL/LF/HF.
- **Expect** the early section to be truncated abruptly at SPAN (gated reverb), and the late section to reuse the S/M/L loops and diffuser.
- **Tests:** EDC cut-off sharpness at SPAN; periodicities in the early and late sections; tap counting; whether tap positions scale with SPAN; which loop set the late part uses; how reversal works.
- References: `algorithm_research.md` §6 and §7.4; data sheet; Studio Magazin 99.

### 2.5 Gaps the library won't fill

- **Low-pass steps** (Tomkin has none). The LP is downstream; type and order are unknown. Look for LP-varied IRs in other libraries:
  - Wave Arts / Impulse Record "Vintage German DSP";
  - Past to Future CLA combos;
  - Altiverb;
  - the Foleytronics video (LP 2.0 / 8.0 / 10 / 12 settings seen).

  Failing that, approximate it.
- **Program E** (the EMT 250 program) and **program D** are out of scope.

---

## 3. Questions for Tomkin (when convenient)

1. Sweep method: type, length, frequency range and level. This is needed to map time-variance to rate.
2. Do they still have **the raw sweep recordings** (before deconvolution)? These would allow direct time-variance analysis.
3. Are the files trimmed (is absolute latency kept?) and/or normalised?
4. File-naming scheme and any settings metadata.
5. Which unit (serial or firmware, if known), and was the same unit used for every capture?

---

## 4. Optional: UAD EMT 250 plugin as a modulation test bench

The aim is to validate the modulation detectors on a reverb that is known to be modulated, and to get a "250-style" fingerprint to compare with any 246 residual. It is **not** a template for 246 S/M/L: their harmonic loops rule out 250-style delay modulation.

| Test signal (to generate) | Render through the plugin (100 % wet, offline) | Reveals |
|---|---|---|
| Long sine tones at 250 Hz, 1 kHz, 4 kHz | One render each | FM sidebands (delay modulation, grow with frequency) against AM sidebands (gain modulation); rate and depth |
| The same tones | Instantaneous-frequency tracking | The modulation waveform (random, LFO or stepped) |
| An impulse every ~10 s in a long file | One render | Snapshot IRs with no sweep smear: depth and which part varies |
| The same sine sweep, several times | Several renders | Calibrates the sweep-pair test used on Tomkin's 48k/96k |
| The same render twice after reloading the plugin | Two renders | Whether its randomness restarts identically |

**To do:** generate the signals (`tools/` script), have the user render them, and add the detectors to the harness (stage E).

---

## 5. Later: plugin implementation notes

- **Core at the native 26,373.6 Hz in integer (16-bit storage, 32-bit accumulate) arithmetic.** It costs a few % of one modern core: the 246 itself ran about 9.6 M instructions/s across two TMS32010s.
- **Host-rate conversion doubles as the analogue emulation:** the 12 kHz band-limit kernel, the sub-20 Hz high-pass, and the L half-sample offset.
- **Noise modelled separately:** analogue hiss (73 dB S/N spec), any L-specific behaviour, and truncation noise (comes free from fixed-point emulation).
- **Structure first.** Every free FIR replaced by true structure makes the model both more accurate and cheaper. The hardware budget was about 50–70 weighted taps per sample.
- **Trademark:** "EMT" belongs to EMT International GmbH. Use a non-infringing product name.
