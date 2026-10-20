# EMT 246: status and next steps

This is the hand-over for the next session. Read it with `EMT246_unit_facts.md` (what we know) and `model_s_results.md` §8–13 (where the model stands).

*Written October 2026, at the end of the single-IR phase.*

---

## 1. Where we are

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

## 2. When the Tomkin files arrive

The library will contain more than 2,400 settings: S/M/L RT plus complete Non-Lin. It comes at **96 kHz and 48 kHz, captured separately** (confirmed by Tomkin), and has **no low-pass settings**.

### 2.1 Intake (do first)

1. Put the files under `irs/tomkin/96k/` and `irs/tomkin/48k/`. They are third-party: check the EULA, and add `irs/tomkin/` to `.gitignore` unless the licence allows it in a private repo.
2. Run `.venv/bin/python tools/harness/run.py irs/tomkin/96k --templates` and add a naming pattern to `tools/harness/config.json` → `naming`.
3. **Harness tweak, not yet done:** add a config switch, e.g. `"separate_rate_captures": true`. With it on, same-label cross-rate pairs are classified as `repeat` rather than `repeat_or_resampled`, so they feed the modulation tests (spec A7 → E1). The code is in `tools/harness/library.py` `classify_pair` and `pair_flags`.
4. Run the harness on each rate (`--name tomkin96`, `--name tomkin48`). Review `analysis/library/<name>/report.md`:
   - flags;
   - normalisation (A6);
   - the latency regression (B3);
   - clock and L/R offset (B1/B4);
   - the unit fingerprint against the CLA unit (B5);
   - step tables (C3);
   - coverage (C5).
5. Add the M and L loop sets to `config.json` → `known_loopsets` once seen. They appear as `U:<lengths>`.

### 2.2 Questions to settle (each needs library IRs)

| # | Question | Test |
|---|---|---|
| Q1 | **Does Path B generalise?** Keep the S structure (loops, diffuser, cluster, injection), change only the loop gains and filters, and predict other RT/HF/LF steps | Held-out null against library IRs at other settings. This is the decisive test |
| Q2 | **Per-loop injection/read-out structure** (the −23 dB plateau) | Joint fit across IRs that differ only in RT/HF/LF, with injection, diffuser and cluster shared. The extra constraint should expose the structure |
| Q3 | **Loop sets and cluster patterns for M and L** | Harness D1/D3, loop scan, cluster fingerprint. The 1986 review says S/M/L differ "mainly in the temporal spacing of the early reflections", so expect different cluster spacing per program. Check the memory budget (~39 k words) |
| Q4 | **Control laws and step tables:** RT (16 steps per program), LF× (16), HF× (10), AMPL (11), DELAY | Harness C3 plus per-harmonic loop-filter fits. Seen so far: RT 1.0/3.6/~5/15; LP 2.0/8.0/10/12; AMPL 20/~50/80; DELAY 10/120 |
| Q5 | **LF reference: 50 Hz or 250 Hz?** (the data sheet conflicts) | LF sweeps: find where the factor applies exactly |
| Q6 | **HF ×1.0: is there residual fixed damping?** | HF sweep including the max step |
| Q7 | **Are the loop filters shared or length-scaled?** First order, or higher? | HF ×0.1 and LF ×2 extremes (`model_s_results` §3.2) |
| Q8 | **Is the 460 stage parallel or series? Does its gain vary with any control?** | Across RT and HF steps |
| Q9 | **Exact coefficients (Q15 snapping)** for bit-exactness | Joint fits across many IRs. One IR resolves loop gains to about ±8 LSB; the gain law across RT steps should pin the table |
| Q10 | **Fixed-point behaviour** (16-bit storage, truncation or rounding) | Only once model error nears −60 dB |
| Q11 | **Absolute timing:** pre-delay = DELAY exactly? | Harness B3 regression, if delays vary |

### 2.3 Modulation tests (hypothesis only; there may be none)

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
