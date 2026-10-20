# EMT 246: status and next steps

Hand-over for the next session. Read it with `EMT246_unit_facts.md` (what we know), `library_notes.md` (the Tomkin library and harness results) and `model_s_results.md` §8–13 (where the model stands).

*Updated October 2026, at the end of the library-intake phase. The model itself has not changed since the single-IR phase; what changed is the data, the hardware knowledge and the analogue chain.*

---

## 1. Where we are

- **Unit topology:** **mono in, stereo out**. One input, two outputs from one multiplexed DAC (L half a sample late). Each stereo IR is the complete 1×2 response.
- **Program S** is a small, deterministic, integer-delay design. Almost all of its structure was identified from one IR (the CLA "Short RT" capture by eugtone):
  - device rate 19.2 MHz / 728 = 26,373.6 Hz (now confirmed by the processor board's timing sequencer);
  - the 8-tap reflection cluster, starting at the DELAY time;
  - 4 integer loops (1321/1389/1461/1587) with gains length-scaled to RT and a first-order HF filter at 7 kHz;
  - an all-pass diffuser AP20 · AP72 · AP230 (g ≈ −0.70) plus a 460-sample comb;
  - left/right read-out lags of 26/332 samples;
  - the analogue chain (§4).
- **Nulls** (from one IR; details in `model_s_results.md`):

  | Window | Null | Note |
  |---|---|---|
  | Cluster / first 50 ms | −50 / −52 dB | In-sample |
  | 400–600 ms | **about −28 dB held out** | |
  | R tail after 0.45 s | within 4–6 dB of the noise floor | Loops + 460 only, no fitted input |
  | L tail | 15–22 dB above noise | Fits about 1 % slow gain jitter on L's read-out (unit or capture chain) |

- **Main gap:** the per-loop injection/read-out structure. Compact fits plateau at about −23 dB where a matching synthetic reaches −46 dB.
- **Ruled out:** coupled FDNs; all-passes inside the loops; delay-length modulation (EMT 250 style); randomised injection taps; a separate late-input generator; a long EMT 250-style loop; a shared input or small tap set across loops.
- **New since the single-IR phase:**
  - **The library is in and analysed** (`library_notes.md`): 2,437 settings, S/M/L/Non-Lin, LP and DELAY fixed. The harness run is done; 116 unknown loop sets are the M/L candidates.
  - **The 48k and 96k sets are the same capture** (`library_notes.md` §5). There are no independent repeat captures, so the time-variance question cannot be answered by pair tests.
  - **The service manual is read** (text pages, processor board, both filters, output stage): `EMT246_unit_facts.md` §2b, page map in `sources/emt-documents/EMT246S-service-manual-page-map.md`.
  - **The analogue chain is modelled** (`tools/analog/`): two filters from the manual plus one S&H pole. This replaces the free band-limit kernel.
  - **Hard hardware limits:** 256 tap-offset entries, 64 K words of delay memory, 182 instructions per sample per DSP. The sparse injection fits (1,719 taps per channel) cannot be the real structure.

---

## 2. Working rules for this phase

- **Work on the 96 kHz set only.** The 48 kHz files are the same capture and give no held-out evidence. Held-out tests must be other settings (or the CLA IR, a different unit), never the other rate.
- **RT comes from the label, not the file** (files are trimmed and faded). Absolute noise floors are not available either.
- **Report nulls** held out, with excess over the noise floor and the parameter count (`CLAUDE.md`). Beware that the library's tails are faded or gated, so "noise floor" there is processing, not the unit.
- **The harness pair classes are unreliable** (`library_notes.md` §3); use `crossrate.csv` for anything cross-rate.
- **Cross-unit check is available:** the CLA IR (S, LP 12, DELAY 10, AMPL 80, LF 1.0, RT 1.0, HF 0.9) is a different unit; 190 library files share its exact S loop set, so the firmware matches. The nearest library settings (S, LF 1.0, RT 1.0, HF 1.0, AMPL 90) are not identical, so use it as a structural check (spec G7), not a null.

---

## 3. Open questions (each needs library IRs)

| # | Question | Test / what is known |
|---|---|---|
| Q1 | **Does Path B generalise?** Keep the S structure (loops, diffuser, cluster, injection), change only the loop gains and filters, and predict other RT/HF/LF steps | Held-out null against library IRs at other settings. The decisive test. S grid: 8 RT × 4 LF × 4 HF × 4 AMPL |
| Q2 | **Per-loop injection/read-out structure** (the −23 dB plateau) | Joint fit across IRs that differ only in RT/HF/LF, with injection, diffuser and cluster shared. **Constraint:** the identified part uses about 34 memory accesses of a 50–70-tap budget, so the unidentified part has about 20–35 accesses |
| Q3 | **Loop sets and cluster patterns for M and L** | The 116 unknown loop sets in the 96k set (`library_notes.md` §3) and the 15 `novel_cluster` files; assign them to programs from the folder labels (the harness can do D1/D3). Expect different cluster spacing per program (Studio Magazin 99). Memory budget: about 39 K words |
| Q4 | **Control laws and step tables:** RT, LF×, HF×, AMPL | The library has a sparse grid, not 16 steps. Fit the laws from the values present; spec says RT 16 steps, LF× 16 (0.1–2), HF× 10 (0.1–1), AMPL 11 steps |
| Q5 | **LF reference frequency** | **S: 50 Hz** (operation card, data-sheet table, Studio Magazin; only one data-sheet diagram says 250). Confirm with the LF sweeps |
| Q6 | **HF ×1.0: residual fixed damping?** | HF sweeps including the max step (S has 1.0, M 0.9, L 1.0) |
| Q7 | **Are the loop filters shared or length-scaled? First order, or higher?** | HF 0.1 and LF 2 extremes (`model_s_results` §3.2) |
| Q8 | **Is the 460 stage parallel or series? Does its gain vary with any control?** | Across RT and HF steps |
| Q9 | **Exact coefficients (Q15 snapping)** for bit-exactness | Joint fits across many IRs; one IR resolves loop gains to about ±8 LSB |
| Q10 | **Fixed-point behaviour** (16-bit storage; truncation or rounding) | **Target: the REVERB program's noise is 15 dB above DELAY's** (73 vs 88 dB re nominal, unweighted rms; nominal = about −6 dBFS). Simulate truncation against rounding in the loops and compare the added noise; the TMS32010 has no rounding mode, so truncation is likely (**H**). Only meaningful once the model nears −60 dB |
| Q11 | **Absolute timing: pre-delay = DELAY exactly?** | Not testable from the library (DELAY is always 0); only the CLA IR has a delay (10 ms) |
| Q12 | **The LP's type and order** (downstream, digital?) | No library data (LP always 12). The through path and the kernel agree at 12 kHz, so the LP at its 12 kHz step adds little. Other libraries (§7) |

---

## 4. The analogue stage (done; use it)

`tools/analog/` simulates the transcribed input filter (7 246 131) and output filter (7 246 130) and compares them with a kernel; `EMT246_unit_facts.md` §2b has the details.

- Chain = input filter (elliptic, 12.2 kHz, with a +5 dB HF boost) × output filter (elliptic, flat to 11 kHz) × **the output S&H lag: a discrete one-pole at the device rate, `y[n] = y[n−1] + (1−p)(x[n] − y[n−1])`, p = 0.325 (CLA) / 0.301 (Tomkin), 0.11 dB rms over 0.5–11 kHz**. p is set by a trimmer (P202), so it is a per-unit parameter.
- Use it as the analogue stage when fitting the digital model, with p (and a small skirt shift, about 0.2 % in frequency) free per unit. This separates algorithm from analogue and makes the L/R and noise-floor questions cleaner.
- Open: the input-side integrate/hold window (board 121, pp. 41–43), the finite-GBW question (ideal op-amps fit better than 3–5 MHz GBW models), and the unit-to-unit meaning of p.

---

## 5. Suggested order of work

1. **M and L structure (Q3).** Assign loop sets and cluster patterns to programs. This is data-driven and unblocks everything for M/L.
2. **Path B generalisation on S (Q1)** at held-out settings, with the analogue stage from §4. If it holds, the same machinery applies to M and L.
3. **Injection/read-out (Q2)** with the access-count constraint; joint fit across settings.
4. **Control laws (Q4–Q8)** from the fits.
5. **Fixed-point noise (Q10)** and exactness (Q9) once nulls reach about −50 dB.
6. **Non-Lin (§6)** last.

---

## 6. Non-Lin (N), the bonus

- **Model:** pre-delay → early section of length SPAN (0–0.5 s; can be reversed) → late reverb.
- **Parameters:** DELAY, SPAN, EARLY RT and LATE RT; no AMPL/LF/HF. LP applies as for the other programs (operation card).
- **From the operation card** (`EMT246_unit_facts.md` §3): the early-RT line **pivots about the middle of the SPAN** (so a negative early RT is a rising ramp that reaches the late-RT start level at the end of SPAN), and the late RT starts at the end of SPAN at the same level as the early start. **H**: it is a schematic.
- **Tests:** EDC cut-off sharpness at SPAN; whether the level at mid-SPAN is independent of EARLY RT (the pivot reading); periodicities in the early and late sections; tap counting; which loop set the late part uses.
- **Library:** 1,046 settings per rate (SPAN 0–0.5, EARLY −3…10, LATE 0–2). Many Non-Lin settings differ only in LATE RT, so their first section is identical, a free consistency check. Reversed (negative-EARLY) files have a rising envelope that the harness onset anchor handles badly.

---

## 7. Gaps the library won't fill

- **Low-pass steps** (none in the library). Look for LP-varied IRs elsewhere: Wave Arts / Impulse Record "Vintage German DSP"; Past to Future CLA combos; Altiverb; the Foleytronics video (LP 2.0 / 8.0 / 10 / 12 seen). Failing that, approximate it.
- **DELAY steps** (always 0 in the library).
- **Independent repeat captures** (none: the rate sets are one capture). The time-variance question (the CLA IR's L-only "1 % gain jitter") can only be approached by joint fits across settings, by the cross-unit comparison with the CLA IR, or by a repeat from another library.
- **Program E** (EMT 250) and **program D** (delay/echo) are out of scope.

---

## 8. Modulation / time-variance (hypothesis only; there may be none)

- **Hardware makes it possible:** the tap-offset memory can be rewritten by processor A while running, so jittered or modulated read-out taps are feasible (`EMT246_unit_facts.md` §2b). Nothing in the manual links the built-in noise generator (Tests 21/22) to the "stochastic control".
- **What can be done without repeat captures:**
  - residuals across similar settings (one RT/HF step apart) share early structure; the scatter gives an upper bound on modulation depth and location (not its waveform);
  - L versus R behaviour in joint fits (the CLA IR's L-only excess against the library);
  - per program: "Grad der Stochastik" differs between S/M/L (1986 review); compare the L excess and residuals between programs.
- **How to read a residual:**

  | Residual | Meaning |
  |---|---|
  | At noise | LTI |
  | At mode frequencies, flat with frequency, mode pitches unchanged | Gain/tap-weight modulation |
  | Growing with frequency and time, modes shifting | Delay modulation |
  | L only | The ~1 % L read-out jitter |
  | First ~300 ms only | Injection |
  | Persisting through the tail | Loops / read-out |

- Detail: `algorithm_research.md` §7.2b and `model_s_results.md` §13.5. The optional UAD EMT 250 test bench (§10) still validates the detectors on a known-modulated reverb.

---

## 9. Questions for Tomkin (when convenient)

1. **What exactly was recorded twice?** The two rate sets are the same capture (`library_notes.md` §5). Was the 48 kHz set exported from the 96 kHz master, or the reverse, or are the sweeps from one session processed twice?
2. Sweep method: type, length, frequency range and level. Needed for any time-variance analysis and for noise-floor interpretation.
3. Do they still have **the raw sweep recordings** (before deconvolution and before the trimming, fading and DC removal)? These would allow direct time-variance analysis and true noise floors.
4. What are the **trimming, fade and `dc` steps** (the files are cut and faded per setting)?
5. Which unit (serial or firmware, if known), and was the same unit used for every capture?

---

## 10. Optional: UAD EMT 250 plugin as a modulation test bench

The aim is to validate modulation detectors on a reverb that is known to be modulated, and to get a "250-style" fingerprint to compare with any 246 residual. It is **not** a template for 246 S/M/L: their harmonic loops rule out 250-style delay modulation.

| Test signal (to generate) | Render through the plugin (100 % wet, offline) | Reveals |
|---|---|---|
| Long sine tones at 250 Hz, 1 kHz, 4 kHz | One render each | FM sidebands (delay modulation, grow with frequency) against AM sidebands (gain modulation); rate and depth |
| The same tones | Instantaneous-frequency tracking | The modulation waveform (random, LFO or stepped) |
| An impulse every ~10 s in a long file | One render | Snapshot IRs with no sweep smear: depth and which part varies |
| The same sine sweep, several times | Several renders | Calibrates a sweep-pair test |
| The same render twice after reloading the plugin | Two renders | Whether its randomness restarts identically |

**To do:** generate the signals (`tools/` script), have the user render them, and add the detectors to the harness (stage E).

---

## 11. Later: plugin implementation notes

- **Core at the native 26,373.6 Hz in integer (16-bit storage, 32-bit accumulate) arithmetic.** It costs a few % of one modern core: the 246 itself ran about 9.6 M instructions/s across two TMS32010s.
- **Host-rate conversion doubles as the analogue emulation:** the two elliptic filters, the S&H one-pole (one IIR per output channel at the device rate), the sub-20 Hz coupling high-pass, and the L half-sample offset (the timing PROM places the two S&H pulses 91–92 steps of 208.3 ns apart).
- **Noise modelled separately:** analogue hiss (73 dB S/N spec), any L-specific behaviour, and truncation noise (comes free from fixed-point emulation).
- **Structure first.** Every free FIR replaced by true structure makes the model both more accurate and cheaper. The hardware budget was about 50–70 weighted taps per sample (256 tap-offset entries at most).
- **Trademark:** "EMT" belongs to EMT International GmbH. Use a non-infringing product name.
