# EMT 246: status and next steps

Hand-over for the next session. Read it with `EMT246_unit_facts.md` (what we know), `library_notes.md` (the Tomkin library and harness results), `model_s_results.md` §8–13 (where the S model stands) and `ml_structure.md` (programs M and L).

*Updated October 2026, after Q2 §15–16 (`model_s_results.md`: fade-free bench; the injection is one fixed, setting-independent excitation; offline it is five equal signed taps per loop, pending the bench) and `period_dsp_practice.md`. Before that: Q1 (Path B across the library, §14) and the M/L structure work (Q3, `ml_structure.md`).*

---

## 1. Where we are

- **Unit topology:** **mono in, stereo out**. One input, two outputs from one multiplexed DAC (L half a sample late). Each stereo IR is the complete 1×2 response.
- **Program S** is a small, deterministic, integer-delay design. Almost all of its structure was identified from one IR (the CLA "Short RT" capture by eugtone):
  - device rate 19.2 MHz / 728 = 26,373.6 Hz (now confirmed by the processor board's timing sequencer);
  - the 8-tap reflection cluster, starting at the DELAY time;
  - 4 integer loops (1321/1389/1461/1587) with gains length-scaled to RT (the Jot law, exact to ±0.1 %), a one-pole HF lowpass per loop (length-scaled, referenced at about 6.4–6.6 kHz) and a first-order LF shelf with one fixed pole (0.9167). The control laws are 13 numbers (`model_s_results.md` §14.2);
  - an all-pass diffuser AP20 · AP72 · AP230 (g ≈ −0.70) plus a 460-sample comb, on the input side. **AP230's loop filter is a two-tap average (1 − b) + b z⁻¹ with b stepped 0 / 0.17 / 0.33 / 0.50 by a table on RT and HF** (short RT with low HF only; §14.3);
  - left/right read-out lags of 26/332 samples;
  - the analogue chain (§4).
- **Nulls** (from one IR; details in `model_s_results.md`):

  | Window | Null | Note |
  |---|---|---|
  | Cluster / first 50 ms | −50 / −52 dB | In-sample |
  | 400–600 ms | **about −28 dB held out** | |
  | R tail after 0.45 s | within 4–6 dB of the noise floor | Loops + 460 only, no fitted input |
  | L tail | 15–22 dB above noise | Fits about 1 % slow gain jitter on L's read-out (unit or capture chain) |

  **Library, held out (Q1, `model_s_results.md` §14.5;** one injection shared by 96 training settings, loops from the laws, clock and gain fitted per capture; 100 Hz–10 kHz, L / R): 0–50 ms −38 / −37 dB; 50–300 ms −29…−32 / −27…−28; 300–600 ms −21…−27 / −21…−25; 600–1000 ms −21…−28 / −21…−26. The repeat-capture ceiling in the same windows is −57 / −46 / −40 dB, so the model is **17–30 dB above it**. Held out equals in-sample (paired captures: +0.0…+0.4 dB). The same injection predicts the CLA unit at −32 / −30 dB (50–300 ms) with nothing digital fitted on CLA.

- **Programs M and L (new, `ml_structure.md`):**
  - **Early section = S's time-scaled by 2 (M) and 4 (L):** same 8-tap cluster and signs (M L 0/108/388/576, R 0/285/388/704; L L 0/216/776/1152, R 0/570/776/1407), reverb onset and L/R gap ×2 / ×4 (**C** for the clusters, **S** for the onsets).
  - **Late reverb is a different structure:** one five-term recursion each, four taps symmetric about N plus one at 2N (L 10826/11384/11618/12176/23002, M 5416/5694/5812/6090/11506; one 2N-word line, **H**). It beats random lag sets by 7–10 dB held out; S's four parallel loops predict S's tail to **−45 dB held out** by the same test.
  - **M and L are time-variant, S is not (S).** AMPL pairs act as repeat captures: S's agree to about −57 dB once each capture's clock ratio (a few ppm) is fitted; M/L tails decorrelate progressively and flat across frequency (L −25 dB at 0.15–0.3 s → −3 dB at 2–4 s). This is the "stochastic control", strong in M/L only.
- **Q2 (October 2026, `model_s_results.md` §15):** fade-free windows and fade compensation; then the main result: **the injection is one fixed excitation per loop and channel, the same for every setting, but about Pᵢ + 1000…1700 samples long, not Q1's Pᵢ + 300.** Trained on 6 flat settings it reaches the repeat-capture null held out (−54 / −54 / −49 dB at 50–150 / 150–300 / 300–600 ms) and predicts all 128 settings at −48.5 / −43 / −39 dB with the plain laws. **What limits S now is the loop-filter laws at strong settings** (LF 0.1, HF 0.2, RT 0.6–1.0): about 1 % per-loop deviations from the length-scaled forms. **Offline (§16, in-sample, pending the bench):** that free injection is **five equal-magnitude signed taps per loop and channel**, with R = L delayed by 26 or 332 per tap, no comb460 on the loop paths, and an empty first-pass FIR. That is 98.7 % of its energy from about 60 discrete parameters plus 2 continuous ones, about 24 memory accesses (the EMT 250 patent's tapped read-out plus output distributor, `period_dsp_practice.md`). The other period candidates (K2–K7) are not supported. **On the bench (§16.6):** the taps hold (equal gains, R = L + 26/332, no comb460), but the tied model reaches only about −32 dB everywhere; a dense, loop-specific feed (1.3 % of the excitation energy, worth about 20 dB of null) is still missing.
- **Main gap (S):** the per-loop injection/read-out structure. Compact fits plateau at about −23 dB where a matching synthetic reaches −46 dB. **Q1 sharpens it:** one free injection shared by every setting reaches only −28…−33 dB (50–300 ms) at long RT and fails at large per-pass damping (RT 0.3 s: −14…−21 dB; HF 0.2: −22…−26; LF 0.1: −24…−28), although the loops themselves follow the laws to < 1 %. The injection needs a branch that passes the loop gain/filter once (x0 + G·x1, **H**: the filter sits partway along each delay line).
- **Ruled out (S):** coupled FDNs; all-passes inside the loops; delay-length modulation (EMT 250 style); randomised injection taps; a separate late-input generator; a long EMT 250-style loop; a shared input or small tap set across loops.
- **New since the single-IR phase:**
  - **The library is in and analysed** (`library_notes.md`): 2,437 settings, S/M/L/Non-Lin, LP and DELAY fixed. The harness run is done; its 116 unknown loop sets are M's and L's recursion seen through the ACF (assigned in `ml_structure.md` §3.1).
  - **The 48k and 96k sets are the same capture** (`library_notes.md` §5), so they are not repeats. **But settings that differ only in AMPL are separate captures of the same reverb** (AMPL changes only the cluster), and they do answer the time-variance question (`ml_structure.md` §4).
  - **The service manual is read** (text pages, processor board, both filters, output stage): `EMT246_unit_facts.md` §2b, page map in `sources/emt-documents/EMT246S-service-manual-page-map.md`.
  - **The analogue chain is modelled** (`tools/analog/`): two filters from the manual plus one S&H pole. This replaces the free band-limit kernel.
  - **Hard hardware limits:** 256 tap-offset entries, 64 K words of delay memory, 182 instructions per sample per DSP. The sparse injection fits (1,719 taps per channel) cannot be the real structure.

---

## 2. Working rules for this phase

- **Work on the 96 kHz set only.** The 48 kHz files are the same capture and give no held-out evidence. Held-out tests must be other settings (or the CLA IR, a different unit), never the other rate.
- **RT comes from the label, not the file** (files are trimmed and faded). The library shows no capture noise: tails reach the 24-bit floor (about −88 dB below the 50–300 ms reverb) unless Tomkin's fade comes first, which it usually does (`ml_structure.md` §4.2). CLA's floor is −71/−77 dB on the same reference.
- **Mind Tomkin's fade.** Every file ends with a linear amplitude ramp of 0.700 s, ending 0.65 ms before the raw file's end. Fit only before it (`analysis/q2/fade2.json`: last non-zero − 0.70 s); report nulls fade-compensated (model × ramp, `tools/q2_injection.py --stage nulls`), never fade the model itself. q1's `fade.json` is wrong for files shorter than about 1.5 s (`model_s_results.md` §15.1–15.2).
- **Smoke-test long runs first** (`--only` a few settings) and check the output before the full run.
- **Fit a clock ratio per capture** when nulling against library files (a few ppm; `ml_structure.md` §4.1). Without it, S repeat captures null only to −34 dB at 0.6–1 s; with it, to about −57 dB. That −57 dB is the most any model can demonstrate against one Tomkin capture.
- **Null target (decided October 2026):** **−50 dB minimum for every program wherever a static IR can support it; otherwise accuracy over null.**

  | Program, part | Target | Why this is the limit |
  |---|---|---|
  | **S**, up to each file's fade | **−50 dB minimum** held out, clock ratio fitted per capture; about −60 dB is "at the measurement ceiling" | Repeat captures agree to about −57 dB with the clock fitted (`ml_structure.md` §4.1); a perfect model against one capture should land about 3 dB deeper |
  | **S against CLA**, R, first ~300 ms | about **−70 dB**, the cross-unit stretch proof (digital part fitted on Tomkin only; only per-unit analogue parameters fitted on CLA, counted) | CLA's noise-limited null is −76 dB there, −54 dB by 400–600 ms; L is noisier |
  | **M and L**, cluster and the gap before the reverb (about 18 ms / 36 ms) | **−50 dB minimum** | The cluster is deterministic |
  | **M and L**, reverb | **Accuracy over null** (below) | Time-variant: two captures agree only to −27…−38 dB at 0.05–0.15 s (M) and −3 to −7 dB after 1 s, so no model can reach −50 dB against one capture |

  **M/L accuracy metrics:** (1) structure: the five-tap recursion, the cluster and onsets, and the mechanism of the stochastic control; (2) per-band decay and spectrum across the library settings; (3) **the capture-to-capture decorrelation curve**: the model, rendered twice with different random states through a simulated sweep, must differ the way Tomkin's AMPL pairs do (−14 → −9 → −7 dB over 0.15–1 s, flat across frequency). Stretch only: if the stochastic control is a deterministic pseudo-random generator, one seed per capture could restore waveform nulls; that also needs the (unknown) sweep simulated.
- **Report nulls** held out, with excess over the noise floor and the parameter count (`CLAUDE.md`). Beware that the library's tails are faded or gated, so "noise floor" there is processing, not the unit.
- **The harness pair classes are unreliable** (`library_notes.md` §3); use `crossrate.csv` for anything cross-rate.
- **Cross-unit check is available:** the CLA IR (S, LP 12, DELAY 10, AMPL 80, LF 1.0, RT 1.0, HF 0.9) is a different unit; 190 library files share its exact S loop set, so the firmware matches. The nearest library settings (S, LF 1.0, RT 1.0, HF 1.0, AMPL 90) are not identical, so use it as a structural check (spec G7), not a null.

---

## 3. Open questions (each needs library IRs)

| # | Question | Test / what is known |
|---|---|---|
| Q1 | ~~**Does Path B generalise?**~~ **Answered (`model_s_results.md` §14): partly.** The loops generalise (laws exact to < 1 %), the first pass does once the AP230 step table is in, and held-out settings null exactly as well as training ones, also across units (CLA). But one shared injection reaches only −28…−33 dB at 50–300 ms and −21…−28 dB later (17–30 dB above the repeat-capture ceiling), worst where the damping per pass is large | Held out by RT step (0.8, 1.8), HF step (0.4) and LF step (0.4); `tools/q1_generalise.py` |
| Q2 | **Per-loop injection/read-out structure** | **Form found (`model_s_results.md` §15.7–15.9):** one fixed, setting-independent excitation per loop and channel, at the measurement ceiling on flat settings, held out. **Structure found offline (§16, in-sample):** 5 equal signed taps per loop (positions in §16.2), R taps = L taps + 26 or 332, no comb460 on the loop paths. **Bench (§16.6), done:** the tied taps are confirmed (equal gains, comb460 removed), but reach only about −32 dB held out and on all 128 settings; a dense loop-specific feed worth about 20 dB is missing (not a common diffuser kernel, not taps inside the diffuser chain). **Found (§16.7):** that feed spans only the taps' range (lags 230–2,230), and is mostly one kernel common to all loops (correlation 0.92–0.94), read through each loop's own taps, with values at the all-pass delays 20 / 72 / 230 and 31 % of its energy before lag 0: most likely the diffuser's all-pass gains slightly off (**H**). **Done (§16.8):** the all-pass gains are **−1/√2**; with them the 9-parameter model reaches −43.5 / −44.4 / −44.1 dB held out (flat) and −42.0 / −38.3 / −35.0 dB on all 128 settings, within 1–6 dB of the free injection. **AP230 steps re-measured with −1/√2 (§16.9):** unchanged within 0.002 (close to 1/6, 1/3, 1/2), and the compact model's score is unchanged. **Checked (§16.10):** a loop refit with the compact injection frozen is contaminated by the injection's error (use the free injection to measure the loops); the tap paths have no fractional delay (+0.0006 samples), no L/R offset and no shared read-path kernel. **Dense remainder with the −1/√2 diffuser (§16.11):** mostly out of band (59 % above 11 kHz); in band small (−29 dB re the taps at 100–400 Hz, −32…−38 elsewhere), not a kernel read through the taps, and not the analogue high-pass (undoing it costs 3.6–3.9 dB). **CLA cross-unit check (§16.12):** the compact model predicts the CLA unit at −36…−37 dB (L / R, 50–600 ms), only 1.4–1.8 (L) / 3.2–4.2 (R) dB short of the free injection, against 8–10 dB on the library: most of the free model's library lead is the library's measurement chain. **R's distributor (§16.13):** binary routing exactly as labelled (other lags ≤ 0.0005), so not the source of R's CLA gap. **CLA ceiling (§16.14):** per-unit corrections (high-pass corners and S&H pole, a 63-tap chain FIR, refitted cluster gains) do not lift it; CLA's cluster-only region is held at about −39 dB by something that is neither a filter nor a gain. On CLA's reverb the compact model is 1.4 (L) / 4.2 (R) dB short of the free injection; whether the free model's 8–10 dB library lead is library-specific stays a hypothesis. **Loop laws (§16.15):** per-loop deviations are negligible; the laws saturate at two limits common to all loops (HF pole ≤ 0.84, per-pass DC gain ≥ 0.18). The HF limit is confirmed (HF 0.2 / RT 0.6 at the null with one constant); the LF limit needs a top-end correction under the cap. **Next:** that LF correction; HF 0.2 at RT 0.8–1.4; comb460 removed outright; the AP230 damping form (first pass only −39…−41 dB where b > 0); the per-loop loop refit with the *free* injection frozen; then (4) per-setting fits with the tied model (the clean K4/K7 check); (5) M/L first recirculation with delays ×2 / ×4. **Then** the loop filters at strong settings (per-loop refit with the injection frozen; per-loop quantised coefficients, **H**) and the short-RT (0.6–1.0) excess. Unit facts §5/§6 (the 460 "fifth recursion", the read-out row) wait for the bench |
| Q3 | ~~**Loop sets and cluster patterns for M and L**~~ **Answered for structure** (`ml_structure.md`) | Clusters and onsets = S's ×2 / ×4. Late reverb: one five-term recursion per program (2N-word line). Open: the integer lags to ±1 (a common −0.8-sample offset, probably a loop filter; §3.3 there), loop gains/filters per setting, the ACF peaks at 2641/5279 and 460/918 (injection or read-out taps?), and the R-only cluster taps being 0.3 % short of ×2 / ×4. The 15 `novel_cluster` files are not yet looked at |
| Q3b | **The M/L time variance ("stochastic control")** | Separate captures decorrelate flat across frequency, progressively (`ml_structure.md` §4). Next: decorrelation rate vs RT/HF/LF and per pass vs per second; which part varies (tap gains, tap positions, read-out) by fitting the five-tap recursion with time-varying parameters; depth and rate bounds. **Without hardware this is the main M/L modelling question** |
| Q4 | **Control laws and step tables:** RT, LF×, HF×, AMPL | **S answered for the library steps (`model_s_results.md` §14.2–14.3):** RT = Jot law exactly; HF and LF laws as below; AMPL scales the cluster linearly and the reverb by ≤ 0.07 dB; a new 4-step AP230 table on RT × HF. Open: the step values between the library's (16 RT, 16 LF, 10 HF steps), the AP230 rule (one candidate in §14.3, **H**), and the same for M/L |
| Q5 | **LF reference frequency** | **S: the LF control sets a first-order shelf (fixed pole 0.9167, corner ≈ 365 Hz)**, so "RT × LF at 50 Hz" holds to ±4 % for LF 0.4 and 2.0 but not for LF 0.1 (0.14 at 50 Hz, 0.13 at DC); the shelf also tilts the top end (RT at Nyquist 1.28 / 1.09 / 1 / 0.96 × RT for LF 0.1 / 0.4 / 1 / 2) |
| Q6 | ~~**HF ×1.0: residual fixed damping?**~~ **S: none** (RT(7 kHz)/RT = 1.000, IQR 1.000–1.002). The HF reference is about 6.4–6.6 kHz, not 7 kHz | |
| Q7 | ~~**Are the loop filters shared or length-scaled? First order, or higher?**~~ **S: first order each.** HF one-pole per loop with a length-scaled pole (a shared pole fits 15 dB worse); LF shelf with one shared pole and length-scaled gains (a length-scaled pole fits 4–14 dB worse) | |
| Q8 | **Is the 460 stage parallel or series? Does its gain vary with any control?** | **Partly:** no control changes it in the first pass; AP230 is the only diffuser stage that steps (§14.3). Series/parallel still open |
| Q9 | **Exact coefficients (Q15 snapping)** for bit-exactness | Joint fits across many IRs; one IR resolves loop gains to about ±8 LSB. **Needed for −50 dB in tails, not only for bit-exactness:** coefficient errors compound every loop pass (±8 LSB at about 50 passes/s ≈ 1.2 % after 1 s, about −38 dB) |
| Q10 | **Fixed-point behaviour** (16-bit storage; truncation or rounding) | **Target: the REVERB program's noise is 15 dB above DELAY's** (73 vs 88 dB re nominal, unweighted rms; nominal = about −6 dBFS). Simulate truncation against rounding in the loops and compare the added noise; the TMS32010 has no rounding mode, so truncation is likely (**H**). Only meaningful once the model nears −60 dB |
| Q11 | **Absolute timing: pre-delay = DELAY exactly?** | Not testable from the library (DELAY is always 0); only the CLA IR has a delay (10 ms) |
| Q12 | **The LP's type and order** (downstream, digital?) | No library data (LP always 12). The through path and the kernel agree at 12 kHz, so the LP at its 12 kHz step adds little. Other libraries (§7) |

---

## 4. The analogue stage (done; use it)

`tools/analog/` simulates the transcribed input filter (7 246 131) and output filter (7 246 130) and compares them with a kernel; `EMT246_unit_facts.md` §2b has the details.

- Chain = input filter (elliptic, 12.2 kHz, with a +5 dB HF boost) × output filter (elliptic, flat to 11 kHz) × **the output S&H lag: a discrete one-pole at the device rate, `y[n] = y[n−1] + (1−p)(x[n] − y[n−1])`, p = 0.325 (CLA) / 0.301 (Tomkin), 0.11 dB rms over 0.5–11 kHz**. p is set by a trimmer (P202), so it is a per-unit parameter.
- Use it as the analogue stage when fitting the digital model, with p free per unit. A skirt shift is optional: it fits **+0.25 % (CLA) and −0.44 % (Tomkin)** and helps 11–12.3 kHz only at a 0.15–0.2 dB cost at 4–11 kHz (`analog_results.md` §2.3). This separates algorithm from analogue and makes the L/R and noise-floor questions cleaner.
- **High-pass (closed, `analog_results.md` §3):** use the unit's coupling corners **3.35, 0.72, 0.23, 0.16 Hz** (fixed; the 3.35 Hz is C106/R117 on the input board) **plus one per-capture-chain term: Tomkin + 1.8 Hz (one corner), CLA + 2 × 6.65 Hz**. Not 15 + 15 Hz: that form is −14 dB (CLA) / −8.6 dB (Tomkin) wrong in the 63 Hz octave, which only a long free FIR hides. Expect the 63 Hz band to improve once a compact model uses the measured form (untested). For the plugin, only the unit's corners belong in the emulation.
- **Input integrate/hold (closed, `analog_results.md` §2):** a leaky integrator, critically damped at τ = 1 µs, in a 9.1 µs window (p. 39 timing plot); −0.05 dB at 12 kHz. Available as `compare_kernel.py --input-sh` / `sample_hold.input_sh`, but not needed (≤ 0.02 dB in any band). It does not explain the skirt.
- Open:
  - the output staircase's zero-order hold is physically present but fits worse than the default (0.18–0.35 dB against 0.11 dB, `--zoh`). So p looks like an effective parameter, and something at 4–11 kHz (pre-emphasis transcription or op-amp model) is 0.2–0.4 dB off (`analog_results.md` §2.5);
  - the per-section skirt shape (needs a longer kernel or a time-domain fit);
  - the unit-to-unit meaning of p (output τ = (R203 + 2·P202)·C205, corrected).
- **When the analogue gaps matter:** not before S's held-out nulls reach about −40 dB. Then they set the band edges: the worst bands now are 63 Hz (−37/−41 dB) and 8 kHz (−42/−40 dB) at 50–300 ms (`model_s_results.md` §9.3). The high-pass and the input window are now closed (above). What is left for the band edges is the ZOH/4–11 kHz question and the skirt. Converter nonlinearity is not needed (THD ≤ 0.03 % = −70 dB).

---

## 5. Suggested order of work

1. ~~**M and L structure (Q3).**~~ Done (`ml_structure.md`).
2. ~~**Decide the null target.**~~ Decided (§2): −50 dB minimum wherever a static IR supports it; accuracy over null for M/L reverb. No hardware work blocks the next step.
3. ~~**Path B generalisation on S (Q1).**~~ Done (`model_s_results.md` §14): the loops, laws and first pass generalise; one shared free injection does not get past −28…−33 dB.
4. **Injection/read-out (Q2).** Form found (§15) and a 5-tap structure found offline (§16). Next, in order (the compact model is §16.8: 5 taps per loop, all-pass gains −1/√2; AP230 steps re-checked in §16.9): then the per-loop loop-filter refit with the tied injection frozen, new laws, and one final shared fit; then update unit facts §5/§6.
5. **M/L time variance (Q3b)**, then the M/L loop gains and filters (Q4–Q7 for M/L) from the deterministic part.
6. ~~**Control laws (Q4–Q8)** from the fits.~~ Done for S at the library's steps (§14.2–14.3); M/L still to do.
7. **Exact coefficients (Q9) and the analogue band edges (§4)** once S's held-out nulls reach about −40 dB; both are needed for −50 dB in every band and in the tails. **Fixed-point noise (Q10)** after that.
   - Alongside Q3b (step 5): read the unread service-manual pages (processor B, the control-CPU side) for how and how often the tap-offset memory is rewritten; that bounds the M/L modulation rate.
8. **Non-Lin (§6)** last.

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

## 8. Modulation / time variance (found in M and L)

- **Hardware makes it possible:** the tap-offset memory can be rewritten by processor A while running, so jittered or modulated read-out taps are feasible (`EMT246_unit_facts.md` §2b). Nothing in the manual links the built-in noise generator (Tests 21/22) to the "stochastic control".
- **Found (`ml_structure.md` §4, S):** AMPL pairs are separate captures of the same reverb. M's and L's tails decorrelate between captures, progressively and flat across frequency (gain or tap-switching type, not delay modulation); S's agree to about −57 dB with a per-capture clock ratio fitted, with no L-only excess. The tail predictor shows the same thing inside single files (a flat, signal-proportional −18 to −22 dB residual in M/L, −45 dB in S).
- **Still open:** depth, rate, which taps vary, and whether the variation is per pass or per second. The captures are sine sweeps, so a time-variant response is seen through an unknown sweep (§9); existence and program dependence don't depend on it, but depth and rate do.
- **Other tools:**
  - AMPL pairs (S: 128 groups of 4, M: 160 of 3, L: 180 pairs), now the main repeat-capture test;
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

## 9. Questions for Tomkin (closed)

Tomkin was asked and does not remember the details: the library was made a long time ago. Treat everything below as **unknown and unrecoverable**, and extrapolate from the data:

- **Sweep method** (type, length, range, level): unknown. The sweep processing gain is therefore unknown, so absolute noise floors and any mapping of time-variance to a rate are out of reach.
- **Raw sweep recordings:** not available (assume gone). No direct time-variance analysis.
- **Trim, fade and `dc` steps:** unknown. The files are cut and faded per setting; work from the labels and the measured structure (RT from the label; fades and zero tails flagged by the harness).
- **Which unit and firmware, and whether one unit made every capture:** unknown. Indirect evidence: all 384 S files share one L/R timing offset (19.19–19.21 µs) and one S loop set, so one unit and one firmware are likely (**H**); the CLA unit differs by about 3 timing steps and a slightly different S&H pole (§4).
- **Recorded twice:** the data say one capture (`library_notes.md` §5). Closed.

Consequence for the working rules (§2): do not rely on tail levels or noise floors in the library as unit noise, expect no repeat captures other than the AMPL pairs (`ml_structure.md` §0), and when a result depends on an unknown above, say so explicitly.

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
