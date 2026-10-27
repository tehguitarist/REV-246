# REV-250 (EMT 250 program) — initial investigation

*Status note (October 2026): a first measurement pass on the Altiverb REV-250 presets, not yet a model. Read with `library_notes.md` §6 (what the files are) and `algorithm_research.md` §1 / `period_dsp_practice.md` §3.1 (the EMT 250 patent and the 1977 algorithm drawing). Every number below is a direct measurement on the deconvolved plugin captures unless marked **[H]** (hypothesis).*

**Scope (confirmed by the user, October 2026):** the REV-250 captures are the **246's Program E — the EMT 250 program running on the 246's hardware** — *not* a real EMT 250. So the target is the 250 algorithm **on the 246 specifically**: device rate 26,373.6 Hz, 64 K-word delay RAM, and the 246's analogue chain. Real-250 captures are welcome later only as a cross-check of the "sounds different from a real 250" claim (`NEXT_STEPS.md` §0); they are not the emulation target. §5 qualifies the shared-chain assumption. (For reference, the real 250 ran at ~24 kHz, not the often-quoted 32 kHz — see §9.)

## 0. Data

- 9 labelled presets `EMT256_REV-250 <RT>s.wav`, RT 0.4 / 1.0 / 1.5 / 2.0 / 2.5 / 3.0 / 3.6 / 4.0 / 4.5 s, plus one unlabelled `REV-250` (a short preset, §6). The RT steps are the EMT 250's own range.
- Deconvolved caches: `analysis/altiverb/cache/REV-250*.npy` (48 kHz, 26 s). All measurements here are on those.

## 1. Decay law

Octave-band T60, L channel (s); the label RT is the 1 kHz value:

| RT label | 125 Hz | 250 | 500 | 1 k | 2 k | 4 k | 8 k |
|---|---|---|---|---|---|---|---|
| 0.4 | 0.49 | 0.53 | 0.54 | 0.45 | 0.47 | 0.44 | 0.41 |
| 1.0 | 1.37 | 1.24 | 1.12 | 1.01 | 0.95 | 0.80 | 0.72 |
| 1.5 | 1.81 | 1.95 | 1.64 | 1.49 | 1.49 | 1.23 | 1.01 |
| 2.0 | 2.66 | 2.64 | 2.27 | 1.96 | 1.85 | 1.63 | 1.29 |
| 2.5 | 3.54 | 3.19 | 2.79 | 2.53 | 2.32 | 2.06 | 1.63 |
| 3.0 | 4.05 | 3.72 | 3.29 | 3.00 | 2.75 | 2.48 | 1.96 |
| 3.6 | 5.03 | 4.36 | 4.07 | 3.55 | 3.38 | 2.97 | 2.30 |
| 4.0 | 5.85 | 5.01 | 4.46 | 3.95 | 3.74 | 3.24 | 2.57 |
| 4.5 | 6.49 | 5.40 | 4.96 | 4.41 | 4.12 | 3.63 | 2.92 |

- **RT(1 kHz) = the label** (0.45 … 4.41 s). **[measured]**
- **LF decays slower, HF faster**, with a fixed per-preset balance: RT(125 Hz)/RT(1 k) ≈ 1.1–1.5, RT(8 k)/RT(1 k) ≈ 0.66–0.91. This is the 250's LF/MF/HF decay set once per preset; they are **not neutral**, the same caveat as the 246 presets (`library_notes.md` §6). **[measured]**

## 2. Onset and stereo

- **R arrives before L:** R first arrival 20.23 ms, L 32.15 ms, so **R leads L by 11.92 ms** (= 314 device samples at 26,373.6 Hz), constant across all nine RT presets. The unlabelled `REV-250` differs (L 19.9 / R 24.9 ms). **[measured]**
- **Early reflections: 3 taps** at 0 / 16.31 / 27.27 ms (= 0 / 430 / 719 device samples) with gains 1.00 / 0.71 / 0.59 (R channel), **identical at every RT step**. L carries the *same* 3 taps delayed 11.92 ms (gains 1.00 / 0.71 / 0.59). This is the 250 patent's 3-tap early-density path (B1–B3). **[measured]**
- **The 3-tap gains are fixed, not RT-dependent.** The patent says their gains are "damped as a function of RT in accordance with an e-function" (`algorithm_research.md` §1.1). Not observed here: either the Altiverb preset freezes the early-reflection level, or the 246's E program uses fixed early taps. **[measured; interpretation H]**
- **The reverb tail is decorrelated between L and R** (ρ ≈ 0 at every lag, including 0 and ±11.92 ms; delaying R by 11.92 ms nulls L to only −2.6 dB). So the 12 ms offset is an **onset-only** L/R delay, and the reverb itself is synthesized per channel — the patent's output discorrelator/distributor (`period_dsp_practice.md` §3.1). **[measured]**

## 3. Loop structure — three combs

Using the harness's decay-flattened ACF method (`tools/harness/core.py: find_loops`, re-implemented at 48 kHz, 200 Hz–8 kHz), which recovers the 246 S loops correctly as a control (50.1 / 52.7 / 55.4 / 60.2 ms):

- **Three loops at 167.98 / 190.12 / 210.52 ms** (48 kHz samples 8063 / 9126 / 10105), present in **both channels**, stable across RT 3.0 / 3.6 / 4.0 / 4.5 s. **[measured]**
- This matches the 1977 algorithm drawing's **"three parallel combs"** (`period_dsp_practice.md` §3.1, P31/P32), not the patent's "about 15 loops".
- **Only resolvable at RT ≥ 3.0 s:** at short RT the feedback gain is too low and the decay-flattening amplifies noise; the ACF at these lags is at the noise floor for RT 0.4–2.0. The lengths themselves are RT-invariant. **[measured]**
- **Integer device samples = 4430 / 5014 / 5552** at 26,373.6 Hz. The measured times sit **+54…58 ppm** above those integers (residuals +0.26 / +0.28 / +0.30 device samples, scaling with the loop length), i.e. a capture-clock ratio, the same order as the known Altiverb clock offset (`library_notes.md` §6, −45 ppm). At 32 kHz the same times give 5375 / 6084 / 6737 — only the middle is near-integer (and 32 kHz is itself now doubtful for a real 250; §9) — so the 26.4 kHz reading is the consistent one, as the confirmed 246 hardware requires. **[measured]**
- **No further loops found** up to 2.5 s by the same method (longer lags are noise). **[measured]**

## 4. Modulation

- **Measurable now, from the single-realisation captures:**
  - **Comb coherence bounds the depth.** Decay-flattened ACF at k × the loop period falls only 0.23 → 0.07 (168 ms loop) from k = 1 to 8, i.e. the three loops stay coherent over ~8 recirculations (1.3 s). That is a little below the damping-only prediction, so a delay modulation is present but **small** — well under one sample, or slow — the same order as 246 M/L's 0.11–0.14 rad rotation (`ml_structure.md` §9). **[measured]**
  - Mode contrast 2.8–3.8 dB vs 5.1 dB for a static 246 S comb (same metric). **[measured; interpretation H]**
  - The sweep-ghost method that gave M/L's rate (`ml_structure.md` §9: ghost lags ±L·f_m/f, halving per octave) needs the loop model first, so `f_m·L` is not readable yet.
- **What one realisation cannot give:** the rate, depth and decorrelation curve. These need repeat captures; the plugin is a convolution of one frozen capture, so re-capturing reproduces the same IR.
- **Data sources, best first:**
  1. **246-E repeat captures at one setting, AMPL varied.** The E program has an AMPL control (remote levers: "pre-delay and ampl for S/M/L/250", `EMT246_unit_facts.md` §3), and AMPL moves only the early section, not the loops. So AMPL pairs are repeat captures of the E's own reverb — exactly the method that found M/L's time variance. **This is the only way to measure the E's modulation without assuming it matches a real 250.**
  2. **Real EMT 250 impulse-train captures** (one impulse every ~10 s, the `NEXT_STEPS.md` §10 method): direct snapshots of the real 250's modulation (depth, rate, which part varies), with no sweep smear. Since the E is "exactly" the 250 program, this is a legitimate **prior** for the E's modulation — but the E's implementation still has to be confirmed against it (the "sounds different from a real 250" caveat).
  3. **Real EMT 250 sweep captures:** ghost-lag analysis as in (1) of the list above, but a sweep smears the modulation, so it is the weakest of the three.
- **A possible shortcut:** if the E reuses the 246's stochastic-control hardware (Tests 21/22) rather than the 250's delay randomisation, M/L's 0.83 Hz / 0.11–0.14 rad may transfer directly. Whether it does is open (§8).

## 5. Chain — shared, but not chain3

- The REV-250 single-tap spectrum matches the **LONG RT 6.0s** preset tap exactly (L and R separately, to ≤ 0.03 dB), and is **darker than chain3 / S / M / CLA**: −2.4 dB at 8 kHz (L), −5.9 dB (R), −10.9 dB at 11 kHz, against chain3's +0.33 / +0.37 dB at 8 kHz. **[measured]**
- So the Altiverb library carries **two tap shapes**: bright (S / M / CLA, = chain3) and dark (L / REV-250). Either the L and REV-250 presets bake in a lower output LP, or they were captured through a different path. The "same hardware chain" assumption therefore holds only within a family: **the reverb algorithm is still separable, but REV-250 needs its own chain kernel** (a "chain4" fitted on the LONG RT cluster, reusing `tools/chain_kernel.py`). **[measured; cause H]**

## 6. The unlabelled `REV-250` preset

- Not in the RT family: T30 ≈ 18 ms (L) / 45 ms (R), first arrival L 19.9 / R 24.9 ms. A short ambience/gate-like preset. Treat separately. **[measured]**

## 7. Discovery plan (reusing the 246 pipeline)

1. **Chain4** for the L / REV-250 family: `tools/chain_kernel.py` on the LONG RT presets' cluster-only windows, so the darker tap is divided out.
2. **Device-rate recovery** (`tools/deconv_grid.py` at 26,373.6 Hz) and a per-capture clock fit (the +54…58 ppm measured in §3) to land the loops on their integer device samples 4430 / 5014 / 5552.
3. **Grey-box fit.** Three uncoupled loops at 168 / 190 / 211 ms (integers once the rate is known); the loop filters carry the LF/HF decay (§1). Loops are **linear in the read-tap coefficients** once the lengths are fixed, so least squares applies — same as 246 S (`model_s_results.md` §17). Priors from the 1977 drawing: **three read taps per comb**, and the output lines use adjacent-tap interpolation (out = f·x[τ] + (1−f)·x[τ+1], `period_dsp_practice.md` §3.1).
4. **Onset model:** 3-tap early path (fixed gains) + the 11.92 ms L/R offset + a per-channel discorrelator; fit after the loops.
5. **Modulation = the accuracy ceiling.** Model the static structure and score accuracy-over-null (as for M/L). Get the modulation from **246-E AMPL repeat captures** (the E's own; §4), and use **real-250 impulse-train captures** as a prior/cross-check. The UAD 250 bench (`NEXT_STEPS.md` §10) validates detectors.
6. **Held-out nulls:** the 9 RT presets are one per RT step, so hold out by RT step (train on 7, predict 2) with per-capture clock/shift/gain, exactly the 246 S protocol.

## 8. Open questions

- The 3-tap early path is fixed in RT here, against the patent's e-function claim — an Altiverb-preset artefact or a real difference in the 246's E program?
- Are there really three read taps per comb (the drawing), and what are their positions/weights? Not yet measured.
- How much of the mode smearing is delay randomisation vs the long-loop mode overlap? §7 step 5.
- **Which modulation mechanism does the E use:** the 250 patent's randomised delay lengths (Fig. 3/4), or the 246's own stochastic control (the 0.83 Hz / 0.11–0.14 rad rotation found in M/L)? If the latter, M/L's numbers transfer. §4.
- **Is the E's modulation the same as a real 250's?** The "sounds different from a real 250" claim needs quantifying; real-250 impulse-train captures answer this directly.
- (Settled) The source is the 246's Program E at 26,373.6 Hz; the loops are 4430 / 5014 / 5552 device samples. A real EMT 250 is not the emulation target, but its modulation is a useful prior.

## 9. Patent and hardware cross-check (review pass, October 2026)

A re-read of the patent PDF and a search for newer sources turned up several details our first summary had compressed, plus one fact that needs correcting.

**9.1 Corrected: the real 250 ran at ~24 kHz, not 32 kHz.** Blesser via Dattorro gives 32 kHz / 8 K words [C: R9], but Morrin's 2024 board-level teardown measures **~24 kHz** (~166 instructions at the 4 MHz clock) and 8,192 × 16-bit words [C: R44]. The patent's "long-loop total delay ≥ 300 ms" decides it: 300 ms = 7,200 words at 24 kHz (fits 8 K) but 9,600 at 32 kHz (does not fit). So the production 250 ran at ~24 kHz; 32 kHz is inconsistent with its own memory. (`algorithm_research.md` §1.1 updated.) Our target is unaffected — the 246's E runs at 26,373.6 Hz — but any future real-250 comparison should assume ~24 kHz.

**9.2 The output stage takes the delayed dry input too [C: R1].** The discorrelator/distributor receives **both** line 12 (delay unit 1, the delayed dry input) **and** line 14 (the reverb), combining/delaying them "such as shown by Schroeder page 225". Our summary mentioned only the reverb feed.

**9.3 The 3-tap early path is per-loop and RT-dependent [C: R1].** Claims 1 and 8 place a second 3-tap delay (t_B, B1–B3, h1–h3) **inside every loop circuit**, fed from the joint input; its gains are "damped as a function of RT". Our measured 3 early taps are **fixed** across RT (§2) — so they are probably **not** the t_B path. Morrin's teardown offers a better candidate: the Reverb program reads **three pre-delay taps** (A/D #1–#3 off the variable +t_d line) [C: R44], which would be fixed. That also explains the 0 / 16.31 / 27.27 ms spacing better than the t_B gains would. Open question (was §8): which path the measured taps are.

**9.4 Randomisation is LFSR address jitter, not a slow sine [C: R44].** The 250's I/O mux has a "Noise" function (an LFSR) and an "8-rotate"; Morrin reads them as generating **pseudo-random memory addresses** for the reverb program. This is the hardware form of the patent's "randomly controlled delay-length variation". It is *not* the 246 M/L's 0.83 Hz rotation — so if the E reuses the 250 scheme, its modulation is fast address jitter, and if it reuses the 246's stochastic control it is a slow rotation. §4's coherence bound (loops coherent over ~8 recirculations) suggests the E's jitter is small either way.

**9.5 Memory check corroborates the 246-E reading.** 8 K words at 24 kHz is 341 ms — enough for one ≥300 ms long loop, but our three loops (168 + 190 + 211 ms = 569 ms) need ~13.7 K words at 24 kHz. They do not fit a real 250; they fit the 246's 64 K words easily. Independent confirmation that REV-250 is the 246's E, not a real 250.

**9.6 Smaller corrections to `algorithm_research.md` §1.1.** Feedback is from "an output of the delay circuit" (drawn as the end), not provably the end; the t_B notation (the PDF OCR showed "t5"); the tapped-loop idea is Kuttruff 1962 [C: R42], with Blesser's contribution the non-harmonic placement + energy e-function + the parallel unreverberated path; the instruction timing is 100 ns (one OCR instance read "100 µs"); claim 5's "300 ns/µs" is garbled — the description's "300 ms" is authoritative.

**Sources checked online:** the patent (Google Patents, full text + claims + citations) [R1]; Morrin's teardown [R44]; DuckDuckGo result set (UA product page — marketing only, no algorithm detail; vintagedigital — 403). No other independent algorithmic analysis of the 250 was found beyond the patent, Dattorro/Blesser and Morrin.

## 10. UAD bench and best-guess specification (25 October 2026)

Detail and methods: `UAD250_test_plan.md` §2 (the plugin is a commercial emulation and a *prior*, not the E; the fits below are to the Altiverb REV-250 captures). Labels: **C**/**S**/**H** as elsewhere.

| Item | Best guess for the 246's Program E | Label | Source |
|---|---|---|---|
| Structure | three parallel combs, no coupling; loop lengths 4430 / 5014 / 5552 device samples | S | §3; plugin has the same design at 24 kHz (4018 / 4560 / 5046) |
| Early taps | 0 / 429 / 718 device samples, gains +1 / −0.709 / +0.59, same in L and R, fixed in RT, LF, HF | S | `UAD250_test_plan.md` §2.9 (fitted integers ±0.1) |
| L/R | R leads L by 314 device samples (11.9 ms) at onset; the tails are decorrelated | S | §2 |
| Decay | T60 = RT × m_b with m_b = 1.40 / 1.26 / 1.11 / 1.00 / 0.92 / 0.82 / 0.67 (125 Hz … 8 kHz), RT ≥ 1.0 s; held-out error 0.5–3 %; RT 0.4 departs | S | §2.9 |
| Modulation type | delay modulation of the combs (not M/L's frequency-flat rotation); early taps and direct path static | H (weak) | plugin; §2.7 comb coherence |
| Modulation rate | slow, ≲ 0.55 Hz (the plugin: 0.552 Hz common to all loops, ≈ 0.74 Hz extra on the longest loop); rates ≥ 0.74 Hz fit worse | H | §2.8 |
| Modulation depth | ≈ ±0.03 ms, ±1 device sample (range 0.01–0.06 ms) on the common line; the plugin's is 10–15× larger | H | §2.8 |
| Not available | read taps and positions per comb, loop-filter shape beyond the per-band multipliers, output discorrelator, the modulation's waveform (sinusoid / drift / random), whether its rate follows the sample clock (0.552 Hz at 24 kHz would be 0.607 Hz at 26,373.6 Hz) | open | needs E repeat captures, which do not exist |

Offline: all UAD captures are archived (`tools/uad250/offline_bundle.py`); the fits here read only the Altiverb caches and synthetic simulation.

## 11. Structure of Program E, decoded (26 October 2026)

*Tool: `tools/e250/e_model.py` (`--stage prep`, `--stage static [--law 246]`); device-rate IRs `analysis/e250/dev_<RT>.npy` (gitignored), each normalised by R's first tap. Supersedes §3's loop lengths, §4's modulation reading and §10's "three parallel combs" / modulation rows. Labels C / S / H as elsewhere.*

**11.1 Each preset is a different realisation; only the direct path is static.** The nine presets' three early taps agree to −53 dB (R) / −61 dB (L) against the mean of the other eight. The dense reverb, from about 20 ms after R's first tap, does not: a leave-one-out mean of the other eight presets nulls it by only about −1 dB in every octave and window from 30 ms to 1.5 s. The decorrelation does not depend on how far apart the RT steps are (4.0 / 4.5 s as poor as 1.0 / 4.5 s), so it is time variance, not RT-dependent structure. An SVD of the nine first passes gives one common component (28 % of the energy) and eight flat, noise-like ones: each capture is **a fixed static IR (≈ 28 % of the reverb energy, ≈ 0.25–0.3 at every time from 30 ms to 1.5 s, more at 125 Hz) plus a fresh random part (≈ 72 %)**. The decorrelation is already complete in the first pass and does not grow with time, unlike M / L's in-loop rotation (captures −27…−38 dB alike early). **[S]** Consequence: for the E's reverb the bar is accuracy over null, as for M / L; the −50 dB null applies to the direct path only.

**11.2 The static network** (both channels; positions in device samples after R's first tap; every number held out by RT step unless stated):
- **Direct path:** taps 0 / 429 / 718, gains +1 / −0.709 / +0.589, the chain kernel (tap 1's own response, 410 samples). L carries the same direct path **314 samples later, inverted, × 1.114**. **[S]**
- **All-pass chain:** three all-passes, **46 / 146 / 392 samples, g = 0.709 / 0.710 / 0.710** (H(z) = (−g + z^−D) / (1 − g z^−D)), outside the loops, shared by the static and the random parts. Found by deconvolution sparsity (L1/L2 0.67 → 0.46), then by least squares on the static first pass, where moving any delay by one sample costs about 2 dB; the same in R and L; deconvolving them also sparsifies each preset's random part (0.77 → 0.69), which the static fit never saw. **[S]**
- **Three loops with negative feedback: 2215 / 2507 / 2776 samples (84.0 / 95.1 / 105.3 ms).** The old 4430 / 5014 / 5552 (§3) are 2 P: a loop with negative feedback has its first positive autocorrelation peak at twice its length. With the chain kernel and the all-passes divided out, the static part is single taps alternating in sign every pass: loop A from 2223 (4438, …), B from 805 (3312, 5819, 8326, 10833), C from 1802 (4578, …). Each loop is read once, its first read scales with the loop gain (the read is after the loop filter: line = G{in − line(t − P)}), gains c = +0.327 / +0.400 / −0.397 (A / B / C, chain-kernel units). **[S]**
- **L's reverb is R's static reverb 315 samples later, same sign** (one sample later than its direct path). **[S]**
- **Loop filter = the 246's own loop law.** `ml_lines.line_coefs` at the E's line lengths with **LF 1.44 / HF 0.72** (one setting for the whole Altiverb family) predicts the octave T60 of the RT ≥ 1.0 presets held out by RT step to 1.6–3.3 % at 500 Hz–8 kHz (6.7 % at 125 Hz, where the measured T60 scatter most), with 2 numbers; a free 7-knot T60 = RT · m(f) curve gets 0.5–1.1 % at 1–8 kHz with 7. In the IR fit the two laws are level (ρ² within 0.01). So the E uses the 246's LF / HF control, plausibly unchanged (**S** for the shape, H for "the same code"). **RT 0.4 departs** (measured T60 8–15 % longer than the label at 500 Hz–2 kHz and up to 34 % at 8 kHz): the 246's HF saturation at short RT would explain 8 kHz (fully saturated −34 → −5 %), not the mid bands; label or T30 bias at a 0.4 s decay, open. **[H]**
- **Count:** 16 numbers (2 direct gains, L gain, 3 loop gains, 3 all-pass gains, 7 law knots; with the 246 law 11), plus the integer delays and the measured chain kernel.

**11.3 The random part: moving reads on the loop lines.** **[S for the observations, H for the mechanism]**
- The random part's energy sits at the same fine-time positions across presets more than by chance (fine-energy correlation of the residuals 0.05–0.27 against ≈ 0.01 for shifted controls) and the strongest static taps are stable across presets (e.g. B's first read −0.10…−0.13). It is not sign or gain jitter on the static taps.
- It recirculates in the loops: its autocorrelation is negative near each loop length, but **split into two lobes about 5–7 samples either side of P** (A −6 / +5, B −7 / +6, C −8 / +6), with a dip at exactly P and nothing at control lags. The reads that make it move by about 6 samples per pass: **≈ ±60 device samples per second, both directions, the same in all three loops** (a delay glide of ≈ ±0.23 %, ≈ ±4 cents). Within one capture this smears each moving read over ~20 samples; between captures the reads are elsewhere, hence the full decorrelation.
- Its energy is spread evenly over the first pass from about 600 samples on (not clustered near the static reads), −7 dB under the static part before B's first read, about equal in the first pass and +2…+3 dB above it after the first recirculation.
- **L's random part is partly R's with the opposite sign** (correlation −0.50 at the static alignment + 1 sample) and partly independent. The output stage adds the moving reads into L and R with opposite polarity, as the patent's output decorrelator / random mixing (§9, Fig. 4) would.
- Hardware reading: Morrin's "Noise" (an LFSR) with the "8-rotate", "to generate pseudo-random memory addresses in the reverb algorithm" (R44), here as read addresses that glide. The motion law (sawtooth, triangle, random walk with reversals), the number of moving reads and their address range are not yet fitted.

**11.4 Held-out scores** (RT ≥ 1.0, each RT step predicted from a fit on the other seven; ρ² = fraction of the held-out capture's 100 Hz–10 kHz energy explained with an optimal gain; the ceiling is the static fraction, ≈ 0.28):

| held-out RT | direct path L / R (dB, to 19 ms) | ρ² 30–165 ms L, R: model / LOO mean | 170–500 ms | 500–1500 ms |
|---|---|---|---|---|
| 1.0 | -52.9 / -49.6 | 0.44, 0.25 / 0.45, 0.26 | 0.21, 0.21 / 0.19, 0.22 | 0.21, 0.22 / 0.09, 0.11 |
| 1.5 | -50.7 / -47.7 | 0.47, 0.25 / 0.42, 0.15 | 0.21, 0.26 / 0.15, 0.16 | 0.23, 0.25 / 0.11, 0.15 |
| 2.0 | -51.5 / -48.5 | 0.48, 0.28 / 0.45, 0.27 | 0.25, 0.25 / 0.22, 0.22 | 0.27, 0.24 / 0.21, 0.23 |
| 2.5 | -49.6 / -46.8 | 0.45, 0.28 / 0.46, 0.30 | 0.26, 0.26 / 0.21, 0.23 | 0.25, 0.27 / 0.19, 0.21 |
| 3.0 | -52.4 / -48.5 | 0.45, 0.28 / 0.37, 0.18 | 0.27, 0.23 / 0.17, 0.13 | 0.28, 0.25 / 0.23, 0.22 |
| 3.6 | -51.7 / -48.4 | 0.47, 0.30 / 0.39, 0.21 | 0.25, 0.26 / 0.16, 0.15 | 0.27, 0.26 / 0.17, 0.16 |
| 4.0 | -51.7 / -47.6 | 0.48, 0.27 / 0.45, 0.25 | 0.27, 0.26 / 0.20, 0.22 | 0.27, 0.26 / 0.17, 0.18 |
| 4.5 | -50.8 / -47.6 | 0.46, 0.26 / 0.43, 0.25 | 0.27, 0.27 / 0.20, 0.19 | 0.28, 0.25 / 0.20, 0.19 |
| **median** | -51.6 / -48.0 | 0.46, 0.27 / 0.44, 0.25 | 0.26, 0.26 / 0.19, 0.20 | 0.27, 0.25 / 0.18, 0.18 |

The direct path is at the capture's own floor there (the leave-one-out mean of the other presets reaches −50.4 / −49.2 dB in the same window at RT 2.5): random reads already start around tap 2. The reverb's static part is explained to its ceiling in the first pass and better than the leave-one-out mean afterwards, with 16 numbers.

**11.5 Literature check (26 October 2026).** No public source gives the 250's delays, gains or LFSR (Morrin's pages and the patent are structural only; nothing documents the 246's Program E). An unsourced forum recollection (Eventide forum, 2008) describes "predelay, 3 all-pass diffusers, a tank, a 4-output decorrelator" with tank times of 80–120 ms; the decoded structure fits it (three all-passes; loops of 84–105 ms).

## 12. The moving reads: a generator fitted to the random part's statistics (26 October 2026)

*Tools: `tools/e250/e_sim.py` (Program E in the time domain at the device rate, integer-address moving reads, simulated sweep captures deconvolved like the library's, statistics), `tools/e250/e_genfit.py` (held-out fit). Data statistics: `analysis/e250/rand_stats_data.json`; fits: `analysis/e250/genfit_*.json`.*

**12.1 Method.** The static network is §11's (the 246 law, all-data fit, `static_all_246.json`). The random part cannot be nulled (one realisation per preset), so the generator is fitted to statistics of r = capture − static model, computed identically on the eight RT ≥ 1.0 presets and on simulated sweep captures at the same RTs: random / static energy in four windows (before the first loop read, first pass, 2215–4400, 4400–9000 samples), r's autocorrelation at each loop length ± 12 samples ("lobes"), the L / R correlation of r at the static alignment ± 2, and the static fraction per octave (250 Hz–4 kHz). Held out: fitted on the statistics of RT 1.0 / 2.0 / 3.0 / 4.0, scored on 1.5 / 2.5 / 3.6 / 4.5, and the reverse. The reference is the data against itself: the two halves' statistics differ by a cost of 41 (energy 1.6, lobes 36, L / R 1.0, octaves 2.1).

**12.2 Result. [H, statistically constrained]**
- **Six moving reads per loop**, each at an integer output delay o(t) **between ≈ 405 and ≈ 3600 samples** after R's tap 1, gain **± 0.81 × the loop's static read** (random sign), on the loop lines (they recirculate with the loops).
- **They move at a peak speed of ≈ 75–90 samples/s** (≈ 0.3 %, ≈ 5 cents of glide). The sweep captures fix the speed but not the turning rate: a sinusoid at 0.3 / 0.55 / 0.8 / 1.2 Hz scores the same at equal peak speed (f · A ≈ 14). With the UAD plugin's rate (0.55 Hz, `UAD250_test_plan.md` §2) the depth is **± 25 samples (± 0.95 ms)**, about 2.5 × the plugin's. An asymmetric triangle (up 60, down 79 samples/s, speeds ± 30 %) is level with the sinusoid (held-out cost 146 / 212 against 131 / 148); the plugin's sinusoid is the better-founded choice and is the recommended form.
- **L:** half the reads are shared with R with the **opposite sign** (share 0.50), the rest independent; L's moving reads sit 314 samples later (the direct path's delay), the static reverb 315.
- **Scores (sine, 0.55 Hz, ± 25):** energy ratios at the data's own level (cost 1.7 vs 1.6); lobes 2 × the data's (71 vs 36: the simulated lobes are a little sharper); static fraction per octave 0.20 vs 0.22–0.28; held-out total 131–148 vs 106 in sample vs 41 data-to-data.
- **Robustness:** the unknown Altiverb sweep length does not bias it (held-out cost 119 / 131 / 100 for 5 / 10 / 20 s sweeps). The read count is **not** identifiable: 3, 4 and 6 reads per loop at equal total energy (M b² = 6 × 0.81²) score 127 / 122 / 127 averaged over three seeds (seed scatter ± 25); only the energy is fixed. Any M from 3 to 6 is consistent; more reads make the sound smoother, fewer grainier, which is a listening decision for the plugin.
- **Not matched:** the L / R lag pattern (data −0.22 / −0.44 / −0.14 / +0.17 at lags −2…+1; an opposite-sign copy gives −0.13 / −0.39 / −0.32 / −0.04 at any delay). The odd shape suggests L's moving reads pass an L-only filter (the patent's all-pass output decorrelator, §9), which the static reverb does not. Open.

**12.3 Cross-check with the UAD 250 plugin.** Its loops (4018 / 4560 / 5046 samples at 24 kHz by the same autocorrelation method, i.e. 2 P) are 2208 / 2506 / 2773 at 26,373.6 Hz when halved: **the E's 2215 / 2507 / 2776 to 0.1–0.3 %**, so the plugin has the same three negative-feedback loops in milliseconds, and its early taps are static too. Its "loop delay modulation" (0.552 Hz, ± 0.4 ms) was read from sidebands, which cannot tell a modulated loop from moving reads on the loop lines; in the E the loops themselves are static (§11.2: the static taps stay sharp to the fifth pass), so the moving reads are the likely mechanism in both.

**12.4 Program E for the plugin (best guess).** Device rate 26,373.6 Hz, integer delays, no interpolation:
1. direct path: x + (−0.709) x[−429] + 0.589 x[−718]; L: −1.114 × the same, 314 later;
2. u = AP46 · AP146 · AP392 (x), g 0.709 / 0.710 / 0.710;
3. loops k = A / B / C, P = 2215 / 2507 / 2776: line_k[t] = G_k{u[t] − line_k[t − P_k]}, G_k = `ml_lines.line_coefs` / `line_resp` at (RT, LF, HF) (the 246's own law; the Altiverb family sits at LF 1.44 / HF 0.72);
4. static reads: c_k · line_k[t − d_k], d = 2223 / 805 / 1802, c = 0.327 / 0.400 / −0.397 (relative to the direct tap, chain-kernel units);
5. moving reads (**superseded by §13.4:** a triangle LFO at 0.553 Hz, ± 32 samples, ± 30 % speed spread): M per loop (3–6; data fix only M b², 6 × 0.81² ≈ 3.9), gain ± b |c_k|, ~~output delay o(t) = centre + 25 sin(2π 0.55 t + φ)~~, centre uniform in [410, 3600], independent phases; half shared into L with the opposite sign, half independent; L's static reverb 315 later, its moving reads 314 later, × 1.114;
6. the 246's analogue output chain as for S / M / L.
About 3 + 3 + 3 + 3 + 18 multiply-adds per sample per channel plus the loop filters: inside the hardware's budget. **Open:** the L-only filtering of the moving reads; RT 0.4 (§11.2); AMPL, pre-delay and LP for E (one setting in the Altiverb family); whether the turning rate is the plugin's 0.55 Hz (needs E repeat captures at one setting, which do not exist).

## 13. The UAD 250 plugin as E's reference, and the closing checks (26 October 2026)

*Data: the UAD 250 bench captures (`UAD250_test_plan.md`; gitignored, used read-only). Tools: `tools/e250/e_sim.py` (`--stage render`, the reference renderer; `lfo` law), `tools/e250/e_genfit.py` (`--law lfo_tri | lfo_sine`). Supersedes §12.2–12.4 where they differ.*

**13.1 The plugin runs the same algorithm.** Its 960 impulse snapshots at Time 0.4 (C_dense, every 0.25 s for 240 s), brought to its 24 kHz grid and averaged, give its static part; the §11 search on it finds **early taps 0 / 389 / 653 (+1 / −0.711 / +0.587), L 285 later and inverted; all-passes 41 / 132 / 356, g 0.70 / 0.69 / 0.67; loops with negative feedback P = 2009 / 2280 / 2523 (the bench's 4018 / 4560 / 5046 are 2 P), first reads 2008 / 724 / 1631**, each recirculating with alternating sign. In milliseconds every element is E's within 0.1–0.4 ms (E's delays are about one sample longer than an exact rescale to 26,373.6 Hz). **[S]** So the 246's E is the 250 design re-timed to the 246's clock ("exactly the EMT 250 program").

**13.2 Its moving reads, measured directly.** Snapshot minus the mean is one realisation of the moving reads, free of sweep smear. Between snapshots 0.25 s apart the moving pattern shifts by **± 18 samples (24 kHz), bimodal, rarely near 0: a constant speed of ≈ 72 samples/s (≈ 79 at 26,373.6 Hz), reversing every half period of 1.81 s (0.553 Hz)**: a **triangle** LFO, ≈ ± 32 samples. From 1200 to 3200 samples every window moves at the same rate and phase (folded velocity ± 15 samples per 0.25 s, flat-topped); the earliest reads (400–1200) follow another phase. This is the bench's 0.552 Hz line (§2.4–2.5 there); its "loop delay modulation" is these reads. **[S for the plugin]**

**13.3 The plugin and E have the same recirculation signature.** The leave-one-out residual's autocorrelation at each loop length (the §12 "lobes") on the plugin's impulse trains (Time 4.5, 24 impulses; Time 2.5, 160 snapshots) has **lobes at ± 5…6 samples, the + side stronger, and an empty centre** (+0.003 at exactly P); the Altiverb E presets, by the same leave-one-out pipeline: lobes at −6 / +5…6, + side stronger, centre −0.007. A simulated triangle LFO gives the empty centre, a sinusoid fills it (−0.02…−0.05). (A first run of the plugin's own sweep captures through a quick deconvolution showed a sharp centre instead; that pipeline was not clean enough for this statistic and the reading is withdrawn.) **[S]**

**13.4 A bias in §12's statistics, fixed.** §12 took E's random part as capture − static model; any static-model error recirculates at exactly P and fills the lobes' centre in the data but not in a simulation, which pushed §12 towards the sinusoid. Lobes and L / R now use the leave-one-out residual for data and simulation alike (`e_sim.stats`). The data's own halves then differ by a cost of 100 (the leave-one-out statistics are noisier). Held out, the plugin's triangle law (A, speed spread, gain, range, share fitted) scores 170 / 206 and a sinusoid 175 / 131: E's statistics cannot tell the two apart; the plugin can. **Adopted for E** (`e_sim.ADOPTED`): **triangle LFO, 0.553 Hz, ± 32 samples (≈ 71 samples/s), per-read speed spread ± 30 %, independent phases; 6 reads per loop, gain ± 0.81 × the static read, centres 410–3600 samples after R's tap 1; half shared into L with the opposite sign (L's moving reads 314 later).** The rate, the waveform and the phase grouping come from the plugin (H for E); the speed, energy, range and L / R share are E's own fit. Shared phases (groups of reads) fit E's statistics no better (cost 161 vs 147); the plugin's grouping is not resolved, independent phases are adopted.

**13.5 RT 0.4 is estimator bias.** The model built exactly to the law at RT 0.4, captured and measured with the same T60 estimator, reads 1.31–1.55 × its label, more than the data's 1.02–1.35 ×: the "RT 0.4 departs" of §11.2 is the T30 estimator at a 0.4 s decay, not the unit. No HF saturation is needed.

**13.6 The decay law's high-band calibration depends on the original capture.** Moving reads lengthen the measured high-band T60 of a sweep capture (the deconvolution of a time-variant system smears energy late): simulated at RT 1.0, 8 kHz reads 0.75 / 0.77 / 0.83 / 1.02 × RT for 2 / 5 / 10 / 20 s sweeps and 0.76 for an impulse. The Altiverb data read 0.72, as the static law itself (0.70; 1 kHz 1.00 vs 1.01, 4 kHz 0.83 vs 0.80): Audio Ease's capture smeared little (a short sweep or an impulse-type method). The law (LF 1.38 / HF 0.74, §11) stands; its high bands are uncertain by about ± 10 % with the capture method unknown. **[H for the HF calibration]**

**13.7 The capture clock: 42.8 ppm (checked).** The 18 static loop taps out to lag 15,682 (passes 1–5, nine-preset mean, kernel and all-passes divided out) drift linearly by 0.19 samples (0.012 rms about the line), which first read as +12.3 ppm (55.1 in all; §3's +54…58 the same way). Re-prepared at 55.1 ppm, the direct path gets 4–5 dB *worse* held out (−51.6 / −48.0 → −46.8 / −45.8 dB) while the reverb is unchanged: the drift is the loop filter's group delay (a few hundredths of a sample per pass), not the clock. Scanned by the model's direct-path null (three folds), the clock peaks at **42.8 ppm** (24 / 33 / 42.8 / 52 ppm: L −43.3 / −47.5 / −50.7 / −47.4, R −44.0 / −46.3 / −47.6 / −46.2 dB), the S presets' value. Kept; §11.4's scores stand. Lesson: time the clock on the direct path, never on recirculations through a loop filter.

**13.8 Status of Program E.** Everything the Altiverb presets and the UAD bench can constrain is fitted and held out: the direct path (a null at the capture's floor), the static network (11 numbers on the 246 law), the decay law, the moving reads' energy, range, speed and L / R share (E's own statistics) and their motion law (the plugin's). **Not resolvable here:** the turning rate and phase grouping *in E* (needs E repeat captures; the plugin's are adopted), the L-only filter on L's moving reads (§12.2), AMPL / pre-delay / LP for E (one setting in the family; the 246's conventions apply, H), the decay law's high bands to better than ± 10 %. Reference renderer: `e_sim.py --stage render` (checked: L 11.92 ms after R, taps at 16.29 / 27.25 ms, octave T60 at RT 2.5 within 1–9 % of the preset from 250 Hz to 8 kHz, 125 Hz 16 % long, the law's known low-band scatter).

## 14. Last checks: L's moving reads, the read inventory, the dark chain, the unlabelled preset (26 October 2026)

*Tools: `tools/e250/uad_moving.py` (plugin snapshots; outputs plugin-derived, kept local), `tools/e250/lp_chain.py`, `tools/e250/rev250_unlabelled.py` (Altiverb; `analysis/e250/lp_chain_*.json`, `rev250_unlabelled.json`), `e_sim.py` (`law="plugin"`). Run by two sub-agents, checked here.*

**14.1 The dark chain is the L program's output kernel, not the LP lever. [S]** Per channel, the dark presets' first-tap response (REV-250 and LONG RT, identical to 0.02 dB / 0.2°, L = R to 0.05 dB) divided by the bright presets' (S, MED) is not a clean causal low-pass (best Butterworth / Bessel / Chebyshev / elliptic of orders 1–6: 0.86 dB / 2.3° rms) but a 2-tap device-rate FIR [0.77, 0.23] (one number, 0.068 dB / 0.50°, held out 0.075 dB / 0.59°; free tap spacing 36.5 µs against the device period 37.9 µs; 44.1 / 48 kHz spacings fail). Its tap a = 0.230 for E / S and LONG / S, 0.047 for MED / S, 0.001 for E / LONG: the M / L output kernels already measured on Tomkin's real unit (`ml_structure.md` §7.5, p ≈ 0.05 M / 0.23 L; magnitudes within 0.1 dB to 12 kHz). So **E has L's output kernel**; the LP lever is at its open setting in every Altiverb family and its response cannot be measured here. Correction to §5: its −2.4 / −5.9 dB at 8 kHz were absolute levels; the dark-over-bright excess is the same in L and R. (`ml_structure.md` §7.5's reading of that kernel as averaged ±1-sample address jitter would fit E's strong time variance; H.)

**14.2 The unlabelled `REV-250` is not E. [C]** Two bare taps, L 20.04 ms and R 24.99 ms (528 / 659 device samples), equal gains, same polarity, S-type chain (kernel a = 0.000); no E taps (< −100 dB at +429 / +718), no reverb (tail ≤ −71 dB from 50 ms, ≤ −92 dB later; L / R tail coherence 0.99), no E network (LS gain on E's static template ≤ −60 dB; the all-passes do not sparsify it). Probably the D (delay) program at feedback 0 (H). It says nothing about E's pre-delay or L / R controls.

**14.3 L's moving reads in the plugin: no L-only filter. [S]** Every strong moving read appears in L at the same address and gain as in R (offset 0.00 ± 0.15 samples, L / R gain 1.00 ± 0.02), 285 samples later, with a fixed sign per read (+R for one set, −R for the others); the static reverb is +R, the direct path −R. A per-address sign map explains 0.39 of L's moving energy held out (first / last 120 s). The plugin's L / R lag pattern is symmetric; **E's asymmetric pattern (−0.22 / −0.44 / −0.14 / +0.17) stays unexplained** (H: some of E's L reads 314 and others 315 samples later, as its direct path and static reverb are).

**14.4 The plugin's read inventory. [S for the plugin]** Seven strong moving reads (single sharp taps), centres 578 / 1091 / 1285 / 1327 / 1710 / 1891 / 2176 at 24 kHz (first / last 120 s agree to 0.2 sample), equal gains (0.17 × the direct tap, 0.42 × B's / C's static read), on loops A / C / B / A / A / B / C (recirculation at c + m P_k with alternating sign); about 45 % of the moving energy; the rest is a population of weaker reads, not resolved. **All reads share one displacement g(t)**: a triangle, ± 32 samples, 71 samples/s, ≈ 0.5535 Hz with a phase that wanders by ≈ 0.1 cycle per 30 s (the "other phase" of §13.2 was intermittency). **The strong reads switch on and off in sets** ({1285, 1710, 2176}, one per loop, present ≈ 62 % of the time; {1091, 1327} ≈ 29 %; {578, 1891} ≈ 9 %; within-set presence correlation 0.3–0.9, confirmed independently here; the sets' mutual exclusion is H), with dwell between ~0.05 and 0.25 s.

**14.5 Tested on E: E prefers more diverse motion.** The plugin-structured law (the seven reads rescaled × 1.0989 with their signs and sets on one common triangle, plus a weak population on the same triangle, its gain fitted on E) reaches E's energy ratios but scores 250 (seed scatter ± 40) against 182 for §13.4's adopted law (independent phases, ± 30 % speed spread) and 100 for the data against itself: lockstep motion makes the recirculation lobes too sharp for E. §13.4's law stays the E generator; `law="plugin"` keeps the plugin-faithful structure as a listening alternative. Whether this is a real 246-versus-250 difference or the Altiverb capture smoothing the lobes cannot be told from one realisation per preset. **[H]**

**14.6 Program E: final status.** Measured and held out: the direct path, the static network, the decay law, the output kernel (L's), the moving reads' energy, range, speed and L / R share. From the plugin (same algorithm, H for E): the triangle motion at ≈ 0.55 Hz, the read inventory and sets. Educated guesses, documented: E's LFO rate (0.553 Hz, the plugin's; 0.608 Hz if the 250 counted it in samples), AMPL (the 246's S / M / L convention), pre-delay range (the 250's 0–315 ms, "extended" on the 246), the LP lever (shared with S / M / L, not measurable from any existing capture), the high bands of the decay law (± 10 % with the original capture method unknown), the L-only detail of E's moving reads (§14.3).

## 15. Program E: the plugin specification (26 October 2026)

**Single source:** `analysis/e250/e_plugin_spec.json`, generated by `tools/e250/e_spec.py` from the reference model's own constants (every number with its label and section). **Reference model:** `tools/e250/e_sim.py` (`render_ir`, `--stage render`), the time-domain structure a plugin mirrors; it agrees with the independent frequency-domain model (`e_model.render`) to −150 dB (static, RT 0.4–2.5; −117 dB at 4.5 where the 4 s render truncates).

**Signal flow at 26,373.6 Hz, integer delays, no interpolation:** input → direct path (taps 0 / 429 / 718, +1 / −0.709 / +0.589, × AMPL / 80) and → three all-passes (46 / 146 / 392, g 0.709 / 0.710 / 0.710) → three loops with negative feedback (P 2215 / 2507 / 2776; line = G{u − line[t − P]}, G the 246's loop law at RT / LF / HF) → static reads (2223 / 805 / 1802, gains 0.327 / 0.400 / −0.397) + moving reads → R; L = −direct (314 later) + static reverb (315 later) + moving reads (314 later; half shared with opposite sign) × 1.114 → L's output kernel [0.77, 0.23] → DELAY (whole output) → LOW PASS → the 246's analogue chain.

**Controls and the educated guesses (H), as implemented in `render_ir`:**
| Control | Implementation | Basis |
|---|---|---|
| RT / LF / HF | the 246's loop law; the Altiverb family at LF 1.38 / HF 0.74 | S, §11.2 / §13.6 (HF bands ± 10 %) |
| Moving reads (the 250-specific modulation) | **superseded by §16:** E's nine fixed reads (`fixed_reads.json`), each gliding ± 32 samples on a triangle at 0.553 Hz (`lfo="time"`; `"samples"` = 0.608 Hz) × (1 ± 0.5) per read, independent phases; L the same addresses with signed amplitudes; §13.4's random reads as `law="random"` | reads, amplitudes, depth: E's low band (S); speed spread: E's statistics (S); triangle and rate: the UAD plugin (H) |
| Alternative motion | `law="plugin"`: the plugin's 7 reads in switching sets on one triangle + weak reads | §14.4–14.5, fits E worse; for listening |
| AMPL | direct taps × AMPL / 80 (the reverb unchanged) | S's linear cluster law; the family taken as AMPL 80 (H) |
| DELAY | a pure delay of the whole output, 0–999 ms | operation card (H for E) |
| LOW PASS | open at 12 kHz; below, a 2nd-order Butterworth at the device rate | never measured (H) |
| Output kernel | [0.77, 0.23] at the device rate | S, §14.1 |

**How to check a plugin's E:** (1) the static network (moving reads off): null `analysis/e250/ref/e_static_rt{0.4,1.0,2.5,4.5}.wav` (device rate, (L, R), R's tap at index 64, AMPL 80, DELAY 0, LP open, family LF / HF) sample for sample; (2) the moving part: render at the eight RT ≥ 1.0 steps, take the statistics with `e_sim.stats` (leave-one-out residuals) and compare with the spec's `targets` (the Altiverb presets'); the adopted generator lands at about 1.3–2 × the data's own half-vs-half cost; (3) octave decay against §1's table (the law: 1.6–3.3 % at 0.5–8 kHz). **Not resolved, carried as guesses:** see the spec's `not_resolved`.

## 16. The full band, E's fixed moving reads, the kernel tail and LF (26 October 2026)

*Prompted by a full-band check (the earlier scores were 100 Hz–10 kHz). Tools: `e_model.py` (`render_mean`, `kernels` with `KLONG` and the tail, `--stage tail`, `E_LF`), `e_sim.py` (`law="fixed"`, now `ADOPTED`), `analysis/e250/fixed_reads.json`, `kernel_tail.json`, `static_cv3.json`. Supersedes §12.4 / §13.4 / §15's moving-read law and §11's LF.*

**16.1 What the full band showed.** Direct path 20–100 Hz −23 dB against a −38 dB floor (the leave-one-out mean); reverb below 200 Hz 92–94 % repeatable between captures while the model explained 17–45 %; 63 Hz decaying 15–20 % too fast. 10–12.4 kHz was fine (decay 0.62 against 0.56–0.64 × RT; direct path −48 dB). Above 12.4 kHz the captures hold −46…−48 dB of their energy (the analogue band limit), left to the plugin's chain. (S / M / L had this check long ago: `ml_structure.md` §11, `model_s_results.md` §23, §25, §27.)

**16.2 E's moving reads are fixed reads with a moving offset. [S]** A ±35-sample glide barely moves the phase below 200 Hz, so a read with a fixed centre and sign is static there. A greedy sparse fit of "a read on loop k at centre c, averaged over its glide (a centred boxcar)" on the presets' mean minus the static model finds **nine reads: 641 A, 1205 C, 1415 B, 1462 A, 1883 A, 2086 B, 2392 C, 2550 B, 2814 C** (output delay after R's direct tap, device samples), amplitudes 0.30–0.42 (direct tap = 1). Seven are the UAD plugin's seven (§14.4: 635 A, 1199 C, 1412 B, 1458 A, 1879 A, 2078 B, 2391 C at E's clock) on the same loops, within 1–8 samples; signs all inverted (the plugin's known polarity, §2.6 of the bench). **Held out** (reads from seven presets, the eighth predicted), R, ρ² 31.5 / 63 / 125 / 250 Hz: static only 0.54 / 0.36 / 0.49 / 0.38 → **0.96 / 0.94 / 0.84 / 0.54**, against the leave-one-out mean's 0.96 / 0.92 / 0.79 / 0.48; 9 reads best held out (7: 0.90 / 0.88; 12 or 16 no better). The best boxcar is 70 samples (50 / 100 worse at 250 Hz): **E's own glide depth is ± 35 samples**, the plugin's ± 32 × 1.0989. **L: every read at the same address with |L / R| = 0.83–1.10 and a fixed sign** (no filter): {1415 B, 1883 A, 2392 C} same sign as R, the other six inverted, the plugin's sets (§14.3). Unlike the plugin, E's reads look always present (similar amplitudes across the plugin's sets). (A first print of the centres was 64 samples short, a tap-offset slip in the fit script; corrected, and the "60 samples earlier" note it caused is withdrawn.)

**16.3 The generator.** `law="fixed"`: the nine reads, each at centre + A tri(f_j t + φ_j), independent phases, f_j = 0.553 Hz × (1 ± fsd). On E's statistics (leave-one-out, held out by RT, three seeds; data against itself 100): **A 32, fsd 0.5: 122.8**; A 35 / fsd 0.3: 161.5; common phase 301; fsd 0: 249; §13.4's random reads 187. The moving / static energy ratios follow from the low-band amplitudes with nothing fitted on energy (−16.0 / +4.0 / +5.5 / +4.7 dB against the data's −15.5 / +4.3 / +5.5 / +4.8). **Adopted** (`e_sim.ADOPTED`). In E the reads move independently (in the plugin together, §14.4). The averaged time-domain model matches the analytic expected response (`render_mean`) to −20.6 dB below 200 Hz over 12 realisations.

**16.4 The chain kernel's length and slow tail. [S]** The kernel stopped at 410 samples (tap 2 at 429). Removing tap 2's own copy recursively, k[n] = d[n] + 0.709 k[n − 429], extends it to 540 (the first reads start at ~606), its tail −61…−77 dB (the first version started the recursion 16 samples late and left tap 2's pre-ringing at −25 dB; fixed). Beyond it a cubic B-spline tail (540–4000 samples, 16 numbers per channel, fitted below 250 Hz on E's expected response) holds the capture chain's slow low-frequency tail (as chain3's B-spline tail does for the library, `ml_structure.md` §13). Held out by RT: direct path 20–100 Hz L −23.9 → −33.7, R −25.6 → −29.5 dB; full band L −47.1 → −52.4, R −45.1 → −47.3 dB. A per-capture DC term, earlier or longer splines did not help R. This tail is the Altiverb capture chain's, not E's; the plugin's analogue chain replaces it.

**16.5 LF from the decay statistics: 1.49. [S]** With the reads and the tail in, the expected response's octave T60 against the data (63 / 125 / 250 / 500 Hz, both channels, the same estimator) puts LF at 1.49 in both RT splits; held out the T60 error is 0.7–3.6 % (at the IR-fit LF 1.383: 5.1–6.5 % at 63 / 125 Hz). The IR least-squares LF was biased by the random part (as §13.5–13.6 found for RT and HF). Adopted (`static_all_246.json`); reads and tail refitted on it (unchanged to a few %).

**16.6 Held-out scores now** (246 law, LF 1.49, the other numbers fitted per fold on seven RT steps; reads and tail from all eight — their own held-out tests are 16.2 / 16.4), medians over the eight RT ≥ 1.0 folds:
- **Direct path (to 19 ms): 100 Hz–10 kHz L −53.3 / R −50.9 dB** (floor −54.1 / −53.5); **20 Hz–12.4 kHz L −52.5 / R −47.4 dB** (floor −54.6 / −53.1).
- **Reverb, R, explained energy, model / leave-one-out mean: 30–165 ms 0.28 / 0.25, 170–500 ms 0.29 / 0.20, 500–1500 ms 0.30 / 0.18**; per octave (170–1500 ms) at or above the data's repeatability in every octave from 31.5 Hz to 11 kHz.
- Decay: 0.7–3.6 % at 63–500 Hz (16.5), 1.6–3.3 % at 0.5–8 kHz (§11.2), 10–12.4 kHz within the data's spread.

**16.7 The two findings left by 16.6, closed. [S]**
- **L / R of the moving part.** The data's L / R correlation (−0.22 / −0.44 / −0.13 / +0.19 at lags −2…+1) needs the **same-sign reads (1415 B, 1883 A, 2392 C) to reach L one sample later than the inverted ones: 315 (the static reverb's delay) against 314 (the direct path's)**. With that (`dl_same=1`, adopted) the generator gives −0.13 / −0.44 / −0.10 / +0.24 and **a held-out statistics cost of 86.3, below the data's own half-against-half level (100)** (all at 314: 122.8; +2 / −1 variants 114–217). So L's output sums the inverted reads with the direct path and the same-sign reads with the static reverb.
- **R's direct path below 100 Hz** was the captures' **non-causal low-band content before the tap** (low-passed, −60…0 samples: model error −62 dB against the leave-one-out mean's −87; after the tap within ~3 dB of the data's level): deterministic across presets, the capture's deconvolution, not the unit's causal chain. As a capture-chain term (`e_model.pretap`, the other presets' mean minus the expected response, low-passed below 150 Hz, kept before the kernel starts; for comparisons with captures only, never part of E or the plugin), held out by RT, the **direct path (to 19 ms) is L / R −32.9 / −32.3 dB at 20–100 Hz, −54.0 / −53.1 at 100 Hz–10 kHz, −49.1 / −50.3 at 10–12.4 kHz and −52.8 / −51.6 dB at 20 Hz–12.4 kHz** (floors −38.0 / −37.4, −54.1 / −53.5, −55.1 / −54.8, −54.6 / −53.1); R full band per fold −48.3…−52.6 dB (6 of 8 at −50 or deeper).

**16.8 Still open (not resolvable from existing data).** The LFO rate (0.553 / 0.608 Hz) and the triangle shape (from the UAD plugin, H); the decay law's high bands ± 10 % (Audio Ease's capture method); AMPL / DELAY / LP for E; the direct path's remaining 5–6 dB excess at 20–100 Hz and 10–12.4 kHz (capture chain; small).

## 17. The reverb-metric audit (G1–G6) and its fixes (27 October 2026)

*Tools: `tools/e250/e_pulse.py` (written by a sub-agent: pulse_check's metric functions on Program E — the model's realisations (5 seeds, impulse renders and 10 s sweep captures) against each Altiverb preset, inside if within 2.5 × the realisations' spread or an estimator floor; both channels, extended bands 20–100 Hz and 10–12.4 kHz), summaries `analysis/e250/e_pulse*.json` (v1 before the fixes, v3 after). Supersedes §16.3's amplitudes and §16.5 / §11.2's HF.*

**17.1 What the first audit found.** Inside the spread: decay shapes 125 Hz–2 kHz, EDT, coherence / IACC, fine structure, ringing. Outside: (1) **the reverb level against the direct path fell with RT in the model** (the capture's early reverb / tap ratio is RT-independent, §11.1 / the first-day check, −33.5 ± 0.4 dB; the model needed +3.3 dB at RT 1.0 … −0.6 at 4.5, +11.8 at 0.4), which the earlier optimal-gain scores (ρ²) could not see; (2) high bands decaying 5–7.5 % too long at 8–12.4 kHz; (3) L's late level 0.5–0.7 dB low; (4) echo density and mixing time fail for impulse renders but pass for sweep captures (the captures are sweep-deconvolved: judge density on sweep captures); a 10 s capture sweep lengthens the late high-band decay far beyond the data, a ~2 s sweep matches it (Audio Ease's capture smeared little, §13.6).

**17.2 The moving reads read the loop line before the loop filter. [S]** Fitted per RT with after-filter atoms, the nine reads' amplitudes rise at short RT exactly as 1 / |G| of the loop filter in the low band (normalised, RT 1.0 … 4.5: fitted 1.22 / 1.04 / 0.98 / 1.00 / 0.97 / 0.91 / 0.91 / 0.97; 1 / |G(125 Hz)| 1.28 / 1.10 / 1.01 / 0.97 / 0.94 / 0.91 / 0.90 / 0.89). So in each loop **v[t] = u[t] − line[t − P], line = G{v}; the static reads take line (after G, their first output ∝ G, §11.2), the moving reads take v (before G)**. With that the gain the data needs over the model's early reverb (30–165 ms, R) is −0.2…−0.6 dB at every RT (was +3.3…−0.6), and the audit's per-RT level gains are 0.995–1.07 (≤ 0.6 dB) including RT 0.4. Both forms fit the low band equally (a normalised measure cannot see a per-RT scale); the level decides.

**17.3 HF from the decay statistics: 0.69. [S]** Octave T60 of impulse renders (two seeds) against the data, fitted on four RT steps and scored on the other four: HF 0.69 in both splits; held out the T60 error is **0.6–2.1 % in every octave from 63 Hz to 8 kHz and in 10–12.4 kHz** (at the IR-fit 0.736: +4.7 / +7.0 % at 8 / 10–12.4 kHz). LF 1.49 confirmed in the same run. As for RT (§13.5) and LF (§16.5), the IR least-squares value was biased by the random part.

**17.4 The moving reads' broadband levels.** The low-band amplitude fit sets R about 1 dB high broadband: a ×0.89 scale puts the moving / static energy ratios on the data (−16.0 / +4.2 / +5.7 / +4.7 against −15.5 / +4.4 / +5.6 / +5.0 dB) and the held-out statistics cost at 91 (data against itself 100); a one-pole low-pass instead is no better (98–134). L's moving reads need ×1.12 relative to that (the L / R energy balance over time, held out by RT: within ± 0.25 dB in every window to 1.6 s; a small early-high / late-low residual remains), i.e. L at its low-band value. The low band stays within 0.3 dB of the data (63–500 Hz: +0.3…−0.1 dB). Per-channel decay slopes already agreed (T60 L / R within 0.02 of the data's).

**17.5 After the fixes** (impulse renders, inside / total, before → after): 

| Family (impulse renders, all nine RTs, both channels) | First audit (v1) | After 17.2–17.6 (v4) |
|---|---|---|
| Octave T60 L / R | 40 / 39 of 56 | **53 / 53** |
| Third-octave T60 L / R | 129 / 117 of 152 | **141 / 138** |
| Octave EDT L / R | 45 / 46 of 56 | **49 / 53** |
| Octave late decay L / R | 52 / 52 of 56 | **55 / 56** |
| C50 L / R | 23 / 32 of 56 | **41 / 41** |
| D50 L / R | 23 / 33 of 56 | **41 / 41** |
| Ts L / R | 19 / 28 of 56 | **44 / 43** |
| Early-window level L / R | 10 / 8 of 32 | **26 / 28** |
| Early level L / R | 103 / 126 of 152 | **131 / 140** |
| Late level L / R | 71 / 104 of 152 | **126 / 133** |
| L / R balance early / late (broadband) | 2 / 0 of 8 | **8 / 5** |
| L / R balance late (per band) | 35 of 56 | **55** |
| IACC (time-resolved) | 48 of 48 | 48 |
| Broadband T60 / EDT | 7 of 8 | **8 of 8** |
| **All families** (echo density and mixing time excluded: judged on sweep captures, ned_b 8 / 8, ned_a 6 / 8, mixing time 4 / 8) | **2986 / 3832 (78 %)** | **3402 / 3832 (89 %)** |

The level gain the model needs per RT went from +3.3 dB (RT 1.0) / +11.8 dB (RT 0.4) to 0.956–1.027 (≤ ± 0.4 dB) at every RT, RT 0.4 included.

**17.6 The reverb's level against the direct tap: +0.35 dB. [S]** After 17.2–17.4 every window after 30 ms sat 0.3–0.45 dB low against the (exactly matched) direct tap, which read as C50 +0.4–0.8 dB and Ts 4–7 ms short; the audit's level gains said the same (1.03–1.07 at every RT). One gain on the whole reverb (the static reads and the moving reads together, `rev_gain` in both files) leaves the energy ratios, the direct path and the decay untouched; on three RTs it took C50 R / L +0.44 / +0.47 → +0.27 / +0.07 dB and Ts R / L −4.4 / −6.8 → −0.1 / −1.4 ms (v4 above includes it). So §17.4's ×0.89 on the moving reads was half the story: part of the moving / static excess was the static part being low against the direct path.

**17.7 A gentle low-pass on the moving reads (one pole, 0.07). [H]** v4's remaining outlier was the top band: early level right (6.4 / 10 kHz −0.08 / +0.24 dB) but late level high (+0.65…0.87 / +1.6…1.8 dB) with the 8 kHz decay slope right, i.e. C50 8 kHz −1.05 dB and Ts +8–9 ms in both channels: the moving reads (before the loop filter, §17.2) carry a little too much top end. Same three RTs, like for like: no low-pass / pole 0.07 / 0.15: late 6.4 kHz +0.80 / +0.35 / −0.22 dB, 10 kHz +1.78 / +0.99 / +0.10; C50 8 kHz −1.08 / −0.46 / +0.35 dB; Ts 8 kHz +9.3 / +3.6 / −3.1 ms; all checks 1254 / 1248 / 1226 of 1437; statistics cost 90.9 / 89.9 / 94.1. Adopted 0.07 (`ADOPTED["mlp"]`; the kind of roll-off interpolation or address stepping gives; below 1 kHz < 0.05 dB, so the low-band fit is untouched). Left: late 10 kHz about +1 dB (the band limit's edge). **The full nine-RT audit with the low-pass alone (v5) was 88 % against v4's 89 %**: the high and extended bands gained (+8 checks), clarity and the early-window level lost (the low-pass removes ~0.15 dB of reverb energy, so the level gains rose to 1.01–1.05). **Adopted with that energy restored: the reverb gain is +0.5 dB in total (17.6's 0.35 + 0.15).** Same three RTs, like for like, in-band / extended checks: no low-pass +0.35 dB 1254 / 438, low-pass +0.35 dB 1248 / 442, **low-pass +0.5 dB 1269 / 449** (of 1437 / 552); early-window level 10 / 5 / 10 of 12 (R). Statistics cost 89.9; static references unchanged in form (−150 dB against the frequency-domain model).

**17.8 The rate and shape of the glide, E's onset and level. [S / H]**
- **Rate from E's own data.** With the depth now E's own (the low band, §16.2), the lobes' speed fixes the rate: held-out cost 104 / 102 / 91 / 98 / 120 at 0.45 / 0.50 / 0.553 / 0.608 / 0.70 Hz (seed scatter ± 15): **E's LFO is at 0.5–0.65 Hz** (S), the plugin's 0.553 Hz applies. 0.553 against 0.608 (the 250's rate kept in seconds or in samples) does not separate (8 seeds: 94.4 ± 5.2 / 97.8 ± 7.0); 0.553 stays the default (the delays were re-timed in milliseconds, H).
- **Shape.** Triangle and sine at equal depth score alike (91 / 92); the plugin's snapshots (constant speed) decide for the triangle (H for E).
- **Onset.** First arrivals in the Altiverb captures (one capture chain): S 10.17, M 20.33, L 40.42 ms (×1 / ×2 / ×4) and **E 20.33 ms (R), with M's**, L 32.25. At DELAY 0 E's direct tap arrives with M's first cluster tap. [S]
- **Level.** The Altiverb normalised every preset to the same total energy (64.33–64.63 dB over all 32 S / M / L / E presets), so the programs' absolute levels are not in the data. Calibration (the spec's `controls.level`): **E's direct tap 1 at the level of an S / M / L cluster tap at the same AMPL** (one converter and headroom for every program). [H]

**17.9 Status.** Everything the existing data constrains is fitted and held out, and the G1–G6 audit is at 89 % inside the model's own realisation spread (S / M / L's pulse check is the precedent). Left, not resolvable from existing data: AMPL / DELAY / LP for E, the absolute level (calibrated by rule), 0.553 vs 0.608 Hz, the glide's shape in E, the 16-bit mode for E (a plugin item, modelled as S's), and small residuals: the direct path 5–6 dB over the data's floor at 20–100 Hz and 10–12.4 kHz (capture chain), the late 10 kHz level (+1 dB), L / R over time ± 0.25 dB.

