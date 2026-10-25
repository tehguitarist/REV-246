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
