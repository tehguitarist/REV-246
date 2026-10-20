# EMT 246 algorithm research: lineage, classic structures, and a recreation plan

> **Status note (Oct 2026):** the hypotheses here were later tested against the CLA IR (`model_s_results.md` §8–13). Since then:
> - a diffuser AP20 · AP72 · AP230 + comb460 has been found;
> - randomised injection is ruled out;
> - the "late input" is explained;
> - the topology ranking in §7.1/§7.2a is out of date.
>
> Current facts are in `EMT246_unit_facts.md`.

*Scope: what the literature, patents and documented hardware say about how the 246's S/M/L programs (and, more briefly, program N) were most likely built, and how to identify and fit a model from IRs at the device rate (fs = 19.2 MHz / 728 = 26,373.6 Hz). This document is a companion to `EMT246_discovery.md` and `analysis_harness_spec.md`.*

**Labels used on every claim:**
- **[C: source]** means confirmed by the cited primary or near-primary source.
- **[M]** means measured by us (from the discovery doc or coordinator notes; a single IR unless stated).
- **[I]** means an inference from confirmed or measured facts.
- **[S]** means speculation: plausible, but nothing points to it directly.

Reference numbers like [R3] point to §9.

---

## 0. Key conclusions

1. **The 246's vocabulary is Blesser's.** In a 1997 letter printed in Dattorro's JAES paper, Barry Blesser frames reverb design entirely in terms of *eigentone* (Eigenton) density, and says the EMT 250 got its eigentone density "entirely by randomizing delay lines" [C: R9 App. 1]. The 246 brochure's "optimised eigentone distribution" therefore reads as the deliberate **time-invariant replacement** for the 250's randomisation [I].
2. **Uncoupled harmonic loops fit Blesser's own advice.** Our measurement of exactly harmonic loop modes [M] matches that letter's warning: an all-pass inside a large loop makes loop time and RT vary with frequency and leaves the result "very colored" [C: R9]. It also rules out coupled FDNs (Stautner–Puckette, Jot), Griesinger/Dattorro figure-8 tanks, and 250-style modulation [I].
3. **The first reverb pass is too dense for read-out taps alone.** The first 1,300 samples after onset have kurtosis 5.3–5.6, and the top 100 samples hold only 53 % of the energy [M]. That implies at least about 200 non-negligible samples per channel. At roughly 4–6 TMS32010 cycles per external tap, 2 × 182 cycles allow only about 50–80 taps in total for both channels and every function [I].
   - So **the most likely topology is a recursive diffuser (short all-pass loops, as in the 250 patent) feeding four uncoupled damped loops, read by a few output taps per channel** [I].
   - A pure 250-style "many exponentially-weighted taps per long loop" read-out is the second hypothesis.
   - Both make sharp, different, testable predictions (§7.2).
4. **"Optimised eigentone distribution" has two readings, and the IRs can tell them apart.**
   - One reading is the choice of loop lengths.
   - The other is equalised *mode amplitudes*: a narrow distribution of modal excitation, which Heldmann & Schlecht showed is what makes a reverb sound colourless [C: R27].
   - The four S lengths are unremarkable on the frequency-spacing and echo-coincidence metrics we tried [M, §4.3], which points towards the amplitude reading. Test: the spread of per-harmonic magnitudes per loop (§7.2, T2).
5. **"Adapted density" should mean more total loop memory for longer-RT programs.** Blesser notes that the required eigentone density depends on the nominal RT [C: R9].
   - Predicted total loop memory: M about 10–12 k words, L about 30–39 k words (memory-capped). Most likely still about 4 loops, but longer ones [I/S].
   - L's loops would then be around 300 ms, i.e. close to the 250 patent's "≥ 300 ms" long loops [S].
6. **Fitting needs grey-box least squares, not black-box optimisation.** With integer loop lengths known, the IR is *linear* in the injection and read-out coefficients, and the loop filters can be read directly from per-mode frequency offsets and decay rates. GA and differentiable-FDN methods are useful only as a final polish (§5).
7. **Update after running T1/T2 on the CLA IR (§7.2a):** mode magnitudes are noise-like, which disfavours a pure all-pass diffuser with single read taps. Every loop's L/R relation uses the same two lags (26 and 332 samples), which points to a fixed post-loop output distributor as in the 250 patent. The ranking now favours multi-tap (C/B) designs.
8. **Program N is most plausibly an LTI construction.** Candidates are an early section that is truncated at SPAN (cancellation taps, i.e. a truncated IIR, or a gain-shaped multi-tap FIR), summed with a separate late loop set.
   - This is the same idea as the Lexicon 480L SHAPE/SPREAD envelope control [C: R31, mechanism I].
   - It is testable from EDC shape, periodicity and tap counting (§6).

---

## 1. EMT lineage

### 1.1 EMT 250 patent: US 4,181,820 / DE 27 19 276 (Blesser & Bäder; priority 29 Apr 1977) [R1, R2]

Every structural statement below is **[C: R1]** unless marked otherwise.

- **Overall chain** (Fig. 1):
  - an input delay unit t_v ("0…60 ms" is legible in the figure);
  - a "reverberation generating means" made of loops t2–t16 with different delays, roughly 15 loops;
  - an "output discorrelator and distributor" with its own delays t17–t19 and four outputs (mono, stereo or quad).
  - The patent does not say how the short and long loops are interconnected (series or parallel) [I from text].
- **Two loop classes:** short loops "having an all-pass characteristic", and long loops "with tappings". The long loops' total delay is "preferably at least 300 ms". Claim 5's "300 µs" is an OCR or typo artefact; the description says ms.
- **Long-loop structure** (Fig. 2):
  - The input feeds the delay t_A via summer 6.
  - The end of the line feeds back through damping element g1, which also goes to the output.
  - Taps A2…An feed damping elements g2…gn into summer 7, then to the loop output.
  - The tap gains are chosen so that the *energy* of successive taps follows an exponential ("e-function"), "and not, as was initially obvious, the voltage distribution". The stated aim is to avoid the **stepwise** decay of a loop that is read only once per period (listed as problem 2).
- **Tap placement:** the tap delays, and their "sums, differences and quotients", are in a "non-harmonic relationship". The preferred rule is **T_n = ⁿ√k**, where n is the tap number and k the total delay. This was confirmed from the scanned PDF, where the root index is n [C: R1 PDF image].
  - **Quantitative reading [I]:** taken literally (k = 300 ms), it gives taps at 300, 17.3, 6.7, 4.2, 3.1 … ms, crowding at the start of the line. That is implausible as a complete read-out pattern.
  - The usable content is the *design criterion*: tap delays, their pairwise sums and differences, and their ratios should not be close to small-integer ratios of one another. In integer-sample terms that roughly means distinct pairwise differences (a Golomb-ruler-like set) and no near-rational ratios with small denominators.
  - The measured 246 cluster positions {0, 54, 143, 194, 288, 353} happen to have all 15 pairwise differences distinct [M/I]. That is weak evidence, since about 60 % of random 6-sets would too [I].
- **Frequency-dependent feedback:** the feedback factors of the long loops are frequency dependent, set by front-panel controls t_L, t_M, t_H (low, mid and high decay). The patent gives no filter topology.
- **Early-density path:** a parallel delay t5 with **3 taps** B1–B3 and gains h1–h3 adds "a delayed and repeated combination of unreverberated input signals" to the loop output. Its gains are "damped as a function of the reverberation time … in accordance with an e-function". The aim is to fix the slow build-up (problem 4).
- **Eigentone density by randomisation:** "Since the ear picks up a time sequence of eigentones at different frequencies as a greater density of eigentones", randomly controlled delay-length variation increases the effective eigentone density and suppresses beats (Fig. 3: random-noise-controlled feedback gains f_z1, f_z2 at several points along the line).
  - Fig. 4 adds **random mixing** of the loop outputs (gains g21–g216 driven by noise generators).
- **Hardware assumptions:** at most 100 ns per instruction, an ADC of at least 10 bits, optional gain-ranging.

**Other facts about the 250 (from Blesser via Dattorro):** the production EMT 250 ran at **32 kHz with an 8 K-word main memory** [C: R9 footnote 5 and App. 1]. By comparison, the 246 has 64 K words at 26.4 kHz, eight times the memory.

### 1.2 Blesser's own account (Dattorro 1997, Appendix 1, "Reverberation Recollections") [R9]

This is the most useful primary text found. All points are **[C: R9]**:
- In the 250, "the required eigentone density was emulated entirely by randomizing delay lines."
- All-pass "colourless" reverberation is *perceptually* not colourless. Many cascaded all-passes sound "like a machine shop", because frequency regions get bunched in time. The value of an all-pass is that every downstream element sees the same spectral density, so comb peaks do not align. "Parallel structures of non-all-pass elements achieve a similar" effect.
- "All-pass elements are more critical for small delay values." An all-pass inside a larger loop "must be used with great care" because its group delay varies with frequency. The effective loop time and RT then vary with frequency, and after many trips the result is very coloured.
- Schroeder's 3 eigentones/Hz criterion (about 3 s of memory) was empirical. Blesser explains it as random beating between modes within the mode bandwidth. He stresses that the required density **depends on the nominal RT**, because the RT sets how many eigentones a narrow-band input excites.
- With early boxes (about 150 ms of memory), envelopes showed a roughly 6 Hz periodicity, which "sounds bad".
- Physical rooms have eigentone density ∝ f², while electronic structures have constant density.

**Implication [I]:** the 246 S/M/L design looks like a direct application of these principles.
- Uncoupled loops with no in-loop all-passes (frequency-independent loop time, so a uniform RT).
- Diffusion kept outside the loops, so that every loop "sees the full spectrum".
- Loop memory scaled with RT ("adapted density").
- A static, optimised mode distribution instead of randomisation.

### 1.3 EMT 244/245/251/252

Very little is documented.
- **EMT 251:** added an LCD, improved HF response and more flexibility, and used "not the same algorithm" as the 250 [C: R5, Mix, secondary]. It also had a Non-Lin program giving a "non-linear reverberant delay-time characteristic" [C: R6 dealer page via search snippet; not fetched, 403].
- **EMT 252:** listed with a "non-linear program" [C: search snippet, R6].
- **EMT 244 (≈1980) / 245 (1981):** Blesser algorithms, 13-bit converters, RT 0.4–4.5 s; the 245 added pre-delay and early-reflection level [C: R7 review; R8, secondary].
  - The PSP 2445 is described as based on Blesser's original algorithms, with Blesser consulting [C: R7, vendor claim].
  - No structural description of the 244/245 was found in any document [I: absence after search].
- **Other credits:** Mix credits the 250's digital hardware to Ralph Zaorski (Dynatron) and the algorithms to Blesser and Bäder [C: R5].
- **Designer of the 246:** no source names one; Bäder remains plausible [S].
- **Patents:** no patent was found for the 244, 245, 251, 252 or 246 algorithms. The 250 family's citing patents are later Blesser "artificial ambiance" filings (US 7,062,337, from 2006) and unrelated others [C: R1 "Cited by" list].

### 1.4 Blesser writings

- **"An interdisciplinary synthesis of reverberation viewpoints"** (JAES 49(10):867–903, Oct 2001) [R10]: paywalled; only the abstract fragment was seen. Its perceptual framing argues that an algorithm's performance can be predicted without extensive ad-hoc listening tests [C: R10 abstract snippet]. The full text is worth obtaining for its discussion of eigentone and echo-density criteria [S that it contains them].
- **Spaces Speak, Are You Listening?** (Blesser & Salter, MIT Press 2007/2009) [R11]: argues that random time variation in algorithmic reverbs emulates thermal fluctuation in real spaces [C: R12 secondary summary].
  - That philosophy fits the 250, but the 246 S/M/L are static [M], which further suggests either a different designer or a deliberate change of approach after Blesser's EMT period [S].
- **Interview with David Griesinger** (Recording Musician, 1992) [R13]: Griesinger says Blesser, then consulting for EMT, pointed him to Schroeder 1962 [C: R13]. He also explains that static designs let some modes ring 2–3× longer than others, and that the Lexicon 224 used random delay changes to shorten them [C: R13].

---

## 2. Classic algorithms and how each would show in our IR

**TMS32010 cost model [now C: TMS32010NL confirmed on the processor board; timings verified in the TI TMS32010 User's Guide 1985: IN/OUT 2 cycles, TBLR/TBLW 3, LT/LTA/LTD/MPY/APAC 1, 144-word data RAM, cycle = CLKIN/4; see EMT246_discovery.md §2.3a].** One instruction is 4 clocks, so 19.2 MHz / 4 / fs = **182 instructions per sample per DSP** [M/C arithmetic]. The 144-word on-chip RAM is too small for any delay line, so all delays sit in external RAM.
- From memory, and to be verified against the TMS320C1x User's Guide: IN/OUT take 2 cycles and TBLR 3 cycles, with no single-cycle MAC (LT/LTA/LTD + MPY = 2 cycles per product).
- So one external delay tap costs about 4 cycles if a hardware address counter supplies the address, or about 6 if the DSP must OUT the address first.
- **Budget:** about 30–45 taps per DSP, or roughly 50–80 for the whole machine, before filters and housekeeping.

| Algorithm | Structure | Static / integer? | IR signature | Est. cycles (mono→stereo) | Fit to 246 S |
|---|---|---|---|---|---|
| Schroeder 1961/62 [R14, R15] | 4–6 parallel combs (30–45 ms, mutually prime, ≈1:1.5 span) + 2–5 series all-passes (≈100/3ⁱ ms) | yes / yes | Harmonic comb modes; all-passes in *series outside* the combs leave modes harmonic. Echo density builds slowly. Modal density = ΣPᵢ/fs | ~60–80 | **Topologically close** (parallel uncoupled combs), but the 246 has matched decays and much faster density |
| Moorer 1979 [R16] | 6 combs with **one-pole lowpass in the feedback path** + early-reflection FIR (~19 taps) + all-pass | yes / yes | As Schroeder; per-loop HF decay; effective loop length P + τ_g(LP) | ~120–160 | **Very close**: lowpass combs explain the P + 0.2 offset; the ER FIR is like the cluster |
| Gerzon 1971/76 [R17] | Unitary (energy-preserving) feedback matrix between delays | yes / yes | Modes are eigenvalues of the coupled system, *not* on individual loop harmonic series | ~80+ | **Ruled out** by exact per-loop harmonic modes [I] |
| Stautner & Puckette 1982 [R18] | 4-delay FDN with orthogonal mixing matrix | yes / yes | Same as Gerzon | ~80 | **Ruled out** [I] |
| Jot & Chaigne 1991/92 [R19, R20] | FDN + per-delay absorptive filters with dB loss ∝ delay length, + tone-correction filter | yes / yes | Matched decay at all frequencies; coupled modes | ~100 | **Its filter-design rule is directly applicable** to our uncoupled loops; its coupling matrix is ruled out |
| Griesinger / Lexicon 224/480, Dattorro 1997 [R9, R13, R31] | Input all-pass diffusers → figure-8 "tank" with **all-passes inside loops**, modulated delays, many output taps | no (modulated) / no (interpolated) | Inharmonic, slowly drifting modes; time-variant | ~200+ | **Ruled out** for S (harmonic, static) [I]; its output multi-tap scheme is relevant |
| 250-style tapped long loops [R1] | Long loops (≥300 ms), many read taps with exponential energy, short all-pass loops, randomisation | no (random) / yes | Time-variant; smooth intra-period decay | taps × 4–6 | **The static subset is hypothesis B** (§7.1) |
| Ursa Major SST-282 [R22] | One 255 ms buffer, ~15 summed taps fed back | partly (taps wander) | Strong multi-tap periodicities | low | Different (multi-tap feedback, not parallel loops) |
| AMS RMX16 [R23] | Undocumented; described as "a static reverb" | yes (per R23) | — | — | No structural information |
| Yamaha REV1, Sony DRE-2000, Quantec QRS (1982), Eventide SP2016 | Undocumented in public sources found. The SP2016 reportedly uses no modulation [C: forum snippet] | ? | — | — | No information; not useful as a template |

**Signature checks against our measurements [I]:**
- Exact harmonic series per loop, an autocorrelation with peaks only at multiples of Pᵢ, and a time-invariant tail together mean **diagonal feedback (uncoupled loops), no in-loop all-pass, and no modulation**.
- That leaves Schroeder/Moorer-type parallel lowpass combs and static 250-type loops as the only families consistent with the tail.
- The differences between those families lie entirely in **how the loops are fed and read**.

---

## 3. Colouration and mode distribution theory

### 3.1 Criteria

- **Mode density** of parallel combs: ΣPᵢ/fs modes per Hz [C: R24 JOS "Mode Density Requirement"].
- **Schroeder & Logan:** 0.15 modes/Hz is adequate for RT = 1 s, and the requirement scales with RT. That gives ΣPᵢ ≥ 0.15·RT·fs [C: R24].
- **Modal overlap:** the mode half-power bandwidth is Δf = 6.91/(π·RT) ≈ 2.2/RT Hz. The room-acoustics criterion (Schroeder frequency) needs about 3 modes per bandwidth. Blesser and Dattorro quote Schroeder's "3 eigentones/Hz" [C: R9]. Dattorro notes his own good-sounding plate violates it (it would need about 88 K words) and that echo density matters as much [C: R9].
- **Mutually prime lengths** maximise the repetition period of the lossless prototype [C: R25].
- **Coloration ∝ spread of modal *excitation* (residue magnitudes).** Narrow distributions are perceived as colourless and as higher modal density [C: R27]. The modal decomposition of FDNs shows that few modes carry most energy unless this spread is controlled [C: R26 abstract].
- **Griesinger:** colouration comes from a spread of *decay rates* across modes [C: R13]. The 246's matched loop decays [M] remove that cause. Only Jot-style per-loop filter scaling keeps that true at LF and HF [I].

### 3.2 Numbers for program S [M arithmetic]

| Quantity | Value |
|---|---|
| ΣPᵢ | 5,758 samples = 0.218 s, i.e. 0.218 modes/Hz (mean spacing 4.6 Hz) |
| Max RT that meets Schroeder–Logan | 0.218/0.15 = **1.46 s** (S's range is 0.2–2.5 s) |
| Modal overlap (modes per 3 dB bandwidth) | 2.4 at RT 0.2 s, 0.48 at RT 1 s, 0.19 at RT 2.5 s |
| Prime factors | 1321 (prime); 1389 = 3·463; 1461 = 3·487; 1587 = 3·23² |
| Exact shared modes | All three ×3 loops share a mode at fs/3 = **8,791 Hz**. No other coincidences below Nyquist |
| Length ratios | 1 : 1.0515 : 1.106 : 1.2014. The first three are nearly geometric (r ≈ 1.0515); 1587 breaks the series (1321·r³ ≈ 1536) [S: possible design trace] |
| Memory left with 999 ms pre-delay | 65,536 − 26,347 = **39,189 words** (1.49 s) for loops, cluster and diffusion |

### 3.3 Are the S lengths "optimal"? [M, quick Monte-Carlo against 5,000 random 4-sets in 1300–1610 with ≥ 20-sample separation]

- **Frequency spacing** (100–8,000 Hz union of harmonic series):
  - Coefficient of variation of nearest-neighbour spacing: 29th percentile (lower is more uniform).
  - Minimum spacing relative to the mean: 92nd percentile.
- **Time-domain echo coincidence** (first time any two loops' echo trains coincide within ±1 sample):
  - The earliest is 1321 vs 1389 at 7.2 s; the median for random sets is 4.9 s. That is the 65th percentile.
  - All pairs coincide later than the longest RT in S [I].
- **Verdict [I]:** the lengths are decent but not extreme on any simple metric. Other constraints probably dominated:
  - memory layout or addressing (the shared factor 3 might reflect interleaved or 3-phase processing) [S];
  - echo-arrival-time distribution;
  - choice by ear (Schroeder himself admitted picking delays "random[ly]" [C: R9 App. 1]).

### 3.4 What "optimierte Eigentonverteilung" might mean

- **(a) Frequency distribution:** loop lengths chosen so that mode frequencies interleave evenly. The weak metrics above slightly argue against this as the main meaning [I].
- **(b) Amplitude distribution:** every mode excited at near-equal magnitude, i.e. narrow modal excitation [R27]. Two things could achieve it:
  - a diffuser or injection whose spectrum is flat at every loop's harmonics (an all-pass is ideal; one read tap per loop preserves flatness);
  - per-loop gains that balance the four series.

  Random sparse multi-tap read-outs give Rayleigh-distributed mode magnitudes, which are *not* narrow [I].
- **(c) Decay distribution:** matched RT across loops and frequencies (Griesinger's criterion) [I].
- **Most likely [I]:** (b) + (c), with (a) as a secondary constraint. This is decidable from the IR (§7.2, T2).

### 3.5 "Adapted density" predictions for M and L [I/S]

- **If loop memory scales with max RT as in S** (5,758 words for 2.5 s, i.e. 0.58× the Schroeder–Logan requirement): M (5 s) needs about 11.5 k words and L (20 s) about 46 k words.
- L exceeds the 39 k words left with full pre-delay. So L is probably **memory-capped at about 30–39 k words**, or the pre-delay range is reduced on L (testable).
- Loop count costs cycles while loop length costs nothing, so **4 loops of about 2.6–3.2 k samples (M) and about 7–9.8 k samples (L)** is the cheapest design [I].
- Longer loops lower echo density per recirculation. So L in particular needs more diffusion or more read taps, which brings it close to the 250's long tapped loops [S].

---

## 4. Loop filters for LF and HF

- **Jot's rule** [C: R19, R20; formulas per R24/R28]: each loop filter's dB attenuation should be proportional to its length: A_i(ω) = −60·Pᵢ/(fs·RT(ω)) dB.
  - For a first-order filter H_i(z) = g_i(1−p_i)/(1−p_i z⁻¹) with per-pass DC and Nyquist gains R₀ and R_π:
    - p_i = (R₀ − R_π)/(R₀ + R_π)
    - g_i = 2R₀R_π/(R₀ + R_π)
    - Rω ≈ 1 − 6.91·Pᵢ/(fs·RT(ω)) [C: R28].
  - The loop gain law gᵢ = 10^(−3Pᵢ/(fs·RT)) that we measured [M] is this rule at mid frequency.
- **What the 246's ranges imply [I, our arithmetic]:**
  - **HF ×0.9 at RT 1 s** needs an extra 6.7 dB/s at 7 kHz, i.e. 0.33–0.40 dB per pass. A one-pole lowpass does this with p ≈ 0.034–0.041, giving a group delay of only **0.035–0.042 samples**.
  - **HF ×0.1** with a plain one-pole lowpass would need p ≈ 0.94–0.97. That would wreck the mid band, so **the HF control must be a shelving section** (first-order high-shelf or equivalent).
  - **LF ×0.1–2 at 50 Hz** likewise needs a low-shelf (for LF > 1, the loop is lossier at mid than at LF).
  - The measured **P + 0.2** offset is therefore *not* explained by the HF damping at the probable setting (HF 0.9). Candidate explanations:
    - the LF shelf's group delay;
    - a fixed 2-tap FIR element (e.g. 0.8 + 0.2z⁻¹);
    - a clock mismatch between the capture interface and the unit (0.2/1440 ≈ 140 ppm).
  - **Discriminating test:** a clock offset scales with Pᵢ; a filter offset is the same for all loops and varies with frequency.
  - **Update [M]:** the harness measured this capture's clock at **−99 ppm** relative to 19.2 MHz/728 (from cluster-tap spacing). That alone stretches the effective lengths by 0.13 samples (P = 1321) to 0.16 samples (P = 1587). So most of the "+0.2" is the clock, and very little group delay is left for the loop filter.
  - **Correction (second CLA image, numeric display):** the IR was taken at **HF ×0.9, LF ×1.0, RT 1.0 s**, not HF ×1.0. The measured RT(7 kHz)/RT(1 kHz) of 0.90–0.91 is therefore exactly the HF setting: the HF factor is applied at 7 kHz as the data sheet says.
    - At this setting it is about 0.37–0.40 dB of extra loss per pass at 7 kHz; a mild one-pole (p ≈ 0.04) fits it.
    - Whether HF ×1.0 still has residual fixed damping is unknown until an HF sweep is available.
- **Whether filters are per-loop (Jot) or shared [S → test]:** EMT might have used one coefficient set for all loops to save table space. Then the per-second HF loss differs between loops in the ratio 1587:1321 (about 20 %). Example: identical p = 0.05 gives 9.9 against 8.3 dB/s at 7 kHz. This is measurable with the per-loop separation.

---

## 5. Modern fitting methods and their applicability

| Method | What it does | Applicability here |
|---|---|---|
| Coggin & Pirkle 2016 [R29] | GA over FDN parameters to match IR features | Low: our topology and integer delays are already known; a GA wastes the structure |
| Ibnyahya & Reiss 2022 [R30] | IR features + GA with MFCC loss | Low, same reason; perceptual losses are useful only as a final check |
| Shen & Duraiswami 2020 [R32] | Data-driven pre-tuned FDNs from room parameters | Not applicable |
| Dal Santo et al. 2023 DAFx; 2024–25 "tiny colorless FDNs" [R33, R34]; Mezza et al. 2024 [R35] | Differentiable FDN; spectral-flatness + sparsity losses; learnable delays | Medium: useful to **polish** loop filters, injection FIRs and output gains jointly by gradient descent once the structure is fixed. Their colourlessness loss is also a way to ask whether EMT's design optimised mode amplitudes |
| Schlecht & Habets 2019 modal decomposition [R26] | Exact poles and residues of an FDN | High for **verification**: compare model poles and residues with measured per-loop harmonic decays and amplitudes |
| Heldmann & Schlecht 2021 [R27] | Coloration ↔ modal-excitation spread | High for interpreting "Eigentonverteilung" (T2) |
| ESPRIT / matrix pencil / FZ-ARMA (Karjalainen et al. 2002) [R36] | High-resolution estimation of mode frequency and decay | High, but simpler here: modes sit on *known* harmonic grids, so per-loop least squares at fixed frequencies (our comb projection) plus a per-harmonic decay fit is better conditioned than blind ESPRIT [I] |
| Karjalainen et al. 2002 decay + noise-floor fit [R37]; Schroeder backward integration | Band-wise T60 with a noise floor | High for LF/HF law fitting from library IRs |
| Abel & Huang 2006 NED [R38] | Echo density profile | Medium: for checking diffusion-model build-up against measurement |
| Välimäki et al. 2012 / 2016 reviews [R39, R40] | Surveys of all of the above | Background |

**Recommended approach [I]:** structured linear least squares.
- With integer Pᵢ and a parametric loop filter, the tail of each channel is y = Σᵢ Cᵢ(z)·uᵢ, where Cᵢ is the known damped comb and uᵢ is the unknown injection and read-out sequence.
- This is linear in uᵢ, so it can be solved exactly (as `combsep.py` already does in the scratchpad).
- The filter parameters are an outer loop of 2–4 nonlinear parameters per loop (Gauss–Newton or a grid).
- Each loop's harmonic family gives several hundred independent (frequency offset, decay) pairs, so |H_i(ω)| and τ_g,i(ω) are heavily over-determined.

---

## 6. Non-linear / gated programs and the 246 N

**What is documented:**
- **AMS RMX16 NonLin 2:** developed to reproduce the Townhouse SSL gated-room drum sound with feedback from Hugh Padgham [C: R23 and search snippets, secondary]. No structural documentation was found. LiquidSonics calls the RMX16 a static reverb [C: R23]. The IR envelope of one Townhouse NonLin preset shows a build-up of about 60 ms, about 100 ms of near-constant level, then a decay [C: R23 description of an IR].
- **Lexicon 480L:** SHAPE and SPREAD control the build-up and sustain of the envelope.
  - SPREAD ≈ sustain duration in ms. At mid SHAPE a secondary plateau appears, simulating a diffuse back-wall reflection; at high SHAPE it exceeds the early part, giving an "inverse" sound.
  - The 480L also offers diffused pre-echo clusters [C: R31 pp. 3-2/3-3].
  - The mechanism is not stated. Our inference is that it is a set of input taps (or a multi-tap FIR) with shaped gains feeding the tank [I].
- **Yamaha gate / reverse:** "reverberation quickly cut off" and "reverse increase in reverb" (MU80 manual) [C: R41]. No REV1 or REV7 internals were found.
- **EMT 251 Non-Lin:** "non-linear reverberant delay-time characteristic" [C: R6 snippet].
- **246 N controls:** DELAY (ms), SPAN (0.1–0.5 s, "Anhalldauer"), EARLY RT and LATE RT. The diagram shows crossing early and late decays [M/C: card and brochure via coordinator].

**Candidate implementations [I], from most to least likely:**
1. **N-A: early section truncated at SPAN plus a separate late loop set.**
   - The early section is either a gain-shaped sparse multi-tap FIR (tap gains following 10^(−3t/RT_early) up to t = SPAN), or a set of loops with RT_early truncated exactly at SPAN by cancellation taps. The truncation works because (1 − gᴷz^(−KP))/(1 − gz^(−P)) is FIR, which costs one extra tap per loop.
   - The late part is loops with RT_late, fed in parallel (with the crossover visible where it exceeds the early curve) or fed from a tap at SPAN.
   - It is cheap: N need not run S/M/L concurrently, so all 364 cycles are available.
2. **N-B: 250-style.** The 3-tap "early density" path, generalised to many taps with RT-dependent e-function gains, followed by the long loops.
3. **N-C: a level-triggered gate on the reverb output** (a truly non-linear system). This is unlikely: the parameter set has no threshold [I].

**Tests:**
- An EDC drop to the noise floor at SPAN within ≤ 1 sample means FIR or truncated IIR.
- Periodicity analysis of the early section: loop lags present means truncated loops; absent means FIR.
- Tap counting in the device domain.
- Whether tap positions scale with SPAN (a stretched pattern) or the pattern is fixed and gains are windowed.
- Whether the late part's loop set matches S, M or L.
- Level dependence (sweep vs impulse, two levels) to rule out N-C.

---

## 7. Synthesis: recreation plan

### 7.1 S/M/L topology, ranked [I]

Common to every hypothesis below [M]:
- input → LP (2–12 kHz) → pre-delay;
- 8-tap fixed cluster → outputs;
- reverb section;
- L/R sum → single DAC (L half a sample late);
- 4 uncoupled loops with gain gᵢ = 10^(−3Pᵢ/(fs·RT)) and LF/HF shelving damping.

The hypotheses differ in how the loops are fed and read:

| Rank | Hypothesis | Why |
|---|---|---|
| **A (≈45 %)** | **Recursive diffuser before the loops** (1–4 short all-pass or "all-pass character" loops, as in the 250 patent's short loops), common to all four loops, with **a few output taps per loop per channel** (different positions per channel, signed) | Explains: dense first pass (≥ ~200 significant samples, kurtosis 5.5); different folded waveforms per loop (a long diffuser response folded mod Pᵢ); fast Gaussianity; cheap (~10–12 cycles per all-pass stage); matches Blesser's "each element sees the same spectrum" |
| **B (≈30 %)** | **Static 250-style:** input written directly into each loop; **10–20 read taps per loop per channel** with exponentially weighted (possibly RT-dependent) gains; no diffuser | Explains: lineage, per-channel differences, smooth intra-period decay. Strained by: cycle budget (80–160 taps) and the 53 %-in-top-100 statistic |
| **C (≈20 %)** | **Hybrid:** short diffuser plus moderate multi-tap read-out (4–8 taps per loop per channel) | The most "EMT-like" mix of both patent elements; fits every statistic |
| D (≈5 %) | Per-loop input multi-tap FIRs (long sparse injection) with single read taps | Possible, but per-channel differences would then have to be pure circular shifts |
| — | Coupled FDN, Dattorro tank, modulated 250 | Ruled out (§2) |

### 7.2 Discriminating tests from existing and library IRs [I]

All tests work in the device domain (deconvolved, fs = 26,373.6 Hz) after per-loop separation of the periodic waveforms wᵢᴸ and wᵢᴿ.

- **T1, circular-shift test.**
  - The input is mono and each loop's injection is shared by L and R. So wᵢᴿ is a pure circular shift of wᵢᴸ **if and only if** each channel reads loop i at one tap (A or D).
  - Under B or C, Wᵢᴿ/Wᵢᴸ (over harmonics) is a ratio of sparse-tap spectra. Sparse deconvolution (OMP or L1) of wᵢᴿ against wᵢᴸ will reveal the relative tap sets.
- **T2, modal-excitation spread.**
  - Take the per-harmonic magnitude |Wᵢ(m)| after removing smooth LP and HF trends, and compute the standard deviation of the dB values.
  - About 5.6 dB (Rayleigh) means random multi-tap or noise-like injection (B).
  - Well under 1 dB means an all-pass diffuser read by one tap (A), and directly supports the "optimised eigentone (amplitude) distribution" reading (§3.4b).
- **T3, common-diffuser test.**
  - Recover the injections xᵢ = wᵢ ⊛ (1 − gᵢHᵢz^(−Pᵢ)) applied to the per-loop *non-periodic* response (first pass plus tail model).
  - Under A, all xᵢ are the same sequence d(n) up to delay and gain.
  - Equivalently, fit one d(n) (~2,000 taps) jointly so that fold_Pᵢ(d ∗ tᵢ) = wᵢ for all i.
- **T4, first-pass sparsity.**
  - Use the 1,300 samples after onset, with no recirculation yet [M: kurtosis 5.3–5.6, top 100 samples = 53 % of energy].
  - Recompute after envelope normalisation (local RMS over ~64 samples). An exponential build-up or decay alone inflates kurtosis.
  - Count isolated spikes above 4σ. B predicts ≤ ~20 spikes per loop per channel, single-sample, with no smooth envelope. A predicts noise with a smooth all-pass-like envelope.
- **T5, intra-period envelope vs RT.**
  - Under the patent's exponential tap weighting (B), the envelope of wᵢ across one period (ordered from the write point) falls by exactly the per-pass loss: 60·Pᵢ/(fs·RT) dB, i.e. ~3.3 dB at RT 1 s and **~16 dB at RT 0.2 s**. Its shape changes with the RT setting.
  - Under A, wᵢ's shape is independent of RT apart from overall scale.
  - Use S at RT 0.2 and 2.5 s.
- **T6, onset geometry.** First arrivals come 408 samples (L) and 576 samples (R) after the first cluster tap [M, coordinator; the discovery doc's earlier estimate was 580–650, so reconcile].
  - Under A, the per-channel onset = diffuser latency + the minimum read-tap offset, and the diffuser latency is common to L and R.
  - Under B, it equals each channel's earliest read tap. Its position relative to the loops' write points must then be consistent across all four loops.
- **T7, cluster and loops.** Difference IRs (AMPL 0 vs 100 %) show whether the cluster feeds the loops or only the outputs. The 250 patent's 3-tap path went only to the output [C: R1].
- **T8, LP position.** If changing LP alters loop decay rates, the LP is inside the loops; if it only tilts the spectrum, it sits outside (probably at the input, since it is applied "proportionally to cluster and reverb").
- **T9, filter offset.** Measure τᵢ(ω) = (effective length − Pᵢ) per loop vs frequency (§4). Constant across loops means a filter; ∝ Pᵢ means a clock offset.
- **T10, M/L memory.** Identify the loop sets of M and L (harmonic-comb search over integer P). Check ΣPᵢ against §3.5, the number of loops, whether lengths are scaled versions of S's, and whether pre-delay is limited on L.

### 7.2a First results of T1 and T2 on the CLA IR (added after the research pass) [M]

Script: `tools/topology_tests.py`. It runs per-loop periodic decomposition below 6 kHz on the CLA IR and on two synthetic controls with the same loops and RT:
- **control A:** a common 6-stage all-pass diffuser → loops → one read tap per channel;
- **control B:** direct injection → loops → 16 exponentially weighted random read taps per channel.

| Test | EMT 246 (CLA IR) | Control A | Control B |
|---|---|---|---|
| **T2** mode-magnitude spread (dB, 150 Hz–5 kHz, detrended) | **4.7–6.0** | 2.8 (the method's estimation floor) | 4.3–5.6 |
| **T1** PHAT L/R lag per loop (2nd/1st peak) | −332, −26, −332, −26 (0.8–0.97) | per-loop lags recovered (≈0.78) | random lags (0.77–0.9) |

**Reading:**
- **T2:** the 246's per-loop mode magnitudes are noise-like (Rayleigh ≈ 5.6 dB), like control B. They are clearly wider than an all-pass diffuser read by single taps.
  - This **argues against hypothesis A in its pure form**.
  - It also argues against the "equalised mode amplitudes" reading of *Eigentonverteilung* (§3.4b). On this evidence the optimisation is more likely about frequencies, decays, or something not yet identified.
  - Multi-tap read-out or injection (B, C or D) is favoured. That moves the ranking toward the user's 250-patent hypothesis.
- **T1:** the PHAT peak ratio does not separate the controls at this estimation quality (control A only reaches 0.78), so it is inconclusive as a single-tap test.
  - Notable instead: **the best L/R lag for every loop is one of the same two values, −26 or −332 samples.** These are also the two peaks of the broadband L/R cross-correlation.
  - Per-loop independent read taps would give unrelated lags (as control B does).
  - So the L/R relation is probably set by a **fixed output network after the loops are summed**. That would be a static version of the 250 patent's "output discorrelator and distributor" (delays t17–t19) **[I]**.
  - This should be modelled explicitly: loop sum(s) → a small multi-tap distributor → L and R. Fitting it is linear given the loop outputs.

**Revised ranking [I]:**
- **C (hybrid: some diffusion plus multi-tap injection/read-out) ≈ 40 %**
- **B (250-style multi-tap) ≈ 35 %**
- **D (multi-tap injection) ≈ 15 %**
- **A (pure diffuser plus single taps) ≈ 10 %**

All are combined with a common L/R output distributor. Tests T4–T6 on library IRs remain the deciders.

### 7.2b "Stochastic control": a 1986 claim that needs reconciling [C: Studio Magazin 99, Nov 1986]

The only period review says S/M/L use "eine neuartige, erstmals realisierte stochastische Steuerung … Eigentonverdichtung … ohne wahrnehmbare Frequenzmodulation" (a novel stochastic control giving eigentone densification without perceptible FM). It adds that S/M/L differ in "dem Grad der eingesetzten Stochastik".

Our CLA IR shows exactly harmonic loop modes and stable spectral fine structure. Three readings are compatible with both [I]:
1. **Slow or sparse randomisation of read-tap positions or gains** (not the loop lengths). Mode *frequencies* stay harmonic, but mode *amplitudes* (residues) drift. A sweep-measured IR would show this as the noise-like per-mode magnitudes we see (T2 ≈ 5 dB), plus excess noise in the tail. L's 6 dB higher floor and eugtone's L-only null error fit this if L carries more of the randomisation.
2. **Random coefficient switching between a few fixed configurations**, slow enough that one sweep mostly sees one configuration.
3. **The production firmware differs from the tested prototype**, or the reviewer described the design philosophy rather than a literal modulation.

**Tests:**
- ~~Repeat captures~~ (not possible: no hardware). Instead, use separate-capture duplicates if any library has them (harness A7 → E1), or joint fits across library IRs that share structure. A residual well above the noise floor that concentrates at mode frequencies means (1) or (2). Within one IR, mode linewidth in excess of what the decay predicts would point to time variance.
- Per-mode amplitude and phase stability across repeats (E2).
- Compare S, M and L for degree, since "Grad der Stochastik" differs per program.
- If confirmed, model it as randomised read-out gains or tap positions on top of the static loops. That is cheap and keeps the measured harmonic structure.

### 7.3 Step-by-step identification and fitting

1. **Data:** follow the harness spec.
   - Minimum per program: all 16 RT steps at neutral LF/HF; LF and HF sweeps at one RT; AMPL 0 and 100 %; two LP settings. Repeat captures are not possible, so the L-channel time-variance question rests on library duplicates and joint fits.
2. **Loop sets** for M and L (T10): autocorrelation and harmonic-comb search, then confirm with the four-periodic-component variance test.
3. **Per-loop separation:** least squares over the tail with known Pᵢ to get wᵢᴸ, wᵢᴿ, per-harmonic complex amplitudes and per-harmonic decay rates.
4. **Loop filter identification:**
   - From mode m of loop i: ω_m(Pᵢ + τᵢ(ω_m)) = 2πm gives the phase and group delay of Hᵢ.
   - The per-second decay σ_m gives |gᵢHᵢ(ω_m)| = exp(−σ_m·(Pᵢ + τᵢ)/fs).
   - Fit gain + one-pole high-shelf + one-pole low-shelf (or a biquad). Check Jot scaling vs shared coefficients.
   - Then fit the laws vs RT, LF and HF step indices (expect lookup tables; snap to likely 16-bit coefficient values).
5. **Feed and read structure:** run T1–T6 and choose A, B or C.
   - Then fit the integer tap positions (sparse recovery on the first pass + folded waveforms).
   - Fit all-pass lengths and coefficients for A/C. Use a grid over integer lengths; coefficients are likely "nice" binary fractions such as 0.5, 0.625 or 0.75 [S].
6. **Fallback (Path B in the discovery doc):** use the measured wᵢ (or the first-pass FIR) directly as injection. RT, LF and HF stay parametric through the identified filters.
7. **Assemble and null-test** at the device rate, then add the remaining non-idealities:
   - 16-bit storage with truncation;
   - the L half-sample DAC offset;
   - the converter/anti-image response;
   - possible L-only modulation.
8. **Polish (optional):** run a differentiable-FDN-style gradient refinement [R33, R35] of the continuous parameters only, with integer delays frozen.

### 7.4 Program N, in brief

- Capture SPAN at 0.1/0.3/0.5 s × EARLY RT at short and long × LATE RT at short and long, plus DELAY 0 and a two-level repeat.
- Apply the §6 tests:
  - drop sharpness at SPAN;
  - periodicities in the early and late sections;
  - spike counting;
  - whether tap positions scale with SPAN;
  - whether the late loop set matches S, M or L.
- Model N-A first. The early section is either a gain-windowed multi-tap FIR or truncated loops; the late section is the identified loops of one program.

---

## 8. Open questions and caveats

- Several facts rest on **secondary sources** (Mix, dealer pages, reviews, a vendor blog) and are flagged as such.
- ~~The TMS32010 cycle counts are from memory~~: verified against the TI User's Guide (October 2026).
- The Jot formulas were taken from JOS's online book, not the original preprint.
- The Monte-Carlo metrics in §3.3 are a quick screen, not proof of the designer's intent.
- **CLA IR settings, read from a second image's numeric display:** SHORT, LP 12, DELAY 10 ms, AMPL 80 % (faulty segments), LF ×1.0, RT 1.0 s, HF ×0.9. These supersede the LED-bar reading that suggested HF ×1.0 (EMT246_discovery §3.2).
- **Worth obtaining:**
  - Blesser 2001 JAES (full text);
  - Kuttruff, Frequenz 16(3) 1962 pp. 91–96 (tapped loops, cited by the 250 patent);
  - any 244/245/251/252 service documentation (algorithm block diagrams);
  - the Lexicon 224X/480L service manuals (lost-manuals.com, R31b), for how SHAPE/SPREAD and inverse programs were built.

---

## 9. References

*"Fetched" means read for this document. "Cited" means bibliographic details from search results or other papers; the full text was not read.*

- **R1** Blesser, B. A. & Bäder, K.-O., *Electric reverberation apparatus*, US 4,181,820 (1980; priority DE 2719276, 29 Apr 1977). https://patents.google.com/patent/US4181820A/en (fetched; local PDF/OCR in scratchpad)
- **R2** DE 2719276 A1, *Elektronisches Nachhallgerät*. https://patents.google.com/patent/DE2719276A1/de (fetched, OCR)
- **R5** Mix, "1976 EMT Model 250 digital reverb". https://www.mixonline.com/technology/1976-emt-model-250-digital-reverb-377973 (fetched)
- **R6** Vintage Digital EMT 251/252/246 pages. https://vintagedigital.com.au/emt-251 , https://www.vintagedigital.com.au/emt-246/ (403; search snippets only)
- **R7** Barry Rudolph, PSP 2445 review. https://and.barryrudolph.com/newtoys/toys16/pspaudioware_psp_2445.html (fetched); PSP product page https://pspaudioware.com/products/psp-2445-emt (cited)
- **R8** Sonicstate, "Five minutes with EMT 245". https://www.sonicstate.com/news/2019/09/27/five-minutes-with-emt-245-reverb/ (fetched)
- **R9** Dattorro, J., "Effect Design Part 1: Reverberator and Other Filters", JAES 45(9):660–684, 1997, incl. Appendix 1 (letter from B. Blesser). https://ccrma.stanford.edu/~dattorro/EffectDesignPart1.pdf (fetched)
- **R10** Blesser, B., "An interdisciplinary synthesis of reverberation viewpoints", JAES 49(10):867–903, 2001. https://aes2.org/publications/elibrary-page/?id=10176 (cited)
- **R11** Blesser, B. & Salter, L.-R., *Spaces Speak, Are You Listening?*, MIT Press (cited)
- **R12** Valhalla DSP, "Modulation in reverbs: reality and unreality". https://valhalladsp.com/?p=1563 (cited, search summary)
- **R13** White, P., "Illusions of Space" (Griesinger interview), Recording Musician, Aug 1992. https://www.muzines.co.uk/articles/illusions-of-space/2902 (fetched)
- **R14** Schroeder, M. R. & Logan, B. F., "Colorless artificial reverberation", JAES 9(3), 1961 (cited)
- **R15** Schroeder, M. R., "Natural sounding artificial reverberation", JAES 10(3):219–223, 1962 (cited); summary at https://ccrma.stanford.edu/~jos/pasp/Schroeder_Reverberators.html (fetched)
- **R16** Moorer, J. A., "About this reverberation business", Computer Music Journal 3(2), 1979. https://freeverb3-vst.sourceforge.io/doc/aboutThisReverberationBusiness.pdf (cited)
- **R17** Gerzon, M. A., "Unitary (energy-preserving) multichannel networks with feedback", Electronics Letters 12(11), 1976 (cited)
- **R18** Stautner, J. & Puckette, M., "Designing multi-channel reverberators", CMJ 6(1), 1982 (cited)
- **R19** Jot, J.-M. & Chaigne, A., "Digital delay networks for designing artificial reverberators", AES 90th Conv., preprint 3030, 1991. https://aes2.org/publications/elibrary-page/?id=5663 (cited)
- **R20** Jot, J.-M., "An analysis/synthesis approach to real-time artificial reverberation", ICASSP 1992 (cited)
- **R22** Valhalla DSP, "Stability through time variation: Ursa Major Space Station". https://valhalladsp.com/2010/05/14/stability-through-time-variation-ursa-major-space-station/ (cited, search summary)
- **R23** LiquidSonics, "The AMS RMX16". https://www.liquidsonics.com/?p=40684 (fetched); Vintage Digital RMX16 article https://www.vintagedigital.com.au/?p=185690 (cited)
- **R24** Smith, J. O., *Physical Audio Signal Processing*, "Mode Density Requirement". https://ccrma.stanford.edu/~jos/pasp/Mode_Density_Requirement.html (fetched)
- **R25** Smith, J. O., PASP, "Choice of Delay Lengths". https://ccrma.stanford.edu/~jos/pasp/Choice_Delay_Lengths.html (fetched)
- **R26** Schlecht, S. J. & Habets, E. A. P., "Modal decomposition of feedback delay networks", IEEE TSP 67(20):5340–5351, 2019. https://arxiv.org/abs/1901.08865 (cited)
- **R27** Heldmann, J. & Schlecht, S. J., "The role of modal excitation in colorless reverberation", DAFx20in21. https://dafx2020.mdw.ac.at/proceedings/papers/DAFx20in21_paper_17.pdf (cited, abstract)
- **R28** Smith, J. O., PASP, "First-Order Delay-Filter Design". https://ccrma.stanford.edu/~jos/pasp/First_Order_Delay_Filter_Design.html (fetched)
- **R29** Coggin, J. & Pirkle, W., "Automatic design of feedback delay network reverb parameters for impulse response matching", AES 141st Conv., 2016. https://aes2.org/publications/elibrary-page/?id=18470 (cited)
- **R30** Ibnyahya, I. & Reiss, J. D., "A method for matching room impulse responses with feedback delay networks", AES 153rd Conv., 2022. https://www.eecs.qmul.ac.uk/~josh/documents/2022/21917.pdf (cited)
- **R31** Lexicon 480L Owner's Manual (1993), §3 "About the Reverberation Algorithm". https://freeverb3-vst.sourceforge.io/doc/Lexicon/lexicon_480l.pdf (fetched). R31b: service manuals at https://www.lost-manuals.com/equipment.php?id=16685 (cited)
- **R32** Shen, J., "Real-time audio reverberation for virtual room acoustics" (MS thesis, advisor R. Duraiswami, UMD, 2020). https://drum.lib.umd.edu/handle/1903/26425 (cited)
- **R33** Dal Santo, G., Prawda, K., Schlecht, S. J. & Välimäki, V., "Differentiable feedback delay network for colorless reverberation", DAFx23. https://dafx.de/paper-archive/details/bwyPN0CBpFcNVpD-ciAPcA (cited)
- **R34** Dal Santo, G. et al., "Optimizing tiny colorless feedback delay networks" (arXiv 2402.11216, v3 2025). https://arxiv.org/abs/2402.11216 (cited)
- **R35** Mezza, A. I. et al., "Data-driven room acoustic modeling via differentiable feedback delay networks with learnable delay lines", EURASIP JASMP 2024. https://asmp-eurasipjournals.springeropen.com/counter/pdf/10.1186/s13636-024-00371-5.pdf (cited)
- **R36** Karjalainen, M., Esquef, P., Antsalo, P., Mäkivirta, A. & Välimäki, V., "Frequency-zooming ARMA modeling of resonant and reverberant systems", JAES 50(12), 2002 (cited from memory; verify)
- **R37** Karjalainen, M. et al., "Estimation of modal decay parameters from noisy response measurements", JAES 2002. https://aes2.org/publications/elibrary-page/?id=11059 (cited)
- **R38** Abel, J. S. & Huang, P., "A simple, robust measure of reverberation echo density", AES 121st Conv., 2006 (cited)
- **R39** Välimäki, V., Parker, J. D., Savioja, L., Smith, J. O. & Abel, J. S., "Fifty years of artificial reverberation", IEEE TASLP 20(5):1421–1448, 2012, doi:10.1109/TASL.2012.2189567 (cited)
- **R40** Välimäki, V. et al., "More than 50 years of artificial reverberation", AES 60th Int. Conf., 2016. https://aes2.org/publications/elibrary-page/?id=18061 (cited)
- **R41** Yamaha MU80 manual (Gate/Reverse Gate descriptions). https://de.yamaha.com/files/download/other_assets/5/317065/MU80F2.pdf (cited, search snippet)
- **R42** Kuttruff, H., "Künstlicher Nachhall", Frequenz 16(3):91–96, 1962. https://www.degruyterbrill.com/journal/key/freq/16/3/html (cited; tapped-loop idea per R1)
- **R43** Texas Instruments, TMS320C1x User's Guide / TMS320C10 data sheet (instruction timings; to verify). Overview table: https://www.komponenten.es.aau.dk/fileadmin/komponenten/Data_Sheet/Microprossor/TMS320C15.pdf (cited)
