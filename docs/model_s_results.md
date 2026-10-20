# Program S grey-box model (Path B): first results on the CLA IR

The fitted model is `tools/model_s.py`. Every number below comes from
`.venv/bin/python tools/model_s.py` (log: `analysis/model_s_run.log`; data:
`analysis/model_s/summary.json`; plot: `analysis/model_s/model_s.png`). The exceptions are numbers
marked *[exploratory]*, which come from one-off scratch runs while the model was being developed.
Labels: **[M]** measured here, **[A]** assumption, **[H]** hypothesis.

![model vs measured](../analysis/model_s/model_s.png)

> **Current status (after round 13):** §1–7 describe the first model (B). The current state is in §8–13:
> - best nulls: cluster −50/−52 dB; held-out 400–600 ms about −28 dB; R tail predicted within 4–6 dB of noise after 0.45 s (§10.1);
> - diffuser found (§11);
> - group-sparse and long-line taps are negative (§12);
> - time variance: R is LTI within about 3 %, L fits about 1 % read-out jitter, and the late input is withdrawn (§13).
>
> Open items are listed in `NEXT_STEPS.md`.

## 1. Summary

- **What fits.**
  - With 4 uncoupled loops, length-scaled gains and one shared first-order HF shelf, the decay matches: T30 is within 0.014 s in every octave from 63 Hz to 8 kHz, including the HF roll-off.
  - Echo-density curves overlap.
  - The sample-domain null is about −30 dB from 50 to 1000 ms in both channels (model B).
- **What doesn't fit.**
  - The late-tail null stops about 15 dB above R's noise floor, and the residual is still structured. Adding a 5th recursion at 460 samples (model C) gains 2–8 dB.
  - The cluster null is limited to −28/−32 dB by how the analysis chain's band-limit kernel is modelled.
  - The loops alone cannot produce the first pass (−19/−16 dB at 9–50 ms). That needs a non-recirculating first-pass FIR.
  - The output-distributor hypothesis fails in its strict form.
  - **Tightening (§8, §9):**
    - In the 20 Hz–11.3 kHz comparison band, the in-sample null reaches −48 to −50 dB at 50–300 ms. The capture-domain cluster reaches −50 to −52 dB.
    - The *held-out* null at 400–600 ms is only −27/−30 dB, about 21–24 dB above noise. The long-tail metrics (EDC, decay rates, fine structure) agree with measurement to within their noise.
    - Loop gains are resolved to about ±8 Q15 LSB, not to the LSB.
  - **Injection (§10):**
    - R's input ends at about 0.45 s. After that, 4 loops + the 460 stage predict R to within 4–6 dB of noise using only the IR's own past.
    - The generative models fail at 400–600 ms because nothing yet generates the dense, per-loop input arriving up to 0.45 s.
    - Not found: a long loop, extra 460-type stages, or one input shared across loops.
  - **Diffuser (§11):**
    - A common all-pass chain AP20 · AP72 · AP230 with g ≈ −0.70, plus the 460 comb (g 0.445), was found on R and confirmed on L. The delays are sharp to ±1 sample.
    - It improves sparse-tap models by 10–13 dB at every tap budget.
    - It does not generate the late input, so the held-out null with a generative model is still about −29 dB in band.
  - **Group-sparse taps (§12):**
    - No small tap set is shared across loops; within a loop, sharing between L and R costs about 4 dB.
    - Long-delay-line (EMT 250-style) taps and loop feedback don't generate the late input.
    - L's excess is in the recirculating read-out, not the injection.
  - **Time-variance through the sweep (§13):**
    - An LTI simulation of the fitted structure reproduces the "late input". It is the diffuser tail plus deconvolution band edges, not an extra generator.
    - R is time-invariant to within about 3 % injection jitter.
    - L's excess fits an L-only read-out gain variation of about 1 % every 0.2 s (unit or capture chain).
    - The non-sparse injection is structural mismatch, not stochastic.
  - What is left over in R sits at loop-mode frequencies; what is left over in L is mostly broadband noise. L's *model* error is not worse than R's (§6.5).
- **HF.**
  - The HF *decay* is right: T30 at 8 kHz is 0.874/0.868 s against 0.878/0.872 s measured. The loop filter's 7 kHz loss is 0.37–0.40 dB per pass. The panel's numeric displays read HF ×0.9, and the measured RT(7 kHz)/RT(1 kHz) ≈ 0.90 is exactly that setting, with the HF factor applied at 7 kHz as the data sheet says. A pure one-pole lowpass with p ≈ 0.037–0.044 fits it as well as a general first-order section. Whether any damping remains at HF ×1.0 is unknown.
  - The HF *fine structure* in the late tail is 1–4 dB worse than mid band (§6.3).
- **The "+0.2 sample" loop offset is a clock offset, not a filter.**
  - It grows with Pᵢ at about 100 ppm, and the cluster taps drift by the same ratio.
  - Once the clock is accounted for, the loop filters' own delay is ≤ ±0.02 samples.
- **Caveat [A].** The excitations and first-pass FIR (about 16.7 k free samples in model B) are fitted *in-sample* to this one IR. Only 5 numbers (RT₀, RT_π, pole, τ₀, ε) set the loop recursion, and those were identified separately (§3). The null depths show the structure *can* represent the IR; they are not a held-out prediction. The real Path-B test is to keep the excitations, change only the loop gains, and predict another RT setting. That needs the library.

## 2. Set-up and assumptions

| Item | Value |
|---|---|
| Data | `analysis/shortRT/h_dev.npy`, regenerated with `tools/analyze_ir.py` on `irs/eugtone_EMT246_CLA_shortRT.wav`; fs_dev = 19.2 MHz/728 = 26,373.6 Hz; L's ½-sample offset already removed by `deconv_grid.recover` |
| Time origin **[A]** | first cluster tap (index 263); all windows are in ms after it |
| Fit range | 0–1.0 s. Weight: decay-flattening 10^(3t/0.98), capped at t = 0.3 s **[A, a trade-off, §6.4]** |
| Band limit **[A]** | the IR only holds information up to 12.3 kHz (the deconvolution window rolls off from 11 kHz). The cluster taps are rendered through that window's analytic device-rate kernel. Everything else is fitted inside the band |
| Settings **[A]** | from a screenshot of CLA's unit showing the numeric displays: SHORT RT, LP 12 kHz, DELAY 10 ms, AMPL 80 % (the display has dead segments, so 60 % is possible; the bar favours 80), **LF ×1.0, RT 1.0 s, HF ×0.9**. This supersedes an earlier reading from the LED bars alone. DELAY 10 ms puts the first cluster tap at the pre-delay itself (index 263 ≈ 9.97 ms from file start; 10 ms = 263.7 device samples), so eugtone appears to have kept absolute timing and the reverb's onset times are meaningful |
| Noise floor **[A]** | from 1.6–1.9 s, assumed stationary. The "noise-limited null" is noise energy over signal energy per window: L is about 9 dB noisier than R |

## 3. Loops: lengths, clock and damping filters

**Method.** Each loop's per-pass response Gᵢ(ω) is identified *without* the excitations. The annihilator e = ∏ᵢ (1 − Gᵢ(z) z^−Pᵢ) · h is zero, up to noise, once every loop's excitation has ended. The fit window is 0.25–0.75 s.

Two estimators are used:
- a non-parametric 13-tap FIR per loop, by alternating least squares;
- a parametric fit by nonlinear least squares.

### 3.1 The +0.2 offset: a clock ratio [M]

| Estimate | Value |
|---|---|
| Cluster taps: sub-sample drift with lag (fit of the band-limit kernel; e.g. the R tap at 353 sits at 353.035) | **104.5 ppm** |
| Loop LF offset τᵢ (100–500 Hz, FIR) for P = 1321/1389/1461/1587 | 0.152 / 0.157 / 0.176 / 0.180 → slope **113 ± 33 ppm**, intercept 0.003 |
| Parametric, clock only (τ₀ = 0) | ε = **98.5 ppm**, cost 0.5448 |
| Parametric, clock + constant offset | ε = 103.3 ppm, τ₀ = −0.007, cost 0.5447 |
| Parametric, constant offset only (ε = 0) | τ₀ = 0.140, cost **0.6192 (+14 %)** |

- **Reading:**
  - The offset scales with Pᵢ, and the cluster shows the same ratio on its own. The IR-library harness independently measured −99 ppm for this capture relative to 19.2 MHz/728, from cluster-tap spacing; the same ratio expressed with the opposite sign. So the capture and the unit's clock differ by about 100 ppm. The analyser's 0.05-sample grid had rounded this to "+0.2".
  - Two explanations fit equally well. Either the interface clock ran fast, or the unit's 19.2 MHz crystal runs about 100 ppm slow, which would put the device rate at about 26,371 Hz **[H; these two cannot be told apart from audio]**.
  - The constant-delay explanation proposed earlier (algorithm_research §4) is rejected.
  - **For the plugin:** use integer Pᵢ, with no fractional delay in the loops.
- **How the model handles it:** to null against this capture, the model renders the loops at Pᵢ(1+ε) and the cluster taps at t(1+ε).

### 3.2 Gains and HF damping [M]

Per-loop RT computed from |Gᵢ| (FIR estimate):

| Hz | 100 | 1k | 2k | 4k | 6k | 7k | 8k | 10k | 11k |
|---|---|---|---|---|---|---|---|---|---|
| 1321 | 0.991 | 0.987 | 0.977 | 0.949 | 0.911 | 0.883 | 0.854 | 0.848 | 0.827 |
| 1389 | 0.998 | 0.991 | 0.977 | 0.943 | 0.911 | 0.891 | 0.865 | 0.844 | 0.818 |
| 1461 | 0.993 | 0.989 | 0.981 | 0.954 | 0.907 | 0.886 | 0.863 | 0.847 | 0.831 |
| 1587 | 0.985 | 0.981 | 0.975 | 0.949 | 0.904 | 0.889 | 0.869 | 0.845 | 0.828 |

- **Length-scaled gains:** the DC loss per pass is 3.034 / 3.168 / 3.346 / 3.667 dB, i.e. 2.30 / 2.28 / 2.29 / 2.31 mdB per sample. That is gᵢ = 10^(−3Pᵢ/(fs·RT)) within ±0.6 %, at RT₀ = 0.991 s.
  - *[exploratory]* Freeing a gain trim per loop changed the cost by < 0.1 % (trims ≤ 0.08 %).
- **HF:** the extra loss at 7 kHz is 0.37 / 0.38 / 0.40 / 0.40 dB per pass, close to the ≈ 0.35 dB in the updated algorithm_research §4. That is the HF ×0.9 setting **[A: from the panel]**: RT(7 kHz)/RT(1 kHz) ≈ 0.90 in every loop. With LF ×1.0 the filter is flat at LF, which matches the 100 Hz values. So at this setting the loop filter needs HF shaping only.
  - Removing the HF damping doubles the cost (1.187 against 0.545).
- **Parametric model (used in B):**
  - One first-order section per loop, (b₀ + b₁z⁻¹)/(1 − p z⁻¹). Its DC and Nyquist gains follow the length-scaling rule from RT₀ = 0.991 s and RT_π = 0.818 s; the pole p = 0.030 is shared.
  - RT_π is an extrapolation, because the data stop at 12.3 kHz.
  - The pole is poorly determined: a zero-pole section (p = 0) gives cost 0.545 as well. First order is enough on this IR.
- **Pure one-pole lowpass** H = R₀(1−p)/(1−p z⁻¹) (fitted with τ₀ = 0, so ε absorbs everything):
  - **Length-scaled one-pole** (pole per loop set by its DC and Nyquist gains): cost **0.5449**, the same as the general section's 0.5447. RT₀ = 0.992 s, RT_π = 0.820 s, ε = 98.6 ppm. The implied poles are p = 0.0365 / 0.0384 / 0.0403 / 0.0438, with extra loss at 7 kHz of 0.36 / 0.38 / 0.40 / 0.43 dB per pass and DC group delay of 0.038–0.046 samples.
  - **One fixed one-pole in every loop:** p = 0.0393, cost **0.5491**, i.e. +0.8 %.
  - **Reading:** a mild one-pole (p ≈ 0.04) is enough for HF ×0.9. The data slightly favour scaling the pole with loop length over one shared coefficient, but 0.8 % is not decisive on one IR.
- **Shared coefficients vs length-scaled: not decided.**
  - Cost is 0.5467 against 0.5447 with free ε, and 0.5702 against 0.5673 with ε fixed at the cluster value. That is a 0.4–0.5 % difference.
  - At HF ×0.9 the two hypotheses differ by ≤ 0.015 s in RT. An HF ×0.1 setting will separate them, and HF ×1.0 will show whether any fixed damping remains.
- **Filter phase:** after removing εPᵢ, the per-loop delay varies by about +0.015 samples (2–3 kHz) and −0.01 (9–11 kHz). That is consistent with a mild first-order section **[H]**.

### 3.3 A fifth recursion at 460 samples [M, interpretation H]

- The annihilator residual of the 4-loop model is not white. *[exploratory]* A scan over single lags 150–3500 found that one extra recursion at lag **460** removes 60 % (L) and 65 % (R) of the residual. The lag-scan peaks also include 230, 920 and 690.
- Fitting 1 − G₅ z^−460 together with the loops:
  - cost 0.5447 → **0.1945**;
  - g₅ = 0.439 per pass (−7.15 dB, i.e. RT ≈ 0.146 s), HF ratio 0.89;
  - the loops refit to RT₀ = 1.000 s.
- **Where it sits:**
  - In the full model, a parallel 460 loop (C) and a 460 recursion in series after the loop sum and the first-pass FIR (Cs) give the same null within 0.5 dB. Its position can't be determined from this IR.
  - A plausible reading is a 17.4 ms all-pass or comb stage with |g| ≈ 0.44 in the diffusion path **[H]**. It is not in the §3.4 topology of the discovery doc.
- **Test on the library:**
  - Is the 460 stage present in M and L?
  - Does g₅ stay fixed when RT changes? Fixed would mean a diffuser.

## 4. Injection

- **Model B** (per channel):
  - 4 excitations xᵢ of Pᵢ + 300 samples, starting 230 samples after the first tap;
  - a first-pass FIR of 1,400 samples (230–1,630, i.e. until just after the shortest loop's first recirculation at 1,551);
  - 8 cluster gains.
  - All of these are solved jointly by LSQR with the loop filters fixed.
- **Reverb onset:** about 240 samples after the first tap in L, and about 570 in R (330 later; compare the 332 L/R lag).
- **The loops cannot carry the first pass by themselves.**
  - Model A has no first-pass FIR. Its null is −19.4/−15.6 dB at 9–50 ms, and fitting that early part costs the tail 8–11 dB.
  - In B, the first-pass FIR carries −10.7 dB (L) and −7.1 dB (R) of the first 1,400 samples' energy. That is content that does not recirculate the way the loop excitations do.
  - This fits injection or read-out at several points along each delay line, where some taps see an input before it wraps **[H]**. It is evidence against a single write and single read per loop.
- **The spec's literal method fails.**
  - Model P takes one period per loop from periodic decomposition of the flattened tail, scales it per loop and channel, and uses the measured first pass as the FIR. It nulls only −7.4/−7.7 dB (50–300 ms) and −9.3/−9.6 dB (300–1000 ms).
  - *[exploratory]* The mapping is correct: on an undamped synthetic it recovers the excitations with correlation ≥ 0.999.
  - On a noise-free synthetic with the fitted damping it again gives only −7 to −10 dB. Periodic decomposition averages over passes, while the HF loss and phase change from pass to pass.
  - **The excitations must be fitted in the time domain with the loop filters**, which is what `Linear` does.

## 5. Output: common loop sum + distributor vs per-channel components

| Test | Result |
|---|---|
| D: common excitations, L = d_L ∗ s, R = d_R ∗ s, with 701-tap d_c (alternating LS, 5 rounds) | null **−10.8/−10.7 dB** (50–300), −8.9/−7.9 (300–1000), against −31/−30 and −29/−29 for B. The ALS cost was still falling (0.73× the start), but D would need about 100× less cost to compete. d_R stayed near its start values (lags 332: 0.55 and 26: 0.53) |
| L/R harmonic-ratio test on B's excitations. A common distributor implies X_iL/X_iR(f) is the same function for every loop | harmonics of *different* loops within 1 Hz: phase score 0.60 (0 = equal, random pairs 1.01); the magnitude ratio differs by a median 4.85 dB (random 4.82 dB) |
| Per-loop sparse L→R relation | lags 26 and 332 explain 47 / 49 / 53 / 54 % of R's periodic component in loops 1321 / 1389 / 1461 / 1587, with **loop-specific** gains (0.45/0.59, 0.55/0.49, 0.34/0.73, 0.67/0.33). *[exploratory]* Greedy OMP over all lags picked 26 and 332 first in every loop. Four more taps add only 7–9 % |

- **Reading:**
  - The strict "loop sum → one fixed distributor" is **rejected**. L/R phases are partly shared across loops, but magnitudes are not.
  - The same two lags (26 and 332) appear in every loop, with different weights per loop. So the output side shares its *delays* across loops but not its *weights*. Examples: two output buses that mix the loops differently, or read taps at fixed offsets from each loop's write point **[H]**.
- **For now:** keep independent per-channel components (B).

## 6. Null test (model − measured) [M]

### 6.1 Null depth, dB (L / R)

| Model | 0–50 ms | 50–300 ms | 300–1000 ms | 0–9 ms (cluster only) | 9–50 ms | 300–600 | 600–1000 |
|---|---|---|---|---|---|---|---|
| A: loops only | −22.6 / −19.8 | −23.6 / −21.3 | −21.3 / −18.7 | −28.3 / −32.5 | −19.4 / −15.6 | −21.7 / −18.9 | −12.3 / −11.6 |
| **B: loops + first-pass FIR (main)** | **−30.5 / −34.5** | **−31.3 / −30.1** | **−29.1 / −29.3** | −28.3 / −32.5 | −68.0 / −66.4 ¹ | −30.6 / −29.7 | −15.1 / −20.2 |
| C: B + parallel 460 loop | −30.5 / −34.5 | −35.8 / −37.7 | −30.9 / −34.4 | −28.3 / −32.5 | −62.2 / −60.0 ¹ | −33.1 / −34.9 | −15.7 / −24.3 |
| Cs: B + 460 in series | −30.4 / −34.4 | −35.7 / −37.2 | −30.8 / −34.2 | −28.3 / −32.5 | −55.3 / −51.0 ¹ | −33.0 / −34.8 | −15.6 / −23.7 |
| D: common sum + distributor | −30.5 / −34.5 | −10.8 / −10.7 | −8.9 / −7.9 | −28.3 / −32.5 | −56.9 / −58.1 ¹ | −9.0 / −8.0 | −7.8 / −6.9 |
| P: periodic decomposition | −30.5 / −34.5 | −7.4 / −7.7 | −9.3 / −9.6 | −28.3 / −32.5 | (copied) ¹ | −9.4 / −9.6 | −9.1 / −9.8 |
| *noise-limited floor* | −66.6 / −74.7 | −56.2 / −64.8 | −35.6 / −44.3 | −71.8 / −80.1 | −63.4 / −71.2 | −39.2 / −48.0 | −18.9 / −27.6 |

¹ The first-pass FIR is free over 9–62 ms, so this column is a fit **by construction**, not a prediction. P copies the measured first pass exactly.

- **0–50 ms** is dominated by the cluster, which limits it to −28/−32 dB. That residual sits *between* the taps, in the ±0.12 band-limit ripple: within ±3 samples of the taps it is −39/−42 dB.
  - *[exploratory]* An empirically fitted kernel did not improve it.
  - This is a limit of modelling the analysis chain, not of the cluster model **[A]**: on the hardware the taps are single samples.
- **Cluster gains** (kernel units): L 1.128, 1.128, 1.119, −1.140; R 1.120, 1.116, −1.125, 1.104. All are equal within ±1.6 %, with signs + + + − / + + − + as expected.
- **300–1000 ms:**
  - L's 600–1000 ms is noise-limited (−15.1 against a floor of −18.9).
  - R's is about 7 dB above its floor (−20.2 against −27.6); its 300–600 ms is 18 dB above (−29.7 against −48.0). So R shows the model error remaining in the tail.
  - *[exploratory]* The residual's autocorrelation still peaks at the loop lengths and their doubles (0.2–0.3), and at 460/920. Part of the loops' periodic content is still not captured.

### 6.2 T30 per octave (s), measured → model B (L | R)

| Octave | 63 | 125 | 250 | 500 | 1k | 2k | 4k | 8k |
|---|---|---|---|---|---|---|---|---|
| L meas | 0.994 | 1.004 | 0.987 | 0.976 | 0.978 | 1.000 | 0.955 | 0.878 |
| L model B | 0.987 | 0.993 | 0.981 | 0.973 | 0.974 | 0.994 | 0.951 | 0.874 |
| R meas | 1.028 | 0.994 | 1.003 | 1.008 | 0.996 | 1.003 | 0.950 | 0.872 |
| R model B | 1.014 | 0.985 | 0.995 | 1.007 | 0.992 | 0.997 | 0.947 | 0.868 |

- Every model B octave is within −0.014 to −0.001 s, slightly short.
- Model C is within ±0.002 s in every octave of both channels: e.g. L 0.993, 1.002, 0.988, 0.977, 0.977, 1.000, 0.954, 0.876.
- The L/R differences in the measured T30 (e.g. 500 Hz: 0.976 against 1.008) come from modal beating, because the loops are common to both channels. The model reproduces them because the excitations are fitted.
- Method: T30 is computed at the device rate with `analyze_ir.edc`; the 8 kHz band is capped at 12.3 kHz.

### 6.3 Residual per band, model B (null dB, L / R)

| Band | 63 | 125 | 250 | 500 | 1k | 2k | 4k | 8k | 10–12.3k |
|---|---|---|---|---|---|---|---|---|---|
| 50–300 ms | −33.0/−30.0 | −33.7/−40.3 | −32.9/−41.3 | −30.9/−26.9 | −32.0/−29.8 | −32.9/−31.6 | −33.9/−31.8 | −35.0/−32.3 | −36.2/−31.6 |
| 300–1000 ms | −33.6/−27.3 | −31.7/−36.4 | −33.6/−40.3 | −30.6/−27.7 | −33.2/−31.2 | −33.5/−32.0 | −33.2/−31.1 | −31.1/−29.5 | −29.8/−27.7 |

- **HF:** in 50–300 ms, HF nulls as well as mid band or better. In 300–1000 ms, 8 kHz and 10–12.3 kHz are 1–4 dB worse than 1–4 kHz. That matches small HF loss or phase errors that accumulate over passes, from the first-order shelf.
  - Closing that gap needs an HF/LF sweep; one IR doesn't constrain the filter shape further. *[exploratory]* A 33-tap free FIR per loop lowered the annihilator cost by only about 1 % against the parametric model (0.538 against 0.545).
- **R at 500 Hz and R at 63 Hz** are the worst bands (about −27 dB). The reason is unknown.

### 6.4 Weighting trade-off *[exploratory]*

- With excitations starting at 0 rather than 230: capping the flattening at 0.1 s gives 50–300 ms −35.0/−34.7 and 300–600 ms −24.9/−23.9.
- Capping at 0.3 s (as used) gives −31.3/−30.2 and −30.6/−29.7.
- With the B configuration, extending the first-pass FIR to 2,500 samples improved the null to −36.8/−38.9 (50–300) and −32.8/−34.5 (300–600). It then covers part of the second pass, so it was not used.

### 6.5 Is the residual at loop modes or broadband? L vs R [M]

This bears on the 1986 "stochastic control" claim (algorithm_research §7.2b). Tests: `residual_mode_test`, on the decay-flattened window.

- **on/off:** spectral density within 0.3 Hz of any loop harmonic, over density more than 1.2 Hz from every harmonic. These cover 12 % and 57 % of the 100 Hz–11 kHz bins. 0 dB means broadband; the measured IR itself gives the fully modal value.
- **Mode-shaped share:** in each 250 Hz band, fit S_res = a·S_meas + b (a, b ≥ 0). The share is the fraction of residual energy that follows the IR's modal fine structure.
- **Control:** the measurement's own noise tail (1.4–1.9 s).

| Signal | Window | on/off L / R (dB) | mode-shaped share L / R |
|---|---|---|---|
| measured IR | 0.1–0.6 s | 11.7 / 11.5 | (reference) |
| residual, B | 0.1–0.6 s | 3.9 / 4.8 | 0.37 / 0.42 |
| residual, B | 0.3–0.8 s | 8.1 / 11.4 | 0.61 / 0.79 |
| residual, C (+460) | 0.1–0.6 s | 0.7 / 4.3 | 0.08 / 0.33 |
| residual, C (+460) | 0.3–0.8 s | 2.6 / 9.5 | 0.17 / 0.66 |
| noise tail (control) | 1.4–1.9 s | 0.4 / 1.6 | 0.06 / 0.13 (the test's false-positive level) |

- **R:** the unexplained residual is **concentrated at loop-mode frequencies**, more so later in the tail. With C, 66 % of R's 0.3–0.8 s residual follows the modal structure, with on/off 9.5 dB against 11.7 for the IR itself. So the leftover error is mainly in the modes' amplitudes and decays, not broadband.
  - This fits the "randomised read-out gains" reading (§7.2b, reading 1).
  - It equally fits an ordinary static mis-fit: a slightly wrong loop-filter shape, which also leaves error at mode frequencies that grows with time.
  - **One IR cannot tell these apart.** §8.6 tests one consequence within this IR: the line widths match an LTI model. Library IRs at other settings can test another: whether the mode-level error moves with RT, which a randomised read-out would not.
- **L:** the residual is **mostly broadband**: 8–17 % mode-shaped with C, close to the control. L carries more measurement noise:
  - L's noise-limited floor is 9 dB above R's.
  - L's late floor is broadband (control row) but **not stationary**. It falls from −88.1 dB at 1.1 s to −90.9 dB at 1.8 s (about 4 dB/s), while R's is flat at −98.4 to −98.9 dB from 1.4 s.
  - A slowly decaying broadband floor is what time-variance smear in a sweep measurement can leave. It could equally come from the capture chain **[H]**.
- **Is L's model error worse than R's? No.** Subtracting the noise floor's energy from the null (arithmetic on the §6.1 numbers):
  - Model B at 300–600 ms: model error is −31.2 dB in L and −29.8 dB in R.
  - Model C: −34.3 dB in L and −35.1 dB in R.
  - So the *model* error is about the same in both channels. **L's extra residual is broadband noise**, not extra mode-level error.

### 6.6 Echo density

The NED curves of models B and C overlap the measured ones over 0–300 ms (plot, middle right). The first crossing of 0.9 falls at:

| | measured | model B | model C |
|---|---|---|---|
| L | 75.9 ms | 75.9 ms | 75.9 ms |
| R | 80.9 ms | 90.8 ms | 80.9 ms |

The crossing time is a noisy statistic, so compare the curves. NED is computed at the device rate (20 ms window) and is not comparable with the 96 kHz figure in the discovery doc.

## 7. What this means for Path B and the library work

1. **Use as-is:**
   - integer loops with length-scaled gains, from RT₀;
   - one first-order HF section per loop;
   - per-channel excitation FIRs (Pᵢ + 300 taps) and a first-pass FIR (1,400 taps per channel);
   - 8 fixed cluster taps.

   Parameters are in `analysis/model_s/model_B_params.npz`.
2. **Held-out test (decisive):** render B with RT changed, keeping the excitations, and null it against the library's RT steps. If the null stays near −30 dB, Path B is viable. If not, the excitations depend on RT, which would mean exponential tap weighting (T5).
3. **HF/LF sweeps:** decide length-scaled against shared coefficients (§3.2), and the filter order and shape (§6.3).
4. **The 460 recursion:** is it in M and L, and does g₅ change with RT? Model C or Cs is the better candidate structure if it holds.
5. **Output stage:** model per-loop L/R weights on shared lags 26 and 332, which would cut the per-channel parameters roughly in half if it holds (§5). A strict common distributor is ruled out.
6. **Clock:** re-grid captures with fs_dev/(1+ε). Expect about 100 ppm on this interface; check whether ε differs between libraries. If it does, the drift is in the capture chain.

## 8. Tightening against one IR: hold-out, sparsity, structure, capture domain

Code: `.venv/bin/python tools/model_s.py --tighten` → `analysis/model_s/tighten.json` (log: `analysis/model_s_tighten.log`). Everything below is from that run unless marked otherwise.

- **Settings [A]:** this section uses the confirmed panel settings (§2): LF ×1.0, RT 1.0 s, HF ×0.9, DELAY 10 ms.
- **Absolute timing:** with DELAY 10 ms, cluster tap 0 at index 263 (9.97 ms from file start) is the pre-delay output. So the reverb onsets are absolute: about 240 device samples (9.1 ms) after tap 0 in L and about 570 (21.6 ms) in R.

### 8.1 Ground rules

- **Hold-out:** excitations, first-pass FIR and cluster gains are fitted on **0–400 ms only**. The loop and 460 parameters are also identified from data before 400 ms (annihilator window 0.25–0.40 s). Windows after 400 ms are predictions.
- **Comparison band 20 Hz–11.3 kHz** (`band_limit`).
  - Below 20 Hz the IR holds a large non-reverb component (sub-20 Hz energy is −22 dB of the signal at 50–300 ms), which is presumably sweep/deconvolution leakage **[A]**. In broadband nulls it dominates the residual: 0–20 Hz carries −36/−38 dB of residual against −44.5/−42.1 dB at 44 Hz–11.3 kHz (§9.3).
  - Above 11.3 kHz the deconvolution window rolls off.
  - Both the broadband and the in-band figures are kept.
- **Excess over noise** = null − noise-limited floor of the same window and band. 0 dB is the target.
  - The in-band floors are much lower than the broadband ones. The broadband "noise floor" is itself mostly below 20 Hz: R's in-band floor is −110 dB against −98 dB broadband.

### 8.2 Hold-out table, device domain, 20 Hz–11.3 kHz (null dB, L / R; excess over noise in brackets)

| Model | Params | 0–50 ms (fit) | 50–300 (fit) | 300–400 (fit) | **400–600 (predicted)** | **600–1000 (predicted)** |
|---|---|---|---|---|---|---|
| B: 4 loops + first-pass FIR | 16,729 | −31.2/−35.2 (+48.8/+50.7) | −35.1/−33.0 (+34.5/+43.2) | −32.6/−30.7 (+23.8/+32.3) | −20.6/−19.7 (+27.3/+34.7) | −12.8/−12.9 (+19.5/+26.1) |
| C: + parallel 460 | 18,252 | −31.2/−35.2 (+48.8/+50.7) | **−48.1/−50.1** (+21.4/+26.0) | −43.5/−46.1 (+12.9/+16.9) | −26.7/−30.1 (+21.2/+24.3) | −17.2/−21.3 (+15.1/+17.6) |
| Cs: + series 460 | 16,732 | −31.2/−35.2 | −47.4/−46.8 (+22.2/+29.3) | −42.6/−43.4 | **−28.4/−27.5** (+19.5/+27.0) | −19.5/−18.5 (+12.8/+20.4) |
| Cs without first-pass FIR | 13,932 | −27.6/−22.3 | −28.9/−22.7 | −26.8/−20.6 | −17.8/−12.1 | −10.4/−6.7 |
| Cs, loops + 460 taken from the **L-only** fit | 16,732 | −31.2/−35.2 | −47.0/−46.1 | −42.5/−42.9 | −28.5/−27.3 | −19.6/−18.5 |
| *noise floor (in band)* | | −80.0/−85.9 | −69.6/−76.1 | −56.4/−63.0 | −47.9/−54.4 | −32.3/−38.9 |

The broadband equivalents are in `tighten.json` (`null`, `excess`). For example, C at 50–300 ms is −36.3/−38.5 broadband.

- **Reading:**
  - The 460 recursion is the one structural addition that clearly helps: +13/+17 dB in-band at 50–300 ms, and +6/+10 dB on the held-out 400–600 ms.
  - At the fit/hold-out boundary the null drops by about 17 dB, from −43/−46 dB (300–400 ms) to −27/−30 dB (400–600 ms). The 16–18 k free excitation samples fit the 0–400 ms window much better than they generalise. **The hold-out numbers are the honest ones.**
- **Cross-channel:** the 9 shared parameters (loops + 460) fitted on L alone predict R as well as R's own parameters.
  - Annihilator cost on R: 0.1048 against 0.1042 (+0.5 %). On L with R's parameters: +0.4 %.
  - R's in-band nulls with L-derived loops differ from the R-derived ones by ≤ 0.7 dB.
  - The early-window loop parameters (RT₀ 0.998 s, ε 98.4 ppm, g₄₆₀ 0.445) match the full-window ones. Without the 460 stage they drift (RT₀ 0.970, ε 113.7).

### 8.3 Sparse integer-tap excitations

Greedy batch OMP picks integer-position excitation and first-pass taps (Cs basis, fit 0–400 ms):

| Taps per channel (+8 cluster) | 103 | 414 | 1,719 |
|---|---|---|---|
| L null 50–300 / 400–600 (broadband) | −1.7 / −1.3 | −4.8 / −3.9 | −12.8 / −8.6 |
| R null 50–300 / 400–600 | −1.8 / −1.5 | −4.9 / −4.0 | −12.4 / −8.9 |

- **The excitations are not a modest set of integer taps.**
  - They are heavy-tailed (kurtosis 11–20; top 1 % of samples hold 25–34 % of the energy), and some positions repeat across loop pairs at offsets of 20, 52 and 72 samples. *[exploratory]*
  - Yet 1,719 taps per channel recover less than the dense model does with 8,400.
  - Something dense, a diffuser we have not identified, sits between a sparse tap set and the loops.
- **The 460 stage as an all-pass:** *[exploratory]* the excitations keep a −0.31 lag-460 regression coefficient after the 460 recursion is modelled. That fits an all-pass numerator (−g + z⁻⁴⁶⁰). But removing that numerator does not make them sparser.
- **Other short recursions:** lags 72, 20, 52, 230 each explain ≤ 0.2 % of the excitation variance *[exploratory]*. So no second short all-pass is visible this way.

### 8.4 What else lowers the annihilator cost (loops + 460, window 0.25–0.75 s, cost L+R = 0.195)

| Added freedom | Cost |
|---|---|
| none (4 loops + 460) | 0.0902 + 0.1043 = 0.1945 |
| + recursions at 230, 690, 920 (gain + HF each) | 0.1825 (−6 %); singly 0.190 / 0.194 / 0.194 |
| per-loop free 81-tap FIR instead of the shelf | 0.0905 + 0.1008 (−2 %) |
| per-loop free shelves (12 parameters) | 0.1944 |
| LF shelf in every loop | 0.1923 (corner 9 Hz, i.e. a DC effect) |
| cross-coupling between loops (12 off-diagonal gains) | 0.1944 (couplings ≤ 0.08, with no gain) |
| *noise-only level of this objective* | 0.0123 + 0.0025 |

- Annihilator residual over its noise-only level, in 50 ms blocks:
  - **R:** +47 dB at 0.25 s, +7 dB at 0.40 s, +1.4 dB at 0.50 s, and ≤ 0 dB from 0.55 s on. **In R, 4 loops + 460 explain the recursion to the noise floor from about 0.5 s.**
  - **L:** +15 dB at 0.40 s, falling to +4 dB at 0.85 s. It never reaches noise, and no LTI change above removes it.
- So two things limit the null, and neither is the loop filters:
  1. **Unmodelled input arriving at 0.25–0.4 s**, in both channels: annihilator residual +47 → +7 dB (R). This is content that neither the loops, the 460 stage, nor the excitation and first-pass supports can carry. A longer free excitation support did not help: *[exploratory]* +1,000 samples moved Cs at 400–600 ms (broadband) from −27.9/−26.8 to −27.8/−26.5 dB. It is probably the unidentified diffuser of §8.3 having a longer tail than the model allows **[H]**.
  2. **An L-only excess in the tail.** Its periodicities are the loop-length *differences* (68, 126, 140, 198, 266 samples) *[exploratory]*, which is what a slow modulation of the loop or read-out gains in L would leave **[H]**.

### 8.5 Capture domain (raw 96 kHz WAV)

- **Method:** the device model uses *integer* Pᵢ, delta cluster taps and the clock ratio ε only in the resampling. It is resampled onto the 96 kHz grid with a Kaiser-sinc interpolator, then passed through the analogue chain. That chain is one shared 1,193-tap FIR at 96 kHz × a 2nd-order high-pass (15 + 15 Hz), with L allowed a fractional offset. The comparison is against the raw WAV low-passed at 11.5 kHz.
- **Cluster-only test** (reverb-free regions, L < 9 ms, R < 21 ms; kernel fitted on R, then **held out on L**):

| Analogue high-pass | R in-sample | L held out | L offset |
|---|---|---|---|
| none | −37.0 | −35.4 | 0.492 device samples |
| 10 Hz | −38.3 | −40.5 | 0.492 |
| 10 + 10 Hz | −46.5 | −41.2 | 0.492 |
| 15 + 15 Hz | −47.5 | −40.6 | 0.492 |
| *noise floor* | −81.5 | −75.6 | |

- **What the cluster shows:**
  - The L offset is **0.492 device samples**: the ½-sample DAC multiplex, fitted freely.
  - The R tap gains are **0.9999 / 1.0010 / −0.9992 / 0.9999** (cluster-only, 15 + 15 Hz), i.e. the taps are unity gain. L's three in-region taps are 1.0065–1.0075, an L/R level offset of about 0.06 dB. L's tap at 288 lies outside the L region and is not determined there.
  - The analogue chain needs a 2nd-order sub-20 Hz high-pass (coupling capacitors **[H]**). A plain FIR cannot represent its long, low-level tail. That tail was the main reason the device-domain cluster stopped at −28/−32 dB.
- **Full model Cs in the capture domain** (kernel fixed from the cluster fit, device parameters solved on 0–400 ms):

| Window | 0–9 ms | 0–50 | 50–300 | 300–400 | 400–600 (pred.) | 600–1000 (pred.) |
|---|---|---|---|---|---|---|
| null, 20 Hz–11.3 kHz | **−49.5/−51.9** | **−51.0/−52.9** | −47.8/−47.7 | −42.7/−44.4 | −27.1/−29.2 | −17.0/−19.4 |
| excess over noise | +36.0/+39.7 | +29.2/+33.2 | +22.0/+28.5 | +13.9/+18.7 | +21.0/+25.3 | +15.6/+19.6 |
| null, broadband (to 11.5 kHz) | −47.7/−50.1 | −49.4/−51.2 | −42.4/−43.2 | −35.6/−38.4 | −26.3/−28.5 | −15.2/−18.7 |

- Refitting the kernel jointly with the reverb changes these by ≤ 1.1 dB (`kernel_joint`).
- Moving to the capture domain gains about 20 dB on 0–50 ms (cluster) and about 6 dB broadband on 50–300 ms, the latter from modelling the sub-20 Hz analogue high-pass. In band, the reverb windows are the same as the device domain.

### 8.6 Line width (the "stochastic control" question within one IR)

Spectral concentration at loop modes, on/off in dB (§6.5 measure), for the measured IR against the LTI model C plus measured noise:
- 0.1–0.4 s: 6.61 / 6.50 measured against 6.61 / 6.50 model.
- 0.4–0.7 s (the model's held-out range): 6.65 / 6.55 against 6.64 / 6.54.

**No line broadening beyond the LTI decay** is visible at this resolution (0.3 s windows, about 3 Hz). Modulation slower than about 1 Hz, or affecting only the L excess of §8.4, is not excluded **[limit]**.

### 8.7 What limits the null now

1. **Unmodelled late input (0.25–0.4 s), both channels.** This is the main limit on the held-out windows (+21 to +25 dB over noise at 400–600 ms). It needs a structural idea, most likely the diffuser of §8.3. It does not need more free taps: those overfit (§8.2).
2. **L-only tail excess**, about +4 to +15 dB over noise in the annihilator, which no time-invariant change removes. Either L-side time variance (the 1986 "stochastic control", if confined to L) or a capture artefact **[H]**.
3. **Analogue-chain model**, for the cluster: −50/−52 dB against a −86/−92 dB floor. The first-order kernel + high-pass is the limit; it is not the digital model.
4. **Sub-20 Hz content** of the capture: broadband nulls only, excluded from the comparison band.
5. **Fixed-point effects** are not reached. They only matter once the model error nears −60 dB, which it does nowhere yet.

So −50 dB is reached only for the cluster and the first 50 ms in the capture domain, and in-sample at 50–300 ms. Nowhere is the null at the noise floor, and **this IR does not yet prove that the true algorithm has been found**.

## 9. Long tail: noise-robust metrics, noise model, loop-parameter precision

Code: `.venv/bin/python tools/model_s.py --tail` → `analysis/model_s/tail.json`.
- **Loop fit:** loops + 460 come from the whole decay (annihilator over 0.25–1.8 s). The residual of a correct model there is filtered noise, so it is fitted unweighted, then iteratively reweighted by the inverse local residual RMS (50 ms blocks, 3 passes). That is the SNR weighting.
- **Model:** model C uses those loops, with the excitations fitted on 0–400 ms. **Everything after 400 ms is a prediction.**

### 9.1 Loop parameters and their precision

| Parameter | Value | 1σ (Jacobian × √τ_int, τ_int = 54) | in Q15 LSB | L-only / R-only | 0.25–0.7 s / 0.7–1.8 s |
|---|---|---|---|---|---|
| g₁₃₂₁ (per-pass DC gain) | 0.707439 | 0.000246 | 8.1 | 0.707138 / 0.707452 | 0.707443 / 0.707425 |
| g₁₃₈₉ | 0.695463 | 0.000231 | 7.6 | 0.695132 / 0.695470 | 0.695486 / 0.695272 |
| g₁₄₆₁ | 0.682337 | 0.000229 | 7.5 | 0.682280 / 0.682338 | 0.682345 / 0.682162 |
| g₁₅₈₇ | 0.660021 | 0.000235 | 7.7 | 0.659779 / 0.660036 | 0.660034 / 0.659852 |
| κ, HF excess at Nyquist (length-scaled) | 12.84 dB/s | 0.10 | | 13.40 / 12.82 | 12.82 / 13.86 |
| shelf pole p | 0.0392 | 0.0049 | 159 | 0.028 / 0.040 | 0.040 / 0.016 |
| clock ε | 98.39 ppm | 0.14 | | 97.86 / 98.41 | 98.41 / 97.12 |
| g₄₆₀ | 0.458 | 0.020 | 656 | 0.397 / 0.501 | 0.531 / 0.016 |

- **Per-loop RT₀ from the gains:** 0.9997 / 1.0017 / 1.0011 / 1.0004 s. That is length-scaled gains exactly at the RT 1.0 s setting, within ±0.2 %.
- **Q15 resolution, the bit-exactness question:**
  - The loop gains are resolved to about **±8 LSB (1σ)**.
  - Early and late halves agree to 0.6–7 LSB, but L-only and R-only differ by 2–11 LSB, systematically lower in L (the L excess again).
  - So this IR is roughly 3–4 bits short of pinning Q15 loop coefficients. The filter pole (±160 LSB) and the 460 gain (±650 LSB; it is unconstrained after 0.7 s, where it has decayed) are far from it.
  - **[A]** The per-pass DC gain equals a hardware coefficient only if the loop filter has unity DC gain. That is not known.
- **What would tighten it:** several library IRs at different RT settings sharing the same loop structure. Each adds an independent decay, and the gain law across RT steps constrains the coefficient table.

### 9.2 Tail metrics, model against measured (both over the predicted range)

- **Noise-compensated EDC per octave**, max |difference| down to −50 dB (L / R):

| Octave | 63 Hz | 125 | 250 | 500 | 1 k | 2 k | 4 k | 8 k |
|---|---|---|---|---|---|---|---|---|
| L / R | 0.19 / 0.68 | 0.59 / 0.17 | 0.15 / 0.15 | 0.08 / 0.13 | 0.07 / 0.07 | 0.05 / 0.05 | 0.03 / 0.03 | 0.06 / 0.05 |

- **Per-band decay rate and frequency** (free 81-tap FIR per loop over the whole decay, against the parametric loops):
  - 125 Hz–8 kHz: decay deviation within −0.51…+0.31 dB/s (on 60 dB/s, i.e. < 1 %); frequency deviation within ±7.4 mHz.
  - 10 kHz: −0.6 to −2.1 dB/s and up to 19.5 mHz. Near the band edge the first-order shelf is slightly off.
- **Spectral fine-structure correlation**, model against measured (ceiling = model against model + a measured-noise realisation):

| Window | 0.3–0.6 s | 0.6–0.9 | 0.9–1.2 | 1.2–1.5 |
|---|---|---|---|---|
| L: model vs measured (ceiling) | 0.998 (1.000) | 0.986 (0.997) | 0.886 (0.942) | 0.51 (0.53) |
| R | 0.998 (1.000) | 0.992 (0.999) | 0.964 (0.972) | 0.66 (0.65) |

- **Reading:** in the tail the *algorithm* matches measurement to within the noise: EDCs within 0.2 dB in most bands to −50 dB; decays within 1 %; mode frequencies within a few mHz; fine structure at the noise ceiling in R. The waveform null degrades because the excitation is not exactly right (§8.7, item 1), not because the loops are wrong. L sits below its ceiling at 0.9–1.2 s, the L excess again.

### 9.3 Noise model and "model error = residual − noise"

- **Noise model** (`noise_model`): an octave-band spectrum from 1.4–1.9 s, and a level that is linear in dB against time (fitted from 1.3/1.4 s).
  - **Broadband:** L −89.0 dB at 1.4 s, falling 4.2 dB/s; R −98.4 dB, falling 1.0 dB/s.
  - **In band (20 Hz–11.3 kHz):** L −103.3 dB falling 1.8 dB/s; R −110.3 dB and flat (+0.4 dB/s).
  - **So L's "decaying floor" is mostly below 20 Hz.** In band, L's floor is 7 dB above R's and drifts only slightly.
  - The noise spectrum rises with frequency in both channels: the 8 kHz octave is −18/−16 dB of the total, the 63 Hz octave −44/−35 dB.
  - **[A]** The model extrapolates the 1.3–1.9 s level back to earlier windows. For L both a drifting and a flat extrapolation are reported.
- **Model error per window, in band, with model C** (dB re signal in the window; "res/noise" = residual over noise model):

| Window | L null | L noise (drift / flat) | **L model error** | R null | R noise | **R model error** | res/noise L (drift) / R |
|---|---|---|---|---|---|---|---|
| 50–300 ms | −43.8 | −66.8 / −69.4 | −43.8 | −42.0 | −76.1 | −42.0 | +23.0 / +34.1 |
| 300–400 | −40.1 | −53.9 / −56.2 | −40.2 | −38.9 | −63.0 | −38.9 | +13.9 / +24.1 |
| 400–600 (pred.) | −25.7 | −45.8 / −47.8 | −25.7 | −26.7 | −54.4 | −26.7 | +20.1 / +27.8 |
| 600–1000 (pred.) | −16.7 | −30.7 / −32.2 | −16.8 | −19.6 | −38.9 | −19.7 | +14.0 / +19.3 |
| 1000–1500 (pred.) | −6.3 | −6.6 / −7.3 | (≈ 0) | −13.2 | −13.1 | (at noise) | +0.3 / −0.1 |

- **Model error, broadband, per octave** (dB re that band's signal, L / R):

| Window | 63 Hz | 125 | 250 | 500 | 1 k | 2 k | 4 k | 8 k |
|---|---|---|---|---|---|---|---|---|
| 50–300 ms | −37.4/−41.2 | −41.7/−46.0 | −45.2/−50.2 | −48.7/−48.6 | −52.7/−50.7 | −52.0/−50.9 | −48.2/−47.0 | −42.3/−39.7 |
| 300–600 | −34.8/−44.4 | −37.7/−41.2 | −37.7/−41.5 | −40.3/−41.2 | −42.1/−40.8 | −41.2/−38.9 | −38.7/−37.2 | −31.8/−30.4 |
| 600–1000 | −40.9/−41.9 | −32.7/−34.4 | −23.1/−26.8 | −28.0/−26.6 | −28.7/−29.2 | −30.4/−27.1 | −26.9/−25.8 | −18.7/−18.4 |

- **Where the residual sits** (50–300 ms, share of total signal energy, L / R):
  - 0–20 Hz: −35.9/−38.0 dB (signal there: −22.2/−23.0);
  - 20–44 Hz: −55.7/−63.1;
  - 44 Hz–11.3 kHz: −44.5/−42.1;
  - 11.3–13.2 kHz: −56.6/−54.4.
- **Reading:**
  - Up to 1 s the residual is model error, not noise, in every window and band (+14 to +34 dB over noise).
  - Mid bands (500 Hz–2 kHz) fit best: −49 to −53 dB at 50–300 ms. The extremes are worse: 63 Hz and 8 kHz.
  - From 1 s on, the residual equals the noise model (res/noise ≈ 0 dB).

## 10. Injection / diffuser structure

Code: `.venv/bin/python tools/model_s.py --diffuser` (needs `tighten.json` and `tail.json`) → `analysis/model_s/diffuser.json` (log: `analysis/model_s_diffuser.log`). All nulls are in the 20 Hz–11.3 kHz band. **Excess** is over the in-band noise floor of the same window.

### 10.1 The key test: where does the input end?

If the input has stopped by time T, D(z)·h = 0 after T, where D is the product of the four loop factors and the 460 factor. So h after T follows from h before T alone:

h_after = −(1/D) · [(D · h_before) restricted to n ≥ T]

This prediction has **no excitation parameters**: just the 9 loop + 460 numbers, with the measured IR up to T as the state.

| Loop + 460 parameters from | T | Predicted window | L null (excess) | R null (excess) |
|---|---|---|---|---|
| 0.25–0.40 s only (fully honest) | 0.30 | 400–600 | −19.5 (+28.4) | −23.9 (+30.5) |
| 〃 | 0.35 | 400–600 | −22.6 (+25.3) | −29.5 (+24.9) |
| 〃 | **0.40** | **400–600** | **−26.0 (+21.9)** | **−38.8 (+15.6)** |
| 〃 | 0.40 | 600–1000 | −15.6 (+16.7) | −28.5 (+10.5) |
| whole decay (§9; uses data after 0.4 s for these 9 numbers) | 0.40 | 400–600 | −26.3 (+21.6) | −41.7 (+12.7) |
| 〃 | 0.45 | 450–650 | −25.3 (+19.5) | **−45.5 (+5.9)** |
| 〃 | 0.50 | 500–700 | −24.4 (+17.1) | −41.8 (+6.2) |
| 〃 | 0.60 | 600–800 | −20.2 (+14.9) | −38.0 (+3.7) |

- **R:** the prediction improves steadily as T moves from 0.30 to 0.45 s, then stays 4–6 dB above noise.
  - **So R's input ends at about 0.45 s, and from there the 4 loops + 460 stage explain R to within a few dB of the noise floor.**
  - With early-window parameters, the R hold-out at 400–600 ms is −38.8 dB, against −27.5 dB for the generative model Cs (§8.2). Cs fails there because its excitations end by about 80 ms (plus the 460 tail) and so cannot carry input that arrives up to 0.45 s.
- **L:** the prediction stays 15–22 dB above noise whatever T is. **L has something that R does not, all the way down the tail.** That is the L excess of §8.4. It is not input timing.

### 10.2 What the late input looks like

The de-recirculated IR D·h is the input as the loops see it, in dB over its noise-only level (10 ms blocks):
- It is high up to 0.24 s (≥ 52/62 dB for L/R). That is the excitations' own support, ending about 230 + P + 300 + ΣP + 460 samples after tap 0.
- Then it **drops by 20–30 dB and decays smoothly, with no spikes**:
  - R: +52 dB at 0.25 s to 0 dB at 0.43 s, about 290 dB/s;
  - L: +35 to +10 dB over 0.25–0.45 s, about 130 dB/s, and L never reaches noise.
- So the late input is a dense, decaying component with its own decay rate: RT ≈ 0.2 s in R. It is neither isolated taps nor a delayed copy of earlier input.

### 10.3 Hypotheses for the late input

| Hypothesis | Free params | In-sample 50–300 ms | Held-out 400–600 ms L / R (excess) | Held-out 600–1000 ms L / R (excess) |
|---|---|---|---|---|
| Cs reference (§8.2): free excitations + first-pass FIR, 4 loops + series 460 | 16,732 | −47.4 / −46.8 | −28.4 / −27.5 (+19.5 / +27.0) | −19.5 / −18.5 (+12.8 / +20.4) |
| Input ended by 0.40 s; loops + 460 only (whole-decay parameters) | 9 + measured state | n/a | −26.3 / −41.7 (+21.6 / +12.7) | −16.0 / −29.8 (+16.3 / +9.1) |
| + a second 460 stage | 11 + state | n/a | −26.1 / −41.3 | −15.7 / −29.2 |
| + a 920 stage | 11 + state | n/a | −26.1 / −41.2 | −15.7 / −29.2 |
| + a 230 stage | 11 + state | n/a | −26.1 / −40.4 | −15.8 / −28.4 |
| 460 stage as a free 49-tap FIR (frequency-shaped gain) | 49 + state | n/a | −26.2 / −41.6 | −16.1 / −29.7 |
| + AR(1200) on the input (fitted 0.25–0.40 s) | 1,209 + state | n/a | −26.2 / −43.2 | −16.2 / −31.3 |
| Long loop, lags 3,500–16,000 (scan of the input, 0.22–0.40 s and 0.30–0.50 s) | n/a | best single lag removes ≤ 9 %, at sums of loop lengths (4,297, 4,437, 5,759), i.e. artefacts | | |
| Shared mono input d (0–0.4 s) × 8 tap sets of 300 (bilinear ALS) | 12,957 | −7.3 / −7.5 | +1.3 / +1.3 (+49 / +56) | +1.6 / +1.7 |

Reading of each hypothesis:
1. **A long, low-gain loop (EMT 250 style)** is **not present**: nothing at lags 3,500–16,000.
2. **A second 460 stage, 920 or 230 stages**: not present. The 460 stage's per-pass gain is **flat across frequency**: 0.45 / 0.42 / 0.51 / 0.44 / 0.52 / 0.45 at 0.1 / 0.5 / 1 / 2 / 4 / 8 kHz (49-tap FIR). The two-parameter stage is enough.
3. **A shared long input** feeding all loops is **rejected**, on two counts:
   - The fine structure of the Cs excitation spectra (detrended at 20 Hz) is uncorrelated across *different* loops: mean |r| 0.03, against 0.02 for random sequences.
   - It is correlated between L and R *of the same loop*: 0.17 / 0.17 / 0.16 / 0.38. **Each loop has its own long injection**, and its two read-outs are related (the 26/332-lag relation of §5).
   - The bilinear model also fails outright. From a naive start the alternating fit converges to −7 dB in-sample, and its prediction is worse than nothing. Better initialisation was not tried.
4. **The 460 stage as a true all-pass:**
   - Its numerator (−g + z⁻⁴⁶⁰) cannot be tested while excitations are free, because they absorb it.
   - *[exploratory, §8.3]* The excitations keep a −0.31 lag-460 regression coefficient, close to the −g/(1+g²) = −0.37 that an all-pass numerator with g = 0.44 would leave.
   - g₄₆₀ = 0.439–0.458 (σ ≈ 0.02). Simple values such as 7/16 = 0.4375 or 15/32 = 0.469 cannot be told apart.
   - Whether it sits before the loops (common) or per loop: the annihilator cannot tell, because a pole is a pole. The parallel (C) and series (Cs) placements null within 0.5–2 dB of each other.
5. **AR model of the late input:** it explains 89–96 % of it in-window, mostly as spectral colouring (alternating lags 1–9, an HF tilt), plus weak lags 230, 72, 30 and 20. As a predictor it adds only about 1.5 dB in R.

### 10.4 Timing and plausibility

- **Timing:** relative to the pre-delay output (cluster tap 0, the 10 ms pre-delay), the reverb starts about 240 samples later in L and about 570 in R. They differ by 330, the 332 L/R lag of §5. No integer injection taps could be isolated (§8.3), so no tap set can be given relative to the pre-delay.
- **Memory budget:** model C uses 5,758 + 460 = 6.2 k words of loop memory. That's comfortably inside the 39 k free words.
- **Tap budget:** model C needs 2 × 8.4 k taps of dense excitation FIR per channel, against about 50–70 gain-weighted taps per sample for the whole machine. **That is implausible as hardware.** The real injection must come from a compact structure.
- **What that structure would have to do:**
  - feed each loop its own dense sequence;
  - keep arriving, decaying at 130–290 dB/s (RT about 0.2–0.5 s), until about 0.45 s;
  - share fine structure between the L and R read-outs of one loop, but not across loops.

  Per-loop diffusers (a short all-pass chain or nested all-pass in front of each loop, each with its own delays) would fit that description **[H]**. The 230 / 72 / 30 / 20 lags that keep appearing are the obvious first candidates for their delays.

### 10.5 What limits the held-out null now

1. *(Withdrawn in §13.2: the apparent late input is the diffuser tail plus deconvolution band edges.)* **No generator for the per-loop input over 0–0.45 s.** In R, everything else is within 4–6 dB of noise once the input has ended (§10.1). Until the injection structure is found, the best *generative* hold-out stays at about −28 dB (Cs). The structural information is there: a predictor that uses the IR's own past reaches −42 to −46 dB in R.
2. **L-only excess**, 15–22 dB over noise at every T. Neither input timing nor any time-invariant recursion explains it (§8.4, §10.1). It is the main obstacle in L.
3. Free excitations overfit (§8.2), and sparse integer taps fail (§8.3). **The next step is a structural search for per-loop diffusers.** For example, fit a small all-pass chain per loop (delays from {20, 30, 72, 230, 460, …} and gains) to each loop's excitation, judged by this held-out predictor and not by in-sample nulls.

## 11. All-pass diffuser search

Code: `.venv/bin/python tools/model_s.py --allpass` (needs `tighten.json`) → `analysis/model_s/allpass.json` (log: `analysis/model_s_allpass.log`). Rows marked *[exploratory]* come from scratch runs of the same functions.

**Scoring.**
- A candidate diffuser is a fixed filter placed in series with the loops (and with the first-pass path).
- Integer taps in front of it, per loop and channel, are chosen greedily by batch OMP on 0–400 ms only.
- The score is the broadband null at a fixed tap budget (215 taps per channel): in-sample at 50–300 ms, **held out at 400–600 ms**.
- Loop and 460 parameters are the early-window ones (data before 0.4 s).
- The search runs on R; L is used afterwards as an independent check.
- If a filter is the real diffuser, the same few taps should explain far more of the IR.

### 11.1 What was found

| Step (R, 215 taps) | in-sample 50–300 | **held-out 400–600** |
|---|---|---|
| 460 stage only (comb or all-pass) | −3.2 | −2.6 |
| + all-pass 20, g = −0.6 / +0.6 | −5.9 / −2.6 | −5.9 / −1.8 |
| + all-pass 72, g = −0.6 / +0.6 | −5.8 / −2.7 | −5.9 / −2.0 |
| other single stages (30, 52, 100, 144, 230, 690, 920; ±0.6) | −1.8 … −3.9 | −1.2 … −3.7 |
| 460 · AP20 · AP72, gains (−0.5 / −0.6 / −0.7 / −0.75 / −0.8 / −0.9) | −6.6 / −9.1 / **−12.1** / −11.3 / −9.1 / −5.3 | −6.5 / −9.2 / **−13.4** / −12.0 / −9.0 / −4.7 |
| delays 18, 19, 21, 22 instead of 20 | −3.2 … −3.7 | −2.6 … −3.3 |
| delays 70, 71, 73, 74 instead of 72 | −3.3 … −3.7 | −2.7 … −3.2 |
| delays 458, 459, 461, 462 instead of 460 | −7.2 … −8.1 | −7.3 … −8.2 |
| 460 as a comb (pole only) instead of an all-pass | −13.3 | −14.4 |
| AP20 · AP72 with no 460 stage | −10.9 | −11.6 |
| + third stage, full scan of every delay 2–460 at g = −0.6 *[exploratory]*: best 230 (−14.3 / −14.9); next 115 (−9.7), 229 (−9.0), 3 (−8.9) | | |
| **comb460 · AP20 · AP72 · AP230, all at g = −0.7** | **−15.7** | **−17.2** |
| AP230 gain −0.45 / −0.6 / −0.8 / +0.6 | −13.6 / −14.3 / −13.2 / −10.9 | −14.6 / −14.9 / −14.7 / −12.1 |
| + a fourth stage (3, 7, 30, 52, 115, 144, 153; ±0.6) | all worse (−9.7 … −11.1) | |
| common gain −0.65 / −0.6875 / −0.7 / −0.72 / −0.75 | −12.8 / −15.3 / −15.7 / −15.3 / −13.2 | −12.9 / −16.3 / −17.2 / −17.0 / −14.1 |
| **L check** (chain chosen on R): 460 only → + AP20 · AP72 → + AP230 | −3.0 → −13.6 → **−15.9** | −2.4 → −14.5 → **−17.6** |

**Structure found [M]:** a common series diffuser of three all-pass sections with integer delays and equal gain:

AP(20, g) · AP(72, g) · AP(230, g), with g ≈ −0.70; each AP(M, g) = (−g + z⁻ᴹ) / (1 − g·z⁻ᴹ),

combined with the 460 recursion (pole gain 0.445, better modelled as a plain comb than as a fourth all-pass).

- **The delays are sharp:** ±1 sample loses 8–10 dB (20, 72), 4–7 dB (230 at −0.6) and 4–6 dB (460).
- **L agrees:** the chain was chosen on R alone, and it gives L the same improvement at the same gain.
- **Gain:** −0.70 to −0.72 is best. −0.6875 (−11/16) is 0.4–0.9 dB worse, so it is not resolved to a Q15 value **[limit]**.
- **Delay set:** 20, 72, 230 and 460 explain the 20/52/72 tap offsets (§8.3) and the 230/460/690/920 lags (§8.4, §10). AP20 · AP72 alone has echoes at 20, 52 (= 72 − 20), 72, 92, … .
- **Cost:** 3 all-pass sections + 1 comb ≈ 8 multiplies and 4.1 k words of delay at most. That's well inside the budget.
- **Input or read-out side?** Because the chain is common to all loops and both channels and is linear, one IR cannot tell input side from read-out side. The usual design, and the per-loop L/R fine-structure sharing of §10.3, favour an **input diffuser in front of the loops [H]**.
- **What it implies for timing:** with the diffuser in front of the loops, the "late input" decay of §10.2 (R about 380 dB/s in band) is close to AP230's own rate at g = −0.7 (356 dB/s) **[H]**.

### 11.2 Effect on the models (broadband null dB, L / R; excess over the broadband noise floor in brackets)

| Model | Free params | in-sample 50–300 | **held-out 400–600** | **held-out 600–1000** |
|---|---|---|---|---|
| sparse taps, 460 only, 3,345 taps/ch | 6.7 k + 16 | −20.6 / −20.3 | −12.7 / −12.4 (+21.8 / +30.7) | −9.5 / −10.2 |
| **sparse taps + diffuser**, 414 taps/ch | 0.8 k + 19 | −19.8 / −20.5 | −21.1 / −20.6 (+13.4 / +22.5) | −14.5 / −19.6 |
| **sparse taps + diffuser**, 1,719 taps/ch | 3.4 k + 19 | −26.7 / −27.2 | −23.2 / −23.6 (+11.3 / +19.5) | −15.2 / −20.9 |
| **sparse taps + diffuser**, 3,345 taps/ch | 6.7 k + 19 | −31.6 / −32.3 | −25.4 / −26.3 (+9.1 / +16.8) | −15.3 / −21.0 |
| Cs: dense excitations + first-pass FIR, no diffuser (§8.2, broadband) | 16,732 | −36.1 / −37.9 | −27.9 / −26.8 (+6.6 / +16.3) | −14.2 / −17.3 |
| dense + diffuser *[exploratory]*, in band (20 Hz–11.3 kHz) | 16,732 | −43.3 / −49.0 | −29.2 / −29.1 (in-band excess +18.7 / +25.3) | −20.3 / −20.0 |
| (same model without diffuser, in band, §8.2) | 16,732 | −47.4 / −46.8 | −28.4 / −27.4 (+19.5 / +27.0) | −19.5 / −18.5 |

- **For a parsimonious model, the diffuser is worth 10–13 dB at every tap budget.** With 414 taps per channel, the held-out null (−21 dB) is 17 dB better than without it. With 3.3 k taps it comes within 1–2.5 dB of the 16.7 k-parameter dense model.
- **With dense free excitations it adds only 0.8–1.7 dB held out.** Dense excitations can already absorb a short diffuser; what they cannot carry is the late input.
- So **the diffuser is real structure, but it is not the missing late input generator.** Even with it, thousands of taps per channel are needed. The taps in front of the diffuser are still far from a "handful": the injection has more structure we have not found.

### 11.3 What still limits

1. **The injection taps are still dense** (thousands per channel). There must be further structure in front of or around the diffuser: per-loop multi-tap writes, a second diffuser branch, or taps read from a long delay line carrying the diffused signal. One test that fits the 50–70 tap budget: look for a small set of tap delays that is **shared across loops and channels** with loop-specific gains, i.e. group-sparse OMP.
2. **The late input** (§10) is still not generated. In R, the input-ended predictor (§10.1) remains the best held-out evidence (−38.8 dB from early-window parameters, −41.7 dB with whole-decay parameters, at 400–600 ms). No generative model reaches it.
3. **L-only excess** (§8.4, §10.1): unchanged by any of this.
4. **Gains** are resolved to about ±0.02 (diffuser) and ±8 Q15 LSB (loops). That is not bit-exact.

## 12. Group-sparse injection taps, and the late-input generator

Code: `.venv/bin/python tools/model_s.py --groups` (needs `tighten.json`) → `analysis/model_s/groups.json` (log: `analysis/model_s_groups.log`).

**Set-up.**
- Fixed parts: the §11 diffuser (comb460 · AP20 · AP72 · AP230, g = −0.7), the early-window loops, and the cluster with free gains.
- Taps sit at integer device-sample delays after cluster tap 0 (the 10 ms pre-delay output), placed with the (1 + ε) clock stretch.
- Each tap feeds a "slot": loop i (4 slots) or the direct, non-recirculating path through the diffuser (1 slot), in each channel.
- Fit: 0–400 ms. Held-out windows: 400–600 and 600–1000 ms. Excess is over the noise floor of the same band.

### 12.1 Tap delays shared across loops and channels (group OMP)

Each group is one delay feeding all 10 slots (5 per channel) with its own gains.

| Shared delays | Params | 50–300 ms (fit) | **400–600 held out**, broadband (excess) | 400–600 held out, 20 Hz–11.3 kHz (excess) | 600–1000 held out, broadband |
|---|---|---|---|---|---|
| 4 | 48 | −1.0 / −0.9 | −0.9 / −0.8 (+33.6 / +42.3) | −0.9 / −0.9 (+47.0 / +53.5) | −0.9 / −0.9 |
| 8 | 88 | −1.5 / −1.4 | −1.4 / −1.3 | −1.5 / −1.4 | −1.4 / −1.5 |
| 16 | 168 | −2.3 / −2.9 | −2.1 / −2.7 | −2.2 / −2.9 | −2.1 / −2.8 |
| 32 | 328 | −4.4 / −5.3 | −4.6 / −5.1 | −5.1 / −5.6 | −4.3 / −5.3 |
| 64 | 648 | −8.9 / −9.8 | −9.8 / −10.6 (+24.7 / +32.5) | −12.0 / −13.1 (+35.9 / +41.3) | −9.3 / −11.5 |

**Negative result.** A small delay set shared by all loops and both channels does not exist: 4–16 shared delays explain almost nothing. The 64 selected groups show no structure either:
- RMS gain against delay: L −201 dB/s (r = −0.27), R +5 dB/s (r = 0.01);
- sign agreement across loops: 0.45 / 0.30, against about 0.375 for random signs;
- their delays bunch just after the reverb onsets (L at 240; R at 570–1,250).

This agrees with §10.3: different loops share no spectral fine structure.

### 12.2 Sharing within a loop only (L and R read the same delay)

Each group is (loop or direct path, delay), feeding both channels with its own gains.

| Groups | Params | 50–300 ms (fit) | **400–600 held out**, broadband | 400–600 held out, in band | 600–1000 held out, broadband |
|---|---|---|---|---|---|
| 50 | 108 | −6.5 / −6.3 | −7.0 / −6.4 | −7.9 / −7.2 | −6.5 / −6.7 |
| 100 | 208 | −10.5 / −10.9 | −12.5 / −12.7 | −17.4 / −18.0 | −11.5 / −13.8 |
| 200 | 408 | −12.4 / −12.7 | −14.0 / −14.0 | −17.5 / −18.3 | −12.3 / −15.1 |
| 400 | 808 | −15.6 / −15.6 | −16.9 / −16.7 (+17.6 / +26.4) | −19.4 / −19.7 (+28.5 / +34.7) | −13.4 / −17.1 |
| *for comparison:* independent taps per slot and channel (§11), 414/ch | about 830 | −19.8 / −20.5 | −21.1 / −20.6 (+13.4 / +22.5) | n/a | −14.5 / −19.6 |

At the same parameter count, sharing a loop's delays between L and R costs about 4 dB. The L and R read-outs of a loop are related (§5, §10.3), but not by identical tap positions.

### 12.3 Late-input generator: long-delay-line taps (EMT 250 style)

Independent per-slot taps allowed out to 8,500 samples beyond each loop's normal support, i.e. up to about 10,300 samples (0.39 s) after tap 0. These are the "taps read from a long delay line".

| Taps per channel | normal support (ext 300), 400–600 held out | long support (ext 8,500), 400–600 held out |
|---|---|---|
| 414 | −21.1 / −20.6 | −20.4 / −20.1 |
| 1,719 | −23.2 / −23.6 | −22.1 / **+88.1 (R unstable)** |

- **Negative result.** Late taps do not improve the held-out null.
- At 1,719 taps, R's prediction blows up: taps near the end of the fit window are almost unconstrained, and greedy selection fits them to it.
- The exponential gain law could not be tested as a parametric model, because the selected taps show no decay with delay at all (§12.1). Fitting a·exp(−σt) to them would describe noise.
- A long multi-tap line would also inject *bursts* (one diffuser response per tap). The late input of §10.2 is smooth and only mildly heavy-tailed: kurtosis 4–10 over 0.27–0.45 s *[exploratory]*.

### 12.4 The two alternatives

- **(a) Feedback from the loop sum or comb460 output into the diffuser input.**
  - That couples the loops, which moves the modes off the harmonic series. §8.4 found that 12 free cross-coupling gains lower the annihilator cost by nothing (0.1944 against 0.1945). The modes fit Pᵢ(1+ε) per loop to a few mHz (§9.2).
  - So any such feedback must be too small to matter. It is not the late-input generator. **[M, indirect: frequency-dependent coupling through the diffuser was not fitted explicitly]**
- **(b) L-only time-varying injection.**
  - The input-ended predictor (§10.1) keeps L 15–22 dB above noise even when started at T = 0.45–0.6 s, after the input has ended. R falls to +4–6 dB.
  - So **L's excess lives in the recirculating part, not in the injection**: L's read-out of the loops differs from what a fixed set of read-out filters predicts, all the way down the tail **[M]**. The same predictor also rules out time-variance in the injection as the cause.
  - A time-varying L read-out (the 1986 "stochastic control" confined to L) remains consistent with this **[H]**. The line-width test (§8.6) bounds it to modulation slower than about 1 Hz, or to small depths.

### 12.5 What still limits

1. The injection is **per loop, dense, and not sparse in any basis tried**: integer taps, loop-shared taps, L/R-shared taps, the §11 diffuser, or long-delay taps. Each loop's excitation behaves like a long, loop-specific diffused sequence. The best generative held-out results remain:
   - dense model + diffuser: −29.2 / −29.1 dB in band at 400–600 ms (§11.2);
   - sparse + diffuser at 0.8 k params: −21 dB.
2. *(Withdrawn in §13.2.)* **The late input (0.25–0.45 s) is still unexplained.** Every candidate generator tried here fails: a long loop, a 460 cascade, a common all-pass chain, long-line taps, and loop feedback. The only description that works is the R predictor from the IR's own past (§10.1).
3. **L's read-out excess** in the tail (12.4b).
4. With only one setting, structural candidates can't be cross-checked against changes in RT, HF, LF, DELAY or AMPL. Several remaining questions are best settled on the library:
   - does the late input scale with RT (recirculating) or not (feed-forward)?
   - does the diffuser gain change with any control?
   - does L's excess change between captures of other settings?

## 13. Is it "stochastic control"? Time-variance seen through a sine sweep

Code: `.venv/bin/python tools/model_s.py --stochastic` (needs `tighten.json` and `tail.json`) → `analysis/model_s/stochastic.json` (log: `analysis/model_s_stochastic.log`).

**Hypothesis under test:** S/M/L randomly switch or drift their tap gains or positions, in the injection and, on L, in the read-out. Measured through a slow ESS, that would make a sparse system look dense and spoil the L predictor.

### 13.1 Sweep parameters from the IR itself

- No distortion images. Searching the WAV for scaled copies of the cluster pattern (15 ms template, normalised matched filter) finds:
  - in the 10 ms before the first tap: best match 0.35 (L) / 0.34 (R);
  - in the last second: 0.16 / 0.22;
  - typical (median) values: 0.009 / 0.005.

  No discrete images sit at consistent negative offsets. With a slow sweep, the k = 2, 3 images would lie hundreds of ms before t = 0, and the file starts only 10.1 ms before the first tap. **So the sweep parameters cannot be recovered from this file.**
- **The pre-onset energy is diffuse.** L: −52.6 dB re its first tap. R: −78.6 dB. That is **26 dB higher in L**. L is also 18 dB above its own last 0.3 s (−70.3 dB). R is 8 dB above its own (−86.2 dB).
- Another L-only anomaly. It is consistent with L-side time-variance smear, which an ESS spreads into negative time too, but this file can't separate that from a capture-chain artefact **[H]**.

### 13.2 Simulator

- **Synthetic program S** (`Synth246`), built from the fitted structure:
  - integer loops with the fitted filters;
  - comb460 · AP20 · AP72 · AP230 (g = −0.7);
  - sparse injection: 8 random taps per loop, common to L and R;
  - sparse read-out: 3 taps per loop and channel, at lags 0 / 26 / 332;
  - the 8-tap cluster.

  That is about 96 effective taps per channel.
- **Measurement:** an exact ESS (20 Hz–12 kHz, 10 s and 30 s) through the system, deconvolved by regularised spectral division (`deconvolve`). Noise is added at −105 dB re tap.
- **Variants** (taps change every *period*):
  - injection gain jitter (relative rms 3–30 %);
  - injection position jitter (±1 sample);
  - full injection re-draws;
  - L-only read-out gain jitter (1–10 %);
  - combinations.
- **Diagnostics (identical for CLA and the simulations):**
  - batch OMP held out at 400–600 ms, after 103/414 taps. The atoms include the measurement's band-limit (`post·B`); see below for why;
  - the state predictor (§10.1) on a *gently* band-limited IR (30–80 Hz and 9–11 kHz raised-cosine edges, `gentle_band`), given as excess over noise;
  - cross-band envelope correlation of the de-recirculated input (`cross_band`).

**Two corrections to earlier sections, both found with the LTI control:**

1. **The "late input" up to about 0.4–0.45 s (§10) is reproduced by the LTI control,** whose taps all lie within the first 82 ms.
   - Its state predictor behaves exactly like CLA's R: excess +21.0 dB at T = 0.35 s, falling to +9.0 / +8.3 dB at 0.40 / 0.45 s.
   - The apparent late input is the **diffuser's own recursive tail** (AP230 at g = −0.7 decays at about 356 dB/s) **plus band-edge ringing of the deconvolution**: a sharp 10–20 Hz edge rings for hundreds of ms.
   - **No extra late-input generator is needed.** §10.5 item 1 and §12.5 item 2 are withdrawn.
2. **OMP atoms must include the measurement band-limit.**
   - Without it, the LTI control (truth about 96 taps per channel) needs more than 400 OMP taps to reach −18 / −20 dB held out.
   - With it, 103 taps reach **−45.5 / −46.0 dB**.
   - So the poor sparse fits in §8.3 and §11–12 were partly a basis artefact. For CLA, adding the deconv_grid band-limit (11–12.3 kHz) to the atoms raises the 103-tap held-out null from −14.5 / −14.6 to −21.8 / −21.0 dB. **But CLA then plateaus at about −23 dB** (−23.8 / −23.1 at 414 taps), where the LTI control reaches −46 to −49 dB.

### 13.3 Results (held-out 400–600 ms; predictor excess over noise at T = 0.35 / 0.40 / 0.45 s)

| System | OMP 103 taps L / R | OMP 414 taps L / R | predictor excess L | predictor excess R |
|---|---|---|---|---|
| **CLA (measured)** | −21.8 / −21.0 | −23.8 / −23.1 | **+23.6 / +21.1 / +19.3** | **+16.5 / +8.1 / +6.2** |
| LTI control (10 s and 30 s sweeps identical) | −45.5 / −46.0 | −46.5 / −48.6 | +18.5 / +7.8 / +7.7 | +21.0 / +9.0 / +8.3 |
| injection gain 3 % / 0.2 s | −30.1 / −30.0 | −32.5 / −33.2 | +19.3 / +13.7 / +11.7 | +21.5 / +15.4 / +13.6 |
| injection gain 10 % / 0.2 s | −20.0 / −20.0 | −23.1 / −23.9 | +22.5 / +19.9 / +17.3 | +22.8 / +20.0 / +18.3 |
| injection gain 30 % / 0.2 s | −11.1 / −11.6 | −14.0 / −14.7 | +24.0 / +21.7 / +19.1 | +23.2 / +21.0 / +19.4 |
| injection gain 10 % / 0.05 s | −19.1 / −19.2 | −21.0 / −21.2 | +22.6 / +21.2 / +19.0 | +23.6 / +22.4 / +20.5 |
| injection gain 10 % / 1 s | −21.7 / −22.2 | −26.9 / −28.4 | +19.5 / +13.8 / +12.6 | +22.4 / +16.2 / +14.5 |
| injection positions ±1 / 0.2 s | −3.9 / −4.8 | −7.4 / −8.2 | +22.8 / +20.5 / +19.3 | +23.1 / +20.9 / +19.3 |
| injection full re-draw / 0.2 s | −1.0 / −1.0 | −2.4 / −2.8 | +24.7 / +22.7 / +21.0 | +23.1 / +21.4 / +21.3 |
| injection full re-draw / 1 s | −2.3 / −2.6 | −5.1 / −5.6 | +27.4 / +23.0 / +17.7 | +26.9 / +25.2 / +22.1 |
| **L read-out 1 % / 0.2 s** | −38.9 / −46.0 | −38.5 / −48.6 | **+26.3 / +22.1 / +18.9** | +21.0 / +9.0 / +8.3 |
| L read-out 2 % / 0.2 s (10 s sweep; 30 s sweep in brackets) | −33.5 / −46.0 | −32.5 / −48.6 | +31.5 / +27.9 / +24.6 (+30.3 / +28.3 / +24.3) | +21.0 / +9.0 / +8.3 |
| **L read-out 2 % / 1 s** | −34.5 / −46.0 | −35.4 / −48.6 | **+27.1 / +21.4 / +17.5** | +21.0 / +9.0 / +8.3 |
| L read-out 3 % / 0.2 s | −30.1 / −46.0 | −29.0 / −48.6 | +34.8 / +31.3 / +27.9 | +21.0 / +9.0 / +8.3 |
| L read-out 10 % / 0.2 s | −19.5 / −46.0 | −18.4 / −48.6 | +42.9 / +39.5 / +36.1 | +21.0 / +9.0 / +8.3 |
| injection 3 % + L read-out 2 % / 0.2 s | −28.5 / −30.0 | −29.4 / −33.2 | +31.1 / +27.5 / +24.3 | +21.5 / +15.4 / +13.6 |

- **Cross-band envelope correlation** of the de-recirculated input is 0.91–0.98 in every variant, and 0.94–0.99 in CLA. It is dominated by the shared decay envelope, so it **does not discriminate**.

### 13.4 Reading

- **R is time-invariant within this test's resolution.**
  - CLA R's predictor excess after the input ends (+8.1 / +6.2 dB) matches the LTI control (+9.0 / +8.3).
  - Injection jitter of 3 % / 0.2 s already gives +15.4 / +13.6, and 10 % / 1 s gives +16.2 / +14.5. So any randomisation of the injection, which R shares with L, is **below about 3 % rms per 0.2 s**. Larger re-draws are ruled out by a wide margin.
- **L's excess fits an L-only read-out variation of about 1 % rms gain every 0.2 s, or about 2 % every 1 s** (+22 / +21 dB at T = 0.40 against CLA's +21.1). 2–3 % at 0.2 s already overshoots.
  - This is the model's best explanation of the L excess. It also fits the 26 dB higher L pre-onset energy (13.1) and eugtone's "L-only null error".
  - L-only gain modulation of about 1 % is far too small to hear as FM or AM, which fits "ohne wahrnehmbare Frequenzmodulation".
  - **[H]** The same signature would come from any slow L-only variation in the capture chain (e.g. an unstable L converter or cable). One IR cannot tell unit from capture chain.
- **The dense, non-sparse injection is not explained by time variance.**
  - The variants that fit the predictors (R LTI, L read-out about 1 %) leave the simulated OMP at −38 to −49 dB.
  - CLA's OMP plateau at −23 dB in **both** channels must therefore be **structural mismatch**: our loops, diffuser or read-out form differs from the real one in a way that sparse taps cannot absorb.
  - Injection jitter of 10 % or more would produce the plateau, but it contradicts R's predictor.
- **Not supported:** "stochastic control" as switched or re-drawn injection taps. Re-draws at 0.2–1 s destroy both the sparse fit (−1 to −6 dB) and the predictor (+18 to +27 dB). CLA shows neither.

### 13.5 What it means for the model, and what the library can confirm

- **For the plugin:**
  - R path: time-invariant.
  - L path: optionally a slow random gain modulation of the L read-out (about 1 % rms, re-drawn every 0.2–1 s). That costs nothing in hardware terms. Whether it belongs to the unit or to eugtone's capture is open.
  - The injection still needs a structural model, which remains unidentified (§11–12). It is not a stochastic one.
- **The library can test this:**
  1. If L's excess is in the unit, it appears in every capture, with the same rough depth, and in L only. If it is the capture chain, it varies between libraries.
  2. Captures that differ only in RT or HF reuse the same injection. If the injection were randomised, the excitations fitted on different captures would disagree beyond noise; in an LTI unit they agree. That is a direct test of injection randomness without repeat captures of the same setting.
  3. Programs M and L ("Grad der Stochastik" differs per program, 1986 review): run the same predictor test. A larger L or R excess in M/L than in S would point to the unit.
