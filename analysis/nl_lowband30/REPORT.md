# NonLin_program §11.6 bullet 2: the low band (100 Hz-1 kHz) in the first 30 ms. Report (27 Oct 2026)

Scripts: `analysis/nl_lowband30/lb_*.py` (read-only on `tools/` and `analysis/nl/bench/`; `lb_common.py` loads the bench alignment and renders the
model with `nl_model.render`). Outputs: `scan*.json`, `oracle*.json`, `*.pkl`. Metric here: window = file sample 20 ... front-frame zero + 30 ms, residual
e = y - m with the bench's own o / clock / gain, 4th-order band-pass, untrimmed captures with EARLY >= 0.3 (steep negative EARLY has the front taps
40-60 dB down, so that window is noise), median over captures of the L/R mean. My absolute numbers differ from the doc's (different window and
filter); the pattern reproduces.

## 1. What the error is (S = several measurements agree)

1. **It is a sub-300 Hz error.** Per band, first 30 ms (all captures): 20-100 Hz -34...-46, 100-1k -51...-58, 150-1k -53...-59, 300-1k -54...-59,
   1-5k -55...-58 dB. Above 300 Hz every group is alike (-53...-59); the spread between groups sits in 20-300 Hz, rising towards LF. "100 Hz-1 kHz -49
   dB" is the upper skirt of it. (In E/M terms, front section [0, 25 ms], 20-32 Hz -25...-33 dB, 50-125 Hz -38...-55, 125-315 Hz -41...-55, >=500 Hz -55...-65.)
2. **Time course: it is a step in the residual at the file start**, not an error on the taps. Cluster-mean full-band e (units 1e-4 of the 410-tap peak)
   is flat in the pre-roll (cache index < 45) and jumps to ~0 at **cache index 45-46, four device samples before file sample 0 (index 50)**, in every
   session (front-shift cluster) and both channels (`lb_step.py`): R -7.5 (fs 0) / -2.9 (fs 20) / -4.0 (fs 23) / -0.5 (fs 16) before, +1 after; L +2 ... -0.5. After the
   step e is -90 dB or lower re the peak (the doc's "after 30 ms all alike"). The pre-tap baseline itself (-12e-4 of the peak) is in data and model alike; only
   the difference steps. Band-limited (100-1k) residual by absolute index (`lb_time_rep2.py`, re the peak): indices 0-100, -71...-90 dB (worst SPAN 0.5 fs 0, -71...-76 dB; SPAN 0.0 / 0.1, -85...-91); after the first tap -88...-100 dB.
3. **Same in L and R, not the same size**: |c| median (1e-4) L / R = 0.4-0.6 / 1.8-1.9 for SPAN 0.0, 0.1; 2.5-3.0 / 3.0-3.8 for SPAN 0.2, 0.3; 2.8 / 5.0 for 0.4; 1.4 / 4.2-7.7 for 0.5. R is bigger and mostly the opposite sign of L in the fs 0 session.
4. **It is deterministic per capture group, not noise and not a common model error**: for captures sharing EARLY (LATE <= 0.6, o 0-1.7), e in [-40, 25 ms] correlates
   0.6-0.86 within a front-shift cluster (different SPAN / LATE) and -0.06...0.37 between clusters (`lb_corr.py`). A fixed transfer error of the model would correlate across clusters; it does not.
5. **Two separate contributions**, which the doc's "SPAN 0.3/0.5 sessions" blends:
   - **A LATE (late-level) dependent part, all sessions.** Held to one (SPAN, EARLY), the first-30-ms null worsens monotonically with LATE: SPAN 0.4 EARLY 10, 100-1k -56 (LATE 0) -> -35 (LATE 0.6) -> -31 dB (LATE 2.0); R's pre-roll step grows with LATE (1.0 -> 5.8e-4). It depends on the late *level*: EARLY 3 / 10 (late at 0 dB) is bad, EARLY 0.3-2 at SPAN 0.5 (late 40-75 dB down) is not. The SPAN 0.5 "fs 0" session holds **only EARLY 3 and 10** (17 captures), the SPAN 0.5 "other" session only EARLY 0.3-2: the doc's session contrast at SPAN 0.5 is mostly EARLY / late level. (Model side: the late reads' signs sum to +-2 per loop, `nlt.READ_SIGNS`, so the model's late has a large DC / sub-20 Hz gain (+12 dB re 1 at 1 Hz at LATE 2) that the processing's acausal low edge turns into a pre-roll baseline; a small DC / LF mismatch of the real late then shows up as exactly this step. H.)
   - **A genuine session / SPAN part at LATE 0.** With no late section at all, fs 0 (SPAN 0.3) is -51.7 dB against -56.4 for the others (100-1k, n 9 / 38; 20-100 Hz -34.5 vs -44...-48). Controlled for late level and LATE, fs 0 is 1-6 dB worse than the other sessions in every cell. SPAN 0.2 (fs 23) is as bad (-52.8); SPAN 0.4 and 0.5 (fs 23) are not (-55.8 / -57.6). So SPAN 0.2 / 0.3 (and 0.5 fs 0) carry a step at the file start even without a late (|c| L / R = 2-3 / 3-4e-4).
6. **Not a front-tap / transfer error:** the three front taps and the first main tap carry less than -90 dB residual at 100 Hz-1 kHz (full-band cluster means at t = 42, 408-412, 460-466, 500 samples are 5-30e-4 of the peak, but the 100 Hz-1 kHz band residual there is <= -86 dB re the peak in every 40-160-sample bin after the first tap, `lb_time_rep.py`); two captures of one session with the same front taps (gain ratio 0.989) can differ 0.2x in e (SPAN 0.5, EARLY 10 vs 3), so e does not scale with the front-tap level (additive, not multiplicative).

## 2. Capture artefact or model error? (what the evidence supports)

- **Capture / processing artefact at the file start, S for the form, H for the cause.** The error is a baseline-step mismatch at the file start (the point where Tomkin's file ends and the deconvolution's zero PAD begins), session-specific (point 4), additive (point 6) and the same sign pattern per session in L and R. It looks like a DC / sub-20 Hz level that the file carries (or lacks) at its start: the zero-phase inverse (3-12 Hz edge) converts it into a constant before the file start and the residual keeps a step. This matches the E 16.7 finding of a non-causal pre-tap low band in the captures' chain.
- **Part of it is a model LF error that the file start merely exposes (H):** the LATE-growth (point 5) and the +-2 sign sums say the late's DC / sub-20 Hz gain is not what the unit has (e.g. a DC-blocking path or truncation in the loop; the repeat captures cannot show it because it is a *static* LF level). The same level mismatch lives in the whole response at 3-20 Hz, in line with the known 20-100 Hz reverb band (-40 dB) and the "20-100 Hz" row of §11.1.
- Ruled out here: front-tap / main-tap error in the band (point 6); leakage of the model's own branch into the window (the branch leak is -64...-69 dB re the peak, e is 20 dB below it; `lb_leak.py`, apart from captures whose branch offset o is -64 / -84, where the first main tap at 685 + o falls inside the window: exclude them from the 30 ms windows, or end the window at 685 + o - 20); a session-independent model error (point 4).

## 3. Test of a compact correction, held out (`lb_oracle.py`, `lb_oracle2.py`)

Correction: subtract c . u(t) from the capture, u = tanh edge at cache index 45.5 (width 1.5), c = per channel constant.
- **Measured in the capture's own pre-roll** (indices 8-40, no content, outside every scored window; 2 numbers per capture, none fitted on the scored window): first-30-ms medians 100-1k **-51.4 -> -54.1 (fs 0 sessions), -53.1 -> -56.4 (other sessions)**; 20-100 Hz -35.2 -> -38.6 and -37.4 -> -47.4; SPAN 0.2 -51.8 -> -56.6, SPAN 0.4 -51.4 -> -54.4, SPAN 0.5 (fs 23) -53.6 -> -57.1, SPAN 0.5 (fs 0, 6 captures) -41.1 -> -45.6. The excess over the noise floor was not recomputed per band (the doc's 0-0.05 s full-band floor is -78.8 dB); the corrected nulls are still well above it and above the repeat level (-66.5 dB).
- **Predicted from other settings (held out by EARLY)**: c_L = group median, c_R = a + b . LATE per (SPAN, session) group, fitted on the other EARLY values of that group (3 numbers per group, 9 groups, 27 parameters for 497 captures): 100-1k -51.6 -> -53.2 (fs 0), -53.1 -> -55.6 (others); median change -0.87 dB, better by > 0.5 dB in 304, worse by > 0.5 dB in 18 of 497. **A held-out fixed correction works for the session / SPAN part and keeps 2/3 of the measured gain, so the effect is largely a per-group constant.** It does not remove the whole gap: fs 0 stays 2.5 dB behind the others.
- No correction is proposed for the plugin: the step lies before file sample 0 and in the metric; the device does not make it. What the plugin needs from this is only the late DC check below.

## 4. Status of the doc's hypothesis H and the labels

- S: the excess is < 300 Hz, an additive step in the residual at the file start (cache index 45-46), per capture group, both channels, R larger.
- S: it is not the front taps' read shape or a transfer error on them, and it is not noise (within-group correlation 0.6-0.86).
- S: a large part of "SPAN 0.3 / 0.5 sessions bad" is late level and LATE (the SPAN 0.5 fs 0 group has only EARLY 3 / 10); the session contrast should be quoted at LATE <= 0.4 / EARLY 0.3-2.
- H: the cause is a DC / sub-20 Hz level at the capture's start (file offset or the processing's LF edge applied to a truncated start), plus a model LF error of the late (sign sums +-2).
- Ruled out: trim, front / branch path difference (already), a session-independent transfer error, the branch leak.

## 5. Recommended next step

1. In `nl_model.py` bench: score the first-30-ms windows only for captures whose o is 0-1.7 (or end the window at 685 + o - 20) and quote the contrast at LATE <= 0.4; report the step-removed null (pre-roll mean subtracted, outside the window) beside the raw null. This takes 2.7-3.3 dB off the "session gap" at no model cost.
2. Test the late's DC gain: refit the late reads with the sign-sum constraint relaxed or a first-order DC block in the loop feedback (a 3-parameter test: one corner frequency, fitted on LATE 1.0-2.0 R pre-roll baselines of SPAN 0.4 / 0.5 EARLY 10 at held-out SPAN). If the R step vs LATE slope (R c = -1.8 -> -7.3e-4 over LATE 0 -> 2 at SPAN 0.4) disappears, the LATE part is a model LF term.
3. Do not spend a model parameter on the session part; document it as a capture feature (c = 0.4-5e-4 of the peak, 1-5 dB in this metric).

---

# Part 2: the 20 Hz-1 kHz question for Program N as a whole (27 Oct 2026, second brief)

Scripts: `lb_rep.py` (repeat level), `lb_model_band.py` (model null in the same window), `lb_low*.py`, `lb_tf*.py`, `lb_biquad.py`, `lb_edge.py`
(spectral analysis, T(f) and parametric corrections). Alignment/gain/clock from `bench_bk_a0.65` (env `NL_BK=1 NL_BK_ALPHA=0.65 LB_BENCH=bench_bk_a0.65`).

## 6. The repeat-capture level, reproduced (`lb_rep.py`)

Method: SPAN 0.5, captures of one EARLY that differ only in LATE (the late starts at 0.465 s, after the window), same front-shift cluster **and the same
branch offset o (|do| < 0.02)**, window front-zero + 0.05...0.45 s, 4th-order band-pass, gain + fractional shift + clock fitted per pair, band and channel,
median of up to 20 pairs per EARLY (L/R mean, raw pair null; subtract 3 dB for one capture's noise).

| EARLY | 20-100 | 100-1k | 1-5k | 5-10k | 10-12k |
|---|---|---|---|---|---|
| 10 | -66.8 | -73.7 | -66.6 | -65.2 | -59.1 |
| 3 | -68.3 | -70.7 | -63.8 | -62.7 | -59.5 |
| 2 | -64.0 | -69.3 | -62.8 | -61.9 | -61.2 |
| 1.5 | -62.0 | -68.1 | -61.6 | -60.5 | -60.1 |
| 1.0 | -58.0 | -66.7 | -59.7 | -58.7 | -58.0 |
| 0.7 | -55.1 | -63.9 | -57.8 | -56.4 | -55.7 |
| 0.5 | -48.5 | -54.7 | -52.3 | -50.9 | -51.5 |
| 0.3 | -44.5 | -51.2 | -49.7 | -49.1 | -48.7 |
| 0.2 | -38.8 | -50.1 | -48.3 | -47.3 | -44.5 |
| **pooled, all EARLY** | **-57.1** | **-66.0** | **-59.1** | **-57.9** | **-56.1** |
| **pooled, EARLY >= 0.7** | **-62.0** | **-68.2** | **-61.7** | **-60.6** | **-59.4** |
| doc §11.1 | | -66.5 | -60.6 | -59.4 | -58.9 |

It reproduces §11.1 within about 1 dB (the pooled figures bracket it). **Why the orchestrator's rough check was 12-18 dB worse:** it used EARLY 0.1-0.3 only. There the
early section decays at 60 dB per 0.3...0.1 s, so by 0.05-0.45 s the signal is within reach of the capture noise: the same method on EARLY 0.2-0.3 gives
-44...-51 dB at 100 Hz-1 kHz, i.e. the repeat level is a signal-to-noise figure that depends on the EARLY ramp. Pairs with a different o (a different branch
offset, even 0.2 samples) are not repeats of one early section and show -5...-40 dB at 5-12 kHz: filter them out. Compare model and repeat on the same EARLY.

## 7. Model versus repeat, same captures and window (`lb_model_band.py`, 130 captures, SPAN 0.5)

Model null, bench alignment (alpha 0.65 branch kernel) fixed | repeat level, median L/R, dB:

| EARLY | 20-100 model / repeat | 100-1k | 1-5k | 5-10k | 10-12k |
|---|---|---|---|---|---|
| 10 | -41.0 / -66.8 | -58.8 / -73.7 | -63.1 / -66.6 | -60.7 / -65.2 | -35.2 / -59.1 |
| 2 | -40.2 / -64.0 | -57.6 / -69.3 | -61.0 / -62.8 | -57.7 / -61.9 | -33.4 / -61.2 |
| 1 | -39.3 / -58.0 | -56.8 / -66.7 | -59.4 / -59.7 | -57.5 / -58.7 | -34.8 / -58.0 |
| 0.5 | -37.4 / -48.5 | -54.3 / -54.7 | -56.8 / -52.3 | -55.1 / -50.9 | -36.8 / -51.5 |
| pooled model | -39.2 | -56.6 | -59.1 | -57.3 | -35.6 |
| pooled repeat | -57.1 | -66.0 | -59.1 | -57.9 | -56.1 |

Gap (model - repeat, pooled): **20-100 Hz +18 dB (+22...+27 at EARLY >= 1), 100 Hz-1 kHz +9 (+13 at EARLY >= 1), 1-5 kHz 0, 5-10 kHz +0.6, 10-12 kHz +20.**
Re-fitting gain, fractional shift and clock in the band (the repeat method's freedom) changes the model null by < 1 dB at 20-1000 Hz: the gap is not an alignment effect.
**Verdict: 20-100 Hz is not capture-limited. It is a model-side error 18-27 dB above the repeat level, and 100 Hz-1 kHz carries about 10 dB of it too (S).** At EARLY <= 0.3 both bands are noise-limited and the gap closes, which is why the weak ramps hide it.

## 8. What the 20-100 Hz error is (S measured, H cause)

1. **Deterministic, setting-independent, not a start effect.** Local null in 25-ms bins from 0.05 to 0.45 s is flat: -36...-44 dB in 20-100 Hz and -55...-59 dB in 100 Hz-1 kHz in every bin, every EARLY, both sessions (`lb_tcourse.py`), apart from the worst-case bins at the file start (the §1-5 step) and, for EARLY 0.7, the last 0.1 s (low signal). So it is proportional to the local signal: a transfer error, not an additive artefact.
2. **It is a fixed LTI error of the same size and sign in both channels and both sessions.** The pooled data/model transfer ratio T(f) = <Y M*>/<|M|^2> (`lb_tfshow.py`, 33 captures, +-12 % bands): 20.8 Hz -0.02...-0.24 dB / -1.5 deg, 24 Hz -0.1...-0.18 dB / -1.1...-0.7 deg, 28 Hz -0.07 dB / -0.5...-0.2 deg, 33-100 Hz |T-1| < 0.05 dB / < 0.6 deg, 11-13 Hz +0.7...+2.0 dB, 8 Hz +6...+10 dB, 6 Hz +14...+19 dB (the three below 15 Hz are in the processing response's 3-12 Hz edge region, where the model itself is small). The numbers agree to 0.1 dB / 1 deg between fs 0 and fs 23, L and R. These are 1-3 % deviations, i.e. -34...-40 dB: exactly the size of the error.
3. **A held-out fixed T(f) removes 10-14 dB.** T estimated on one group, applied to the other (40 log-spaced bands, 2 numbers per band and channel, so 80 numbers; spectral null, window 0.05-0.45 s): train fs 23 (EARLY 0.7-2) -> test fs 0 (EARLY 3, 10): **20-100 Hz -41.8 -> -53.0 dB**; train fs 0 -> test fs 23: **-40.0 -> -50.5**; leave-one-EARLY-out (4 folds): -39.2...-41.0 -> -48.6...-54.5. By sub-band (fs 23 -> fs 0): 20-25 Hz -31.7 -> -46.9, 25-32 -38.6 -> -48.3, 32-50 -42.3 -> -54.2, 50-100 -50.0 -> -55.7; 100-300 and 300-1000 Hz unchanged (T applied only inside the band; those bands keep their -55 / -59). The repeat level is -57...-62, so about 4-8 dB above it remains in 20-100 Hz (not closed; and 80 numbers is not a compact model).
4. **It is not a compact filter I could find.** Held out the same way: a 4-number biquad T(f) fitted on 16-200 Hz gives only +1...+2 dB (-40 -> -41...-43; fitted on 8-200 Hz it is worse than nothing); re-shaping the 3-12 Hz sin^2 edge (a, b free) gives 0...+1.6 dB; a first-order high-pass corner gives +0.8...-5 dB; edge x HP gives up to +3 dB. The structure of T (non-monotone: +17 dB at 6 Hz, about 0 at 9.6 Hz, +2 dB at 11 Hz, -0.2 dB at 21 Hz) is sharper than those forms, which looks like a *window / edge* effect (a 3-12 Hz edge that differs in shape, or the deconvolution's regularised inverse at the chain's own LF corners) rather than a device filter: the S / M / L chain kernel (`chain3`) was fitted on other programs.
5. **Ruled out:** independent per-tap relative gain errors (the taps do not cancel at LF here: the amplification of a per-tap relative error is 0.0...+0.8 dB in every band, `lb_amp.py`, so tap-gain errors cannot give a 17-dB-bigger error at 20-100 Hz); a noisy band (repeat floor is -57...-81 dB at 20-100 Hz in the spectral comparison, `lb_low_an.py`); the 20-100 Hz errors of the file-start step (part 1) - the step explains only the first 30 ms.
6. **The bench's 20-100 Hz number also sees the sub-20 Hz mismatch**: the 4th-order band-pass has a 24 dB/oct skirt, T is +2...+19 dB at 6-13 Hz, so a share of "20-100 Hz -40 dB" is the skirt of the 6-15 Hz mismatch. (My spectral null uses only bins >= 20 Hz and gives the same -40, so the 20-32 Hz bins themselves are bad.)

## 9. Verdict and next step

- **20-100 Hz: model error (processing-response / chain-kernel LF behaviour), not capture-limited; the doc's "-40 dB, 8 dB above the repeat level" should read 18-27 dB above (EARLY >= 1; the wording "100 Hz-1 kHz 8 dB above" also stands at +9...+13).** Label: S (measurement), H (cause: the LF shape of the chain kernel x 3-12 Hz edge in `nl_structure.processing_response`, or the S/M/L kernel's LF end, not the device).
- **A fixed correction exists and generalises** (80 numbers, 10-14 dB held out by session and EARLY, SPAN 0.5 only, 33 captures) but no compact parametric form (biquad, edge, HP) reproduces it.
- **Recommended next step:** (a) take T(f) = data/model on the 15-60 Hz grid from *all* SPANs and sessions (the full 1,044 bench, not only SPAN 0.5) and check it is SPAN-independent; (b) if so, tabulate it as an additional zero-phase+phase correction of `nl_structure.processing_response` (bench side, not the device - the plugin does not need it) and re-run the bench: expected 20-100 Hz -40 -> -50 dB; (c) look for the physical term by fitting T(f) with the chain kernel's own HP corners (`k.hp_hz`, first-order analog corners at the deconvolution's inverse) and the F_LO window as free, using the *full* kernel structure rather than a stand-alone biquad; (d) the same method on 100 Hz-1 kHz (gap +9...+13 dB) with a time-domain, not LTI, ansatz (T(f) did not help there: -58.1 -> -56.6).
- Open: whether T differs between L and R at < 15 Hz (the estimates differ by up to 5 dB at 6-8 Hz); whether it follows from the late section's DC gain (part 1, sign sums +-2).

---

# Part 3: closing the low band: what fixes it and what does not (27 Oct 2026, third brief)

Scripts: `lb_tf_all*.py`, `lb_tf_eval*.py`, `lb_fineT*.py`, `lb_final.py`, `lb_stage.py`, `lb_apg.py`, `lb_envel.py`, `lb_additive.py`, `lb_edge.py`, `lb_biquad.py`.
Deliverables: `fineT_all_sm6.npz` (recommended table) and `lowband_T.json` (coarse log grid), loader `lowband_fix.py`.

## 10. The error is a fixed linear transfer error of the whole pipeline (S)

- **Independent of setting.** The pooled data/model ratio T(f) at 20-60 Hz is the same for SPAN 0.2, 0.3, 0.4, 0.5, for both sessions (fs 0 / 20 / 23), both channels, EARLY 0.3...10, and for LATE 0.0...1.0 (file lengths 0.9...5.7 s): -0.1...-0.2 dB and -1...-2.6 deg at 19-25 Hz, +0.03...+0.06 dB around 40-55 Hz, bootstrap sd < 0.01 dB above 25 Hz (`lb_tf_tab.py`, `lb_tf_late.py`). Below ~18 Hz T varies by session (the file-start step of Part 1 plus the 3-12 Hz edge) and is not used. So it is **not** a file-end/truncation effect (T does not change with file length), not a tap-gain, stage-gain or all-pass-coefficient error (below), and not capture noise.
- **Size by band** (rms |T-1|, 70-point table): 20-100 Hz -41 dB (0.9 %), 100-300 Hz -54, 300-1000 Hz -59...-63, 1-3 kHz -64 dB. The deviation falls smoothly with frequency, the signature of a **chain-kernel fit error** (`chain3` was fitted on S/M/L cluster windows: 3,731-tap FIR + 4 analytic high-pass corners): T = G_true / G_kernel. H.
- **Ruled out for the 100 Hz-1 kHz gap (S):** free per-stage gains (15 numbers per capture, in-sample: -56.4 -> -56.6, -58.6 -> -58.7 dB); all-pass coefficients g(AP124), g(AP180/188) (fitted dg 1e-5...1e-3; leave-EARLY-out change < 0.2 dB); per-tap relative gain errors (no LF amplification, `lb_amp.py`); an error that concentrates at the taps (3 % of the error energy within +-8 samples of tap positions, which are 3.5 % of the window); a smooth or low-order T (biquad 4 numbers +1...+2 dB; re-shaped 3-12 Hz edge 0...+1.6; HP corner -5...+1; rational fits of order 2/4/6 +0.6/+2.3/+5.5 dB at 20-100 Hz, none helps 100 Hz-1 kHz).
- **What does remove it is fine frequency structure.** A pooled per-bin T (0.4 Hz bins, smoothed +-2.4 Hz, 51 captures, SPAN >= 0.2, LATE <= 0.2, EARLY >= 0.3, window front-zero + 0.05 s ... min(0.45 s, late feed - 10 ms)) applied as an LTI filter to the model record.

## 11. Held-out result (leave-one-SPAN-out; time-domain 4th-order band nulls as in the bench; median of 102 channel-captures; correction estimated without the SPAN under test)

| Band | raw | corrected | SPAN 0.2 | 0.3 | 0.4 | 0.5 (raw -> corrected) |
|---|---|---|---|---|---|---|
| 20-100 Hz | **-39.8** | **-48.2** | -38.5 -> -47.6 | -38.7 -> -46.5 | -40.3 -> -49.0 | -40.1 -> -48.7 |
| 100-300 Hz | **-52.8** | **-60.0** | -52.4 -> -61.5 | -52.3 -> -58.0 | -53.3 -> -57.2 | -53.5 -> -62.3 |
| 300 Hz-1 kHz | -58.1 | -60.8 | -57.8 -> -61.5 | -58.4 -> -58.8 | -57.1 -> -57.7 | -58.5 -> -64.1 |
| 1-5 kHz | -59.6 | -59.8 | | | | |
| 5-10 kHz | -56.6 | -56.6 | | | | |

(+-2.4 Hz smoothing; 0 of 102 channel-captures worse than +0.5 dB at 20-100 Hz and 100-300 Hz; 101 of 102 better than 3 dB at 20-100 Hz, 84 at 100-300 Hz.) **With the late section present** (204 captures with LATE >= 0.4, window 0.05-0.6 s, T estimated from the other SPANs): 20-100 Hz -38.6...-39.6 -> -46.5...-47.6, 100-300 Hz -51.7...-52.8 -> -54.7...-60.5, 300 Hz-1 kHz +0.7...+3.2 dB, above 1 kHz unchanged. A first estimate from only 12 captures of one SPAN (SPAN 0.5, one session) transferred to the other SPANs and sessions at 20-100 Hz -40 -> -44...-47.5 and 100-300 Hz -53 -> -58...-60, so it is not a within-session artefact. Plus the file-start step (Part 1, 2 numbers per capture measured in the pre-roll): 20-100 Hz a further -49.8 -> -50.7 dB.
**Parameter count:** the fine table is 2 channels x about 2,500 complex bins (14 Hz-1 kHz at 0.4 Hz, smoothed): it is a correction filter, not a compact model. The 70-point log-grid version (280 numbers) reaches only 20-100 -48.5, 100-300 -56.2, 300 Hz-1 kHz -57.8 (28 of 102 worse in 300 Hz-1 kHz): the structure is finer than 6 %. The 36-point band-averaged LF table (17-150 Hz) reaches 20-100 Hz -49.8 alone.
**Against the repeat level** (pooled SPAN 0.5: 20-100 -57...-62, 100 Hz-1 kHz -66...-68): 20-100 Hz is now about 9-14 dB above it, 100 Hz-1 kHz about 6-8 dB (repeat level per band measured only for 100 Hz-1 kHz as a whole). The remaining error is model-side but no longer a smooth transfer error.

## 12. Verdict

- **20-100 Hz: closed to -48 dB held out (-50.7 with the file-start step removed), against -40 before.** The -50 dB minimum is met only with both corrections; as a pure LTI correction it stops 2 dB short.
- **100 Hz-1 kHz: 100-300 Hz now -60 dB, 300 Hz-1 kHz -60.8 dB (from -52.8 / -58.1); both within 6-8 dB of the repeat level.** (The bench quotes the whole 100 Hz-1 kHz as -55.2 / -56.8 at alpha 0.65; I did not recompute that exact figure with the correction in the bench pipeline.)
- **What the correction is.** An empirical estimate of the chain kernel's residual error for Program N's captures, not a device property: it does not belong in the plugin. It belongs in `nl_structure.processing_response` (bench side) as `R *= T`, per channel, `lowband_fix.fine(f, ch)`. If T is a property of `chain3`, **the same table should improve S, M and L**, which use the same kernel (not tested here: check it first on S, where §37 puts the Altiverb-chain-limited nulls at -39...-41 dB).

## 13. Recommended next steps

1. Fold `fineT_all_sm6.npz` into `processing_response` behind a switch (bench only), re-run the N bench and the per-band table: expected from the held-out test above: 20-100 Hz -40 -> -48 dB, 100-300 Hz -53 -> -60 dB, 300 Hz-1 kHz -58 -> -61 dB.
2. Test the same T on S, M, L captures (SPAN-like leave-one-setting-out). If it transfers, the right fix is to refit `chain3` once (a kernel `chain3b` with T folded in) rather than carrying a Program-N table.
3. If the file-start step (Part 1) is wanted in the bench, subtract the pre-roll mean (indices 8-40 of the cache, 2 numbers per capture) from the capture before band-limiting.
4. Not resolved: why T has structure finer than 6 % (smoothing +-2.4 Hz is the best of 0 / 1 / 3 / 6 / 12 bins: 20-100 -48.3 / -48.3 / -48.3 / -48.2 / -47.8, 100-300 -59.3 / -59.4 / -59.7 / -60.0 / -60.1); the sub-18 Hz region (T +2...+19 dB, session-dependent); the remaining 5-8 dB at 100 Hz-1 kHz to the repeat level.

---

# Part 4: is T(f) a chain3 kernel error? Test on Program S (27 Oct 2026, `lb_s_test.py`, `s_test.log`)

Material: the stored S accumulators of `sml_fullband.py --stage s` (`analysis/fullband/s/*.npz`, read only): per setting, channel and 2 Hz bin the sums Y M*, |M|^2, |Y|^2 over S's reverb (150 ms to the compensated end) for S's bench model (own clock/shift/gain fit per capture, held out by RT fold), all 127 scored settings. T_N is my `fineT_all_sm6` (interpolated to the 2 Hz bins). Error energies with a given T follow from the sums: |Y - T M|^2 = |Y|^2 + |T|^2 |M|^2 - 2 Re(T conj(Y M*)). S's own fixed transfer T_S is `fold_T` pooled over LF >= 1.0 (what §41.3 calls S's loop-path transfer).

**(a) T_S and T_N are the same kind of object at low frequency.** rms |T - 1| (S / N): 20-40 Hz 2.1-2.3 % / 1.5-1.6 %, 40-80 Hz 0.34-0.60 / 0.37-0.63, 80-160 0.17-0.25 / 0.10-0.19, 160-270 0.23-0.32 / 0.21-0.23, above that S 0.2 % / N 0.05-0.09 %. Complex correlation of (T_S - 1) with (T_N - 1), L / R: **20-40 Hz 0.94 / 0.89, 40-80 Hz 0.88 / 0.47, 80-160 Hz 0.82 / 0.40**, 160-270 Hz 0.59 / 0.55, 270-800 Hz 0.12-0.45, 0.8-1.6 kHz 0.72 / 0.71. (The correlation phase is -31...-40 deg at 20-40 Hz: unresolved, see below.)

**(b) N's table applied to S's model, nothing refitted, S's whole 127-setting bench (median band error energy over settings, dB re |Y|^2; L / R):**

| Band | T = 1 | T = T_N | settings better / worse by > 0.5 dB (L) | T_S (other RT fold, held out) |
|---|---|---|---|---|
| 20-40 Hz | -35.3 / -35.7 | **-38.5 / -36.7** | 83 / 5 | -40.2 / -39.4 |
| 40-80 Hz | -44.6 / -47.8 | **-48.4 / -48.0** | 91 / 4 | -49.6 / -49.0 |
| 80-160 Hz | -47.2 / -48.5 | -48.6 / -47.6 | 93 / 17 (R: 14 / 63) | -51.2 / -49.5 |
| 160-270 Hz | -48.2 / -47.8 | **-50.9 / -49.7** | 114 / 1 | -53.4 / -53.5 |
| 270-400 Hz | -50.1 / -51.7 | -51.5 / -51.6 | 91 / 0 (R: 16 / 72) | -54.7 / -55.5 |
| 400-800 Hz | -52.3 / -51.6 | -52.9 / -51.9 | 81 / 0 | -55.8 / -55.5 |
| 0.8-1.6 kHz | -51.4 / -50.7 | -52.5 / -51.5 | 83 / 3 | -55.2 / -54.7 |
| 1.6-3.2 kHz | -51.1 / -50.7 | -51.2 / -50.5 | 6 / 0 | -52.8 / -52.8 |

**Reading (S for the numbers, H for the cause).** A transfer measured on a different program, a different part of the unit and different capture sessions improves S at 20-270 Hz in the left channel by 3-4 dB at 20-80 Hz and 2.7 dB at 160-270 Hz, in 80-90 % of the settings, with almost none worse; S's own held-out transfer gains roughly twice as much (5-5.5 dB at 20-80 Hz). So:
- **A large common component exists** (correlation 0.9 at 20-40 Hz, 0.8-0.9 at 40-160 Hz in L): it is shared by S and N, i.e. it belongs to what S and N have in common, the capture chain, the deconvolution kernel `chain3` with its analytic high-pass corners and the processing response. That supports the kernel/processing explanation for the sub-160 Hz part.
- **A program-specific part exists on top** (S: the loop-path/AP230 structure of §41.3, which §41.4 shows does not transfer to M/L; N: whatever makes its 100-300 Hz structure), and the right channel agrees much less (correlation 0.40-0.47 at 40-160 Hz; N's table makes S's R worse at 80-160 and 270-400 Hz in 63 / 72 settings). The R channel's sub-300 Hz T differs between programs, which a pure kernel error would not allow unless the kernel's R side is also program-sensitive (the R DAC phase, the L half-sample offset or the R pinned gain).
- **Not settled:** a pure `chain3` error would give identical T for S and N in both channels. It does not (L high, R mixed), so the kernel is *a* contributor, not the whole story. M and L cannot arbitrate (their time-variant reverb has repeat levels of -27...-38 dB, far above a 0.1-1 % transfer).

**Recommended next step.** Estimate one joint T from N and S together (they are the only programs with a deterministic reverb) and check it against both held out (N by SPAN, S by RT fold). If both sit within about 1 dB of their own T, refit `chain3` once as `chain3b` (L first) with it folded in; if S's R still disagrees, keep the correction per program (N: this table; S: its loop-path transfer) and call the sub-160 Hz common part the kernel's.

---

# Part 5: joint N + S low-band transfer, on the current models (27 Oct 2026, `lb_joint_sums.py`, `lb_joint.py`, `lb_joint2.py`, `lowband_joint_T.json`)

State re-checked after the S/M/L close-out: S's baseline unchanged (`analysis/fullband/s/*.npz` valid), N's default is now the long-kernel warp with the adopted branch kernel at 0.7 (`bench_warp2_bk_kf_a0.7`, 1,044 captures). The N sums were recomputed on that model (51 captures: SPAN >= 0.2, LATE <= 0.2, EARLY >= 0.3, window 0.05 s ... min(0.45 s, late feed)). Both programs are scored with the same device (S's stored per-2-Hz-bin sums; N's recomputed in the same convention), every number held out: **N leave-one-SPAN-out, S with the other RT fold**. Joint T = N and S with equal program weight (each normalised by its own 20-1600 Hz |M|^2).

Median band error energy, dB re |Y|^2, L / R; T estimated from: none | N | S | joint N+S.

Scored on **N**:

| Band | none | T_N | T_S | joint |
|---|---|---|---|---|
| 20-40 Hz | -33.5 / -31.5 | -40.9 / -37.5 | -35.9 / -32.8 | **-40.8 / -38.4** |
| 40-80 Hz | -42.3 / -46.8 | -53.3 / -52.2 | -47.9 / -48.1 | **-51.6 / -52.3** |
| 80-160 Hz | -53.7 / -56.1 | -59.1 / -56.2 | -53.6 / -53.7 | -56.2 / -54.4 |
| 160-270 Hz | -52.5 / -51.3 | -58.7 / -55.5 | -52.6 / -50.9 | -55.2 / -53.4 |
| 270-400 Hz | -52.7 / -57.0 | -58.0 / -57.5 | -52.1 / -50.9 | -55.9 / -54.9 |
| 400-800 Hz | -58.6 / -58.8 | -61.9 / -58.1 | -52.9 / -51.0 | -57.7 / -54.8 |

Scored on **S**:

| Band | none | T_N | T_S | joint |
|---|---|---|---|---|
| 20-40 Hz | -35.3 / -35.7 | -38.2 / -35.6 | -40.2 / -39.4 | **-39.6 / -38.1** |
| 40-80 Hz | -44.6 / -47.8 | -48.0 / -47.7 | -49.6 / -49.0 | **-50.0 / -48.6** |
| 80-160 Hz | -47.2 / -48.5 | -48.8 / -48.4 | -51.2 / -49.5 | -50.6 / -49.6 |
| 160-270 Hz | -48.2 / -47.8 | -50.9 / -50.2 | -53.4 / -53.5 | -53.2 / -53.3 |
| 270-400 Hz | -50.1 / -51.7 | -51.9 / -52.2 | -54.7 / -55.5 | -54.7 / -55.0 |
| 400-800 Hz | -52.3 / -51.6 | -53.0 / -51.8 | -55.8 / -55.5 | -55.7 / -54.9 |

**Reading.**
- **Below about 100 Hz there is one common transfer** (S for the numbers; the cause is H: the capture chain / `chain3` kernel / processing response at the low end). The joint T is within 0-1.7 dB of each program's own held-out transfer in 20-40 and 40-80 Hz on both programs, both channels (N 20-40: -40.8 / -38.4 against own -40.9 / -37.5; 40-80: -51.6 / -52.3 against -53.3 / -52.2; S 20-40: -39.6 / -38.1 against -40.2 / -39.4; 40-80: -50.0 / -48.6 against -49.6 / -49.0). T_N and T_S themselves agree to 0.05 dB and 0.5 deg at 33-61 Hz in both channels (`lowband_joint_T.json`); at 21-29 Hz the phases agree within about 1 deg but the magnitudes differ by 0.1-0.2 dB (R at 29 Hz: -0.10 vs +0.03 dB).
- **Above about 100 Hz the programs differ** and the joint transfer is the wrong compromise: at 80-160 Hz it keeps about 1/2 of N's gain (-59.1 -> -56.2) and nearly all of S's (-51.2 -> -50.6), from 270 Hz up it costs N 2-4 dB and S 0-1.4 dB against their own. T_S on N is harmful above 270 Hz (-58.6 -> -52.9 at 400-800 Hz L); T_N on S is harmless. So 100 Hz-1 kHz needs per-program transfers (S: the loop-path transfer of model_s_results 41.3; N: this work's fine table), the low end a shared one.
- **So the right fix is two-level:** (1) a **shared correction below 100-130 Hz** applied to the kernel (or to `processing_response`, which every program's model uses), `joint` columns of `lowband_joint_T.json`; (2) the per-program fine transfers above, which for N is `fineT_all_sm6.npz` blended in above 130 Hz.
- **Compactness.** The joint transfer on 15-100 Hz is 43 two-Hz bins per channel; rational fits of it (real coefficients, s-domain, T(inf) = 1, SK iterations; `lb_joint2.py`) reach the table's quality only at order 10 (20 numbers per channel): S 20-40 / 40-80 Hz -39.2 / -48.8 (table -38.8 / -48.9), N -39.4 / -49.8 (table -39.9 / -52.4). Orders 2-8 leave 1-6 dB. So it is a measured table; the ripple (+-0.05 dB, +-0.5 deg, period about 25-30 Hz) is the whole content.

**Verdict on "is it a chain3 kernel error?":** the sub-100 Hz part, yes (H, on two programs and two RT/SPAN families held out); the 100 Hz-1 kHz part, no: it is program-specific (N's table is N's, S's is S's loop paths), so the chain is not the whole story there.

**Recommended action.** Give `processing_response` (nl_structure; the S/M/L equivalents) a shared `lowband_joint_T` correction below 100 Hz (blend to 1 at 130 Hz), then keep N's fine table above; re-run the N bench (expected 20-100 Hz about -40 -> -48 dB) and S's gate. Decision for the user: whether the 20-100 Hz band is worth a measured table in the bench (it is not a device property and does not enter the plugin), or whether the headline numbers should simply be quoted with the explained 8-10 dB at 20-100 Hz as "kernel limit".

---

# Part 6: the shared sub-100 Hz table folded in, behind switches, and re-run (27 Oct 2026)

**Code (uncommitted):**
- new `tools/lowband_joint.py`: `response(f, mode)`; tables `analysis/lowband/chain3_lowband_joint.json` (shared N+S transfer, 12-130 Hz, 2 Hz bins, tapered 1 -> T over 12-16 Hz and T -> 1 over 100-130 Hz) and `analysis/nl/lowband_N_fine.npz` (N's own fine table, re-estimated on the current model, <= 3.5 kHz); built by `analysis/nl_lowband30/lb_make_tables.py` from `lowband_joint_T.json` and `fineT_all_sm6.npz`.
- `tools/nl_structure.py` `processing_response`: multiplies the (now complex) response by `lowband_joint.response(f, NL_LOWBAND)`; **`NL_LOWBAND=0` (default) is byte-identical to before**. `NL_LOWBAND=1`: the shared table. `NL_LOWBAND=2`: the shared table below 100 Hz and N's own table above 130 Hz (linear blend over 100-130 Hz).
- `tools/nl_model.py`: bench folder suffix `_lb<mode>` when the switch is on.
- `tools/sml_fullband.py`: `SML_LOWBAND=1` multiplies the shared table onto S's `full_model` pipeline response (cluster and reverb); S results go to `s_lowband/` (stage `s`) or the gate folder.

**N bench** (`analysis/nl/bench_warp2_bk_kf_a0.7_lb1`, `_lb2`, 1,044 captures each, default alignment fits redone with the correction in the model; baseline `bench_warp2_bk_kf_a0.7` = the committed default). Medians over captures and channels; `lb_bench_cmp.py`:

| | 20-100 Hz | 100 Hz-1 kHz | 1-5 kHz | 5-10 kHz | windows 0.05-0.15 / 0.15-0.3 / 0.3-0.6 s |
|---|---|---|---|---|---|
| default | -39.5 | -55.8 | -59.1 | -56.6 | -58.7 / -57.6 / -52.2 |
| `NL_LOWBAND=1` | **-46.8** | -55.9 | -59.1 | -56.6 | -58.7 / -57.6 / -52.2 |
| `NL_LOWBAND=2` | **-46.8** | **-60.8** | -59.5 | -56.7 | **-59.4 / -58.1** / -52.3 |

- The 983 captures NOT used to estimate the tables (the table came from 61 captures: SPAN >= 0.2, LATE <= 0.2, EARLY >= 0.3, untrimmed): 20-100 Hz -39.6 -> -46.8 (mode 1 and 2), 100 Hz-1 kHz -55.9 -> -56.0 (mode 1) / **-61.2** (mode 2): the gain holds on SPAN 0 / 0.1, steep negative EARLY and LATE > 0.2 captures, which the estimation never saw. The 61 estimation captures: -38.9 -> -46.8 and -53.0 -> -56.4.
- Strict rule (channel-windows worse by more than 0.5 dB, 0.0-0.05 ... 0.6-1.0 s and the 20-100 Hz band): **mode 1: 0 of 11,000+ worse, worst +0.03 dB. Mode 2: 1 of 1,966 at 0.15-0.3 s worse by +0.52 dB (the strict rule's 0.5 is exceeded by 0.02 dB), nothing else above +0.32 dB.**
- The 20-100 Hz band at -46.8 dB is still 3 dB short of the -50 dB bar; the remaining pieces are the file-start step (Part 1, about 1 dB at 20-100 Hz, a capture feature) and the joint table's loss against N's own (about 1.3 dB: own -48.2 held out).
- The windowed nulls (100 Hz-10 kHz) move by 0.0 (mode 1) and 0.2-0.7 dB (mode 2) because the long windows are dominated by higher bands.

**S full-band gate** (`sml_fullband.py --stage gate --tag lowband`, `SML_LOWBAND=1`): **PASS** (16 S settings in 12 bands at 50-150 / 150-300 / 300-600 ms; M / L json identical to the control: the switch is wired into S's `full_model` only). The control run with the switch off reproduces `gate_ref` (PASS, identical). S medians, L / R, dB (change against the reference):

| Window | 20-40 Hz | 40-80 Hz | 80-160 Hz | 160 Hz-3.2 kHz |
|---|---|---|---|---|
| 50-150 ms | -26.8 / -32.0 -> -30.5 / -38.5 (-3.6 / -6.6) | -41.7 / -47.1 -> -53.2 / -51.5 (-11.4 / -4.5) | -0.6 / -0.3 | unchanged |
| 150-300 ms | -35.6 / -36.1 -> -44.9 / -45.0 (-9.3 / -8.8) | -44.7 / -49.3 -> -51.5 / -50.9 (-6.8 / -1.7) | -0.35 / -0.28 | unchanged |
| 300-600 ms | -35.5 / -36.6 -> -39.8 / -43.2 (-4.3 / -6.6) | -40.1 / -44.2 -> -47.3 / -46.8 (-7.3 / -2.6) | -0.29 / -0.33 | unchanged |

(The 16 gate settings are inside the 127 LF >= 1.0 settings from which S's half of the joint table was estimated, so these S figures are partly in-sample; the held-out S figure from Part 5 is 20-40 Hz -35.3 -> -39.6 and 40-80 Hz -44.6 -> -50.0 dB, L. The full 127-setting S run with the switch (`sml_fullband.py --stage s`, 13 min, 7 GB) was not done.)

**Recommendation.** Adopt `NL_LOWBAND=1` as N's default (gain 7 dB at 20-100 Hz, nothing else moves, nothing worse); `NL_LOWBAND=2` is the user's call under the strict rule (+5 dB at 100 Hz-1 kHz, +0.7 dB at 50-150 ms, one channel-window 0.02 dB past the 0.5 dB limit). M / L: not wired (their time-variant reverb cannot show it); S: switch available for the full-band tools. Neither enters the plugin (a capture-processing correction). The docs (`NonLin_program.md` 11.2 / 11.6, `model_s_results.md` 41) are not edited; this report and `lowband_joint_T.json` carry the numbers.
