# Program N (Non-Lin) — initial investigation

*Status note (October 2026): a first measurement pass on the Tomkin Non-Lin set (1,046 settings at 96 kHz), read against the printed operation card. Not yet a model. Confidence as elsewhere: **[measured]** = direct measurement here; **[C]** = printed in a primary source; **[H]** = hypothesis.*

**Scope:** program N only. S/M/L are covered in `EMT246_unit_facts.md` §5–§7 and `model_s_results.md`; E in `EMT250_program.md`.

## 0. The operation card — transcription (`sources/device-images/emt246_operation_and_function.jpg`)

*Until now only §3 of `EMT246_unit_facts.md` summarised the N panel; the card itself was not transcribed in full. The text below is read from the image.*

**Header** — "EMT 246 S / OPERATION AND FUNCTION". Program entry fields:
- **Prog.: N** — `Delay (ms)`, `Span (s)`, `early —— Reverb Time —— late`
- **Prog.: D** — `left —— Delay —— right`, `left —— Echo (%) —— right`
- **Prog.: S, M, L, E** — `Delay (ms)`, `Amplitude (%)`, `low. freq. —— Reverb —— high freq.`

**Left panel — "Prog. S, M, L: Short, Medium, Long Reverb. / Prog. E: EMT 250 Reverb Program"** — "Time Characteristics of Reverberation": level axis "Pegel LEVEL", time axis "Zeit/TIME". A dashed "Verzög. DELAY" precedes a burst at "Rückwurfamplitude REFLECTION AMPL.", then a straight decay line "Nachhall REVERB", then an arrow "Tiefpaß LOW PASS". Text:
- "Reverberation Time 1000 Hz = R_T"
- "Reverberation Time 50 Hz = R_T × LF"
- "Reverberation Time 7000 Hz = R_T × HF"

**Middle top — "Optimaler Nachhallzeitbereich / OPTIMUM REVERB TIME RANGE"** — SHORT 0.2–1 s (full to ~2.5), MEDIUM 0.8–4 s, LONG 3–10 s (full to 20), on a 1…20 s axis.

**Middle bottom — "Prog. D: Delay/Echo"** — input splits to `L Verzög. L DELAY` and `R Verzög. R DELAY`, each with an `ECHO (%)` gain and a feedback summing node; outputs `Links LEFT`, `Rechts RIGHT`.

**Right panel — "Prog. N: Non Lin"** (the source for everything below):
- axes "Pegel LEVEL" / "Zeit/TIME"; "Verzög. DELAY" at the start;
- "Anhalldauer SPAN" marked `0…0,5 s`;
- **"Frühe Rt / EARLY RT"** at the **start** of the span, **"Späte Rt / LATE RT"** at its **end**, both drawn at the **same top level**;
- two straight decay lines **crossing in an X**, with a **"+" above and "−" below** the crossing (the early RT can be set positive or negative, i.e. falling or rising), a **dashed vertical line at the span end**, and a **thin dashed tail from the early line continuing past SPAN**;
- "LOW PASS: 2 kHz… 12 kHz".

**Reading [C/H]:** the early section runs for SPAN at a rate set by EARLY RT (falling for +, rising for −), then the late section runs at a rate set by LATE RT. The card draws both sections starting at the same level; §5 shows the hardware instead continues from the early section's **end** level.

## 1. Data

- **Grid:** SPAN {0, 0.1, 0.2, 0.3, 0.4, 0.5} × EARLY {−3.0, −2.0, −1.0, −0.7, −0.5, −0.2, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 10.0} × LATE {0 … 2.0 in 0.2 steps} = 1,056, of which **1,046** exist. LP 12 kHz and DELAY 0 throughout.
- **Files:** 0.42–5.76 s, median 1.63 s; trimmed and faded per setting, so **RT must come from the label, not the file** (`library_notes.md` §2). Peak −33…−45 dBFS. First arrival at ≈0 ms; the forward-decay case peaks at ≈16.5 ms (the early build-up).
- Naming `LP_DELAY_SPAN_EARLY_LATE_dc`, implied decimal on the last three fields (`library_notes.md` §2).

## 2. The two-section envelope — confirmed

Broadband envelope (20 ms RMS, dBFS re the first window), SPAN 0.5, LATE 2.0:

```
EARLY      0ms   50   100   200   300   400   500   600   800  1000  1500
 +0.1      -58  -79  -90  -93  -94 -105  -126  -126  -126  -126  -126   (file 0.48 s)
 +1.0      -55  -65  -67  -75  -78  -84   -88   -89   -90   -97  -123   (file 1.25 s)
 +3.0      -55  -64  -64  -68  -68  -69   -72   -75   -80   -85   -90   (file 1.79 s)
+10.0      -55  -63  -62  -64  -62  -61   -62   -65   -71   -77   -88   (file 2.02 s)
 -1.0      -84  -89  -85  -82  -74  -66   -62   -65   -71   -77   -85   (file 3.49 s)
 -3.0      -65  -72  -70  -71  -66  -63   -62   -65   -71   -77   -88   (file 2.02 s)
```

- **Positive EARLY = falling** early section, steeper for smaller EARLY; the tail is short because the late section inherits a low level at SPAN (file lengths 0.48–2.02 s). **[measured]**
- **Negative EARLY = rising** early section (time-reversed), peaking at ≈SPAN; the tail is long because the late section inherits the peak (file 3.49 s). **[measured]** — matches the card's "+/−" and "inverse" reading.
- The two sections are visible as one slope to SPAN and a different one after; see §3–§4.

## 3. Early section — rate and reversal

- **Rate follows EARLY.** §2: the early slope steepens monotonically from EARLY +0.1 to +10.0. A rough 20–380 ms slope fit gives RT ≈ 1.0 / 1.2 / 1.4 / 1.9 / 2.8 / 4.8 s for EARLY 0.1 / 1.0 / 1.5 / 1.5… (the fits are biased by the short window and the trimmed files; the card and the label are the authority). **[measured trend]**
- **Reversal is a time-reversal of the early section, not a negative gain.** Negative EARLY gives a rising envelope; its starting level rises with |EARLY| (a slowly-decaying early is reversed into a slowly-rising one) and it peaks just before SPAN. §2 and §5. **[measured]**
- **Start level:** the forward case starts high (≈−55 dBFS) and the reversed case starts low (−65…−88 dBFS, lower for shorter |EARLY|) — the same envelope mirrored. **[measured]**

## 4. Late section — rate and loops

- **Rate follows LATE.** SPAN 0.5, EARLY −1.0 (late loud), envelope after 500 ms:

```
LATE      500   600   800  1000  1500  2000
 0.0      -62   -90   -98  -116  -116  -116   (gate)
 0.8      -62   -69   -83   -87   -89   -91
 1.2      -62   -67   -77   -84   -88   -90
 2.0      -62   -65   -71   -77   -85   -84
```

  **LATE 0 cuts the tail off at SPAN (gated); larger LATE gives a slower late decay.** (Exact RT fits are unreliable — the files are trimmed and faded; label is the authority.) **[measured]**
- **Loops of the late section.** Decay-flattened ACF (the `harness/core.py: find_loops` method, band 100 Hz–8 kHz, n≈8σ):
  - medium/long LATE: **1108 / 1254 / 1386 device samples** (42.0 / 47.6 / 52.6 ms), the harness's `1108-1254-1386` set (310 files);
  - short LATE (e.g. 0.4): **628 device samples** (23.8 ms), the harness's `628` set (35 files).
  So the late loop set is **LATE-dependent** — the "adapted density" idea, as in S/M/L. **[measured]**
- A few detections also return ~898 device samples; it is a spurious peak in the same scan, not part of the harness's sets. **[measured]**

## 5. Level at SPAN — the card and the data differ

Card reading: the late RT "starts at the end of SPAN at the same level as the early start" (both drawn at the top level; `EMT246_unit_facts.md` §3 marks this **H**). **The data does not support that:**

```
EARLY +1.0, LATE 2.0      start    SPAN−30 ms   SPAN+30 ms
 SPAN 0.1                 -55.5      -66.4        -68.2
 SPAN 0.3                 -55.4      -75.1        -78.3
 SPAN 0.5                 -55.4      -86.6        -89.5
```

The level is **continuous across SPAN** (jump −0.1…−3.9 dB), and the late section starts at the early section's **end** level, not its start. So the hardware **concatenates with continuity**; the card's same-level drawing is schematic. **[measured; card reading H]**

## 6. Consistency check — the early section is invariant across LATE

SPAN 0.5, EARLY −1.0, the full 0–500 ms envelope is **identical to the displayed precision** for every LATE 0.0…2.0 (−84, −87, −84, −82, −75, −66, −62 at 0/50/100/200/300/400/500 ms). The late section alone carries LATE. This is the library's "settings that differ only in LATE are the same capture" note (`library_notes.md` §3, §5) confirmed section by section. **[measured]**

## 7. Gated and inverse cases

- **Gated:** LATE 0 with a long EARLY (the early runs ~flat, then LATE 0 cuts the tail at SPAN). The envelope falls ≈20 dB over ~2 ms near SPAN and stays there — a gate, not a hard digital zero. **[measured]**
- **Inverse:** negative EARLY (rising early to a peak at SPAN), as in §2. **[measured]**
- **Short SPAN** starts the late section earlier and at a higher level, so the whole tail is longer; **long SPAN** hands the late section a low level and the file is short. **[measured]**

## 8. Open questions and discovery plan (October 2026; answered in §9–11)

- **Is the early section a reverb or a shaped FIR?** The card's dashed tail past SPAN hints at a process that keeps running. In the 0–SPAN window I find **no robust loop structure** at the late-loop lags (ACF ≈ 0.00 at 1108/1254/1386), even when the early is nearly flat (EARLY +10). Either the 0.5 s window is too short for the ACF, or the early section is a shaped multi-tap FIR (the `algorithm_research.md` §6 **N-A** option: "a gain-shaped sparse multi-tap FIR … up to t = SPAN"). **This is the first thing to settle.** Test: device-domain tap counting and EDC cut sharpness at SPAN (the doc's Q for §6).
- **Which loop set does the late section use, and does it switch with LATE?** {1108,1254,1386} vs {628} — map the crossover LATE and compare with S/M/L's "adapted density".
- **How is the late section fed?** Continuity at the early's end level (§5) suggests the late is released at SPAN from the early's state (the N-A "fed from a tap at SPAN" option), not a parallel sum. Test: exact phase continuity across SPAN in the device domain.
- **Exact laws:** RT_early(EARLY) and RT_late(LATE) need device-domain EDC fits on un-trimmed renders (or fresh captures), since the library trims. The card gives the intended straight-line-in-dB decay.
- **Lineage:** the 251/252 Non-Lin is "two concatenated processes" (`EMT246_discovery.md` §10), which fits the measured splice. The RMX16 / 480L precedents (`algorithm_research.md` §6) suggest the early section is a shaped-gain FIR — consistent with the missing early loops.

## 9. The early section, decoded (26 October 2026)

*Device-rate IRs of all 1,046 settings (`tools/nl_structure.py`, chain3 kernel, file sample 0 at K0; the files are trimmed, so no onset search), fits in `tools/nl_early.py` on the 96 LATE 0 files (nothing but the early section is there). Labels as in `EMT246_unit_facts.md`: **C** confirmed, **S** strong, **H** hypothesis.*

### 9.1 What changes with EARLY and SPAN (S)
- At fixed SPAN, every EARLY gives **the same waveform times a smooth gain**: window correlation with EARLY 10 is 0.95–0.98 (positive) and 0.88–0.89 (negative) in every 20 ms window. A negative EARLY is **a rising gain ramp, not a time-reversed signal**.
- At fixed EARLY, every SPAN gives **the first part of SPAN 0.5's waveform** at 0.0 dB (correlation 0.89–1.00): one fixed process, cut by SPAN. After the cut it dies at ≈ −440 dB/s (the diffuser's tail, §9.2): the card's thin line continuing past SPAN.

### 9.2 Structure (S)
- **Front taps:** three plain taps per channel, no diffuser, the place of S's cluster: **L 41 / 410 / 463, R 41 / 278 / 463**, amplitudes **1 : 2 : 1** in units of a main tap (R's 463 negative).
- **Main taps:** the input through **AP(124)**, then **27 taps** (`TAPS`: 685, 1230, 1767, 2523, 2817, 3747, 3835, 4821, 4902, 5775, 5988, 6697, 7005, 7587, 7953, 8445, 8832, 9271, 9642, 10065, 10383, 10827, 11055, 11557, 11658, 12192, 12255), then **AP(180) in L and AP(188) in R**. All three all-passes at **g = −1/√2** (S's diffuser gain; free feed-forward and feedback gains land within 2·10⁻⁴ of it, ±1 sample on any delay costs 13–16 dB). Undoing the diffuser leaves 94 % of the early energy in 50 taps.
- **R's taps are L's 221 samples earlier**, each with its own sign (`SIGNS`). R's first tap is the exception: its direct (feed-forward) term sits at 465 in the front taps' frame, its delayed terms in the main taps' frame (§9.4); with that, R's first-tap residual goes from −18 to −45 dB.
- **Stages:** the taps share gains in pairs, 14 stages (the first a single tap): {685}, {1230, 1767}, {2523, 2817}, … {12192, 12255}. The tap intervals shrink through the train (two interleaved families whose spacing falls by 32 and 69 samples per step), so the density rises with time.
- **SPAN** switches whole stages: SPAN 0 / 0.1 / 0.2 / 0.3 / 0.4 / 0.5 runs 0 / 2 / 4 / 7 / 10 / 14 stages (last tap 1767, 3835, 7005, 9642, 12255). Each cut falls 480–870 samples before SPAN × fs; consistent with "a stage runs if it ends at least 15–33 ms before SPAN" (H for the steps between the library's).

### 9.3 The EARLY law (S; every stage step and front tap to ±0.01 dB)
Law time τ in 4 ms ticks: front taps at L 0 / 4 / 12 ms, R 0 / 8 / 12 ms; one stage = 36 ms (9 ticks).
- **EARLY > 0:** stage j at τ = 16 + 36 j ms, gain **−60 τ / EARLY dB** (stage steps 0.72 / 2.16 / 4.32 / 7.2 dB at EARLY 3 / 1 / 0.5 / 0.3).
- **EARLY < 0:** the last stage at τ_end = SPAN + 12 ms, the others 36 ms apart before it, gain **+60 (τ − τ_end) / |EARLY| dB**: the ramp rises over the whole SPAN (front taps included) to the flat level. This places the front taps at −2.11 / −2.89 / −2.67 / −2.44 / −1.23 stages before stage 0 at SPAN 0.1 … 0.5, as measured.
- **EARLY 10 is flat** (0 dB per stage, not 6 dB/s): the top step means "no slope".
- **The late section is fed at the law's gain of stage min(ns, 13)** for EARLY > 0 (T = 16 / 88 / 160 / 268 / 376 / 484 ms for SPAN 0 … 0.5; measured from the late level 16 / 88 / 160 / 267 / 377 / 487) and at 0 dB for EARLY < 0.

### 9.4 Per-capture branch offset (S as data; H as cause)
- The AP(124) branch (main taps **and** the late section) sits o samples from the front taps, o per capture: modal +1 (§12), else about −63…−64, −84…−87 or −163, mixed at random among captures of one (SPAN, EARLY) that differ only in LATE. R's first-tap direct term (§9.2) does not move with o.
- **H:** the front taps and the branch run on different processors joined by a buffer whose pointer phase is set at program load (2 × TMS32010). For the model, o is one alignment per capture; the plugin uses o = +1.
- **Corrected (§11.4):** the large "offsets" (−486 … −9,316 samples) of steep negative EARLY at long SPAN are **not** device behaviour. Those responses start 40–60 dB down (the front taps and the ramp's first stages), so Tomkin's trim (the file starts at its first loud sample) cut the quiet start away; the front taps are then gone and the "front frame" the alignment finds is arbitrary. Every tap inside the file is at its law gain with the full diffuser, and the tails of the trimmed taps are present after the cut. Only the small offsets are the unit's.

### 9.5 Held-out score (`nl_early.py --stage law`, `analysis/nl/early_law.json`)
Per capture only the alignment (front shift, o, clock) and one gain per channel; every structural number is shared and no setting's gains are fitted. On the 96 LATE 0 captures, 100 Hz–10 kHz, from 20 samples before the first front tap to 3000 samples past the last main tap: **median L −45.3, R −45.1 dB** (the free per-tap fit, ≈ 31 gains per channel per file, in-sample: −44.2 / −44.4). Clock ratios +1…+14 ppm.
- *Superseded by §11 (final bench):* the "dense residual" below was mostly two things found later: the deconvolution's own response near the band edge, which the model did not carry (§11.3), and a small read shape on every main tap (§11.5). The paragraph is kept for the record.
- **Floor:** captures that share their early section (settings differing only in LATE) agree to −52…−62 dB, so ≈ 10 dB of model error is left. The residual is dense, HF-heavy (70 % of it at 6–10 kHz while the data is white), follows the signal's envelope, repeats across captures (correlation 0.99 in L) and is not a linear function of the main taps (a shared 63-tap output kernel, per-tap 5-tap shapes, all-pass damping or gain changes each gain ≤ 1 dB; a 1,300-tap shared kernel 3–4 dB held out). Local timing errors of ±0.003 samples would explain its size and spectrum (H: sample-level timing in the unit's output); open.
- Steep negative EARLY at long SPAN puts the front taps below the noise (−60 dB re the ramp's end), so their windows are noise-only.

### 9.6 Model size
Structure: 27 tap positions, 2 × 27 signs, 3 front positions per channel with 1 : 2 : 1 amplitudes, R's offset (−221) and first-tap direct position (465), three all-pass delays at one gain, 14 stages and the SPAN → stage table; law: 4 ms tick, 36 ms stage, front ticks, the 16 ms / 12 ms offsets. Per sample: 3 + 27 taps per channel and 3 all-passes (≈ 36 multiplies per channel), within the 2 × TMS32010 budget.

### 9.7 Defects
`defects.NL_DEFECTS`:
- (SPAN 0, EARLY −0.5, LATE 0): a 0.27-sample slip between the front taps 410 and 463 (−15.7 dB against the model; every other file with visible front taps holds 369.00 / 53.00 samples to ±0.2);
- (SPAN 0.1, EARLY −1.0, LATE 0): a foreign tail, a −52 dB plateau (re the early section) decaying about 17 dB over 5 s with an ACF peak at 537 samples, in a gated capture whose neighbours end at 0.55 s (§11.1).

## 10. The late section, decoded (26 October 2026)

*`tools/nl_late.py` (per-frequency excitation solve, renderer), `tools/nl_model.py` (whole program), on the clock-corrected captures.*

### 10.1 What it depends on (S)
- **LATE only.** At fixed LATE every (SPAN, EARLY) gives the same late waveform (correlation 0.96–1.00), shifted in time by SPAN and scaled by EARLY. **RT = LATE exactly and flat in frequency** (125 Hz–8 kHz within ±3 % from LATE 0.4 up; LATE 0.2 reads 0.23–0.29 s, the short-RT estimator bias). No LF / HF damping. **LATE 0 = no late section at all**, not even a first pass (the gate).
- **Feed:** the AP(124) line, one sample before a main tap: **684 / 1766 / 3834 / 6696 / 9270 / 12254** for SPAN 0 / 0.1 / … / 0.5 (to ±0.002 samples, −50 dB nulls with the late shift free); R's late is fed 220 samples earlier. The late moves with the per-capture branch offset (§9.4). (The earlier note in §4 that the late set "switches with LATE" was the excitation's structure, see 10.2.)
- **Level:** EARLY > 0: the early law's gain at stage min(ns, 13), τ = 16 + 36 min(ns, 13) ms (measured 16 / 88 / 160 / 267 / 377 / 487 ms for SPAN 0 … 0.5; the late's level is within 0.2 dB of the law down to −50 dB); EARLY < 0 and EARLY 10: 0 dB. Where the law puts the late below about −55 dB (EARLY 0.3 / 0.5 at long SPAN) the captures are 6–12 dB louder than the law (H: a floor or a small leak; ≥ 55 dB down, left open).

### 10.2 Structure (S; equal-gain read taps C-level consistent: all 12 within 0.4 %)
- **Input path:** AP(124) (shared with the early section) → the feed → **AP(314)**, g = −1/√2 (the excitation's 314-sample chains fall 3 dB per step with an all-pass's sign pattern +, +, −, +, −).
- **Three uncoupled flat loops, 1108 / 1254 / 1386 samples** (ACF: each loop's multiples, no sum lags), gain 10^(−3 P / (fs · LATE)). A joint per-frequency solve over ten LATE values (one LATE-independent part + one excitation per loop, nothing else) fits every SPAN 0 capture to **−59…−70 dB** (the noise); the LATE-independent part is empty (−22…−29 dB re the loops).
- **Four equal read taps per loop**, the same positions in both channels, per-channel signs; the last one sample short of the loop length: **1108: 301 / 712 / 923 / 1107; 1254: 398 / 689 / 1022 / 1253; 1386: 582 / 897 / 1179 / 1385**. Read amplitude 0.5 main-tap units (0.2494–0.2500 of the 410 front tap). S has five equal signed taps per loop; this is the same idea.
- **R = L's network 220 samples earlier** (the early's R taps are 221 earlier). H: one TMS32010 per channel running the same program with its own offsets and signs.

### 10.3 What is left in the late section
With the session clock (11.6 ppm) fixed, the late model nulls SPAN 0 captures at −47…−55 dB at 50–300 ms; captures sharing a setting agree to −57…−64 / −50…−57 dB there. Free per-loop gain corrections are −11…−15 millidB per pass (an ≈ 1 % RT shortfall, as in M / L) and gain 0.3 dB. The later windows are noise-limited.

## 11. Bench: the whole program held out on the library (26–27 October 2026)

`tools/nl_model.py --stage noise | bench | report` → `analysis/nl/bench_warp2/` (default, no kernel) and `analysis/nl/bench_warp2_bk_kf_a0.7/` (`NL_BK=1`, the adopted branch kernel, §11.5); the figures of 26–27 October on the old clock warp are kept in `bench/` and `bench_bk_a0.65/` (`NL_WARP2=0`). Since 27 October (§11.6) the clock correction uses a long interpolation kernel (`nl_late.unwarp`) and the model carries the harness's own resampler filter instead of a fitted EQ. 1,044 captures (`NL_DEFECTS` out), 96 kHz set only. **Per capture only:** the front taps' trim shift (measured, bias-corrected), the branch offset o, the clock ratio (11.6 ± 3 ppm) and one gain per channel; every structural number and law is shared, no setting's gains are fitted. The model goes through the captures' own processing response and the harness's resampler filter (§11.3, `nl_structure.resampler_gain`, no fitted numbers). Nulls per window from the input, each window from where the capture's content starts (§11.4) to the fade guard; windows within 10 dB of the noise are not scored. Floor = the capture noise (median of the LATE 0 tails) over the window's energy. Also per band over 0.05–0.6 s from the full-band capture and model, and an in-band (< 9.5 kHz) null.

### 11.1 Results

| Window (s) | n | Null L / R | In-band L | Noise floor | Excess | With the branch kernel (`NL_BK=1`, refit, α 0.7, §11.5) |
|---|---|---|---|---|---|---|
| 0–0.05 | 1036 | **−52.9 / −53.0** | −53.6 | −78.8 | 24.5 | −53.0 / −53.3 (excess 24.3) |
| 0.05–0.15 | 1026 | **−55.5 / −55.3** | −56.1 | −70.7 | 14.9 | **−58.7 / −58.8** (11.8) |
| 0.15–0.3 | 983 | **−55.1 / −54.6** | −55.6 | −67.9 | 13.4 | **−57.8 / −57.6** (10.5) |
| 0.3–0.6 | 904 | **−51.4 / −51.4** | −51.5 | −60.5 | 11.4 | −52.2 / −52.1 (9.2) |
| 0.6–1.0 | 645 | −41.9 / −41.7 | −41.9 | −48.6 | 6.9 | −42.0 / −41.8 (6.5) |
| 1.0–2.0 | 344 | −31.2 / −31.0 | −31.3 | −36.6 | 5.6 | −31.3 / −31.1 (5.6) |

(27 October: with the long-kernel clock warp. On the old warp, `bench/` and `bench_bk_a0.65/`: −52.6 / −54.8 / −54.2 / −50.9 dB without the kernel and −52.8 / −58.4 / −57.3 / −52.0 with it at 0–0.05 / 0.05–0.15 / 0.15–0.3 / 0.3–0.6 s, 10–12 kHz −34.5.)

Per band, 0.05–0.6 s (L / R), no kernel: 20–100 Hz −40.1 / −39.1; 100 Hz–1 kHz −54.4 / −56.5; 1–5 kHz −56.5 / −56.3; 5–10 kHz −53.6 / −52.9; **10–12 kHz −40.6 / −39.3** (old warp −34.5 / −34.1). With the kernel: −40.2 / −39.1, −55.2 / −56.9, **−59.2 / −59.0, −56.9 / −56.4**, −40.8 / −39.9. By group (L, 0.05–0.15 / 0.15–0.3 / 0.3–0.6 s, without the kernel): EARLY > 0 −55.2 / −53.9 / −44.8; EARLY < 0 −55.8 / −55.7 / −54.7; EARLY 10 −56.1 / −55.4 / −53.5. Trimmed-start captures (202) score like the rest.

- **The −50 dB target is met in every window up to 0.6 s** (median, both channels), held out, with 5 numbers per capture.
- **Repeat-capture level (the ceiling), per band**, from captures that differ only in LATE (identical early section, SPAN 0.5, 0.05–0.45 s, alignment fitted in band; measured on the old warp): 100 Hz–1 kHz −66.5, 1–5 kHz −60.6, 5–10 kHz −59.4, 10–12 kHz −58.9 dB; the model on the same captures −57.9 / −56.3 / −53.6 / −33.7, with the first kernel −58.4 / −61.3 / −59.2 / −33.8. A second, quick estimate on the new warp (`analysis/nl_hfdrift/repeat_ceiling.py`: 20 pairs of adjacent LATE, SPAN 0.5, EARLY 10 / 0.5 / 1 / −1, one delay and gain fitted on 8–11 kHz, so 10–12 kHz is held out of the fit) gives 10–12 kHz −53.2 / −48.9 dB and 11.8–12.8 kHz −46.6 / −46.3 dB, the same on the old warp (the warp artefact is common to repeats); the −58.9 figure was not reproduced, and against the quick estimate the model at 10–12 kHz (−40.8 / −39.9) is 8–12 dB short, not 25.
- **Positive EARLY at 0.3–0.6 s (−44.7 dB)** is data-limited: there the window is dominated by the late section at the law's −40…−75 dB, 9–20 dB above the noise; the late's shape is right (correlation 0.83–1.00) and its level matches the law within 0.5 dB wherever the late is ≥ 25 dB above the noise (§10.1).
- **After 0.6 s** the excess is 5–7 dB with the late section's timing and level exact (local shift 0.000 samples, gain within 0.02 dB): noise and Tomkin's fade.
- **Defects:** two captures (§9.7), one found here: (0.1, −1.0, 0.0) carries a foreign −52 dB tail (ACF 537) in a gated capture.

### 11.2 What the first bench (§11 of 26 October) got wrong, and the fixes (measurement side, not the model)
- **0–50 ms (−42 → −52.6 dB):** the front-tap position from the peak picker is biased by −0.0055 samples (its neighbours pull it), the same in every capture and channel, in both analysis bands; corrected as a constant (`nl_early.FRONT_BIAS`). A free per-capture front shift adds only 0.3 dB (median) and is not used.
- **10–12 kHz (−8 → −34.5 dB):** the deconvolution's band window (11.0–12.3 kHz raised cosine in the harness) and its regularised inverse shape every cached IR; the model now carries the same zero-phase response (`nl_structure.processing_response`, high and low band edges). The band window is widened to 12.0–12.8 kHz for N (`NL_WIDE`, default; own cache), since the unit's chain is only 4 dB down at 12 kHz and 26 dB at 12.5 kHz.
- **Low end:** with the deconvolution's 3–12 Hz edge matched, level and phase agree within ≈ 0.5 dB / a few degrees from 2 to 100 Hz; no extra corner is needed.

### 11.3 The top-end roll-off is the harness's own resampler (S; replaces the fitted "top-end EQ" and its H "LOW PASS")
The first bench pooled a zero-phase EQ of −0.24 / −1.03 / −1.85 / −3.0 dB at 11.5 / 12 / 12.3 / 12.6 kHz (`analysis/nl/out_eq.json`, 8 captures, later 141) and took it for the unit's LOW PASS at its 12 kHz step. It is not the unit's. Every cached IR went through `core.to_device_grid`'s `resample_poly(25, 91)` (96 kHz → device rate), whose Kaiser(5) FIR has a 3.9 kHz transition band around the device Nyquist (13.19 kHz): gain **−0.22 / −0.90 / −1.50 / −2.33 / −3.65 dB at 11.5 / 12 / 12.25 / 12.5 / 12.8 kHz**. The EQ refitted after the clock warp was fixed (§11.6; 120 captures, the fit on 60 of them, 12 knots, pooled, L and R alike) is −0.18 / −0.86 / −1.36 / −1.92 / −2.56 dB: within 0.05 dB of the filter up to 12 kHz, +0.15 / +0.4 / +1.1 dB above it (`analysis/nl_hfdrift/eq_vs_resampler.py`). The old numbers were that filter plus the old warp's own droop. The model now multiplies the exact filter gain in (`nl_structure.resampler_gain`, **0 fitted numbers**). **The unit shows no 12 kHz roll-off of its own here**; the plugin's N output needs none, and the spec's `output_eq` says so.

### 11.4 Trimmed starts (S)
The large "branch offsets" of steep negative EARLY at long SPAN are Tomkin's trim (§9.4): their files start mid-response. Every tap inside the file is at its law gain with the full diffuser, the trimmed taps' tails are present after the cut, and before the content starts the file is silent. Scoring therefore starts where the capture's content starts (within 45 dB of its peak, less 20 samples) plus a 40-sample guard where a file is loud from its first sample; 202 captures are flagged `trimmed` and reported apart (they score like the rest).

### 11.5 The main taps' read shape (S measured / H cause; adopted, refit on the long-kernel warp, strength 0.7)
- Averaging the residual of captures that share an early section (correlation 0.98 across captures) shows **the same small shape on every main tap**: arrivals at −13…+12 samples around each tap, 59–73 dB down (largest at −3, ±1, +3, ±11), in both channels; the front taps do not have it.
- As a kernel on the whole branch (29 taps, deviations ≤ 0.0007, band-limited below 11 kHz so it cannot touch the band edge; fitted on 18 LATE 0 captures at SPAN 0.2 / 0.4, `nl_model --stage bkernel`; the bench used the first fit on 10 of them, which agrees to 1·10⁻⁵), held out: 0.05–0.15 s −54.8 → −59.5 dB, 0.15–0.3 s −54.0 → −58.1, 0.3–0.6 s −50.9 → −52.2 (medians), 1–5 kHz +3.7 dB, 5–10 kHz +4.8 dB, at the repeat-capture level.
- **At full strength it fails the strict rule:** 41 of 2,052 channel-windows at 0.05–0.15 s and 3 at 0.15–0.3 s get worse by more than 0.5 dB (worst +0.96 dB), all steep ramps (EARLY ±0.1…0.7) at SPAN > 0, mostly R; those 41 captures also null worse without the kernel (−51.7 dB against −54.9 dB for the rest at 0.05–0.15 s, R, noise floor −62 against −72 dB), so the kernel's shape is not the main error there. Training on steep settings too gives the same kernel; giving R's first-tap direct term the kernel changes nothing.
- **Adopted at strength 0.65 (user decision, 27 October 2026).** `NL_BK_ALPHA` scales the kernel's deviation from a unit tap (k = δ + α(k₁ − δ)); sweep on the 41 losers plus 61 captures sampled evenly: α 0.5 / 0.65 / 0.75 / 1.0 keep −2.9 / −3.5 / −3.9 / −4.1 dB of the 0.05–0.15 s gain and leave the worst loss at +0.27 / +0.43 / +0.56 / +0.96 dB (0 / 0 / 3 / 41 channel-windows over 0.5 dB). **Full bench at α = 0.65 (all 1,044 captures, `analysis/nl/bench_bk_a0.65/`, held out as before, 5 numbers per capture; the kernel's 29 taps per channel are shared):** against the baseline, no channel-window is worse than +0.43 dB at any window (0.05–0.15 s median −3.6 dB; 0.15–0.3 s −3.0, worst +0.16; 0.3–0.6 s −1.1, worst +0.03; 0–0.05 s −0.15, worst +0.23; later windows ≤ 0.14), so **the strict rule holds**. Medians L / R: 0–0.05 −52.8 / −53.1, 0.05–0.15 −58.4 / −58.4, 0.15–0.3 −57.3 / −57.1, 0.3–0.6 −52.0 / −52.0, 0.6–1.0 −42.0 / −41.8, 1.0–2.0 −31.2 / −31.0 dB; per band 100 Hz–1 kHz −55.2 / −56.8, 1–5 kHz −59.0 / −58.7, 5–10 kHz −56.4 / −55.8, 10–12 kHz −34.6 / −34.3. `NL_BK=1` selects it (`BK_ALPHA_DEFAULT` 0.65; `NL_BK_ALPHA=1` the fitted kernel, `0` off); the spec's `branch_read_kernel` carries the fitted taps and `strength`. Inaudible either way (−63 dB deviations). A gain-dependent or per-tap shape was not pursued: the remaining headroom is about 0.6 dB median.

- **Refit under the long-kernel warp (27 October 2026, `tools/nl_model.py --stage bkernel` with `NL_WARP2` default; `analysis/nl/branch_kernel_warp2.json`; the first kernel is `branch_kernel.json`, legacy `NL_WARP2=0`).** About a quarter of the first kernel was the old warp's interpolation error (its kernel is ±10 input samples long, the shape spans −13…+12): the refit has the same shape at 74 / 79 % of the amplitude (L / R; rms deviation 1.4 / 1.7·10⁻³ against 1.7 / 1.9·10⁻³), largest taps −3, ±1, +3 again. The rest is real. The warp fix alone, without any kernel, already takes the medians 0.05–0.15 / 0.15–0.3 / 0.3–0.6 s from −54.8 / −54.0 / −50.9 to −55.5 / −54.8 / −51.4 dB; the old 0.65 kernel on the new warp is too strong (strict rule fails, worst +1.56 dB against the 0.65 bench, 1,119 channel-windows over 0.5 dB at 0.05–0.15 s). **Full bench of the refit kernel at α = 0.7 (all 1,044 captures, `analysis/nl/bench_warp2_bk_kf_a0.7/`, held out, 5 numbers per capture, the kernel's 29 taps per channel shared):** medians L / R 0–0.05 −53.0 / −53.3, 0.05–0.15 −58.7 / −58.8, 0.15–0.3 −57.8 / −57.6, 0.3–0.6 −52.2 / −52.1, 0.6–1.0 −42.0 / −41.8, 1.0–2.0 −31.3 / −31.1 dB (excess over the noise floor 24.3 / 11.8 / 10.5 / 9.2 / 6.5 / 5.6 dB). **Against the previous adopted bench (`bench_bk_a0.65`) no channel-window is worse than +0.07 dB at 0–0.6 s (+0.15 at 0.6–1.0, +0.36 at 1.0–2.0, noise-limited), and no band from 20 Hz to 10 kHz is worse than +0.39 dB;** medians are 0.2–0.3 dB better at 0.05–0.6 s and 1–10 kHz 0.2–0.3 dB better. The refit kernel at α = 1 fails the strict rule as the first one did (128 channel-windows over +0.5 dB against `bench_bk_a0.65`, worst +0.98). **Against the same bench without the kernel** (the rule's original form) α = 0.7 loses up to +0.81 dB (22 channel-windows over 0.5 dB at 0.05–0.15 s, 2 at 0.15–0.3 s and 2 at 0.3–0.6 s, all steep ramps such as 0.2 / −0.5 / 2.0, R); α = 0.6 on the 151 worst-case-plus-sampled captures leaves 10, α = 0.5 none (worst +0.49), with less gain (`NL_BK_ALPHA=0.5`; screened on that subset, not benched in full). 0.7 is the default (`BK_ALPHA_DEFAULT`) because it is never worse than what was adopted; 0.5 is the setting that also holds against the no-kernel bench. 10–12 kHz is unchanged by the kernel (it is band-limited below 11 kHz).

### 11.6 Open, with what is known
- **11–12.3 kHz timing drift: found, a harness artefact (S; fixed 27 October 2026).** It grew with file time (0.01–0.15 samples at 0.1–0.9 s, 11–12.8 kHz only, same in repeats, early and late, both channels, SPAN / EARLY / LATE alike) because `nl_late.unwarp`, the per-capture clock correction, upsampled ×32 with scipy's default Kaiser(5) FIR (10 input samples per side). Its transition band, 3.9 kHz wide around the device Nyquist (gain −0.22 / −0.90 / −1.65 / −2.72 / −3.65 dB at 11.5 / 12 / 12.3 / 12.6 / 12.8 kHz), leaves the image of 11–13 kHz in the 32× signal, so resampling at the slowly advancing fractional position δ = ε·n gives `G(f) + G(fs − f)·e^{i2πδ}` instead of a pure delay: a phase error ∝ δ and an amplitude loss ∝ δ² at the band edge, turning round when δ passes about 0.2 (about 0.9 s at 12 ppm); below 11 kHz G(f) + G(fs − f) = 1 and nothing happens. Evidence: `unwarp` applied to the model alone, no capture, reproduces it (+0.079 samples at 12.3–12.8 kHz and +0.037 at 11.6–12.3 kHz at 0.86 s, against +0.09 / +0.04 measured; measured minus predicted is flat in time), the drift follows file time (trimmed files start near 0), and a long kernel (150 input samples per side, Kaiser β = 9, flat to 12.8 kHz) gives a common delay in every band to ±0.001 samples. Not the clock ratio's precision (refitting ppm on 8–11 kHz moves it by sd 0.7 ppm and gains ≤ 0.5 dB at 9–12 kHz), not aliasing, not the unit, not the capture chain. **Fixed** (`NL_WARP2`, default; `NL_WARP2=0` the old warp): 10–12 kHz −34.5 / −34.1 → **−40.6 / −39.3 dB** (median, 1,044 captures, no fitted EQ), 1–5 kHz 0.5 and 5–10 kHz 0.7 dB deeper (median, without the kernel), and no window worse than +0.1 dB. What is left at 11.6–12.8 kHz (−27 … −36 dB per sub-band) is a static phase of 0 / 4 / 11 / 17 / 28 / 40 mrad at 11.5 / 11.75 / 12 / 12.25 / 12.5 / 12.8 kHz (about 0.01 samples at the edge; H: the real chain against the chain kernel's regularised inverse, which dominates above about 12.1 kHz) and a 0.03–0.07 dB level step; a pooled complex residual of 22 numbers gains 0.5–2.3 dB at 10–12 kHz held out by capture halves but 30 % of the captures get worse by more than 0.5 dB (late-dominated settings), so it is **not kept** (`analysis/nl_hfdrift/resid_eq.py`). Report: `analysis/nl_hfdrift/REPORT.md`.
- **100 Hz–1 kHz, first 30 ms:** the SPAN 0.3 / 0.5 capture sessions (front shift ≈ 0) null at −49 dB there against −57…−58.5 for the others; after 30 ms all sessions are alike (−55 dB). Not the trim (a simulation of the cut through the chain and its inverse changes nothing), not a front / branch path difference (separate gains and low-frequency terms gain ≤ 0.7 dB), not timing or level per window (≤ 2 dB). H: capture-session-specific low-band behaviour at the file start. The low band stays 8 dB above the repeat level.
- **Late levels below about −55 dB** (steep positive EARLY): 3–5 dB below the law, at 9–20 dB SNR; a 16-bit recursive gain computation explains only ≈ 0.1 dB there. Not resolvable from the library.
- **Nothing in this library tests:** SPAN steps between 0.1 multiples (H rule §9.2), LATE beyond 2.0, DELAY, the LOW PASS below 12 kHz.
- **11–12.8 kHz after the warp fix: a term of the capture chain, per capture run (S; opt-in `NL_HF=1`, 27 October 2026, `analysis/nl_hfdrift/REPORT.md`).** Headroom: the repeat ceiling (576 LATE-only pairs, alignment fitted on 8–9 kHz so every sub-band is held out) is −56.5 / −56.0 / −52.4 / −40.4 dB at 9–11 / 11–11.6 / 11.6–12.3 / 12.3–12.8 kHz against the model's −52.6 / −42.2 / −36.0 / −28.4: 15–19 dB above 11 kHz, model-side, not noise (62–76 dB down) and not repeat-limited; the 10–12 kHz ceiling is −55.8 dB. What is left is the same in every time window after 0.05 s, SPAN, LATE, channel and trimming (the chain and its regularised inverse, not the model), and is measured from the front taps alone (R: −1.5 / −5.7 / −12.1 / −20.8 / −30.3 / −40.0 mrad at 11.5 / 11.75 / 12 / 12.25 / 12.5 / 12.8 kHz, +0.03–0.09 dB). It follows the **capture run**, not SPAN: the WAV `bext` chunks carry the origination times (1,038 of 1,046 files, 9–12 April 2016), which split the library into 12 runs; one scalar per run on a rank-1 shape of 9 complex knots explains 98 % of the run-to-run variance (a band-edge filter shift: not a delay, the best per-run delay leaves 28 % of the phase variance; not the clock, ppm correlates +0.5 across runs and with changing sign within them). Re-fitting o, the front shift and ppm on a high band fails: it costs 5–10 kHz 4–9 dB for at most 2 dB at 11–12.8 kHz; the bench's alignment is optimal to 1·10⁻⁴ samples there, and the +0.004-sample high-band delay is the same on the front taps and the branch (a chain phase). **`NL_HF=1`** (`analysis/nl/hf_corr.json`, built by `analysis/nl_hfdrift/f_table.py`) multiplies the model by exp(z) above 10 kHz, z = shared mean + scalar × shape at 9 knots (10.5–12.8 kHz; shape fitted on the other runs, scalar from the R front taps of the 15 captures nearest in time, the capture's own front taps excluded; 36 shared numbers + one scalar per capture run; a measurement-side term like the clock and the gain, not the unit). Full 1,044-capture bench against the current default (`bench_warp2_bk_kf_a0.7_lb2_hf` vs `…_lb2`): **10–12 kHz −40.8 / −39.9 → −48.1 / −48.4 dB (median)**; every window and every band 20 Hz–10 kHz unchanged (worst channel-window +0.49 dB, worst band +0.20 dB). **But 320 of 2,088 channel-bands at 10–12 kHz get worse by more than 0.5 dB (worst +13.9 dB)**, mostly captures that already nulled well (−47 → −42 dB), clustered in SPAN 0.3 (115 / 332), SPAN 0.1 (86 / 350), the trimmed captures (123 / 404, gain only 1.2 dB) and the (0.2, −0.2, *) family; so it is left **opt-in, default off**. A pooled single correction (22 numbers from the reverb, or 20 from the front taps) gains only 1–4 dB with 22–30 % of captures worse.

### 11.7 Second capture: the Altiverb N presets (held out)
The Altiverb set (a third-party capture of the same design) holds two N presets; its `NON-LIN SHORT RT` files are S (S's cluster, no N taps). **NON-LIN = SPAN 0.3, EARLY 10, LATE 0.4; PRESET54 GATED = SPAN 0.2, EARLY −3, LATE 0.8** (stage count, stage-gain slope and tail, best fit over the library steps). Front taps at 41 / 410 / 463 (R 278) with 1 : 2 : 1 and every main tap at the table position + 36 (the Altiverb onset), clock +43 ppm (its own, as for E). Nulls with alignment, clock, gain and a 7-tap chain FIR per channel fitted on 0.05–0.15 s, held out on the rest: L −32…−34, R −29…−31 dB, by band −43…−46 dB at 1–6 kHz and −18 dB at 100–250 Hz / −36 dB at 8 kHz: the Altiverb chain's known weak spots (S reached −39…−41 there with a 29-number chain model, `model_s_results.md` §37). **The structure is confirmed on a second unit's capture.**

### 11.8 The low band: a measured transfer of the captures' processing (27 October 2026; S measured / H cause; **adopted as the default, `NL_LOWBAND=2`**)

Full record: `analysis/nl_lowband30/REPORT.md` (parts 1–6). Switch: `tools/lowband_joint.py`, wired into `nl_structure.processing_response`; `NL_LOWBAND=0` off, `1` the shared table, `2` (default) the shared table below 100 Hz plus N's own table above 130 Hz (blend 100–130 Hz). Bench-side only: it is a property of the captures' processing, not of the unit, and does not enter the plugin or the spec.
- **20–100 Hz was a model error, not capture-limited.** The repeat-capture level (SPAN 0.5, captures differing only in LATE, window 0.05–0.45 s, |Δo| < 0.02, same method as §11.1: pooled over EARLY 20–100 Hz −57.1, 100 Hz–1 kHz −66.0, 1–5 kHz −59.1, 5–10 kHz −57.9, 10–12 kHz −56.1 dB; with EARLY ≥ 0.7 −62.0 / −68.2 / −61.7 / −60.6 / −59.4) sits 18–27 dB below the model's −39 dB at 20–100 Hz and 9–13 dB below it at 100 Hz–1 kHz. Repeat levels depend on the EARLY ramp (a decaying early section is noise-limited: EARLY 0.2–0.3 gives −44…−51 dB at 100 Hz–1 kHz).
- **The error is a fixed linear transfer** T(f) = data/model (complex, per channel): identical for SPAN 0.2–0.5, all sessions, L and R, EARLY 0.3–10 and every file length; 0.9 % rms at 20–100 Hz falling to 0.06 % at 1–3 kHz, fine structure (±0.05 dB, ±0.5° with a period of about 25 Hz). Not stage gains, all-pass coefficients, tap gains or file ends (each tested). Below 100 Hz **the same transfer fits Program S** (joint N+S table within 0–1.7 dB of each program's own, held out on both: N by SPAN, S by RT fold; correlation of T_S and T_N 0.9 at 20–40 Hz): the chain kernel's / processing response's low-band error (H). Above 100 Hz it is program-specific (S: its loop-path transfer, `model_s_results.md` 41.3).
- **Result (bench, 1,044 captures, alpha 0.7 branch kernel, medians L / R):** 20–100 Hz −39.5 → **−48.9 / −48.0**, 100 Hz–1 kHz −55.8 → **−61.5 / −60.2**, 1–5 kHz −59.1 → −59.7 / −59.2, windows 0.05–0.15 s −58.7 → −59.7 / −59.2, 0.15–0.3 s −57.6 → −58.3 / −57.9. Held out (983 captures the 61-capture estimation set never contained): 20–100 Hz −39.6 → −48.6 (pooled L/R medians; the 61 estimation captures −38.9 → −47.9), 100 Hz–1 kHz −55.9 → −61.2. Strict rule: 1 of 1,966 channel-windows at 0.15–0.3 s worse by +0.52 dB (limit 0.5; the user accepted it), the next worst +0.32 dB. The shared table has no extra smoothing (2 Hz bins; ±2-bin smoothing was 1.4 dB worse at 20–40 Hz on N, S unchanged; the first version of the table used it and gave 20–100 Hz −46.8 / −45.9 with the same strict-rule result). Mode 1 alone (that first version): 20–100 Hz −46.8, nothing worse than +0.03 dB.
- **S:** `SML_LOWBAND=1` (off by default) in `sml_fullband.full_model`: full-band gate PASS, S 20–40 Hz +3.7…+12.1 dB, 40–80 Hz +1.8…+13.2 dB, 80–160 Hz ≤ +0.7 dB, nothing above 160 Hz changes (partly in-sample for S; held out 20–40 Hz −35.3 → −39.6).
- **Not closed:** 20–100 Hz is still 1–2 dB short of −50 dB (the file-start step of the next bullet is worth about 1 dB there) and 3–5 dB above the repeat level; no low-order form (biquad, edge, high-pass, rational ≤ 8) reproduces the table (rational order 10 does); the sub-18 Hz region (session-dependent) is not corrected.
- **The first 30 ms (§11.6 bullet 2) is the pipeline's response to a DC offset in each file's pre-roll (S; mechanism found 27 October).** A baseline step in the residual at cache index 45–46 (four device samples before file sample 0), flat before it and about −90 dB after, 0.4–8·10⁻⁴ of the 410 tap (R larger). A forward replica of the whole pipeline (model → clock warp → 96 kHz → chain kernel → file length → `core.deconvolve` → `to_device_grid`, which reproduces the cache to 1.3·10⁻⁸) gives no step from the model's finite record (file end, fade, PAD, FFT size, kernel high-pass tails: below 10⁻⁴), but feeding the real files' own pre-roll through it reproduces the step exactly (L −0.46 / R −1.46 against measured −0.44 / −1.45). The files carry a DC of 0–4.5·10⁻⁵ of full scale (plus a Nyquist tone) before the first arrival; the step is that DC through the deconvolution's low edge (unit-DC response about −9.6 before index 45, +2.4 after, fading over about 0.1 s). c = −9 · DC/peak with DC the mean of the file's first 30 samples: correlation −0.98 / −0.98 (L / R) over 86 captures, residual 12–17 % of c's rms. The R excess, the session differences and the growth with LATE are changes of the file's DC (a converter or session offset, H), not of the unit. Adding the DC's pipeline response to the model (no fitted parameter; two numbers read from each file) takes 20–100 Hz in 0–50 ms −7 dB (−12 dB with NL_LOWBAND=2) and 100 Hz–1 kHz −2…−2.5 dB, but gains nothing after 50 ms and worsens 20–100 Hz over 0.05–0.6 s in some captures (worst +3.75 dB): a correction for scoring the first window, not adopted in the bench (`analysis/nl_lowband30/cand_pipeline/`).
- **Candidates tried and not adopted (27 October, `analysis/nl_lowband30/cand_*/`):** a 38-number kernel-tail form (cubic B-spline correction to the chain kernel over 0–100 ms, knots every 6 ms) reproduces the 2 Hz shared table within 0.9 dB (0.16 dB mean) at 20–80 Hz held out, but is no more compact than a rational order-10 fit and does not carry to 130–1000 Hz; the energy of the kernel error sits at 25–40 ms in the time domain, past the 38.9 ms FIR only in R. After the shipped transfer, 100 Hz–1 kHz is −63.7 dB against a like-for-like repeat level of −66.7 (the bench median −61.5 is noise-limited: a perfect fix would gain 2 dB at most), with 126 derivative columns per channel (tap, stage and front gains and timings, all-pass coefficients, branch kernel, read shape, tap signs) gaining at most 1.6 dB held out and a pooled fixed transfer at any resolution −2…+3 dB; the 20–100 Hz remainder is an EARLY-dependent low-frequency model term (sub-18 Hz transfer varies with the EARLY class, not with session or file length; a pooled sub-18 Hz table gains about 1 dB and R never reaches −50 dB held out).
- **Remaining ideas, tested 28 October and left out:** (1) an EARLY-class transfer (EARLY ≤ 1 against ≥ 1.5, held out by SPAN, with S in the pool): 0.1–0.4 dB over the pooled table (`analysis/nl_lowband30/lb_earlyclass.py`), so the EARLY dependence is not a transfer; (2) a pooled sub-18 Hz table (about 1 dB at 20–100 Hz, R never −50 dB; not shipped, see above); (3) the front-tap terms (3–6 numbers, −1.2…−2.2 dB held out but unstable, ill-conditioned: not tried further); (4) the 38-number kernel-tail form and the DC scoring term (above). Nothing is worth adopting; 20–100 Hz stays at −48.9 / −48.0 dB (1–2 dB short of the −50 dB bar).
- **S and M / L (28 October).** The same table goes into S's full-band model (`SML_LOWBAND=1`) and M / L's (`ML_LOWBAND=1`, per channel on the pipeline response; `sml_fullband.py --stage ml`); the full-band gate (`--stage gate`) now runs S and M / L with both switches forced on and `gate_ref` was rewritten with them (S medians as in the bullet above; M / L change by −0.9…+0.3 dB, the best at 20–40 Hz: L cluster/gap −0.9, first pass −0.8, R/L pooled M −0.1…−0.4; they are at their time-variant repeat level there, so the table is neutral, not a gain). The exploratory stages keep the switches off by default, because the fold-transfer machinery (`fold_T`, `analysis/fullband/s/*.npz`) needs the uncorrected sums.

## 12. For the plugin

`tools/nl_spec.py` → `analysis/nl/nl_plugin_spec.json` (every constant with label and section) and `analysis/nl/ref/nl_<SPAN>_<EARLY>_<LATE>.wav` (device rate, float32, input impulse at sample 0, main-tap units). `simulate()` is the sample-by-sample reference (delay lines, all-passes, loops); it equals the bench model (device only, `render(proc=False)`) to −136 dB or better.

- **Signal flow (per channel):** front taps on the input; the input through AP(124) into a delay line read by the main taps (stage gains from the EARLY law), summed into AP(180) L / AP(188) R; the same AP(124) line at FEED[SPAN] (R 220 earlier) through AP(314) into three flat loops (Jot gain for RT = LATE) read by four equal taps each, at the late gain; LATE 0 drops the late. The branch sits **one sample** after the table positions (the captures' modal offset; then R's first tap needs no special case).
- **Cost:** per channel 3 front + up to 27 main taps, 4 all-passes (AP124 shared), 3 loops with 12 reads: well inside 2 × TMS32010.
- **Controls (N's panel slots):** DELAY (H: a pure pre-delay); SPAN in the AMPL slot, library values 0–0.5 in 0.1 (finer steps H by the stage rule); EARLY in the LF slot, **16 values −3, −2, −1, −0.7, −0.5, −0.2, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 10 (flat)** = the LF column's 16 steps; RT slot blank; LATE in the HF slot, **0–2.0 in 0.2 (11 values), 0 = gated**; LOW PASS as the other programs (its 12 kHz step most likely the §11.3 EQ). The plugin's current N tables (`plugin/src/params/Parameters.cpp`) are guesses and differ: to be changed in the plugin session.
- **Optional:** the branch read kernel (adopted at strength 0.7, §11.5; −63 dB, inaudible). No top-end EQ: the captures' 12 kHz roll-off is the harness's resampler (§11.3). The plugin's `CaptureEq` is the separate low-band transfer of `model_s_results.md` §27, and the plugin does not read the spec's `output_eq`.
