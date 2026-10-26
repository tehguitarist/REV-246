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

## 8. Open questions and discovery plan

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
- The AP(124) branch (main taps **and** the late section) sits o samples from the front taps, o per capture: 0 in most, else −63…−64, −84…−87 or −163, mixed at random among captures of one (SPAN, EARLY) that differ only in LATE. For steep negative EARLY at long SPAN (−0.2…−0.7 at SPAN 0.3–0.5) o reaches −486 … −9,316 samples; the train keeps its exact shape and law gains, and the stages that would land before the file start are missing (Tomkin's trim cut them).
- R's first-tap direct term (§9.2) does not move with o.
- **H:** the front taps and the branch run on different processors joined by a buffer whose pointer phase is set at program load (2 × TMS32010). For the model, o is one alignment per capture; the plugin uses o = 0.

### 9.5 Held-out score (`nl_early.py --stage law`, `analysis/nl/early_law.json`)
Per capture only the alignment (front shift, o, clock) and one gain per channel; every structural number is shared and no setting's gains are fitted. On the 96 LATE 0 captures, 100 Hz–10 kHz, from 20 samples before the first front tap to 3000 samples past the last main tap: **median L −45.3, R −45.1 dB** (the free per-tap fit, ≈ 31 gains per channel per file, in-sample: −44.2 / −44.4). Clock ratios +1…+14 ppm.
- **Floor:** captures that share their early section (settings differing only in LATE) agree to −52…−62 dB, so ≈ 10 dB of model error is left. The residual is dense, HF-heavy (70 % of it at 6–10 kHz while the data is white), follows the signal's envelope, repeats across captures (correlation 0.99 in L) and is not a linear function of the main taps (a shared 63-tap output kernel, per-tap 5-tap shapes, all-pass damping or gain changes each gain ≤ 1 dB; a 1,300-tap shared kernel 3–4 dB held out). Local timing errors of ±0.003 samples would explain its size and spectrum (H: sample-level timing in the unit's output); open.
- Steep negative EARLY at long SPAN puts the front taps below the noise (−60 dB re the ramp's end), so their windows are noise-only.

### 9.6 Model size
Structure: 27 tap positions, 2 × 27 signs, 3 front positions per channel with 1 : 2 : 1 amplitudes, R's offset (−221) and first-tap direct position (465), three all-pass delays at one gain, 14 stages and the SPAN → stage table; law: 4 ms tick, 36 ms stage, front ticks, the 16 ms / 12 ms offsets. Per sample: 3 + 27 taps per channel and 3 all-passes (≈ 36 multiplies per channel), within the 2 × TMS32010 budget.

### 9.7 Defects
`defects.NL_DEFECTS`: (SPAN 0, EARLY −0.5, LATE 0) carries a 0.27-sample slip between the front taps 410 and 463 (−15.7 dB against the model; every other file with visible front taps holds 369.00 / 53.00 samples to ±0.2).

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
