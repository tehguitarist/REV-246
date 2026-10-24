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
