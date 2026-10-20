# EMT 246 reverb: discovery and viability study

> **Status note (Oct 2026):** this is the original study, with corrections added inline. Current findings live in `EMT246_unit_facts.md`, and current status and next steps in `NEXT_STEPS.md`. Where they differ, those two win.

*Scope: programs S/M/L (A/B/C in the brochure): RT, LF and HF decay, reflection cluster (delay and amplitude), and lowpass. Program E (the EMT 250 algorithm) and program D (delay/echo) are out of scope. Program N (Non-Lin) is an open scope question; see §9.*

*Confidence labels used throughout:* **[confirmed]** means measured or printed in a primary source. **[strong]** means several independent measurements agree. **[hypothesis]** means it is consistent with the data but needs more IRs or documents. Every number from the single IR should be re-checked once more captures are available.

---

## 1. Verdict

**Viable, and more viable than expected.** Analysing the one IR you supplied turned up the following:

1. **The unit runs at about 26,373.6 Hz, which fits 19.2 MHz ÷ 728.** It does not run at exactly the 26.4 kHz the brochure states. Once the converter filter is deconvolved, the cluster taps land on exact integer sample positions on this grid **[strong]**. That means every delay inside the algorithm can be read off as an integer.
2. **The reverb is essentially time-invariant [strong for R, probable for L].** The spectral fine structure stays the same across the decay in both channels.
   - *Caveat added later:* a 1986 review describes "stochastic control" in S/M/L (§10a). Any randomisation must therefore be slight and must leave the loop modes harmonic. See algorithm_research §7.2b.
   - The person who made this IR (eugtone) null-tested it against the hardware. They reported that R "nulls the HW down to quantization noise". L had about 50× more error, and they suspected "a very subtle modulation on the left channel".
   - This IR shows no L-only modulation: mode-phase stability is the same in both channels. However, L's noise floor is 6 dB higher (−104 against −110 dB).
   - For practical purposes an IR captures nearly everything. A possible slight L-only time variation can't be checked with repeat captures (none are possible). It can only be checked by joint multi-IR fitting or separate-capture duplicates in libraries.
3. **The late reverb in program S is four parallel recirculating loops** of 1321, 1389, 1461 and 1587 samples (50.1 / 52.7 / 55.4 / 60.2 ms) **[strong]**. The evidence:
   - The tail's autocorrelation has peaks at exactly these lags and at their multiples, with nothing else above the noise.
   - The spectral modes sit on these four harmonic series, about 22 dB above random frequencies.
   - A model made only of four periodic components explains 98.7 % of the tail's variance below 1 kHz and 97 % below 3 kHz. A control set of four similar but wrong lengths explains 36 %.
4. **The "Reflexionscluster" is 8 single-sample taps, four per channel, at fixed integer positions** **[strong]**:
   - L: 0 (+), 54 (+), 194 (+), 288 (−)
   - R: 0 (+), 143 (+), 194 (−), 353 (+)
5. **The left output lags the right by half a sample (18.6 µs) [strong].** The two channels otherwise have identical analogue chains. The cause is **one DAC multiplexed into two S&H outputs** **[confirmed: printed block diagram, and one PCM53 on the output board]**.

This is a small, deterministic, integer-delay architecture. That is the best case for grey-box reverse engineering: a structural model whose topology and integer constants are pinned down from IRs, with the parameter laws fitted from IRs swept across settings. The main risk is **getting enough well-labelled IRs**, not the analysis itself.

No patent covers the 246 algorithm (§2.4), and no manual is publicly available (§10), so the IRs are the main source of truth.

There is a large commercial IR set: Tomkin Mastering's library of more than 2,400 settings across Short, Med, Long and Non-Lin (§10). That makes the data problem far more tractable than first thought.

---

## 2. What the unit is (from your brochure, manual photos and parameter card)

### 2.1 Specifications relevant to modelling [confirmed, brochure *Technische Daten*]

| Item | Value | Modelling implication |
|---|---|---|
| Programs | S/M/L (= A/B/C), E (= EMT 250), D (delay/echo), N (Non-Lin) | S/M/L are the "new algorithm" with "adapted density" for each program, so expect different loop sets per program |
| RT @ 1 kHz | S 0.2–2.5 s, M 0.5–5 s, L 2–20 s, **16 steps** each | Discrete, so each step can be sampled exhaustively |
| LF factor @ 50 Hz | ×0.1 … ×2, **16 steps** | Loop filter gain law at LF |
| HF factor @ 7 kHz | ×0.1 … ×1, **10 steps** | Loop filter gain law at HF |
| Pre-delay | 0–999 ms, 1 ms steps | Integer delay, which is trivial to model |
| Reflection cluster | **8 reflections, fixed pattern**, amplitude 0–100 % in **11 steps** | Pattern measured (§4.2) |
| Lowpass | 2–12 kHz, **8 steps**, "proportional control of cluster and reverb" | Applies to both paths |
| Converters | 16-bit linear SAR, quoted fs 26.4 kHz | Measured as 26,373.6 Hz (§4.1) |
| DSP | 2 × TMS32010NL (confirmed from photos, §2.3a), 16/32-bit, 208 ns per instruction, 19.2 MHz (38.4 MHz ÷ 2) | 182 instructions per sample per DSP (§2.3) |
| Memory | 128 kB RAM (64 k words = 2.49 s at fs), 12 kB ROM | Total delay-memory budget |
| Internal arithmetic | 32-bit accumulate, **16-bit storage of intermediate results** | Quantisation in delay lines, worth emulating for grain/noise |
| Frequency response | 30 Hz–12 kHz | Matches the measured band-limit |
| Reverb S/N | 73 dB unweighted | Sets the noise floor to emulate |
| Presets | 90 factory (F1–F90), 99 user, F91–F99 test routines | Preset table lives in the remote (see below). **Out of scope for now:** preset settings are unknown |

### 2.2 Front panel and remote (manual photos pp. 65–74, plus panel photos)

**Panel layout (from panel photos):**
- **PROGRAM:** keys S, M, L (red LEDs), D (green), E (red), N (yellow).
- **MEMORY:** keypad 0–9, −, +, FIX/USER, RCL (EXC), STO. An EMT 246 M memory-module slot sits on the left.
- **Input level meter:** REG, +3, 0, −3, −6, −9, −12 dB.
- **Displays and controls (left to right):**
  - an 8-character program name (e.g. "SHORT RT"), then the memory display;
  - **LOWPASS:** one parameter, `f (kHz)`;
  - **REFLECTIONS:** two parameters, `DELAY (ms)` and `AMPL. (%)`;
  - **DECAY:** three parameters, `LF×`, `RT (s)` and `HF×`;
  - each section has its own LED wedge bar per parameter and **one** trough controller. The "|", "|–|" and "|–|–|" symbols show how many parameters that controller addresses.
- **Legend colour = program LED colour [strong]:**

  | Legends | Colour | Program LEDs | Programs |
  |---|---|---|---|
  | `f (kHz)`, `DELAY (ms)`, `AMPL. (%)`, `LF×`, `RT (s)`, `HF×` (below the displays) | red/orange | S, M, L, E are red | reverb programs |
  | `L — DELAY (ms) — R`, `L — ECHO (%) — R` (above) | green | D is green | Delay/Echo |
  | `DELAY (ms)`, `SPAN (s)`, `EARLY — RT (s) — LATE` (above) | yellow/cream | N is yellow | Non-Lin |

  The parameter card groups the controls the same way: "Progr.: S, M, L, E", "Progr.: D" and "Progr.: N" rows.
  - **Program N therefore has four parameters** [inference]: DELAY and SPAN (0.1–0.5 s, "Anhalldauer") on the REFLECTIONS controller, plus EARLY RT and LATE RT on the DECAY controller.
  - That probably leaves the middle DECAY display unused, as D likely does for its two-sided legends.
  - N appears to have **no cluster amplitude and no LF/HF factors**.
  - Whether LOWPASS (red legend only) also acts in D and N is unknown.

**Block diagram (printed on the mainframe front) [confirmed]:**
- **Input:** a single (mono) balanced input → amp → anti-alias low-pass → S&H → A/D.
- **Digital section:** CPU, ROM, RAM.
- **Output:** **one D/A feeding two separate S&H stages** (L and R), each followed by a reconstruction low-pass and an output amp.

  This confirms that the measured 18.6 µs (half-sample) L/R offset is a property of the design. An emulation should reproduce it.

**Remote electronics (manual photos):**

- **The remote (EMT 246 S) is a complete Z80 computer.** It contains:
  - a Z80A CPU (Z601) and a 27C256 EPROM labelled "BGR 1A" (firmware revision);
  - 6116 RAM;
  - a Z80 CTC (Z603), plus an optional MIDI interface on J607;
  - an ICM7218B 8-digit LED display driver (Z612), driving the display board 7 246 124;
  - an MK48Z02 battery-backed RAM on the plug-in memory module (EMT 246 M, board 7 246 128).

  So the **factory preset table and the parameter-to-display mapping live in the remote's 27C256**. An EPROM dump would give exact preset values and the display scaling, but not the DSP algorithm.
- **The controls are incremental, not potentiometers.** The "trough-shaped slide elements" are read by opto-switch pairs (H601–H610). Those go through LM339 comparators and CD4066 switches into 74HCT244/LS244 buffers, which is the classic quadrature-encoder arrangement. Combined with the step counts above, **every parameter is a small integer index**. There are no continuous values to interpolate, and a plugin can reproduce the stepping exactly.
- **Keyboard (board 7 246 127):**
  - program keys S, D, E, M, L, N;
  - 0–9, −, +, STO, F/U and RCL;
  - an 8-digit, 7-segment display, plus LED bar "analogue light columns" for each parameter.
- **Link to the mainframe:** 5-pin DIN carrying +25 V, GND, DATA, R/W and ACT/IN. The data runs as 32-bit pulse packets (one pulse = 0, three pulses = 1), with gaps that let the mainframe reply. *The DSP program is in the mainframe's 12 kB ROM, not in the remote.*
- **Parameter card wording** (from the studiorepair.com photo):
  - Programs S/M/L/E share the REFLECTIONS controls "Verzögerung (ms) / Amplitude (%)" and the DECAY controls "tiefe Freq. / Nachhall / hohe Freq.".
  - The card's "optimum RT range" bars are: S 0.2–1 s, M 0.8–4 s, L 3–10 s, with L extending to 20 s.
  - The diagram puts the low-pass on both cluster and reverb, with the reverb starting after the cluster.

### 2.3 Compute budget: why the algorithm must be small [hypothesis, derived]

- 19.2 MHz ÷ 4 gives a 4.8 MHz instruction rate, and 4.8 MHz ÷ 26,373.6 Hz is **exactly 182 instruction cycles per sample per DSP**. That is 364 across both DSPs.
- The TMS32010 has only 144 words of on-chip data RAM. Delay memory must therefore be reached through I/O or table instructions at 2–3 cycles each.
- Each delay tap (address out, data in, multiply, accumulate) costs roughly 5–8 cycles. That leaves room for a few dozen delay-memory accesses per sample across both DSPs, for:
  - 8 cluster taps;
  - 4 loops with their damping filters;
  - input diffusion;
  - L/R output taps;
  - the lowpass;
  - the pre-delay.

  The four-loop finding fits this budget. A big FDN with dense matrices would not.
- The exact integer ratio (728 = 4 × 182) is good circumstantial support for the 19.2 MHz ÷ 728 clock.

### 2.3a Mainframe internals (studiorepair.com photo gallery, October 2026)

Source: http://www.studiorepair.com/gallery/EMT/EMT_246/index.html (Nikolaus Riehm; 822×622 photos). These are part numbers read off the photos, so treat small labels as **[strong]** rather than certain.

**Processor board 7 246 123:**

| Item | Parts | Notes |
|---|---|---|
| DSPs | **2 × TMS32010NL** (TI, date codes 8540 and 8618) | **Confirmed: TMS32010, not TMS32020** |
| Clock | **MCO-5100 38.400 MHz** oscillator (Musashino) | The TMS32010 tops out at 20.5 MHz, so this must be halved to the brochure's 19.2 MHz. Its instruction cycle is CLKIN/4 = 208.3 ns, which gives **182 cycles per sample per DSP at fs = 38.4 MHz / 1456 = 26,373.6 Hz**. That matches the measured rate exactly |
| DSP program memory | **6 × AM27S191 bipolar PROM** (2 K × 8), labelled **BGC 1A–6A** | 6 × 2 KB = **12 KB, matching the brochure's ROM**. BGC 1–4 sit by the left DSP (4 K × 16 words, the TMS32010's maximum program space); BGC 5–6 by the right one (2 K words) |
| Empty PROM sockets | **2 empty 24-pin sockets** beside BGC 5/6 | Room for 2 K more words on the second DSP. Possibly the optional **"EMT 250 program"** (a 2020 listing mentions an "EMT 250 Option") **[hypothesis]** |
| Delay memory | **4 × NEC μPD41464C-12** (64 K × 4 DRAM, 120 ns) | = **64 K × 16-bit words, 128 KB, matching the brochure**. One 16-bit word per access cycle |
| Spare DRAM | One Fujitsu MB81464-12 in a position silk-screened "Spare" | Seven unpopulated 18-pin footprints sit along the top edge |
| Other RAM | 2 × UM6116-3L (2 K × 8 SRAM) | Likely the coefficient/parameter RAM written by the control Z80, or inter-DSP data **[hypothesis]** |
| Glue logic | Small bipolar PROMs / PALs labelled **BGP 1A, BGP 2A, BGD 1A**; 74S188 (32 × 8) PROMs; 74LS193 counters; 74LS374/245 latches and buffers; 2 × Am2966 DRAM drivers | Hardware address generation for the delay RAM **[hypothesis]**. The BGP/BGD parts may hold delay-offset tables. That would fit fixed loop lengths per program |

**Interface board 7 246 125 (in the mainframe):**
- Z80A CPU, Z80 CTC, Z80 PIO;
- EPROM **"BGX 1B"**;
- UM6116 SRAM;
- 74LS374 latches to the processor board.

This is the brochure's control processor. **The step-index → coefficient tables (RT → loop gains, LF/HF → filter coefficients, and so on) most likely live in BGX 1B** **[hypothesis]**.

**Analog input board 7 246 121:**
- An input amp and an anti-alias filter on a daughterboard (TL071/TL074 active filter);
- an S&H, then a **discrete 16-bit SAR ADC** built from a **PCM53KP-I** DAC, **2 × 74LS503 successive-approximation registers** and an LT1016 fast comparator. This matches the brochure's "successive approximation, 16 bit".
- Trimmers: LEVEL, HF CMR, LF CMR, OFFSET.

**Analog output board 7 246 122:**
- **One PCM53KP-I DAC**;
- **LF13333 analog switches** feeding **two hold stages**;
- **two identical active-filter daughterboards** (TL074-based, high order), one per channel;
- trimmers LEVEL 1 and LEVEL 2.

**This confirms in hardware that one DAC is multiplexed into two S&H channels**, which is the source of the measured 18.6 µs (half-sample) L/R offset.

**Other details:**
- The front badge reads **"BARCO-EMT"**, so units were still built after Barco's 1989 purchase.
- The inside of the front flap carries a printed board map with trimmer names. It is not a schematic.
- There is no separate memory board. The whole digital section is the processor board plus the interface board.

**Firmware family:** BGC (DSP code), BGP/BGD (logic tables), BGX (mainframe Z80), BGR (remote Z80). The DSP code sits in **socketed bipolar PROMs**, which any owner with a PROM reader could dump.

**A lead for you:** your manual photos are pages 65–74, covering boards **7 246 125–129**. The board numbers run in order: 121 input, 122 output, **123 processor**, 124 display, 125 interface, and so on. That makes it likely the **processor-board schematic (7 246 123) and the analog boards are on the pages just before 65** of the same service document. Asking whoever posted pages 65–74 for the earlier pages is the most valuable request you can make.
- The processor schematic would not contain the algorithm.
- It *would* show the delay-memory addressing (how many address counters and offset tables, the I/O port map). That sets the hard limit on delay taps per sample.

**Cycle-cost check (TMS32010 User's Guide, TI 1985; the one-cycle MAC instructions are covered in its §3.4.3):**
- IN/OUT take 2 cycles; TBLR/TBLW take 3; LT, LTA, LTD, MPY, MPYK and APAC take 1 each; on-chip data RAM is 144 words.
- **A gain-weighted delay tap costs about 5–7 cycles:** address OUT (unless hardware supplies it), data IN, LT/MPY/APAC.
- That gives roughly 25–35 such taps per DSP, or **about 50–70 per sample for the whole machine**, shared among the cluster, loops, diffusion, output taps and filters.

### 2.4 Patents [confirmed by search; the absence of a patent is an inference]

- The only digital-era EMT/Franz reverb patent family is the EMT 250's:
  - **DE 27 19 276 / US 4,181,820** "Electric reverberation apparatus";
  - Bäder & Blesser; priority 1977-04-29.
- It describes:
  - pre-delay;
  - about 15 loops (short all-pass loops plus long tapped loops);
  - tap gains with exponentially decaying energy, and non-harmonic tap placement;
  - frequency-dependent feedback;
  - random delay modulation;
  - a 3-tap "early density" path.

  It almost certainly describes **program E**, not S/M/L.
- **No patent was found for the 246's "optimierte Eigentonverteilung" algorithm.** The search covered:
  - DEPATISnet, by assignee, inventor and G10K 15/12 classification, 1980–1995;
  - Google Patents.

  Espacenet full text was not reachable, so a filing under an individual engineer's name with an unusual title could have been missed.
- Useful contrast: our IR shows **no measurable modulation and four harmonic loops** (but see the 1986 "stochastic control" claim, §10a). That makes the 246 S program a different design philosophy from the 250 patent.
- Related: DE 30 20 247 / US 4,370,523 (Bäder, 1980) covers A/D conversion, not reverb.
- Local copy of the 250 patent text (scratchpad, not in the repo): retrieved by the research agent; the source URL is https://patents.google.com/patent/US4181820A/en

### 2.5 History, documentation and existing captures

See §10.

---

## 3. The supplied IR

`eugtone_EMT246_CLA_shortRT.wav`: 96 kHz, 24-bit, stereo, 2.0 s. Your notes on its settings: program S, LP max, delay "≈10 %", amplitude "≈90 %".

Full numbers are in `analysis/shortRT/summary.json`, and the plots are in `analysis/shortRT/report.png`.

### 3.1 Sample clock and converter chain

- **Out-of-band check:**
  - The response is flat to 9 kHz, −1.9 dB at 11 kHz and −7.6 dB at 12 kHz.
  - It falls to −57 dB at 13.2 kHz, the device Nyquist.
  - Energy above 13.2 kHz is −75 dB.

  So the capture holds the device's full useful band. Nothing is lost by working at the device rate.
- **Device clock:**
  - Cluster tap spacings measured sub-sample at 96 kHz are 54.055, 143.155, 194.219, 288.320 and 353.392 "26.4 kHz samples". The fractional parts grow in proportion to the lag.
  - Refitting to integers gives fs ≈ 26,370–26,374 Hz.
  - 19.2 MHz ÷ 728 = 26,373.63 Hz fits. **[strong]**
- **Deconvolution:**
  - The first R cluster tap is used as the kernel of the analogue chain.
  - After deconvolving with it, every cluster tap in both channels collapses to a single device sample. The only residue is ±0.12 Gibbs ripple from the 12.3 kHz analysis band-limit. That validates the grid recovery.
- **Channel offset:** L lags R by 18.6 µs (0.49 device samples). The magnitude matches within 0.1 dB.
  - A repair thread names a **PCM53** (16-bit DAC) in the output path. One DAC time-multiplexed between L and R produces exactly this half-sample offset. This is **[confirmed]** by the block diagram and the studiorepair photos (§2.3a).
- **Noise floor:** L sits about 6 dB above R. This may be related to eugtone's L-channel null error.

### 3.2 Reflection cluster

| Channel | Tap offsets (device samples from first tap) | Signs | ms |
|---|---|---|---|
| L | 0, 54, 194, 288 | + + + − | 0, 2.05, 7.36, 10.92 |
| R | 0, 143, 194, 353 | + + − + | 0, 5.42, 7.36, 13.38 |

- All taps have equal magnitude, within 1 %.
- The first tap sits 263 samples (9.97 ms) after the file start.
- **CLA setting, read from two images of CLA's unit (October 2026).** The second image shows the numeric displays. Several 7-segment digits in the REFLECTIONS display have dead right-hand segments (b and c), so "10" shows as "'C" and "80" as "EC".

  | Parameter | Display | LED bar | Value | Confidence |
  |---|---|---|---|---|
  | Program | SHORT RT | — | S | **[confirmed]** |
  | LOWPASS | 12 | top | **12 kHz** | **[confirmed]** |
  | DELAY | 'C (faulty) | near the bottom | **10 ms** | **[strong]** |
  | AMPL | EC (faulty) | near the top | **80 %** (60 is possible from the glyph, but the bar favours 80) | **[strong]** |
  | LF× | 1.0 | about two-thirds up | **1.0** | **[confirmed]** |
  | RT | 1.0 | about two-thirds up | **1.0 s** | **[confirmed]**; measured 0.98–1.00 s |
  | HF× | 0.9 | near the top | **0.9** | **[confirmed]** |

  The first screenshot's "HF bar at the very top" was HF 0.9 (one step below max), not 1.0. The bar can't separate adjacent steps.
- **These are the first measurements of the control laws, and they match the data sheet:**
  - **HF ×0.9 → measured RT(7 kHz)/RT(1 kHz) = 0.90–0.91** (per-loop estimates in model_s_results §3.2: 0.883–0.891 against 0.981–0.991). So **the HF factor is applied at 7 kHz as specified.** Whether any fixed HF damping remains at HF ×1.0 is still unknown.
  - **LF ×1.0 → measured RT(100 Hz)/RT(1 kHz) ≈ 1.00.** It is flat at LF. This doesn't settle the 50 Hz against 250 Hz reference question.
  - **DELAY 10 ms → the first cluster tap is at 9.97–10.08 ms from the file start.** So eugtone's file **keeps absolute timing**: the latency was compensated, with only about 0.1 ms of converter delay left. The pre-delay puts tap 0 of the cluster **at the DELAY time itself** **[strong]**. In samples, 10 ms = 263.7, so the unit uses 263 or 264 samples plus about 2 samples of converter latency.
  - **AMPL 80 %:** cluster taps sit about 29 dB above the early reverb level. The amplitude law needs the library's AMPL sweep.
- Two open questions for further IRs:
  - Whether the pattern is the same in M and L ("feste Zuordnung", fixed assignment, suggests yes).
  - How the 11 amplitude steps map to gain.

### 3.3 Late reverb

- **Decay:**
  - T30 is 0.98–1.00 s from 63 Hz to 2 kHz, 0.95 s at 4 kHz and 0.88 s at 8 kHz, matching closely in L and R.
  - The likely setting is RT ≈ 1.0 s, LF ×1.0, and HF around ×0.9 or ×1.0. The last two cannot be separated from the converter roll-off with one IR.
  - The EDC is very straight over more than 60 dB.
- **Onset:**
  - The reverb's first arrivals come about 580–650 samples (22–25 ms) after the first cluster tap.
  - Echo density reaches Gaussian (NED ≈ 0.9) about 35 ms after the cluster. That is very fast, which implies **dense diffusion before or after the loops**.
- **Loops:** as described in §1, item 3.
  - **Effective loop lengths fitted from the modes are P + 0.2 samples.** *Resolved later:* this is a **~100 ppm clock offset** between the capture interface and the unit. It scales with Pᵢ, and the cluster taps drift by the same ratio. The loop filters add ≤ 0.02 samples (model_s_results §3.1). The loops are pure integer delays.
  - Because the modes are exactly harmonic, **there are no all-passes inside the loops**. The diffusion has to sit outside them.
  - **The four loops decay at matched rates.** The flattened autocorrelation stays about 0.2 at every multiple. That is consistent with loop gains gᵢ = 10^(−3·Pᵢ / (fs·RT)).
- **Excitation:**
  - The waveform injected into each loop is noise-like and spread evenly across the whole period. Something long and dense feeds each loop, or each loop has many input taps.
  - Each loop gets a different pattern.
- **Stereo:**
  - L and R are the same four loops read with a different time offset per loop. L−R correlation of each loop's component is 0.6–0.7 at lags −568, −131, −27 and −330.
  - The broadband L/R correlation is 0.004.
- **Mode distribution:**
  - The modal density of 5,758 samples ÷ fs is 0.22 modes/Hz. At RT 1 s that gives a modal overlap of about 0.5. That is below Schroeder's "≥3" criterion, so colouration control must come from something else.
  - The four lengths are *not* obviously optimal on simple spacing metrics. Three of them share a factor of 3, so some harmonics coincide, for example near 8.79 kHz.
  - *Update:* test T2 found the per-loop mode amplitudes are noise-like (about 5 dB spread), not equalised. So "optimierte Eigentonverteilung" is **not** about equal mode amplitudes. It remains unexplained; the 1986 review links it to "stochastic control" (§10a).

### 3.4 Working topology (program S) [hypothesis, testable; updated with prototype results]

```
in (mono) ─► A/D (SAR, 16 bit) ─► pre-delay = DELAY setting (integer samples)
                                    │
                                    ├─► cluster: tap 0 at the DELAY time; 4 taps → L (0,54,194,288 / + + + −),
                                    │            4 taps → R (0,143,194,353 / + + − +), × AMPL
                                    │
                                    └─► injection / diffusion (dense first pass; possible 460-sample stage, g ≈ 0.44)
                                           │      │      │      │
                                        loop1  loop2  loop3  loop4      each: z^-P · gᵢ · Hᵢ(z), integer P
                                        1321   1389   1461   1587       gᵢ = 10^(−3Pᵢ/(fs·RT)); Hᵢ = first-order, HF factor at 7 kHz
                                           │      │      │      │
                          read-out per channel: shared lags 26 / 332 samples, per-loop weights
                                           └──────┴──┬───┴──────┘
                                                 L, R (+ cluster) ─► LP (2–12 kHz) ─► one D/A ─► 2 × S&H (L ½-sample late)
```

What is still unknown:
- ~~Where the LP sits~~: **downstream** of the reverb and cluster. Studio Magazin 96 says "nachgeschaltet"; the data sheet diagram and the video narration agree.
- **How the injection is built.** Free FIRs work; the true sparse structure is not yet found.
- **Whether the 460 stage is real.**
- **What the "stochastic control" is.**

The live version of these facts is kept in **`docs/EMT246_unit_facts.md`**.

---

## 4. Analysis toolkit (in this repo)

```
tools/deconv_grid.py   # 96 kHz capture → device-rate IR (deconvolve analogue chain, resample on 19.2 MHz/728 grid)
tools/analyze_ir.py    # full pipeline → analysis/<name>/{summary.json, report.png, h_dev.npy}
```

The pipeline stages, each a reusable function:

1. **Basics:** onset, peak, noise floor, out-of-band energy.
2. **Device-grid recovery:** analogue chain deconvolved, fractional L/R offset removed.
3. **Cluster taps** as integers.
4. **Decay:** Schroeder EDC with noise compensation; T20, T30 and EDT per octave.
5. **Echo density:** Abel & Huang NED.
6. **Loop identification:**
   - autocorrelation of the decay-flattened tail, keeping fundamentals only;
   - harmonic-series validation on the high-resolution spectrum.
7. **Periodic decomposition:** a least-squares split into one periodic component per loop, giving the variance explained and the L/R lag per loop.
8. **Time-invariance test:** spectral fine-structure correlation between decay windows. Reference points:
   - 0.67 measured here;
   - 0.60 for a synthetic static modal reverb;
   - 0.02 for exponentially decaying noise.

Run it on any new IR:

```bash
.venv/bin/python tools/analyze_ir.py irs/<file>.wav --out analysis/<name>
```

### 4.1 Analyses to add as more IRs arrive

| Technique | What it nails down |
|---|---|
| **Per-harmonic decay fitting.** Each loop's harmonics are known exactly, so track the decay of every harmonic (short-time DFT at the known frequencies, or a matrix-pencil/ESPRIT fit per harmonic) | Each loop filter's magnitude response \|Dᵢ(f)\| at hundreds of frequencies. Fit 1st/2nd-order shelves and map them to LF/HF steps |
| **Group delay of the harmonics** (the effective-length offset against frequency) | Damping filter phase, which tells a one-pole apart from a shelf apart from an FIR |
| **Differencing IRs across one parameter** (e.g. amp 0 % against 100 %, RT step n against n+1) | Isolates the cluster from the reverb exactly. Shows whether the cluster feeds the loops. Gives gain-law slopes |
| **Sparse diffusion identification** on the per-loop excitation: greedy search over integer all-pass delays and Q15 gains, or a cascaded-all-pass fit | Diffuser topology, if it is an all-pass chain. If it is not, fall back to measured FIR injection filters (§5, path B) |
| **First-pass analysis:** the first P samples of the reverb, before any recirculation | Injection taps and diffusion structure without loop wrap-around |
| **Cross-program comparison** (S/M/L) | Whether M/L use different loop lengths or counts ("adapted density"), and memory usage against the 64 k-word budget |
| **Level series** (the same setting at −30/−18/−6 dBFS) | Nonlinearity, 16-bit storage truncation, and the noise floor |
| **Repeat captures** of the same setting | Confirms time-invariance and gives the noise floor of the measurement chain |
| **Settings inference.** Measured RT, cluster gain, pre-delay and LP corner snap to the quantised steps | Lets unlabelled or mislabelled library IRs be used |

---

## 5. Path to a plugin

**Path A: full structural model (target).**

1. Run the core at fs_dev = 26,373.6 Hz. Resample in and out, and emulate the 12 kHz analogue band-limit and the ½-sample L offset.
2. Build the network from the identified topology: integer delays, 4 loops for each program, damping filters, cluster taps, diffusion and LP. Use the measured constants.
3. Fit the parameter laws. Map RT, LF, HF, amp, LP and delay step indices to coefficients, as lookup tables or closed forms (e.g. gᵢ from RT).
4. Optionally store state in 16-bit with the original truncation, for noise-floor authenticity.

The CPU cost is tiny, so it would run comfortably as a plugin.

**Path B: grey-box hybrid (fallback, also a stepping stone).** Keep the identified loops, but feed each one with a *measured* injection FIR of length about Pᵢ (estimated by the periodic decomposition) instead of a reverse-engineered diffuser. RT, LF and HF stay fully parametric through the loop gains and filters. This works even if the diffuser never yields to analysis. It costs four partitioned convolutions of about 1.5 k taps per channel, which is trivial.

**Path C: interpolated convolution (baseline only).** This is not recommended as the product. The quantised parameter space (S alone is 16 RT × 16 LF × 10 HF × 11 amp × 8 LP) is far too large to capture, and convolution gives up the loops' parametric behaviour. It is still useful as a reference renderer for A/B checks.

**Validation:**
- A sample-domain null against device-rate IRs (residual dB against time);
- EDC per band within ±0.5 dB;
- NED curves;
- the fine-structure correlation;
- the modal-frequency match;
- blind ABX listening on program material.

---

## 6. IR sampling plan

Because the parameters are discrete and the structure looks separable, a few hundred IRs cover everything. About 40 already pin down the core.

| Set | Program(s) | Sweep (everything else held) | # IRs | Purpose |
|---|---|---|---|---|
| Core RT | S, M, L | RT all 16 steps; LF=1, HF=1, amp 0 %, LP 12 k, delay 0 | 48 | Loop lengths per program; gain law; check whether lengths change with RT |
| LF law | S (spot-check M, L) | LF all 16 steps at mid RT | 16 (+8) | LF shelf mapping |
| HF law | S (spot-check M, L) | HF all 10 steps at mid RT | 10 (+6) | HF shelf mapping |
| Cluster | S, M, L | amp 0, 50, 100 % (+ all 11 on S) | 9 (+8) | Pattern for each program, gain law, whether the cluster feeds the reverb |
| Lowpass | S | all 8 steps, amp 100 %, plus amp 0 % at 2 k | 9 | LP response and where it sits |
| Pre-delay | S | 0, 1, 10, 100, 999 ms | 5 | Confirms the delay is pure and applies to both paths |

*New captures are not possible (no hardware access, ever). The sets above are a checklist to match against existing libraries (harness C5 coverage map). Repeat and level tests are replaced by **joint multi-IR fitting** across library IRs, plus any separate-capture duplicates the harness finds (A7).*

**Capture spec** (kept for reference only; no new captures are possible):
- 96 kHz/24-bit ESS sweeps of at least 10 s;
- ≥ 3 s tail, or 25 s for L at long RT;
- one mono input, two outputs;
- a **loopback channel recorded at the same time**, for absolute latency and clock-drift measurement;
- 100 % wet;
- all settings **photographed from the display** (step values, not slider "%").

**Sources, given that you have no hardware:**
1. **Tomkin Mastering's EMT 246 library: buy this first.**
   - More than 2,400 settings across Short/Med/Long RT and Non-Lin.
   - 96 kHz/24-bit, captured through a Prism Orpheus; about €69.
   - It has no low-pass settings, which is a gap for the LP set. It may well cover most of the Core RT, LF, HF and Cluster sets above.
   - Check its naming and metadata scheme and its EULA (analysis use is fine; redistribution is not).

   Other sources to check:
   - **Impulse Record / Wave Arts "Vintage German DSP":** includes the 246; IR count unknown.
   - **Past to Future "CLA" combos:** 36 IRs across four units.
   - **Altiverb:** the 246 is probably in the library (unverified).
   - **eugtone's Gearspace thread:** includes a null-test file.

   Run settings inference on all of them.
2. ~~A paid capture session~~: **not available** (no hardware access). The LP steps will have to come from other libraries, if any have them, or stay approximated.
3. ~~Mainframe ROM dump~~: **not available.** Ruled out (October 2026). The project is IR-only, and the analysis harness (`docs/analysis_harness_spec.md`) is designed for that.

---

## 7. Risks and unknowns

| Risk | Likelihood | Mitigation |
|---|---|---|
| Not enough labelled IRs | Medium (High before Tomkin was found) | Buy the Tomkin library; settings inference; a paid capture day for the LP and robustness sets |
| **"Stochastic control" (Studio Magazin 1986, §10a)** plus slight time variation on L (eugtone's null test) | **Medium–High** | The 1986 test describes S/M/L as using stochastic control to raise eigentone density "without audible FM". Our single IR looks near-static: harmonic loops, stable fine structure, though L's noise floor is 6 dB higher. Possibilities: the randomisation is too slight to show in one IR; it acts on parts that average out in a sweep; it affects mainly L; or the production firmware differs from the tested prototype. **No repeat captures are possible.** Use joint multi-IR fitting and any separate-capture duplicates in the libraries (harness A7 → E1). If confirmed, characterise rate, depth and location from sweep-to-sweep differences and model it |
| The digital-I/O variant (44.1/48 kHz sync) may run the algorithm at a different rate | Low | Label each capture with the unit variant; the loop lengths in ms will show it |
| M/L programs use a different topology, not just different lengths | Medium | Same toolkit; the core RT sets will show it immediately |
| Injection diffusion can't be identified structurally | Medium | Path B measured FIR injection (no audible loss) |
| Converter and analogue chain colour the "sound" | Medium | Already measured; emulate the kernel. Check the analogue input stage with level series |
| Unit-to-unit variation, or captures from a faulty unit | Low–medium | Compare IRs from different sources; check the half-sample offset and the clock |
| A clock offset between the capture and the unit biases loop lengths at large lags | Low | Loopback channel; integer constraints self-calibrate |
| Single-IR conclusions overfit | Medium | Every claim here is labelled; re-run on each new IR |
| Trademark: "EMT" and "246" naming | — | Use a non-infringing product name; check who holds the EMT marks (§10) |
| IR licence terms (eugtone / library EULAs) | — | Analysis use is usually fine. Do not ship third-party IRs in the product; check each EULA |

---

## 8. Recommended next steps

> Superseded by `NEXT_STEPS.md`. Tomkin will supply **separately captured 96 kHz and 48 kHz** files, including complete Non-Lin. The list below is kept for history.

1. **Buy the Tomkin library** (about €69, 2,400+ settings), plus any other 246 IRs available. Drop them in `irs/` and run `analyze_ir.py` on all of them.
   - The first questions to answer are: are the loop lengths the same across RT steps, and what are the loop sets for M and L?
2. **Write a batch comparator.** It should line up loop sets, cluster patterns and decay laws across IRs and flag outliers.
3. **Extend the toolkit with per-harmonic decay and group-delay fitting** to identify the loop filters.
4. ~~Build a Path-B prototype for program S~~: **done**, and being tightened (`tools/model_s.py`, `docs/model_s_results.md`).
5. **Documents:** ask whoever posted manual pages 65–74 for the earlier pages (the processor board, 7 246 123); also email emt-archiv.de.
6. **Decide on program N.**

## 9. Open questions for you

- ~~Is the panel photo the CLA setting?~~ No. A screenshot of CLA's unit (§3.2) shows a different setting: AMPL high, HF at max. Only the LED bars are readable, not the digits.

- ~~Is program N in scope?~~ **Yes, as a bonus** after S/M/L (`NEXT_STEPS.md` §2.4).
- ~~How was the eugtone capture aligned?~~ Unknown. This doesn't block anything: the harness recovers latency from the library itself by regressing first-tap time against the pre-delay setting (harness spec, B3).
- Are there other EMT 246 IRs you own, and what metadata comes with them?

## 10a. Period sources found (October 2026 search)

**Primary documents** (local copies in the session scratchpad; the URLs below work):
- **EMT 246 data sheet** (the same document as your PDF): https://emt-archiv.de/wp-content/uploads/EMT-246-Digital-Reverberator-DE.pdf. Also https://emt-archiv.de/emt-246/ and the EMT Hauptkatalog 1988/89, catalogue pp. 24–27.
  - **emt-archiv.de** is an EMT archive project. It is worth contacting directly about the 246 operating and service manuals.
- **Studio Magazin Nr. 96 (Jul 1986) pp. 14–15**, "Neue Rechnerstruktur EMT 246": https://bibliothek.echoschall.de/studio-magazin-1986-juli-heft-96/page/14-15
- **Studio Magazin Nr. 99 (Nov 1986) pp. 35–38**, "testreport: EMT 246 – Hohe Eigentondichte": https://bibliothek.echoschall.de/studio-magazin-1986-november-heft-99/page/34-35
  - This was a listening test of a **prototype**; it has no bench measurements.
- **EMT Courier Special Issue No. 1 (July 1985)**, covering the 251 and 252; editorial responsibility K.O. Bäder and Dr. Frank Hirsch: https://tile.loc.gov/storage-services/master/mbrs/recording_preservation/manuals/EMT%20Courier%20(July%201985).pdf
- Ads and listings:
  - Studio Sound 1986–88: "optimized eigentone distribution … 0.25 to 20 seconds … 190 program memories".
  - Mix (Oct 1987): US introduction March 1987, US $9,460.
  - db (Jul/Aug 1987).

**New facts about S/M/L** (quoted sources, German):
- **"Stochastic control" [confirmed claim, mechanism unknown].**
  - Studio Magazin 99: S/M/L share "ein neuer einheitlicher Algorithmus … der eine neuartige, erstmals realisierte stochastische Steuerung verwendet, die den Effekt einer Eigentonverdichtung bewirkt … ohne wahrnehmbare Frequenzmodulation". In English: a new unified algorithm using a novel stochastic control that densifies eigentones without perceptible FM.
  - The programs differ in "dem Grad der eingesetzten Stochastik, also der Eigenton- und Rückwurfdichte, aber auch in der Geschwindigkeit des Energieaufbaus … und in der Verteilung der frühen Reflexionen". In English: the degree of stochastic control (eigentone and echo density), the speed of energy build-up, and the distribution of early reflections.
  - **This must be reconciled with our near-static measurement.** See the §7 risk table and the research doc.
  - Note that "the distribution of early reflections" differs per program. So **the cluster pattern may differ between S, M and L.**
- **Cluster:** four reflections per output, "deren zeitlicher Abstand zueinander nicht verändert werden kann. Möglich ist allerdings ein zeitlicher Versatz der beiden Reflexionsgruppen insgesamt". In English: their spacing is fixed, but the delay setting shifts both groups as a whole. That matches the measured patterns.
- **Low-pass:** the timing diagram places it after the cluster. In the video, the narrator says it is "across the whole reverb mix".
- **Stereo:** a correlation meter showed "leichte Gegenphasigkeit (100–110 Grad) im Ausklang" for S and M, i.e. slight anti-phase in the tail. That is consistent with our near-zero L/R correlation.
- **LF reference frequency discrepancy:** the data sheet table and Studio Magazin say **50 Hz**. The data sheet diagram and the remote's parameter card say **250 Hz**. Treat this as unknown until LF sweeps are measured.
- **Non-Lin (N):**
  - "Durch Aneinanderreihen zweier Nachhallverläufe", i.e. two decay curves joined end to end: pre-delay, then SPAN (the early decay), then the late reverb.
  - **The early section can be reversed** (a negative sign). In that case the reverb reaches full amplitude only at the end of the span.
  - SPAN is 0.1–0.5 s per the data sheet; the prototype allowed 0–2 s. A short SPAN cuts off abruptly, which gives gated or inverse effects.
- **Further details from a second read of Studio Magazin 96 and 99:**
  - The low-pass is "nachgeschaltet" (downstream of the reverb).
  - S/M/L differ audibly "mainly in the temporal spacing of the early reflections".
  - L has "comparatively low echo density".
  - The SPAN in N truncates the early decay abruptly.
  - The 246 is the direct successor of the 245; the 245's service manual describes 19 loops in mixed series/parallel.
  - The prototype's display read `EMT 246 | F-86 | 10 | 950 0 | 2 12 05`, i.e. no decimals on the prototype.
- **Remote behaviour:**
  - lever 1 = low-pass;
  - lever 2 = delay and amplitude;
  - lever 3 = RT, LF and HF (or early/late RT in N);
  - holding a lever for more than 1 s speeds up the change.
- **Variants:** a digital-I/O version has **4 × TMS320**, fs locked to 44.1/48 kHz, 256 KB RAM and 32 KB ROM. IRs from such a unit would have a different sample grid.
- **Lineage:**
  - The 251 and 252 used a bit-slice processor at about 32.6–32.75 kHz with different reverb programs.
  - Studio Magazin says the 246 software "mit seinem Vorgänger in keinerlei Beziehung mehr zu vergleichen", i.e. it no longer bears any comparison with its predecessor.
  - Only the 250 program (E) is shared. The 252 also had a Non-Lin with two concatenated processes, and a 6-reflection cluster; the 246 has 8.
- **Launch:** shown at AES Montreux in March 1986, priced at DM 13,268 + VAT with the 246 S remote.

**Foleytronics video (YouTube PZj4ZSpi9P8, Jul 2025): reviewed in full** (captions plus frame inspection).
- It is a listening comparison of the 252 and 246. **The 246's internals are not shown.**
- Usable display readings:
  - preset **F‑51 "S SNARE1"**: LP **12**, DELAY **10**, AMPL about 50 (blurred);
  - the presenter later set RT to **1.0** and LP to **2.0** in the same preset;
  - preset **F‑40 "SPACE"** (program L; matches the parameter card): LP **8.0**, DELAY **120**, AMPL **20**;
  - the STRINGS preset (L) at RT **3.6 s**; Cathedral at about **5 s**; SPACE at about **15 s**; L goes down to **2 s**.
- **Display formats:** the program display shows the letter and the name ("S SNARE1"). LP shows one decimal below 10 kHz ("2.0", "8.0") and "12" at the top. RT shows one decimal ("1.0", "3.6").
- The narrator's specs ("44.1 kHz, 20 kHz bandwidth" for the 246; "24 kHz" for the 252) contradict the data sheet and our measurements. **Treat them as unreliable.**
- **Presets go beyond the card's F1–F50** (SNARE1 = F51).

## 10. History, documentation and capture sources

*Sources are dealers, forums and Wayback captures unless stated otherwise. Verify before relying on them.*

**History**
- **Introduced around 1986.** Vintage Digital calls it the successor to the EMT 244/245.
  - Timeline: 250 (1976), 251 (1980/81, the brochure's "second generation"), 244 (1980), 245 (1982), 252 (1985), 246 (1986), 248 (1993).
  - Vintage Digital's marketing says it carried "EMT 250 algorithms". Your brochure says S/M/L are a *new* algorithm and the 250 is program E. A 2020 Reverb listing mentions an "EMT 250 Option", so program E may have been an optional extra.
- **Designers:** no source names them. Karl-Otto Bäder (EMT technical director; co-inventor of the 250 with Barry Blesser) is plausible but unconfirmed.
- **Variants:**
  - mainframe with remote (9 246 000);
  - built-in panel (9 246 001);
  - an optional digital-I/O version (SONY 1610/1630 or AES/EBU, synced to 44.1/48 kHz);
  - the EMT 246 F fader unit.
- **Corporate:**
  - Barco bought EMT in 1989.
  - In 2003 it sold "EMT including trademarks" to Walter Derrer, and the company continued as EMT Studiotechnik GmbH in Mahlberg.
  - The trademarks are now held by **EMT International GmbH, Switzerland** (EMT Tontechnik imprint).
- **Users and character:**
  - Chris Lord-Alge uses one as a drum reverb, alongside a Sony DRE-2000, a Yamaha REV-1 and a Bricasti M7. Your IR is his "Short RT" drum setting.
  - Owners describe it as "dark (not dull but natural) and rich".
  - A 2005 owner of both the 252 and the 246 preferred the 246.
- **Market:** used units sold for $1.2 k–$7.5 k (2016–2022). A 2020 listing with remote and service manual asked about $14 k.

**Documentation**
- **No public operating or service manual exists** on any archive the agent could reach: audiocircuit.dk, archive.org, Vintage Digital, elektrotanya, hifiengine, synthmanuals or radiomuseum.
- **GroupDIY "EMT 246 schematic" thread (2017):** https://groupdiy.com/threads/emt-246-schematic.65233/
  - It linked a full service manual on mega.nz, which is now dead.
  - It also points to schematics posted by "Mr. Harrison" on Gearspace: https://gearspace.com/threads/whats-a-good-price-for-an-emt-246-with-remote-these-days.982513/
  - **Resolved:** the manual photos used in this study are the ones Mr. Harrison posted there. They cover the remote/panel boards (pp. 65–74) only. No mainframe DSP-board schematics are known to exist publicly.
- **Other leads:**
  - The 2020 Reverb seller "MadSparrow Shop" (Guillaume) had a copy of the service manual.
  - The EMT 250 service manual is public: https://audiocircuit.dk/downloads/emt/EMT-250-rev-sm.pdf. It is useful for the analogue and converter design lineage.
- **No ROM or firmware dumps** are online.

**Existing captures and emulations**
- **No plugin emulates the 246.** The nearest is PSP 2445 EMT, which covers the 244/245.
- **IR libraries:**
  - **Tomkin Mastering:** more than 2,400 settings (Short/Med/Long/Non-Lin, no LP settings), 96/48/44.1 kHz. https://tomkinmastering.com/content.php?c=emt-246-ir&l=e
  - **Wave Arts / Impulse Record "Vintage German DSP":** https://wavearts.com/products/libraries/vintage-german-dsp
  - **Past to Future "CLA Ultimate Drum Reverb Combos":** https://pasttofuturereverbs.gumroad.com/l/fjrkr
  - **eugtone's free CLA IR and null-test thread:** https://web.archive.org/web/20230630181220/https://gearspace.com/board/high-end/1352418-study-cla-drum-verbs-emt-246-dre-2000-rev-1-free-irs.html
  - A HISE forum user also captured a 246 and noted it "really only goes to 12k", which matches the measured band-limit: https://forum.hise.audio/topic/2712/capturing-impulses-software-techniques/30
- **People who might help** with a capture session:
  - Tomkin Mastering;
  - Foleytronics FX Repair (YouTube 252/246 teardown, July 2025);
  - Gearspace users t-hiho (owner of the CLA-setting unit) and eugtone;
  - MadSparrow Shop;
  - lalalandsynth (HISE forum);
  - Nikolaus Riehm (studiorepair.de).

**Legal notes (not legal advice)**
- "EMT" is a registered trademark of EMT International GmbH. A product name or marketing that uses "EMT"/"EMT 246" needs a licence; Universal Audio's EMT 140/250 plugins are licensed. Otherwise use a descriptive, non-infringing name.
- The 1977 EMT 250 patent has long expired.
- Reimplementing an algorithm identified from measurements is generally a different matter from copying firmware code. Take advice before using anything derived from a ROM dump.
