# EMT 246: consolidated unit facts

The single place for *what we know about how the unit works*. Details and evidence live in the linked docs; this sheet is updated whenever a finding changes.
- `discovery` = EMT246_discovery.md
- `research` = algorithm_research.md
- `model` = model_s_results.md

**Confidence:**
- **C** = confirmed by a primary source or a direct measurement;
- **S** = strong: several independent measurements agree;
- **H** = hypothesis.

Constraints: there is no hardware access, ROM dump or new capture, ever. Evidence comes from existing IRs, documents and photos only.

*Last consolidated: October 2026.*

---

**At a glance:** **mono in, stereo out** (C). A single balanced input; L and R outputs come from one multiplexed DAC, with L half a sample late. Each stereo IR is the complete 1-in/2-out response of a setting.

## 1. Hardware

| Fact | Conf. | Source |
|---|---|---|
| **2 × TMS32010NL** DSPs (TI, date codes 8540/8618); processor board 7 246 123 | C | studiorepair photos (discovery §2.3a) |
| **38.400 MHz oscillator ÷ 2 → 19.2 MHz CLKIN → 208.3 ns per instruction → 182 instructions per sample per DSP** | C (÷2 is S) | photos; TI TMS32010 User's Guide (max CLKIN 20.5 MHz, cycle = CLKIN/4) |
| **fs = 38.4 MHz / 1456 = 26,373.63 Hz**, not the brochure's 26.4 kHz | S | integer cluster-tap grid; loop harmonics; the 728 = 4 × 182 arithmetic |
| TMS32010 timing: IN/OUT 2 cycles; TBLR/TBLW 3; LT/LTA/LTD/MPY/APAC 1; 144-word on-chip RAM. A gain-weighted delay tap costs about 5–7 cycles, i.e. **about 50–70 taps per sample for the whole machine** | C (timings) / H (budget) | TI User's Guide |
| DSP program: **6 × AM27S191 bipolar PROM (2 K × 8), "BGC 1A–6A" = 12 KB** (4 K words on DSP 1, 2 K on DSP 2) | C | photos; data sheet ("ROM 12 kByte") |
| **2 empty 24-pin PROM sockets** beside BGC 5/6 | C (sockets) / H (purpose: the EMT 250 option) | photos; a 2020 listing mentions an "EMT 250 Option" |
| Delay RAM: **4 × NEC μPD41464C-12 (64 K × 4 DRAM) = 64 K × 16-bit words** (128 KB = 2.48 s at fs) | C | photos; data sheet |
| A spare MB81464 DRAM sits in a "Spare" position; 7 unpopulated 18-pin footprints | C | photos |
| 2 × UM6116 2 K × 8 SRAM on the DSP board | C | photos. Coefficient/parameter RAM is **H** |
| Address-generation logic (74LS193 counters, 74S188 PROMs, "BGP 1A/2A", "BGD 1A"), so delay addressing is partly in hardware. **Now confirmed: base counter + tap offset memory (§2b)** | C | photos; service manual p. 53 |
| Mainframe control: **Z80A + CTC + PIO, EPROM "BGX 1B"** (interface board 7 246 125). The step → coefficient tables are probably here | C / H | photos; data sheet |
| Remote 246 S: Z80A, 27256 EPROM "BGR 1A", 6116 RAM, ICM7218B display driver, MK48Z02 memory module, quadrature opto controllers | C | manual pp. 65–74 |
| Remote link: 5-pin DIN, 32-bit pulse packets (one pulse = 0, three = 1). **The unit won't run without the remote** | C / reported | manual p. 65; Gearspace |
| Arithmetic: 32-bit accumulate, **16-bit storage of intermediate results**, no floating point | C | data sheet |

## 2. Converters and analogue path

| Fact | Conf. | Source |
|---|---|---|
| **Mono input** → amp → anti-alias LP (TL07x daughterboard) → S&H → **discrete 16-bit SAR ADC**: PCM53 DAC + 2 × 74LS503 SAR + LT1016 comparator | C | block diagram; photos (board 7 246 121) |
| Output: **one PCM53 DAC → LF13333 switches → 2 × S&H → 2 identical active LP filters** (TL074 daughterboards) → amps | C | block diagram; photos (board 7 246 122) |
| **L output is 18.6 µs (0.49 sample) later than R**, with identical magnitude response (±0.1 dB). Caused by the multiplexed DAC | C | IR measurement + hardware |
| Overall band-limit: flat to 9 kHz, −1.9 dB at 11 k, −7.6 dB at 12 k, −57 dB at 13.2 k; energy above 13.2 k is −75 dB | S | IR |
| Spec: 30 Hz–12 kHz (+1/−3 dB, DELAY program); reverb S/N 73 dB unweighted; THD ≤ 0.03 % | C | data sheet |
| **Noise:** L's floor is about 6–9 dB above R's, broadband and slowly decaying (−88 → −91 dB over 1.1–1.8 s); R's floor is flat. Other reports also show L/R mismatches | S | IR; model §6.5; user |
| Trimmers: input LEVEL, HF CMR, LF CMR, OFFSET; output LEVEL 1/2 | C | photo of the front flap |

## 2b. Service manual findings (EMT 246S service manual, read October 2026)

Source: `sources/emt-documents/EMT 246S service manual 1/2.pdf` (OCR of the English column in `sources/emt-documents/ocr/`, gitignored). Page numbers are the printed ones. Only the text pages are read so far; the schematic sheets are not.

| Fact | Conf. | Source |
|---|---|---|
| **Delay addressing is in hardware: `address = Base Address Counter + Tap Offset Memory[tap number]`**. A 16-bit base counter (Z333–Z336) counts up every sample; a 16-bit adder (Z337–Z340) adds the offset from the TOM (Z331/Z332, buffered by Z329/Z330), which has its own **Tap Number Counter**. So every delay element is an **integer-sample tap of one 64 K-word circular memory**, and the tap offsets are loaded per setting. The two TOM chips are probably the two UM6116 2 K × 8 SRAMs seen in the photos (so up to 2 K tap entries, far above the time budget) | C (structure) / H (TOM = UM6116) | p. 53; SP self-test error list (Test 12) |
| Two TMS320 processors: 16-bit data bus, 12-bit program address bus, **32-bit internal registers**, run synchronously under "main processor A". Data memory: 4 × 4-bit chips (Z347–Z350) = 16-bit words, shared through buffers and transceivers; arbitration in two PALs (Z325, Z326). ROM: processor A has four ROM chips (Z304–Z307), processor B two (Z314/Z315): matches 4 K / 2 K words | C | p. 53; Test 12 |
| **Instruction time about 200 ns; sample rate 26.4 kHz; RAM 128 kB; ROM 12 kB (the DSP PROMs); control CPU Z80** | C | technical data |
| **Reverb S/N is 15 dB worse than the DELAY (pass-through) program:** relative to nominal level, unweighted rms 88 dB (DELAY) against 73 dB (REVERB); CCIR 468-2 peak 82 against 69 dB; weighted 76 against 64 dB. The extra noise in REVERB is the algorithm's own (16-bit storage, truncation or rounding noise through the loops). It is a **quantitative target for the fixed-point noise model (Q10)**. The setting used for the spec is unknown | C (spec) / H (cause) | technical data |
| **Levels:** 0 dB = 0.775 V. Factory nominal level +6 dBu in and out (adjustable −10…+15 dB). The 0 dB LED marks nominal, which is **6 dB below the A/D overload point (REG LED)**: nominal = about −6 dBFS, headroom 6 dB (stated for f < 2 kHz only, so the reserve is smaller at HF, which fits the input pre-emphasis below) | C | pp. 31, 38 (calibration); technical data |
| **Input chain:** balanced active input → **10th-order elliptic low-pass, fc 12.2 kHz, stop-band better than 60 dB, phase-compensated, with HF pre-emphasis** (daughterboard 7 246 131, Z103) → integrate/hold S&H → 16-step successive-approximation A/D (PCM-type DAC Z107, comparator, 2 × 8-bit registers Z108/Z109), **16 comparisons per sample period**, MSB inverted for the processor (two's complement). P104 sets the A/D offset by searching for the **noise minimum** inside a small zero window | C | pp. 37–39; Test 27 |
| **Output chain:** two 16-bit buffers per channel → **one current-output D/A converter, two conversions per sample period** (R then L) → two S&Hs triggered at staggered times → output low-pass (Z210/Z210', daughterboard 7 246 130) → balanced amp. This confirms the half-sample L/R offset by design. A **"Charge Current Adjust" (P202)** trims the HF response; P201 sets the converter gain (+1.5 dB at TP3 with the 1 kHz internal tone) | C | pp. 46, 31–32 |
| **THROUGH test (Test 25, A/D → D/A, no reverb):** the HF trim is set so that **12 kHz in gives −2.5 dB** at the output, and 10 kHz −0.5 to −1 dB (re 1 kHz). **This agrees with the reference kernel** (the cluster-tap response of the CLA IR: −0.1 dB at 10 kHz, −2.0 dB at 12.0 kHz, −7.6 dB at 12.2 kHz; the old 'band-limit' row in §2 quotes the same skirt by rounded frequencies). So the through path and the reverb path have the same band-limit; **the earlier suggestion that the program's own LP adds about 5 dB at 12 kHz is withdrawn** | C (spec) / S (match) | p. 32; `analysis/analog/summary.json` |
| Spec: DELAY response 30 Hz–12 kHz (+1/−3 dB); THD ≤ 0.03 % (nominal, 1 kHz); output CMR ≥ 41 dB (30 Hz–11.5 kHz) | C | technical data; p. 32 |
| **Built-in test signals:** Test 21 white noise and 22 pink noise (−6 dB, **uncorrelated between outputs**), 23 sine 1 kHz, 24 sine 10 kHz, 25 THROUGH, 26 digital zero. The noise generator exists in the DSP code, but nothing in the text links it to the "stochastic control" | C / H | service manual, test-program pages (about pp. 26–31) |
| Spec step counts repeated: RT 16 steps (S 0.2–2.5, M 0.5–5, L 2–20 s at 1 kHz); **LF× 0.1–2 at 50 Hz, 16 steps; HF× 0.1–1 at 7 kHz, 10 steps**; delay 0–999 ms; **reflection cluster = "a group of eight individual reflections with a fixed relationship to one another", amplitude 0–100 % in 11 steps**; LP 2–12 kHz in 8 steps | C | technical data |
| **Processor board schematics (7 246 123, sheets 1–5, printed pp. 54–58):** the board is the whole DSP engine. Two TMS320 (A, B; CLKIN about 20 MHz, CLKOUT about 5 MHz) with instruction ROMs of **4 K × 16 (A)** and **2 K or 4 K × 16 (B)**; a 19.2 MHz crystal through a clock-sync stage; **data memory 64 K × 16 DRAM** (4 × NEC D41464C-12) addressed through an **address adder = base-address counter (16 bit, 4 × 74LS193) + tap-offset memory (TOM)**; the DRAM sits behind a write register (**MIR**) and a read register (**MOR**), with PAL1 (BGP 1A) as the DRAM controller (RAS/CAS/MUX/W) and PAL2 (BGP 2A) as the memory bus arbiter between A and B | C | pp. 54–58 |
| **The TOM is 256 × 16 bit** (two HM6116LP3 with only A0–A7 used), addressed by an **8-bit tap-number counter** (2 × 74LS193) that steps on every memory access. So the machine has **at most 256 programmed tap offsets per sample**, and every memory access (read or write) goes through `base + TOM[n]`. The TOM is written from **processor A's data bus** (strobe LTOM) and the tap counter is loaded by A (LTNC), so A can reprogram offsets while running; whether it does so (tap jitter) is unknown | C (structure) / H (run-time use) | p. 58; p. 57 |
| **Processor A's port map (74S138 decode):** OUT 0 = write memory (MIR), 1 = load tap counter, 2 = load TOM, 3 = OUT FP (to the control CPU / panel via the external bus J301), 4 = to processor B, 5 = DAC 1 data, 6 = DAC 2 data; IN 0 = read memory (MOR), **1 = increment the base-address counter (IN BAZ)**, 2 = from B, 4/5 = front-panel data (FP1/FP2), 6 = A/D data. **The delay-line pointer is therefore advanced by the DSP program (one IN BAZ per sample), and the parameters arrive from the Z80 over the external bus**, which fits the step-to-coefficient tables living in the Z80 EPROM. Processor B has no external ports fitted (its decoders Z360/Z361 and J302 are "not mounted"): it talks to A and to memory only | C | p. 57 |
| **Sample timing is a hardware sequencer:** an 8-bit timing counter (2 × 74LS193) clocked by CLKOUT A (about 208 ns) addresses a timing PROM (MMI 63S281, "BGD 1A") and a register that drive **DA enable 1/2, DA S&H 1/2, A/D S&H A/B, A/D clock (SAR) and BIO** (the DSP's sync input). The jumpers JP306 set the counter preload, so the period is **182 CLKOUT cycles = 728 CLKIN cycles = 19.2 MHz / 728 = 26,373.6 Hz**, exactly the device rate found from the IRs | C | p. 57 |
| **L/R offset is a timing-ROM step:** the two D/A S&H pulses are placed by the PROM in 208.3 ns steps, so half a sample is 91 steps (18.96 µs). Measured offsets: CLA unit 18.6 µs ≈ 89 steps; the Tomkin library gives **19.19–19.21 µs in every one of the 768 S files with a cluster** (92.1 steps; 92 steps = 19.17 µs), i.e. no per-setting variation at all. The step quantisation is C; reading the two units as different firmware/analogue delay (89 vs 92 steps) is **H** (the S&H aperture and filter delays also enter the measurement) | C (quantisation) / H | p. 57; harness B4 |
| **Input filter 7 246 131 (Z103, p. 44) and output filter 7 246 130 (Z210, p. 51), transcribed and simulated** (`tools/analog/`: netlists, MNA solver, `compare_kernel.py`; ideal op-amps or single-pole TL07x). Both are all-op-amp chains of state-variable biquads from 68.1 k + 220 pF integrators (10.6 kHz). **Input:** 14 op-amps; passband **boosts HF by +5.1 dB at 10–11 kHz** (the 'pre-emphasis'), then an elliptic skirt (+0.8 dB at 12.2 kHz, −19 dB at 13 kHz, −52 dB at 14 kHz, all re 1 kHz) and a stop-band below −74 dB. **Output:** 10 op-amps + follower, flat to ±0.5 dB to 11 kHz, −1.7 dB at 12.0 kHz, −3.5 dB at 12.2 kHz, −20 dB at 13 kHz, stop-band about −50 dB. A coupling capacitor (C3109 with R3134, C3009 with R3031) gives high-pass corners near 1 Hz or below. The transcription is validated: the product of the two responses reproduces the kernel's band edge shape | C (parts) / S (simulation) | pp. 44, 51 |
| **Analogue chain = input filter × output filter × a gentle one-pole lag.** Simulated H_in·H_out has a +5 dB hump at 8–11 kHz that the kernel does not have. What is left over (kernel minus the simulated filters) is a smooth roll-off (−1 dB at 3 kHz, −2.7 at 5, −5.2 at 10 kHz) that a **single pole at about 5.3–5.9 kHz** reproduces to 0.2–0.3 dB rms over 0.5–11 kHz (a squared sinc does worse: 0.9 dB rms). That lag is **H = the output S&H / charge-current trim (P202) and the input integrate-and-hold**, not a hold droop. With finite op-amp GBW (3–5 MHz) the skirt falls at 12.0/12.2/12.4/12.6 kHz to within about 0.5/0.5/2/3 dB of the kernel (deeper stop-band values differ by 5–10 dB, tolerance-limited). Nominal component values only | S (fit) / H (lag origin) | `analysis/analog/summary.json`, `chain.png` |
| The unit's sample period is 37.92 µs, so half a period is 18.96 µs. Measured L/R offsets: CLA unit 18.6 µs, Tomkin library median 19.20 µs (MAD 0) | S | IRs; harness B4 |

## 3. Controls and parameters

| Fact | Conf. | Source |
|---|---|---|
| Sections and controls:<br>• PROGRAM (S M L D E N);<br>• MEMORY;<br>• LOWPASS `f (kHz)`;<br>• REFLECTIONS `DELAY (ms)` / `AMPL. (%)`;<br>• DECAY `LF×` / `RT (s)` / `HF×`.<br>One trough controller per section; one LED wedge bar per parameter | C | panel photos |
| **Legend colour = program LED colour:**<br>• red legends = S/M/L/E;<br>• green = D (`L — DELAY — R`, `L — ECHO (%) — R`);<br>• yellow = N (`DELAY`, `SPAN (s)`, `EARLY — RT — LATE`) | S | panel; parameter card |
| Step counts:<br>• RT: 16 steps per program (S 0.2–2.5 s, M 0.5–5, L 2–20);<br>• LF×: 0.1–2, 16 steps;<br>• HF×: 0.1–1, 10 steps;<br>• AMPL: 0–100 %, 11 steps;<br>• LP: 2–12 kHz, 8 steps;<br>• DELAY: 0–999 ms in 1 ms steps | C | data sheet. Exact step values are **unknown** |
| Observed display values:<br>• LP 2.0 / 8.0 / 10 / 12;<br>• RT 1.0 (S), 3.6 and about 5 (L), 15 (L), L minimum 2;<br>• DELAY 10 / 120;<br>• AMPL 20 / about 50 / 80 | C (as read) | video; CLA images |
| Display formats: program shows letter + name ("S SNARE1"); LP one decimal below 10 kHz; RT one decimal | C | video |
| **HF factor is applied at 7 kHz as specified:** HF ×0.9 gives a measured RT(7k)/RT(1k) of 0.90–0.91 | S | CLA IR + confirmed settings |
| LF ×1.0 gives a flat LF decay. **The operation card (high-res scan, `emt246_operation_and_function.jpg`) defines the three decay references explicitly: RT at 1000 Hz = RT, at 50 Hz = RT × LF, at 7000 Hz = RT × HF.** That agrees with the data-sheet table and Studio Magazin (50 Hz) and with the measured 7 kHz HF point. Only the data-sheet diagram says 250 Hz, and the older parameter-card reading of 250 Hz is superseded. Still to be confirmed by LF sweeps in the library | S (card, table and magazine say 50 Hz; one diagram dissents) | operation card; data sheet; Studio Magazin |
| Remote levers work like a **gear shift**: horizontal selects the parameter, vertical steps it.<br>• lever 1 = LP (one position);<br>• lever 2 = 6 parameters: pre-delay and ampl for S/M/L/250, L/R delay for D;<br>• lever 3 = 9 parameters: RT/LF/HF for S/M/L/250, L/R feedback for D, early/late RT for N.<br>Holding a lever more than 1 s accelerates the change. The display shows a coarse LED bar plus the exact value | C | Studio Magazin 99 |
| **The low-pass is downstream** ("Nachgeschaltet ist ein Tiefpassfilter"). The data sheet diagram also places it after the cluster. So model the LP at the output, after the reverb and cluster sum | C (claim) / S | Studio Magazin 96; data sheet |
| **Operation card, optimum RT ranges** (black = optimum, outline = full range): S optimum 0.2–1 s (full 0.2–2.5); M optimum 0.8–4 s (full 0.5–5); L optimum 3–10 s (full 2–20). Full ranges match the data sheet and the library's RT values | C | operation card |
| **Operation card, S/M/L/E model:** DELAY, then a burst of reflections (REFLECTION AMPL. is the burst level against the start of the reverb), then a straight-line (in dB) decay; LOW PASS is drawn as acting on the whole output after the reflections. Entry fields on the card header: Delay (ms), Amplitude (%), "low freq. — Reverb — high freq." (LF× / RT / HF×) | C | operation card |
| **Operation card, N (Non-Lin):**<br>• shown with LOW PASS 2–12 kHz, so LP applies to N as well;<br>• the DELAY precedes the early section; SPAN ("Anhalldauer") is 0…0.5 s;<br>• the early-RT line **pivots about the middle of the SPAN** (a "+/−" arrow around the crossing point at mid-level), so a negative early RT is a rising ramp that reaches the late-RT start level at the end of SPAN, and a positive one is a fall from the same level at the start;<br>• the late RT starts at the end of SPAN at the same level as the early start, and its slope is set by LATE RT (arrow swings the line);<br>• a thin tail from the early line is drawn continuing past SPAN's end (meaning unclear) | C (what is drawn) / H (the pivot reading, and a schematic) | operation card |
| **Front panel in N** (photo `emt246_front.jpg`):<br>• display 1 is an 8-character program/preset name (`NON-LIN.`); display 2 is a 4-digit memory number (blank); display 3 is LP `f (kHz)` (12);<br>• display 4 shows DELAY (0) and SPAN (0.4) in the DELAY and AMPL slots; display 5 shows **EARLY RT in the LF slot (−3.0) and LATE RT in the HF slot (1.2)**, with the RT slot blank;<br>• N legends are amber, D legends green, S/M/L/E legends red;<br>• bars: DELAY bar empty at 0, SPAN about 80 % up, early and late bars mid-low, the middle decay bar dark;<br>• lever symbols: 1, 2 and 3 columns for LOWPASS, REFLECTIONS, DECAY;<br>• level ladder REG +3 0 −3 −6 −9 −12 dB; a memory-card slot at the left | C | photo |
| **The library's Non-Lin file order `LP_DELAY_SPAN_EARLY_LATE` follows the panel's display order** (DELAY, SPAN, EARLY, LATE); the photo's −3.0 early confirms that **negative EARLY RT is a valid setting** and that the library's negative values are those | S | photo; library names |
| **The Tomkin Mastering library** (user-confirmed; README credit Thomas Kienast): 2,437 settings per rate at 48 kHz and 96 kHz (separate captures), S/M/L RT plus Non-Lin; LP always 12, DELAY always 0. **Licence (user): commercial use is allowed; the only restriction is on redistributing the files as they are.** The folder stays gitignored | C (user) | user; README |
| Presets: 90 factory (F1–F90; card shows F1–F50; SNARE1 = F51), F91–F99 service, 99 user | C | data sheet; card; video |

## 4. Programs

| Program | What is known | Conf. |
|---|---|---|
| **S / M / L** (A/B/C), sound | Studio Magazin 99:<br>• audibly, S/M/L differ **mainly in the temporal spacing of the early reflections** ("zeitliche Staffelung der frühen Reflexionen"), so expect a different cluster spacing per program;<br>• S has the highest echo density (small rooms, densification);<br>• M is the universal "concert hall";<br>• L is a "cathedral" with **comparatively low echo density** (interesting at short RT);<br>• all have even decays;<br>• the HF absorption simulation and the LF behaviour with factor > 1 are praised | C (review) |
| **S / M / L** (A/B/C) | "Completely new algorithm, colouration-free through optimised eigentone distribution". They share one unified algorithm with a novel **"stochastic control"** that densifies eigentones "without perceptible FM". They differ in the degree of stochastic control (eigentone/echo density), the speed of energy build-up, and the **early-reflection distribution**. S has the highest density and fastest build-up | C (claims) |
| **E** | "Exactly" the EMT 250 program, with extended pre-delay; likely an option (empty PROM sockets) | C / H |
| **D** | Delay/echo: 2 taps (L, R), 0–999 ms, feedback 16 steps (~0 to −60 dB per repeat) | C |
| **N** (Non-Lin) | SPAN is a "superordinate time value" for the early decay, after which the separate late decay begins. A shorter SPAN **cuts the early decay off abruptly** at its end, giving gated reverb; reversal gives a rising curve, i.e. inverse reverb. Two decay curves joined end to end: pre-delay → SPAN (0.1–0.5 s; prototype 0–2 s) early decay → late reverb. The **early section can be reversed**. Probably 4 parameters (DELAY, SPAN, EARLY RT, LATE RT) with no AMPL/LF/HF | C / H (parameter set) |
| Digital-I/O variant | 4 × TMS320, fs locked to 44.1/48 kHz, 256 KB RAM, 32 KB ROM, so a different sample grid | C (catalogue 1988/89) |

## 5. Program S algorithm (measured from the CLA IR)

The CLA IR settings, from a numeric display with some dead segments: **SHORT RT, LP 12, DELAY 10 ms, AMPL 80 %, LF ×1.0, RT 1.0 s, HF ×0.9** (S).

| Fact | Conf. | Source |
|---|---|---|
| **Time-invariant to measurement precision.** Fine-structure correlation is 0.67, against 0.60 for a static modal reference and 0.02 for noise. Loop modes are exactly harmonic. This must be reconciled with the 1986 "stochastic control" claim | S / open | discovery §1; research §7.2b |
| **Pre-delay:** cluster tap 0 lands at the DELAY time (10 ms → 9.97–10.08 ms in the file; eugtone preserved absolute timing) | S | IR + settings |
| **Cluster:**<br>• L taps at 0, 54, 194, 288 samples, signs + + + −;<br>• R taps at 0, 143, 194, 353, signs + + − +;<br>• all 8 equal in magnitude (±1.6 %);<br>• single integer samples;<br>• spacing fixed; DELAY shifts both groups together | S (pattern) / C (fixed spacing) | IR; Studio Magazin 99 |
| Cluster level at AMPL 80 %: about +29 dB relative to early reverb level. Amplitude law unknown | S | IR |
| **Late reverb = 4 uncoupled recirculating loops of 1321, 1389, 1461 and 1587 samples** (50.1–60.2 ms). They are integer and pure delays: the "+0.2" was a ~100 ppm capture clock offset. No all-pass inside the loops | S | autocorrelation, harmonic test, model §3.1 |
| **Loop gains are length-scaled exactly at the RT setting:** gᵢ = 0.707439, 0.695463, 0.682337, 0.660021, giving per-loop RT₀ 0.9997–1.0017 s at RT 1.0. They are resolved to about ±8 Q15 LSB (1σ), so 3–4 bits short of pinning exact coefficients. Library RT sweeps would close the gap | S | model §9.1 |
| **Loop and 460 parameters are shared by L and R:** parameters fitted on L alone predict R within 0.5 % in annihilator cost | S | model §8.2 |
| **In the tail the algorithm matches to within the noise:**<br>• noise-compensated EDC within 0.2 dB in most octaves to −50 dB;<br>• decay rates within 1 %;<br>• mode frequencies within a few mHz;<br>• R's fine structure at its noise ceiling.<br>The *waveform* null is limited by the excitation (injection), not the loops | S | model §9.2 |
| **Loop filters:** first order is sufficient here. At HF ×0.9 the extra loss is 0.37–0.40 dB per pass at 7 kHz (one-pole p ≈ 0.04). Shared or length-scaled coefficients is undecided | S / open | model §3.2 |
| Modal density 0.218 modes/Hz. The lengths are not optimal on simple spacing metrics; three share a factor of 3 | C (arith.) | discovery §3.3 |
| **Per-loop mode amplitudes are noise-like (spread ≈ 5 dB, Rayleigh-like)**, so not equalised. This argues against a pure all-pass diffuser read by single taps | S | topology T2 |
| **First pass is dense** and carries non-recirculating content (−11/−7 dB of the first 1,400 samples). The reverb starts about 240 samples (L) and 570 (R) after tap 0. So there is multi-point injection or read-out, not a single write and read per loop | S / H | model §4 |
| **5th recursion at 460 samples** (g ≈ 0.45, i.e. RT ≈ 0.15 s). It is the one structural addition that clearly helps: +13/+17 dB in-sample, +6/+10 dB held out. Its excitations show an all-pass-like −g + z⁻⁴⁶⁰ signature. Whether it sits in parallel or in series is unresolved | S (exists) / H (form) | model §3.3, §8.2–8.3 |
| **R: once the input ends (about 0.45 s), the 4 loops + 460 stage predict the rest to within 4–6 dB of the noise floor.** This prediction uses only the 9 loop/460 numbers and the IR's own past, e.g. −45.5 dB at 450–650 ms. So the loops and the 460 stage are right | S | model §10.1 |
| **Correction (§13):** the "late input up to 0.45 s" is **not** a separate generator. A fully LTI synthetic with all taps inside 82 ms reproduces it. It is the **diffuser's own tail** (AP230 at g = −0.7, about 356 dB/s) **plus band-edge ringing from the IR deconvolution**. The row below is kept for history | S | model §13.2 |
| *(superseded)* **Late input:** each loop receives its **own dense injection** that keeps arriving and decays smoothly (R ≈ 290 dB/s, RT ≈ 0.2 s; L ≈ 130 dB/s) until about 0.45 s. There are no spikes. Fine structure is uncorrelated across loops but correlated between L and R of the same loop. **Ruled out:** a shared input, a long EMT 250-style loop (lags 3.5–16 k), and extra 230/460/920 stages. The 460 stage's gain is flat across frequency, g ≈ 0.44–0.46, consistent with an all-pass | S | model §10.2–10.3 |
| **Diffuser found: a common series all-pass chain AP(20) · AP(72) · AP(230), all at g ≈ −0.70**, with the 460 recursion as a plain comb (pole gain 0.445). The delays are sharp (±1 sample costs 4–10 dB). It was chosen on R and confirmed on L at the same gain, and it explains the recurring 20/52/72 and 230/460/690/920 lags. Its cost is about 8 multiplies and ≤ 4.1 k words. For a sparse model it is worth 10–13 dB on the held-out null. Input side vs read-out side is undecidable from one IR; input side is favoured | S (chain) / H (placement) | model §11 |
| **Round 12 (all negative):**<br>• no small tap-delay set is shared across loops/channels (4–64 shared delays explain ≤ 13 dB);<br>• sharing tap positions between L and R of one loop costs about 4 dB;<br>• long-delay-line taps (EMT 250 style) don't help, and the selected taps show no exponential gain decay;<br>• loop-sum feedback is ruled out indirectly, because the modes stay harmonic.<br>The injection is **per loop, dense, and not sparse in any basis tried** | S | model §12 |
| **Sparse-fit basis artefact:** sparse tap fits must include the measurement band-limit. With it, CLA reaches −22 dB at 103 taps but **plateaus at about −23 dB in both channels**, where a matching LTI synthetic reaches −46 to −49 dB. So **the remaining injection/read-out structure differs from our model structurally**. It is not explained by time variance | S | model §13.2–13.4 |
| **R is time-invariant within about 3 % rms per 0.2 s.** Injection randomisation, switched taps or re-draws ("stochastic control" in that form) is **ruled out** | S | model §13.4 |
| **L's excess fits an L-only read-out gain jitter of about 1 % rms per 0.2 s (or 2 % per 1 s).** It is also consistent with L's 26 dB higher pre-onset energy and eugtone's L-only null error. It may be the unit (a mild "stochastic control" on L) or the capture chain; one IR can't tell. The library test is whether it appears in L only, at the same depth, across captures | S (fit) / H (cause) | model §13.4–13.5 |
| **L's excess lives in the recirculating part, not the injection.** L's read-out departs from any fixed read-out all the way down the tail. A time-varying L read-out (the "stochastic control", modulation slower than about 1 Hz or of small depth) is consistent | S / H | model §12.4b |
| The diffuser does **not** generate the late input. The best generative held-out null is still about −29 dB in band (400–600 ms), and the taps in front of the diffuser are still thousands per channel. Next candidate: a small tap set shared across loops and channels (group-sparse) | S | model §11.2–11.3 |
| **Hardware plausibility:** free excitation FIRs (2 × 8.4 k taps) are far beyond the ~50–70 taps per sample budget, so the real injection must be compact. The best guess is **per-loop diffusers** (short all-pass chains or nested all-passes), with lags 20/30/72/230 as candidate delays | H | model §10.4 |
| **Injection is not a sparse tap set.** 1,719 integer taps per channel recover less than the dense model. Something dense (an unidentified diffuser) feeds the loops and **keeps delivering new input at 0.25–0.4 s**. This is now the main limit on the null | S | model §8.3–8.4 |
| **Output:** L and R share lags **26 and 332 samples** in every loop, with **per-loop weights**. A strict common output distributor is rejected | S | model §5; T1 |
| L/R correlation of the tail ≈ 0; the 1986 review's correlation meter showed slight anti-phase (100–110°) | S / C | IR; Studio Magazin 99 |
| Echo density goes Gaussian at about 35 ms (96 kHz NED) after the cluster | S | IR |
| **Best nulls so far (model Cs, raw 96 kHz domain, 20 Hz–11.3 kHz):**<br>• cluster / 0–50 ms: −50/−52 dB;<br>• 50–300 ms: −48 dB (in-sample);<br>• **held-out 400–600 ms: −27/−29 dB** (+21/+25 dB over noise);<br>• 600–1000 ms: −17/−19 dB.<br>The null is not yet at the noise floor anywhere, so this is **not yet proof of the true algorithm** | S | model §8 |
| Cluster taps are equal to 0.1 % (R: 0.9999 / 1.0010 / −0.9992 / 0.9999 in kernel units). L sits about 0.06 dB above R. The freely fitted L offset is 0.492 samples | S | model §8.5 |
| **Analogue chain** = band-limit FIR × a **2nd-order sub-20 Hz high-pass** (coupling capacitors **H**). That high-pass is what had limited the device-domain cluster null | S | model §8.5 |
| **L-only tail excess:** +4 to +15 dB over noise (and 15–22 dB in the input-ended predictor at every cut-off time T), which no time-invariant change removes. It shows periodicities at loop-length *differences* (68, 126, 140, 198, 266 samples). This fits a **slow modulation of loop/read-out gains on L**, a candidate for the "stochastic control" | H | model §8.4 |
| **No mode-line broadening** beyond the LTI decay at about 3 Hz resolution; modulation slower than about 1 Hz is not excluded | S (limit) | model §8.6 |
| Noise model (in band): L −103 dB at 1.4 s, drifting −1.8 dB/s; R −110 dB, flat. The noise spectrum rises with frequency. Most of L's "decaying floor" is below 20 Hz | S | model §9.3 |

## 6. Memory and compute budget (derived)

- 64 K words of delay RAM. Maximum pre-delay is 26,347 words, leaving **about 39 K words** for the loops, cluster and diffusion.
- S's loops use 5,758 words. M is expected at about 11.5 K; L is memory-capped at about 30–39 K.
- 364 instructions per sample in total, i.e. about 50–70 weighted delay taps for the whole machine.

All of this is **H**, derived from the parts and TI timings.

## 7. Open questions (ranked)

All of these need library IRs. The tests for each are in `NEXT_STEPS.md` §2.2.

1. **The remaining injection/read-out structure.** Sparse fits plateau at about −23 dB where a matching synthetic reaches −46 dB, so something structural is still missing. The diffuser (AP20/72/230 + comb460) is found; the late input is explained.
2. **What the "stochastic control" is.** Injection randomisation is ruled out (R is LTI to about 3 %). The remaining candidate is about 1 % slow gain jitter on L's read-out, or the capture chain. Library test: is it L-only and constant across captures?
3. **Exact step tables** for RT, LF, HF, AMPL and LP; confirming the LF reference (the card says 50 Hz); the HF ×1.0 behaviour.
4. **Loop sets and cluster patterns for M and L.** The 1986 review says the early reflections differ per program.
5. **The LP's filter type and order** (it is downstream per Studio Magazin 96), and the loop filters' type and order (from LF/HF sweeps).
6. **Fixed-point details** (rounding/truncation, coefficient grid), needed for bit-exactness.
7. **N's structure** (FIR against truncated loops; how reversal works).

## 8. Source index

- **Documents** (local copies of the key sources are in `sources/`, which is gitignored):
  - **the EMT 246S service manual** (2 PDFs, scans, bilingual DE/EN; local in `sources/emt-documents/`): schematics and descriptions for every board incl. the processor board 7 246 123, self-test programs, calibration, technical data. Not yet read in detail (`NEXT_STEPS.md` §2.6);
  - the EMT 245 service manual: lineage only. It describes **19 loops in mixed series/parallel** on several parallel microprocessors; the 246 is the 245's direct successor, but its software is "in keinerlei Beziehung" to it;
  - data sheet (your PDF; emt-archiv.de);
  - manual pp. 65–74 (remote boards);
  - parameter card (studiorepair.com);
  - operation card scan and front-panel photo (`irs/EMT 246 captures/`, from the Tomkin library; gitignored);
  - Studio Magazin 96 and 99 (1986);
  - EMT Courier (Jul 1985);
  - EMT 252 brochure;
  - TI TMS32010 User's Guide (1985);
  - US 4,181,820 (the EMT 250 patent).
- **Photos and video:**
  - studiorepair.com internals gallery;
  - panel photos;
  - two CLA setting images;
  - Foleytronics video PZj4ZSpi9P8.
- **Measurements:**
  - `analysis/shortRT/summary.json`;
  - `tools/topology_tests.py`;
  - `docs/model_s_results.md`;
  - harness on the synthetic library.
