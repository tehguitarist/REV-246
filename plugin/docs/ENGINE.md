# The reverb engine

How the real algorithm runs in the plugin (`src/dsp/ReverbEngine.h`). The research side (what each program
is and why) is in `../../docs/`; this file is the plugin side's contract: the signal path, who owns which file, and
the tests every part must pass. The rules in `DEVELOPMENT.md` apply throughout (cross-platform C++, no `M_PI`,
real-time code allocation- and lock-free, warnings are errors in spirit).

## Signal path

```
host rate                         device rate (19.2 MHz / 728 = 26,373.6 Hz)                       host rate
mono in -> input analogue chain -> [A/D clip, 16-bit mode] -> alignment -> ProgramCore -> DELAY -> LOW PASS
           (input filter,          sampled at the device      delay        S / M / L / E / N
            coupling high-passes)  instants                                (two in a program crossfade)
        -> (L, R) -> output chain -> [hiss, 16-bit mode] -> wet L, R
                     (per channel: correction, S&H lag, output filter; L 0.5037 device samples late)
```

- **Device rate inside, host rate outside.** Each program core runs at the device rate with integer delays, exactly
  as its spec and reference IRs. `DeviceBridge` converts host -> device and device -> host and carries the analogue
  chain (`../../docs/NEXT_STEPS.md` §4 and §11, `../../tools/analog/`), so the cores never see the host rate.
- **Mono in, stereo out.** One core per engine; L and R come out at the same device instants. L's half-sample lag
  (the timing PROM's S&H pulses) is the output side's job: evaluate L's interpolation 0.5037 device samples later (the
  library's L - R offset, `../../analysis/plugin_check/REPORT.md` §4: 0.5 plus a fixed L / R chain-delay difference).
- **No oversampling anywhere** (31 October 2026). The host -> device conversion is exact for the band-limited signal
  the host samples stand for, and the unit's own input filter is the anti-alias filter: whatever it passes above
  13.19 kHz folds, as in the unit. The output side reads the analogue output filter's response after the captures'
  anti-alias sinc, so it adds no images. The A/D clip and the 16-bit cut are non-linear at the device rate, as the
  unit's were (their folding is part of the sound). The one non-linear stage at the host rate, XFMR, distorts the core
  flux (the integral of the signal), so its harmonics fall steeply with frequency (0.0001 % at 1 kHz): what could
  fold from the top octave is far below -120 dB.
- **The chain correction** (`../../analysis/plugin_check/REPORT_chainfix.md`): the netlists miss the measured chain
  (chain3) in the output filter's skirt above 8 kHz (up to 70°) and by about a degree in the mid band. Per channel, on
  the device samples before the S&H lag, the bridge runs a first-order 1.8 Hz high-pass (H: Tomkin's 6.25 Hz total over
  the schematic's 4.46 Hz, firmware DC blocker or capture rig, `../../docs/analog_results.md` §3) and a 65-tap FIR
  (`chain_fix.json` via `tools/plugin_codegen/bridge.py`) with a bulk delay of 17 device samples, compensated in the
  latency. Always on; no parameter.
- **The A/D clips at the device rate**, on the input filter's output at the device instants (`DeviceBridge`): in
  16-bit mode a hard clip at full scale, in both modes the overload flag that lights OVL. Full scale follows CAL
  (`src/dsp/AdcCalibration.h`: a sine 6 dB over 0 VU), and the same number scales the cores' 16-bit storage and the
  hiss. (The host-rate `AdcClipper` and its reinjection scheme, which stood in for this before the bridge, are
  retired.)
- **Engine-level controls**, not the cores': DELAY (a pure delay of the whole output, device samples), LOW PASS (a
  guessed device-rate low-pass on the reverb path, H, with 12 kHz an exact bypass), program switching (crossfade) and
  the latency report (bridge + the common alignment). The mono fold-down `(L + R) / sqrt(2)` is the processor's,
  after each channel's output stage.
- **Always on:** the analogue chain and the time variance (M / L's theta, E's moving reads). Only the null tests
  switch the time variance off (`ProgramCore::setModulation(false)`).

## Files and owners

| Path | What | Built in |
|---|---|---|
| `src/dsp/engine/ProgramCore.h` | `Settings`, `deviceRate`, the core interface | phase 0 |
| `src/dsp/engine/DeviceBridge.h` | host <-> device conversion, the analogue chain, the A/D clip and overload flag | phase 1 |
| `src/dsp/engine/BridgeFilter.h` | the bridge's exact band-limited discretisation of one analogue filter (tables + one recursion per pole) | phase 1 |
| `src/dsp/generated/bridge_spec.h` | the analogue chain's poles and residues (from `tools/plugin_codegen/bridge.py`) | phase 1 |
| `src/dsp/engine/AnalogueHiss.h` | the 16-bit mode's analogue hiss on the wet outputs (host rate) | phase 3 |
| `src/dsp/programs/SCore.h` | program S | phase 2 |
| `src/dsp/programs/MLCore.h` | programs M and L (one class, a program flag) | phase 2 |
| `src/dsp/programs/ECore.h` | program E (the EMT 250 algorithm) | phase 2 |
| `src/dsp/programs/NCore.h` | program N (non-linear) | phase 2 |
| `src/dsp/generated/<prog>_spec.h` | constant tables from the spec JSON (committed) | phase 2 |
| `tools/plugin_codegen/<prog>.py` | the generator for that header; `common.py` holds the shared writer | phase 2 |
| `tests/dsp/<Prog>Tests.cpp` | that part's tests; `RefIr.h` holds the shared helpers | phases 1, 2 |
| `src/dsp/ReverbEngine.h` | the engine: bridge, the five cores, alignment, program crossfade, DELAY (`OutputDelay`), LOW PASS (`DeviceLowPass`), hiss | phase 4 |
| `src/dsp/AdcCalibration.h` | the A/D full scale for a CAL setting (`adcFullScaleForZeroVu`) | phase 4 |
| `tests/dsp/EngineTests.cpp` | the engine end to end at host rates | phase 4 |
| `tools/render/RenderIr.cpp` | offline renders for the hardware check (`tools/plugin_check.py`), through the engine's offline hooks (`setModulation`, `processDeviceRate`) | phase 5 |

Cores are header-only (or header + one `.cpp` listed in both CMake source lists if one is truly needed). Each part
touches only its own files; shared files (`ProgramCore.h`, `RefIr.h`, `CMakeLists.txt`, `common.py`) change only in
the integration step, so the parts can be written in parallel.

## The core contract (`ProgramCore`)

- `prepare()` allocates for the worst case of the program (the longest delays at any setting). `reset()`,
  `setSettings()` and `process()` never allocate or lock.
- `setSettings()` maps the panel values (`Settings`: LF x, RT, HF x, AMPL; N's EARLY / LATE / SPAN in their slots,
  see `params::joyValue`) to coefficients with the spec's own law, **for any value in range**, not only the library's
  steps. The resolved settings in the spec (S's 128, M's 160, L's 180, ...) are test vectors for that law.
- `process()` is mono in, (L, R) out, float I/O. Internal state may be float or double, whichever the null bar needs;
  say which and why in the header.
- Every write into recirculating storage (loop and all-pass delay lines) goes through `ProgramCore`'s storage cut
  (`store()`, or `storeAcc()` / `runShelfPole16()` / `cut16()` for a value held in the accumulator), where the
  16-bit mode quantises. Non-recirculating taps need not.
- `latency()` reports any look-ahead (N's centred branch read kernel) left after `setAdvance()`, which shifts the
  whole output earlier where every path is deep enough (`maxAdvance()`; E and N, see The engine). Time variance starts from a fixed phase at
  `reset()`, so renders are repeatable.
- Code reads like the surrounding plugin code: short comments that say what and why, with the spec section in
  brackets where a number comes from (e.g. `// model_s_results §38.1`).

## Generated data

`tools/plugin_codegen/common.py` writes `src/dsp/generated/*.h` with doubles that round-trip exactly (`repr`). Each
program's generator reads only its committed spec JSON (`analysis/s/s_plugin_spec.json`,
`analysis/ml/ml_plugin_spec.json`, `analysis/e250/e_plugin_spec.json`, `analysis/nl/nl_plugin_spec.json`). What a law needs from outside (S's RT folds, cells and HF pole curve; E's loop
law, test vectors and moving-read draw) is embedded in the spec by its tool (`s_spec.py` / `e_spec.py`, `embed_codegen_data`).
Laws stay laws: generate their constants and evaluate them in C++, rather than tabulating every panel step.

## Tests

Build the light runner on its own (about 1 GB peak, under a minute clean):

```bash
cmake -S plugin -B plugin/build-dsp -DCMAKE_BUILD_TYPE=Release -DRV_BUILD_DSP_TESTS=ON -DRV_DSP_TESTS="S"
cmake --build plugin/build-dsp --target REV246DspTests --parallel 2
plugin/build-dsp/REV246DspTests_artefacts/Release/REV246DspTests
```

`RV_DSP_TESTS` picks `tests/dsp/<prefix>*Tests.cpp` (empty = all). `RefIr.h` loads a reference IR, renders a core's
impulse response with its latency removed, and nulls them in the spec tools' windows (0-50, 50-150, 150-300,
300-600, 600-1000 ms; error energy over the reference's, per channel).

**Pass bars (phase 2).** The reference IRs are the Python models' own output, so a correct float port nulls far
below any audible or measurable level:

- every reference IR, every window, both channels: **-100 dB or lower** (report the actual numbers; expect -120 or
  better in double, and say where float storage sets the floor);
- the coefficient law reproduces every resolved setting in the spec (filters to about 1e-9 relative, or say why not);
- M / L and E with time variance on: the statistics the spec gives (theta's rms and rate; E's read motion) are met,
  and theta = 0 / frozen reads null as above.

The full-chain check against the Tomkin hardware IRs (-50 dB held out, excess over the noise floor, parameter
count; `../../CLAUDE.md`) is phase 5 and runs through the bridge at host rates.

## 16-bit mode

Float by default. With `Settings::sixteenBit` the cores' `store()` quantises as the unit's 16-bit storage did; the
reference is `../../tools/loop_read/c/s_loops_fx.c` and `../../docs/model_s_results.md` §36.5. In 16-bit mode the
engine also adds analogue hiss at 91 dB S/N re nominal (the user's choice, 3 dB better than the data sheet's 88 dB
pass-through figure: a very good unit), about -109 dBFS rms at the default CAL. The REVERB programs' total (spec
73 dB) is a check, not a target: its extra 15 dB is the 16-bit storage noise the cores make by themselves.

The cut is **magnitude truncation (toward zero)**, the user's choice of 30 October (H for the unit). The scheme is
bit-exact with `s_loops_fx.c` in its floor and rounding modes (`ProgramCore::setStorageCut`), but with those a loop
can fall into a dead band that holds a loop-periodic pattern for ever after the input stops (up to -28 dB re nominal
at L RT 20, -41...-52 dB in S at HF 1.0; idle S/N 47-72 dB, noisier than the spec). Toward zero keeps the decays'
16-bit grit (-11...-48 dB against float at 50-300 ms) and always decays to silence, so idle noise is the hiss alone.
Measurements: `tests/dsp/SixteenBitTests.cpp`.

## The engine (`ReverbEngine`)

One per mono input; it owns a `DeviceBridge`, the five cores (S, M, L as two `MLCore`s, E, N), all prepared up front
for their worst case, so nothing allocates after `prepare()` (EngineTests counts allocations in `process()`: zero).
Per device block the input goes into a ring; each running core reads it through its alignment delay; the output goes
through DELAY and LOW PASS; the bridge turns it into host-rate L and R; in 16-bit mode the hiss is added.

**One time grid, one latency.** Every program comes out on the same grid, so the reported latency is the same for all
of them, both modes and both engines, and a program switch never shifts the timing. Two cores are late against the
common origin: N by its branch kernel's look-ahead (40 device samples), and E because the spec puts E's R direct tap
(k0 = 64) with M's first cluster tap (at M's kernel centre, 10) at DELAY 0 (`e_plugin_spec.json` controls.onset,
`../../docs/EMT250_program.md` §17.8, S): 54 samples. Both are pure delays inside the core, so instead of delaying the
other programs by 54 (the 30 October build), E and N read their input that much sooner (`ProgramCore::setAdvance`,
31 October): E's every path starts at its front delay k0; N's front taps start at 41 and its branch reads at 685, so
the kernel's look-ahead fits inside them (`static_assert`s in `NCore::prepare`). An exact time shift
(EngineTests "E's and N's advance is an exact time shift": bit-identical; with E's moving reads on, their fixed
start phase now meets the input 54 samples later, an equally valid realisation). The alignment is 0, and

  latency = bridge latency + round(17 x host rate / device rate) host samples

(17: the bridge's chain correction, `DeviceBridge::correctionDelay()`; the fraction left, under half a host sample,
is not compensated). On R: L is another half device sample later.

| Host rate | Bridge | Correction | Reported latency (30 October: with the 54-sample alignment) |
|---|---|---|---|
| 44.1 kHz | 21 | 28 | 49 samples, 1.11 ms (140, 3.17 ms) |
| 48 kHz | 19 | 31 | 50 samples, 1.04 ms (148, 3.08 ms) |
| 96 kHz | 14 | 62 | 76 samples, 0.79 ms (272, 2.83 ms) |

What is left: the correction FIR's bulk delay (17 device samples, 0.64 ms; a refit with a shorter one is the
candidate) and the bridge's band-limiting sinc (0.15-0.48 ms), which stays.

The modelled analogue chain adds its own peak delay on top (about 0.17 ms: S's cluster tap 0 peaks 8 samples after
the latency at 44.1 kHz), as the unit does; it is part of the sound, not of the latency.

**Program switch.** The new core is reset and its *input* fades in over 30 ms while the old core keeps running and its
*output* fades out (equal power, the gains a sin / cos of a sin^2 ramp so neither has a kink). Fading the new core's
input rather than its output matters: a reset core fed an abrupt signal clicks at each of its internal taps as the
onset reaches them (S -> M under a 200 Hz sine: 2-8 kHz peak 7.6x the steady level with an output fade, 1.2x with the
smooth input fade, i.e. at the floor). A switch requested during a fade waits for it to end (at most 30 ms).

**Settings within a program apply at once**, as the unit's coefficient loads did: the cores' filters are first- and
second-order sections with fixed poles inside the unit circle at every setting, so a coefficient step cannot blow up,
and the unit had no smoothing either. Stepping LF, RT, HF and AMPL through every step every 5.3 ms under noise
peaks at 1.0-1.5x the default's level (N and E most, from their louder settings), never anything like a blow-up. N's
SPAN and LATE switch whole stages and clear the late loops at LATE 0, which clicks as the unit's would.

**DELAY** is an integer number of device samples, `round(ms x 26.3736)`, after the core. A change crossfades from
the old read position to the new over 10 ms (sin^2 / cos^2): the unit's read pointer would jump, a click per step in
a running tail, and a held lever steps every few milliseconds. Changes during a fade are taken in turn.

**LOW PASS** (H, `e_plugin_spec.json` controls.LOW_PASS): a 2nd-order Butterworth at the device rate, bilinear with
the cutoff prewarped, on the whole reverb output (the operation card draws it after the reflections). 12 kHz is open:
an exact bypass. Into and out of open it fades over 10 ms; between two cutoffs the coefficients switch at once.

**Cost** (31 October 2026, after the optimisation pass; Apple-silicon Mac, Release, block 512, 10 s of noise, best of
3, % of one core in real time). Hidden test cases, run on demand: `REV246DspTests "[.profile]"` (per stage) and
`REV246Tests "[.profile]"` (the whole `processBlock`).

Whole processor (one engine: stereo in / stereo out on the original summed path; TRUE ST doubles the engine part),
before the pass in brackets:

| Program | 44.1 kHz float / 16-bit | 48 kHz float / 16-bit | 96 kHz float / 16-bit |
|---|---|---|---|
| S | 0.52 / 0.65 (1.12 / 1.20) | 0.54 / 0.67 (1.12 / 1.21) | 0.81 / 0.99 (1.48 / 1.61) |
| M | 0.58 / 0.67 (0.92 / 1.01) | 0.60 / 0.70 (0.93 / 1.02) | 0.88 / 1.02 (1.29 / 1.42) |
| L | 0.64 / 0.72 (1.10 / 1.19) | 0.66 / 0.74 (1.12 / 1.21) | 0.94 / 1.06 (1.47 / 1.60) |
| E | 0.53 / 0.62 (0.82 / 0.90) | 0.55 / 0.64 (0.82 / 0.92) | 0.82 / 0.96 (1.17 / 1.32) |
| N | 0.69 / 0.78 (1.12 / 1.20) | 0.71 / 0.80 (1.13 / 1.22) | 0.98 / 1.12 (1.48 / 1.62) |

XFMR on and off are the same within 0.01 %: the stage always runs (warm state for the toggle's crossfade), so XFMR off
costs as much as on. Per stage (before the pass in brackets):

| Stage | Cost |
|---|---|
| bridge alone (analogue chain, correction, rate conversion; identity core) | 0.36 % at 44.1 / 48 kHz, 0.54 % at 96 kHz (0.64, 0.89) |
| core S / M / L / E / N, float (device rate) | 0.06 / 0.13 / 0.18 / 0.07 / 0.26 % (0.37 / 0.18 / 0.36 / 0.07 / 0.40) |
| the same, 16-bit | 0.15 / 0.17 / 0.22 / 0.12 / 0.31 % (0.41 / 0.22 / 0.40 / 0.12 / 0.45) |
| processor around the engine (gains, XFMR, width, meters, dry delay) | about 0.1 % (0.2 % at 96 kHz) |

The bridge is now 50-70 % of an engine. A program crossfade runs two cores for 30 ms.

**The optimisation pass** (31 October 2026; each step nulled against renders of every program at 44.1 / 48 / 96 kHz,
float and 16-bit, from the build before it: worst -193 dB, i.e. float rounding; every DSP test bar unchanged):
- Long dot products in four (two) partial sums instead of one running sum, which the compiler may not reorder without
  fast-math (`dotProduct`, `ProgramCore.h`): the correction FIR, the bridge's head rows, M / L's and N's output
  kernels (their histories doubled so the kernel reads contiguously).
- S's read-out taps summed per channel in registers: the old `acc[tapChannel[t]] +=` sent every add through memory
  in one dependency chain (S's core 0.37 -> 0.06 %).
- The bridge's L read-out: its tail weights are R's times a fixed per-mode factor (L is a constant 0.5037 device
  samples later), one complex product per mode instead of a second table lookup and series.
- Not done: the 2-biquad + 33-tap correction (`REPORT_chainfix.md` §2) would save under 0.02 % now that the 65-tap FIR
  costs about 0.035 %, for an accuracy loss (-51 dB at 11-12.3 kHz); skipping XFMR while off saves nothing measurable;
  the 16-bit cut sits in S's loop recursion (its serial chain), so its 0.09 % stays.

## Decisions taken for the build (30 October 2026)

- Panel steps: the library's RT values are steps (C); the range ends are the data sheet's; between them (H) the widest
  gaps in log RT are split at the geometric mean, rounded to the display's 0.1 s, until there are 16 (tables in
  `src/params/Parameters.cpp`). Each spec's law runs between the library steps. The default RT step is exactly 1.0 s.
  N's tables come from its spec (EARLY: 16 values, LATE: 0-2.0 s by 0.2, SPAN: 0-0.5 s); LATE's 11 values make the HF
  parameter 0-10, with S / M / L / E stopping at step 9.
- LOW PASS: a guessed device-rate low-pass on the reverb path (H), 12 kHz = exact bypass. Kept by the user's decision (31 October 2026) after the Altiverb presets gave no lever law: their darker M / L / E presets differ in the output kernel's 2-tap coefficient, not by a filter after the program; PRESET23's steep linear-phase cut is one file, either the lever at 8 kHz or the Altiverb's processing (`../../docs/library_notes.md` 7.x). No preset-only kernel parameter (user, 31 October).
- N's branch read kernel: on, as the spec adopts it; NL_BK=1 (the main taps' read kernel) stays off.
- E's moving reads: as the spec adopts them (`e_sim.ADOPTED`): rates 0.553 Hz x (1 +- 0.5) per read, a triangle
  (E's data cannot tell it from a sine; the plan's "sine" was a misreading, corrected 30 October).
- 16-bit mode adds the analogue hiss at 91 dB S/N re nominal (about -109 dBFS rms at the default CAL); FLOAT stays
  clean. The REVERB programs' 73 dB total is a check, not a target.
- Resampling: the bridge adds no error above -90 dB to the modelled chain over 20 Hz-12.8 kHz; its latency is
  reported.
- (Integration, 30 October) One latency for every program: all cores aligned to the latest (E's onset, 54 device
  samples; since 31 October E and N advance instead, so the alignment is 0). Program switch: 30 ms, the new core's input faded in, the old core's output out. Settings within a
  program: at once. DELAY changes: a 10 ms crossfade. The A/D clip moved into the bridge; `AdcClipper` retired.
