# The reverb engine

How the real algorithm replaces the placeholder in `src/dsp/ReverbEngine.h`. The research side (what each program
is and why) is in `../../docs/`; this file is the plugin side's contract: the signal path, who owns which file, and
the tests every part must pass. The rules in `DEVELOPMENT.md` apply throughout (cross-platform C++, no `M_PI`,
real-time code allocation- and lock-free, warnings are errors in spirit).

## Signal path

```
host rate                         device rate (19.2 MHz / 728 = 26,373.6 Hz)                       host rate
mono in -> input analogue chain -> [A/D clip, 16-bit mode] -> ProgramCore -> (L, R) -> output chain -> wet L, R
           (input filter,          sampled at the device       S / M / L / E / N     S&H lag per channel,
            coupling high-passes)  instants                    + engine: DELAY,      output filter, L half a
                                                                LOW PASS             device sample late
```

- **Device rate inside, host rate outside.** Each program core runs at the device rate with integer delays, exactly
  as its spec and reference IRs. `DeviceBridge` converts host -> device and device -> host and carries the analogue
  chain (`../../docs/NEXT_STEPS.md` §4 and §11, `../../tools/analog/`), so the cores never see the host rate.
- **Mono in, stereo out.** One core per engine; L and R come out at the same device instants. L's half-sample lag
  (the timing PROM's S&H pulses) is the output side's job: evaluate L's interpolation half a device sample later.
- **The A/D clip moves to the device rate.** Once the bridge produces the input filter's output at the device
  instants, the 16-bit mode clips it there (full scale from CAL, as `AdcClipper` does now) and the reinjection scheme
  in `AdcClipper.h` is no longer needed. The bridge exposes the clip and an overload flag; the processor wiring is the
  integration step's.
- **Engine-level controls**, not the cores': DELAY (a pure delay of the whole output, device samples), LOW PASS (a
  guessed device-rate low-pass on the reverb path, H, with 12 kHz an exact bypass), program switching (crossfade),
  the mono fold-down `(L + R) / sqrt(2)`, and the latency report (bridge + core look-ahead + the A/D stage).
- **Always on:** the analogue chain and the time variance (M / L's theta, E's moving reads). Only the null tests
  switch the time variance off (`ProgramCore::setModulation(false)`).

## Files and owners

| Path | What | Built in |
|---|---|---|
| `src/dsp/engine/ProgramCore.h` | `Settings`, `deviceRate`, the core interface | phase 0 |
| `src/dsp/engine/DeviceBridge.h` | host <-> device conversion and the analogue chain | phase 1 |
| `src/dsp/programs/SCore.h` | program S | phase 2 |
| `src/dsp/programs/MLCore.h` | programs M and L (one class, a program flag) | phase 2 |
| `src/dsp/programs/ECore.h` | program E (the EMT 250 algorithm) | phase 2 |
| `src/dsp/programs/NCore.h` | program N (non-linear) | phase 2 |
| `src/dsp/generated/<prog>_spec.h` | constant tables from the spec JSON (committed) | phase 2 |
| `tools/plugin_codegen/<prog>.py` | the generator for that header; `common.py` holds the shared writer | phase 2 |
| `tests/dsp/<Prog>Tests.cpp` | that part's tests; `RefIr.h` holds the shared helpers | phases 1, 2 |
| `src/dsp/ReverbEngine.h` | dispatch, engine-level controls, 16-bit wiring | phase 4 |

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
- Every write into recirculating storage (loop and all-pass delay lines) goes through `ProgramCore::store()`, the
  one place the 16-bit mode will quantise (phase 3). Non-recirculating taps need not.
- `latency()` reports any look-ahead (N's centred branch read kernel). Time variance starts from a fixed phase at
  `reset()`, so renders are repeatable.
- Code reads like the surrounding plugin code: short comments that say what and why, with the spec section in
  brackets where a number comes from (e.g. `// model_s_results §38.1`).

## Generated data

`tools/plugin_codegen/common.py` writes `src/dsp/generated/*.h` with doubles that round-trip exactly (`repr`). Each
program's generator reads only its committed spec JSON (`analysis/s/s_plugin_spec.json`,
`analysis/ml/ml_plugin_spec.json`, `analysis/e250/e_plugin_spec.json`, `analysis/nl/nl_plugin_spec.json`). Laws stay
laws: generate their constants and evaluate them in C++, rather than tabulating every panel step.

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
engine also adds analogue hiss calibrated to the spec's 73 dB S/N (H: the internal scaling is unknown). Phase 3.

## Decisions taken for the build (30 October 2026)

- Panel steps between the library's: the data sheet's 16 geometric RT steps, with each spec's law between the
  library steps (H). The default RT step is exactly 1.0 s. N's tables come from its spec (EARLY: 16 values, LATE:
  0-2.0 s by 0.2, SPAN: 0-0.5 s).
- LOW PASS: a guessed device-rate low-pass on the reverb path (H), 12 kHz = exact bypass.
- N's branch read kernel: on, as the spec adopts it; NL_BK=1 (the main taps' read kernel) stays off.
- E's moving reads: as the spec adopts them (`e_sim.ADOPTED`): rates 0.553 Hz x (1 +- 0.5) per read, a triangle
  (E's data cannot tell it from a sine; the plan's "sine" was a misreading, corrected 30 October).
- 16-bit mode adds the analogue hiss (73 dB S/N); FLOAT stays clean.
- Resampling: the bridge adds no error above -90 dB to the modelled chain over 20 Hz-12.8 kHz; its latency is
  reported.
