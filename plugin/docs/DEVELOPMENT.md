# Development rules (plugin)

## Cross-platform: every line must build on macOS, Windows (MSVC) and Linux (GCC)
- **No `M_PI`, `M_E`, `M_SQRT2` and friends.** MSVC only defines them if `_USE_MATH_DEFINES` is set before `<cmath>`,
  which is fragile. Use `juce::MathConstants<float>::pi` (or `std::numbers` once C++20). The test targets define
  `_USE_MATH_DEFINES` as a safety net, not as a licence to use the macros.
- Include what you use (`<cmath>`, `<array>`, `<tuple>`, `<functional>`, `<iterator>`): libc++ often leaks them
  transitively, MSVC and libstdc++ do not.
- No `#pragma` or compiler-specific attributes without a `JUCE_` / `#if` guard; no `__attribute__`, no VLAs, no
  `ssize_t`, no `M_`-anything, no POSIX-only calls. Use `juce::File` for paths (no hard-coded `/` or `\\`).
- Warnings are on (`-Wall -Wextra`, `/W4`): treat a new warning on any platform as a build break.
- Fonts: only system fonts the code asks for by name with a fallback (`Arial`, then `Liberation Sans` on Linux), or
  embedded ones. The panel's digits are drawn as shapes, not a font.
- File-name case matters on Linux: the embedded images are looked up by their exact file name
  (`ui-layout.csv` -> `Layout.h` -> `assets/images`).
- Real-time code (`ReverbEngine::process`, `processBlock`) is allocation- and lock-free; no `std::string`,
  `std::vector` growth, or `juce::String` there.

## Testing in a DAW
- Configure a build with `-DRV_INSTALL_AFTER_BUILD=ON` and every build copies the AU and VST3 into the user's plug-in
  folders (`~/Library/Audio/Plug-Ins/Components` and `.../VST3` on macOS), which DAWs scan. The system folders
  (`/Library/Audio/Plug-Ins/...`) are root-owned; copy there by hand with `sudo` if a host only looks there. After an AU
  update, `killall -9 AudioComponentRegistrar` makes macOS rescan it.

## CI/CD
- Workflows are in `../.github/workflows/plugin-ci.yml` and `plugin-release.yml`. They are **`workflow_dispatch` only**
  until the reverb is ready, so they use no GitHub minutes. To switch CI on, uncomment its `push` / `pull_request`
  triggers.
- Every `uses:` action must run on the current Node runtime (all four currently declare `node24`: `actions/checkout@v6`,
  `actions/upload-artifact@v6`, `actions/download-artifact@v8`, `softprops/action-gh-release@v3`). When touching a
  workflow, check the action's `action.yml` (`runs.using`) and bump the major version rather than pinning an old
  Node (a stipulated `node20` stops working once the runners drop it).
- Signing uses the repo secrets already set (`APPLE_*`); the release job is `workflow_dispatch` and needs a green CI
  run on the same commit unless `skip_ci_check` is ticked.

## Assets
- `../ui/` holds the masters and is never modified. `tools/build_assets.sh` writes `assets/images/` (1:1 with the
  3495 x 1530 design space = 2.5x the default 1398-wide editor, pngquant-compressed) and `src/ui/Layout.h`.
- Fader / funnel / knob behaviour described in the designer's notes lives in `src/ui/Design.h` with the numbers
  measured from `main.png` (tick positions, the funnel's curved edges).

## Parameters
- IDs in `src/params/Parameters.h` are permanent. The joystick step tables (`src/params/Parameters.cpp`): step counts
  and ranges are the data sheet's (`../docs/EMT246_unit_facts.md`); the Tomkin library's values are steps (C: RT
  S 0.3 0.6 0.8 1.0 1.4 1.8 2.2 2.5, M 0.6 0.8 1.2 1.8 2.6 3.0 4.0 5.0, L 2.2 2.6 3.6 4.5 6 8 10 15 20; E's the
  EMT 250 panel's 0.4-4.5 (16 steps; H that the 246's E uses them); LF x 0.1-1.2, 1.3, 1.5, 1.8, 2.0 (above 1.2 H); every library LF and HF); the rest are guesses (H): RT fills
  the widest log gaps at their geometric means, rounded to 0.1 s. N's columns are its spec's: EARLY RT's 16 values
  in the LF slot (10 = flat), LATE RT 0-2.0 s by 0.2 in the HF slot, SPAN 0-0.5 s by 0.05 in the AMPL slot. LATE has
  11 values, so the HF parameter runs 0-10 and S / M / L / E stop at step 9 (`joyMaxStep(Program, Joy)`; a higher step
  left over from N reads as their 1.0).
- **Pre-release range change (flagged 31 October 2026):** the HF parameter's range grew from 0-9 to 0-10 when N's
  LATE steps went in. Saved state and presets are unaffected (the APVTS stores the step itself), but a host stores
  automation as 0-1, so an HF lane recorded with an earlier build lands on a different step (old step k at k/9 now reads
  round(10k/9)). No build has been released, so nothing is migrated; after the first release a step table's range must
  not change again (add steps only through a new parameter ID or a versioned state migration).
- The buses allow mono in / mono out, mono in / stereo out and stereo in / stereo out. The reverb engine is mono in /
  stereo out, so a mono output takes its left side. For stereo in / stereo out, `trueStereo` (default off, the
  ORIGINAL / TRUE ST button) selects between the original summed path (0.5 x L+R into one engine) and two engines,
  one per input channel. Presets skip it, like `power`, `link` and `transformer`.
- `transformer` (default on, the XFMR button) switches on `src/dsp/Transformer.h`: an input transformer core
  on the mono input, then per output the low-band EQ (`src/dsp/CaptureEq.h`) and an output transformer core. Presets
  skip it, like `power` and `link`. The stage always runs (warm state, a 20 ms crossfade when toggled); off, its
  weight is exactly 0 and the signal is unchanged.
- The EQ is a 13-number cascade fitted to CLA's capture transfer (`../analysis/unit_chain/cla13.json`; a first-order
  high-pass, an RBJ low shelf and three RBJ peaking sections, the formulas of `../tools/unit_transfer.py`: `hp1_sos`,
  `lowshelf_sos`, `peak_sos`). It replaced `bq.json`'s 10-number "CLA" entry on 31 October 2026 by the user's decision
  (`../docs/model_s_results.md` §37.1-37.3): 6-12 dB deeper nulls at 160-400 Hz and about 2 dB at 400-800 Hz, for
  about 1 dB at 80-160 Hz and 0.8-1.6 kHz; the two differ by at most 0.17 dB. A joint CLA + Altiverb cascade was tried
  and loses 6-7 dB at 20-80 Hz (the units differ in phase there). It was fitted at the device rate; the plugin designs
  the same sections at the host rate from the same parameters, which matches the Python response within 0.001 dB at
  50 Hz-2 kHz from 22.05 to 192 kHz (no prewarping needed; the fit is undefined below ~50 Hz). If the fit changes,
  update `CaptureEq::cla` and the reference table in `tests/TransformerTests.cpp` together.
- `stereoWidth` (default 100, the WIDTH halo) is a sound setting and part of presets. `src/dsp/StereoWidth.h` holds
  the level-matched M/S `StereoWidener` (100 is a bit-exact bypass, 0 mono, 200 wide, bass mono above 100). It runs
  on the whole wet: pure M/S gains move nothing in time, so the early reflections are never smeared (the design
  intent). The design's late-tail decorrelator was dropped (31 October 2026): the unit's late L / R are already
  uncorrelated (|r| <= 0.04 over 0.1-1 s in every program; up to 0.39 only at 100-400 Hz, which the bass-mono
  crossover handles), so it adds no width; and its Schroeder all-passes became pure delays at amount 0, so engaged it
  put the tail ~8 ms behind the early reflections. An early / late split for it would have cost about 0.25 % of a core.
- `sixteenBit` (default on, the FLOAT / 16-BIT button) is a sound setting and part of presets. It reaches the engine as
  `Settings::sixteenBit`, and switches three things at once: the cores' 16-bit storage (`ProgramCore::store()`), the
  A/D's hard clip, and the analogue hiss (`src/dsp/engine/AnalogueHiss.h`, 91 dB S/N re nominal, about -109 dBFS rms at
  the default CAL). The A/D is inside the engine (`src/dsp/engine/DeviceBridge.h`): the unit's converter is a SAR,
  which sticks at the end of its range, so in 16-bit mode each held sample of the input filter's output is clipped at
  the device rate (the input filter lifts 10-11.5 kHz by up to 5.1 dB, so HF clips early; the clip's harmonics fold
  around 13.2 kHz). Switching it is instant, like the unit's coefficient loads; it never changes the latency.
- The A/D's full scale follows the meter calibration (`Rev246Processor::setMeterCal`, the CAL button; UI state, kept
  with the session, not in presets): the calibration stands for the unit's input alignment, and its full scale is a
  sine 6 dB over 0 VU (unit facts §2; `src/dsp/AdcCalibration.h`). The engine passes it to the bridge's clip, the
  cores' 16-bit storage and the hiss (`ReverbEngine::setFullScale`). OVL lights at +6 VU, on a full-scale input, or
  when the A/D is out of range (`ReverbEngine::takeOverload`; in float mode too: where it would clip), as the unit's
  REG LED did.
- Latency: the engine's (`ReverbEngine::latency`, `docs/ENGINE.md`): the same for every program, both modes and both
  engines (49 / 50 / 76 samples at 44.1 / 48 / 96 kHz, about 1 ms). The dry and power-off signals are delayed to match.
- Display fields other than the program display open a typed-entry dialog on double-click. `params::joyStepFromValue`
  validates the text for the current program and returns a step; `joyInputHint` supplies the unit/range/list shown in
  the dialog. Discrete controls (N's EARLY and LATE included) require an exact step, while RT (S/M/L/E) accepts any
  value in range and snaps to the nearest step (in log). RT in N is rejected with an information dialog. The same `setJoyStep` path is
  used, so the levers, funnels and displays refresh.
- Initial/factory defaults: S, the CLA reference levers (LP 12 kHz, DELAY 10 ms, AMPL 80 %, LF x1.0, RT 1.0 s, HF
  x0.9), MIX 100 %, INPUT/OUTPUT 0 dB, LINK on, 16-bit on, XFMR on, WIDTH 100, ORIGINAL stereo mode and power on.
  Factory presets ("Init S" ... "Init N", one per programme) reset the stored settings to these defaults, with the
  programme set (and E's LF / HF at the Altiverb family's x1.49 / x0.69, the nearest steps: 1.5 / 0.7), but leave the
  preferences (`power`, `link`, `transformer`, `trueStereo`) as they are. Every value is a step of its programme's
  table (`tests/StateTests.cpp`).
- Preset saving: SAVE always opens the Save As dialog prefilled with the current name. Factory names are rejected
  with a message; a chosen name that matches an existing user preset asks Yes/No before overwriting. No returns to the
  name dialog. A modified preset's name is shown italic with a trailing `*`.
- The transformer cores are a guess (H): flux = a leaky integral of the signal (5 Hz corner, normalised at 40 Hz),
  distortion voltage = a·φ|φ| + b·φ³ subtracted from the signal (the nonlinear magnetising current across the source
  resistance). a and b are set for 0.08 % THD at 40 Hz and a core peak of 1 (= +12 dBu), three quarters from the a
  term; the input core runs at half the drive. The stage maps the A/D's full scale (which follows CAL) to that peak
  (`Transformer::setFullScale`), so the cores are driven as in a unit at the factory alignment (nominal +6 dBu, A/D
  clip at +12 dBu) whatever CAL is. A studio running the unit at another nominal would drive the input transformer,
  which sits before the level trim, by that difference. The datasheet points it was checked against are in `Transformer.h`. The flux
  term is low-pass by construction, so no oversampling is needed. `tests/TransformerTests.cpp` holds the calibration.


## Source layout and behaviour notes

(Moved from the plugin README when it became the user-facing page.)

### Layout

| Path | What it is |
|---|---|
| `src/PluginProcessor.*` | Parameters, channel handling (mono/mono, mono/stereo, stereo/stereo), power (bypass), pan, mix, gains, input level for the meter |
| `src/dsp/ReverbEngine.h` | The reverb engine: device bridge, the five programme cores, crossfade, DELAY, LOW PASS (`docs/ENGINE.md`) |
| `src/dsp/engine/`, `src/dsp/programs/`, `src/dsp/generated/` | The bridge and analogue chain, the cores, and their generated tables |
| `src/dsp/AdcCalibration.h` | The A/D's full scale for a CAL setting |
| `src/dsp/StereoWidth.h` | The wet stereo width: level-matched M/S, bass mono above 100 |
| `src/dsp/Transformer.h`, `CaptureEq.h` | The optional TRANSFORMER stage: low-band EQ and guessed transformer distortion (see below) |
| `src/params/` | Parameter IDs (permanent), the joystick step tables, the fader scales |
| `src/presets/` | Preset manager: factory presets (one per programme) and user presets as `.r246preset` XML |
| `src/PluginEditor.*`, `src/ui/` | The panel; `ui/Layout.h` is generated, `ui/Design.h` holds measured constants, `ui/HelpText.h` the tooltip texts (program- and layout-aware) |
| `ui-layout.csv` | The cleaned-up positions from `../ui/Asset positioning.csv` (centres, design space 3495 x 1530) |
| `tools/build_assets.sh` | Resamples `../ui/*.png` (never modified) to `assets/images/` and regenerates `src/ui/Layout.h` |
| `installer/` | macOS `.pkg`, Windows NSIS `.exe`, Linux `.deb` builders |
| `../.github/workflows/plugin-*.yml` | CI and release: manual-trigger only for now |
| `docs/DEVELOPMENT.md` | Rules: cross-platform, assets, CI, parameters |

### Channel modes

The plugin accepts mono in / mono out, mono in / stereo out and stereo in / stereo out. The unit's own algorithm is mono
in / stereo out, so:

- **Mono out** folds the two read-outs down as `(L + R) / sqrt(2)` after each channel's output stage (the transformer
  stage runs on both). In the Tomkin library L and R are uncorrelated (|rho| < 0.05 in every window from 0 to 1.5 s,
  S/M/L/N) and carry the same energy to within about 1 dB, so either side alone is as loud as a side of the stereo
  image, but the fold-down is what the unit gives on a mono bus, with both read-outs' reflections, at the same level
  (power sum; the engine's fold-down is within 0.1 dB of the sides' mean in every programme, `EngineTests`). A plain
  `(L + R) / 2` would be 3 dB down. The second read-out costs a few taps, because the loops are shared, so a separate
  mono engine is not worth building.
- **Mono in / stereo out** is the unit's native path.
- **Stereo in / stereo out** defaults to **ORIGINAL**: the input is summed to the unit's mono input at -6 dB (0.5 x L+R,
  so correlated channels keep their level rather than gaining 6 dB), then the mono reverb feeds both outputs. The
  **ORIGINAL / TRUE ST** button in the bottom bar switches to true stereo, which runs one mono engine per input
  channel (left through engine 0 and right through engine 1). The button is greyed out for the mono layouts. Like LINK
  and XFMR it is a routing preference, so presets neither store nor change it.

### TRANSFORMER (optional, on by default)

The service manual offers the unit with optional balanced input and output transformers ("harmonic distortion
... increased slightly", pp. 37 and 46). The XFMR button in the bottom bar (bool parameter `transformer`)
models such a unit on the wet path:

- **Low-band EQ (measured, S):** the colour that CLA's and the Altiverb's IRs share but Tomkin's library and the stock
  unit do not (`../docs/model_s_results.md` §27): about +3 dB at 40–60 Hz, −0.6 dB at 160 Hz, ripple up to 400 Hz,
  minimum phase, on both outputs. The transformer option is a candidate explanation (H), not a demonstrated one: with the stock coupling capacitor a
  line transformer's resonance would sit at 2–10 Hz, not 40–60 Hz.
- **Transformer distortion (guessed, H):** a flux model of an input transformer (mono, before the reverb) and an output
  transformer per channel. Odd harmonics, mostly 3rd, growing with level and falling steeply with frequency: at
  0 dBFS (= +12 dBu, the unit's maximum) 0.35 % at 20 Hz, 0.08 % at 40 Hz, 0.01 % at 100 Hz, 0.0001 % at 1 kHz;
  calibrated against Jensen and Lundahl line-transformer datasheets, not measured on a 246 (IRs can't show it).

Off, the output is exactly the stock unit's, faithful to Tomkin and the schematic. It is on by default, and a preference
like LINK, so presets neither store nor change it.

## Open items

- [x] Late-tail decorrelator on an early / late split: closed 31 October 2026, not needed (see `stereoWidth` above:
  the late tail is already uncorrelated; the widener alone keeps the early reflections unsmeared).
