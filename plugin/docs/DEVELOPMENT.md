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
- IDs in `src/params/Parameters.h` are permanent. The joystick step tables are placeholders (H) until the unit's exact
  step tables are known; only step counts and ranges are from the data sheet (`../docs/EMT246_unit_facts.md`).
- The buses allow mono in / mono out, mono in / stereo out and stereo in / stereo out. The reverb engine is mono in /
  stereo out, so a mono output takes its left side. For stereo in / stereo out, `trueStereo` (default off, the
  ORIGINAL / TRUE ST button) selects between the original summed path (0.5 x L+R into one engine) and two engines,
  one per input channel. Presets skip it, like `power`, `link` and `transformer`.
- `transformer` (default on, the XFMR button) switches on `src/dsp/Transformer.h`: an input transformer core
  on the mono input, then per output the low-band EQ (`src/dsp/CaptureEq.h`) and an output transformer core. Presets
  skip it, like `power` and `link`. The stage always runs (warm state, a 20 ms crossfade when toggled); off, its
  weight is exactly 0 and the signal is unchanged.
- The EQ is the 10-number cascade of the "CLA" entry in `../analysis/unit_transfer/bq.json`: a first-order
  high-pass, an RBJ low shelf and two RBJ peaking sections, with the formulas of `../tools/unit_transfer.py`
  (`hp1_sos`, `lowshelf_sos`, `peak_sos`). It was fitted at the device rate; the plugin designs the same sections at
  the host rate from the same parameters, which matches `bq_resp` within 0.001 dB at 50 Hz–2 kHz from 22.05 to
  192 kHz (no prewarping needed; the fit is undefined below ~50 Hz). If `bq.json` is refitted, update
  `CaptureEq::cla` and the reference table in `tests/TransformerTests.cpp` together.
- `stereoWidth` (default 100, the WIDTH halo) is a sound setting and part of presets. `src/dsp/StereoWidth.h` holds
  the level-matched M/S `StereoWidener` (100 is a bit-exact bypass, 0 mono, 200 wide, bass mono above 100) and the
  `LateDecorrelator` for the late tail. The widener runs on the whole wet (pure M/S gains do not smear transients); the
  decorrelator is not in the chain yet because the placeholder engine has no early/late split. See the open item below.
- `sixteenBit` (default on, the FLOAT / 16-BIT button) is a sound setting and part of presets. It is passed to the
  engine as `Settings::sixteenBit`; the placeholder has no loops yet, so the real algorithm applies the 16-bit storage
  quantisation when it lands. It also switches on the A/D's clipping (`src/dsp/AdcClipper.h`, one per engine, after
  the input transformer; always running, 20 ms crossfade). The unit's converter is a SAR, which sticks at the end of
  its range: a hard clip of each held sample at the device rate, on the input filter's output (+5.1 dB at 10-11.5 kHz,
  so HF clips early). The stage runs the input filter beside the signal (matched-z from the transcribed netlist's
  poles and zeros, one extra zero tuned per rate: within 0.17 dB to 12.6 kHz), samples it at the device instants
  (Kaiser sinc), clips, and reinjects each error through g0 = IFT{taper / H}, the input filter's band-limited inverse.
  Its output is the input delayed by `latency()` (about 1.4 ms) plus those corrections (times 0 in float mode), so
  below full scale it is an exact delay. That delay is the plugin's reported latency, the same in both modes (so
  16-BIT never changes it and its crossfade has no comb); the dry and power-off signals are delayed to match. When the
  real engine lands, add its own latency (resampling) to the report. **The engine must apply the input filter after this stage**
  (as part of its analogue chain); once the engine runs that chain at the device rate, clip there and drop the
  reinjection. What survives the scheme: the clip error to the taper's floor (-35 dB of it for a sine; error above
  11.5 kHz is lost). Cost 0.15 % of a core idle, 0.8 % clipping hard at 48 kHz (2 % at 192 kHz). Design and checks:
  `../tools/analog/adc_clip.py`; tests: `tests/AdcClipperTests.cpp`.
- The A/D's full scale follows the meter calibration (`Rev246Processor::setMeterCal`, the CAL button; UI state, kept
  with the session, not in presets): the calibration stands for the unit's input alignment, and its full scale is a
  sine 6 dB over 0 VU (unit facts §2). OVL lights at +6 VU, on a full-scale input, or when the A/D clips (in float
  mode too: where it would clip), as the unit's REG LED did.
- Display fields other than the program display open a typed-entry dialog on double-click. `params::joyStepFromValue`
  validates the text for the current program and returns a step; `joyInputHint` supplies the unit/range/list shown in
  the dialog. Discrete controls require an exact step, while the geometric RT (S/M/L/E) and HF (N) accept any value in
  range and snap to the nearest step. RT in N is rejected with an information dialog. The same `setJoyStep` path is
  used, so the levers, funnels and displays refresh.
- Initial/factory defaults: S, the CLA reference levers (LP 12 kHz, DELAY 10 ms, AMPL 80 %, LF x1.0, RT 1.0 s, HF
  x0.9), MIX 100 %, INPUT/OUTPUT 0 dB, LINK on, 16-bit on, XFMR on, WIDTH 100, ORIGINAL stereo mode and power on.
  Factory presets reset the stored settings to these defaults but leave the preferences (`power`, `link`, `transformer`,
  `trueStereo`) as they are.
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
| `src/dsp/ReverbEngine.h` | The model's interface (placeholder body) |
| `src/dsp/AdcClipper.h` | The A/D's clipping in 16-bit mode (see `sixteenBit` above) |
| `src/dsp/StereoWidth.h` | The wet stereo width: level-matched M/S, bass mono above 100, and the late-tail decorrelator |
| `src/dsp/Transformer.h`, `CaptureEq.h` | The optional TRANSFORMER stage: low-band EQ and guessed transformer distortion (see below) |
| `src/params/` | Parameter IDs (permanent), the joystick step tables (placeholders), the fader scales |
| `src/presets/` | Preset manager: factory presets (one per programme) and user presets as `.r246preset` XML |
| `src/PluginEditor.*`, `src/ui/` | The panel; `ui/Layout.h` is generated, `ui/Design.h` holds measured constants |
| `ui-layout.csv` | The cleaned-up positions from `../ui/Asset positioning.csv` (centres, design space 3495 x 1530) |
| `tools/build_assets.sh` | Resamples `../ui/*.png` (never modified) to `assets/images/` and regenerates `src/ui/Layout.h` |
| `installer/` | macOS `.pkg`, Windows NSIS `.exe`, Linux `.deb` builders |
| `../.github/workflows/plugin-*.yml` | CI and release: manual-trigger only for now |
| `docs/DEVELOPMENT.md` | Rules: cross-platform, assets, CI, parameters |

### Channel modes

The plugin accepts mono in / mono out, mono in / stereo out and stereo in / stereo out. The unit's own algorithm is mono
in / stereo out, so:

- **Mono out** uses the left side of the reverb (one read-out is enough; the real output taps can revisit this).
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

- [ ] When `ReverbEngine` exposes separate early and late outputs, route `LateDecorrelator` on the late tail only,
  before it is summed with the early reflections, then run `StereoWidener` on the sum. Check that the early reflections
  are not smeared, the decorrelator's ~8 ms delay is confined to the tail, and the width sweep through 100 is
  click-free. If the engine cannot expose the split, keep the widener only and reconsider, preserving the design
  intent: width changes must not smear the early reflections.
