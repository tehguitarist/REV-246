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
- `captureEq` (default off) switches on `src/dsp/CaptureEq.h`, the 10-number cascade of the "CLA" entry in
  `../analysis/unit_transfer/bq.json`: a first-order high-pass, an RBJ low shelf and two RBJ peaking sections, with
  the formulas of `../tools/unit_transfer.py` (`hp1_sos`, `lowshelf_sos`, `peak_sos`). It was fitted at the device
  rate; the plugin designs the same sections at the host rate from the same parameters, which matches `bq_resp`
  within 0.001 dB at 50 Hz–2 kHz from 22.05 to 192 kHz (no prewarping needed; the fit is undefined below ~50 Hz).
  `tests/CaptureEqTests.cpp` checks the coefficients and the processed impulse response against `bq_resp` values;
  if `bq.json` is refitted, update `CaptureEq::cla` and the test's reference table together. The stage always runs
  (warm state, a 20 ms crossfade when toggled); off, its weight is exactly 0 and the wet signal is unchanged.
  Presets skip it, like `power` and `link`.
