# REV-246 (the plugin)

The VST3 / AU plugin that will carry the reverb model worked out in `../docs` and `../tools`. This folder is the
scaffold: the panel, parameters, presets, installers and (switched-off) CI are in place; the reverb itself is a
placeholder (`src/dsp/ReverbEngine.h` passes the mono input to both outputs).

## Quick start

Build the plugin with CMake (JUCE is the submodule `libs/JUCE`):

    git submodule update --init --depth 1 plugin/libs/JUCE
    cmake -S plugin -B plugin/build -DCMAKE_BUILD_TYPE=Release -DRV_STANDALONE=ON
    cmake --build plugin/build --config Release --parallel 3

The Standalone app is `plugin/build/REV246_artefacts/Release/Standalone/REV-246.app` (macOS); the VST3 and AU
bundles sit beside it. Add `-DRV_BUILD_TESTS=ON` and run `ctest --test-dir plugin/build -C Release` for the tests.

## Layout

| Path | What it is |
|---|---|
| `src/PluginProcessor.*` | Parameters, mono-in / stereo-out path, power (bypass), pan, mix, gains, input level for the meter |
| `src/dsp/ReverbEngine.h` | The model's interface (placeholder body) |
| `src/params/` | Parameter IDs (permanent), the joystick step tables (placeholders), the fader scales |
| `src/presets/` | Preset manager: factory presets (one per programme) and user presets as `.r246preset` XML |
| `src/PluginEditor.*`, `src/ui/` | The panel; `ui/Layout.h` is generated, `ui/Design.h` holds measured constants |
| `ui-layout.csv` | The cleaned-up positions from `../ui/Asset positioning.csv` (centres, design space 3495 x 1530) |
| `tools/build_assets.sh` | Resamples `../ui/*.png` (never modified) to `assets/images/` and regenerates `src/ui/Layout.h` |
| `installer/` | macOS `.pkg`, Windows NSIS `.exe`, Linux `.deb` builders |
| `../.github/workflows/plugin-*.yml` | CI and release: manual-trigger only for now |
| `docs/DEVELOPMENT.md` | Rules: cross-platform, assets, CI, parameters |
