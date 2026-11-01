# REV-246 (plugin repository)

A JUCE VST3/AU/standalone reverb plugin in `plugin/`; the repository root holds the README, licence and the listening
page (`docs/listen`, served by GitHub Pages). The measurement and analysis work behind the reverb is in the history,
at the tag `research-archive`.

Read `plugin/docs/DEVELOPMENT.md` (rules, parameters, source layout) and `plugin/docs/ENGINE.md` (the reverb engine) first.

## Rules

- **All C++ must be cross-platform** (macOS, Windows/MSVC, Linux/GCC): no `M_PI` (use `juce::MathConstants`), no
  POSIX-only calls, `juce::File` for paths, `juce::exactlyEqual` for float equality.
- **Workflows are manual-only** (`workflow_dispatch`): CI first (`plugin-ci.yml`), then `plugin-release.yml`, which refuses a
  commit whose CI has not passed. Every GitHub Action must run on the current Node.
- The version is the `project(REV246 VERSION x.y.z)` line of `plugin/CMakeLists.txt` (release.yml greps it; keep the
  format). Bump it when an Audio Unit changes, or Logic keeps the cached one.
- Don't commit `plugin/build*/`. `ui/` holds the artwork masters and is never modified; `plugin/assets` is built from it.
- Parameter IDs are permanent. After the first release a step table's range must not change (see DEVELOPMENT.md).
- Generated headers (`plugin/src/dsp/generated/`, `plugin/src/presets/FactoryPresets.h`) came from generators that
  now live only in the history; edit them by hand.
- Commit messages end with the attribution line from the session's system reminder.
