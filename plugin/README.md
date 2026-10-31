# REV-246

![Build](https://github.com/tehguitarist/REV-246/actions/workflows/plugin-ci.yml/badge.svg?branch=main)
[![License](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://opensource.org/license/agpl-v3)
[![Downloads](https://img.shields.io/github/downloads/tehguitarist/REV-246/total)](https://somsubhra.github.io/github-release-stats/?username=tehguitarist&repository=REV-246&page=1&per_page=30)

![REV-246](docs/images/screenshot.png)

**A free, open-source digital reverb with a vintage heart.** REV-246 is inspired by the EMT 246, the rare late-1970s
studio reverb prized for its dense, smooth tails and its distinctly period character: a little grit, a little weight,
and a sound that sits in a mix instead of on top of it. It runs as AU, VST3 and standalone, on macOS, Windows and Linux.

Five programmes cover everything from a tight room to a twenty-second cathedral, including a non-linear programme for
gated and reverse drum sounds and the older EMT 250 algorithm. The panel keeps the original's six levers, glowing readouts and
fader pack. Beneath it are modern conveniences that the original made hard or impossible: an extra pre-delay, a ducker,
a stereo-width control, true-stereo operation, presets, and switches for the 16-bit and transformer character that
were fixed in the hardware (or, for the transformers, a factory option you had to order).

> REV-246 is an independent project. It is not affiliated with or endorsed by EMT or any other manufacturer.

**Status:** pre-release. The reverb engine, panel, presets and installers are in place and are being checked and
listened to before the first release, so there is no download yet. When there is, the latest one will be here:
**[Download the latest release](https://github.com/tehguitarist/REV-246/releases/latest)**.

## Highlights

- **Five programmes**, each with the full set of levers.
- **Period character, switchable.** 16-bit processing and an optional transformer stage add the grit and weight of the
  original; switch them off for a clean, modern tail.
- **Six levers:** low-pass, reflections delay and amplitude, and decay: low-frequency multiplier,
  reverberation time and high-frequency multiplier.
- **Modern conveniences the original lacked:** pre-delay (0 to 200 ms), an input-keyed ducker, a WIDTH control that
  keeps the bass mono as it widens, a mono-sum or true-stereo input mode, and one-click switches for the 16-bit and
  transformer character.
- **64 factory presets** from small rooms and chambers to halls and cathedrals, plus your own user presets.
- **Light on the CPU:** about 0.5 to 0.8 % of one Apple-silicon core at 44.1 / 48 kHz and about 1 % at 96 kHz, with about
  1 ms of latency that your DAW compensates for.

## Programmes

| Programme | Character | Reverberation time |
|---|---|---|
| **S** | Short reverb: rooms, booths, live spaces (best between 0.2 and 1 s) | 0.2 to 2.5 s |
| **M** | Medium reverb: the all-round concert hall (best between 0.8 and 4 s) | 0.5 to 5 s |
| **L** | Long reverb: the cathedral, with widely spaced early reflections (best between 3 and 10 s) | 2 to 20 s |
| **E** | The EMT 250 algorithm: recirculating delay loops with slowly moving taps, with extended pre-delay | 0.4 to 4.5 s |
| **N** | Non-linear: an early section of fixed length (SPAN) with its own decay, then the late reverb; short SPAN gives gated reverb and a negative EARLY time gives reverse reverb | SPAN, early and late times |

Every programme has the same six levers. The readouts above them show the value for the programme you are in; in N they
change to the hardware's own labels (DELAY, SPAN, EARLY and LATE).

## Features

The first group is the original's; the second is what REV-246 adds.

**The original's:**

- **Six joystick levers** for the reverb's parameters, stepping through the same kinds of values as the hardware.
  Click, hold, drag, use the mouse wheel or double-click to type a value.
- **Mix, input and output faders** with a PAN control for the dry signal, and an optional LINK to move input and output
  together.
- **Level meter** with an overload light and a calibration button (**CAL**) for the reference level.
- **16-bit mode:** reproduces the original's 16-bit reverb storage, for its faint, signal-dependent noise floor and grit in
  the decays, and clips the input where the original's A/D converter did. Set it to FLOAT for full precision.
- **XFMR:** the weight in the low end and gentle low-frequency harmonic distortion of the unit's optional input and
  output transformers. Off, you hear the standard unit.

**Added for modern use:**

- **Pre-delay:** up to 200 ms before the reverb, on top of the programme's own reflections delay.
- **Ducker:** lowers the reverb while the input is over a threshold, so it blooms in the gaps. Range, threshold,
  attack and release are adjustable.
- **WIDTH:** sets how wide the wet signal sits, from mono through the original's width to extra wide, with the bass kept
  mono as it widens.
- **16-BIT / FLOAT and XFMR switches:** the original's storage depth was fixed and its transformers were an ordered
  option; here both are one click, so you can A/B the period character against a clean tail.
- **Factory and user presets**, with a bypass (power) button.
- **Mono in / mono out, mono in / stereo out and stereo in / stereo out.** With a stereo input the reverb can sum to
  mono first (**ORIGINAL**, the default) or run a reverb per channel (**TRUE ST**).

## Install

Download the installer for your system from the [latest release](https://github.com/tehguitarist/REV-246/releases/latest):

| System | Installer | Formats |
|---|---|---|
| macOS | `.pkg` | AU, VST3 |
| Windows | `.exe` | VST3 |
| Linux | `.deb` | VST3 |

Then rescan plugins in your DAW. Logic caches Audio Units, so if REV-246 does not show up, restart Logic or run
`auval -v aufx Rv24 Lprc`.

## Using it

The initial state is fully wet, with INPUT and OUTPUT at 0 dB and LINK on, a short programme with a one-second decay, WIDTH at
100 and both character options on.

**Pick a programme** (S, M, L, E or N) and set the six joysticks. Each joystick is the original's lever: press the
upper half of the round button to step the value up and the lower half to step it down. Hold it to repeat; after a
second it speeds up. For convenience you can also drag or click the wedge above it, use the mouse wheel for one step per
notch, and Option/Alt-click to return to the default.

| Control | What it does |
|---|---|
| **LOWPASS** (f, kHz) | Low-pass filter on the reverb, 2 to 12 kHz |
| **REFLECTIONS** (DELAY, AMPL.) | Delay before the early reflections, 0 to 999 ms, and their level, 0 to 100 % |
| **DECAY** (LF ×, RT, HF ×) | Reverberation time in seconds, with multipliers for how much faster or slower the low and high frequencies decay |
| **WIDTH** | Wet stereo width: 0 = mono, 100 = unchanged, 200 = extra wide. Above 100 the bass stays mono |
| **PRE-DELAY** | Extra delay before the reverb, 0 to 200 ms (after the ducker). Option/Alt-click for 0 |
| **DUCK** | Lowers the reverb while the input is over a threshold. Click for its settings: range (0 to -60 dB), threshold, attack and release |
| **MIX** | Wet/dry balance. PAN places the dry signal |
| **INPUT / OUTPUT** | Levels, with LINK to move them together |

In the non-linear programme (N) the readouts change meaning to match the hardware's own labels (DELAY, SPAN, and early
and late RT). Hover over any control for a description.

## Width

The WIDTH halo is the first control in the right-hand group. Its arc runs from 0 to 200 with the value in the middle;
100 leaves the wet signal unchanged, 0 sums it to mono and 200 widens it. Above 100 the low band is kept mono so the
bass does not wander. Drag up to widen and down to narrow (Shift for fine), or use the wheel; double-click to type a
value, or Option/Alt-click to return to 100. Width only affects the wet path, never the dry signal.

## Typed values

Double-click any display field except the program display to type a value instead of using its lever. The dialog lists
the valid values or range for the current program and unit; an invalid entry leaves the current value untouched. S/M/L/E
and N have different meanings and ranges (for example amplitude in % vs span in s), and RT is not used in N.

The bottom bar holds the preset browser (previous, next, name, **SAVE**) and **CAL** for the meter reference on the
left, then the **WIDTH** halo and the routing and character buttons (**ORIGINAL / TRUE ST**, **FLOAT / 16-BIT**,
**XFMR**), **LINK** and **?** on the right. A changed preset's name turns italic with a `*`. Saving under an existing user-preset name asks
before replacing it; factory names are read-only. ORIGINAL / TRUE ST, XFMR and LINK are preferences, so presets
neither store nor change them.

## Character options

16-bit and XFMR are independent, both default on, and both only change the reverb path.

**16-bit.** The original stored the intermediate results in its reverb loops as 16-bit words, and truncating them at
that depth costs a little noise and grit: its service data put the reverb around 15 dB noisier than the unit's own
straight-through path. With **16-BIT** selected, REV-246 quantises its loops the same way, and adds a little analogue
hiss: very quiet, below the original's own specification, and only in 16-BIT mode. **16-BIT** is the default; switch
to **FLOAT** for a cleaner tail.

In 16-bit mode the input also clips where the original's A/D converter did. That converter left only 6 dB of headroom
above the 0 mark on its meter, and it clipped hard: past full scale it simply holds the last code. Its input filter
lifts the treble by up to 5 dB first, so cymbals and sibilance clip sooner than the meter suggests. And because it
clipped at its own sample rate of 26.4 kHz, the clipping comes out as inharmonic grit rather than plain harmonics.
Where it clips follows **CAL**: CAL sets what 0 VU means in dBFS, and the converter clips 6 dB above that for a low
sine (at the default -18 dBFS, that's a sine peaking at -9 dBFS). **OVL** lights whenever the converter clips, and in
FLOAT mode it lights where it would have clipped.

REV-246 reports about 1 ms of latency (0.8-1.1 ms: the modelled output filter's correction and the conversion to and
from the original's own sample rate), and your DAW compensates for it. It is the same in every programme and in FLOAT and 16-BIT, so switching
never shifts the timing.

**XFMR.** The original could be ordered with balanced input and output transformers, which the manual says raise the
harmonic distortion slightly. **XFMR** is on by default; with it on, the wet path gains the low-frequency colour and
gentle odd-harmonic distortion of that option. The low-frequency colour comes from measurement; the distortion is an informed
estimate, since measurements of a unit's response cannot show it.

## Hear it

**[Listen in your browser](https://tehguitarist.github.io/REV-246/listen/)**: looping players that switch between dry, IR and REV-246 in sync.

Four recordings, each with the **dry** take, a convolution with a reference **impulse response (IR)** and **REV-246** set to the
same settings. All at a **35 % wet mix** (0.65 dry + 0.35 wet); the two wets are level-matched, so what you hear is the reverb's
character and not its level. The recordings are summed to mono for the reverb (the original is mono in), the plugin's
~1 ms latency is removed as a DAW would, and the three files of an example share one gain. The IR sources are numbered
1 to 3 and stay unnamed.

| Example | Settings | Dry | IR source | REV-246 |
|---|---|---|---|---|
| Room mic | Programme S, RT 1.0 s, LF 1.0, HF 0.9, AMPL 80 %, DELAY 10 ms | [dry](../docs/listen/audio/room-mic_dry.mp3) | 1: [IR](../docs/listen/audio/room-mic_ir.mp3) | [REV-246](../docs/listen/audio/room-mic_plugin.mp3) |
| Piano | Programme M, RT 1.5 s, LF 1.2, HF 0.7, AMPL 60 %, DELAY 40 ms (factory preset *Altiverb P1 Piano*) | [dry](../docs/listen/audio/piano_dry.mp3) | 2: [IR](../docs/listen/audio/piano_ir.mp3) | [REV-246](../docs/listen/audio/piano_plugin.mp3) |
| Acoustic guitar | Programme S, RT 1.4 s, LF 1.0, HF 0.6, AMPL 60 %, DELAY 0 | [dry](../docs/listen/audio/acoustic-guitar_dry.mp3) | 3: [IR](../docs/listen/audio/acoustic-guitar_ir.mp3) | [REV-246](../docs/listen/audio/acoustic-guitar_plugin.mp3) |
| Snare (gated) | Programme N, SPAN 0.2 s, EARLY RT -3.0 s, LATE RT 0.8 s, DELAY 39 ms (factory preset *Altiverb P54 Gated*) | [dry](../docs/listen/audio/snare_dry.mp3) | 2: [IR](../docs/listen/audio/snare-ir2_ir.mp3) | [REV-246](../docs/listen/audio/snare-ir2_plugin.mp3) |
| Snare (gated) | The same, DELAY 0 | [dry](../docs/listen/audio/snare_dry.mp3) | 3: [IR](../docs/listen/audio/snare-ir3_ir.mp3) | [REV-246](../docs/listen/audio/snare-ir3_plugin.mp3) |

The reverb is the plugin's own engine set to FLOAT with XFMR off: an impulse response is linear and noise-free, so 16-BIT's noise and grit and XFMR's added low-end weight would only move the plugin away from it. Both are there to switch on. The IRs are single measurements of a
unit, so a time-varying programme (M here) and the original's analogue chain will differ from them in detail: see
[analysis/ab_samples/report.md](../analysis/ab_samples/report.md) for how close the wets are, band by band.

## Build

You need CMake 3.15 or newer and a C++17 compiler (Xcode on macOS, Visual Studio 2022 on Windows, GCC or Clang on
Linux). JUCE is a git submodule.

```bash
git clone --recurse-submodules https://github.com/tehguitarist/REV-246.git
cd REV-246
cmake -S plugin -B plugin/build -DCMAKE_BUILD_TYPE=Release -DRV_STANDALONE=ON
cmake --build plugin/build --config Release --parallel 3
```

The plugins land in `plugin/build/REV246_artefacts/Release/{AU,VST3,Standalone}` (AU on macOS only). Add
`-DRV_BUILD_TESTS=ON` and run `ctest --test-dir plugin/build -C Release` to build and run the tests.

Development rules (cross-platform code, assets, CI) and the source layout are in
[docs/DEVELOPMENT.md](docs/DEVELOPMENT.md). The plugin is one half of this repository: the other half, in `docs/` and
`tools/` at the top level, is the research behind the reverb.

## Licence

REV-246 is free software under the [GNU AGPLv3](../LICENSE). It uses [JUCE](https://juce.com) under its AGPLv3 option,
and the notices for it and the other third-party code are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Credits

Built by Leigh Pierce, using [JUCE](https://juce.com/) for the plugin framework and interface.
