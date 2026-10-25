# UAD EMT 250 plugin: headless capture and test plan

*Status (25 October 2026): harness built and smoke-tested; first full capture run in `analysis/uad250/captures/` (gitignored). No analysis of the captures yet except the headline checks in §2.*

**Why.** `docs/EMT250_program.md` §4 and `NEXT_STEPS.md` §10: the Altiverb REV-250 presets are single frozen captures, so the 250 program's modulation, level behaviour and parameter-change behaviour cannot be read from them. The UAD EMT 250 is a software emulation we can drive with any stimulus and any parameter state, so it is a **test bench** for (a) detectors (known-modulated reverb), (b) a *prior* for what the E program's modulation and non-linearities look like, and (c) a broad cross-check of the static structure (loop lengths, early taps, band RTs). It is **not** the emulation target: the target is still the 246's Program E (`EMT250_program.md` scope). Claim labels below follow the project convention (C / S / H).

## 1. How the capture works (no DAW)

`tools/uad250/au_render.swift` hosts the AU **in-process** with `AVAudioEngine` in offline manual-rendering mode: a source node feeds the file, the plugin renders, the result is written as 32-bit float WAV. Parameters are set by name (`--set Time="2.5 s"`) and can change **mid-render on a sample boundary** (`--auto events.json`: `{t, param, value}` or `{t, action: "reset"}`); blocks are split at event times.

What was learned getting it to work (all **C**, measured here):
- **Out-of-process hosting (`.loadOutOfProcess`) does not work**: the sandboxed host cannot reach the UAD driver and the plugin returns only a constant −98 dBFS noise floor, ignoring input and parameters. In-process works.
- Calling the AU render block directly fails with −50 in-process; the engine route works. The plugin segfaults in its atexit handlers, so the host ends with `_exit(0)` after the file is closed.
- It runs the **UAD-2 DSP path** (log: "UAD2 mode", Apollo Twin X / UAD-2 Satellite present). Reported latency is **523 samples at 48 kHz**; every IR carries it (measure it per rate, do not assume it).
- **Renders are faster than real time** (a 30 s file takes about a second), so the full suite is minutes, not hours.
- **Bit-identical repeatability across fresh instances**, with modulation both on and off (null −300 dB over 0–4 s). The plugin's randomness is a seeded generator that restarts the same way each load.
- **WetSolo = 1 forces 100 % wet whatever Mix says.** Mix-law and dry-path cases therefore set WetSolo = 0.
- Parameters (15): Program (Reverb / Delay / Phase / Chorus / Echo / Space), Time (16 steps 0.4–4.5 s), LFDecay (×0.5 / 1 / 1.5 / 2), HFDecay (×0.25 / 0.33 / 0.50 / max), LDelayTime and RDelayTime (0–315 ms in 5 ms steps), Phase (16 steps to 15 ms), ChorusMode (I–IV), Echo (0–315 ms), Predelay (0 / 20 / 40 / 60 ms), Outputs (Rear / Front), Mix, WetSolo, **Noise (on/off)**, Power. **Noise is the modulation switch**: Noise 1 versus 0 differs by −34 dB (0–150 ms) growing to 0 dB (fully decorrelated) by 2–4 s at Time 2.0.

Run it: `.venv/bin/python tools/uad250/uad250_suite.py build | list | run [--only 'C_mod/*'] [--smoke]`. One process at a time under `tools/memguard.py`; finished cases are skipped, so a run can be stopped and resumed.

## 2. Headline checks already done (smoke / bring-up)

| Check | Result |
|---|---|
| Same file, two fresh instances, Noise 1 | identical, −300 dB (0–4 s) **[C]** |
| Same, Noise 0 | identical, −300 dB **[C]** |
| Noise 1 vs Noise 0, Time 2.0 | −34.5 / −25.2 / −18.4 / 0.0 dB at 0–0.15 / 0.15–0.5 / 0.5–2 / 2–4 s **[C]** |
| First arrival, 48 kHz, Time 1.0, predelay 0 | R 13.0 ms, L 24.9 ms (R leads L by ≈ 11.9 ms, as in the Altiverb REV-250) **[C]** |

## 3. The tests

Groups, case counts and audio length from `uad250_suite.py list` (447 cases, about 107 min of audio). Inputs are mono, fed to L and R identically unless stated.

### A. Static structure, Noise off (260 cases): the algorithm, broadly
- **A1** an impulse at every one of the 16 Time × 4 LF × 4 HF settings (256). Gives per-band T60 (RT, LF/HF multipliers), the loop lengths by the harness's decay-flattened ACF, the early taps, the L/R onset offset, and how each loop's feedback gain and loop filter move with the controls. Directly comparable with `EMT250_program.md` §1–§3 (167.98 / 190.12 / 210.52 ms loops, 3 early taps at 0 / 16.31 / 27.27 ms, R before L by 11.92 ms) and with the 246 E model fit.
- **A2** log sweeps (ESS, 20 Hz–20 kHz, −12 dBFS) at four Times: the same IR with ~40 dB more SNR than an impulse, plus a linearity cross-check against A1.
- Questions it answers: do the plugin's three loops match 4430 / 5014 / 5552 device samples, or the real 250's (≈ 24 kHz)? Is the early-tap gain fixed in RT (Altiverb says yes; the patent says no)? What is the LF/HF decay law versus the 246 S/M/L law? Is the loop count/length RT-dependent?

### B. Non-linear and level behaviour (62 cases): *the main target, part 1*
- **B1 impulse ladder**, +12 … −90 dBFS in 6 dB steps at Time 0.4 and 4.5, Noise off. Normalise each IR by its input gain and compare: any residual is non-linearity. Reads input clipping/saturation (the +12 and +6 dBFS cases), internal overflow, **quantisation / truncation noise in the loops** (the tail floor versus level), and a **dead zone or limit cycles** in very quiet tails (does the decay reach the floor or stall or hold a pattern). The 246's own 16-bit word length is the question this informs.
- **B2 sine ladder**, 250 Hz / 1 kHz / 4 kHz at −60 … +6 dBFS, Time 0.4, Noise off: steady-state harmonic content versus level (THD and its level slope), plus a two-tone (1000 + 1100 Hz) intermodulation case.
- **B3 floor**: silence for 20 s at Noise 0/1 × Power on/off (idle noise, dither, spectrum, any idle tones), and white noise through Power off (bypass path) and Mix 0 (dry path) for the converter/analogue chain's bandwidth and gain.
- Caution: float inputs above 0 dBFS may be clipped by the host/driver rather than the plugin; compare the +12 / +6 / 0 steps before attributing it to the emulation.

### C. Modulation and time variance (35 cases, 39 min): *the main target, part 2*
- **C1 impulse trains**, 24 impulses 12 s apart at Time 1.0 / 2.5 / 4.5, Noise on, plus Noise-off controls and a train with +1 sample drift per impulse. If the generator free-runs, each impulse sees a different modulation state, so the 24 are **24 independent snapshots of the same setting**: the capture-to-capture decorrelation curve per band (the `ml_structure.md` §4–5 metric), the depth against the Noise-off IR, and which part of the IR varies. If the 24 are identical the state restarts on onset, which is itself the finding. The Noise-off train tests **time invariance** including block phase.
- **C2 long tones**, 250 Hz / 1 kHz / 4 kHz, 60 s, Noise on and off (plus 1 kHz at Time 0.4 and 4.5). Complex demodulation of the reverb output gives phase (FM, delay modulation) versus amplitude (AM, gain modulation) sidebands: rate, depth, waveform (random, stepped or smooth), and whether depth grows with frequency as delay modulation must. The Noise-off twin is the control.
- **C3 onset shift**: the same impulse at 0 / 1 / 2 / 3 / 100 / 511 / 512 / 513 / 1000 samples, 0.5 / 1 / 2 / 3.7 s, each in a fresh instance, plus an exact-repeat control. Tells whether the random state is tied to stream time (and at what granularity: 1 sample, one block of 512, or the device clock) or to the instance.
- **C4 reset and toggle**: plugin reset between two impulses (does it re-seed), Noise switched off/on during a tone, Noise switched off mid-tail (does the modulation stop at once, hold, or ramp).
- **C5 repeated sweeps**, six ESS 8 s apart, Noise on and off: the sweep-pair calibration of `NEXT_STEPS.md` §10, for the detectors that run on library sweeps.
- Mapping to the 246: `EMT250_program.md` §8's open question is **which modulation the E uses** (the 250's LFSR address jitter, or the 246's stochastic control with M/L's 0.83 Hz / 0.15 rad rms rotation). The UAD answers the first half: what a modelled 250 looks like (fast jitter versus slow rotation, depth, per-loop or global). It does not settle which one the E uses; that still needs the 246-E AMPL pairs.

### D. Other programs, rate and routing (77 cases)
- **D_programs**: Chorus I–IV and Space (impulse train + 1 kHz tone, Noise on/off), the other 250 time-based programs with randomised delay; Delay at several L/R times, Echo, Phase and Predelay at **96 kHz**. The 5 ms delay step measured in 96 kHz samples against the stated step recovers the plugin's internal rate and whether delays are integer samples of an internal clock (the same method as the 246's 26,373.6 Hz).
- **D_rate**: the same impulse and white noise at 44.1 / 48 / 88.2 / 96 / 192 kHz, Noise on and off. Is the emulation rate-independent (an internal rate plus SRC, with its own filter) or run natively; sets latency, bandwidth, and whether Noise-state timing depends on host rate.
- **D_route**: Mix law (0 … 1 at WetSolo 0), WetSolo override, Outputs Rear/Front, L-only / R-only / both input (the input summing law and channel gain).

### F. Parameter changes while ringing (13 cases)
Automation lands on a sample: Time 4.5 → 0.4 / 0.4 → 4.5 / 4.5 → 2.5 mid-tail, LF and HF decay jumps, Predelay jumps under a 1 kHz tone (dry-only and wet), L delay 0 → 315 ms under a tone, a Mix ramp (zipper noise or smoothing), Power toggled mid-tail, Program switched mid-tail. These reveal how the feedback gains are applied (instant, ramped or recomputed), the **delay-change implementation** (a click, a crossfade or a pitch ramp = the interpolation quality, which is what the 250's address randomisation exercises), and state retention on bypass and program switch.

## 4. Analysis to build next (not done)
1. A1 harness: T60 per octave band, `find_loops` per setting, early taps, L/R offset, grid versus `EMT250_program.md` and the Altiverb REV-250 presets (and a fit of the plugin as a prior for the 246 E).
2. B1 normalised-difference ladder; noise-floor versus level; deep-tail dead-zone check. B2 harmonic ladder.
3. C1 decorrelation curve per band versus onset time (and the repeat-period question). C2 demodulation: the rate, depth and waveform of the modulation. C3/C4 as tables.
4. F: the instantaneous-frequency tracks (predelay, delay) and the gain envelope across each Time/LF/HF jump.

## 5. Limits (state these with any result)
- The UAD plugin is a **commercial emulation of the 250**, not the 246's E and not an oracle. Everything it says about modulation and non-linearity is a *prior* until the 246's own evidence agrees. Keep it labelled **H** for the 246 until then.
- **Licensing:** the captures are renders of a licensed plugin and stay local (`analysis/uad250/captures/` is gitignored); do not commit them or derived assets that could stand in for the plugin.
- The plugin runs on the UAD-2 DSP; its sample-rate handling and latency (523 samples at 48 kHz) are part of what is measured, so align each IR by its measured first arrival, not by assuming zero latency.
