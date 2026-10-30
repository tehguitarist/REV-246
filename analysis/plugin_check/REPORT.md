# Phase 5: the plugin's own output against the Tomkin 96 kHz IRs

*30 October 2026. Tool: `tools/plugin_check.py`. Renders: `plugin/tools/render/RenderIr.cpp` (target `REV246Render`, option `RV_BUILD_RENDER`, built in `plugin/build-render`). All values are medians over captures and are nulls (error energy over capture energy, dB) unless stated. Tracked summaries: `summary_{s,ml,n,e}.json`, `chain.json`, `chain_ratio.npz`, `chain_correction_fir.json`, `loff.json`, `rates.json`, `align.json`, `sim/`. Renders and per-setting files are gitignored.*

## 0. Verdict
- **The plugin fails the bar in every program, and the analogue chain alone is the cause.** The bridge's netlist chain (input filter × output filter × S&H lag × coupling corners) differs from the measured chain (`chain3`):
  - phase: about 0.7° across the mid band, +3 to +13° from 8 to 11 kHz, up to +70° at 12.3 kHz;
  - magnitude: within 0.1 dB up to 6.4 kHz, −1 to +3 dB at 11–12.3 kHz.
- That limits any plugin-vs-capture null to about −35 to −38 dB in the mid band, −25 dB at 6.4–10 kHz and −14 to −3 dB above 10 kHz. The result is −28 to −33 dB broadband (100 Hz–10 kHz) in S, M, L and N, and 0 scored captures at −50 dB.
- **The cores match the benches.** With chain3 swapped in for the bridge's chain:
  - S: −52.6 dB at 50–150 ms (bench −53.0);
  - M / L cluster and gap: −51.4 / −51.8 (bench −53.7 / −52.6);
  - E: 4,588 of 5,358 cells inside (bench 4,601).
- **Host rates agree.** 44.1 and 48 kHz match 96 kHz to −88 dB.

## 1. Method
- **Render.** `REV246Render` runs `rv::dsp::ReverbEngine` on a unit impulse: DELAY 0, LOW PASS 12 kHz (bypass), float mode. Options: modulation off (`mod=0`), device-rate core path (`domain=device`), any mono input (`input=`).
- **Capture path.** Each 96 kHz render is deconvolved with chain3, resampled to the device grid as the benches' cache builders do, and aligned per channel on cluster tap 0 as each bench aligns its captures.
- **Per-capture fit.** Each bench's own fit runs with the plugin in the model's place:

  | Program | Per-capture fit |
  |---|---|
  | S | `q1._clock_job`: clock + gain, shift only for AMPL 0; variant: shift free (still 3 numbers) |
  | M / L | `ml_lines.fit_warp`: shift, clock and 2 gains per channel on cluster + gap + first pass |
  | N | shift, clock, 2 gains; branch offset fixed at +1 |
  | E | e_pulse statistics |

- **Held out.** S is held out by RT step, as in the bench: its library-step tables are the fold constants that left that RT step out. M / L: everything after the first pass is predicted.
- **Noise floor.** S: the cross-take repeat pairs. M / L: the law files' repeat captures. N: the bench's noise floor. Excess = model null − floor.
- **Runs.** All under memguard with at most 4 workers. Peak memory: S 2.25 GB, M / L 5.5 GB, N 1.1 GB.

## 2. S (127 settings, 491 captures)

### Nulls, 100 Hz–10 kHz, R (L in brackets)

| Window (ms) | Plugin, bench fit | Plugin, shift free | Plugin + chain3 | Bench | Repeat | Excess: shift free / + chain3 / bench | ≤ −50 dB: shift free / + chain3 / bench |
|---|---|---|---|---|---|---|---|
| 0–50 | −25.8 | −28.2 (−27.5) | −49.7 (−51.8) | −51.5 (−52.1) | n/a | – | 0 / 185 / 329 of 491 |
| 50–150 | −28.7 | −29.4 (−28.6) | −52.6 (−52.2) | −53.0 (−53.2) | −57.2 | +27.8 / +4.6 / +4.2 | 0 / 467 / 421 |
| 150–300 | −32.8 | −32.9 (−32.1) | −51.7 (−51.0) | −51.4 (−51.2) | −54.5 | +21.7 / +2.8 / +3.2 | 0 / 362 / 334 |
| 300–600 | −32.8 | −35.1 (−34.8) | −48.8 (−48.4) | −46.7 (−46.5) | −47.6 | +12.5 / −1.2 / +0.8 | 0 / 187 / 105 of 484 |
| 600–1000 | −28.4 | −31.4 (−31.3) | −42.6 (−42.4) | −39.6 (−39.4) | −40.0 | +8.6 / −2.5 / +0.4 | 0 / 37 / 21 of 414 |
| 1000–1500 | −25.1 | −29.0 | −37.2 | −35.2 | −34.6 | +5.6 / −2.6 / −0.6 | – |

The S bar fails by 21–25 dB at 0–300 ms. With chain3, S matches the bench and beats it after 300 ms; the shift-free fit accounts for that late gain. The CLA stretch was not run: the chain mismatch rules it out, and the chain3 variant carries Tomkin's capture chain, not CLA's.

### Full band, R, 50–150 ms (0–50 ms in brackets)

| Band (Hz) | Plugin, shift free | Plugin + chain3 | Bench | Bench + sub-100 Hz T | Repeat | Plugin vs bench model |
|---|---|---|---|---|---|---|
| 20–40 | −24.5 (−21.3) | −35.2 (−32.2) | −32.1 (−25.7) | −41.2 (−27.3) | −28.0 | −23.8 |
| 40–80 | −29.4 (−31.2) | −49.6 (−39.4) | −47.2 | −51.6 | −37.8 | −29.2 |
| 80–160 | −33.8 | −55.6 | −51.0 | −51.2 | −49.4 | −33.5 |
| 160–270 | −36.0 | −48.2 | −46.0 | −46.0 | −61.6 | −34.8 |
| 400–800 | −36.6 | −57.6 | −52.4 | −52.4 | −61.4 | −34.4 |
| 1.6–3.2 k | −32.5 | −59.1 | −54.6 | −54.6 | −60.2 | −28.5 |
| 3.2–6.4 k | −29.3 | −53.1 | −56.0 | −56.0 | −58.7 | −24.0 |
| 6.4–10 k | −26.9 (−26.5) | −49.6 (−48.3) | −49.9 (−52.9) | −49.9 | −53.6 | −22.7 |
| 10–11 k | −16.4 (−17.3) | −48.6 (−46.3) | −49.5 (−49.9) | −49.5 | −47.4 | −20.8 |
| 11–12.3 k | −12.9 (−13.5) | −44.2 (−43.2) | −44.0 (−38.1) | −44.0 | −42.5 | −15.5 |

The last column is the plugin against the bench's own model under the same measurement assumptions. It is the chain difference of §3 plus the tap-0 alignment of a chain-distorted tap, not a port error: the cores null their reference IRs at −136 dB.

## 3. Chain3 over the bridge's chain (best gain and delay; R and L within 0.6 dB)

| Band (Hz) | Magnitude (dB) | Phase (°) | Null-equivalent (dB) |
|---|---|---|---|
| 20–40 | −0.06 to +0.12 | −7.5 to −2.9 | −22.0 |
| 40–80 | +0.02 to +0.12 | −2.9 to −1.6 | −28.3 |
| 80–160 | +0.03 to +0.05 | −1.5 to −0.8 | −33.9 |
| 160–800 | +0.01 to +0.06 | −0.9 to −0.5 | −37.2 to −38.7 |
| 0.8–6.4 k | −0.10 to +0.10 | −1.3 to −0.4 | −35.3 to −36.3 |
| 6.4–10 k | −0.22 to +0.08 | −1.3 to +8.8 | −24.9 |
| 10–11 k | −0.04 to +0.29 | +8.8 to +12.8 | −14.1 |
| 11–12.3 k | −1.45 to +3.13 | +12.8 to +69.5 | −3.2 |

**Verdict:** the bridge's chain hurts against Tomkin compared with chain3, by 21–25 dB in S, about 21 dB in the M / L cluster and about 20 dB in N. The mismatch has three parts:
- the filter-skirt phase above 8 kHz;
- a near-constant −0.7° in the mid band;
- below 80 Hz, Tomkin's own capture chain (the ~1.8 Hz corner; H).

A compact device-rate FIR correction does not close the gap:

| FIR length | Above 400 Hz | 6.4–10 kHz | 10–11 kHz | 11–12.3 kHz |
|---|---|---|---|---|
| 9 taps | −39 to −51 dB | −28 | −21 | −14 |
| 33 taps | −40 to −49 dB | −41 | −34 | −28 |
| 65 taps | ≤ −42 dB | −51 | −51 | −47 |

At every length 20–80 Hz stays at −23 / −31 dB, so the low band needs a pole-zero form.

## 4. The two H items (no default changed)

### Coupling corners: 0.72 / 0.16 Hz vs the netlists' 0.617 / 0.157 Hz
An exact LTI swap on every S render, with the bench fit kept. With the netlist corners the error grows in almost every channel-window (2,932 at 50–1000 ms):

| Band (Hz) | Error change, median | Netlist better in |
|---|---|---|
| 20–40 | +0.46 dB (+0.005 to +0.65) | 0 / 2,932 |
| 40–80 | +0.42 dB | 23 / 2,932 |
| 80–160 | +0.35 dB | 111 / 2,932 |
| 160–270 | +0.26 dB | 38 / 2,932 |

The data favour the bridge's corners, weakly and not independently: chain3's high-pass is exactly these corners, and below 80 Hz the capture-side ~1.8 Hz term (3–4° at 30 Hz) outweighs the ≤ 0.3° between the two sets. **Keep 0.72 / 0.16.**

### L's output offset: 0.5 device samples vs the harness's 19.19 µs (0.506)
One fractional shift per channel, plugin onto file at 96 kHz, over cluster tap 0 −1 ms to +60 ms:

| Program | Captures | Library L−R offset (device samples) | IQR | µs |
|---|---|---|---|---|
| S | 368 | 0.5037 | 0.50366–0.50383 | 19.10 |
| M | 316 | 0.5035 | – | 19.09 |
| L | 180 | 0.5021 | – | 19.04 |

Neither 0.5 nor 0.506 fits. 0.5037 sits between 91 and 92 PROM steps, so a fixed L/R chain-delay difference is folded in. Keeping 0.5, with everything else exact, would cap the L channel:

| Band | L-channel limit |
|---|---|
| 1.6–3.2 kHz | −53 dB |
| 3.2–6.4 kHz | −47 dB |
| 6.4–10 kHz | −43 dB |
| 10–12.3 kHz | −40 dB |

The benches cannot see this because they align L and R separately. **Recommend 0.5037 (S), applied together with the chain fix.**

## 5. M, L, N, E, host rates, parameter counts

### M and L (modulation on; fit: 6 numbers per capture on cluster + gap + first pass; R)

| Window | M plugin | M + chain3 | M bench | M repeat | L plugin | L + chain3 | L bench | L repeat |
|---|---|---|---|---|---|---|---|---|
| Cluster + gap | −30.5 (0 / 458) | −51.4 (312 / 458) | −53.7 (403) | −54.8 | −32.0 (0 / 179) | −51.8 (174 / 179) | −52.6 | −52.6 |
| First pass | −29.9 | −48.6 | −48.2 | −49.0 | −31.6 | −51.6 | −45.6 (old law file?) | −54.1 |
| Recirc 1 / 2 / 3 | −18.1 / −12.8 / −8.9 | same | −20.8 / −16.0 / −12.6 (static) | −19.7 / −13.7 / −9.8 | −14.5 / −12.0 / −9.5 | same | −19.1 / −15.9 / −12.7 | −17.6 / −13.4 / −10.0 |

- **Reverb statistics match as well as the static bench does:**
  - decay slope medians ≤ 0.25 dB/s, 90th percentile 0.4–2.7 dB/s (M) and 0.2–0.5 dB/s (L);
  - level medians ≤ 0.08 dB;
  - EDC decay-time ratio median 0 (90th percentile 1.2–7.6 %).
- **Decorrelation through the sweep:** 6 start phases, 15 pairs per group. At recirc 1:

  | Group | Tomkin | Python sim, θ 0.15 rms | Plugin |
  |---|---|---|---|
  | M 1.2_5.0_0.9 | −17.0 | −14.7 | −13.8 |
  | L 1.2_4.5_0.8 | −20.4 | −16.9 | −14.2 |
  | L 1.2_20.0_1.0 | −18.0 | −14.3 | −13.3 |

  The plugin decorrelates 2–6 dB more than Tomkin in 6 of 8 groups; two groups are at Tomkin's level. About 0.8 dB comes from the evenly spaced phases. Tomkin sits between the θ 0.11 and 0.15 rms simulations, so the adopted 0.15 rms looks high (H).

### N (every 4th capture, 261)
- 186 of 261 have a bench-fitted branch offset other than +1 (82 fractional at 0.5–0.7, 104 at −64 to −165 or trims). The plugin cannot follow them and nulls at 0 to −10 dB there.
- On the 75 captures with offset ≈ +1:

  | Window (s) | Plugin (excess) | + chain3 (excess) | Bench (excess) |
  |---|---|---|---|
  | 0–0.05 | −29.8 (+50) | −50.3 (+30) | −56.2 (+24) |
  | 0.05–0.15 | −29.9 (+44) | −50.6 (+23) | −61.6 (+12) |
  | 0.15–0.3 | −29.9 (+39) | −50.1 (+19) | −58.2 (+9) |
  | 0.3–0.6 | −29.8 (+30) | −48.0 (+12) | −51.6 (+7) |
  | 0.6–1.0 | −29.4 | −38.7 | −39.0 |

- Per band, 0.05–0.6 s:

  | | 20–100 Hz | 0.1–1 kHz | 1–5 kHz | 5–10 kHz | 10–12 kHz |
  |---|---|---|---|---|---|
  | Plugin | −27.4 | −37.1 | −35.0 | −27.3 | −9.6 |
  | + chain3 | −39.2 | −57.7 | −54.4 | −47.8 | −41.9 |
  | Bench | −50.9 | −61.3 | −60.6 | −58.2 | −49.5 |

- The 5–10 dB left between + chain3 and the bench is the bench's measurement-side terms (EQ, HF correction, sub-100 Hz transfer) and its per-capture branch offset.

### E (e_pulse against the Altiverb REV-250 family)
- The plugin's device-rate core, with its own moving reads and 5 start phases. Its [0.77, 0.23] output kernel is divided out and e_pulse's Altiverb chain kernel applied.
- Result: 4,588 of 5,358 cells inside, against the bench's 4,601.
- Families differ both ways:

  | Family | Plugin | Bench |
  |---|---|---|
  | levs_e L [ext] | 41 / 56 | 52 / 56 |
  | ts L | 34 / 56 | 45 / 56 |
  | lr_l | 2 / 8 | 5 / 8 |
  | d50 R | 48 / 56 | 39 / 56 |

- The host-rate chain cannot be tested here (the Altiverb chain is another unit's).

### Host rates
- At 44.1 and 48 kHz, every program, modulation on and off, matches 96 kHz to −88.3 to −88.9 dB over 20 Hz–12.8 kHz (worst band −84.7 dB).
- A constant 10–11 µs offset, the uncompensated latency fraction, is absorbed by the per-capture fit.
- The time variance is locked to the device rate, so it is identical at every host rate.

### Parameter counts
- S: 371 fitted numbers, plus 3 per capture.
- M / L: S's shared laws, output kernels of 21 / 61 taps, a 3-value b table and θ (2 numbers), plus 6 per capture.
- N: 4 per capture here (bench 5).
- Bridge: not fitted except the S&H pole 0.301.
- Constants stored in the generated headers:

  | Header | Constants |
  |---|---|
  | `s_spec.h` | 3,728 |
  | `ml_spec.h` | 597 |
  | `n_spec.h` | 2,546 |
  | `e_spec.h` | 430 |
  | `bridge_spec.h` | 228 |

- A chain3 correction as stored would add about 3,700 taps per channel at 96 kHz, or a compact form still to be fitted.

## 6. Failures and proposed fixes (none applied)
1. **All programs fail −50 dB against Tomkin.**
   - Cause: the netlist chain vs chain3.
   - Fix: a per-channel correction on the bridge's output side fitted to chain3 ÷ bridge, held out by RT step.
   - Split it into the unit part (the skirt above 8 kHz, the 0.7° mid-band term) and the capture part (below 80 Hz, ~1.8 Hz corner), with the capture part off or behind a switch.
   - A 33–65-tap device-rate FIR covers 400 Hz–12.3 kHz; the low band needs pole-zero sections.
2. **L offset:** use 0.5037, together with fix 1.
3. **M / L θ depth:** re-check 0.15 rad rms; the sweep pairs suggest 0.11–0.13.
4. **N branch offset:** a fixed +1 is right if the offset is a load-time state. The 82 sub-sample offsets in the sample are not explained by a load phase alone; see NonLin_program.md §9.4.

## 7. Reproduce
```bash
cmake -S plugin -B plugin/build-render -DCMAKE_BUILD_TYPE=Release -DRV_BUILD_RENDER=ON
cmake --build plugin/build-render --target REV246Render --parallel 3
# each line: .venv/bin/python tools/memguard.py --limit-gb 6 -- .venv/bin/python tools/plugin_check.py <args>
#   --stage render --set all        (about 1,900 renders, about 6 GB, about 2 min)
#   --stage align ; --stage chain
#   --stage s --workers 4           (about 57 min)
#   --stage loff --workers 4 ; --stage rates
#   --stage ml --workers 4          (about 23 min)
#   --stage sim --workers 4 ; --stage n --workers 4 ; --stage e ; --stage report
```
