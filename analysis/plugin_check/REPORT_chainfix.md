# The bridge's chain correction: unit part of chain3, the 1.8 Hz corner, L at 0.5037 (smoke test)

*30 October 2026. Follows `REPORT.md` (phase 5).*
- *Tools:*
  - `tools/plugin_check.py --stage chainfit` (the fit, writes `chain_fix.json`);
  - `tools/plugin_codegen/bridge.py` (to `bridge_spec.h`);
  - `PLUGIN_CHECK_TAG=fix` for the smoke runs (results go to `*_fix`).
- ***This is a smoke test, not the full validation.*** *Sample: 6 S settings (24 captures), 2 M and 2 L groups (10 captures), and 3 N captures with branch offset +1.*
- *Values are medians of R (L in brackets) and are nulls in dB (error energy over capture energy).*

## 0. Verdict
- **The correction closes the chain gap.** On the smoke sample every program reaches its bench or does better:
  - **S, 50–150 ms:** −55.9 (−55.8) dB, 24 / 24 captures at −50 or lower. Before: −29.0. The repeat floor is −56.5, so the excess is +0.6 dB. The bench reaches −52.7 on the same captures.
  - **M / L, cluster + gap:** −60.6 to −62.7 dB. Before: −30.5 to −32.0. Bench: −52.5 to −53.5.
  - **N, 0.05–0.15 s:** −59.9 to −62.0 dB. Bench: −61.6 to −63.7.
- **The 1.8 Hz corner belongs in the chain.**
  - S at 20–40 Hz, 50–150 ms: −49.2 dB with the corner, −24.8 without it, −35.4 with chain3 swapped in.
  - So chain3's own shape below 40 Hz is a kernel artefact, and no capture residual is left.
  - Applying chain3's leftover low band in post makes 20–40 Hz worse (−33.9). Run the full validation with nothing in post.
- **Cost:**
  - 2 × 65 FIR taps plus one first-order section per channel at the device rate;
  - latency +17 device samples (+29 host samples at 48 kHz);
  - CPU +0.2–0.3 % of one core.

## 1. The split (unit part / capture part)
How the ratio was formed:
- chain3 over the netlist chain, per channel.
- The bridge's L is taken 0.5037 device samples late; chain3's L is read at the cache's 19.19 µs.
- With that, the ratio's L − R delay is 0.0002 device samples, which confirms 0.5037.
- R's best delay (0.788 device samples, a time-origin convention) and R's level are taken out of both channels.

| Part | What | Where it goes | Label |
|---|---|---|---|
| Output-filter skirt above 8 kHz | phase −1 to −70°; magnitude −0.4 to +1.0 dB at 9.5–12 kHz, −3 to −4 dB at 12.3 kHz | FIR, always on | S |
| Mid-band ~0.7° term | Mostly the 1.8 Hz corner's 1/f tail traded against the best delay. An HP alone plus a delay fits 25 Hz–3.2 kHz with fc 1.9–2.1 Hz to 0.20–0.26° rms. Adding a constant phase gives 0.31–0.41°, fc 1.5–1.7 Hz and 0.15–0.21° rms, a gain not separable from chain3's ±0.2° / ±0.1 dB ripple. | corner + FIR | S (corner tail) / H (any constant-phase rest) |
| ~1.8 Hz corner | Tomkin's 6.25 Hz total over the schematic's 4.46 Hz | device-rate first-order HP, always on (user, 30 October) | H: firmware DC blocker or Tomkin's rig |
| L / R difference | L +0.10 dB re R in the mid band; −0.1 to −0.9 dB and −0.3 to −2° above 10 kHz | per-channel FIR | S (shape) / H (unit's output filters or the interface) |
| chain3 below ~40 Hz beyond the corner | the 39 ms kernel's truncation / DC term | nowhere: the captures reject it (§3) | S |
| CLA / Altiverb 40–60 Hz bump | not Tomkin's | the separate XFMR stage, unchanged | – |

Where the corner is labelled in the docs:
- `docs/analog_results.md` lines 17, 20, 185 and 192;
- `docs/NEXT_STEPS.md` line 295, which calls it per-capture-chain;
- `docs/EMT246_unit_facts.md` line 84, plus the new row.

## 2. The correction and its fit error against the unit-part target
- **Form:** HP 1.8 Hz (bilinear, prewarped), then a 65-tap FIR per channel, bulk delay 17 device samples, unity at 1 kHz on R.
- **Fit:** linear least squares on the complex response, with equal weight per scoring band from 20 Hz to 12.3 kHz and weight 0.05 at 12.3–12.8 kHz.
- **Parameters:** 130 FIR coefficients fitted to the chain ratio. The corner comes from `analog_results.md`. The delay and the length were chosen from a scan; there are no per-capture numbers.
- **Top of the band:** maximum gain from 12.3 kHz to the device Nyquist is −3.9 / −3.2 dB (L / R).

| Band (Hz) | No correction | Corner + FIR (C++, 48 and 96 kHz) | FIR alone, no corner |
|---|---|---|---|
| 20–40 | −22.1 / −21.9 | −36.5 / −35.7 (target is chain3's artefact) | −22.6 / −22.4 |
| 40–80 | −29.0 / −28.0 | −43.9 / −45.8 | −31.3 / −30.5 |
| 80–160 | −34.6 / −33.2 | −56.0 / −59.4 | −44.0 / −43.3 |
| 160–270 | −37.1 / −35.8 | −53.3 / −52.2 | −40.9 / −40.1 |
| 270–400 | −38.6 / −37.0 | −54.6 / −57.1 | −39.5 / −39.2 |
| 400–800 | −38.0 / −37.2 | −58.4 / −61.1 | −43.5 / −43.3 |
| 0.8–1.6 k | −35.1 / −36.1 | −63.3 / −64.9 | −45.7 / −45.7 |
| 1.6–3.2 k | −33.3 / −36.4 | −66.4 / −67.1 | −48.5 / −48.4 |
| 3.2–6.4 k | −35.2 / −34.5 | −63.9 / −66.6 | −50.9 / −50.9 |
| 6.4–10 k | −24.5 / −24.9 | −60.2 / −62.0 | −53.8 / −53.5 |
| 10–11 k | −13.7 / −14.1 | −60.7 / −63.2 | −58.9 / −58.2 |
| 11–12.3 k | −3.7 / −4.0 | −51.1 / −54.5 | −50.0 / −50.0 |

FIR length scan (R, corner in, best delay):

| Taps | 80 Hz–10 kHz | 11–12.3 kHz |
|---|---|---|
| 17 | ≥ −29 | −21 |
| 33 | ≥ −42 | −37.5 |
| 49 | ≥ −53 | −46.6 |
| 65 | ≥ −52 | −52.9 |

A lower-latency alternative: 2 biquads + a 33-tap FIR with only 2 samples of delay reached −51 dB at 11–12.3 kHz but −56 to −59 elsewhere. I kept the 65-tap FIR for accuracy.

## 3. Smoke-test nulls
The plugin goes through the benches' capture path and their per-capture fits.

### S: 6 settings across RT 0.3–2.5, 24 captures
- **Fit:** shift-free, 3 numbers per capture.
- **Held out:**
  - The correction comes from the chain ratio, not from captures.
  - chain3 was fitted on cluster-only windows, so the reverb windows are not in any fit.
  - The S laws stay on the RT-step fold tables.

| Window (ms) | Before | Before + chain3 | **After** | After, no corner | Bench | Repeat | Excess | ≤ −50 dB |
|---|---|---|---|---|---|---|---|---|
| 0–50 | −28.3 (−27.7) | −49.6 (−51.3) | **−52.3 (−52.4)** | −49.4 | −51.2 | n/a | – | 22 / 24 |
| 50–150 | −29.0 (−28.3) | −52.4 (−51.9) | **−55.9 (−55.8)** | −51.4 | −52.7 | −56.5 | +0.6 | 24 / 24 |
| 150–300 | −30.5 (−29.9) | −51.5 (−50.8) | **−52.9 (−53.6)** | −49.2 | −49.6 | −53.0 | +0.1 | 16 / 24 |
| 300–600 | −28.7 | −48.6 | **−49.3 (−49.5)** | −42.9 | −44.0 | −42.9 | −6.4 | 10 / 24 |
| 600–1000 | −33.0 | −40.0 | **−40.2 (−40.0)** | −37.2 | −38.8 | −36.8 | −3.4 | 0 / 16 |

- "No corner" divides the plugin's HP out exactly and uses the bench's refit, not the shift-free fit.
- The capture residual in post changes no window here by more than 0.1 dB.

### S, full band, R, 50–150 ms (0–50 ms in brackets)

| Band (Hz) | Before | Before + chain3 | **After** | After + residual | No corner | Bench + low band | Repeat |
|---|---|---|---|---|---|---|---|
| 20–40 | −24.5 (−21.1) | −35.4 | **−49.2 (−47.3)** | −33.9 | −24.8 | −41.4 | −27.4 |
| 40–80 | −29.4 | −49.8 | **−53.1 (−49.4)** | −53.1 | −30.4 | −52.3 | −37.5 |
| 80–160 | −33.8 | −56.9 | **−59.7** | −59.9 | −36.4 | −51.7 | −49.3 |
| 160–270 | −36.0 | −48.5 | **−60.6** | −60.7 | −41.6 | −46.4 | −59.4 |
| 270–400 | −37.0 | −58.3 | **−60.0** | −60.1 | −46.4 | −53.5 | −58.4 |
| 400–800 | −36.4 | −58.3 | **−57.9** | −58.0 | −51.8 | −52.9 | −58.9 |
| 0.8–1.6 k | −34.1 | −59.0 | **−58.1** | −58.1 | −53.6 | −52.1 | −58.6 |
| 1.6–3.2 k | −32.2 | −59.5 | **−59.9** | −60.0 | −57.2 | −55.0 | −58.4 |
| 3.2–6.4 k | −28.9 | −53.0 | **−55.4** | −55.4 | −54.9 | −55.0 | −57.7 |
| 6.4–10 k | −26.8 | −49.4 | **−54.2 (−53.8)** | −54.2 | −53.0 | −49.2 | −53.6 |
| 10–11 k | −16.5 | −48.6 | **−50.9 (−48.9)** | −50.9 | −48.8 | −49.7 | −47.7 |
| 11–12.3 k | −12.9 | −44.3 | **−44.2 (−42.5)** | −44.2 | −44.8 | −44.3 | −43.0 |

- From 400 Hz to 3.2 kHz, after matches before + chain3 within 1.2 dB.
- 11–12.3 kHz sits at the repeat floor.

### M / L: 4 groups, 10 captures
- **Fit:** 6 numbers per capture, on cluster + gap + first pass.
- **Held out:** cluster + gap overlaps the windows chain3 was fitted on, so it is not fully held out from the chain fit. Everything after the first pass is predicted.
- Values are R (L).

| Group | Window | Before | Before + chain3 | **After** | Bench law | Bench repeat | ≤ −50 dB |
|---|---|---|---|---|---|---|---|
| M 0.1_1.8_0.4 | cluster + gap | −30.5 | −51.4 | **−62.7 (−63.8)** | −53.5 | −54.7 | 3 / 3 |
| | first pass | −29.9 | −51.0 | **−54.0 (−58.6)** | −49.3 | −49.8 | 3 / 3 |
| M 1.2_5.0_0.9 | cluster + gap | −30.5 | −51.3 | **−62.0 (−63.5)** | −53.4 | −54.5 | 3 / 3 |
| | first pass | −29.9 | −46.6 | **−47.6 (−52.7)** | −46.1 | −42.8 | 0 / 3 |
| L 0.1_10.0_0.3 | cluster + gap | −32.0 | −51.7 | **−60.6 (−61.2)** | −52.5 | −52.5 | 1 / 1 |
| | first pass | −31.6 | −52.1 | **−56.4 (−58.0)** | −45.5 | −54.8 | 2 / 2 |
| L 1.2_4.5_0.8 | cluster + gap | −32.0 | −51.7 | **−60.6 (−61.3)** | −52.5 | −52.5 | 1 / 1 |
| | first pass | −31.6 | −51.6 | **−55.0 (−56.8)** | −45.4 | −58.3 | 2 / 2 |

- **Recirculation passes are unchanged by the chain fix:** M −16.0 to −18.2 and L −13.9 to −14.3 at recirc 1. They are time-variant (`ml_structure.md` §4–5), so no null bar applies there.
- **Beating the bench by about 9 dB at cluster + gap is unexplained (H).** The bench model assumes chain3 exactly, while the plugin now has the analytic corner where chain3 has its artefact. The full run must confirm it.

### N: 3 captures with branch offset +1, R
Fit: 4 numbers per capture.

| Window (s) | Before | Before + chain3 | **After** (excess) | Bench (excess) |
|---|---|---|---|---|
| 0–0.05 | −29.8 | −50.0 to −50.7 | **−54.4 to −57.4** (+22 to +27) | −54.2 to −56.7 |
| 0.05–0.15 | −29.8 to −29.9 | −50.5 to −51.7 | **−59.9 to −62.0** (+14) | −61.6 to −63.7 (+12) |
| 0.15–0.3 | −29.7 to −29.9 | −48.0 to −51.5 | **−52.4 to −62.0** (+6 to +14) | −52.6 to −64.4 |
| 0.3–0.6 | −27.4 to −29.9 | −30.9 to −51.4 | **−30.9 to −62.2** | −30.9 to −64.0 |

Per band, 0.05–0.6 s (blank cells were not tabulated in the smoke summary):

| Band | Before | Before + chain3 | After | After + residual |
|---|---|---|---|---|
| 20–100 Hz | −27.6 | −40 | **−49.7 to −52.6** | −40 |
| 0.1–1 kHz | | | −59.8 to −61.8 | |
| 1–5 kHz | | | −59.9 to −63.7 | |
| 5–10 kHz | | | −57.2 to −61.2 | |
| 10–12 kHz | −9.6 | −43 to −44 | **−43.2 to −47.9** | |

### Host rates
At 44.1 and 48 kHz, all 6 programs (modulation on and off) match 96 kHz to −88.5 to −88.9 dB broadband; the worst band is −85.6 dB.

## 4. Code, latency, tests
- **`DeviceBridge`:**
  - `correct()` runs per channel on the device samples before the S&H lag: the HP, then the 65-tap FIR from a doubled ring, all in double, with no allocation.
  - `lLateDeviceSamples` comes from the generated header (0.5037).
  - `correctionDelay()` = 17.
  - The stale AdcClipper comment is fixed.
- **`ReverbEngine`:** `hostAlign = round((54 + 17) × host / device)`.

  | Host rate | Before | After |
  |---|---|---|
  | 44.1 kHz | 111 | 140 |
  | 48 kHz | 117 | 148 |
  | 96 kHz | 211 | 272 |

- **Codegen:**
  - `bridge.py` reads `chain_fix.json` and writes `lLateDeviceSamples`, `fixHighPassHz`, `fixDelayDeviceSamples`, `fixFirL` and `fixFirR`.
  - The R and L test references and the alias reference include the correction.
  - The unit-part target is written at 6 points per band.
- **Tests:**
  - REV246DspTests (55 cases) and REV246Tests (52 cases) pass.
  - The bridge matches its reference to −110 to −130 dB, and L its own reference to −110 to −122 dB.
  - New test of the corrected chain against the target: the §2 table. Bars: −50 dB from 80 Hz up, −42 at 40–80 Hz, −32 at 20–40 Hz.
  - The L test checks L / R against the references to −181 dB.
  - At 1 kHz, L reads 19.05 µs after R against 19.10 µs for 0.5037 device samples; the L FIR's own phase accounts for the difference.
- **Cost:** 0.9–1.1 % of one core per engine at 48 kHz and 1.1–1.4 % at 96 kHz (EngineTests). Before: 0.5–0.9 % and 0.8–1.1 %.

## 5. Not done / next
- **Full s / ml / n runs** are waiting for approval. Estimates are in the summary above. Run them with `PLUGIN_CHECK_TAG=fix`, with the post term dropped.
- auval was not run.
- The E stage was not rerun: its device-rate renders don't pass through the bridge.
- loff was not rerun. The new alignment shows L − R after the cache read at −0.00005 device samples (before: −0.0067).

## 6. Full validation (30 October 2026, `PLUGIN_CHECK_TAG=fix`, nothing in post)

Whole library, same fits and hold-outs as `REPORT.md` (S: RT-step folds; M / L: all after the first pass predicted). Summaries: `summary_{s,ml,n}_fix.json`, `loff_fix.json`, `rates_fix.json`. Medians, dB, R (L).

**S (127 settings, 491 captures), shift-free fit, 3 numbers per capture**

| Window (ms) | Plugin | Bench | Repeat | Excess | ≤ −50 dB |
|---|---|---|---|---|---|
| 0–50 | −52.4 (−52.9) | −51.5 | n/a | – | 442 / 491 |
| 50–150 | −55.8 (−55.8) | −53.0 | −57.2 | +1.4 | 472 / 491 |
| 150–300 | −53.4 (−53.4) | −51.4 | −54.5 | +1.1 | 381 / 491 |
| 300–600 | −49.4 (−49.4) | −46.7 | −47.6 | −1.8 | 213 / 484 |
| 600–1000 | −42.6 (−42.3) | −39.6 | −40.0 | −2.6 | 59 / 414 |

**M / L (modulation on; 6 numbers per capture fitted on cluster + gap + first pass)**

| | Plugin | Bench | Repeat | ≤ −50 dB |
|---|---|---|---|---|
| M cluster + gap | −62.5 (−63.7) | −53.7 | −54.8 | 417 / 458 |
| M first pass | −50.0 (−55.4) | −48.2 | −49.0 | 226 / 457 |
| L cluster + gap | −57.3 (−56.7) | −52.6 | −52.6 | 178 / 179 |
| L first pass | −55.0 (−56.9) | −45.6 | −54.1 | 359 / 359 |
| M recirc 1 / 2 / 3 | −18.3 / −12.8 / −8.9 | static −20.8 / −16.0 / −12.6 | −19.7 / −13.7 / −9.8 | time-variant: statistics |
| L recirc 1 / 2 / 3 | −14.6 / −12.1 / −9.5 | static −19.1 / −15.9 / −12.7 | −17.6 / −13.4 / −10.0 | time-variant: statistics |

- Cluster + gap and first pass are inside the per-capture fit window, and chain3 itself was fitted on cluster-only regions, so these windows are **not held out**; beating the repeat floor by 5–8 dB there reflects that (the repeat pair carries two captures' noise and the capture-side slow term). Not evidence on its own.
- EDC decay-time ratio median 0.00 in every octave (90th percentile ≤ 7.6 % M, ≤ 4.8 % L).
- Decorrelation: the plugin decorrelates 0.8–5.6 dB (median ≈ 2.6) more than the Python sweep simulation at the same θ = 0.15 rad rms, which itself matches Tomkin's pairs (median −0.4 dB here; `ml_structure.md` §9.7). A port difference, under investigation.

**N (every 4th capture, 261):** 86 / 256 captures reach −50 dB at 0.05–0.15 s; 75 of the 261 have the bench's branch offset +1, which the plugin fixes, so the rest are the open branch-offset item. All-capture medians are −32 dB for that reason.

**L offset after the fix:** residual library-vs-plugin offset leaves an L-channel limit of −61 dB at 10–12.3 kHz (−64 at 6.4–10 kHz), below the bar.

## 7. M / L at θ = 0.13 rad rms (31 October 2026)

Section 6 scored M / L with renders made before the θ refit (0.15 → 0.13, `ml_structure.md` §9.8). Rerun with the
same fits and windows: `--stage render --set ml / sim`, `--stage ml`, `--stage sim`. The old renders and results are in
`renders_fix_theta015/` and `ml_fix_theta015/` (gitignored); `summary_ml_fix.json` and `sim/` now hold θ = 0.13.
S, N, E and the host rates were not rerun: a 1-in-25 sample of their render lists, rendered by the 31 October build,
equals `renders_fix/` to −197 dB or better.

Medians, dB, R (L), the plugin against each capture; repeat = Tomkin's own capture pairs.

| | θ 0.15 | θ 0.13 | Repeat |
|---|---|---|---|
| M cluster + gap | −62.5 (−63.7) | −62.5 (−63.7) | −54.8 |
| M first pass | −50.0 (−55.4) | −50.0 (−55.4) | −49.0 |
| M recirc 1 / 2 / 3 | −18.3 / −12.8 / −8.9 | −18.8 / −13.4 / −9.6 | −19.7 / −13.7 / −9.8 |
| L cluster + gap | −57.3 (−56.7) | −57.3 (−56.7) | −52.6 |
| L first pass | −55.0 (−56.9) | −55.0 (−56.9) | −54.1 |
| L recirc 1 / 2 / 3 | −14.6 / −12.1 / −9.5 | −15.3 / −12.8 / −10.1 | −17.6 / −13.4 / −10.0 |

- The static windows are unchanged, as they must be (θ only acts on the recirculation).
- The recirculation nulls are 0.5–0.8 dB closer to the captures. M sits 0.2–0.9 dB from Tomkin's repeat floor. L is
  0.0–2.3 dB from it, and nearest by recirc 3.
- EDC decay-time ratio: median 0.00 in every octave for both programs; 90th percentile ≤ 7.2 % (M) and ≤ 4.7 % (L),
  as before.
- Decorrelation (the sweep pairs, `sim/`; the plugin's renders at different start phases against Tomkin's capture
  pairs, plugin minus Tomkin over the 8 settings × 6 windows): median **+1.7 dB** (|p90| 4.5) at 0.13, against
  +2.8 dB (|p90| 5.7) at 0.15. θ 0.13 is closer in every setting and window. Five settings are within ±2.4 dB (M
  −2.4…+1.9, L 2.0_2.6_0.5 −0.8…−0.6). Three L settings still decorrelate 3.2–5.3 dB more than Tomkin's pairs (L
  0.1_10.0_0.3, 1.2_20.0_1.0, 1.2_4.5_0.8), each by a nearly constant amount across windows. Open (H): a
  setting-dependent depth, or those pairs sharing more than the modulation (e.g. close in time). The θ fit's range
  (0.11–0.144, §9.8) brackets it; not changed here.
- Run note: 6 workers passed 8 GB on the long L groups (LF 2.0, RT ≥ 2.6), and memguard stopped the run at 314 / 340.
  The resume with 4 workers peaked at 7.2 GB.

## 8. Audit metrics on the plugin's renders

The pulse check's metrics (G1–G6 of `tools/pulse_check.py`) are now run on the plugin's renders in place of the research
model. The new code is `tools/plugin_pulse.py`, and the results are in `analysis/plugin_pulse/` (`s/`, `ml/`,
`report.log`, `summary.json`). Inputs:

- **Renders:** the existing ones in `renders_fix/96000/`, unchanged and not re-rendered. S is current; M / L are at
  θ = 0.13 rad rms with modulation on.
- **Into the bench's frame:** exactly as `plugin_check`'s `s` / `ml` stages do it. The fixed anchor comes from
  `align_fix.json`, followed by chain3 deconvolution and the device grid.
  - **S:** the per-channel tap-0 shift, then the 3-number refit (clock, shift, gain) per capture.
  - **M / L:** cluster part plus reverb part, with 6 numbers per capture fitted on cluster, gap and first pass. These
    are the same freedoms as the static model's fit.
- **Model and capture pairs:** pulse_check's own full run (`analysis/pulse/s`, `analysis/pulse/ml`), on the same
  captures.

**Checks:**
- The S nulls equal `s_fix/`'s `null_shift` to 0.1 dB.
- One M group's metrics and nulls equal `ml_fix/`'s plugin rows bit for bit.

**Run:**
- S: 128 settings, 512 captures, 20 min with 4 workers, peak 1.74 GB.
- M / L: 340 groups. 280 groups ran with 3 workers (15 min, 1.57 GB). The 60 long L groups (RT ≥ 10) ran with 1 worker
  (20 min, 2.06 GB; a single RT 20 group peaks at 3.3 GB).

**How to read the columns:** plugin vs capture / model vs capture [capture vs capture], medians, R.
- **S:** the brackets are the cross-take pairs, 0–30 and 60–90.
- **M / L:** the brackets are every pair in the group, with slipped captures left out.
- **Early metrics' brackets are not a tolerance.** The capture pairs differ in AMPL, so their bracket values for
  C50 / C80 / D50 / Ts, early levels, early IACC and early L/R mostly measure the cluster-level change between AMPLs.

### S

The plugin and the model are indistinguishable on every G2–G6 metric. Each is about 0.1 % from the capture on the
decay metrics (octave and third-octave decay, late decay, EDT) and about 0.02 dB on the level and clarity metrics
(EDC shape, C50 / C80 / D50 / Ts, tonal balance early and late, coherence, fine-structure sd, IACC).

- **Inside tolerance:** every capture in every band (octave decay 512 / 512, third-octave 511–512 / 512, for both the
  plugin and the model) is within its setting's largest cross-take difference (3 % floor).
- **Echo density** |diff| is 0.002 for the plugin, the model and the pairs alike.
- **The G5 fine-structure correlation is 1.000 everywhere**, so it does not discriminate here.
- **Common bias (H):** the ~0.1 % decay bias that plugin and model share, but the pairs do not, is the capture's fade
  compensation (the capture is divided by its ramp, the render is not). It is a floor of the measurement, not of the
  plugin.
- **EDR worst band:** 2.6 % for the plugin and 2.9 % for the model, against 0.9 % for the repeat (RT ≤ 0.6: 5.4 / 5.4
  [0.5]). Plugin and model share it, so it is a model-level item. The 90th percentile over bands is 0.2 / 0.2 [0.1].
- **G1 nulls**, held out (fold constants) apart from the 3-number fit:

| Window (ms) | Plugin | Model | Repeat |
|---|---|---|---|
| 50–150 | −55.7 | −53.0 | −56.5 |
| 150–300 | −53.2 | −51.2 | −54.2 |
| 300–600 | −48.9 | −46.4 | −47.2 |
| 600–1000 | −42.5 | −39.5 | −40.0 |

The plugin is 1–3 dB better than the model in every group. With the clean captures only, this is the same as §6.

### M and L: which metrics mean something

The plugin's render is time-variant and so are the captures. Plugin vs capture should therefore be read against capture
vs capture, not against the static model, which is a mean response and can sit closer to a capture than any other
realisation can.

- **Meaningful:** metrics computed on each IR and then compared. These are the decay times, EDT, EDC shape,
  C50 / C80 / D50 / Ts, tonal balance, echo density, coherence values, fine-structure sd (colouration), mixing time,
  IACC and L/R level.
- **Measure decorrelation, not accuracy:** the waveform-matching metrics. These are the nulls past the first pass, the
  fine-structure correlation (except L's window, which falls inside the deterministic first pass: 1.000), the
  narrow-band envelope and EDR per 1/24 octave.

**Results inside the capture tolerance:**

- **Octave decay at 1–8 kHz:**
  - M: 0.9 / 0.6 / 0.5 / 0.3 % against [0.8 / 0.6 / 0.5 / 0.2].
  - L: 0.9 / 0.5 / 0.4 / 0.3 % against [0.7 / 0.4 / 0.2 / 0.2].
  - Signed medians are within ±0.4 %.
  - Within the group's largest pair difference: 341–348 of 347–350 for M, 344–357 of 358 for L.
- **Late decay** at every octave (M 4.4–0.6 % against [4.6–0.6]; L 4.5–0.7 % against [4.9–0.6]).
- **Third-octave decay** at every band, at the repeat or within about 0.6 % of it.
- **Colouration:** fine-structure sd is within 0.01 dB of the capture (about 5 dB).
- **Late tonal balance** is within 0.1 dB at every band.
- **L/R:** late IACC within 0.001, L/R level within 0.01 dB, coherence within 0.01 (except M at 400 Hz, +0.014).
- **Early metrics:** C50 / C80 / D50 / Ts, early tonal balance and early IACC match the capture to 0.0–0.08 dB, because
  the cluster and first pass are deterministic and fitted. Here the plugin is better than the static model at
  100–160 Hz (M early level 0.01 against 1.29 dB; L C50 at 125 Hz 0.04 against 0.53 dB).
- **M echo density:** 0.022 / 0.020 [0.022].
- **M mixing time:** 2.7 % / 2.6 % [2.7 %].
- **Decorrelation:** M's nulls sit on the pairs' values (recirc 1–3: −18.8 / −13.5 / −9.7 against [−19.3 / −13.5 /
  −9.9]). The static model is 2–3 dB "better" there, as a mean response must be.

**Results where the plugin is clearly worse than the model and outside tolerance (all in the low band):**

- **M, 125 Hz decay (−5/−25 dB):**
  - Median error 2.6 % against the model's 1.9 % [1.9 %]; signed −1.2 %.
  - EDC deviation 0.81 dB against 0.62 [0.64].
  - Only 253 of 346 captures are within their group's pair spread, against 305 for the model.
  - The excess is spread over RT, worst at RT 0.8 (6.2 % against 2.7 %).
- **M EDT at 125–500 Hz:** 2.4 / 1.6 / 1.5 % against the model's 1.6 / 1.1 / 0.7 % [2.0 / 1.1 / 0.7]. At 500 Hz the
  excess grows with RT, up to 2.2 % against 0.7 % at RT 4–5.
- **L, 250 Hz decay:** 1.3 % against 0.9 % [0.9 %], worst at RT ≤ 3.6 (1.5–1.8 % against 0.7–0.8 %). Within the pair
  spread: 318 of 358, against 339 for the model.
- **L time structure:**
  - Mixing time is 3.6 % from the capture, against 1.3 % for the model [1.4 %].
  - Echo density |diff| is 0.011 / 0.008 [0.009].
  - Recirc 1 null is −15.4 against [−17.6], 2.2 dB more decorrelated than Tomkin's pairs, as already seen in §7's
    sweep pairs.

**Likely cause (H):** the plugin's time variance, not its static structure. The static windows still null at −50 to
−63 dB, and every excess is in metrics a realisation-to-realisation spread can move: few modes per band below 500 Hz,
and L's early recirculation. They sit where the plugin decorrelates more than Tomkin's pairs: L by 2 dB at recirc 1
here, and the +1.7 dB median of §7. The candidates are the same as §7's open item:

- θ too deep for some settings;
- the C++ port's modulation differing from the Python simulation, which matches Tomkin.

A check that fits within the no-hardware rule: the same metrics on two plugin renders of one setting at different
start phases (`renders_fix/*_sweep*` style). If plugin-vs-plugin shows the same low-band spread, the excess is
modulation depth and not a decay-law error.

## 9. Audit metrics on the plugin's N renders

The pulse check's per-IR metrics (G1–G6 of `tools/pulse_check.py`) were run on the plugin's N renders. The code is the
new `tools/plugin_pulse_n.py`. Results are in `analysis/plugin_pulse/n/` (one JSON per capture),
`analysis/plugin_pulse/report_n.log`, `summary_n.json` and `n_run.log`. `plugin_pulse.py` was not changed.

**Inputs:**
- **Renders:** the existing `renders_fix/96000/N_*.wav`, unchanged.
- **Into the bench frame:** `plugin_check.n_device`, with the plugin's fit (shift, clock, 2 gains) read from
  plugin_check's own result files: `n_fix/` for n = 0 and `analysis/nl/branch_slip/deslip/n_fix/` for the de-slipped
  captures.
- **Captures:** `plugin_check.n_keys()` (N_EVERY = 4).
- **Slips:** the library's slip is undone as in `n_deslip.py`, with n from `defects.NL_SLIPPED`, checked against the
  bench's o.
- **Bench model:** `nl_model.render` with the adopted bench's own per-capture o (moved by the slip), clock and front shift,
  and least-squares gains, against the same capture clock-corrected as the bench does it.

**Checks:** the recomputed window nulls equal plugin_check's files (plugin) and the bench files (n = 0 captures) to
0.0000 dB.

**Run:** 257 captures, 2 workers, 15 min. Peak 0.9 GB for one capture and 2.04 GB for the run.

**Excluded or apart:**
- The 4 SPAN 0 / LATE 0 captures were skipped (no branch).
- 18 captures whose bench o is not +1 plus a whole slip in −400…−1 are out ("unresolved": Tomkin's trims of §11.4 and the
  two off-grid offsets).
- 0.0_2.0_0.4 is out: the bench itself nulls at +5…+21 dB, so the file is not its label's response (a new defect
  candidate, not added to `defects`).
- The 33 trimmed-start captures are reported apart.
- 205 untrimmed captures remain.

**Tolerance:** the library has no repeat captures of an N setting (all 1,046 labels are distinct), so two stand-ins are
used.
- **Twin:** the capture against itself plus an independent stretch of the library's own capture noise (LATE 0 tails),
  3 realisations; the tolerance is the largest |error|. This is what a perfect, noise-free model would score: a lower
  bound.
- **LATE pair:** a real capture vs capture comparison on the early section. The capture is set against the one that
  differs only in LATE (SPAN ≥ 0.2, 127 pairs, 1 shift and 2 gains). It covers early levels, early IACC and L/R, and the
  0–50 ms null.

**Which metrics mean something on N:**
- **Defined and reported:** EDC shape (largest EDC difference to −30 dB, from the input), the late section's tail decay
  (EDC from 40 ms after the last main tap, −5/−25 dB, LATE > 0), C50 / C80 / D50 / Ts as energy partitions, early and late
  tonal balance, echo density, coherence, fine-structure correlation and sd, narrow-band envelope, IACC, L/R level, and
  nulls (N is time-invariant).
- **Restricted:** the EDC decay fits (octave or third-octave −5/−25, EDT, −25/−40) are kept only where the capture's EDC
  is straight within 1 dB over the fit range. That holds for 20–60 % of captures per band; gated, reverse and
  stepped-ramp envelopes fail it.
- **Dropped:** EDR per 1/24 octave (the same fit on non-exponential envelopes) and the mixing time (NED never reaches 0.9
  within 1 s on N; defined for 1 capture).

**Results** (205 untrimmed captures). Cells are plugin |err| / bench |err| [twin], medians, R.

| Family | Plugin / bench [tolerance] |
|---|---|
| Octave decay −5/−25, where straight | 0.04–0.05 % / 0.04–0.05 % [0.01] |
| Third-octave decay, where straight | 0.03–0.07 % / same [0.01–0.02] |
| EDT, where straight | 0.03–0.05 % / same [0.00] |
| Late decay −25/−40 | 0.08–0.15 % / same [0.02–0.08] |
| Tail decay (late section), broadband and octaves | 0.05–0.06 % / 0.05–0.06 % [0.01–0.02] |
| EDC shape, octaves, R and L | 0.01–0.03 dB / same [0.00–0.01] |
| Third-octave EDC shape | 0.01–0.05 dB / same [0.01] |
| C50 / C80 | 125 Hz 0.036 / 0.031 dB; elsewhere 0.003–0.009 / same [0.000–0.001] |
| D50 | ≤ 0.0009 / same |
| Ts | 125 Hz 0.22 ms; elsewhere 0.02–0.06 / same [≤ 0.005] |
| Early tonal balance, 0–80 ms | 100–250 Hz 0.024–0.044 / 0.017–0.046 dB; above 400 Hz ≤ 0.01 [0.000–0.001] |
| Late tonal balance | 160–400 Hz 0.008–0.010 / 0.004–0.006 dB; elsewhere ≤ 0.006 / same [0.000] |
| Echo density | 0.0009 / 0.0009 [0.0005] |
| Coherence, IACC, L/R | ≤ 0.0003, ≤ 0.0002, 0.000 dB, plugin = bench |
| Fine-structure correlation | 1.0000 / 1.0000 |
| Fine-structure sd | 0.001 dB / same [0.001], inside |
| Narrow-band envelope (p90 / worst band) | 0.17 / 0.77 dB, bench 0.16 / 0.78 [0.06 / 0.17] |
| LATE pairs (real capture vs capture) | early levels 100–250 Hz 0.03–0.06 / 0.03–0.06 dB [0.004–0.010]; 0–50 ms null −52.7 / −52.8 [pair −60.7] |

**G1 nulls (held out as §9.9),** R (L), dB, medians:

| Window (s) | Plugin | Bench | Excess (plugin) | ≤ −50 dB, plugin / bench |
|---|---|---|---|---|
| 0–0.05 | −53.3 (−51.7) | −53.4 | +26.5 | 148 / 147 of 205 |
| 0.05–0.15 | −58.0 (−60.0) | −59.2 | +14.7 | 184 / 183 |
| 0.15–0.3 | −57.1 (−58.6) | −56.2 | +11.4 | 158 / 158 of 202 |
| 0.3–0.6 | −50.9 (−51.5) | −50.7 | +8.5 | 97 / 97 of 192 |
| 0.6–1.0 | −40.2 | −40.1 | +6.4 | 17 / 14 of 134 |

Per band, 0.05–0.6 s (plugin / bench):

| Band | Plugin | Bench |
|---|---|---|
| 20–100 Hz | −49.4 | −49.4 |
| 100 Hz–1 kHz | −58.3 | −60.3 |
| 1–5 kHz | −58.6 | −58.5 |
| 5–10 kHz | −55.6 | −55.1 |
| 10–12 kHz | −39.7 | −48.4 |

**Reading the results:**
- **The plugin and the bench model are the same on every metric family.** The medians agree to within 0.005 dB, 0.01 %
  or 0.001 in every band and group: EARLY > 0, EARLY 10, EARLY < 0, LATE 0, SPAN 0, slip undone and trimmed. The few
  cases where the plugin is worse by more than twice the tolerance are listed in the next bullet.
- **Where the plugin is clearly worse than the bench: 100 Hz–1 kHz.**
  - Early and late tonal balance at 160–400 Hz: 0.010 against 0.006 dB late; in 100–150 of 205 captures the plugin's
    error exceeds the bench's by more than twice the tolerance.
  - The 100 Hz–1 kHz null: −58.3 against −60.3 dB.
  - The 10–12 kHz null: −39.7 against −48.4 dB.
  - **Cause (S):** two bench-side capture terms the plugin does not and must not carry. They are the measured low-band
    transfer of the captures' processing (`NL_LOWBAND=2`, N's own table above 130 Hz, NonLin §11.8) and the per-run
    band-edge term above 10 kHz (`NL_HF`, §11.6). Both are measurement-side, not the unit. 20 Hz–100 Hz and 1–10 kHz are
    equal.
- **Against the tolerance:** plugin and bench alike are outside the noise twin, and outside the real LATE-pair spread of
  the early section (0–50 ms −52.7 against the pair's −60.7 dB; 100–250 Hz early levels 0.03–0.06 against 0.004–0.010 dB).
  This is the known model-level gap of §11.1 (8–13 dB above the repeat level), shared by plugin and bench. Inside or at
  tolerance: fine-structure sd and correlation, L/R level, and coherence in about half the captures.
  The absolute errors are tiny: decay 0.05 %, EDC 0.02 dB, levels ≤ 0.06 dB.
- **One clear shared failure: LATE 0.2** (20 captures).
  - The late section's tail decays 1.8 % fast: RT 0.244 against the capture's 0.248 s (all 16 measurable captures between
    −0.6 % and −2.0 %; tolerance 0.01 %).
  - Nulls collapse after 0.15 s: −38 / −21 / −16 dB at 0.15–0.3 / 0.3–0.6 / 0.6–1.0 s, every band about −41 dB.
  - Narrow-band envelope worst band 2.5 dB; echo density |diff| 0.010.
  - Plugin and bench are identical, so this is a model item, not the port. It is the "LATE 0.2 reads 0.23–0.29 s" of
    NonLin §10.1, attributed there to estimator bias, but it is the unit's response. H: the LATE 0.2 step's loop gains
    are not Jot gains for RT = 0.2 s (the unit's real RT at that step is about 0.25 s, or the gain is quantised).
    LATE ≥ 0.4 is at the 0.05 % level.
- **Trimmed starts:** they score like the rest after 0.15 s. 0–0.15 s and the bands are limited by Tomkin's trim (20–100
  Hz −33 dB for plugin and bench alike).
