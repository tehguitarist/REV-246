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
