# Programs M and L: clusters, loops and time variance (Q3)

Everything below comes from `.venv/bin/python tools/ml_structure.py` on the Tomkin **96 kHz** set
(log: `analysis/ml_structure/run.log`; data: `analysis/ml_structure/summary.json`). Program S is run
through the same code as the control and the calibration. Labels as in `EMT246_unit_facts.md`
(**C** confirmed, **S** strong, **H** hypothesis).

*October 2026. §12 replaces the library's deconvolution kernel with a per-channel chain kernel (`EMT_KERNEL=chain2`). Answers `NEXT_STEPS.md` Q3 and changes the picture of the "stochastic control" (§8 there). §8 below adds the line gains and filters (Q4–Q7 for M/L), §9 the time variance itself (Q3b).*

> **Summary.**
> - **Early section = program S's, time-scaled by 2 (M) and 4 (L).** Same 8-tap cluster with the same signs, the reverb onset and the L/R onset gap all scale (§1, §2).
> - **The late reverb is S's four loops, ×2 / ×4 and coupled (§7, superseding §3's reading).** S is four parallel loops: the product of four combs predicts its tail to **−45 dB held out**. M and L follow one five-term recursion each (L 10826/11384/11618/12176/23002, M 5416/5694/5812/6090/11506, §3). §7 shows why: the same four lines scaled (M 2641/2775/2919/3171, L 5279/5547/5837/6339) are cross-coupled by an orthogonal 2 + 2 matrix, whose characteristic polynomial has exactly these five lags (pair sums and the sum of all four) and their signs.
> - **The carry-over test (§7):** M's and L's diffuser, 20 read taps per channel, R routing and AP damping steps are S's, time-scaled (with small exact integer offsets); M and L add an output smoothing kernel. Held out, the first pass nulls as well as S's own compact model does.
> - **M and L are time-variant; S is not.** Separate captures of the same M or L setting decorrelate progressively (L: −25 dB at 0.15–0.3 s → −3 dB at 2–4 s), flat across frequency. S's repeat captures agree to −52…−59 dB once a per-capture clock ratio is fitted (§4.1). This is the first direct evidence of the "stochastic control" and its program dependence (§4).
> - **What varies (§9):** one common rotation angle of the coupling's two 2 × 2 blocks, periodic at 0.83 Hz with a random phase per capture, about 0.11–0.14 rad rms. Rendered through a simulated sweep, the rotating model decorrelates between renders the way Tomkin's AMPL pairs do, and it explains §8.3's small RT shortfall.
> - **Consequence for the null target:** past the first recirculation, M and L can't be nulled waveform-for-waveform against one capture; two captures of the same unit don't null each other (§5).
> - **Noise floor vs the CLA capture (§4.2):** not the same. CLA levels off at −71/−77 dB (L/R) below its 50–300 ms reverb level; the library shows no capture noise down to its 24-bit floor, about −88 dB below. What limits a comparison against the library is the capture clock, not noise: with it fitted, the S repeat-capture ceiling is about −57 dB.

## 0. Method notes

- **Device-rate IRs:** each file deconvolved with the library kernel and resampled to 19.2 MHz/728 (harness `core.deconvolve` / `to_device_grid`), onset at index 50, L shifted by the library L/R offset (19.19 µs).
- **AMPL pairs are repeat captures of the reverb.** The library has settings that differ only in AMPL (S: 128 groups of four, M: 160 groups of three, L: 180 pairs). AMPL scales the cluster only, so for an LTI unit the reverb tails of a pair are identical after alignment and one gain. This is the repeat-capture test that `NEXT_STEPS.md` §7–9 said was missing. It also gives the cluster directly, as the pair difference.
- **Tail predictor (characteristic polynomial).** Once the input has ended, every output of an LTI recursive structure obeys its characteristic polynomial, whatever the read-out taps. So the right lag set predicts the tail from its own past down to the noise, and a wrong one doesn't. Predictor: 300 short lags (spectral colour and band-limit; on their own they give −10…−11 dB) plus a 7-tap FIR at each candidate lag, fitted on one time window and scored on the next (held out). Controls: random lag sets of the same size within ±8 %.

## 1. Cluster patterns

From 40 AMPL pairs per program (pair difference, four largest peaks per channel in the first 2,000 samples, sub-sample peak positions):

| Program | L taps (samples after tap 0), signs | R taps, signs | Scatter over 40 pairs |
|---|---|---|---|
| S | 0, 54, 194, 288; + + + − | 0, 143, 194, 353; + + − + | ≤ 0.001 sample |
| M | 0, 108, 388, 576; + + + − | 0, 285, 388, 704; + + − + | ≤ 0.001 |
| L | 0, 216, 776, 1152; + + + − | 0, 570, 776, 1407; + + − + | ≤ 0.001 |

- **C (direct measurement):** the patterns, integer positions and signs.
- **The L-channel taps and the shared tap are exactly ×2 (M) and ×4 (L) of S.** The two R-only taps are not: S's 143/353 become 285/704 in M (2× − 1, 2× − 2) and 570/1407 in L (4× − 2, 4× − 5). The R-only taps come out 0.28–0.35 % short of 2× and 4×. Why is open (**H**: they are computed from a shared table with rounding, not stored per program).
- All 8 taps are equal in magnitude within each program (spread 0.11–0.13 dB). The pair difference measures AMPL-hi minus AMPL-lo, so the AMPL law is not measured here.
- The harness's `cluster_templates` hold S only and its search window is 16 ms, too short for L's 55 ms cluster; that is why M and L showed `no_cluster`.

## 2. Reverb onset after cluster tap 0

First sample of the low-AMPL file above −30 dB of its 0–0.3 s peak, cluster masked; median over 40 pairs (IQR within ±2 samples):

| | L | R | R − L |
|---|---|---|---|
| S | 240.5 | 572.5 | 332 |
| M | 474.4 | 1,133.5 | 659 |
| L | 948.4 | 2,268.4 | 1,320 |

S reproduces the CLA unit's 240/570 (`model_s_results.md` §4). M is 1.97–1.98 × S and L 3.94–3.96 × S in both channels, and the L/R gap is S's 332-sample read-out lag scaled by about 2 and 4 (**S**; a threshold onset is good to a few samples, so the exact factors are not resolved). So the whole early section (cluster, injection onset, L/R read-out lag) is S's, time-scaled.

## 3. Late reverb: the recursion

### 3.1 Library autocorrelation

Averaged decay-normalised tail ACF (300 Hz–6 kHz) over all files of a program: S gives exactly its loop set 1321/1389/1461/1587 and their multiples, which validates the method. M peaks at 2641, 5416, 6090, 8057 (= 2641 + 5416) and 8336; L at 5279, 10826, 12176, 34063, 35413, 23237, 16106 (= 5279 + 10826) and 17456 (= 5279 + 12176) *[exploratory]*. Sums of lags can't come from parallel loops, which give multiples only, so M and L are not S-style loop sets. These are the harness's "unknown loop sets" (`library_notes.md` §3): `5279-10826-12176…` is L and `2641-5416-6090…` is M.

### 3.2 The characteristic polynomial (held out)

| Program, file, windows | Lag set | Held-out null L / R (dB) | Params |
|---|---|---|---|
| **S** `12_0_0_10_25_10` (LF 1.0, HF 1.0, RT 2.5), fit 0.4–1.0 s, test 1.0–1.4 s | short block only | −9.9 / −10.6 | 300 |
| | 4 single lags | −10.9 / −11.6 | 328 |
| | single lags + pair sums | −11.2 / −12.2 | 370 |
| | **4 parallel loops (all 15 subset sums)** | **−45.9 / −45.0** (in-sample −52.5 / −53.4) | 405 |
| | 3 controls (15 random lags) | −9.4…−10.0 / −10.3…−10.8 | 405 |
| **M** `12_0_10_12_50_09` (LF 1.2, HF 0.9, RT 5), fit 1.0–2.3 s, test 2.3–3.3 s | short block only | −11.1 / −11.2 | 300 |
| | 4 mid taps | −11.9 / −12.3 | 328 |
| | **5 taps** | **−17.5 / −20.8** (in-sample −15.6 / −20.2) | 335 |
| | 4 mid taps as parallel loops (15 sums) | −21.3 / −21.1 | 398 |
| | 3 controls (5 random lags) | −9.7…−10.1 / −11.0 | 335 |
| **L** `12_0_0_12_200_10` (LF 1.2, HF 1.0, RT 20), fit 2.0–4.5 s, test 4.5–6.5 s | short block only | −10.5 / −10.5 | 300 |
| | 4 mid taps | −11.5 / −11.4 | 328 |
| | **5 taps** | **−17.5 / −17.6** (in-sample −17.6 / −17.6) | 335 |
| | 4 mid taps as parallel loops (15 sums) | −18.1 / −17.9 | 398 |
| | 3 controls (5 random lags) | −10.5 / −10.4 | 335 |

The parameter counts are dominated by the 300-tap short block, which only models spectral colour; the structural content is the lag set (4 or 15 lags for S, 5 for M/L). Noise floors are not available in the faded library files (`NEXT_STEPS.md` §2), so the nulls are given as ratios to the signal in the window; S's −45 dB shows how far the method reaches when the structure is right.

- **S is four parallel loops (S, now held out to about −45 dB).** Only the full product of four combs predicts the tail; single lags give 1 dB over the short block. This confirms the S loop model and shows S's tail is LTI to about −45 dB in the library unit, at LF 1.0/HF 1.0 and also at LF 2.0/HF 0.6 (−42 to −44 dB *[exploratory]*).
- **M and L: five lags; no further lag found.** The lags are 4 taps placed symmetrically about a centre N, plus one at 2N:

  | | N | Taps | 2N tap | Signs of the fitted coefficients |
  |---|---|---|---|---|
  | L | 11,501 | N ± 117, N ± 675 → 10826, 11384, 11618, 12176 | 23002 | + − − + and − |
  | M | 5,753 | N ± 59, N ± 337 → 5416, 5694, 5812, 6090 | 11506 | + − − + and − |

  Equivalently a + {0, p, q, p + q} with p = 558, q = 792 (L) and p = 278, q = 396 (M), and 2a + p + q. L's lags are 2× M's within 4–10 samples, not exactly.
- Random lag sets of the same size add nothing over the short block, and a greedy search for further lags (up to 2.2 × the largest) finds none that helps held out *[exploratory]*. Treating the four mid taps as parallel loops (which adds their pair sums, including 23002) gives the same null as the five taps, so the pair terms carry no extra information.
- **Reading (H; superseded by §7.3):** one recirculating delay line of 2N words (L 23 K, M 11.5 K, within the ~39 K budget) with five feedback taps, or the equivalent two lines of N with a coupling. The four mid taps cancel at DC (+ − − +), so at low frequencies the loop is the 2N tap alone; the mid taps shape the loop gain as −4β sin(ωp/2) sin(ωq/2) on the half-loop. With uniform decay that is a 2 × 2 rotation-like system in w = z^−N: modes decay by the same factor per N samples wherever |trace| ≤ 2√det.
- **The plateau (−18 to −22 dB) is not a missing lag.** The residual is flat across octaves (within about 2 dB from 125 Hz to 4 kHz, 3–8 dB higher in the 8 kHz octave next to the 10 kHz band edge; S: −56 to −43 dB), the same in-sample as held out, and still there later in the tail (L fitted 8–10.5 s, tested 10.5–12.5 s: −12.3 / −18.9 dB): a fixed fraction of the decaying reverb that no fixed recursion predicts. §4 shows why.

### 3.3 Lag precision

Sub-sample position of each lag FIR's peak (low-passed), over 12 settings each with HF ≥ 0.9 (M) or 1.0 (L), against S as the calibration:

| | Offset from the integer lag (mean, L / R) | Scatter across settings |
|---|---|---|
| S (4 loops 1321…1587, and the 4-fold sum 5758), LF 1.0 and HF 1.0 only | +0.04 / −0.03 (5758: +0.07 / +0.02) | 0.26 / 0.13 |
| M (5416, 5694, 5812, 6090; 11506) | −0.80 / −0.79 (11506: −0.75 / −0.74) | ≤ 0.04 |
| L (10826, 11384, 11618, 12176; 23002) | −0.82 / −0.83 (23002: −0.73 / −0.76) | ≤ 0.04 (23002 on L: 0.09) |

S sits within 0.05 sample of its integer lags (the remainder is the capture clock, `model_s_results.md` §3.1), so the method has no offset of its own. **M and L carry a common −0.8 sample on every lag**, stable to ±0.04 across 12 settings each (RT, LF and AMPL varied): each lag is d − 1 plus about 0.2 sample of group delay that all taps share. **H:** a filter in the common feedback path (e.g. a one-pole with p ≈ 0.17), and the integer lags are d − 1 (L 10825 … 23001). The sums between lags (§3.2) are exact either way, but the absolute integers are uncertain by one sample until the loop filter is modelled. With heavy HF damping (HF 0.1–0.2) the offsets spread to ±1 sample, as the loop low-pass adds frequency-dependent delay *[exploratory]*.

*Update (§7.3):* the lags are sums of four integer lines, and a single recirculation lands on tap + line within +0.05…0.08 sample (the capture clock), so the lines carry no −0.8-sample offset or loop-filter group delay. Why the predictor's lag FIRs peak 0.8 sample early is open.

## 4. Time variance: AMPL pairs as repeat captures

Median pair null after alignment and one gain per window (100 Hz–10 kHz):

| Window (s) | S, AMPL 90 vs 60 (40 pairs) L / R | M, 99 vs 60 (40) | L, 80 vs 0 (40) |
|---|---|---|---|
| 0.05–0.15 | −47.0 / −46.9 | −26.6 / −37.9 | −34.8 / (−1.4)¹ |
| 0.15–0.3 | −41.1 / −41.1 | −13.9 / −14.7 | −25.1 / −32.8 |
| 0.3–0.6 | −38.1 / −38.4 | −9.2 / −9.3 | −13.7 / −14.4 |
| 0.6–1.0 | −34.8 / −35.0 | −7.0 / −7.1 | −9.4 / −9.7 |
| 1.0–2.0 | −22.7 / −22.5 (9 pairs; ends of short files) | −5.9 / −5.1 | −5.3 / −5.4 |
| 2.0–4.0 | – | −7.1 / −7.1 (11) | −2.9 / −2.9 |

¹ L's R cluster runs to 1,407 samples (53 ms), so this window still holds R cluster taps that only one file of the pair has.

Per octave (R channel, median):

| Program, window (s) | 125 | 250 | 500 | 1k | 2k | 4k | 8k |
|---|---|---|---|---|---|---|---|
| S 0.15–0.3 | −64.7 | −62.7 | −58.8 | −53.7 | −47.8 | −41.7 | −36.5 |
| S 0.6–1.0 | −48.3 | −50.1 | −45.8 | −41.7 | −36.2 | −30.8 | −25.1 |
| M 0.15–0.3 | −16.1 | −15.5 | −14.1 | −14.3 | −14.5 | −15.5 | −15.9 |
| M 0.6–1.0 | −5.7 | −8.1 | −6.3 | −6.8 | −7.0 | −7.0 | −5.9 |
| L 0.15–0.3 | −32.6 | −33.6 | −34.4 | −31.2 | −34.6 | −33.7 | −32.9 |
| L 0.6–1.0 | −10.9 | −9.3 | −9.2 | −9.5 | −9.6 | −10.5 | −10.5 |

All three M pairings agree (R, 0.15–0.3 / 0.3–0.6 / 0.6–1.0 / 1–2 s): 60 vs 10 −15.2 / −10.5 / −8.0 / −6.4; 99 vs 60 −14.0 / −8.7 / −6.1 / −4.8; 99 vs 10 −15.4 / −9.4 / −6.6 / −4.9.

- **S: LTI within the capture.** Its pairs agree to −47 dB early and −35 dB at 0.6–1 s with a sub-sample shift only. The difference rises about 5 dB per octave and grows with time: a clock-rate difference between captures, not the unit. With a per-capture clock ratio fitted, S pairs agree to −52…−59 dB out to 1 s (§4.1). L and R agree within 0.3 dB, so the library shows **no L-only excess** in S (the CLA unit's L read-out jitter, `model_s_results.md` §13, is not seen here at this depth).
- **M and L: the reverb decorrelates between captures, progressively and flat across frequency** (**S**). Both agree early (M −27/−38 dB at 0.05–0.15 s, L −25/−33 dB at 0.15–0.3 s, before its first recirculation at 0.41 s), then lose coherence as the reverb recirculates until they are uncorrelated (−3 dB) by 2–4 s (L) or nearly so (M, −5 to −7 dB after 1 s).
- **It is per-capture randomness, not AMPL.** M has three AMPL values; all three pairings (60–10, 99–60, 99–10) decorrelate equally after 0.15 s.
- **Not delay modulation, not clock drift.** Both would grow about 6 dB per octave, as S's small difference does; M's and L's decorrelation is flat from 125 Hz to 8 kHz in every window. (Without the sub-sample alignment step, L's 0.15–0.3 s window looks as if it rises with frequency; that was an alignment artefact *[exploratory]*.)
- **Caveat:** the library IRs are sine-sweep measurements with an unknown sweep (`NEXT_STEPS.md` §9). A time-variant unit's sweep IR depends on the sweep, so depth and rate can't be read off these numbers directly. The existence of the variation and its program dependence don't depend on the sweep.
- **The tail predictor sees the same thing within one file** (§3.2): a flat, signal-proportional unpredictable part at about −18 dB (L) and −18 to −22 dB (M).
- **Reading (H):** the "stochastic control" of the 1986 review (`algorithm_research.md` §7.2b) is real, **strong in M and L and not seen in S** (S's pair difference is timing-type and stays at −46 to −50 dB in the 125–500 Hz octaves at 0.6–1 s, so any flat gain-type variation in S is below about −45 dB there). Its flat spectrum points to randomised gains or switched tap positions (jumps of many samples), not small delay modulation, and it runs free (not restarted per capture). The review's "Grad der Stochastik differs per program" fits, though the review also calls S the densest. Rate, depth and which taps vary are open; the library can bound them (§6).

### 4.1 S repeat captures with the capture clock fitted

`--stage clock`: the 40 S AMPL pairs (90 vs 60 %) in the capture domain (raw 96 kHz, 100 Hz–10 kHz), aligned by cross-correlation, then either a sub-sample shift only or a shift plus a clock ratio, fitted on R over 0.05–0.6 s and applied to both channels. One extra parameter per capture.

| Window (s) | Shift only, L / R | Shift + clock ratio, L / R | R, 10–90 % |
|---|---|---|---|
| 0.05–0.15 | −51.0 / −50.4 | **−57.6 / −56.2** | −58.4 … −49.3 |
| 0.15–0.3 | −49.5 / −49.7 | **−59.1 / −58.5** | −61.9 … −51.3 |
| 0.3–0.6 | −40.9 / −41.2 | **−58.2 / −57.8** | −60.1 … −52.9 |
| 0.6–1.0 | −33.9 / −34.0 | **−51.8 / −51.7** | −54.4 … −36.6 |

- **Clock ratio between the two captures of a pair:** −1.7 to +4.3 ppm (median −1.4), and the same within each LF group: LF 0.1 −1.6…−1.3, LF 0.4 −1.7…−1.5, LF 1.0 +2.0…+4.3, LF 2.0 −0.6…−0.2. That looks like one clock per capture session (**H**). The fit window is in-sample for the first three windows; 0.6–1.0 s is outside it.
- **So the library repeats itself to about −57 dB** (S, until the fade takes over), as long as each capture gets its own clock ratio. This is the ceiling for any model nulled against one Tomkin capture: −50 dB is reachable and checkable, −60 dB can't be demonstrated with this library. What limits it below −57 dB is not resolved (residual timing, sweep artefacts or a tiny time variance); it is not noise (§4.2).
- These numbers supersede the device-domain S pair nulls in the table above for S (which have no clock fit). M's and L's decorrelation is flat across frequency, so a clock can't explain it.

### 4.2 Noise floor: library vs the CLA capture

`--stage floors`: in-band (20 Hz–11.3 kHz) level per 25 ms frame, relative to each file's 50–300 ms reverb level, capture domain. "Natural decay holds to" is how far the tail keeps its own decay rate (local slope 0.7–1.4 × the rate measured between 0.3 s and −30 dB) before it flattens (a floor) or steepens (a fade).

| | CLA (eugtone) L / R | Library S L / R (512 files) | Library M (519) | Library L (360) |
|---|---|---|---|---|
| Natural decay holds to | −65.2 / −71.6 dB, then **flattens** | median −38.9; 5 % −85; deepest −91 / −90 | median −34; 5 % −77 / −78; deepest −91 | median −31; deepest −40 (all faded early) |
| Floor (plateau) | **−71.5 / −77.5 dB** re reverb, **−105.2 / −110.8 dBFS** (slope before it 0.06–0.07 × the decay rate: flat) | 8 files level off, at **−88.3 / −88.8 dB**, **−148 dBFS** | 26 / 20 files level off at −80 / −83 dB, −140 / −142 dBFS | none visible |
| How the quietest frame is reached | floor | fade (slope 2.6–2.7 × the decay rate) | fade (2.8–2.9 ×) | fade (4.6–4.7 ×) |

- **The floors are not consistent.** CLA has a real noise floor, L 6 dB above R. The library shows no capture noise at all: S tails keep their natural decay to about −85…−91 dB re the reverb and stop at −148 dBFS, the 24-bit export floor (rounding noise about −149 dBFS full band). L and R are equal. So the library is at least 11–17 dB cleaner relative to the signal than CLA. A longer or louder sweep (more processing gain) or denoising would both do this; with the sweep unknown (`NEXT_STEPS.md` §9) the two can't be told apart.
- M's plateaus (−80 dB, −140 dBFS) are 8 dB above S's and may be where the time-variant tail meets the format floor; they are not separable as noise.
- **Tomkin's fades usually come first:** the median file is at its natural decay only to −31…−39 dB re reverb, and every Long RT file is faded by about −40 dB. In most library files the usable null window ends where the fade starts, not at a noise floor.
- **Pre-roll is not noise in either source:** it sits 31–48 dB below the reverb, far above both tail floors *[exploratory]*.
- **CLA's noise-limited null per window** (from `model_s_results.md` §9.3, in band, L / R): −67 / −76 dB at 50–300 ms, −54 / −63 at 300–400 ms, −46 / −54 at 400–600 ms, −31 / −39 at 0.6–1 s. So CLA can show −60 dB or better only in roughly the first 300 ms, and mainly on R.

## 5. What this means

- **Q3 is answered for structure:** M and L have their own cluster scale (×2, ×4), onset scale and a different late recursion; the harness's unknown loop sets are assigned (§3.1). §7 explains the recursion as S's loops, scaled and coupled.
- **The null target needs a decision for M and L.** Two captures of the same M/L setting from the same unit null each other only to −9 to −14 dB at 0.3–0.6 s and −3 to −7 dB after 1 s. No deterministic model can beat that against one capture. For M and L the waveform target can only apply to the deterministic part (cluster and first recirculation, about the first 0.1–0.2 s); beyond it, the right target is statistical (decay, spectrum, decorrelation rate between captures). S is unaffected: its repeat captures agree to about −57 dB with the clock fitted (§4.1).
- **What can be demonstrated for S:** against Tomkin, about −57 dB at best (repeat-capture ceiling, one clock parameter per capture), in held-out windows up to each file's fade; the noise floor is not the limit there (§4.2). Against CLA, −60 dB or better only in roughly the first 300 ms, and it is a different unit (S&H pole, L/R offset).
- **Budget:** M and L use the same four lines as S (scaled, 11.5 K / 23 K words) and the same 20 + 20 read taps; the coupling adds a few multiply-adds per sample. The randomisation must also fit the budget.

## 6. Next tests

*Done since: the S carry-over (§7) and the line laws (§8). With the deterministic model complete, items 1–2 can now be tested inside it (which matrix entries, line gains or taps vary).*

1. *(Done, §9: per second, 0.83 Hz, no dependence on RT, LF or HF.)* **Decorrelation rate vs settings:** does the pair-null decay rate depend on RT, HF or LF? Does it scale per recirculation (per N samples) or per second? This separates "random change per pass" from a free-running clock.
2. *(Done, §9: the coupling's rotation angle.)* **Which part varies:** compare the pair difference with the five-tap recursion: is it consistent with randomised tap gains (β, γ), switched tap positions, or a varying read-out?
3. ~~**Loop gains and filters from the recursion** at RT/HF/LF steps (Q4–Q7 for M/L), using the deterministic part.~~ Done (§8): S's laws, length-scaled.
4. **Harness:** per-program cluster templates (§1) with a 60 ms window, and the M/L lag sets as known loop sets.
5. **Non-Lin:** its late part (`1108-1254-1386` and `231-628-…` in the harness) through the same predictor.

## 7. The S carry-over test (October 2026)

`.venv/bin/python tools/ml_carryover.py [--stage diffuser|taps|routing|firstpass]` → `analysis/ml_carryover/` (`*.json`, `*.log`). `NEXT_STEPS.md` §5 item 5: do S's diffuser, read taps and R routing, scaled ×2 / ×4, explain M's and L's injection and first recirculation?

**Method.** Per AMPL group, a cluster-free reverb (y_lo − A_lo/(A_hi − A_lo)·(y_hi − y_lo), the low-AMPL capture aligned onto the high one, cluster tap 0 at a fixed index). The diffuser is all-pass, so correlating the first pass with the diffuser's impulse response (a matched filter) turns it back into the read taps as isolated spikes when the diffuser is right, and smears them when it is not. S is the control: with its own diffuser (AP20 · AP72 · AP230, −1/√2) the matched filter puts 90 / 91 % of the first-pass energy into 20 spikes per channel, exactly at the Q2 taps (`model_s_results.md` §16.2).

> **Summary.** Yes, with three exceptions in the detail and one structural change:
> - **The diffuser carries over** with exact, slightly irregular lengths: **M AP42 · AP144 · AP460, L AP82 · AP286 · AP918**, all at **−1/√2**. (S × 2 would be 40 / 144 / 460; L = 2 × M − 2 for every stage.) Plain ×2 / ×4 smears the taps (24 / 9 % in 20 spikes against 91–92 %). **C**
> - **The read taps carry over:** 20 equal, sign-only taps per channel in every M and L group (160 and 180 groups), with **S's sign pattern tap for tap**, at **2 × S − 8…−12 (M) and 4 × S − 24…−35 (L)**, integers. R = L through **51 / 661 (M) and 101 / 1320 (L)**, the same binary 26 / 332 routing as S, per tap. **C**
> - **The damping step carries over:** the last all-pass takes S's AP230 two-tap damping with **the same four step values** (0 / 0.17 / 0.33 / 0.50) at short RT × low HF. **S**
> - **New in M and L: an output smoothing kernel** S does not have, on the cluster and the reverb alike: symmetric, about **[0.024, 0.952, 0.028] (M)** and **[0.11, 0.78, 0.12] (L)**, setting-independent (−0.8 / −4.2 dB at 10 kHz). **S**
> - **The structural change: M's and L's four loops are coupled.** The lines are S's ×2 / ×4 (**M 2641 / 2775 / 2919 / 3171, L 5279 / 5547 / 5837 / 6339**, integers), but each tap recirculates through **two** lines at ±g/√2, through an orthogonal 4 × 4 matrix that couples loops 1/3 with lines 2/4 and loops 2/4 with lines 1/3. Its characteristic polynomial is exactly the five-term recursion of §3 (lags, and signs + − − + −). S's matrix is diagonal. **C**
> - **Held out** (positions, signs, kernel from 6 flat long-RT settings; per capture one shift, one cluster gain, one tap gain; b per setting): first pass **S −42.8 / −41.5, M −40.7 / −39.2, L −36.6 / −35.9 dB**, against repeat captures **−48.2 / −54.6, −47.1 / −48.8, −41.3 / −50.8 dB (§7.6)**. Cluster and gap **S −49.1 / −44.1, M −46.6 / −42.0, L −43.0 / −38.8 dB** (at the repeat-capture level).

### 7.1 Diffuser

Matched-filter sparsity (share of first-pass energy in the 20 largest isolated spikes per channel; median of 13 settings per program, a ±3-sample grid per stage on 4 of them):

| Program | Best stage lengths (all 4 grid settings) | Share L / R | Plain S × scale | Runner-up (grid) |
|---|---|---|---|---|
| S | 20 / 72 / 230 | 0.90 / 0.91 | (same) | 0.54 |
| M | **42 / 144 / 460** | 0.91 / 0.91 | 0.25 / 0.23 | 0.54 |
| L | **82 / 286 / 918** | 0.92 / 0.92 | 0.10 / 0.08 | 0.58 |

- One sample off in any stage costs a third or more of the share (runner-up 0.54–0.58), so the lengths are exact. Where a grid setting chose the last stage +1 (one per program), it is a damped setting (b > 0, §7.4), whose two-tap filter adds delay.
- **All-pass gain:** the share peaks at **−0.70…−0.71** in all three programs (−0.707 best of the grid −0.80…−0.62). The −1/√2 of S (§16.8 there) carries over.
- **Lengths:** M = 2 × S except the first stage (42, not 40); L = 2 × M − 2 for all three stages (82, 286, 918). The cluster's R-only taps show the same kind of "×2 minus a little" (§1), so the program tables are probably computed with an offset rather than stored as clean multiples (**H**).

### 7.2 Read taps and R routing

Greedy pursuit (band-limited pulses shaped by the program's output kernel, §7.5) on the matched-filter output of every AMPL group, matched to S's tap table:

| | S (128 groups) | M (160) | L (180) |
|---|---|---|---|
| All 20 taps found per channel | 100 % | 100 % | 100 % |
| Sign agreement with S's pattern | 1.000 | 1.000 | 1.000 |
| Magnitude spread over the 20 taps | 0.01 dB | 0.01–0.02 dB | 0.01 dB |
| Position IQR over groups (max) | 0.15 | 0.13 | 0.004 |
| Integer part re cluster tap 0 | the Q2 table (+0.00…0.03) | +0.01…0.06 | +0.03…0.12 |

The fractional parts grow with the lag (L: +0.03 at 952 to +0.12 at 7379): the capture clock against the nominal device rate, not the unit. Tap positions (L channel; R = L + route):

| S | 244 | 276 | 338 | 624 | 677 | 768 | 797 | 826 | 913 | 971 | 984 | 1069 | 1129 | 1256 | 1319 | 1368 | 1449 | 1514 | 1521 | 1780 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sign | + | − | + | − | − | − | − | − | − | + | + | + | − | − | + | + | − | + | + | − |
| M | 480 | 543 | 667 | 1239 | 1344 | 1526 | 1584 | 1644 | 1816 | 1934 | 1960 | 2128 | 2249 | 2504 | 2629 | 2728 | 2889 | 3020 | 3033 | 3550 |
| M − 2S | −8 | −9 | −9 | −9 | −10 | −10 | −10 | −8 | −10 | −8 | −8 | −10 | −9 | −8 | −9 | −8 | −9 | −8 | −9 | −10 |
| L | 952 | 1078 | 1326 | 2469 | 2679 | 3043 | 3159 | 3280 | 3623 | 3860 | 3912 | 4247 | 4489 | 5000 | 5249 | 5448 | 5769 | 6030 | 6059 | 7090 |
| L − 4S | −24 | −26 | −26 | −27 | −29 | −29 | −29 | −24 | −29 | −24 | −24 | −29 | −27 | −24 | −27 | −24 | −27 | −26 | −25 | −30 |

- **R routing:** each R tap is its L tap delayed by **S 26 / 332, M 51 / 661, L 101 / 1320**, the same choice per tap in all three programs (σ ≤ 0.01 sample). M = 2 × S − 1 / −3, L = 2 × M − 1 / −2. R taps sit at 2 × S − 9…−12 (M) and 4 × S − 29…−35 (L).
- The offsets from 2 × S / 4 × S take three or four values per program, not one, and don't follow the loop a tap belongs to. They are exact integers, so the tap table is stored or computed per program with its own rounding (**H**); the rule is open.
- The onsets of §2 (474 / 1,134 and 948 / 2,268) are these first taps (480 / 1,141 and 952 / 2,272) seen through a threshold.

### 7.3 Routing through the lines: the four loops are coupled

For each tap, the matched-filter value at tap + Pⱼ for each line (flattest HF settings: S 32, M 32, L 36 groups):

- **Line lengths are exact integers:** the recirculated taps arrive at tap + Pⱼ + 0.00…0.03 (S), +0.05…0.06 (M), +0.07…0.08 (L): the capture clock again (§7.2), no fractional delay or loop-filter group delay. So **M's lines are 2641 / 2775 / 2919 / 3171** (2 × S − 1, −3, −3, −3) and **L's 5279 / 5547 / 5837 / 6339** (2 × M − 3, −3, −1, −3). **C**
- **Routing matrix** (line ← the tap's S loop; median gain re the tap):

  | | S: loop 1 | 2 | 3 | 4 | M: loop 1 | 2 | 3 | 4 | L: loop 1 | 2 | 3 | 4 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | line 1 | **+0.72** | 0 | 0 | 0 | 0 | **+0.49** | 0 | **+0.53** | 0 | **+0.55** | 0 | **+0.60** |
  | line 2 | 0 | **+0.72** | 0 | 0 | **+0.50** | 0 | **−0.52** | 0 | **+0.55** | 0 | **−0.58** | 0 |
  | line 3 | 0 | 0 | **+0.70** | 0 | 0 | **+0.52** | 0 | **−0.48** | 0 | **+0.58** | 0 | **−0.54** |
  | line 4 | 0 | 0 | 0 | **+0.67** | **−0.51** | 0 | **−0.46** | 0 | **−0.57** | 0 | **−0.53** | 0 |

  Off-pattern entries are ≤ 0.011 in magnitude. Gains are medians over RT, so only their pattern and ratios matter here.
- **M and L share one matrix:** F = (g/√2) · [[0, 1, 0, 1], [1, 0, −1, 0], [0, 1, 0, −1], [−1, 0, −1, 0]] (rows: lines; columns: loops), times each line's own gain and filter. F/g is **orthogonal**, so the coupling itself is lossless. It is bipartite: loops 1 and 3 feed lines 2 and 4 through [[1, −1], [−1, −1]]/√2, and loops 2 and 4 feed lines 1 and 3 through [[1, 1], [1, −1]]/√2. **C** (pattern; σ of the sign over groups 0 where no other tap collides).
- **This is the five-term recursion of §3.** With xᵢ = g·z^−Pᵢ (one gain for all lines, filters left out), det(I − X·F) = 1 − (x₁x₂ − x₂x₃ − x₁x₄ + x₃x₄)/2 + x₁x₂x₃x₄. Its lags are P₁ + P₂ = 5416, P₂ + P₃ = 5694, P₁ + P₄ = 5812, P₃ + P₄ = 6090 and P₁ + P₂ + P₃ + P₄ = 11506 (L: 10826 / 11384 / 11618 / 12176 / 23002), with signs **+ − − + −**, exactly the lags and signs the tail predictor found blind (§3.2). The "2N line" reading of §3.2 is superseded. The tail ACF peaks at 2641 (M) and 5279 (L) are line 1 (§3.1).
- **§3.3's −0.8-sample offset is not a property of the lines:** single recirculations land on the integers. Why the predictor's lag FIRs peak 0.8 sample early is open (a predictor artefact with the time-variant part, **H**).

### 7.4 The damping step

b of the last all-pass (two-tap loop filter (1 − b) + b z⁻¹, as S's AP230), fitted per setting on the first pass of the high-AMPL capture:

| S: RT \ HF | 0.2 | 0.4 | 0.6 | 1.0 |
|---|---|---|---|---|
| 0.3 | 0.50 | 0.50 | 0.33 | 0 |
| 0.6 | 0.50 | 0.33 | 0.17 | 0 |
| 0.8 | 0.33 | 0.17 | 0 | 0 |
| 1.0 | 0.17 | 0 | 0 | 0 |
| ≥ 1.4 | 0 | 0 | 0 | 0 |

| M: RT \ HF | 0.1 | 0.2 | 0.4 | 0.7 | 0.9 |
|---|---|---|---|---|---|
| 0.6 | 0.50 | 0.50 | 0.50 | 0.33 | 0.17 |
| 0.8 | 0.50 | 0.50 | 0.33 | 0.17 | 0 |
| 1.2 | 0.50 | 0.33 | 0.17 | 0 | 0 |
| 1.8 | 0.17 | 0.17 | 0 | 0 | 0 |
| ≥ 2.6 | 0 | 0 | 0 | 0 | 0 |

| L: RT \ HF | 0.1 | 0.3 | 0.5 | 0.8 | 1.0 |
|---|---|---|---|---|---|
| 2.2 | 0.50 | 0.50 | 0.50 | 0.17 | 0 |
| 2.6 | 0.50 | 0.50 | 0.33 | 0 | 0 |
| 3.6 | 0.33 | 0.17 | 0 | 0 | 0 |
| 4.5 | 0.17 | 0 | 0 | 0 | 0 |
| ≥ 6.0 | 0 | 0 | 0 | 0 | 0 |

Medians over LF (the four LF values agree in every cell). Fitted values: 0.500–0.502, 0.330–0.332, 0.171–0.173 and ≤ 0.004.

- **S's fit reproduces its Q2 table exactly** (`model_s_results.md` §16.9), so the bench measures the same thing.
- **M and L use S's four step values** and the same kind of staircase: b rises as RT and HF fall. **S**
- The cells don't map onto S's by a plain RT rescaling (RT / 2, RT / 4) or by RT × HF; the rule that picks the step is open for all three programs (`model_s_results.md` §14.3, §16.17).
- Fitting b per setting matters only where it is on: held out it changes the median first-pass null by 0.3–0.9 dB, but the 90th percentile by 30 dB (S −10.6 → −42.1 dB on L, M −9.1 → −40.4, L −10.5 → −35.5).

### 7.5 The output kernel (M and L only)

A least-squares FIR (−5…5) on the cluster window of the high-AMPL captures (the cluster is 8 known taps, so the kernel is all that is fitted), unit DC, median of the 6 training settings × 2 channels:

| | −2 | −1 | 0 | +1 | +2 | spread (max σ) | at 10 kHz |
|---|---|---|---|---|---|---|---|
| S | −0.0000 | +0.0007 | 1.0001 | +0.0022 | −0.0022 | 0.003 | 0 dB |
| M | −0.0093 | +0.0235 | 0.9518 | +0.0276 | +0.0063 | 0.005 | −0.8 dB |
| L | −0.0392 | +0.1134 | 0.7834 | +0.1217 | +0.0173 | 0.004 | −4.2 dB |

- Symmetric (zero phase) and the same in every setting; it acts on the reverb as well as the cluster (without it, M's first pass is held at about −30 dB, §7.6). In band, L's is close to [1/8, 3/4, 1/8]; the alternating outer taps are out of band (the fit is in 100 Hz–10 kHz).
- **Not in S** (S's kernel is δ within 0.003). So it is program-specific processing, not the analogue chain or the library's deconvolution.
- **Cause (H):** a ±1-sample jitter of a common read or write address, averaged by the sweep measurement, gives exactly a symmetric [p/2, 1 − p, p/2] (p ≈ 0.05 in M, 0.23 in L). That would tie it to the time variance of §4 (stronger in L). A deliberate output low-pass is the alternative. The library alone may not separate them: a jitter averaged over a 10 s sweep looks the same in every capture.

### 7.6 Held-out first pass

Positions, signs and kernel from 6 training settings per program (the flattest HF, the longest RTs); held out on every other group. Per capture: one sub-sample shift, one cluster gain, one tap gain (both windows, per channel). Per setting: b (§7.4). 100 Hz–10 kHz, device domain, no clock fit (the windows end by 0.07 / 0.14 / 0.29 s, where a few ppm don't matter). Repeat-capture null: another capture of the same group as the model, with the same freedoms (shift, cluster gain, scale). Medians, L / R [repeat L / R]:

| Program (captures held out) | Cluster and gap, [0, onset) | First pass, [onset, onset + shortest line) | Excess R, first pass | First pass without the kernel |
|---|---|---|---|---|
| S (366 / 488) | **−49.1 / −44.1** [−47.2 / −45.3] | **−42.8 / −41.5** [−48.2 / −54.6] | +13.7 | −42.5 / −41.6 |
| M (461) | **−46.6 / −42.0** [−46.4 / −44.6] | **−40.7 / −39.2** [−47.1 / −48.8] | +9.6 | −30.3 / −30.8 |
| L (174 / 348) | **−43.0 / −38.8** [−41.3 / −40.6] | **−36.6 / −35.9** [−41.3 / −50.8] | +15.3 | −17.2 / −17.4 |

- **Held out equals training** within 0.5 dB in every row (training: S −42.6 / −40.7, M −40.7 / −39.1, L −36.8 / −35.8), and the 90th percentile is within 0.3–1.1 dB of the median: the model is the same for every setting, HF, LF and RT (per-HF medians vary by ≤ 0.3 dB).
- **The carry-over is complete at S's own level.** M's first pass is predicted as well as S's own compact model predicts S (−40.7 / −39.2 against −42.8 / −41.5; excess over repeat captures +9.6 dB against +13.7). L is 4–6 dB behind (+15.3 on R).
- **Cluster and gap are at the repeat-capture level** in all three (excess +0.9 / +1.5 / +1.7 dB on R). This is the window the null target puts at −50 dB for M and L; in this bench the repeat captures themselves reach only −41…−47 dB there (device domain, no clock fit, 100 Hz–10 kHz), so −50 dB can't be shown with it. A capture-domain comparison with fitted clocks (as `model_s_results.md` §17) is the next step for that window.
- **What the first pass still misses** (the 10–15 dB to the repeat level) is the same gap S's compact model has (`model_s_results.md` §16.11, §17: a dense remainder, mostly out of band). L's extra 4–6 dB is probably the kernel: it is fitted on the cluster window only, and L's kernel is the strongest (**H**).
- **Parameters:** per program the diffuser (3 lengths, one gain), the 40 tap positions and signs, the 11-tap kernel and the step table; per capture 3 (shift, two gains). The kernel and taps are trained on 6 settings, never on the held-out ones.
- **The first recirculation** (the coupling) is measured in §7.3 but not yet nulled: its gains are the line gains and filters per setting (Q4–Q7 for M/L), and from there on the time variance of §4 caps any waveform null.

## 8. Line gains and filters: S's laws, length-scaled (October 2026)

`.venv/bin/python tools/ml_lines.py [--stage law|report|fit]` → `analysis/ml_lines/` (`law/` one file per AMPL group, `report.log`, `fit.json` / `fit.log`). `NEXT_STEPS.md` §5 item 5, the step after §7: what are the line gains and filters per setting (Q4–Q7 for M/L)?

**Method.** The §7 model is closed through the four lines: y = h ∗ [g_c C + g_t D · 1ᵀ(I − ΛF)⁻¹T], Λ = diag(G_j(z) z^−P_j), F the §7.3 coupling matrix, T the read taps grouped by their S loop.
- **Line responses:** each G_j is **S's control laws** (`model_s_results.md` §14.2, the 13 numbers) evaluated at line j's own length: Jot gain, LF shelf with the fixed pole, HF one-pole. S's three limits apply: HF pole ≤ 0.84, per-pass DC gain ≥ 0.18, β = 0.019 (§16.15–16.16 there). M's and L's LF steps that S lacks (0.7, 0.8, 1.2) are interpolated in log LF. **No loop number is fitted on M or L**, so everything after the first pass is a held-out prediction.
- **Per capture** (fitted only on [0, onset + shortest line), which no line reaches): one sub-sample shift, one clock ratio and two gains per channel (cluster, reverb). The fitted clocks are tight: M +12.8 ppm (10–90 % +12.2…+13.2), L +13.4 (+12.4…+13.5).
- **Windows** (100 Hz–10 kHz, up to each capture's fade): cluster and gap, first pass, then three recirculation windows one shortest line long (M R channel: 0.14–0.24 / 0.24–0.34 / 0.34–0.44 s; L: 0.29–0.49 / 0.49–0.69 / 0.69–0.89 s).
- **Repeat-capture null:** another capture of the same AMPL group as the model, with the same freedoms (shift, clock, cluster column, gain).
- **Rendering:** frequency domain on an exponentially weighted grid (no time aliasing at RT 20 s) with integer delays; clock and shift are applied after, as a Kaiser-sinc warp. The tap positions' fractional parts are the capture clock (§7.2), so the taps render as integers.
- **Decay per octave:** 10 log10 of the band energy ratio measured / model per 39 ms frame, from the onset to the fade. Its slope is the decay error, and it is energy-based, so the time variance does not bias it the way it biases a waveform null.

> **Summary.**
> - **S's laws, length-scaled, are M's and L's line laws.** With no parameter fitted on M or L, the model predicts the recirculations **better than a repeat capture does**: recirc1 (R) M −20.8 vs repeat −19.8 dB, L −18.9 vs −17.5. The model beats the repeat in 69 % (M) and 65 % (L) of captures, by a median 2.2 dB. **S**
> - **Coupling confirmed:** with the lines uncoupled (F = I) every recirculation window is at about 0 dB (−1.2…+3.2 dB). **C**
> - **S's coefficient limits carry over unchanged:** the HF-pole limit scan peaks at **0.84** in both programs (recirc1 M −26.9, L −23.6 dB at 0.84; 2–3 dB worse at 0.80 / 0.88; without the limit −16 / −11). The DC floor is best at 0.18 (shallow: ±0.3 dB over 0.14–0.22), and β = 0.019 is at or near the best. **S** (HF), **S/weak** (DC, β)
> - **Decay matches to under 1 %:** measured RT / law RT over all octaves is 0.996 (M) and 0.995 (L), IQR 0.985–1.001. Fitted per setting, M's and L's LF law numbers reproduce S's (interpolated) ones to 0.1–1.5 %, except L at LF 2.0 (DC RT 3 % short) and LF 0.1 (capped settings). **S**
> - **The model is at the time-variance limit:** splitting each window's null into capture-to-capture variance v and the model's own error e (repeat = 2v, model = v + e), e is 2–9 dB below v in every recirculation window. Past the first pass, waveform nulls against one capture can't show more (§4, `NEXT_STEPS.md` §2: accuracy over null). **S** (given the independence assumption)

### 8.1 Nulls against repeat captures

Median over captures, L / R [repeat L / R], dB; v / e for R; held out in time and setting (loops from S's laws):

| Window | M (160 groups, 479 captures) | v / e | L (180 groups, 360 captures) | v / e |
|---|---|---|---|---|
| Cluster and gap | −48.5 / −45.8 [−49.9 / −47.3] | −50.3 / −47.6 | −44.9 / −45.5 [−44.8 / −45.5] | −48.5 / −48.5 |
| First pass | −44.5 / −44.1 [−48.6 / −49.2] | −52.2 / −44.8 | −44.3 / −42.7 [−46.2 / −54.1] | −57.1 / −42.9 |
| Recirculation 1 | **−22.0 / −20.8** [−21.2 / −19.8] | −22.8 / −25.1 | **−19.7 / −18.9** [−18.4 / −17.5] | −20.5 / −23.8 |
| Recirculation 2 | −16.5 / −16.0 [−14.4 / −14.1] | −17.1 / −22.2 | −16.1 / −15.7 [−13.7 / −13.4] | −16.5 / −23.6 |
| Recirculation 3 | −13.1 / −12.6 [−10.9 / −10.2] | −13.3 / −20.8 | −12.8 / −12.6 [−10.4 / −10.1] | −13.1 / −22.3 |
| Recirc. 1, lines uncoupled (R) | −1.2 | | +0.7 | |
| Recirc. 1, no limits (R) | −19.3 | | −18.5 | |

- **The model beats the repeat captures past the first pass, at every RT, LF and HF step** (by-control tables in `report.log`; the only exceptions are three cells within 0.3 dB: M RT 4.0, L RT 2.6, L LF 0.1). That is what a deterministic model of the common response should do when each capture adds its own random part.
- **The first pass is where the model's own error dominates** (e 7–14 dB above v): the same dense remainder S's compact model has (`model_s_results.md` §16.11, §17), about −44 dB. It enters every later window too, but at −44 dB that is far below the time variance there.
- **Cluster and gap with fitted clocks:** M at the repeat-capture level (−48.5 / −45.8 against −49.9 / −47.3), L exactly at it. With fitted clocks the repeat captures reach −50 / −47 dB (M) and −45 dB (L) in this window, so **−50 dB is not demonstrable here against this library**; the model is within 1.5 dB of what is.
- **Recirculation nulls get worse with RT** (M recirc1 −26 dB at RT 0.6, −19 at RT 5.0), and so do the repeat nulls: longer RT keeps more of the earlier, already-decorrelated passes in each window.

### 8.2 The coefficient limits

Recirculation nulls (median over captures and channels, recirc1 / recirc2 / recirc3, dB) in the groups where each limit can act; everything else at the law:

| HF pole limit | M (41 groups) | L (21 groups) |
|---|---|---|
| 0.80 | −22.4 / −11.6 / −7.8 | −18.7 / −10.7 / −7.1 |
| 0.82 | −24.8 / −14.6 / −10.2 | −21.1 / −13.6 / −9.8 |
| **0.84** | **−26.9 / −17.1 / −12.5** | **−23.6 / −16.5 / −12.4** |
| 0.86 | −26.1 / −16.8 / −12.3 | −23.4 / −16.3 / −12.4 |
| 0.88 | −24.0 / −13.8 / −9.8 | −21.3 / −12.9 / −9.4 |
| 0.90 | −21.9 / −10.5 / −6.9 | −18.3 / −9.4 / −6.2 |
| none | −16.2 / −3.0 / −0.7 | −11.4 / −2.5 / −0.5 |

| DC floor | M (45) | L (30) | | β | M (45) | L (30) |
|---|---|---|---|---|---|---|
| 0.10 | −22.0 / −15.7 | −19.5 / −15.7 | | 0 | −21.3 / −15.6 | −19.3 / −15.6 |
| 0.14 | −22.2 / −16.1 | −19.6 / −16.0 | | **0.019** | **−22.3 / −16.4** | **−19.7 / −16.2** |
| **0.18** | **−22.3 / −16.4** | **−19.7 / −16.2** | | 0.04 | −22.2 / −16.4 | −19.7 / −16.2 |
| 0.22 | −22.4 / −16.3 | −19.7 / −16.2 | | | | |
| 0.26 | −22.2 / −15.9 | −19.5 / −15.9 | | | | |
| none | −20.3 / −15.3 | −18.9 / −14.9 | | | | |

- **The HF limit is sharp and the same as S's** (0.84, between 0.84 and 0.86 in both programs; S's refit gave 0.835–0.845). It acts per line, so it bites at longer RT in M and L than in S (more damping per pass for the same RT): 41 M groups (HF 0.1–0.2) and 21 L groups (HF 0.1). **S**
- **The DC floor matters (0.8–2 dB against none) but its value is weakly fixed** by these windows: 0.14–0.22 are within 0.3 dB, 0.18 is in the middle. Its effect is mostly below 100 Hz, outside the null band. **S** that a floor exists; its value carries over from S as the best estimate.
- **β = 0.019 is at the best** (0 is 0.4–1.0 dB worse; 0.04 equal). Consistent with S's value, not an independent measurement of it.
- The decay-based per-setting fit (§8.3) cannot fix the HF pole at HF 0.1: the top octaves decay before 0.15 s, and the pole then trades off against the broadband RT (fitted 0.86–0.87 with an equal RT curve). The waveform scan above decides it.

### 8.3 Decay per octave, and the LF law per step

Measured RT / law RT (median over captures and channels):

| | 125 | 250 | 500 | 1 k | 2 k | 4 k | 8 k Hz |
|---|---|---|---|---|---|---|---|
| M, LF 0.1 / 0.7 / 1.2 / 2.0 | 0.993 / 0.996 / 0.988 / 0.991 | 0.997 / 0.996 / 0.992 / 0.996 | 0.987 / 0.997 / 0.995 / 1.001 | 0.980 / 0.994 / 0.993 / 0.997 | 0.986 / 0.996 / 0.997 / 0.998 | 0.993 / 0.998 / 0.997 / 0.998 | 0.997 / 0.999 / 0.999 / 1.000 |
| L, LF 0.1 / 0.8 / 1.2 / 2.0 | 0.996 / 0.982 / 0.972 / **0.958** | 0.996 / 0.992 / 0.987 / 0.982 | 0.988 / 0.991 / 0.990 / 0.991 | 0.986 / 0.994 / 0.993 / 0.995 | 0.988 / 0.996 / 0.995 / 0.997 | 0.994 / 0.997 / 0.997 / 0.998 | 0.999 / 1.000 / 0.999 / 1.000 |

Band levels at the onset agree with the model within ±0.15 dB in every octave and HF step (`report.log`), so the spectrum is right, not only the decay.

**The law's numbers per LF step, fitted per setting** (`fit` stage: Nyquist-end RT / label and DC RT / Nyquist-end RT, from the band decay, settings at the top HF step where no limit is active; median over RT):

| LF | M: hi | M: dc / hi | L: hi | L: dc / hi | S law (interpolated where S has no step): hi, dc / hi |
|---|---|---|---|---|---|
| 0.1 | 1.281 | 0.103 | 1.321 | 0.098 | 1.278, 0.103 |
| 0.7 / 0.8 | 1.032 | 0.684 | 1.021 | 0.781 | 1.035, 0.684 / 1.022, 0.788 |
| 1.2 | 0.986 | 1.197 | 0.986 | 1.187 | 0.989, 1.205 |
| 2.0 | 0.958 | 2.019 | 0.958 | 1.967 | 0.958, 2.032 |

- **M's LF steps are S's law exactly**, including the two steps S doesn't have (0.7, 1.2), from log-LF interpolation of S's four steps. So the LF control is one continuous law shared by the programs, not a per-program table (**S**). L agrees as well except at the low end of LF 2.0 / 1.2 (DC RT 3 % / 1.5 % short, the 125 Hz shortfall above) and LF 0.1 (only 4 of 9 settings uncapped).
- **A small common shortfall:** measured RT is 0.2–2 % below the law in the mid bands, most at HF 0.1 (2 %) and in L at low frequencies. Two readings, not separated yet: a real slightly stronger damping in M/L, or the sweep measurement of a time-variant system, which loses coherent energy as the response decorrelates, so the measured IR decays a little faster than the unit's energy (**H**). The second predicts a shortfall that grows with the decorrelation rate; testing it needs the time-variance model (Q3b).
- The HF steps behave the same way: fitted HF poles equal the law's at HF 0.2–0.9 in M and 0.3–1.0 in L (M 0.720 / 0.392 / 0.143 / 0.041 against 0.719 / 0.391 / 0.143 / 0.041).

### 8.4 What this answers, and what's next

- **Q4–Q7 for M/L (line gains and filters) are answered at the library's steps:** S's laws at the M/L line lengths with S's limits. The per-program control laws are one set: what changes between programs is the line lengths, the coupling and the kernel, not the control arithmetic. **Parameters added for M/L: none** (13 law numbers and 3 limits from S, the §7 structure trained on 6 settings, 6 per capture: shift, clock, and per channel a cluster gain and one reverb gain shared by all lines).
- **For the null target:** the reverb windows are now limited by the time variance, not the model (e below v). By the M/L rule in `NEXT_STEPS.md` §2 the next proof is statistical, the decorrelation curve, which needs a model of the variation itself.
- **Next (Q3b):** with the deterministic model in hand, the model-minus-capture residual in each recirculation window is the variation itself (plus e, 2–9 dB lower). `ml_lines.fitted(prog, key)` returns each capture of a setting with its model and usable end, from the stored clock and gains (no refit; reproduces the stored nulls exactly). Test which part varies (line gains, matrix signs, tap positions, the output kernel as averaged address jitter) by how that residual grows per pass, and check whether the 0.2–2 % RT shortfall follows it.

### 8.5 The law stage rerun (24 October 2026) **[S: the line laws stand on chain3 IRs, the M / L fades and S's pole curve and HF floor; H: L's shortest RT at LF 0.1]**

`EMT_KERNEL=chain3 tools/ml_lines.py --stage law --workers 4` (all 340 groups, 20 min, 4.0 GB peak), then `--stage fit` and `--stage report`; `analysis/ml_lines/{law,fit.json,fit.log,report.log}` are the new files. The old ones held the 0.72 s fade, the old IRs and the capped pole law. Nothing is fitted in the law stage. Median, R, new [old]:

| Window | M | L |
|---|---|---|
| cluster and gap | −53.5 [−45.8] | −52.6 [−45.5] |
| first pass | −48.0 [−44.1] | −45.6 [−42.7] |
| recirculation 1 / 2 / 3 | −20.7 / −16.0 / −12.6 [−20.8 / −16.0 / −12.6] | −18.9 / −15.8 / −12.6 [−18.9 / −15.7 / −12.6] |

- **Recirculations unchanged and still better than a repeat capture** (repeat: M −19.7 / −13.7 / −9.8, L −17.6 / −13.4 / −10.0); the model's own error is 5–7 dB under the capture-to-capture variance in every window. Cluster and gap are at the repeat level, so −50 dB can't be shown further there.
- **The HF 0.1 decay shortfall is gone** (measured RT / law at 500 Hz 0.980 → 0.999 in M, 0.977 → 0.995 in L; at 1 kHz 0.980 → 0.999 and 0.984 → 1.000): the pole curve and HF floor (§16). Level at the onset within 0.04 dB.
- **Low bands still short** (L LF 2.0 at 125 Hz 0.958 → 0.976; M LF 0.1 at 1 kHz 0.986): the sweep's under-reading of long LF decays (§11.2).
- **Per-LF fit unchanged:** Nyquist-end RT / label 1.280 / 1.033 / 0.986 / 0.958 (M), DC ratio 0.103 / 0.684 / 1.197 / 2.024; S's laws to 0.5 %. §8.3's per-LF numbers stand.
- **Sample changes:** the M fade is a 1.0 s ramp (§15), so 15 of M's 160 groups (short RT) have no usable window and are out of the fit (145); recirculations 2 and 3 have fewer captures (M 384 / 337 of 443 / 436). Several changes came together (IRs, fades, pole curve, floor); they are not separated here.
- **L's shortest RT** was the one setting group the law missed (recirculation 1 at −9 dB at HF 0.1 / 0.3); found and fixed in §8.6, which also updates the L table above (recirculation 1 / 2 / 3: −19.1 / −15.9 / −12.7 dB R, repeat −17.6 / −13.4 / −10.0).

### 8.6 L at RT 2.2: the HF control's floor is 0.40 there **[S: the floor at L's shortest RT step; H: why only that step]**

Recirculation nulls (R, median) at L RT 2.2 were −9 dB at HF 0.1 and −13 dB at HF 0.3 (all four LFs alike), against −20 dB at HF ≥ 0.5 and at RT 2.6 and above (−25 dB at RT 2.6 / HF 0.1). The first pass is fine (−50.8 dB), `b` is the same at HF 0.5 (which nulls well), and the stock pole-limit scan only reaches −11.4.

- **A per-setting refit does not repair it** (the `fit` stage's three numbers with no limits: recirculation 1 −9.6 → −10.7 dB, recirculations 2 / 3 stay at −5 / −3 dB), so it is not decay-rate coefficients. It shows only as measured decay 7–17 % longer than the law at 250 Hz – 8 kHz (fit / law 1.17 at 500 Hz).
- **Scanning the model's HF control** (floor off) at RT 2.2 puts the best HF at **0.40 for HF 0.1, 0.3 and 0.5 alike** (recirculation 1 / 2 −21.5 / −16.8 at HF 0.1, at the repeat level), a sharp optimum (0.35 and 0.45 are 3 dB worse). So HF acts as max(HF, 0.40) there, not the stock floor 0.12 s × 4 / RT = 0.218. Line poles at HF 0.40 are 0.55–0.61, where the stock cells sit at the 0.8475 ceiling.
- **Held out by LF:** the value chosen on LF 1.2 only (4 captures), scored on LF 0.1 / 0.8 / 2.0 (12 captures): recirculation 1 / 2 / 3 −11.4 / −6.0 / −3.5 → **−21.4 / −17.1 / −13.2 dB**, repeat −18.2 / −13.6 / −10.2. HF 0.5 – 1.0 cells unchanged to 0.1 dB. One number; it can't be held out by RT, since only RT 2.2 changes (RT 2.6 and 3.6 at HF 0.1 fit the stock floor: best HF ≤ the ceiling-pinned 0.2, and 0.25 is worse).
- **Adopted** in `ml_lines.line_coefs` (`L22_HF = 0.40`, `ML_L22=0` off); the 20 RT 2.2 law files, `fit.json` and `report.log` are rerun. L overall: recirculation 1 / 2 / 3 −19.9 / −16.3 / −12.9 (L), −19.1 / −15.9 / −12.7 (R), repeat −17.6 / −13.4 / −10.0; the RT 2.2 row −20.8 [−18.6], HF 0.1 decay ratios 0.993–1.002.
- **Altiverb LONG RT 2.0 preset (`tools/altiverb_floor.py`): uninformative about the floor.** The presets give no HF (inferred 0.4); scanning the effective HF with the floors off (L carry-over model, b from `model_s_results.md` §18.2, LF 1.2, floors off), the best recirculation-1 HF per preset is 0.5 at RT 2.0, 0.4 at 3.0 and 4.0, 0.45 at 5.0 – 8.0, 0.5 at 10.0 (nulls −14.5…−16.6 dB, broad dips; recirculation 2 −13…−14.4). So the Altiverb's L HF is about 0.45 at every RT, above any floor candidate (0.24 stock, 0.40 here) even at RT 2.0, and the floor cannot act there. What it does show: RT 2.0 sits inside the other presets' scatter (no jump, no extra anomaly at L's shortest RT on a second unit). A floor test needs an L setting with HF below 0.4 at RT ≤ 2.2; the presets have none.
- **Open (H):** why only this step. Two lines of thought, neither tested: (1) at S RT 0.3 (law 0.12 / 0.3) and M RT 0.6 (0.24 / 0.6) the floor is also 0.40, so each program's shortest RT step may share one HF limit, but L's RT 2.6 (step 2) shows no such floor, and L's RT axis is not 4 × S's; (2) a table in the unit indexed by the RT step. The Altiverb LONG RT 2.0 preset could not test it (above). The plugin must carry this entry when the M / L laws go in.

### 8.7 M / L on S's pole curve of `model_s_results.md` §32.6 (25 October 2026) **[S: neutral to slightly better on M / L's static metrics; adopted so that all three programs share one HF pole map]**

S's HF pole turned out to be one curve of the law's pole with a drop-back (`model_s_results.md` §32.4–32.6): 12 knot values, piecewise linear, flattening at ≈ 0.852, and S's HF 0.2 / RT 0.6 value 0.8465 past a law pole of 0.856. `ml_lines._pole_curve` now uses that map by default (the all-data fit, `analysis/lowband/chain3_pole_cells_allf.json`; `DROP_BACK` 0.8465). `ML_POLE_CURVE=25` gives §25.17's curve (law, ramp, ceiling 0.8475), `0` the capped law. Between law poles 0.62 and 0.83 the two curves agree to ≈ 1·10⁻³. They differ at the top: up to 0.8518 instead of the 0.8475 ceiling on lines with law poles 0.835–0.856, and 0.8465 instead of 0.8475 past the drop-back point. That is a change of up to 4.3·10⁻³ in 108 of M's 160 groups and 88 of L's 180. M / L's own best capped poles (§25.17, one shift shared by the four lines) were 0.8443–0.8517, consistent with either.

**Pulse check on the changed groups** (`pulse_check._ml_job`, nothing fitted on M / L), §25.17's curve → §32.6's, medians:
- M (108 groups): recirculation 1 / 2 / 3 −21.68 / −16.49 / −12.68 → −21.69 / −16.50 / −12.68 dB; octave decay |error| 2 / 4 / 8 kHz 0.56 / 0.49 / 0.39 → 0.53 / 0.49 / 0.34 %, late decay 2 kHz 1.20 → 1.04 %, 500 Hz 1.20 → 1.27 %.
- L (88 groups): octave decay |error| 250 Hz – 8 kHz 1.07 / 0.88 / 0.93 / 0.63 / 0.48 / 0.42 → 0.98 / 0.79 / 0.83 / 0.59 / 0.40 / 0.36 %; late decay 8 kHz 0.82 → 0.57 %, 125 Hz 7.18 → 7.58 % (Tomkin's fades dominate there, §14).
- All groups (report): M recirculation 1 −20.7 dB unchanged, decay |error| 8 kHz 0.5 → 0.4 %; L recirculation 1 −18.9 → −19.1 dB, decay |error| 1–8 kHz 1.0 / 0.6 / 0.5 / 0.5 → 0.9 / 0.5 / 0.4 / 0.4 %.

Everything moves by ≤ 0.2 dB / 0.1–0.4 %, inside the capture-pair spread, and mostly towards the captures at 2–8 kHz, where the HF pole acts. The static M / L metrics cannot distinguish the two maps strongly, so the reason to adopt is that S shows the curve directly, and one map for all three programs is what the plugin needs. Gate (`--tag mlcurve32`): PASS, M / L early decay |error| 1.12 → 1.11 / 0.61 → 0.60 %. The pulse check's 196 changed groups are replaced in `analysis/pulse/ml/`. The `ml_lines` law / fit / report stages have not been rerun.

## 9. The time variance: a slowly rotating coupling (Q3b, October 2026)

`.venv/bin/python tools/ml_variance.py [--stage rot|report|sim|all]` → `analysis/ml_variance/` (`rot/` one file per AMPL group, `report.log`, `sim/`, `sim.log`). Works on the residual capture − model of §8 (`ml_lines.fitted`: stored clocks and gains, nothing refitted).

> **Summary.**
> - **What varies is the coupling, and only in one way:** both 2 × 2 blocks of F (§7.3) rotate by **one common angle θ(t)**, F(θ) = F_b R(θ). The coupling stays orthogonal at every instant, so the decay is untouched, the variation is flat across frequency, and S, which has no coupling, has none. This is the "stochastic control … without perceptible FM" of the 1986 review. **S**
> - **θ is periodic, 0.83 Hz, with a random phase per capture**, free-running in seconds (one cosine fits M's and L's autocorrelations at 0.11, 0.22 and 0.44 s within 0.1); depth about **0.11–0.14 rad rms (6–8°, ±9–11° if sinusoidal)**. **S** (mechanism, rate), **H** (sine shape, exact depth)
> - **Accuracy test passed (`NEXT_STEPS.md` §2, metric 3):** the model with θ = 0.155 sin(2π · 0.83 t + φ), rendered through a simulated exponential sweep with random φ, decorrelates between renders the way Tomkin's AMPL pairs do, from the first recirculation to 2–4 s; the library's spread of model and repeat nulls is the periodic θ's signature. **S**
> - **The 0.2–2 % RT shortfall of §8.3 follows from it:** the simulated sweep of the rotating model shows the same shortfall (largest in L at 125 Hz), so the unit's decay is the law's and the sweep loses the decorrelating part. **S**
> - **Order:** the rotation acts on the line outputs and the 20 taps read after it (the read-out form, the transpose of §8's algebra). **S**

### 9.1 Method: kernels on derivatives, held out across channels

The residual past the first pass is mostly each capture's own random part (§8.1). Any candidate for what varies is a derivative of the model. A parameter that wanders during a sine sweep shows up, after deconvolution, as a **kernel on that derivative's response**: a slow modulation at f_m becomes a pair of ghosts at lags ±L·f_m/f for an exponential sweep with time constant L (lags shrinking with frequency). So per capture and window the residual is fitted as kernel ∗ derivative (ridge, lags ±300 in M, ±600 in L) on **one channel** and scored on **the other**: L and R read the same four lines through different taps, so only a line-side variation predicts across channels.

### 9.2 What is ruled out, and what fits

Recirculation 1 (`M 1.2_5.0_0.9`, first prototype; library-wide numbers in §9.4):

| Candidate (per capture) | Cross-channel null, from −18.5 dB |
|---|---|
| Static change of the 4 line gains, 4 line lengths or the 8 coupling entries (also per octave) | −18.5…−18.6 (nothing) |
| Kernels per line, per destination, per connection (±4…±128 lags) | ≤ 1 dB, or worse |
| **One kernel on the rotation derivative F_b J (J = [[0, −1], [1, 0]]), both blocks** | **−30.5…−33.7** |
| the two blocks with separate kernels | the same (one angle) |
| the blocks rotated in opposite senses; a common gain change of F | −18.3…−18.7 (nothing) |

- **Unique:** over all 128 sign patterns on F's support (one kernel each), the rotation pattern (+ on line 1 ← loop 2 and line 4 ← loop 3, − on the other six) is best by 5.7 dB; patterns on F's zero entries do nothing. Algebraically Fᵀ·dF is antisymmetric, the same J in both blocks: the tangent of a rotation.
- **Per capture, not common:** a capture's kernel applied to the other captures of its group makes them worse, and recirc1 kernels of different captures are uncorrelated (§9.4). One group (`M 1.2_2.6_0.9`) has a residual 90 % common to its captures (model error 8 dB above the variation there); its common part is also fitted by a rotation kernel, so the model's error and the variation are not separable by this test in that group.

### 9.3 Rate, timing and order

- **Ghost lags halve per octave:** the recirc1 kernel's median |lag| is 77 / 37 / 19 / 12 samples at 1 / 2 / 4 / 8 kHz in M and 73 / 38 / 19 / 12 in L (§9.4), lag × f ≈ 2.8–2.9 s·Hz (3.6 at 8 kHz, where 12 samples is coarse), the same in both programs, and it is symmetric (98–99.6 % of its energy): the exponential-sweep signature of a slow real modulation (ghost pairs at ±Δt(f) combine into a zero-phase kernel for any modulation phase). This fixes only f_m · L.
- **Rate in seconds, from the passes:** in recirculation k the rotation acts at up to k coupling events, one line length apart. One kernel per event, each correlated with recirculation 1's (1–9 kHz), gives θ's autocorrelation at a known lag with no sweep needed. Calibrated on simulations with a known θ (§9.5): at 0.83 Hz the estimator returns 0.83 / 0.42 for a true 0.84 / 0.42 (M) and 0.42 / −0.58 for 0.42 / −0.65 (L), against the library's 0.82 / 0.38 and 0.42 / −0.56. Library medians (§9.4): M ρ(0.11 s) = 0.78–0.82, ρ(0.22 s) = 0.38; L ρ(0.22 s) = 0.31–0.42, ρ(0.44 s) = −0.56. **One cosine, f_m = 0.83 Hz, fits all six** (0.84 / 0.41 / −0.66); per pass the programs disagree (one pass: M 0.8, L 0.4), so **θ runs in time, not per pass**. With f_m · L ≈ 2.85 this puts Tomkin's sweep time constant at L ≈ 3.4 s (≈ 24 s for a 20 Hz–20 kHz sweep, **H**).
- **Order:** grouping each recirculation's events by when they happen in the read-out form (rotation on the line outputs, taps after) beats the opposite grouping held out (L recirc2 −25.0 against −22.6 dB, recirc3 −20.0 against −16.6). So the lines are written with the diffused input plus the rotated line outputs, and the 20 taps read the lines, as in the EMT 250 patent's tapped read-out (`period_dsp_practice.md`).
- **Not tied to the output kernel (§7.5):** the cluster does not pass through the lines, so the rotation cannot make it; the kernel's cause stays open.

### 9.4 Library-wide (rot stage)

Every M and L group (M 160 groups / 479 captures, L 180 / 360), per capture, held out across channels (mean of L → R and R → L), dB re the capture, median [10–90 %]:

| Window | M base | M rotation | L base | L rotation |
|---|---|---|---|---|
| Recirculation 1 | −21.3 [−27.1 … −18.9] | **−33.1** [−35.2 … −29.9] | −19.3 [−22.3 … −18.2] | **−31.4** [−32.8 … −30.4] |
| Recirculation 2 | −16.2 [−17.6 … −15.5] | **−27.3** [−28.9 … −19.6] | −15.9 [−16.8 … −15.3] | **−25.7** [−26.2 … −22.6] |
| Recirculation 3 | −12.8 [−13.4 … −11.9] | **−20.2** [−22.9 … −12.2] | −12.7 [−13.2 … −12.4] | **−20.7** [−21.2 … −17.4] |

| Recirculation 1 control (median) | M | L |
|---|---|---|
| two blocks, separate angles | −32.6 | −30.8 |
| blocks in opposite senses | −20.9 | −18.9 |
| common gain change of F | −21.6 | −19.7 |
| **recirc2 / recirc3: one kernel for all events** | −25.1 / −18.8 | −21.3 / −15.4 |
| **recirc2 / recirc3: events grouped in the opposite time order** | −25.8 / −19.4 | −22.9 / −17.2 |

- **The rotation removes 12 dB of recirculation 1 in both programs** (8–11 dB in recirc 2, 7–8 in recirc 3), at every RT, LF and HF step (recirc1 rotation −30.8…−34.5 dB in M, −30.9…−32.4 in L per step; least at HF 0.1, where little is left above 1 kHz). The controls stay at the base level everywhere.
- **θ's size is the same everywhere:** kernel norm 0.110 [0.107 … 0.113] rad (M), 0.109 [0.106 … 0.112] (L); per RT, LF and HF step 0.106–0.113. Depth does not depend on the program or a control. (Calibrated, §9.5: about 0.14 rad true.)
- **Per capture:** recirc1 kernels of different captures of one group correlate +0.02 (M) / +0.16 (L) in the median, [−0.8 … +0.8]: each capture has its own phase.
- **Ghosts per octave** (recirc1, median |lag| holding 50 / 90 % of the kernel's energy, samples; lag50 × f; symmetric fraction): M 77 / 273, 2.92, 0.989 (1 kHz); 37 / 65, 2.81, 0.994 (2 kHz); 19 / 32, 2.88, 0.995 (4 kHz); 12 / 18, 3.64, 0.996 (8 kHz). L 73 / 145, 2.77, 0.984; 38 / 68, 2.88, 0.987; 19 / 32, 2.88, 0.994; 12 / 17, 3.64, 0.995.
- **θ's autocorrelation** (median correlation with recirc1's kernel, divided by the same-instant one): M 0.82 at 0.11 s (recirc2), 0.78 / 0.38 at 0.11 / 0.22 s (recirc3); L 0.42 at 0.22 s, 0.31 / −0.56 at 0.22 / 0.44 s. Same-instant correlations 0.94–0.98: recirculations 1–3 see one and the same θ.

### 9.5 The accuracy test: the rotating model through a simulated sweep (sim stage)

The model in the time domain in the read-out form (blocks of one shortest line; with θ = 0 it reproduces §8's model to about −100 dB), θ(t) = 0.155 sin(2π · 0.83 t + φ), measured with a 10 Hz–13 kHz exponential sweep of time constant 3.4 s at the device rate and deconvolved; 6 renders with random φ per setting. Pair nulls (one gain per window and channel; recirc windows R, time windows mean of L and R), median [range]:

| Setting (LF / RT / HF) | | recirc1 | recirc2 | recirc3 | 0.3–0.6 s | 0.6–1.0 s | 1.0–2.0 s | 2.0–4.0 s |
|---|---|---|---|---|---|---|---|---|
| M 1.2 / 5.0 / 0.9 | Tomkin | −17.0 [−15.0 … −19.7] | −13.8 [−12.0 … −16.3] | −11.0 [−9.1 … −13.7] | −10.7 [−9.0 … −13.7] | −8.2 [−6.6 … −11.1] | −6.5 [−5.1 … −8.6] | −6.5 [−5.6 … −7.6] |
|  | model | −17.4 [−15.3 … −60.0] | −14.3 [−11.9 … −56.9] | −11.6 [−9.2 … −53.9] | −11.5 [−9.3 … −53.7] | −8.5 [−6.6 … −50.9] | −6.7 [−4.6 … −49.1] | −6.4 [−4.6 … −48.8] |
| M 2.0 / 0.8 / 0.9 | Tomkin | −19.6 [−17.5 … −23.6] | −14.9 [−13.1 … −19.0] | −11.4 [−9.3 … −15.1] | −11.9 [−10.0 … −15.8] | −10.0 [−8.4 … −14.0] | – | – |
|  | model | −22.9 [−18.9 … −46.4] | −18.2 [−14.4 … −41.5] | −14.3 [−10.3 … −37.0] | −14.9 [−10.9 … −37.9] | −11.0 [−7.5 … −35.2] | −8.9 [−5.3 … −33.1] | −10.8 [−6.6 … −35.0] |
| M 0.7 / 3.0 / 0.2 | Tomkin | −18.4 [−15.9 … −18.6] | −14.1 [−11.2 … −14.9] | −11.0 [−8.2 … −12.7] | −10.7 [−8.1 … −13.0] | −8.0 [−5.5 … −11.5] | −6.3 [−4.3 … −10.1] | – |
|  | model | −21.3 [−17.4 … −49.3] | −17.0 [−12.8 … −44.9] | −13.6 [−9.3 … −41.7] | −13.8 [−9.6 … −41.4] | −10.3 [−7.1 … −37.6] | −8.2 [−4.4 … −36.2] | −7.0 [−3.5 … −34.1] |
| M 0.1 / 1.8 / 0.4 | Tomkin | −15.5 [−15.1 … −26.6] | −10.8 [−10.4 … −22.9] | −7.7 [−7.5 … −20.1] | −7.7 [−7.4 … −20.2] | −5.2 [−5.1 … −17.8] | – | – |
|  | model | −24.3 [−17.8 … −36.1] | −20.0 [−14.0 … −32.1] | −16.7 [−10.8 … −28.9] | −17.0 [−10.8 … −29.0] | −14.0 [−7.3 … −25.8] | −11.6 [−5.6 … −23.4] | −10.6 [−4.2 … −22.6] |
| L 1.2 / 20.0 / 1.0 | Tomkin | −18.0 | −15.0 | −12.1 | −16.7 | −12.3 | −7.9 | −5.0 |
|  | model | −17.0 [−14.6 … −59.7] | −14.2 [−12.0 … −56.6] | −11.4 [−9.2 … −53.8] | −15.8 [−13.6 … −58.4] | −11.6 [−9.4 … −54.0] | −7.2 [−5.2 … −49.6] | −4.5 [−2.6 … −46.6] |
| L 1.2 / 4.5 / 0.8 | Tomkin | −20.4 | −17.3 | −13.7 | −19.2 | −14.2 | −10.9 | −8.0 |
|  | model | −19.6 [−15.7 … −43.0] | −16.5 [−12.7 … −40.0] | −13.3 [−9.5 … −36.9] | −18.3 [−14.4 … −41.8] | −13.6 [−9.9 … −37.1] | −9.8 [−6.1 … −33.2] | −6.7 [−3.2 … −30.0] |
| L 2.0 / 2.6 / 0.5 | Tomkin | −15.9 | −12.6 | −9.0 | −14.4 | −9.4 | −5.9 | −2.5 |
|  | model | −20.8 [−16.7 … −48.6] | −17.3 [−13.5 … −44.8] | −13.7 [−9.9 … −41.2] | −19.1 [−15.2 … −46.9] | −14.0 [−10.2 … −41.6] | −10.1 [−6.5 … −37.6] | −6.8 [−3.3 … −33.9] |
| L 0.1 / 10.0 / 0.3 | Tomkin | −19.3 | −15.9 | −13.1 | −17.8 | −13.2 | −8.9 | −6.2 |
|  | model | −22.3 [−16.1 … −34.3] | −19.4 [−12.9 … −31.2] | −16.2 [−9.9 … −28.1] | −21.1 [−14.7 … −33.0] | −16.5 [−10.2 … −28.4] | −12.3 [−6.2 … −24.2] | −9.3 [−3.4 … −21.1] |

Tomkin: its AMPL pairs (M 3, L 1 per setting); model: 15 pairs of 6 renders; ranges least to most similar. A pure sine makes some render pairs nearly identical (Δφ ≈ 0), hence the deep model minima.

- **The library's null spreads are a periodic θ's signature** (`report`, from §8's law files, R): against the θ = 0 model every capture is nearly the same distance (IQR 4.6 / 1.3 / 0.6 dB in M recirc1–3, 1.7 / 0.9 / 0.4 in L), while repeat pairs spread 5.4–7.5 dB. A sine with a random phase per capture predicts exactly that: a pair differs by 2(1 − cos Δφ) times the model residual, Δφ uniform, so the repeat median is the model median + 3.0 dB and the 90th percentile + 5.9 dB. Observed M recirc2 −14.1 / −10.4 against −13.0 / −10.1 predicted, L −13.4 / −10.4 against −12.7 / −9.8. Only the deep tail is shallower (10th percentile −20.7 against −26.1 in M recirc2; 1st −25 against −46): about −25 dB of each capture's variation is not periodic, or the repeat fit's freedoms floor it (**H**). Narrowband noise in place of the sine would spread the model nulls as much as the pairs.
- **Tomkin's pairs lie inside the model's range** (between its least and most similar render pair), or within 1.5 dB of its least similar end (M 2.0 / 0.8 / 0.9, L 2.0 / 2.6 / 0.5), in every window out to 2–4 s, **except M 0.1 / 1.8 / 0.4** (LF 0.1, short RT), whose two less similar pairs are 2.2–3.6 dB less similar than any render pair (θ larger there, a capped-LF effect, or chance with 6 phases: open). The medians differ by −1…+9 dB (the model's more often deeper), because 6 random phases sample 2(1 − cos Δφ) coarsely and Tomkin has 1–3 pairs per setting; the library-wide distribution above is the stronger test. Nothing in the model is fitted to these pairs: θ's depth, rate and the sweep constant come from the rot stage.
- **The RT shortfall (§8.3):** decay per octave of the rotating renders against the θ = 0 model: M 0.992 / 1.000 / 1.000 / 0.996 / 0.998 / 0.997 / 0.996 (125 Hz–8 kHz), L 0.958 / 0.992 / 0.984 / 0.985 / 0.992 / 0.994 / 0.992; measured (§8.3, LF 1.2) M 0.988–0.997, L 0.972 at 125 Hz and 0.987–0.997 above. The sweep of a time-variant unit decays faster than its energy; the line laws need no correction.
- **Depth calibration:** on the renders the recirc1 kernel norm reads 0.081–0.087 for a true 0.110 rad rms (ridge shrinkage and the ±LK window), so the library's kernel norm (§9.4) corresponds to about 0.14 rad; the pair nulls favour 0.11–0.12. Depth 0.11–0.14 rad rms (**H** within that range).

### 9.6 What this answers, and what's next

- **Q3b, which part varies: answered.** One common rotation angle of the two coupling blocks; nothing else in the lines, taps or read-out varies at a level this test sees (the leftover after the rotation kernel is −30…−34 dB in recirc1, near the kernel method's floor on simulations, −35…−37 dB).
- **Per pass or per second: per second** (the same 0.83 Hz in M and L), periodic with a random phase per capture (**H**: a sine; a triangle or other smooth periodic shape is not excluded).
- **Budget:** 8 coefficient multiplies per sample as before (the two blocks), with (cos θ, sin θ) rewritten by the control side at a rate well above 0.83 Hz; a table or recursive oscillator suffices.
- **Next:** depth done in §9.7 (0.15 rad rms); the shape is not resolvable from this library (§9.7: only a square-like θ is disfavoured); the non-periodic remainder (the shallow deep tail), and whether θ's rate or depth depends on the program or a control (§9.4 shows none at this precision). For the plugin: the read-out form with F(θ(t)), θ = a sin(2π f_m t + φ₀), a ≈ 0.155–0.2 rad, f_m ≈ 0.83 Hz.

### 9.7 θ's depth: 0.15 rad rms by two methods (25 October 2026) **[S: Tomkin's AMPL pairs match sine-θ renders at 0.15 rad rms, and the kernel estimator's reading on renders matches the library's at a true 0.146; H: the sine shape]**

```
ML_THETA_RMS=<rms> EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 8 -- .venv/bin/python tools/ml_variance.py --stage sim --prog ML --workers 8
                                     # 22 s, 6.8 GB -> analysis/ml_variance/sim_<rms>/, sim_<rms>.log (same seeds, so the same phases, as sim/ at 0.110)
```

§9.5 left the depth at 0.11–0.14 rad rms. The pair nulls were read as favouring 0.11–0.12, the calibrated kernel norm as about 0.14. The sim stage now renders the same 8 groups and 6 phases per group at 0.13, 0.15 and 0.17 rad rms beside §9.5's 0.110. Each Tomkin AMPL pair is compared with its group's median render-pair null over recirculations 1–2 and 0.15–0.3 s (48 pair-windows: M 36, L 12; later windows are past the small-θ regime). Median gap (Tomkin − renders), with bootstrap 16–84 %, and the kernel estimator's reading on the renders:

| θ rms (rad) | gap, dB | M | L | kernel reading |
|---|---|---|---|---|
| 0.110 (§9.5) | +2.88 [+2.40, +3.09] | +2.88 | +1.15 | 0.082 |
| 0.130 | +1.39 [+0.94, +1.62] | +1.39 | −0.31 | 0.098 |
| **0.150** | **+0.16 [−0.28, +0.38]** | +0.16 | −1.54 | **0.113** |
| 0.170 | −0.91 [−1.35, −0.70] | −0.91 | −2.60 | 0.129 |

The gap moves by −1.49 / −1.23 / −1.07 dB per step, the 20 log₁₀ of the depth ratios (1.45 / 1.24 / 1.09), so the renders are in the small-θ regime and the scaling is exact. Tomkin's pairs match at **0.15 rad rms (0.14–0.16)**. The kernel estimator reads 0.75 × the true depth on renders, and the library's own reading (§9.4: 0.110 M, 0.109 L) corresponds to a true **0.146**. The two methods agree. L's four pairs alone point lower (≈ 0.125), but that is 4 pairs, one per group; L's kernel reading says 0.145 like M's. **Depth: θ = a sin(2π · 0.83 t + φ), a ≈ 0.21 rad (0.15 rms, ±12°), the same in M and L (H for a program difference: none shown).** The `sim` default stays at §9.5's 0.155 amplitude for its stored output; `ML_THETA_RMS` sets another. §9.5's "pair nulls favour 0.11–0.12" was read by eye per setting; pooled, they don't. Unchanged: the shape (sine, H), the sweep (unknown, only f_m × SWEEP_L is set), Tomkin's pairs inside the renders' range. Note: these renders use the M / L pole map of §8.7, §9.5's the older one; the static model is common to both renders of a pair, so the pair nulls don't depend on it.

**The shape is not resolvable from this library** (H stays: a sine). The evidence for periodicity is θ's autocorrelation across passes (§9.4), fitted by a cosine at 0.11 / 0.22 / 0.44 s to ±0.1. A triangle wave's autocorrelation differs from a cosine by at most 0.020 and a square wave's by 0.21, so only a square-like θ is disfavoured. Harmonics would add ghosts at 3, 5, … × the fundamental's lag. A triangle puts 1.2 % of θ's power in its 3rd harmonic, so those ghosts would sit about 19 dB under the main ones, near −39 dB re the recirculation. That is below the kernel method's floor on simulations (−35…−37 dB, §9.6). For the plugin a sine is the simplest choice consistent with every test, and a smooth non-sinusoidal shape would not be audibly different at ±12°.

## 10. The first-pass gap was the deconvolution kernel, not the algorithm (October 2026)

`.venv/bin/python tools/firstpass_gap.py --stage collect|report|chain [--model carry|lines]` → `analysis/firstpass_gap/` (`report_*.log`, `chain.{json,log}`; the per-group residual caches `groups_*/` are not committed). `NEXT_STEPS.md` §5 item 5: the first pass (cluster + diffuser + 20 taps, §7.6 / §8.1) sat 4–15 dB above the repeat-capture null in S, M and L alike. What is the residual?

> **Summary.** It is a fixed, setting- and program-independent error of the library's device-domain IRs: the 5.2 ms deconvolution kernel (`analysis/library/tomkin96/kernels/kernel_96000.npz`, the mean of 8 R tap-0 waveforms, used for both channels) leaves a small low-frequency tail 2–6 ms after every impulse. **A per-channel correction kernel fitted only on the cluster-only regions (no reverb sample in the fit) predicts the first pass to the repeat-capture level, held out**: S −58.4 / −57.8, M −53.9 / −49.6, L −49.5 / −46.5 dB (L / R), against −44.7 / −45.7, −44.5 / −44.1, −44.3 / −42.7 without it. The reverb structure of §7–8 needs no change. **S** (cause: **H**, strongly indicated). *Update (§12): cause confirmed. A per-channel chain kernel (19 ms FIR × the unit's high-pass, fitted on the same cluster-only windows) holds exactly this ripple. With the caches rebuilt from it (`EMT_KERNEL=chain2`), the correction is no longer needed: a correction refitted on the new IRs is ≤ 0.04 dB / 0.4° in band.*

### 10.1 Without a clock fit, part of the gap is the clock

The §7.6 bench fits no capture clock and renders the taps at their measured fractional positions, so the diffuser's echoes stay on the nominal clock while the taps carry the capture's (10–15 ppm). Its residual's kernel (one signal read through each channel's taps) is antisymmetric ±1-sample pairs at lag 0 and at every all-pass delay (S 20 / 72, M 42 / 144, L 82 / 286): a time-scale error. With `ml_lines.Model` (integer taps, per capture a shift and clock fitted on [0, onset + shortest line), as §8.1) the first pass is S −43.6 / −45.7, M −42.6 / −44.1, L −42.2 / −42.7 dB (medians over all captures; repeat −48.2 / −62.3, −47.3 / −49.2, −46.0 / −54.1). Everything below uses this model.

### 10.2 What is left is one waveform, the same everywhere

Residual of the first pass (capture minus its fitted cluster, minus the fitted tap part), b = 0 settings (S 368, M 311, L 280 captures):
- **The mean residual, normalised by the tap gain and learned on the even-indexed settings, removes the gap on the odd ones:** S −43.8 / −45.8 → −49.9 / −62.9, M −42.8 / −43.9 → −49.2 / −49.4, L −42.9 / −42.7 → −51.0 / −56.0 dB (repeat as above). Its two halves correlate 0.993–1.000. Short FIRs on the tap part (±4…±48 taps, ridge) gain only 0.3–1.7 dB.
- **60–83 % of it lies at 100–400 Hz,** at the same level re the data in all three programs (−29 dB at 100–200 Hz, −35 dB at 200–400 Hz).
- **Two parts, by AMPL:** in R it does not depend on AMPL (scale 0.92–1.00, correlation 0.94–1.00); in L a further part grows with AMPL and sits just after the last cluster tap inside L's first-pass window (S tap 288, M 576, L 1152): the cluster's own low-frequency tail.
- **The data/model transfer function is one curve for S, M and L, both channels** (pooled, ≈25 Hz smoothing): −0.2 dB at 160 Hz, +0.13 at 280, −0.08 at 400, +0.05 at 520 Hz, phase ±1–2°, flat (≤ 0.03 dB) above 600 Hz. A ripple fixed in device hertz, not scaled ×2 / ×4 with the program, so not in the structure. Its time kernel is causal, at lags 50–150 (peak 93), and each program's correlates 0.88–0.98 with the pooled one.

### 10.3 A chain-correction kernel, fitted on the clusters only (`chain` stage)

One causal kernel g (lags −32…400) applied to the whole fitted model, m → m + g ∗ m; per held-out capture only the two gains per channel refitted (shift and clock as fitted without g). Medians L / R, dB, over all captures of the program; "first pass" held out in every row except where marked:

| Kernel (fitted on) | Free numbers | S cluster/gap | S first pass | M cluster/gap | M first pass | L cluster/gap | L first pass |
|---|---|---|---|---|---|---|---|
| none | 0 | −49.4 / −46.3 | −44.7 / −45.7 | −48.5 / −45.7 | −44.5 / −44.1 | −44.9 / −45.5 | −44.3 / −42.7 |
| joint L+R, on another program (S row: from M) | 433 | −52.3 / −54.0 | −50.7 / −53.6 | −54.9 / −52.5 (from S) | −51.4 / −47.8 | −52.8 / −51.3 (from S) | −48.4 / −45.1 |
| per channel, even settings, all programs (odd settings scored) | 866 | −59.7 / −58.5 | −59.1 / −57.7 | −63.3 / −64.0 | −54.4 / −49.7 | −57.4 / −57.2 | −49.8 / −46.7 |
| **per channel, cluster-only regions of all captures** | 866 | −61.6 / −59.7 (in-sample) | **−58.4 / −57.8** | −63.5 / −64.1 (in-sample) | **−53.9 / −49.6** | −57.1 / −57.0 (in-sample) | **−49.5 / −46.5** |
| per channel, cluster-only regions of M and L only | 866 | **−57.9 / −57.1** | **−56.6 / −56.0** | (in-sample) | −53.7 / −49.6 | (in-sample) | −49.8 / −46.7 |
| *repeat capture* | | | *−48.2 / −62.3* | | *−47.3 / −49.2* | | *−46.0 / −54.1* |

- **Excess over the repeat-capture null, first pass, cluster-only kernel:** S −10.2 / +4.5, M −6.6 / −0.4, L −3.5 / +7.6 dB (negative: the model nulls deeper than another capture does; the L channel's repeat null is limited by capture-to-capture differences, unit facts §2: L's noise floor is 6–9 dB above R's). The L program's R channel keeps +7.6 dB.
- **Program S held out completely** (kernel from M's and L's clusters only) reaches −56…−58 dB in both windows: the correction is not per program.
- Per channel matters (S cluster/gap −53 → −60 dB): the library kernel was built from R tap-0 waveforms and is applied to L too, and each channel has its own output S&H and trimmer (unit facts §2b).
- **In band the correction is small:** −47 dB of energy; 64–67 % in lags 40–200 (peak 93–98), 20–23 % within ±10 samples. Response L / R: +0.18 / +0.17 dB at 100 Hz (−2°), −0.23 / −0.25 at 150, −0.17 / −0.18 at 200, +0.15 / +0.13 at 300 Hz, within ±0.06 dB from 500 Hz to 10 kHz (L −0.06 dB at 10 kHz). Its out-of-band part (−13 dB, lags −12…−8) is unconstrained and is cut off wherever the kernel is applied.
- **Reading (H):** the analogue and capture chain's response continues past the 5.2 ms kernel (its tail is −52 dB at 4 ms, and the high-passes at 0.2–3.4 Hz, `analog_results.md` §3.5, are far longer), so the truncated kernel's deconvolution leaves this low-frequency ripple on every impulse, cluster and reverb alike. It belongs to the measurement pipeline, not to the unit's algorithm, and affects every device-domain comparison made so far by about this much. The principled fix is a longer, per-channel chain kernel in `tools/harness` (rebuild of the IR caches); until then the correction is `chain.json:cluster_perch`. *Done (§12): `chain2`; confirmed **S**.*
- **What it explains elsewhere (H):** the free 28 k-coefficient S injection's 8–10 dB lead over the compact model on the library shrank to 1–2.5 dB on CLA and the Altiverb (`model_s_results.md` §16.12, §18.1): a free injection can absorb a library-wide kernel error; on another chain there was none to absorb.

## 11. Out of band: 25–100 Hz and 10–12 kHz (October 2026)

`.venv/bin/python tools/outofband.py --stage score|report|sweep` → `analysis/outofband/` (`report.{json,log}`, `sweep.{json,log}`; the per-group files `groups/` are not committed). Every null so far is 100 Hz–10 kHz. Do the existing S, M and L models (`ml_lines.Model`, all three programs) hold outside it? **Nothing is refitted:** each capture keeps the law stage's in-band shift, clock and two gains per channel; the model is rendered full band and multiplied by the pipeline's measured response P(f) (flat to 11 kHz, −3.6 dB at 11.5, −19 dB at 12 kHz, nothing from 12.3 kHz: the deconvolution's band window, `tools/harness/core.py`). The data are the uncorrected device-domain IRs (the §10 correction is defined only in band). Bands: 25–50, 50–100, 10–11 k, 11–12 kHz; 100–200 Hz and 5–10 kHz as in-band controls.

> **Summary.** No sign that the band limit hid a structural error. **10–12 kHz:** the decay matches in every program and HF step (RT capture / model 1.000 ± 0.003, level within 0.1 dB), so the HF laws and the 0.84 pole limit hold to 12 kHz; the nulls expose two measurement-side items, S's L-channel chain (the R-only deconvolution kernel, §10) and **the M/L output kernel's phase above 9 kHz** (L +7° at 10 kHz, +18° at 12 kHz; magnitude within 0.08 dB of a symmetric kernel). **25–100 Hz:** S's decay matches (RT ratio 1.00–1.01, LF 0.1 +2…4 %, near the measurement's noise); **M and L fall short by up to 16 % at high LF and long RT, and the rotating model through a simulated sweep reproduces this** (static: no shortfall), so it is how a sweep measures the time-variant unit, not an LF-law error. Low-band nulls and levels are limited by the pipeline (repeat captures null only −14…−36 dB there), not by the model.

### 11.1 High band

Median nulls, strict (in-band gains) L / R [repeat capture], dB:

| Program | Window | 5–10 kHz (control) | 10–11 kHz | 11–12 kHz |
|---|---|---|---|---|
| S | first pass | −50.0 / −52.1 [−57.9 / −59.9] | −41.2 / −47.4 [−46.6 / −48.6] | −32.7 / −43.0 [−36.7 / −41.7] |
| S | recirc 1 | −46.5 / −47.9 [−61.5 / −59.1] | −40.4 / −43.8 [−48.4 / −47.7] | −32.8 / −39.8 [−41.2 / −40.0] |
| M | cluster and gap | −53.5 / −50.4 | −47.3 / −40.1 | −40.7 / −33.2 |
| M | first pass | −51.3 / −48.9 [−58.1 / −67.1] | −47.1 / −39.4 [−51.6 / −60.7] | −40.8 / −33.1 [−45.5 / −54.2] |
| L | cluster and gap | −50.9 / −48.8 | −30.8 / −29.4 | −21.3 / −19.5 |
| L | first pass | −46.0 / −42.3 [−52.6 / −65.2] | −31.3 / −29.2 [−32.8 / −55.8] | −21.3 / −19.8 [−22.9 / −48.2] |

(The cluster-and-gap repeat uses the model's cluster, as in `ml_lines.fit_repeat`, so it is no reference there. M/L recirculations are time-variance-limited, model ≈ repeat.)

- **Decay, 10–11 and 11–12 kHz:** RT capture / model 1.000 / 1.000 (S), 1.000 / 0.999 (M), 1.001 / 1.001 (L), every HF step within 0.997–1.004; level at onset −0.08…+0.08 dB. **S**
- **S:** R at the repeat-capture level to 12 kHz; L 4–8 dB worse at 11–12 kHz in every window including the cluster's, i.e. the chain: the deconvolution kernel was built from R alone (§10). **S**
- **The M/L output kernel, measured** in the R channel's cluster-only window (`--stage kernel` → `kernel.{json,log}`; 4 known taps; 161 M / 60 L captures, gain matched at 2–6 kHz):

  | | 8 kHz | 10 kHz | 11 kHz | 12 kHz |
  |---|---|---|---|---|
  | M measured |K| / phase | −0.54 dB / −0.2° | −0.72 / +1.2° | −0.76 / +2.2° | −0.75 / +2.7° |
  | M stored 11-tap / 5-tap symmetric fitted in band | −0.55 / −0.54 | −0.76 / −0.72 | −0.85 / −0.78 | −0.92 / −0.83 |
  | L measured |K| / phase | −2.75 dB / −0.1° | −4.07 / +6.8° | −4.66 / +12.4° | −5.04 / +18.4° |
  | L stored 11-tap / 5-tap symmetric fitted in band | −2.76 / −2.76 | −4.13 / −4.05 | −4.79 / −4.59 | −5.38 / −4.98 |

  The magnitude is a symmetric kernel's (5 taps fitted in band predict it to 0.08 dB at 12 kHz; the stored kernel's unconstrained outer taps are 0.1–0.35 dB off). **The phase is not** (+18° at 12 kHz in L, −3° at 4 kHz): a non-symmetric part that a ±1-sample address jitter averaged by the sweep (§7.5, **H**) does not give. It is what limits L's high-band nulls (cluster window −21 dB at 11–12 kHz). **S** (measurement) / open (cause). Fix: refit the kernel over 100 Hz–12 kHz with P(f) on the cluster windows.

### 11.2 Low band

- **Nulls are pipeline-limited:** strict model nulls −13…−23 dB at 25–100 Hz in every program and window, but the repeat captures reach only −14…−36 dB in the early windows too (S recirc1 at 50–100 Hz: repeat −51 / −56, model −20 / −20). Adding the unit's sub-20 Hz high-passes to the model makes it 2–4 dB worse (as `model_s_results.md` §16.11 found in band). Level at onset: captures +1.2 dB (25–50 Hz) and +0.5…0.7 dB (50–100 Hz) above the model in all three programs alike: the deconvolution kernel's low-frequency error growing below 100 Hz (§10.3: +0.18 dB at 100 Hz). So waveform and level below 100 Hz can't be tested with these IRs until the chain kernel is fixed (**S**).
- **Decay, S:** RT capture / model 1.007 (25–50 Hz) and 1.008 (50–100 Hz); LF 1.0 / 2.0 within 1.000–1.006 at every RT, LF 0.4 1.005–1.045, **LF 0.1 +2…4 % at every RT** (repeat-pair spread 2–3 % at 25–50 Hz, 0.5–0.8 % at 50–100 Hz, so LF 0.1 is just above it). **S** (match); the LF 0.1 excess is small (**H**: the DC floor or β regime).
- **Decay, M and L: the captures decay faster than the model at high LF, increasingly with the model's LF RT:**

  | | LF 0.1 | LF 0.7 / 0.8 | LF 1.2 | LF 2.0 |
  |---|---|---|---|---|
  | M, 25–50 Hz (RT 5.0) | 1.011 | 0.899 | 0.876 | 0.840 |
  | M, 50–100 Hz (RT 5.0) | 0.988 | 0.995 | 0.998 | 1.000 |
  | L, 50–100 Hz (RT 20) | 1.004 | 0.883 | 0.859 | 0.879 |
  | L, 100–200 Hz (RT 20, in band) | 0.997 | 0.962 | 0.962 | 0.969 |

  (RT capture / model, median over HF; the full RT × LF grids per band are in `analysis/outofband/report.log`.) In L at 50–100 Hz: 0.99 at RT_law 4–5 s, 0.92–0.95 at 7–8 s, 0.83–0.88 beyond 11 s. S, time-invariant, has none at the same per-pass loss.
- **The rotating model reproduces it through a simulated sweep** (`outofband.py --stage sweep` → `sweep.{json,log}`: `ml_variance.sweep_measure` with the sim's band widened to 15 Hz–12 kHz; θ = 0.155 sin(2π · 0.83 t + φ), 3 phases, 12 s tail; RT rotating / static, same decay measure):

  | Setting | 25–50 Hz | 50–100 Hz | 100–200 Hz | 5–12 kHz | captures / model at 25–50 / 50–100 / 100–200 |
  |---|---|---|---|---|---|
  | M 2.0_5.0_0.9 | 0.860 | 1.004 | 0.994 | 0.996–1.005 | 0.840 / 1.000 / 0.999 |
  | M 0.1_5.0_0.9 | 0.993 | 1.003 | 1.005 | 0.997–1.004 | 1.011 / 0.988 / 1.003 |
  | L 2.0_20.0_1.0 | 0.922 | 0.839 | 0.977 | 0.991–0.996 | 0.874 / 0.879 / 0.969 |
  | L 2.0_8.0_1.0 | 0.945 | 0.869 | 0.971 | 0.995–0.998 | 0.834 / 0.827 / 0.936 |
  | L 0.1_20.0_1.0 | 0.996 | 0.984 | 0.997 | 0.989–0.997 | 0.993 / 1.004 / 0.997 |

  Including the program-dependent pattern (M short at 25–50 Hz only, L most at 50–100 Hz). The captures' shortfall is somewhat larger in L 2.0_8.0 (the sweep length and θ's depth are not known exactly, §9.5). **Reading (S):** a sine sweep of a unit whose coupling rotates under-reads long low-frequency decays; the rotation keeps the coupling orthogonal, so the unit's own decay follows the law. The M/L LF laws stand; for the plugin, use the law, not the captures' low-band RT.
- **Not L/R timing:** the L channel's poorer repeat nulls (in band and out) are not a per-capture L/R timing difference: per-channel shifts of repeat-capture pairs differ by 0.0003–0.0005 samples (std; S, L programs).

### 11.3 What this answers, and what's next

- **The algorithm was not mis-identified by the band limit:** structure, HF laws and pole limit hold to 12 kHz; S's LF decay holds to 25 Hz; M/L's LF decay holds once the sweep's time-variance artefact is accounted for.
- **Measurement-side items:** the deconvolution kernel (§10; also the low-band level and S's L channel at 11–12 kHz) and the M/L output kernel's phase above 9 kHz.
- *Update (§12.5), with the chain kernel: the low-band level offset is gone (+0.2 / −0.0 / +0.1 dB at 25–50 Hz); S's first pass reaches the repeat-capture null at 50–100 Hz (−37 / −42 dB) and its recirculations stay 15–20 dB above it (the model's LF side, now testable); S's L channel at 11–12 kHz matches R.*
- **Next:** (1) ~~a per-channel, longer chain kernel and a cache rebuild, then repeat this check below 100 Hz for level and waveform~~ done (§12); (2) ~~refit the M/L output kernel over 100 Hz–12 kHz on the cluster windows (with P(f)), then rescore the high band~~ refitted alongside the chain kernel (§12.5: L's 11–12 kHz first pass −21 → −35 dB); not yet adopted in the model; (3) S at LF 0.1 (+2…4 % RT at 25–100 Hz).

## 12. A per-channel chain kernel replaces the correction (October 2026)

`.venv/bin/python tools/chain_kernel.py --stage data|align|fit|write` → `analysis/chain_kernel/` (`fit.{json,log}`, `chain2.json`, `shifts.json`; the windows `data_*.npz` are not committed) and the kernel `analysis/library/tomkin96/kernels/chain2.npz`. Any tool that reads the device-domain IRs takes it with **`EMT_KERNEL=chain2`** in the environment (`tools/ml_structure.py`): the caches go to `analysis/ml_structure/cache_chain2/` and `analysis/q1/cache_chain2/` (+ `data_chain2.json`), and every result file gets the suffix `_chain2`. Without it everything reads the library kernel as before, so old and new compare directly. `ML_HK=analysis/chain_kernel/chain2.json` (`tools/ml_lines.py`) also swaps the model's stored output kernels for the ones refitted with the chain kernel (S: a unit tap); results get `_hk`. `NEXT_STEPS.md` §5 item 5, §10–11 above.

> **Summary.** The library's deconvolution kernel is replaced by a per-channel chain kernel: a **19 ms FIR per channel times the unit's analytic sub-20 Hz high-pass** (the schematic's 3.35 / 0.72 / 0.23 / 0.16 Hz, nothing fitted), fitted only on the cluster-only windows of 884 S, M and L captures with a capture clock per capture. **Without any correction the rebuilt IRs now carry no measurable transfer error** (S first pass ≤ 0.04 dB, ≤ 0.5° per band, was 0.14 dB / 1°); the low-band level offset is gone (25–50 Hz +1.2 → +0.2 / −0.0 / +0.1 dB in S / M / L) and S's L channel at 11–12 kHz now matches R. The §10 correction is no longer needed: on the compact S model (`model_s_results.md` §19–20) the new IRs beat the corrected old ones by up to 4 dB at 50–300 ms and match them later. The first pass with nothing fitted on the reverb reaches S −52.0 / −53.3, M −50.0 / −48.2, L −47.5 / −45.6 dB (repeat −54.7 / −63.6, −49.9 / −49.1, −50.9 / −54.1); S stays 1.4–6.3 dB short of the corrected old figure. That comes from model parts fitted on the old IRs (the stored output kernels and a residual waveform, §12.3), not the chain. **S** (pipeline); **H** (what limits S's first pass now).

### 12.1 The fit (`fit` stage; held out on the odd settings)

Before the reverb starts a capture holds only the cluster. Per channel there are 4 sign-only taps at known device offsets (§1), for M and L through the program's output kernel (§7.5), all with one gain per capture. The model is the chain kernel c_ch = h_ch ∗ HP:
- h_ch, a free FIR per channel on a 4× device grid (the taps sit on the device grid, so the problem splits into 4 polyphase problems with one Gram matrix);
- HP, the schematic's high-pass, divided out of the capture first (causal, from the file start), so its slow tail, far longer than any window, is known rather than truncated.

The windows run from 25 device samples before tap 0 to 6 before the reverb onset (`ml_carryover.ONSET`: 9 / 18 / 36 ms in L, 22 / 43 / 86 ms in R for S / M / L), so **no reverb sample enters the fit**. Per capture the fit has:
- a gain;
- a sub-sample shift and **a capture clock** (`align` stage, from 13 ppm; fitted 10–14 ppm). Without the clock, M's and L's long windows can't share a kernel with S: 13 ppm moves L's last cluster tap by 0.02 samples, −28 dB at 10 kHz;
- a slow additive term {1, t, t²} per channel, projected out of model and data. It has the same shape in every capture but does not scale with AMPL (M: 0.54 / 0.59 at AMPL 60 / 99), so it is capture-side (an offset, or a wrapped tail), not the chain.

The M and L output kernels (taps −8…8) are refitted alongside: the stored ones (fitted on the old IRs) leak into a joint kernel and cost 3–15 dB. S is weighted ×3 so S alone sets the short lags. The fit uses AMPL ≥ 30 captures: M's AMPL 10 clusters are 19 dB down, and 20 of them have none.

Nulls over the window, median over the AMPL ≥ 60 captures of the odd settings (kernel fitted on the even ones), L / R, dB, with a gain, shift, clock and slow term refitted per capture. Bands: 100 Hz–10 kHz / 25–100 Hz / 10–12 kHz:

| Kernel | Numbers | S | M | L |
|---|---|---|---|---|
| library kernel (5.2 ms, R only) | 0 | −41.2 / −45.4 · −40.4 / −24.4 · −33.6 / −49.0 | −43.8 / −42.4 · −27.5 / −20.5 · −41.5 / −36.7 | −38.7 / −37.0 · −24.2 / −20.3 · −24.7 / −22.6 |
| **chosen: FIR 500 × schematic high-pass** | 4200 + 34 | **−58.8 / −57.3 · −72.1 / −57.1 · −47.0 / −46.3** | **−65.6 / −61.9 · −59.3 / −41.8 · −56.7 / −52.6** | **−51.3 / −50.7 · −37.5 / −38.5 · −39.1 / −38.4** |
| output kernels stored, not refitted | 4200 | −55.8 / −51.9 · −71.4 / −34.6 · −43.7 / −42.6 | −52.0 / −47.4 · −59.7 / −26.2 · −38.8 / −36.4 | −39.3 / −38.7 · −37.7 / −23.8 · −23.8 / −22.9 |
| no high-pass (FIR alone) | 4234 | −58.8 / −56.2 · −72.3 / −44.9 · … | −65.8 / −55.5 · −59.6 / −31.5 · … | −48.1 / −49.5 · −29.2 / −29.5 · … |
| schematic + Tomkin's 1.8 Hz (`analog_results.md` §3.5) | 4234 | −58.8 / −56.6 · −72.1 / −44.6 | −65.4 / −57.6 · −59.3 / −35.8 | −51.4 / −50.1 · −39.1 / −33.7 |
| FIR 2300 (87 ms, the longest window) | 18634 | −58.4 / −56.5 · −69.4 / −58.5 | −63.3 / −61.2 · −60.1 / −50.0 | −52.9 / −52.7 · −51.3 / −45.6 |
| *repeat capture (AMPL 60 from the group's highest)* | | *−53.8 / −52.1 · −66.8 / −55.5 · −41.8 / −40.5* | *−66.0 / −66.3 · −56.6 / −46.4 · −55.9 / −55.6* | (no second cluster capture) |
| S fitted on M and L only (stored output kernels) | 4200 | −43.3 / −43.1 · −63.9 / −47.6 · −31.0 / −32.5 | | |

(`analysis/chain_kernel/fit.log` has the full set: FIR 300 / 800 / 1200, S weight 1, the slow term as a DC offset only.)

- **Excess over the repeat-capture null (100 Hz–10 kHz):** S −5.0 / −5.2 dB (the kernel nulls deeper than another, noisier capture does), M +0.4 / +4.4 dB; at 10–12 kHz S −5.2 / −5.8, M −0.8 / +3.0. L has no repeat cluster capture.
- **The FIR length:** with the high-pass prior the FIR decays to −47 / −53 dB (L / R) at 300–600 lags, and its DC gain stays at 5.1–5.4 × its peak (old kernel 4.7). Longer FIRs score a little better in the R windows' low band, but their long lags are only weakly identified next to the slow term (energy −22 dB at 600–1200 lags in L, DC gain swinging between −6.7 and +3.9 × the peak): unphysical, so not used. 300 lags is too short (M −53, L −49 dB).
- **The high-pass:** the schematic's corners alone (4.46 Hz) score best in the low band. Adding Tomkin's 1.8 Hz (the 6.25 Hz that `analog_results.md` §3.5 fitted with a high-pass-only model on 20 ms windows) is 4.8–12.5 dB worse in the R windows' low band once the FIR is free (L within ±1.6 dB); half the schematic + 1.8 Hz (3.5 Hz) is as good as the schematic. So the data fix the total near 3.5–4.5 Hz, not its split (**S** for the total, **H** for the split).
- **Cross-program:** S from M and L only stays at −43 dB, held back by the stored M/L output kernels (with refitted ones S and the M/L kernels can't be separated). The even/odd split is the held-out test.

### 12.2 The kernel (`write` stage: all AMPL ≥ 30 captures)

R's mid-band gain and delay are pinned to the library kernel at 300–3000 Hz (+0.18 dB and −0.001 samples removed), so levels and timing in the caches keep their convention. L keeps its own relative to R: **+0.10 dB, −0.008 samples (96 kHz)**. Response relative to the library kernel, L / R (`analysis/chain_kernel/chain2.json`):

| | 25 Hz | 50 | 100 | 160 | 280 | 400 | 1 k | 4 k | 8 k | 11 k | 11.5 k | 12 kHz |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dB | +1.12 / +0.83 | +1.18 / +1.04 | +0.36 / +0.21 | −0.16 / −0.25 | +0.20 / +0.13 | +0.02 / −0.09 | +0.10 / −0.00 | +0.10 / +0.02 | +0.05 / −0.01 | −0.07 / −0.02 | −0.13 / +0.01 | −0.36 / +0.04 |
| deg | +13.8 / +11.9 | +0.8 / +1.8 | −3.1 / −2.6 | −1.0 / −0.7 | +0.4 / +0.4 | 0.0 / −0.1 | +0.1 / 0.0 | +0.2 / 0.0 | +0.2 / −0.1 | −0.6 / 0.0 | −1.3 / +0.1 | −1.1 / +0.3 |

It holds exactly the §10 ripple (−0.2 dB at 160 Hz, +0.13 at 280 Hz), the 25–100 Hz excess of §11.2, and an L channel 0.36 dB down at 12 kHz (§11.1: S's L channel 4–8 dB worse there). The refitted output kernels: M 0.955 centre, +0.024 / +0.022 at ±1; L 0.782, +0.122 / +0.115 at ±1 (stored: 0.952 / 0.783 and similar sides; the outer taps differ most). In the deconvolution (`tools/harness/core.py`, `Kernel` with `lead` and `hp_hz`) each channel uses its own kernel, and a kernel with high-pass corners also gets a zero-phase low band edge at 3–12 Hz, so nothing below 3 Hz is inverted. Regularisation costs ≤ 0.02 dB below 12 kHz for both kernels, so `outofband`'s P(f) (the band window) still holds. The library kernel deconvolves bit-identically as before.

### 12.3 Validation, first pass (`firstpass_gap.py --stage collect|chain`, `EMT_KERNEL=chain2`)

ml_lines.Model as in §10 (clock and shift per capture, two gains per channel); first pass held out (nothing in the kernel fit or the correction saw a reverb sample). Medians over all captures, L / R, dB:

| IRs, model | Free numbers | S cluster/gap | S first pass | M cluster/gap | M first pass | L cluster/gap | L first pass |
|---|---|---|---|---|---|---|---|
| old, none | 0 | −49.4 / −46.3 | −44.7 / −45.7 | −48.5 / −45.7 | −44.5 / −44.1 | −44.9 / −45.5 | −44.3 / −42.7 |
| old + `cluster_perch` (§10.3) | 866 | −61.6 / −59.7 | −58.3 / −57.8 | −63.5 / −64.1 | −53.9 / −49.7 | −57.1 / −57.1 | −49.5 / −46.5 |
| **new, none** | 0 (kernel: 4200) | **−54.4 / −54.2** | **−52.0 / −53.3** | **−54.0 / −53.6** | **−50.0 / −48.2** | **−52.3 / −52.5** | **−47.5 / −45.6** |
| new, none, refitted output kernels (S: unit tap) | 34 | −58.1 / −58.0 | −54.1 / −56.4 | −57.5 / −57.0 | −51.5 / −48.9 | −51.5 / −51.1 | −47.3 / −45.4 |
| new + a correction refitted on it | 866 | −62.1 / −60.0 | −56.5 / −58.0 | −63.6 / −63.3 | −52.4 / −49.6 | −56.4 / −58.0 | −49.0 / −46.5 |
| *repeat capture (new / new, refitted kernels)* | | | *−54.7 / −63.6 · −55.8 / −63.7* | | *−49.9 / −49.1 · −50.9 / −49.1* | | *−50.9 / −54.1 · −50.4 / −54.0* |

- **Without any correction the new IRs gain 7.3 / 7.6 (S), 5.5 / 4.1 (M), 3.2 / 2.9 (L) dB in the first pass.** M is at the repeat-capture level (excess −0.1 / +0.9 dB), L +3.4 / +8.5, S +2.7 / +10.3 (+1.7 / +7.3 with the refitted kernels).
- **A correction refitted on the new IRs is tiny:** ≤ 0.04 dB and ≤ 0.4° in every band from 100 Hz to 10 kHz (was 0.25 dB / 2.4°). With the refitted output kernels it gains only 0.6–1.2 dB in the first pass. The first-pass transfer error H = Y/T on S's b = 0 captures (`report` stage) is now within ±0.011 dB and ±0.02° in R in every band, and in L within 0.04 dB / 0.5° (100–200 Hz, at the repeat level there); it was −0.11…−0.14 dB / −1° at 100–200 Hz.
- **What stops S at −54…−56 instead of −58:** model parts fitted on the old IRs. S's stored output kernel is a unit tap with ±0.002 side taps fitted on the old IRs: swapping it for a unit tap gains 2–4 dB. The rest is a setting-independent residual waveform (the mean residual of the even settings still removes 4 / 8 dB on the odd ones, L / R) that no FIR on the tap part removes (±48 taps: 1.8 / 3.2 dB), 44 % of it at 7–10 kHz in R. So it is no longer a transfer error of the pipeline; the stored first-pass tap parameters (positions, signs, b steps; fitted on the old IRs) are the candidates (**H**). *Ruled out (§12.7): refitted on the new IRs they don't move.*
- **Repeat references moved** because the repeat model includes the model's cluster column (`ml_lines.fit_repeat`): S L −48.2 → −54.7.

### 12.4 Validation, compact S model (`q2_injection.py --stage final`, `EMT_KERNEL=chain2`, no `Q2_CHAIN`)

`analysis/q2/final_s_chain2.{json,log}`; nothing in the §17 model refitted (9 + 13 + 3 numbers and the AP230 table, all fitted on the old IRs), per capture clock, shift and gain. All 128 settings, L / R [excess R over the repeat-capture null of `nulls_free`, which is from the old IRs]:

| Window (ms) | old (§17) | old, corrected (§19) | **new** | Flat, corrected | **Flat, new** |
|---|---|---|---|---|---|
| 0–50 | −42.9 / −44.3 | −46.2 / −46.3 | **−44.8 / −47.1** | −46.8 / −46.2 | **−43.7 / −50.1** |
| 50–150 | −41.0 / −42.3 [+13.5] | −44.9 / −44.9 [+11.0] | **−48.7 / −49.1 [+5.9]** | −45.6 / −45.7 [+5.9] | **−50.5 / −50.4 [+0.6]** |
| 150–300 | −39.6 / −38.8 [+16.1] | −43.4 / −42.6 [+13.2] | **−44.3 / −44.6 [+10.8]** | −44.5 / −46.8 [+6.6] | **−49.4 / −49.1 [+3.2]** |
| 300–600 | −36.4 / −35.3 [+13.0] | −38.7 / −40.2 [+10.4] | **−39.2 / −39.7 [+9.2]** | −42.2 / −45.3 [+0.5] | **−46.1 / −45.3 [+0.6]** |
| 600–1000 | −32.0 / −32.6 [+7.1] | −34.5 / −35.5 [+4.8] | **−35.0 / −35.2 [+4.4]** | −37.2 / −37.5 | −37.6 / −37.4 |

- **The new IRs beat the corrected old ones** by 3.8 / 4.2 dB at 50–150 ms, 0.9 / 2.0 at 150–300, and match them from 300 ms on (L +0.5, R −0.5 dB). On flat settings: 4.7–4.9 dB better at 50–150 ms, 4.9 / 2.3 at 150–300, 3.9 / 0.0 at 300–600, equal after.
- **0–50 ms in L is 1.4 dB worse** (flat: 3.1 dB; R 3.9 dB better). The likely reason is that the per-channel kernel takes out the chain's 0.10 dB L/R level difference, which the 8 cluster gains (fitted on the old IRs, per channel) had absorbed. 0.10 dB alone is a −39 dB error (**H**; refitting the 8 gains on the new IRs would show it).
- By control (R excess, 50–150 / 150–300 / 300–600 ms): HF 0.2 +8.8 / +15.2 / +13.9, LF 0.1 +12.7 / +17.6 / +17.1, RT 0.6 +14.5 / +17.6 / +10.4, but LF 1.0 +3.4 / +6.1 / +3.2, RT 2.2–2.5 +2.1…+7.0. The remainder is where the loop laws are weakest, as §19 found.
- CLA unchanged (its own IR, not deconvolved here).

### 12.5 Validation, out of band (`outofband.py --stage score|report`, `EMT_KERNEL=chain2`)

`analysis/outofband/report_chain2.{json,log}` (stored output kernels) and `report_chain2_hk.{json,log}` (refitted). Strict model nulls L / R, dB [repeat capture]:

- **Level at onset (capture re model):** 25–50 Hz +1.18 / +1.23 / +1.17 → **+0.16 / −0.02 / +0.05 dB** (S / M / L); 50–100 Hz +0.53 / +0.68 / +0.68 → **+0.03 / +0.06 / +0.07**; 100–200 Hz −0.15…−0.09 → +0.02…+0.04.
- **Low-band nulls now test the model**, S:

  | S | 25–50 Hz, old | **25–50 Hz, new** | 50–100 Hz, old | **50–100 Hz, new** |
  |---|---|---|---|---|
  | first pass | −16.4 / −14.9 [−17.7 / −31.3] | **−25.9 / −28.2** [−29.4 / −35.8] | −21.9 / −23.1 [−25.9 / −36.0] | **−37.3 / −41.7** [−36.1 / −45.4] |
  | recirc 1 | −16.1 / −14.3 [−24.9 / −38.6] | **−25.8 / −27.3** [−40.2 / −43.3] | −20.1 / −20.4 [−51.2 / −56.1] | **−35.3 / −37.9** [−53.9 / −58.1] |
  | recirc 2 | | −33.3 / −31.3 [−51.4 / −47.9] | | −32.6 / −40.8 [−62.0 / −60.8] |

  S's first pass is at the repeat-capture level at 50–100 Hz; its recirculations stay 15–20 dB above it. That is the model's low band (the LF laws, the DC floor, β; S LF 0.1), now measurable (**S**). The model with the high-passes added is now 10–20 dB *worse* (the IRs no longer carry them). M and L: first pass at the repeat level in both low bands (M 25–50 Hz −24.5 / −26.2 [−24.8 / −26.0], L −25.1 / −25.8 [−21.6 / −24.5]); their recirculations are time-variance-limited (model ≈ repeat).
- **Decay, 25–100 Hz:** S's RT ratio 1.007 / 1.008 → 0.991 / 1.003 (25–50 / 50–100 Hz). M and L unchanged (0.990 / 0.995, 0.963 / 0.972): the sweep artefact of §11.2 stays.
- **S's L channel at 11–12 kHz now matches R:** cluster/gap −32.8 / −43.4 → **−41.3 / −41.7**, first pass −32.7 / −43.0 → **−40.1 / −41.4** [repeat −39.9 / −40.5]; 10–11 kHz L −41.2 → −43.5. **S** (the L-channel chain was the cause).
- **With the stored output kernels, R loses a little at 10–11 kHz** (S first pass −47.4 → −44.4; M and L also lose 1–8 dB at 10–12 kHz): those kernels were fitted on the old IRs. **With the refitted ones (`ML_HK`, 34 numbers, cluster windows only) the high band improves everywhere**, first pass L / R [repeat]:

  | | 5–10 kHz | 10–11 kHz | 11–12 kHz |
  |---|---|---|---|
  | S, old (§11.1) | −50.0 / −52.1 | −41.2 / −47.4 | −32.7 / −43.0 |
  | **S, chain2 + refitted** | **−54.9 / −56.1** [−57.6 / −58.6] | **−51.1 / −49.7** [−48.6 / −47.6] | **−42.5 / −42.5** [−40.8 / −40.5] |
  | M, old | −51.3 / −48.9 | −47.1 / −39.4 | −40.8 / −33.1 |
  | **M, chain2 + refitted** | **−53.5 / −51.0** [−59.5 / −66.3] | **−50.8 / −46.9** [−56.7 / −60.7] | **−47.1 / −45.5** [−51.2 / −54.4] |
  | L, old | −46.0 / −42.3 | −31.3 / −29.2 | −21.3 / −19.8 |
  | **L, chain2 + refitted** | **−44.8 / −42.4** [−49.3 / −65.3] | **−40.2 / −38.7** [−44.4 / −55.8] | **−35.7 / −35.0** [−38.4 / −48.2] |

  So most of §11.1's L-program phase error above 9 kHz (cluster window −21 dB at 11–12 kHz) was the old IRs plus a kernel fitted on them: refitted with the chain kernel, L's 11–12 kHz first pass gains 14–15 dB. S at 10–12 kHz is now at or below the repeat-capture null. `analysis/outofband/report_chain2_hk.{json,log}`.

### 12.6 What this settles, and what's left

- **The pipeline:** the device-domain IRs no longer need `chain.json:cluster_perch` or `Q2_CHAIN`. Use `EMT_KERNEL=chain2` for new work; the old caches stay for comparison. **S**
- **Unit facts:** the library chain's sub-20 Hz high-pass totals about 3.5–4.5 Hz in this fit (the schematic's 4.46 Hz alone fits best), and L's chain is 0.10 dB above R's with 0.36 dB less at 12 kHz.
- **Left (model side, not the chain):**
  1. ~~refit the parameters that were fitted on the old IRs: the 8 cluster gains (L/R), S's output kernel (a unit tap), the M/L output kernels (refitted here), the first-pass taps and the b steps;~~ **done (§12.7, `model_s_results.md` §21)**: the cluster gains become one (0–50 ms +3.9 dB in L); everything else holds;
  2. adopt the refitted M/L output kernels in the model (`ML_HK`; 10–12 kHz first pass +5…+15 dB, §12.5) after checking them on the recirculations;
  3. S's first-pass residual waveform at 7–10 kHz (§12.3);
  4. S's low band in the recirculations (§12.5);
  5. the repeat-capture references (`nulls_free`) recomputed on the new IRs.

### 12.7 The first-pass parameters refitted on the new IRs (`ml_carryover.py --stage firstpass`, `EMT_KERNEL=chain2`) **[S]**

`EMT_KERNEL=chain2 .venv/bin/python tools/ml_carryover.py --stage firstpass` → `analysis/ml_carryover/firstpass_chain2.{json,log}`. `tools/ml_lines.py` (and so `firstpass_gap.py` and `outofband.py`) reads it with `ML_FP=analysis/ml_carryover/firstpass_chain2.json`; results get `_fp`. Same training settings and procedure as §7.6: positions and signs from 6 training settings, the output kernel from their clusters, b per setting.

> **Summary.** Refitted on chain2, the first-pass parameters are the old ones. In S, M and L the integer tap positions and every sign are identical (fractional parts within 0.002 samples); every b is within 0.0008 of the old value; S's output kernel comes out as a unit tap (side taps ≤ 0.0016). So §12.3's last candidate for S's residual first-pass waveform (the stored tap parameters) is ruled out (**S**). That residual is a property of the model's structure or of the library, not of a stale fit.

| | Positions (rounded) | Signs | b, max change | Output kernel |
|---|---|---|---|---|
| S | identical (max fractional change 0.002) | identical | 0.0004 | unit tap, side taps ≤ 0.0016 (old ≤ 0.0022) |
| M | identical (0.002) | identical | 0.0005 | 0.952 centre, +0.021 / +0.027 at ∓1; changes ≤ 0.0024 |
| L | identical (0.002) | identical | 0.0008 | 0.783 centre, +0.112 / +0.121; changes ≤ 0.0017 |

- **The M/L output kernels agree with §12.2's** (refitted with the chain kernel, clock per capture) in magnitude to 0.02 dB up to 10 kHz (0.10 / 0.19 dB at 12 kHz, M / L). They differ in their asymmetry (M ±1 taps 0.021 / 0.027 here against 0.024 / 0.022 there), i.e. in a fractional delay. This bench fits no capture clock, so §12.2's kernels are the better-identified ones: use `ML_HK=analysis/chain_kernel/chain2.json`. With S as a unit tap and positions and b unchanged, §12.3's row "new, none, refitted output kernels" is the held-out validation of this refit (S −54.1 / −56.4, M −51.5 / −48.9, L −47.3 / −45.4 dB first pass).
- **This bench's own held-out nulls** (no clock fit, so below the `lines` bench; first pass, b fitted, L / R, old IRs → chain2): S −42.8 / −41.5 → −45.0 / −43.0, M −40.7 / −39.2 → −41.6 / −39.9, L −36.6 / −35.9 → −36.8 / −36.1 dB. Cluster and gap: S −49.1 / −44.1 → −52.9 / −48.8, M −46.6 / −42.0 → −50.9 / −45.1, L −43.0 / −38.8 → −46.0 / −40.4 dB.
- **The b table is S's four steps in every program** (0.001 / 0.17 / 0.33 / 0.50), unchanged; so is the compact S model's AP230 table (`model_s_results.md` §21).


## 13. chain3: the chain kernel's long tail (October 2026) **[S]**

`.venv/bin/python tools/chain_kernel.py --stage tail` → `analysis/chain_kernel/tail.{json,log}`; `--stage write --name chain3` → `analysis/library/tomkin96/kernels/chain3.npz`, `analysis/chain_kernel/chain3.json`. Caches: `EMT_KERNEL=chain3` (`ml_structure.IRs.build` for S / M / L, then `q1_generalise.py --stage data`). Motivated by `model_s_results.md` §23: S's reverb against the S baseline on chain2 showed a fixed transfer error below 200 Hz whose L part peaks at lag 501, just past chain2's last FIR lag.

> **Summary.** chain2's FIR stops at lag 500, and its per-capture slow term {1, t, t²} absorbs whatever lies beyond: that is why §12.1 found longer FIRs "unphysical". With the slow term as a **DC offset only** (passed through the known high-pass inversion, the physical form of a capture-side offset) and the lags from 300 to 1000 as a **smooth cubic B-spline tail** (knots every 40 samples), the kernel fitted on the cluster-only windows alone (no reverb, no model) **reproduces the transfer error measured independently over S's reverb**: 86 % / 58 % of its energy (L / R) at 20–200 Hz. chain3 has 2,768 numbers (chain2: 4,200). Held out, its cluster windows gain 6–12 dB at 25–100 Hz in M and L and lose ≤ 0.3 dB anywhere in S. On the rebuilt IRs, with nothing in any model refitted, M's and L's first pass gain 1.8–4.3 dB in L at 50–200 Hz, and S's reverb loses 82 % / 85 % of the transfer error's energy (L / R). **Use `EMT_KERNEL=chain3` for new device-domain work.**

### 13.1 Which tail form (`tail` stage)

Fitted on all captures; the new kernel's response over chain2's (gain and delay matched at 300–3000 Hz) against T(f), the capture / model transfer over S's reverb on chain2 (`model_s_results.md` §23.2; rms |T − 1| at 20–200 Hz L 0.0206, R 0.0151). Residual = rms |ratio − T|:

| Slow term, FIR | Numbers | L: change / residual | R: change / residual |
|---|---|---|---|
| quad, 500 (chain2) | 4200 | 0 / 0.0206 | 0 / 0.0151 |
| quad, L 800 free | 5400 | 0.082 / 0.069 | |
| quad, L 1000, tail from 300 | 3484 | 0.133 / 0.127 | |
| dc, 500 | 4200 | 0.003 / 0.0196 | 0.018 / 0.018 |
| dc, L 800 / 1000 / 1200, tail from 300 (knots 40 or 80) | 3448–3504 | 0.018–0.019 / **0.0077–0.0083** | |
| dc, L 1000 free | 6200 | 0.018 / 0.0079 | |
| dc, L 1000, tail from 200 / 400 | 3092 / 3872 | 0.028 / 0.0206, 0.024 / 0.0163 | |
| **chain3: dc, both 1000, tail from 300 / 40** | **2768** | **0.018 / 0.0077** | **0.018 / 0.0098** |

- **With the quadratic slow term the long lags run away** (residual 3–6 × T itself): the quadratic and the tail trade off, as §12.1 saw. With the DC term the result doesn't depend on the FIR length (800–1200) or the knot spacing, and a free 1000-lag FIR gives the same as the spline (the spline only removes numbers).
- The match, chain3 / chain2 against T (dB / deg), L: 25 Hz +0.29 / −1.8 against +0.31 / −2.5; 30 Hz +0.06 / −2.3 against +0.05 / −2.2; 35 Hz −0.13 / −1.9 against −0.15 / −1.3; 40 Hz −0.23 / −0.9 against −0.19 / −1.3; 50 Hz −0.15 / +0.8 against −0.08 / +0.8; 70 Hz +0.15 / −0.2 against +0.14 / −1.2. Correlation over 20–200 Hz +0.93.
- **Two independent data sets agree**: the cluster-only windows (no reverb sample) and S's reverb (through the S model). So the tail is the chain's, not a model error fitted into the kernel.
- Tail energy (re total): L −46.6 dB at lags 300–600, −54.7 at 600–1200; R −45.4 / −51.1. chain3's R mid-band gain is pinned to the library kernel as chain2's (+0.244 dB removed); L / R difference +0.099 dB, −0.008 samples, as chain2.

### 13.2 Held out on the cluster windows (fit on even settings, score odd)

Median nulls L / R, dB, 100 Hz–10 kHz | 25–100 Hz | 10–12 kHz, per-capture freedoms as §12.1 (slow term {1, t, t²} in the scoring):

| Kernel | S | M | L |
|---|---|---|---|
| chain2 | −58.8 / −57.3 · −72.1 / −57.1 · −47.0 / −46.3 | −65.6 / −61.9 · −59.3 / −41.8 · −56.7 / −52.6 | −51.3 / −50.7 · −37.5 / −38.5 · −39.1 / −38.4 |
| dc, 500 | −58.8 / −56.6 · −71.3 / −45.9 · −47.0 / −46.3 | −65.7 / −58.4 · −60.0 / −36.2 · −56.6 / −52.5 | −51.1 / −50.2 · −39.3 / −34.3 · −39.1 / −38.4 |
| **chain3** | −58.5 / −56.8 · −71.8 / −58.7 · −46.7 / −46.1 | −64.9 / **−64.1** · −60.5 / **−50.5** · −56.0 / **−55.0** | **−52.0 / −51.3 · −49.3 / −43.7** · −39.3 / −38.9 |

With a DC-only slow term in the scoring the ranking is the same (`tail.log`). The DC-only slow term alone (dc, 500) is worse in R's low band; the tail is what helps.

### 13.3 On the rebuilt IRs (nothing in any model refitted)

- **S, M, L first pass** (`lowband.py --stage mlcheck`, ml_lines' law model with each kernel's refitted output kernels, `ML_HK=analysis/chain_kernel/chain3.json`; every 4th group, AMPL ≥ 30; median change per capture L / R, dB; held out in time, no kernel fit sees a first-pass sample):

  | | S first pass | S recirc 1 | M first pass | L first pass |
  |---|---|---|---|---|
  | 25–50 Hz | −0.9 / −4.8 | −2.8 / −6.9 | −0.7 / −0.9 | −1.1 / −1.8 |
  | 50–100 Hz | −1.9 / −3.2 | −8.6 / −6.8 | **−2.9** / 0.0 | **−2.8** / −0.1 |
  | 100–200 Hz | −5.8 / −1.6 | −4.6 / −1.3 | **−1.8** / −0.1 | **−4.3** / 0.0 |
  | 100 Hz–10 kHz | −1.9 / +0.1 | −0.3 / 0.0 | −1.4 / −0.1 | −0.4 / −0.1 |
  | cluster / gap, 25–200 Hz | +0.4…+1.4 | | −0.2…+3.4 | −0.0…+0.3 |

  M's and L's R channels don't move: their R first pass is at −32…−34 dB at 50–100 Hz, held by something else there (open). The cluster / gap windows (9–86 ms) barely resolve the low bands; their strict nulls lose a little because chain3 no longer fits a free quadratic per capture.
- **S's reverb**: `model_s_results.md` §24.

### 13.4 What this settles

- **The chain's response runs on past 19 ms with a smooth low-frequency tail** in both channels (−45…−47 dB at lags 300–600), and the library's capture-side slow term is a DC offset, not a free quadratic. **S**
- §12.1's "longer FIRs are unphysical" was the quadratic slow term, not the data.
- Left in S's reverb after chain3: a phase of −1.1…−1.4° at 25 Hz (both channels), ±0.07 dB above 30 Hz (`model_s_results.md` §24).

## 14. The pulse check on M and L (`tools/pulse_check.py --stage ml`, October 2026) **[S: the static model passes its static metrics; S: the apparent late-decay excess is Tomkin's longer fades in M / L; S: the rotation through a sweep does not change the late decay]**

`EMT_KERNEL=chain3 .venv/bin/python tools/pulse_check.py --stage ml --workers 8` → `analysis/pulse/ml/`, report in `analysis/pulse/report.log` (`model_s_results.md` §26 for S and the method). All 340 groups (M 160, 479 captures; L 180, 360 captures): the static model (§8's line law with S's HF pole curve, θ = 0) rebuilt per capture by `ml_lines.fitted` (clock, shift and two gains per channel fitted before any line; everything after predicted), against each capture and every capture pair. EDCs from the R onset + 20 ms; coherence, fine structure and echo density over the 250 ms after it; up to the capture's usable end. Median model [capture pair], R:

| Metric | M | L |
|---|---|---|
| octave decay (−5 / −25 dB), signed, 125 Hz–8 kHz | +0.3…+0.6 % [−0.2…+0.2] | +0.3…+1.0 % [+0.0…+0.3] |
| octave decay, \|error\| | 0.5–2.1 % [0.2–2.3] | 0.5–1.7 % [0.2–2.1] |
| EDT, signed | −0.3…+0.5 % | −0.1…+0.4 % |
| largest EDC difference to −30 dB | 0.18–0.90 dB [0.13–0.88] | 0.20–0.77 dB [0.13–0.98] |
| echo density \|difference\| | 0.021 [0.023] | 0.009 [0.009] |
| L/R coherence, model − capture | ≤ 0.005 [≤ 0.002] | ≤ 0.001 [≤ 0.001] |
| fine structure 100 Hz–1 kHz / 1–4 kHz / 4–10 kHz | 0.980 / 0.983 / 0.984 [0.970 / 0.977 / 0.982] | 1.000 / 1.000 / 1.000 [1.000] |
| late decay (−25 / −40 dB), signed, ≤ 1 kHz | +4.3…+6.0 % [−0.2…+0.4] | +5.2…+8.8 % [+0.4…+1.6] |

- **Early decay, EDCs, coherence, echo density and fine structure are at the capture-to-capture spread** (the EDC differences of 0.2–0.9 dB are the time variance: the pairs differ as much). Subgroups: short RT (M ≤ 1.0 s) +1.1…+1.7 % below 500 Hz, LF 0.1 in M +1.2…+2.5 % at 250 Hz–2 kHz, L at HF ≤ 0.2 +1.3…+2.1 % at 125–250 Hz: all within ±3 %, and the short-RT / LF 0.1 regime is where S's own errors sit too.
- **Decorrelation (null per window, model [pair]):** M first pass −47.7 [−48.3], recirc1 −20.8 [−19.3], recirc2 −16.2 [−13.8], recirc3 −12.7 [−10.2], 0.3–0.6 s −12.9 [−10.5], 0.6–1 s −9.9 [−7.8], 1–2 s −7.3 [−6.2]; L first pass −45.3 [−51.9], recirc1 −18.9 [−17.5], recirc2 −15.9 [−13.5], recirc3 −12.7 [−10.2], 0.3–0.6 s −17.6 [−15.9], 0.6–1 s −13.0 [−10.4], 1–2 s −8.9 [−6.5]. The static model sits 1.5–3 dB closer to each capture than the captures are to each other, as §9.5's periodic θ predicts (a pair differs by twice a capture's own time variance); the decorrelation curve itself is §9.5's sim (rotating model through a sweep), not rerun here. L's first pass is 6.6 dB short of the pair value (−45.3 against −51.9; M 0.6 dB): the first-pass L gap of §7.6 / §12.
- **The late-decay excess is the file fade, not the model:**
  - Through §9.5's simulated sweep, the rotating model's late decay (−25 / −40 dB) differs from the static model's by only −5.5…+4.3 % per octave, scattered around zero (mean ≈ −0.3 %; 8 SIM_KEYS settings × 4 renders), while model vs capture reads +2…+24 % at the same settings. The rotation does not explain it.
  - Kept further from the file end, it goes away: late model / capture at 1 kHz with the fit ending 0.15 / 0.6 / 1.2 s before the usable end: M +2.9 / +0.7 / −0.8 %, L +5.6 / +5.4 / +4.6 %; the slope of the capture / model energy ratio, frames ending 1 / 2.5 / 4 s before the last sample: **L 250 Hz +6.1 / +0.7 / −0.0 %, 1 kHz +4.7 / +0.6 / +0.0 %**, M ±0.5…0.9 % at 1 s.
  - The ratio against time before the end (250 Hz–4 kHz, median): M falls from ~1.0 s before the last sample (−2.2 dB at 0.7 s, where a 0.70 s fade would give 0), L slopes from ~4 s (−0.5 dB at 2 s, −4.2 dB at 1 s). **M's fade is ~1.0 s and L's ~3–4 s, not the 0.72 s `ml_lines.usable_end` assumes** (`library_notes.md` §2). Away from the fades the static model's late decay matches to under 1 %.
- **Consequences:** (1) §8.3's per-octave RT reading (from the onset to `usable_end`) is partly fade-biased in L, so the "RT shortfall" of 0.97 at 125 Hz there needs re-reading with program-specific fade ends before anything is fitted to it; §9.5's sweep effect (L 0.958–0.992) is independent of the fade. (2) `ml_lines`' time windows that end near a file's end (1–2 s, 2–4 s in long-RT L groups) carry fade bias; the recirculation windows don't.

## 15. Tomkin's fade in M and L, and the usable end (`tools/ml_fade.py`, October 2026) **[S: M a linear 1.0 s ramp; S: L a linear ramp of per-file length, about 0.3 × the file; H: a level-guided trim per file]**

`EMT_KERNEL=chain3 .venv/bin/python tools/ml_fade.py --stage measure|report|rt` → `analysis/ml_fade/` (`curves/` one file per AMPL group, `fits.json`, `fade.json`, `report.log`, `rt/`, `rt.log`). All 340 groups (M 160, L 180). Per capture the static model of §8 (`ml_lines.fitted`: stored clock and gains, nothing refitted) is the unfaded reference: octave energies 250 Hz–4 kHz of capture and model in 25 ms frames from the R onset + 0.1 s to the file's end (the cached IR's end is the raw file's end to one device sample), frames above the floor only (model within 90 dB of its first frame, capture within 20 dB of its early level against the model). The level difference d(t) (mean over octaves) is fitted per file and channel as a gain plus 10 log10 of the fade's mean g² over each frame. The gain absorbs everything else; a free slope (the static model's late decay is within ~1 %) trades against the fade length on short files and is not used.

**Result.**
- **Shape: linear in amplitude**, as S's (§15.2 of `model_s_results.md`). With the ramp's start and its zero point both free, linear fits to 0.37 / 0.44 dB rms per frame (M / L, median), below a ramp forced to end at the file's end (0.55 / 0.52), quarter-sine (0.50 / 0.53), raised cosine (1.5) and dB-linear (2.0–2.8). A free power g = u^p gives p 0.82–0.98 (M) / 0.97 (L).
- **M: a 1.0 s ramp.** The fade starts 0.98 / 1.01 / 1.09 s before the file's end (quartiles, 556 file-channels where it is seen). Ramps seen whole are 0.97–1.03 s long (149); the longer fits (1.1–3 s, 260) are 1.0 s ramps the trim cut short (the ramp's zero point lies up to 0.3 s after the file's end; T and that point trade on a partial ramp, the start does not). Some files under 1.5 s carry a 0.3 s ramp instead (33 file-channels at 0.25–0.35 s).
- **L: a linear ramp whose length follows the file:** the start lies 0.31 × the file's length before its end (IQR 0.25–0.40, 10–90 % 0.19–0.57; T 0.8 / 1.6 / 3.3 s, up to ~12 s on 40 s files; Spearman with length +0.83). No fixed length and no fixed fraction: the ratio scatters as much between the captures of one setting as between settings, so it is per file.
- **Looking for the rule (H):** re the reverb 0.1–0.2 s after the peak, the full-band level where the fade starts is −28…−24 dB (M) and −27…−22 dB (L) (IQR; peak-referenced it spreads twice as wide because the cluster follows AMPL), with an LF trend in L (−27 dB at LF 0.1 to −21 dB at LF 2.0) that a full-band trim would produce. So the fades look hand-set per file near a level about 25 dB down. That is ±2.5 dB, about ±1 s at RT 20 s: not precise enough to place a fade, so the per-file measurement is used.
- **Precision:** the start agrees between L and R to −0.01…+0.01 s (IQR; 90 % within 0.05 s). Against an AMPL partner that is still unfaded there (no model): 12 M pairs agree to −0.002 / +0.002 / +0.024 s, 2 L pairs to +0.08 s.
- **Where the fade is not seen** (272 M / 218 L file-channels, mostly long RT): the reverb meets the floor (the model 80–120 dB down) before the fade, and the frames show no fade up to that point.

**`ml_lines.usable_end(y, f, prog)`** now reads `analysis/ml_fade/fade.json`: per file the fade start (mean over the channels that see it) or, where it is not seen, the latest floor point, plus a 0.05 s guard; with no usable frames at all M 1.1 s, L 0.5 × the file's length. S and calls without a file keep S's 0.72 s. All callers pass the file (`ml_lines` law stage and `fitted`, `lowband`, `outofband`, and through `fitted` `pulse_check` and `ml_variance`). The usable end moves earlier by 0.32 / 0.41 / 0.46 s (M) and 0.13 / 1.04 / 3.12 s (L) (quartiles).

**§8.3 re-read** (`rt` stage: `ml_lines.band_decay` from the onset to the old and the new usable end, current model with the HF pole curve; RT / law, median over captures and channels):

| | 125 | 250 | 500 | 1 k | 2 k | 4 k | 8 k Hz |
|---|---|---|---|---|---|---|---|
| M, LF 0.1 / 0.7 / 1.2 / 2.0, new | 0.997 / 0.999 / 0.992 / 0.994 | 0.994 / 0.997 / 0.996 / 1.003 | 0.991 / 1.001 / 1.000 / 1.005 | 0.985 / 0.997 / 0.996 / 0.999 | 0.988 / 0.998 / 0.998 / 1.000 | 0.992 / 0.999 / 0.997 / 1.000 | 0.995 / 0.999 / 0.999 / 1.000 |
| L, LF 0.1 / 0.8 / 1.2 / 2.0, old | 0.999 / 0.982 / 0.971 / **0.957** | 0.997 / 0.995 / 0.989 / 0.985 | 0.989 / 0.993 / 0.992 / 0.993 | 0.995 / 0.996 / 0.995 / 0.996 | 0.993 / 0.997 / 0.997 / 0.998 | 0.997 / 0.998 / 0.997 / 0.999 | 0.999 / 1.000 / 0.999 / 1.000 |
| L, LF 0.1 / 0.8 / 1.2 / 2.0, new | 0.999 / 0.983 / 0.976 / **0.976** | 0.997 / 0.997 / 0.992 / 0.993 | 0.989 / 0.994 / 0.993 / 0.996 | 0.992 / 0.997 / 0.996 / 0.997 | 0.992 / 0.997 / 0.997 / 0.998 | 0.996 / 0.998 / 0.997 / 0.999 | 0.999 / 1.000 / 0.999 / 1.000 |

M moves by at most 0.002 (its old 0.72 s guard was close to its 1.0 s fade). **L's 125 Hz "shortfall" was half fade:** LF 2.0 0.957 → 0.976, LF 1.2 0.971 → 0.976; 250 Hz LF 2.0 0.985 → 0.993. What is left (2–2.5 % at 125 Hz, LF ≥ 1.2) lies inside §9.5's sweep effect (L 0.958–0.992) and the low-band unit transfer still open (`NEXT_STEPS.md` §0 item 2); nothing should be fitted to it until those are separated.

**Checks.** The pulse gate (`pulse_check.py --stage gate --tag mlfade`) passes: S, CLA and the Altiverb unchanged; M / L nulls unchanged; M / L early decay 1.12 / 0.61 % (was 1.06 / 0.66). On the gate's 8 M / L groups the late decay (−25 / −40 dB, model vs capture, R) moves from M +0.8…+10.8 % to −2.4…+2.0 %, L +0.7…+11.9 % to −0.9…+1.0 % above 125 Hz (L 125 Hz −5.1 % from few points).

**Full pulse rerun** (`pulse_check.py --stage ml`, all 340 groups, the old run kept in git history; `analysis/pulse/report.log`), model [capture pair], R, medians: **late decay (−25 / −40 dB) M −0.3…+0.9 % [−0.7…0.0], L 0.0…+0.5 % [−0.3…+1.0]** per octave 125 Hz–8 kHz (was M +4.3…+6.0 %, L +5.2…+8.8 % at ≤ 1 kHz, §14); early decay (−5 / −25 dB) M −0.3…+0.4 %, L 0.0…+0.3 % (was +0.3…+0.6 / +0.3…+1.0); EDC, coherence, echo density and fine structure as §14; nulls unchanged up to 1–2 s, M 2–4 s −4.9 [−6.2], L −5.4 [−3.5] dB. **The static model's decay now matches the library to within 1 % in every octave and window, at the capture-to-capture spread.**

**Not redone:** the `ml_lines` law files (`analysis/ml_lines/law/`) keep their old `end` and band decays; their recirculation nulls end long before any fade and stand. The law stage's `fit` (§8.3's per-LF law numbers) reads those band decays, so it needs the law stage rerun (or the `rt` stage's decays) before it is quoted again.

## 16. The HF floor in M and L (`model_s_results.md` §29.4)

`ml_lines.line_coefs` now floors the HF control at max(HF, 0.12 s × scale / RT), scale 2 (M) / 4 (L) — S's floor in program time (`ML_HFRT=0` switches it off). In M it is the best value at every RT 0.6 setting with HF 0.1 / 0.2 (+2…+3 dB residual; pulse check recirculation 1 −27.2 / −27.6 → −29.0 / −29.7 dB at HF 0.1 / 0.2); elsewhere it changes nothing. In L the line-length poles meet the 0.8475 ceiling first, so it never acts. L's RT 2.2 residual (−12 dB at HF 0.1 against −27 dB at RT 2.6) is not the HF pole's cap; it is that step's own HF floor of 0.40 (§8.6, fixed).

## 17. Reverb-specific metrics on M and L, and the Altiverb MED / LONG presets (G6; `tools/altiverb_ml.py`, October 2026) **[S: the static M / L model matches the library in every G6 metric at the capture-pair level; S: against another unit (the Altiverb MED / LONG presets, a frozen θ each) it matches as well as against Tomkin's own captures in the reverb, and the first pass sits at a cross-unit floor of −35…−36 dB]**

```
EMT_KERNEL=chain3 .venv/bin/python tools/memguard.py --limit-gb 8 -- .venv/bin/python tools/altiverb_ml.py     # 110 s, 0.4 GB -> analysis/altiverb/ml_pulse.{json,log}
```

**Library (`pulse_check --stage ml`, static θ = 0 model vs capture [capture pairs], R, median):** C80 ≤ 0.02 dB per octave at 250 Hz – 8 kHz in M (L 0.01–0.04, 125 Hz 0.13 / 0.34); early level ≤ 0.02 dB above 200 Hz (L 100 Hz +0.65 dB); late tonal balance 0.01–0.13 dB [0.02–0.19]; EDR 90th percentile 6.9 / 4.5 % [7.7 / 5.4 %] (M / L; θ makes captures differ at fine resolution, and the model is inside that); tail coloration 0.00–0.01 dB; mixing time 2.7 % [2.7 %] (M, 185 ms), 0.8 % [1.4 %] (L, 615 ms); IACC 0.000–0.008 [0.001–0.010]. At L AMPL 0 the 0–80 ms window holds no sound (L's reverb starts 86 ms after tap 0 on R), so early metrics are not computed there; at L AMPL 80 they are ≤ 0.04 dB.

**The Altiverb MED (10) and LONG (7) presets** (`altiverb_ml.py`). The model is unchanged; per preset only the chain is fitted (11-tap FIR per channel on the cluster-only region, cluster and reverb gains, the clock refined around the plugin's +42.8 ppm), plus AP230's damping step b from S's four values, chosen on the first pass. The presets are not library steps, and the short ones use b > 0: MED RT 0.5 b 0.33, LONG RT 2.0 b 0.50, LONG RT 3.0 b 0.17; with b = 0 their first pass was −6…−14 dB. Settings MED LF 1.1 / HF 0.7, LONG LF 1.2 / HF 0.4 (library_notes §6, read from these decays, so the decay metrics are not held out).
- **Nulls (R):** cluster and gap −40.5 (M) / −39.6 dB (L), in-sample with the FIR. **First pass −35.2…−35.3 (M), −36.1…−36.2 dB (L) in every preset at b = 0** (LONG RT 2.0 −31.8, b fitted), against −48.0 / −45.6 dB for the library. A constant cross-unit floor, the counterpart of S's −35…−38 dB against CLA and the Altiverb (`model_s_results.md` §26–27: unit / chain differences, mostly below 500 Hz). Recirculation 1 / 2 / 3 −16.4…−25.5 / −13.5…−18.8 / −11.4…−14.0 dB (M), −15.7…−22.8 / −11.4…−17.5 / −8.1…−13.4 (L); library model vs capture −20.7 / −16.1 / −12.7 (M), −19.1 / −15.9 / −12.8 (L), capture pairs −19.0 / −13.4 / −9.9, −17.6 / −13.4 / −10.2. With θ ≠ 0 nothing deeper can be shown against one realisation; the model is as close to another unit's M / L as to Tomkin's.
- **Decay (not held out):** octave |error| 0.4–1.5 % (M), 0.3–1.0 % (L), largest EDC difference 0.20–0.88 dB, within 0.3 dB of the library's model-vs-capture values.
- **G6:** C80 0.03–0.21 dB per octave (M), 0.12–0.15 dB (L); tail coloration within 0.01 dB; mixing time 2.7 % (M) and 0.8 % (L), the library's values; IACC within 0.006 early and late; echo density |difference| 0.022 / 0.013 [library 0.020 / 0.009]; EDR 90th percentile 6.1 % (M) / 4.0 % (L) [library pairs 7.7 / 5.4 %]. Late level differs by 0.3–0.6 dB at 100–250 Hz [library pairs 0.07–0.19], the low-band unit transfer seen in S (`model_s_results.md` §27, §33).

