# Programs M and L: clusters, loops and time variance (Q3)

Everything below comes from `.venv/bin/python tools/ml_structure.py` on the Tomkin **96 kHz** set
(log: `analysis/ml_structure/run.log`; data: `analysis/ml_structure/summary.json`). Program S is run
through the same code as the control and the calibration. Labels as in `EMT246_unit_facts.md`
(**C** confirmed, **S** strong, **H** hypothesis).

*October 2026. Answers `NEXT_STEPS.md` Q3 and changes the picture of the "stochastic control" (§8 there).*

> **Summary.**
> - **Early section = program S's, time-scaled by 2 (M) and 4 (L).** Same 8-tap cluster with the same signs, the reverb onset and the L/R onset gap all scale (§1, §2).
> - **The late reverb is a different structure from S.** S is four parallel loops: the product of four combs predicts its tail to **−45 dB held out**. M and L follow one five-term recursion each, four taps symmetric about a centre plus one at twice the centre: L 10826/11384/11618/12176/23002, M 5416/5694/5812/6090/11506 (§3).
> - **M and L are time-variant; S is not.** Separate captures of the same M or L setting decorrelate progressively (L: −25 dB at 0.15–0.3 s → −3 dB at 2–4 s), flat across frequency. S's tails agree to −35…−47 dB. This is the first direct evidence of the "stochastic control" and its program dependence (§4).
> - **Consequence for the null target:** past the first recirculation, M and L can't be nulled waveform-for-waveform against one capture; two captures of the same unit don't null each other (§5).

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
- **Reading (H):** one recirculating delay line of 2N words (L 23 K, M 11.5 K, within the ~39 K budget) with five feedback taps, or the equivalent two lines of N with a coupling. The four mid taps cancel at DC (+ − − +), so at low frequencies the loop is the 2N tap alone; the mid taps shape the loop gain as −4β sin(ωp/2) sin(ωq/2) on the half-loop. With uniform decay that is a 2 × 2 rotation-like system in w = z^−N: modes decay by the same factor per N samples wherever |trace| ≤ 2√det.
- **The plateau (−18 to −22 dB) is not a missing lag.** The residual is flat across octaves (within about 2 dB from 125 Hz to 4 kHz, 3–8 dB higher in the 8 kHz octave next to the 10 kHz band edge; S: −56 to −43 dB), the same in-sample as held out, and still there later in the tail (L fitted 8–10.5 s, tested 10.5–12.5 s: −12.3 / −18.9 dB): a fixed fraction of the decaying reverb that no fixed recursion predicts. §4 shows why.

### 3.3 Lag precision

Sub-sample position of each lag FIR's peak (low-passed), over 12 settings each with HF ≥ 0.9 (M) or 1.0 (L), against S as the calibration:

| | Offset from the integer lag (mean, L / R) | Scatter across settings |
|---|---|---|
| S (4 loops 1321…1587, and the 4-fold sum 5758), LF 1.0 and HF 1.0 only | +0.04 / −0.03 (5758: +0.07 / +0.02) | 0.26 / 0.13 |
| M (5416, 5694, 5812, 6090; 11506) | −0.80 / −0.79 (11506: −0.75 / −0.74) | ≤ 0.04 |
| L (10826, 11384, 11618, 12176; 23002) | −0.82 / −0.83 (23002: −0.73 / −0.76) | ≤ 0.04 (23002 on L: 0.09) |

S sits within 0.05 sample of its integer lags (the remainder is the capture clock, `model_s_results.md` §3.1), so the method has no offset of its own. **M and L carry a common −0.8 sample on every lag**, stable to ±0.04 across 12 settings each (RT, LF and AMPL varied): each lag is d − 1 plus about 0.2 sample of group delay that all taps share. **H:** a filter in the common feedback path (e.g. a one-pole with p ≈ 0.17), and the integer lags are d − 1 (L 10825 … 23001). The sums between lags (§3.2) are exact either way, but the absolute integers are uncertain by one sample until the loop filter is modelled. With heavy HF damping (HF 0.1–0.2) the offsets spread to ±1 sample, as the loop low-pass adds frequency-dependent delay *[exploratory]*.

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

- **S: LTI within the capture.** Its pairs agree to −47 dB early and still −35 dB at 0.6–1 s, where the decaying level and the fades take over. The S difference rises about 5 dB per octave, the signature of a small timing difference between two captures, not of the unit. L and R agree within 0.3 dB, so the library shows **no L-only excess** in S (the CLA unit's L read-out jitter, `model_s_results.md` §13, is not seen here at this depth).
- **M and L: the reverb decorrelates between captures, progressively and flat across frequency** (**S**). Both agree early (M −27/−38 dB at 0.05–0.15 s, L −25/−33 dB at 0.15–0.3 s, before its first recirculation at 0.41 s), then lose coherence as the reverb recirculates until they are uncorrelated (−3 dB) by 2–4 s (L) or nearly so (M, −5 to −7 dB after 1 s).
- **It is per-capture randomness, not AMPL.** M has three AMPL values; all three pairings (60–10, 99–60, 99–10) decorrelate equally after 0.15 s.
- **Not delay modulation, not clock drift.** Both would grow about 6 dB per octave, as S's small difference does; M's and L's decorrelation is flat from 125 Hz to 8 kHz in every window. (Without the sub-sample alignment step, L's 0.15–0.3 s window looks as if it rises with frequency; that was an alignment artefact *[exploratory]*.)
- **Caveat:** the library IRs are sine-sweep measurements with an unknown sweep (`NEXT_STEPS.md` §9). A time-variant unit's sweep IR depends on the sweep, so depth and rate can't be read off these numbers directly. The existence of the variation and its program dependence don't depend on the sweep.
- **The tail predictor sees the same thing within one file** (§3.2): a flat, signal-proportional unpredictable part at about −18 dB (L) and −18 to −22 dB (M).
- **Reading (H):** the "stochastic control" of the 1986 review (`algorithm_research.md` §7.2b) is real, **strong in M and L and not seen in S** (S's pair difference is timing-type and stays at −46 to −50 dB in the 125–500 Hz octaves at 0.6–1 s, so any flat gain-type variation in S is below about −45 dB there). Its flat spectrum points to randomised gains or switched tap positions (jumps of many samples), not small delay modulation, and it runs free (not restarted per capture). The review's "Grad der Stochastik differs per program" fits, though the review also calls S the densest. Rate, depth and which taps vary are open; the library can bound them (§6).

## 5. What this means

- **Q3 is answered for structure:** M and L have their own cluster scale (×2, ×4), onset scale and a different late recursion; the harness's unknown loop sets are assigned (§3.1).
- **The null target needs a decision for M and L.** Two captures of the same M/L setting from the same unit null each other only to −9 to −14 dB at 0.3–0.6 s and −3 to −7 dB after 1 s. No deterministic model can beat that against one capture. For M and L the waveform target can only apply to the deterministic part (cluster and first recirculation, about the first 0.1–0.2 s); beyond it, the right target is statistical (decay, spectrum, decorrelation rate between captures). S is unaffected: its tail is LTI to −35…−47 dB here.
- **Budget:** the M/L recursion is about 6 accesses (one write, five reads) on a 2N-word line; S's four loops are 8. The randomisation must also fit the budget.

## 6. Next tests

1. **Decorrelation rate vs settings:** does the pair-null decay rate depend on RT, HF or LF? Does it scale per recirculation (per N samples) or per second? This separates "random change per pass" from a free-running clock.
2. **Which part varies:** compare the pair difference with the five-tap recursion: is it consistent with randomised tap gains (β, γ), switched tap positions, or a varying read-out?
3. **Loop gains and filters from the recursion** at RT/HF/LF steps (Q4–Q7 for M/L), using the deterministic part.
4. **Harness:** per-program cluster templates (§1) with a 60 ms window, and the M/L lag sets as known loop sets.
5. **Non-Lin:** its late part (`1108-1254-1386` and `231-628-…` in the harness) through the same predictor.
