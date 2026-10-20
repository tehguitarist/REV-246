# Analogue chain: input S&H and the sub-20 Hz high-pass

*October 2026. Closes the two open analogue items of `NEXT_STEPS.md` §4. Tools: `tools/analog/sample_hold.py` (new), `tools/analog/lf_highpass.py` (new), `tools/analog/compare_kernel.py` (new options; the default is unchanged). Results: `analysis/analog/summary*.json`, `analysis/analog/lf_highpass.json`.*

Labels: **C** confirmed (primary source or direct measurement), **S** strong, **H** hypothesis.

---

## 1. Answers

1. **Input integrate/hold window: read, modelled, and negligible in band.** It is a leaky integrator with a critically damped 1 µs response, and the manual's timing plot gives it a 9.1 µs window. It settles to about 10⁻³ inside the window, so its discrete pole is p ≈ 0.0006. Its in-band effect is −0.04 dB at 10 kHz and −0.05 dB at 12 kHz, plus a 1.5 µs delay. It changes the kernel fits by ≤ 0.02 dB and **does not explain the skirt misfit or the skirt shift**. The fitted skirt shift is per unit with opposite signs: **+0.25 % (CLA), −0.44 % (Tomkin)**, so it is filter component tolerance (**H**), not a missing element.
2. **The sub-20 Hz high-pass is two things:**
   - the unit's own coupling, which includes a **3.35 Hz corner on the input board (C106 10 µF / R117 4.75 k) that had not been transcribed**: schematic sum 4.46 Hz;
   - a capture-chain part that differs between the two captures.

   Fitted on the cluster region of the IRs, the corners sum to:
   - **Tomkin: 6.25 Hz** (IQR 6.24–6.28 over 24 separate captures), i.e. the schematic + **1.8 Hz**;
   - **CLA: 17.75 Hz**, i.e. the schematic + **2 × 6.65 Hz** (second order).

   The CLA excess is 4× larger than any tolerance on the unit's corners, so it is **not the unit's design (S)**; most likely the CLA capture chain (**H**). S&H droop is ruled out by the circuit. A digital DC blocker, if present, is below about 1.8 Hz. The `model_s.py` form (15 + 15 Hz, sum 30 Hz) is wrong for both captures.

Side findings:
- the output S&H's time constant is (R203 + 2·P202)·C205, not (R202 + P202)·C205 (§2.4);
- the output staircase's zero-order hold is not in the current analogue model, and adding it makes the fit worse; this is open (§2.5).

---

## 2. Item 1: input integrate/hold

### 2.1 Circuit (input board 7 246 121, sheet 2/3, printed p. 41; text p. 37)

| Part | Function |
|---|---|
| Z102 (LF353, ½) | buffer after the input filter Z103 |
| R114 787 Ω | input resistor into node A |
| Z105 TL191 (quad switch), control AD SH A | **integrate:** 9-8 and 16-1 closed (node A → Z104 inverting input); **hold:** 3-4 and 6-5 closed (node A → ground), integrator input open |
| Z104 LF356, C105 1 nF in feedback | the integrator |
| R116 1 k → node B (C104 1 nF to ground) → R115 1 k → node A | feedback T from Z104's output to node A, so the integrator is **leaky**: it tracks −2·R115/R114 · Vin = −2.54 · Vin |
| R129 1.2 M from R128/R127 (−15 V × 8.2 Ω / 4.23 kΩ) | constant −24 nA into the integrator node: a DC offset, not signal-dependent |
| **C106 10 µF (tantalum) → R117 4.75 k → Z112 (LF356 follower)** | **AC coupling into the A/D summing point: 3.35 Hz first-order high-pass** |

Sheet 3/3 (p. 42): SH goes DC-coupled into the PCM53 summing node (R118 23.7 k) and the MAT-01/LT1016 comparator; P104/C109 set only the comparator offset. So there is no further coupling after C106.

The input amplifier (sheet 2/3) has C101/C102 680 nF into R105/R106 1 MΩ: 0.23 Hz.

**Confidence:** C for the parts. S for the topology reading: the switch pin functions agree with the text on p. 37.

### 2.2 Timing and response

- **Window (C, ±1 step):** the logic-analyser plot on p. 39 shows ADSHA high for **362 of 1,504 px of one period = 9.1 µs ≈ 44 timing steps** (208.3 ns each). The SAR clocks (ADCLK) run only while ADSHA is low, and the manual puts the conversion in the hold phase, so **ADSHA high = integrate**. The hold is 28.8 µs.
- **Dynamics (S):** with node A at virtual ground, the integrator's characteristic equation is s·C105·(R115 + R116)·(1 + s·C104·(R115‖R116)) + 1 = 0, which is (1 + s·1 µs)² for these values. It is **critically damped with τ = 1 µs** (R114 sets only the gain). A 60 Ω switch on-resistance moves this by about 7 %.
- **Model:** `sample_hold.input_sh(f, Ti, Ron)` solves the window exactly each period. It assumes VB has settled to Vo/2 during the 28.8 µs hold (R115‖R116·C104 = 0.5 µs) and returns the response from continuous input to the held sequence:

  | Window | Discrete pole | 5 kHz | 10 kHz | 12 kHz | Delay |
  |---|---|---|---|---|---|
  | **9.1 µs, R_on 0** | 6.2·10⁻⁴ | −0.010 dB | −0.038 dB | −0.052 dB | 1.5 µs |
  | 9.1 µs, R_on 60 Ω | 2.3·10⁻³ | −0.021 | −0.071 | −0.092 | 1.75 µs |
  | (4 µs, for scale) | 0.055 | −0.31 | −0.84 | −0.95 | 2.8 µs |

So at the real window, the input S&H is a near-transparent sampler: a slight aperture low-pass plus a fixed delay. **It is a discrete-time effect only in principle** (pole ≈ 0).

### 2.3 Kernel comparison, with and without (ideal op-amps)

`compare_kernel.py` options:
- `--input-sh` adds the stage at the schematic window;
- `--fit-ti` frees the window;
- `--skirt-shift` frees a common frequency scale on both filters;
- `--fit-hi 12300` fits over 0.5–12.3 kHz, skirt included (default 0.5–11 kHz).

rms dB error (kernel − model). Per-band values are on a 100 Hz grid, the skirt on 25 Hz. **Parameters** = fitted analogue parameters (each fit also has the 1 kHz normalisation).

**CLA kernel** (`ref/kernel_eugtone_96k.npz`):

| Model | Par. | p | 0.5–11 k | 0.5–1 | 1–2 | 2–4 | 4–8 | 8–11 | **11–12.3** | Other |
|---|---|---|---|---|---|---|---|---|---|---|
| **default** (filters × output pole) | 1 | 0.325 | **0.110** | 0.06 | 0.05 | 0.10 | 0.05 | 0.18 | 1.16 | |
| + input S&H (9.1 µs) | 1 | 0.324 | 0.109 | 0.06 | 0.04 | 0.10 | 0.05 | 0.18 | 1.17 | |
| + input S&H, R_on 60 Ω | 1 | 0.322 | 0.108 | 0.06 | 0.04 | 0.09 | 0.06 | 0.18 | 1.18 | |
| + input S&H, window free | 2 | 0.317 | 0.107 | 0.06 | 0.04 | 0.08 | 0.07 | 0.17 | 1.20 | Ti = 5.9 µs |
| default, fitted to 12.3 kHz | 1 | 0.297 | 0.359 | 0.05 | 0.01 | 0.08 | 0.39 | 0.51 | 0.82 | |
| + skirt shift, fitted to 12.3 kHz | 2 | 0.305 | 0.264 | 0.05 | 0.02 | 0.04 | 0.28 | 0.38 | 0.75 | scale **1.0025** |
| + input S&H + skirt shift | 2 | 0.303 | 0.269 | 0.05 | 0.02 | 0.05 | 0.29 | 0.38 | 0.76 | scale 1.0026 |

**Tomkin kernel** (`analysis/library/tomkin96/kernels/kernel_96000.npz`, 96 kHz set):

| Model | Par. | p | 0.5–11 k | 0.5–1 | 1–2 | 2–4 | 4–8 | 8–11 | **11–12.3** | Other |
|---|---|---|---|---|---|---|---|---|---|---|
| **default** | 1 | 0.301 | **0.109** | 0.03 | 0.04 | 0.07 | 0.07 | 0.17 | 1.10 | |
| + input S&H (9.1 µs) | 1 | 0.300 | 0.109 | 0.03 | 0.04 | 0.07 | 0.08 | 0.17 | 1.10 | |
| + input S&H, window free | 2 | 0.300 | 0.109 | 0.03 | 0.04 | 0.07 | 0.08 | 0.17 | 1.10 | Ti = 23.5 µs (no preference) |
| default, fitted to 12.3 kHz | 1 | 0.294 | 0.144 | 0.03 | 0.03 | 0.05 | 0.16 | 0.19 | 1.06 | |
| + skirt shift, fitted to 12.3 kHz | 2 | 0.277 | 0.318 | 0.02 | 0.01 | 0.09 | 0.38 | 0.41 | 0.84 | scale **0.9956** |
| + input S&H + skirt shift | 2 | 0.275 | 0.324 | 0.02 | 0.01 | 0.10 | 0.39 | 0.41 | 0.84 | scale 0.9956 |

What this shows:
- **The input S&H changes nothing that matters** (≤ 0.02 dB in any band). A free window either runs long (Tomkin 23.5 µs, equivalent to "fully settled") or takes a 6 µs value worth 0.003 dB (CLA). The kernels cannot see the window; the schematic and timing plot fix it. **S.**
- **The skirt (11–12.3 kHz) is the weak spot in both units, at 1.1–1.2 dB rms**, with the passband at 0.11 dB.
  - A common frequency scale on both filters helps the skirt (to 0.75 / 0.84 dB) but costs 0.15–0.2 dB in the 4–11 kHz bands. The best scales are **+0.25 % (CLA) and −0.44 % (Tomkin)**: the "about 0.2 %" in `NEXT_STEPS.md` was one unit's value, and the two units shift in opposite directions.
  - So the skirt misfit is a **per-unit shape difference** (each section's 220 pF / 68.1 k tolerance moves its own pole or zero; a single scale cannot follow that) (**H**). The input S&H plays no part (**S**).
  - Caution: the kernels are only 5.2 ms long (499 taps at 96 kHz, ≈ 190 Hz resolution, 1.5 ms Hann taper; `tools/harness/library.py:build_kernel`). The skirt (≈ 60 dB/kHz) is partly smeared by the truncation. A per-section tolerance fit is only worth doing against a longer kernel or in the time domain.

### 2.4 Output S&H time constant (correction to `EMT246_unit_facts.md` §2b)

On p. 48, R203 (1.21 k) runs from Z209's output back to the node **between R202 and P202**, not to the integrator input. While integrating, node X sits at (Vd + Vo)·P202/(2·P202 + R), so the charge current is (Vd + Vo)/(R + 2·P202):
- **τ = (R203 + 2·P202)·C205**, 1.2–11.2 µs over the trimmer's range;
- gain is −1 (R202 = R203), independent of P202, as a "charge current" trim should be.

The window length for the fitted p is then Tw = τ·ln(1/p): about 7.5 µs at mid-trim for p = 0.3. The DA SH pulse width is not given in the manual (no output timing plot), so **H**.

### 2.5 Open: the output staircase's zero-order hold

Between its windows each output S&H holds its value for the rest of the period (the one DAC serves R and L in turn). So the output is a staircase, and its spectrum carries the zero-order-hold sinc, sinc(f/26.37 kHz): −0.5 dB at 5 kHz, −2.2 dB at 10 kHz. This is certain from the circuit (**C**). The default model has no ZOH; its p carries all the roll-off. `compare_kernel.py --zoh` adds it:

| Kernel | Model | p | 0.5–11 kHz rms | 8–11 k | 11–12.3 k |
|---|---|---|---|---|---|
| CLA | ZOH + pole, ideal op-amps | 0.233 | 0.349 | 0.43 | 2.40 |
| CLA | ZOH + pole, TL07x 3 MHz | 0.272 | 0.184 | 0.28 | **0.66** |
| Tomkin | ZOH + pole, ideal | 0.206 | 0.366 | 0.45 | 1.91 |
| Tomkin | ZOH + pole, 2.5 MHz | 0.255 | 0.191 | 0.32 | 2.27 |

- A GBW scan from 1.5 to 4 MHz has its best passband at 2.5–3 MHz (0.18–0.19 dB), still worse than the default's 0.11 dB.
- An exponential charge pulse of free length on top adds nothing: it fits to its lower bound.
- So the physically complete model (staircase + charge pole + real op-amps) fits the kernels worse in the 4–11 kHz bands than the default does. The likely reading (**H**): the transcribed input-filter pre-emphasis, or the op-amp model, is off by about 0.2–0.4 dB at 4–11 kHz, and the default's p is an *effective* parameter that absorbs the ZOH along with that error. Then p = 0.325 / 0.301 is not the physical charge pole. The physical one would be about 0.25–0.27 with the ZOH.
- **The default stays as is** (best fit, one parameter). Settle this before using p to infer trimmer positions or window lengths. It does not affect using the default chain for nulls.

---

## 3. Item 2: the sub-20 Hz high-pass

### 3.1 Coupling corners in the transcribed schematics (C)

| Where | Parts | Corner |
|---|---|---|
| **Input board, S&H → A/D summing point** | **C106 10 µF tantalum, R117 4.75 k** | **3.35 Hz** (newly transcribed) |
| Input filter Z103 | C3109 1 µF, R3134 221 k | 0.72 Hz |
| Input amplifier | C101/C102 680 nF, R105/R106 1 M (+ R101‖R102) | 0.23 Hz |
| Output filter Z210 | C3009 1 µF, R3031 1 M | 0.16 Hz |
| Output drivers | 100 µF into the load | ≪ 1 Hz |
| **Sum** | | **4.46 Hz** |

### 3.2 Method (`tools/analog/lf_highpass.py`)

Before the reverb starts (R: 570 device samples after tap 0; L: 240), an S IR is just the cluster taps through the analogue and capture chains. The cluster's positions and signs are known and its gains are unity.
- Both the IR and the model are low-passed at 1 kHz (4th-order Butterworth, zero phase); below that the band-limiting chain is flat.
- The free parameters are then the high-pass corners plus gain, delay and DC offset (3 nuisance parameters).
- A high-pass leaves a slow negative tail after each tap. Its start level is −Σωᵢ and its curvature sets the order.
- Fit on R (tap 0 − 1 ms … + 20 ms); **held out on L** (corners fixed from R; gain/delay/DC refitted; to + 8 ms).
- Errors are given as residual energy re the window's signal energy.
- Data: the CLA IR (one capture), and **24 Tomkin S captures** at 96 kHz (highest cluster-to-reverb, AMPL 90, all separate captures; `files.csv`).

Candidates:
- **schematic**: the 4.46 Hz set, fixed;
- **sch+1**: schematic + one free corner;
- **sch+2**: schematic + two free corners;
- **free2**: two free corners alone, the `model_s.py` form.

### 3.3 Results

| Capture | Candidate | Free par. | Extra corners (Hz) | Σ corners (Hz) | R in-sample | **L held out** |
|---|---|---|---|---|---|---|
| CLA | schematic | 0 (+3) | — | 4.46 | −25.4 dB | −23.8 dB |
| CLA | sch+1 | 1 (+3) | 12.86 | 17.32 | −38.6 | −36.4 |
| CLA | **sch+2** | 2 (+3) | **6.65, 6.65** | **17.75** | **−40.9** | **−40.2** |
| CLA | free2 | 2 (+3) | 8.72, 8.72 | 17.45 | −39.0 | −37.4 |
| CLA | fixed 15 + 15 (`model_s.py`) | 0 (+3) | — | 30 | −25.2 | −28.1 |
| Tomkin (median of 24) | schematic | 0 (+3) | — | 4.46 | −39.9 | −38.1 |
| Tomkin | **sch+1** | 1 (+3) | **1.80** | **6.25** (IQR 6.24–6.28) | **−43.7** | **−44.0** |
| Tomkin | sch+2 | 2 (+3) | 1.37, 0.42 | 6.24 | −43.7 | −44.0 |
| Tomkin | free2 | 2 (+3) | 3.15, 3.15 | 6.29 | −43.6 | −43.9 |
| Tomkin (one file) | fixed 15 + 15 | 0 (+3) | — | 30 | −20.3 | −21.8 |
| Tomkin (one file) | CLA's high-pass (sch + 6.65 + 6.65) | 0 (+3) | — | 17.75 | −26.9 | −26.3 |

Noise and robustness:
- **CLA:** the LF-band floor from the 10 ms pre-delay is −44.5 dB, so sch+2 is +3.6 dB over noise held out.
- **Tomkin:** the files have no usable pre-roll, so the floor isn't measured. Their spread across 24 independent captures is the check: Σ corners IQR 0.04 Hz.
- **Robustness to the analysis low-pass:** at 300 Hz / 1 kHz / 3 kHz, CLA Σ = 16.4–17.8 Hz and Tomkin Σ = 5.5–6.4 Hz.
- The kernels' own tails (1.5–3.6 ms) give the same ordering but read high (28 / 12.6 Hz), because the short window mixes in the chain's other slow parts.

### 3.4 Verdict on the candidates

| Candidate | Verdict | Evidence | Conf. |
|---|---|---|---|
| **Input-board coupling not yet transcribed** | **Found:** C106/R117, 3.35 Hz, the largest unit corner. With the rest of the schematic (4.46 Hz in total) it accounts for **about 70 % of Tomkin's high-pass** and a quarter of CLA's | §2.1, §3.3 | C (parts) / S (share) |
| **S&H droop** | **Ruled out.** In hold, the input integrator's input is open (only R129's constant 24 nA flows in, a DC offset); the output integrator's switch is open. Neither has a leakage path proportional to the signal, so neither can act as a high-pass | p. 41, p. 48 | S |
| **Digital DC blocker in the DSP** | **Bounded at about 1.8 Hz or below.** It would be the same in both units (same firmware), and Tomkin's total leaves only 1.8 Hz beyond the schematic | §3.3 | S |
| **Capture chain** | **The CLA excess (≈ 13 Hz, second order) is the capture chain** most likely: for example the interface's D/A and A/D DC blocking (two first-order stages fit better than one by 2.3 dB in-sample and 3.8 dB held out), or the sweep's low edge in the deconvolution. Tomkin's chain adds ≤ 1.8 Hz. The two units share the design; tolerance on the unit's corners (10 µF tantalum ±20 % → 3.35 ± 0.7 Hz) can't produce a 13 Hz difference. Not excluded: a modified CLA unit or different board revision | §3.3 | S (not the unit's design) / H (capture chain) |

### 3.5 What to use

- **The unit's high-pass, fixed, no free parameters:** first-order corners 3.35, 0.72, 0.23, 0.16 Hz (`lf_highpass.SCHEMATIC_HZ`).
- **Per capture chain, counted as parameters:**
  - Tomkin: + one corner at 1.8 Hz (1 parameter, shared by all Tomkin files);
  - CLA: + two corners at 6.65 Hz (1–2 parameters).
- **Not 15 + 15 Hz.** In the 63 Hz octave, the complex error of 15 + 15 against the measured high-pass is −14 dB (CLA) and −8.6 dB (Tomkin); in the 31.5 Hz octave it is −9 / −4 dB. A long free FIR absorbs part of that, which is why `model_s.py` still reached −37/−41 dB at 63 Hz. But a compact model without the FIR cannot, and **the measured form should lift the 63 Hz band** (not yet tested in the reverb windows; `model_s.py` is not edited here).
- For the plugin (`NEXT_STEPS.md` §11): only the unit's 4.46 Hz set belongs in the emulation; the capture-chain corners do not.

---

## 4. Reproduce

```bash
.venv/bin/python tools/analog/compare_kernel.py
.venv/bin/python tools/analog/compare_kernel.py --input-sh --tag insh
.venv/bin/python tools/analog/compare_kernel.py --skirt-shift --fit-hi 12300 --tag skirt
.venv/bin/python tools/analog/compare_kernel.py --zoh --tag zoh
.venv/bin/python tools/analog/compare_kernel.py --kernel analysis/library/tomkin96/kernels/kernel_96000.npz --tag tomkin96000_insh --input-sh
.venv/bin/python tools/analog/lf_highpass.py --n-tomkin 24
```

The other variants in `analysis/analog/summary_*.json` use the tags in their names (`fitti`, `insh_ron60`, `skirtfit`, `insh_skirt`, `insh_skirtfit`, `fitti_skirt`, `zoh_insh`, with and without the `tomkin96000_` prefix).

## 5. Limits

- All fits are magnitude-only against 5.2 ms kernels (the existing `compare_kernel.py` method). Phase, and the skirt finer than about 190 Hz, are not tested.
- The p. 39 timing plot is a logic-analyser screen at about 25 ns/px; the window is ±1 step (±0.2 µs), which changes the input S&H response by < 0.01 dB.
- The high-pass fit uses the cluster region only (≤ 20 ms). It measures Σ corners well and the order moderately (CLA), and cannot split corners below about 2 Hz from each other (Tomkin). The 100 µF output coupling depends on the capture interface's input impedance, so it is not separable from the capture chain.
