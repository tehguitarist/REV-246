# Period DSP practice and the program-S injection: hypotheses for the Q2 bench

*Research note, October 2026. §0–§6 are reading and reasoning: they collect how 1978–1990 reverbs and the TMS32010 were actually used, and turn that into candidate structures for the dense per-loop injection/read-out that Q2 found (`tools/q2_injection.py`, `branches` stage), each with a test recipe. §7 gives the first offline results (`tools/q2_candidates.py`), in-sample only and pending the held-out bench.*

> **Update (October 2026, §7, in-sample):** K1 holds. The frozen flat injection is **five equal-magnitude signed taps per loop**. R uses the same taps shifted by 26 or 332, and the loop paths carry a factor that cancels the bench's comb460. One tied model explains **98.7 %** of the in-band energy of all eight excitations, with about 60 discrete and 2 continuous parameters (against 27.9 k free coefficients). K2–K7 are not needed for flat settings or are not supported. Not yet confirmed on the held-out bench.

**Labels** (as in `EMT246_unit_facts.md`):
- **C** = confirmed by a primary source (document read) or a direct measurement in this repo;
- **S** = strong: several measurements or several independent sources agree;
- **H** = hypothesis.

Sources are numbered [Pn] and listed in §8 with URLs. Measured facts cite the repo doc and section.

---

## 0. Summary

- **First results (§7, in-sample):**
  - K1 is confirmed offline: 5 sign-only taps of one magnitude per loop. Each R tap is an L tap delayed by 26 or 332. The loop feeds bypass comb460, and the first-pass FIR is almost empty.
  - In the hardware's terms that is about 20 reads, two shared output delay lines and no multiplies: about 24 accesses, inside the budget.
  - The bench test (held-out nulls with the tied model) is next.
- **What has to be generated:** one fixed excitation per loop and output channel, about Pᵢ + 1000 samples long, starting 230 samples after cluster tap 0. It is the same at every RT step (flat filters), carries no Gᵢ factor, differs for every loop, and relates L to R by the 26/332 lags (§1). The hardware has about 20–35 memory accesses and about 130–170 instruction cycles left for it (§2).
- **Two cheap ways a fixed overhang past Pᵢ arises (§1.2):**
  - **Read taps past the feedback point.** In a base + offset memory a loop line can be longer than its loop. Feedback is read at Pᵢ, and output taps read older samples of the same stream out to Pᵢ + 1000. Each such tap is a fixed numerator term, does not depend on RT, and costs one read. This is the EMT 250 patent's tapped long loop and the BBC rule "feedback from the last tap only, density from extra taps" [P3, P10]. **H** for the 246.
  - **A recursive pre-stage per loop** (an all-pass or comb with a fixed gain). In the bench's coordinates its ringing appears as a long dense excitation. AP230 and comb460 ring to −13 and −15 dB after 1000 samples, the same order as the overhang. **H**
- **Ranked candidates** (§4; plausibility first, bench cost second):
  1. **K1, tapped read-out along an extended loop line**, optionally with multi-point writes (EMT 250, BBC, Stautner–Puckette, Dattorro). Cost 1 access per tap. **Test:** in-band sparsity of the saved long-window excitations. This is offline, takes minutes, and is decisive on budget grounds: about ≤ 4–8 taps per loop and channel can be afforded.
  2. **K2, a recursive pre-stage per loop**, or loops fed from different points of the diffuser chain (Schroeder 1962, the EMT 250's short all-pass loops, the EMT 245's series/parallel loops). Cost 2 accesses per stage. **Test:** deconvolve by a dictionary of all-passes and combs, then the K1 sparsity score. Offline.
  3. **K3, a tapped injection line per loop, convolved with the tapped read-out** (Matsushita 1982, Moorer 1979, Kendall–Martens 1984). Support = injection span + read span. **Test:** a shared convolutive factor between L and R of a loop.
  4. **K4, read-out taps on both sides of a loop filter placed partway along the line** (Dattorro/Lexicon tank). Explains why one fixed excitation fails once HF/LF move (the `all` subset sits at −29 dB). It does not explain the overhang. **Test:** one `branches` run per split with a filter-only branch.
  5. **K5, a nested all-pass per loop** (Gerzon 1972). **K6, a per-channel decorrelator** (Schroeder–Logan). **K7, HF/LF tap tables.** Lower; tests given.
- **Sources** (§3, §8): TI's TMS32010 documents; period patents (EMT, Quad Eight, Yamaha, Matsushita, Kendall–Martens, and German DE/AT filings from DEPATISnet); the literature from Schroeder to Dattorro; EMT's own brochures and manuals; Valhalla DSP. Together they give EMT-lineage precedents for every top candidate, but **nothing describes the 246's algorithm itself**: EMT never patented it.
- **What this rules out** (§5):
  - chained or series long loops and loop-sum feedback (RT-independence);
  - a free dense FIR (budget);
  - one tapped line shared by all loops (already negative);
  - noise-like ±1 tap sequences;
  - on-chip FIRs;
  - a *dense* K1 (if the excitations stay dense in every basis, tapped read-out alone is over budget).

---

## 1. What has to be explained

### 1.1 Measured facts the candidates must meet

| # | Fact | Conf. | Source |
|---|---|---|---|
| F1 | With flat filters (LF 1.0, HF 1.0), **one excitation per loop and channel, the same for every RT**, predicts held-out RT steps (0.8 and 1.8 s, 8 captures) at −52.7 / −49.8 / −45.3 dB (R, 50–150 / 150–300 / 300–600 ms, fade compensated), at the repeat-capture ceiling. 27,924 free parameters, a single branch (α = 0) | S (narrow held-out set: 2 settings) | `analysis/q2/branches_flat.log` (`branches_flat_0_ext1700`) |
| F2 | **Support ≈ Pᵢ + 1000 samples**, starting 230 samples after cluster tap 0: ext 1000, 1300, 1700, 2000 and 2600 all reach −50…−51 dB at 50–300 ms; ext 300 needs a branch through Gᵢ and still plateaus at −35 dB; ext 3300 is 2–3 dB worse (over-parameterised). So the content ends somewhere between Pᵢ + 300 and Pᵢ + 1000, or is too weak to matter after that | S | same log |
| F3 | **The overhang beyond Pᵢ is a fixed numerator term:** with the long window the G branch adds nothing (`0_1_ext1700` ≈ `0_ext1700`). With the short window the best branch set was α = −1, 0, 1 (−43.4 / −42.9 dB held out at 50–150 / 150–300 ms): a fixed tap at Pᵢ + d, written with taps shorter than Pᵢ, needs a Gᵢ⁻¹ term, so the short window was approximating exactly such content. Either a read past the feedback point of loop i's own line (§1.2) or a feed from outside the loop gives it | S (F1–F2, the α = −1 run) / H (the two readings) | `branches_flat.log` (`-1_0_1`, `0_1_ext1700`, `0_ext1700`) |
| F4 | **With HF/LF moving, the same single-branch form fails:** held out −35.0 / −29.2 / −27.9 dB on the `all` subset (128 captures). Part is the loop-filter laws at LF 0.1/0.4 and HF 0.2/0.4 (the per-setting free injection is +5…+18 dB there too, `model_s_results.md` §15.3), part may be filter placement (K4) | S (fact) / H (split) | `analysis/q2/branches_all.log`; model §15.3 |
| F5 | Loops share no fine structure (no common tap set, no shared long input); L and R of one loop are related: lags **26 and 332** carry about half of R's periodic part with **per-loop weights**; R's first arrivals are L's delayed by exactly 332 samples | S | model §5, §10.3, §12.1, §14.6 |
| F6 | Common chain comb460 (g 0.445) · AP20 · AP72 · AP230 (g −0.70), delays sharp to ±1 sample; AP230 is on the input side and carries a damping step b from a 4-value table on RT × HF | S | model §11, §14.3 |
| F7 | Per-loop mode magnitudes are Rayleigh-like (spread about 5 dB), i.e. not a pure all-pass read by one tap | S | algorithm_research §7.2a (T2) |
| F8 | M and L reproduce S's cluster and onsets at ×2 / ×4 (integer, ±0.001 sample): the early section is a tap table scaled per program | S | ml_structure §1–2 |
| F9 | Budget left for injection and read-out: about **20–35 memory accesses** per sample (the TOM holds 256 entries, so time is the limit, not the table) | C (TOM, cycles) / H (split of the budget) | unit facts §2b, §6 |

### 1.2 What the bench's coordinates hide

The bench models each channel as

> y = cluster + B · post · (F + Σᵢ Xᵢ / (1 − Gᵢ z^−Pᵢ)),  post = comb460 · AP20 · AP72 · AP230

(`q1_generalise.py`, `Joint`, `group_spectra`). Each Xᵢ is therefore **the loop's feed and read-out, divided by the full common chain**, with all of the loop's poles in the 1/(1 − Gᵢ z^−Pᵢ) factor.

**(a) Reads past the feedback point.** In loop i's line the stored stream is wᵢ = uᵢ + Gᵢ · wᵢ(n − Pᵢ). A read at any lag τ gives z^−τ · Uᵢ / (1 − Gᵢ z^−Pᵢ), whether τ is below or above Pᵢ, as long as no other element has overwritten those words. In a base + offset memory, that only means leaving a gap of about 1000 words after the loop's feedback tap: 4 × 1000 words, against about 39 K free words. So read taps at Pᵢ + d are fixed numerator terms in xᵢ. They are RT-independent and cost one access each (§2.3). The same holds on the write side: several writes of the input into the loop line at different offsets (Stautner–Puckette [P14]) add terms at lags below Pᵢ. **H** (for the 246) / **C** (the arithmetic)

**(b) Recursive stages outside the loop.** Suppose loop i's true feed is Qᵢ · tᵢ, with tᵢ a few integer taps and Qᵢ the recursive stages its signal passed. Then

> xᵢ = tᵢ · Qᵢ / post.

- If Qᵢ = post · Eᵢ (one extra stage of its own), xᵢ = tᵢ · Eᵢ: **causal, dense, and as long as Eᵢ rings**. An all-pass or comb with a fixed gain is RT-independent and does not pass the loop gain, so it meets F1–F3 (candidate K2). **H**
- If loop i is fed *before* one of the common stages (Qᵢ = post / Eᵢ), xᵢ = tᵢ / Eᵢ. For a comb that is a short FIR (tᵢ · (1 − g z^−M)); for an all-pass it is the time-reversed all-pass, i.e. anticausal, which a causal window starting at 230 can only approximate. So feeding *earlier* in the chain does not explain a long causal overhang; feeding through *more* stages does. **H**
- If a loop is fed from *inside* a stage's delay memory (a tap on the AP230 line rather than its output), Qᵢ has that stage's comb denominator without its all-pass numerator: xᵢ = tᵢ / (−g + z^−M), again anticausal. **H**

Ringing of the stages already identified, for scale (arithmetic, C):

| Stage | Decay | After 1000 samples | After 1700 samples |
|---|---|---|---|
| comb460, g 0.445 | 1.53 dB per 100 samples (403 dB/s) | −15.3 dB | −26.0 dB |
| AP230, g 0.70 | 1.35 dB per 100 samples (355 dB/s) | −13.5 dB | −22.9 dB |
| AP72, g 0.70 | 4.3 dB per 100 samples | −43 dB | −73 dB |
| AP20, g 0.70 | 15.5 dB per 100 samples | −155 dB | – |

A per-loop stage with a delay of a few hundred samples and a gain of 0.4–0.7 would leave an overhang of about 1000 samples at the −15…−25 dB level, which fits F2. **H**

**So (a) and (b) give the same long, fixed xᵢ on the bench and differ only in what xᵢ is sparse in:** integer taps for (a), taps after a recursive stage for (b). That is the first thing to test (§6).

---

## 2. TMS32010 and 246 hardware constraints that matter

### 2.1 The processor (TI TMS32010, 1983)

All **C** from the TI *TMS32010 User's Guide* (1985) [P1] unless marked; page numbers are positions in the downloaded PDF, check against the scan before quoting.

| Fact | What it means for the 246 |
|---|---|
| Instruction cycle = CLKIN / 4 (p. 43); 19.2 MHz → 208.3 ns → **182 instructions per sample** per DSP | 364 per sample for the pair (arithmetic) |
| **1 cycle:** LT, LTA, LTD, MPY, MPYK, PAC, APAC, SPAC, DMOV, ZALH, ADD/SUB (with shift), SACL, SACH, LARP, LARK, MAR, ZAC. **2 cycles:** IN, OUT, all branches (B, BANZ, BIOZ, CALL). **3 cycles:** TBLR, TBLW (Table 3-2, pp. 54–56) | A memory access costs an IN or OUT (2 cycles) on top of the arithmetic |
| **Multiply-accumulate is a pair:** LTA (or LTD) then MPY/MPYK, "pipelined multiply and accumulate operations at 400-ns rates", i.e. 2 cycles per product (p. 28, p. 90) | A weighted product costs 2 cycles even with the data on chip |
| **LTD = LT + APAC + DMOV** in one cycle: the FIR idiom keeps the delay line in on-chip RAM and shifts it as it multiplies, oldest sample first (§5.3.2, pp. 156–157; Lovrich & Simar, TI appnote SPRA003 = [P2] ch. 3, with straight-line and BANZ-looped FIRs; straight-line is faster) | An on-chip FIR costs 2 cycles per tap and is capped by the RAM: about 85–90 taps per chip per sample at most, with nothing else done |
| **No repeat instruction** (RPT/RPTK are TMS32020 only); loops cost a 2-cycle BANZ per pass. Only **two auxiliary registers** (AR0, AR1), 9-bit auto-increment/decrement wrap (pp. 20, 30–32). **LTP does not exist** on the 32010 | Production code at this rate is straight-line. 182 instructions per DSP per sample is the whole program for one sample |
| **On-chip data RAM 144 × 16**; direct addressing in two pages (0–127, 128–143); DMOV works on 0–127 only (pp. 30–31, 77) | No delay line beyond a few dozen words fits on chip; all reverb delays must be external |
| 16 × 16 multiplier, 32-bit P and accumulator; input barrel shift 0–15 (left), output shift on **SACH only, by 0, 1 or 4** (pp. 20, 28, 30, 102) | Power-of-two gains are free; no right shift on store |
| **MPYK: 13-bit signed immediate, Q12** (step 2⁻¹² ≈ 2.4 × 10⁻⁴), stored with SACH …,4 for full precision (pp. 93, 144). MPY from data RAM gives Q15 | A fixed gain (diffuser, per-loop pre-stage) is most cheaply an MPYK constant on a 1/4096 grid; a control-dependent gain (loop G, filters) must sit in RAM. **Test (H):** whether fitted fixed gains snap to 1/4096 better than to 1/32768. Not resolvable yet: g_AP is known only to ±0.02 |
| **External data memory only through I/O ports** (8 in, 8 out, PA0–PA2) or TBLR/TBLW to program memory. TI's own suggestion for large memory is an **external address counter** on a port: OUT the start address, then sequential IN/OUT auto-increment (§6.1.3, p. 163; J. Bradley, "Interfacing External Memory to the TMS32010", [P2] ch. 11, p. 319: random access 4 cycles, sequential 2) | Random-access taps on a TI-style interface cost about 6 cycles (OUT address + IN + LT + MPY). **The 246's base + TOM hardware removes the address cost entirely:** every access is "sequential" in the TOM, yet lands at an arbitrary offset. That is EMT's improvement on TI's counter scheme, and it is what makes many cheap taps possible (§2.2) |
| Coefficients: TBLR from program memory (3 cycles), MPYK immediates, or a host link on an I/O port with BIO/interrupt handshake (§5.4.1, §6.1.1, pp. 157–158, 162); IN + TBLW (≥ 5 cycles) can write downloaded words into program RAM | The Z80's parameters arrive over a port into data RAM. Writing them into program memory as MPYK immediates would be self-modifying code; the 246's program memory is PROM, so **control-dependent gains are MPY from RAM** (C for PROM, H for the conclusion) |
| **No period TI document describes reverb, comb or all-pass code** on the TMS32010/20: the 1986 applications book [P2] has no such chapter, and no 1983–88 AES/ICASSP paper on reverb on a TMS320 was found | The 246 team worked without a TI template; the structures come from reverb practice (§3), the coding idioms from the FIR appnotes |

### 2.2 How the 246 wraps it (service manual, `EMT246_unit_facts.md` §2b)

- **Every delay access goes through `base + TOM[n]`:** a 16-bit base counter advanced once per sample by the program (IN port 1, "IN BAZ"), plus a 256 × 16 tap-offset memory whose index counter steps on **every** memory access. Reads are `IN 0` (from the read register MOR), writes `OUT 0` (to the write register MIR). **C**
- **Consequences for a designer** (**H** unless marked):
  - The program never computes an address. A read costs the IN instruction and nothing else, so **a tap anywhere in the 64 K words costs the same**: on a loop line, inside an all-pass's line, on the pre-delay line, or at an offset that lands in a neighbouring element's memory. A tap on the diffuser chain is as cheap as a tap on a loop.
  - The **order** of accesses is fixed by the TOM sequence and the straight-line program; the offsets are data. S, M and L could share code and differ in TOM contents, which fits the exact ×2 / ×4 scaling of the cluster and onsets (F8). (M and L have a different late recursion, so at least part of the code differs.)
  - In one circular memory, "line length" is only the gap between one element's write offset and the next element's. Reading past it reads the neighbour. Designers of the period laid elements out contiguously (see §3 on Lexicon); a tap table with a few such cross-reads would look, to us, like per-loop content at lags longer than the loop. **H**
  - The TOM is written by processor A (OUT 2) while running, so per-setting offsets are possible; the parameters come from the Z80 over the external bus (OUT 3 / IN 4–5). **C** (ports) / **H** (use).
  - Processor B has no external I/O ports fitted; it reaches memory and A only. **C** (schematic) — how B's memory accesses are addressed (through the same TOM sequence or its own) is not read yet.
  - The DRAM (µPD41464-12, 120 ns access) behind the MOR/MIR registers is shared by both processors through an arbiter, so the total number of accesses per sample is also bounded by the memory cycle, not only by instructions. **C** (parts) / **H** (bound).

### 2.3 Cost model per element (cycles per sample; instruction counts **H**, timings from §2.1)

| Element | Accesses | Cycles (approx.) | Notes |
|---|---|---|---|
| Signed, unweighted read tap (cluster style) | 1 | 3 (IN 2 + ADD/SUB 1) | The 8 cluster taps are equal in magnitude with signs only (unit facts §5): consistent with no multiply |
| Weighted read tap | 1 | 4 (IN 2 + LTA 1 + MPY/MPYK 1) | MPYK if the gain is a constant in ROM; MPY from RAM if it comes from the Z80 |
| Write (store a result into a line) | 1 | 3 (SACH 1 + OUT 2) | |
| Comb or all-pass section with a fixed gain | 2 | ~10–14 | read, two products, two adds, write |
| Nested all-pass (two lines) | 4 | ~20–26 | |
| Loop gain + HF one-pole + LF shelf, state on chip | 2 (read, write) | ~15–20 | per loop |
| Power-of-two gain | 0 | 0–1 | shifts are free on load (0–15) and on SACH (0, 1, 4) |

**Budget arithmetic (H):** 2 × 182 = 364 cycles per sample. The identified S structure (about 34 accesses, §6 of the unit facts: cluster 8 × 3, loops 4 × ~18, diffuser and 460 ~45, read-outs 8 × 4, I/O and housekeeping ~20) comes to roughly 190–230 cycles. That leaves about **130–170 cycles**: room for 8–12 more all-pass or comb sections, or 30–40 more weighted taps, but nowhere near a dense FIR.

---

## 3. Period practice

### 3.1 Units and patents: how loops were fed and read

| Source | Structure (as documented) | Conf. | Relevance to the 246 |
|---|---|---|---|
| **EMT 250**, US 4,181,820 (Blesser & Bäder, priority 1977) [P3] | Pre-delay → about 15 loops (long tapped loops "preferably at least 300 ms" plus short loops "with all-pass characteristic") → output "discorrelator and distributor" (delays t17–t19, 4 outputs). Long loop: **feedback from the end of the line only** (through the LF/MF/HF decay element), **many read taps along the line** summed with exponentially falling *energy*; a parallel 3-tap delay adds early density; random gain/delay variation for eigentone density | C | Single feedback point plus read taps along the line is the 250's answer to dense read-out. Short all-pass loops next to long loops = per-loop recursive stages (K2). Long tapped loops are already negative for the 246 (§4 table) |
| **EMT 250 algorithm drawing** ("Blesser Associates, June 77"), visible about 2:32–2:36 in Universal Audio's 2009 EMT 250 video, as described on Gearspace by S. Costello and others [P31, P32] | **Three parallel combs with three read taps per comb**; the taps are combined (partly hidden by a graphic) and sent into 3 or 4 short delay lines. Costello's reading of those short lines: each sums **two adjacent output taps, one scaled by +f, the other by −f and also added unscaled**, i.e. out = f·x[τ] + (1 − f)·x[τ + 1], **linear interpolation done on the DSP with one coefficient**. A forum user who owns a 250 adds: 3 output taps for L and 3 for R (the first R tap a straight pass-through), rear outputs from the comb matrix only | S (forum reading of a blurry video still; not seen by us) | (1) **Several read taps per comb**, then a small output network: EMT's own version of K1 (and of K6 for the output lines). (2) **Adjacent-tap interpolation in EMT's lineage**: the 246's AP230 damping step F = (1 − b) + b z⁻¹ (b = 0 / 0.17 / 0.33 / 0.50) is exactly that two-tap form, so reading it as a fractional delay 230 + b selected by a table (`model_s_results.md` §14.3) gains a period precedent. M and L's common −0.80-sample lag offset (`ml_structure.md` §3.3) could be the same idiom. **H** |
| **EMT 245** service manual (1981–82) [P4] | "Der Programmablauf … entspricht einer gemischten Serien- und Parallelschaltung von insgesamt neunzehn solcher Schleifen" (19 loops of delay + multiplier + adder, mixed series/parallel); several microprocessors in parallel because one is not fast enough for two uncorrelated reverb signals; 15-bit A/D, 16-bit internal words, **rounding at every multiply** | C (read locally) | EMT's own immediate precedent: loops **in series** with other loops. In the 246 the four long loops cannot be in series (§5), but short loops in series *in front of* each long loop are exactly K2. The rounding remark is a lead for the fixed-point model (Q6) |
| **Quad Eight**, US 4,219,880 (Nichols, filed 1978) [P5] | About 8,200-word delay RAM at 24 kHz; an address counter steps once per program cycle and **a per-instruction offset from the operand ROM is added**, next to the multiply coefficients | C | The 246's base counter + TOM, documented six years earlier. A designer using it thinks in "taps at offsets", each costing one instruction |
| **Yamaha (Nippon Gakki)**, US 4,472,993 (priority 1981) [P6] | Read address = counter + ΔADR; **10 early-echo taps with gains summed**, then parallel combs, then series all-passes; tap delays and gains from coefficient and delay-length memories; 8 selectable reverb blocks, 16 programs time-shared | C | Base + offset again; a tapped early section **feeding** the combs (Moorer's scheme in hardware); per-program tap tables |
| **Matsushita**, US 4,584,701 (priority 1982) [P7] | One RAM with an address generator; 6 taps (9.9–83.2 ms) in two groups of 3, **each group summed and fed into 2 of the 4 recirculating delays** (83.5 / 74.5 / 63.3 / 58.9 ms, ratio 1 : 0.9 : 0.8 : 0.7) | C | **Per-loop tapped injection** into four parallel loops of about the 246's size (the 246's are 50–60 ms). Candidate K3 |
| **Yamaha** JPH0683515B2 (1985) [P8] | Tapped register for reflections plus a second long tapped register for the tail; per-hall tap and gain tables in ROM, loaded into RAM | C (machine translation) | Per-program tap tables (F8, K7) |
| **Lexicon / Griesinger**: Dattorro 1997, "Effect Design Part 1" (tank "in the style of Griesinger", from conversations with Blesser and Griesinger) [P9] | Input: pre-delay, one-pole band-limit, **4 series all-passes** (142, 107 at g 0.75; 379, 277 at 0.625) at 29,761 Hz. Tank: figure-eight of two legs, each = modulated all-pass (g −0.70) → delay → **one-pole damping** → decay gain → all-pass (g 0.5) → delay. **Output: 7 taps per channel, all at gain 0.6 with signs ±, taken from inside the tank delays and the all-pass delay lines**, 4 from the opposite leg and 3 from its own; mono input, stereo only from the taps | C | Read-out taps on **both sides of the damping filter** (K4) and **inside all-pass lines** (K2b); **equal-magnitude signed taps** like the 246's cluster (equal magnitude, signs only). The figure-eight coupling and in-loop all-passes do not apply (harmonic modes) |
| **Griesinger 1989**, "Practical processors and programs for digital reverberation" (AES 7th Int. Conf.) [P10] | About 100 delay taps per sample adequate at 8 kHz bandwidth; some coefficients need only 3 bits plus sign; 32–64 K words of memory, managed so that it "appears as a continuous shift register" (base + offset); from the BBC (Axon/Gilford) work: **feedback from more than one tap is unstable or coloured, so feed back from a single tap and add density with non-feedback taps**; Lexicon "tailors the shape and spread of the beginning of the diffused reverberation"; most structural detail unpublished | C | The 246's single-feedback harmonic loops are this rule. The "tailored beginning" is what the 246's long excitation might be. 100 taps at 8 kHz is the same order as the 246's 50–70 at 12 kHz |
| Lexicon program sizes (224: 100 steps; PCM60/70: 128; 480L: 80 × 4) and "room = tank + one more diffuser and delay per leg, more output taps, some in the pre-delay" | anonymous reverse-engineering page via freeverb3 [P11] | lore | Order of magnitude only: 100–130 steps per sample, like the 246's 182 per DSP |
| **Ursa Major SST-282** [P21] | One ~255 ms buffer, about 15 taps summed as feedback (moved in steps), a separate set of about 8 output taps | S (secondary) | Multi-tap feedback: excluded for the 246 by the harmonic modes |
| AMS RMX16, Quantec QRS, Eventide SP2016, Sony DRE-2000, Yamaha REV1 | No structural documentation found ("over 1,000 taps" for the RMX16 is marketing) | – | None |

### 3.2 Algorithms from the literature

| Source | Idea | Conf. | How it could give a dense, fixed, per-loop, per-channel feed |
|---|---|---|---|
| **BBC (Axon, Gilford, Shorter, 1955)**, as told by Griesinger 1989 [P10] | Taps added along a delay line, **feedback confined to the last tap**: density goes up with no instability, and the tap pattern repeats unchanged on every pass. Parallel multi-tap combs of incommensurate lengths | C (Griesinger's account) | Exactly K1: a fixed tap pattern per loop on a single-feedback comb. The 246's harmonic, uncoupled loops are what this rule produces |
| **Kuttruff 1962**, "Künstlicher Nachhall", Frequenz 16(3) [P23] | Extra tappings on recirculating loops so that the decay is continuous, not stepped | S (via the EMT 250 patent, not read) | The German root of K1. The 250 patent says EMT simulated it and heard repetition, stepped decay and colouration, fixed with non-harmonic tap spacing, energy-exponential tap gains and randomisation [P3] |
| **Schroeder & Logan 1961; Schroeder 1962** [P15, P16] | Series all-passes; 4 parallel combs (30–45 ms) + 2 all-passes (g 0.7). **Fig. 7 ("ambiophonic"): all-passes first, then the parallel combs**, and a resistor matrix forms up to 16 outputs, each taking every comb output once with sign ±. Quasi-stereo by a pair of different all-passes | C | Schroeder's Fig. 7 is the 246-S topology (common input all-passes → uncoupled parallel loops → per-output signed read-out); the 246 replaces the ±1 read-out with something denser (K1–K3). A different all-pass per output channel is K6 |
| **Gerzon 1971–72**, "Synthetic stereo reverberation" I–II, Studio Sound [P17] | Part I: incommensurate lines cross-coupled by an orthogonal matrix. Part II: the all-pass generalised to any unitary network, and **feedback modifications nested around cascades of unitary networks** (nested all-pass, "delay inside an all-pass") | C | A nested all-pass in front of each loop (K5). His stereo comes from coupling the loops, which the 246's harmonic modes exclude |
| **Moorer 1979**, "About this reverberation business", CMJ 3(2) [P13] | An N-tap early-reflection FIR (7- and 19-tap sets over 0–80 ms, Boston Symphony Hall) whose **whole output feeds the recirculating section**, so "all the density of the first N reflections is then repeated many times"; 6 lowpass combs + 1 all-pass; warns that random tap placement colours | C | Tapped injection (K3), but one FIR for all combs, i.e. not per loop. 40–80 ms is 1,050–2,100 samples at 26.4 kHz, the order of the 246's overhang |
| **Stautner & Puckette 1982**, CMJ 6(1) [P14] | Early reflections **written into each loop's delay at several positions**, 6–8 differently delayed, low-passed copies per pipe, a different set per channel; "the length of each delay determines the amount of the early response that can be simulated"; one read per pipe; unitary 4 × 4 coupling and slow random modulation | C | **Multi-point writes per loop**: a per-loop, RT-independent FIR of span ≤ Pᵢ for a few writes. With reads past the feedback point (K1) the total span reaches Pᵢ + d. The coupling and modulation do not apply |
| **Kendall & Martens**, US 4,731,848 (filed 1984) [P18] | One multi-tap buffer: 6 first-order taps to the outputs, **12 second-order taps feeding 18 recirculating units** through scalers, each unit directionalised to L/R | C | Different taps of one line feed different loops (K3), sparsely |
| **J. O. Smith 1985**, closed waveguide networks (ICMC) [P19] | Lossless scattering junctions, inputs and outputs anywhere; early reflections from a tapped delay line | C | Needs scattering (coupling): excluded by the harmonic modes |
| **Gardner 1992** (just outside the period; credits nested all-passes to Vercoe & Puckette 1985 and to the instrument industry) [P12] | Nested all-passes on one delay line with moving tap addresses; output taps **between** the all-passes; a tap **inside** a nested all-pass gives (1 − g²)/(1 − g·G(z)), resonant and "sounds bad" | C | Where a designer would *not* tap. Dattorro does tap inside all-pass lines, so practice differed |
| Norlin, US 4,215,242 (1978) [P20] | BBD with irregular taps, several summed back as feedback; multiple output summers with different weights on the same taps | C | Per-channel weights on shared taps: one way to get the 246's shared 26 / 332 lags with per-loop, per-channel weights |

### 3.3 German-language sources (EMT's own context)

EMT is a German firm (Lahr), so German patents, brochures and the Kuttruff/Schroeder line were searched separately. Patents were searched in DEPATISnet (DPMA): inventors Bäder/Baeder and Blesser, applicant EMT and Franz Vertriebs, titles `Nachhall?`, and class G10K 15/12, for DE/EP/AT/CH documents published 1978–95. Google Patents rate-limited the searches, and Espacenet and De Gruyter sat behind challenges, so coverage is good but not complete.

| Source | Content | Conf. | Relevance |
|---|---|---|---|
| **EMT patents after the 250** [P24] | None on reverb. Bäder's only other filing is DE 3020247 (1980), a 1-bit DPCM converter; Blesser's only other Franz filing is DE 2842324 (1978), FM sideband suppression | S (search; absence) | The 251, 252 and 246 algorithms were never patented. The 250 patent [P3] stays the only EMT primary text on per-loop taps |
| **EMT 251 brochure** "Digitales Hallgerät" (German, 4 pp.) [P25] | Reverb onset = "eine Kombination von sechs fest zueinander eingestellten Reflexionen", their level sets the "Hallradius", a delay control sets the reverb start; the Delay program reuses "eine Gruppe von sechs Abgriffen mit zueinander fester Zuordnung". **Echo program diagram:** one main feedback loop (t1, A1) feeding **three parallel feedback loops** (t2–t4), each read straight to output bus I or II (one via a pan pot) | C | The fixed-relative tap group is the direct ancestor of the 246's 8-tap cluster (F8: a tap table). The loops feeding loops are in the *Echo* program, not the reverb; for the 246's reverb a series link is excluded by F1 (§5). EMT drew per-channel read-outs of each loop |
| **EMT 245 brochure** (German, 2 pp.) [P26] | "vereinfachte Version des … EMT 250 … nach dem gleichen Programmprinzip"; 18 kHz, 128 kbit RAM, 15-bit A/D, **one time-shared D/A for both outputs** | C | The 246's multiplexed DAC (L half a sample late) is inherited. With the 245 manual's "19 loops in mixed series/parallel" [P4], the 250 → 245 line is tapped loops plus short series stages |
| **DE 3231925 / EP 0101895** (Dittmar, R. Franz, Scheidegger; 1982) [P27] | Microprocessor reverb: **four series all-passes** (1285 / 1120 / 550 / 450 samples at about 10 kHz, k = 0.75 built from shifts); one RAM addressed as a **cyclic base counter minus fixed offsets**, using incomplete address decoding | C | A German design of the same years with the 246's addressing idea, and all-pass gains chosen to be shift-only. No connection to EMT is known |
| **AKG (Frese), AT 378289 / US 4,475,229 / WO 81/03566** (1980–81) [P28] | The prior-art section lists three standard cures for comb periodicity: (1) **tap the loop at several points with different coefficients and sum them**; (2) several different resonators in parallel or series; (3) feedback around those resonators (said to work only in analogue because of digital limit cycles). The invention is a wave-digital ring of ≥ 50 sections of statistically distributed length with input and taps placed symmetrically | C | A period German-language statement that **multi-point read-out of a loop was standard practice**: supports K1. Cure (2) is K2 |
| **Philips (Dijkmans), DE 3015324 and DE 3015357** (1979) [P29] | All-pass sections with the delay split and a second section nested between the sub-delays | C | **Nested all-passes in 1979**, years before Gardner: K5 was available to a 1983 designer |
| Deutsche ITT, DE 2360983 / DE 2360984 (Keller, Adam; 1973) [P30] | MOS delay with feedback, split into sub-lines; claim 5 adds feedback from an intermediate sub-line to the input | C | An early "tap inside the loop fed back", which the 246's harmonic modes exclude |
| Kuttruff, "Künstlicher Nachhall", Frequenz 16(3), 1962 [P23] | Paywalled, not read; its tapped-loop proposal (p. 96) is known only through the 250 patent | – | K1's German root, unread |
| Quantec QRS (W. Schwarz, 1982) | No algorithm patent (his patents concern Quancor synchronisation); reviews only mention "room resonances", a fixed 25 % density and decorrelated outputs | S (secondary) | None |
| Tonmeistertagung 1984/1986, Rundfunktechnische Mitteilungen (IRT), DAGA, Funkschau | No EMT or Bäder reverb paper found online; proceedings exist but no tables of contents are online | – | Library check only |
| Dynacord DRP 20 (1989), AKG ADR 68K (1986) | Secondary only; the ADR 68K is said to be co-developed with Ursa Major | S | None |

**Worth getting through a library (not online):**
- Kuttruff's 1962 Frequenz article;
- Dickreiter, *Handbuch der Tonstudiotechnik* (5th ed., 1987, IRT-linked);
- the Tonmeistertagung 1984 and 1986 tables of contents (for any EMT talk on the 246).

### 3.4 Valhalla DSP (Sean Costello): a modern designer's reading of the period

valhalladsp.com refuses scripted fetches, so the blog posts were read from Wayback Machine captures and the KVR threads directly. Costello's statements about the EMT 250 are his own understanding and listening impressions, not documented sources (**S** at best). Nothing on the blog, in the slides or on KVR concerns the EMT 246, 251 or 252.

| Source | Claim | Conf. | Bearing |
|---|---|---|---|
| "RIP Keith Barr" (2010) [P34], quoting Barr (Alesis/MXR) | "the big loop is great; you inject input everywhere but take it out in only two places"; "Everything I write today is single loop, with modifications as to input signal injection and output tap selection". The **Midiverb had no address counter**: each instruction carried a relative offset, and the offsets summed to one more than the memory size so that a counter "would appear" to exist; a 128-operation machine, most coefficients 0.5 | C (designer testimony, quoted) | **Tapped injection plus selected output taps on one loop was ordinary 1980s practice**: K3 and K1. The Midiverb's offsets are the same idea as the 246's TOM. Coefficients of 0.5 are shift-only (§2.1) |
| KVR post (2012) [P35] | The EMT 250 is unusual among classic algorithms in having "separate modulation processes for the output taps, in addition to the recursive loops" | S (Costello's view) | In EMT's lineage the **output taps are a structure of their own**, separate from the feedback (K1/K6). The modulation does not apply to a static S |
| "Reverbs: diffusion, allpass delays, and metallic artifacts" (2011) [P36]; "Diffusion, vocals and real rooms" (2010) [P37] | The EMT 250 and Lexicon 224 put several series all-passes at the algorithm input; early structures placed the all-passes at the comb inputs, and Gerzon and Griesinger later moved them inside the feedback. The Quantec QRS seems to have no input all-passes, and its all-passes sit inside the recursive network (from impulses he has heard) | S | A common input diffuser is EMT house style (the 246's AP20 · AP72 · AP230). All-passes at the **comb inputs** in early designs fit K2. In-loop all-passes (excluded for S) belong to a later style |
| "Stability through time variation: Ursa Major Space Station" (2010) [P21], with a comment from Chris Moore | One line, 15 taps summed for feedback, a **separate set of 8 output taps**; stable gain limited to about 1 / N without modulation | S (analysis of Moore's patent) | Output taps kept apart from feedback taps (K1-like). Multi-tap feedback is excluded for the 246 |
| AES 2015 reverb presentation slides [P38] | EMT 250 about 300 ms of memory, Lexicon 224 about 1 s; 12-bit "floating point" converters and 16-bit fixed-point arithmetic; algorithms have been trade secrets since the EMT 250; AMS reportedly studied an early 224 before the RMX16 | S | Context. It agrees with Blesser's 8 K words at 32 kHz (256 ms) [P9] |
| "Getting started with reverb design, part 2" (2021) [P39] | Schroeder 1962's transposed form (comb outputs through a matrix for decorrelated outputs); Gerzon's rotation angles for L/R spread; Moorer's sparse FIR; Gardner builds all-passes, nested all-passes and loops in **one memory buffer** | S (summary of papers) | Per-channel output matrix (K6); single-buffer addressing as in the 246 |
| ValhallaVintageVerb modes (2023) [P40] | "Sanctuary", after a "classic German digital reverberator from the 1970s" (the EMT 250), is loosely based on constant-density algorithms with discrete early reflections and an echo density that rises over time. "Palace" follows the QRS, a 1980s room simulator that "didn't have modulated delays" | S | A static German design of the period existed, so a time-invariant S is not unusual. These are inspired-by designs, not structure |
| "Chamber1979" (2023) [P41] | Built from a late-1970s unit's schematics; "nested allpass reverb algorithms" with an "aperture filter"; the unit is not named | S | Nested all-passes were in production by the late 1970s (K5 available), for another maker |

**Not found at Valhalla:** delay lines longer than their loop or feedback from the last tap only (K1 in that exact form), taps on both sides of a damping filter (K4), per-setting tap tables (K7), and the cluster-to-reverb transition in EMT designs.

**Common thread (H):** every period scheme that gives a dense but fixed excitation does it with **taps**: along the loop line (BBC, Kuttruff, EMT 250), into the loop line (Stautner–Puckette), or on a feeding line (Moorer, Matsushita, Kendall–Martens). Recursive pre-stages (Schroeder, the 250's short loops) are used for diffusion in front of the loops. None of the documented schemes needs more than a few dozen taps per sample in total, which fits the 246's budget.

---

## 4. Candidates, ranked

Ranking is by plausibility given §1–§3 and by how cheaply the bench can test it.

**Data for the offline tests.** The saved excitations are in `analysis/q2/branches_flat_0_ext1700.npy`. The layout is as `q1.Joint.split`: 8 cluster gains, the 2 × 1,400 first-pass FIR, then per channel the four xᵢ of Pᵢ + 1700 samples starting at PRE + 230. Every offline test below reads that file and nothing else. A result counts only once the `branches` bench confirms it.

**Common rule for every sparsity test** (from `model_s_results.md` §13.2):
- The excitations are band-limited by B (100 Hz–10 kHz on the bench), so they are undetermined out of band. Score sparsity **in band only**: fit integer taps convolved with the band kernel to the in-band excitation (OMP, or group OMP over L and R), and report the in-band energy explained at fixed tap counts K = 2, 4, 8, 16, 32.
- Include **adjacent-pair atoms** (taps at τ and τ + 1 with gains f and 1 − f) besides single integer taps. The EMT 250's output lines appear to interpolate that way (§3.1), and the AP230 step has the same form. A fractional tap fitted as two free integer taps costs an extra parameter and can hide sparsity.
- Always run the same score on **two controls**:
  - a synthetic with known sparse taps through the known chain (the taps should be found);
  - phase-randomised excitations with the same spectra (should stay dense).
- A candidate counts only if it beats the controls by a wide margin, and only once the bench gives a held-out null with a small parameter count, reported with its excess over the repeat-capture null.

### K1. Tapped read-out along an extended loop line (with or without multi-point writes): rank 1

- **Idea.**
  - Each loop's line is allocated longer than the loop. Feedback is read at Pᵢ, the one feedback point, through the gain and filters. Output taps for each channel read the same line at a few lags, some of them **past Pᵢ** (up to about Pᵢ + 1000).
  - Optionally, the diffused input is written into the line at more than one offset (Stautner–Puckette).
  - Every term is z^−τ · Uᵢ / (1 − Gᵢ z^−Pᵢ) with a fixed gain (§1.2a).
- **Period precedent.**
  - The EMT 250's long loops: feedback only from the end of the line, many read taps along it [P3]. Its 1977 algorithm drawing reportedly shows three read taps per comb feeding short output lines [P31] (forum reading, S).
  - The BBC rule as told by Griesinger: "feedback confined only to the last tap", extra taps for density [P10].
  - Kuttruff's tapped loops [P23].
  - Dattorro's 7 signed taps per channel inside the tank lines [P9].
  - Stautner–Puckette's multi-point writes [P14].
  - In a base + offset memory every such tap is one TOM entry and one IN (§2.2).
  - **C** (precedents) / **H** (for the 246).
- **Mapping.**
  - F1–F3: fixed gains, so RT-independent, and lags past Pᵢ come for free.
  - F2: the overhang length is the longest read lag minus Pᵢ.
  - F5: per-loop, per-channel tap sets. The shared 26 / 332 relation would be the R tap set being, in part, the L tap set offset by 332 (and 26). R's first arrivals are exactly L's delayed by 332 (model §14.6).
  - F7: a few taps give a Rayleigh-like mode spread.
  - F8: the TOM offsets scale ×2 / ×4 in M and L.
- **Budget.** 4 loops × 2 channels × K reads, plus any extra writes. K = 3–4 gives 24–32 accesses, the top of F9. Taps of equal magnitude with signs only (like the 246's cluster and Dattorro's ±0.6) cost 3 cycles each and no multiply. **So K1 predicts very few taps: at most about 4 per loop and channel if they are the whole story.**
- **Why earlier sparse fits do not rule it out.** The sparse-tap results (`model_s_results.md` §8.3, §11.2, §12) all used the short window (Pᵢ + 300). That window cannot hold a tap at Pᵢ + d for d > 300, so OMP had to approximate such taps with many others. **H**: re-run them in the long window before believing "not sparse".
- **Prediction.**
  - In band, each xᵢᶜ is a handful of integer taps on the diffused signal (after `post`), spread over 0 … Pᵢ + 1000, with at least one strong tap past Pᵢ.
  - The R taps of a loop are largely the L taps shifted by 332 or 26.
  - Gains cluster at a few values, possibly one magnitude with signs.
- **Test recipe.**
  1. *Offline, minutes:* the in-band OMP score above on each xᵢᶜ separately, then group OMP over (L, R) with the shifts {0, 26, 332} tied.
  2. *Bench:* fix the selected tap positions and fit their gains jointly across settings. The `branches` stage needs a sparse-tap mode (positions fixed, gains free) in place of the free xᵢ. Run `--subset flat`, held out by RT. Parameters: 8 × K gains plus the first-pass terms.
  3. Then `--subset all`, after K4.
- **Would falsify it.** More than about 8 taps per loop and channel are needed to explain about 99 % of the in-band energy, or the excitations stay no sparser than the phase-randomised control. Then a tapped read-out alone is over budget, and K2/K3 must supply the density.

### K2. A recursive pre-stage per loop, or loops fed from different points of the diffuser chain: rank 2

- **Idea.** Each loop gets its signal through one short recursive section of its own, read and written through the TOM like everything else, after the common diffuser. The section is an all-pass or a comb with a fixed gain and its own delay. Variant K2b: the loops are fed from different points of one longer chain.
- **Period precedent.**
  - Schroeder's all-passes in front of parallel combs [P16].
  - The EMT 250's "short loops with all-pass characteristic" beside its long loops [P3].
  - The EMT 245's 19 loops in "mixed series and parallel" [P4].
  - Lexicon/Dattorro taps inside the input diffusers and tank all-passes [P9].
  - **C** (precedents) / **H** (for the 246).
- **Mapping.**
  - F1/F3: a fixed gain is RT-independent and outside the loop.
  - F2: one stage of delay M and gain g falls to L dB after L·M / (20 log₁₀|g|) samples. M ≈ 200–500 with |g| ≈ 0.45–0.7 gives about 1000 samples at −15…−25 dB (§1.2).
  - F5: different delays per loop give unrelated fine structure. L and R of a loop share the stage.
  - F8: stage delays scale ×2 / ×4 in M and L.
  - Cost: 2 accesses and about 10–14 cycles per stage (§2.3). 4 loops × 1 stage = 8 accesses.
- **Prediction.**
  - For each loop i there is one (M, g), the same for L and R, such that xᵢᶜ deconvolved by that stage is sparse in band.
  - After its last strong tap, the envelope of xᵢ decays in a straight line at 20 log₁₀|g| / M dB per sample, the same in every setting.
  - The autocorrelation of xᵢ peaks at lag M (and 2M, … for a comb).
- **Test recipe (offline, then bench).**
  1. *Envelope and lag screen:* the in-band envelope slope over Pᵢ … Pᵢ + 1000, and the normalised autocorrelation at lags 2–1700. Look for a lag peak shared by L and R of one loop and by no other loop.
  2. *Dictionary deconvolution:*
     - scan M = 2…1700 and g on a grid (±0.3…±0.8, step 0.05);
     - a comb's inverse is the FIR 1 − g z^−M; an all-pass's inverse is the same all-pass run on the time-reversed sequence;
     - score in-band sparsity (above), with the common stages (20, 72, 230, 460 at their fitted gains) among the candidates (a second pass through one of them is K2b);
     - then greedily remove the best stage and repeat for a second one.
  3. *Bench:* add a per-loop factor Eᵢ to `group_spectra` (Rᵢ = PB · Eᵢ / (1 − Gᵢ z^−Pᵢ)), fit tᵢᶜ as short tap sets, run `branches --subset flat` held out by RT. Report taps + 2 parameters per stage.
- **Would falsify it.** No (M, g) beats the phase-randomised control by more than a few dB at K ≤ 16, or the best stages differ between L and R of one loop.

### K3. A tapped injection line per loop, convolved with the tapped read-out: rank 3

- **Idea.** Each loop is written from a short tap set aᵢ on the line that carries the diffused input (Moorer's early FIR feeding the combs, but per loop as in Matsushita's two tap groups and Kendall–Martens' second-order taps). Each channel then reads it with a short tap set bᵢᶜ. So xᵢᶜ = aᵢ ∗ bᵢᶜ. With injection and read sets spanning about 1000 and about Pᵢ, the support is about Pᵢ + 1000 with K_a·K_b impulses from K_a + K_b accesses.
- **Period precedent.** Moorer 1979 [P13]; Matsushita US 4,584,701 [P7]; Yamaha US 4,472,993 [P6]; Kendall–Martens US 4,731,848 [P18]; Keith Barr's "inject input everywhere but take it out in only two places" [P34]. **C** / **H**.
- **Mapping.** F1–F3: all fixed. F5: aᵢ differs per loop (so no shared fine structure across loops) and is common to L and R of a loop, which gives the partial L/R relation. Earlier tests do not cover it: the bilinear model of §10.3 shared one long input across all loops and failed (`model_s_results.md` §10.3). Cost: K_a + 2K_b per loop, e.g. 3 + 2 × 2 = 7, so 28 for four loops.
- **Prediction.**
  - K1's OMP finds a moderate tap count (about 6–20 per loop and channel) whose positions form a **sumset** A + Bᶜ, with A common to L and R of the loop.
  - Equivalently, Xᵢᴸ / Xᵢᴿ = Bᵢᴸ / Bᵢᴿ is a ratio of two very short tap polynomials.
- **Test recipe.**
  - On K1's selected taps, look for repeated difference vectors: a histogram of pairwise lag differences, shared between L and R.
  - Then a joint sparse factorisation of the (L, R) pair with a common aᵢ: alternating OMP on aᵢ and bᵢᶜ, scored in band against the same controls.
- **Would falsify it.** The tap positions show no sumset structure and the L/R ratio is no sparser than K1's taps themselves.

### K4. Read-out taps on both sides of a loop filter placed partway along the line: rank 4 (for HF/LF, not for the overhang)

- **Idea.** As in the Lexicon/Dattorro tank, the damping filter sits partway along each loop line, and output taps sit on both segments [P9]. Taps downstream of the filter hear it once more than those upstream. On flat settings the filter is unity, which is why F1 holds there; once HF/LF move it does not (F4).
- **Mapping.**
  - It explains F4 and the earlier short-window need for "a branch through the loop gain/filter" (`model_s_results.md` §14.4).
  - It does **not** explain the overhang past Pᵢ, so it combines with K1–K3.
  - Cost: no extra accesses (the read taps exist anyway).
- **Prediction.** On the `all` subset with the long window, a second branch through the **filter only** (Hᵢ / Hᵢ(DC), no gain) closes most of the gap between F4 (−29 dB) and the per-setting free injection. A gain-only branch does not.
- **Test recipe (one bench run per split).**
  - Run `branches --subset all --ext 1000` with a filter-only branch.
  - `branch_spectra` currently offers only powers Gᵢ^α of the whole loop factor. The needed forms exist in q1's `group_spectra` (`through_g="filter"`, `"gain"`, `"all"`) and have to be passed through.
  - Hold out by the RT, HF and LF splits. Judge against the per-setting free injection's excess (the LF/HF law errors set a floor there), not against the repeat-capture ceiling.
  - At the strongest settings, part of the residual may be fixed-point error in the loop filters (§5, [P33]). If the filter branch closes the gap except at LF 0.1 / HF 0.2, that is not by itself evidence against K4.
- **Would falsify it.** The filter branch adds less than the gain branch, or nothing.
- **Status:** not supported. A filters-only branch gains about 1 dB on the full library (`model_s_results.md` §15.8; §7.6).

### K5. A nested all-pass per loop: rank 5

- **Idea.** An all-pass with a second all-pass inside its delay (Gerzon 1972's nested feedback modification [P17]; later Gardner [P12]). Two lines, 4 accesses, a dense response with a smooth envelope.
- **Mapping and prediction.** As K2 but denser for the cost. The structure is still all-pass, so its inverse is the time-reversed structure and deconvolution is exact. The autocorrelation shows M₁, M₂ and M₁ ± M₂.
- **Test recipe.** Extend K2's dictionary greedily: take the best single stage as the outer one and scan the inner over the lags the autocorrelation flags. It counts only if it beats the best single stage by clearly more than two extra parameters buy.

### K6. A per-channel output decorrelator shared by all loops: rank 6

- **Idea.** Stereo from one mono reverberator through a different all-pass per output (Schroeder & Logan 1961 [P15]), or the EMT 250's "discorrelator and distributor" [P3].
- **Mapping.**
  - Fits: the 26 / 332 lags are common to every loop.
  - Against: their weights differ per loop, and the strict "loop sum → one distributor" was rejected with the short window (`model_s_results.md` §5).
  - So at most a partial decorrelator next to per-loop read taps.
- **Prediction.** If the channel difference is a common filter Dᶜ, then X_iᴸ · X_jᴿ = X_iᴿ · X_jᴸ for every pair of loops (cross-ratio).
- **Test recipe (offline).** The cross-ratio residual over all six loop pairs on the long-window excitations, against random pairings across files.

### K7. Per-setting tap or coefficient tables for HF and LF: rank 7

- **Idea.** The Z80 loads offsets or gains per setting, as it already does for the AP230 damping step [model §14.3]. Yamaha's 1981 and 1985 patents keep per-program tap and gain tables [P6, P8].
- **Mapping.** F1 says nothing moves with RT on flat settings, so tables could only act on HF/LF (F4).
- **Prediction.** The per-setting free injections (`analysis/q2/v/`) change in discrete steps that line up with the HF or LF step index, concentrated in a few taps or one stage gain.
- **Test recipe.** Add a step-index basis to the `regress` stage. Do this only after K4, since the LF/HF law errors confound it.
- **Status:** not supported. The flat-trained frozen injection predicts all 128 settings with the plain laws (`model_s_results.md` §15.9; §7.6).

### Not ranked (period practice, already tested here)

| Idea | Period precedent | Status here |
|---|---|---|
| A **long** line (≥ 300 ms) with many read taps, EMT 250 style | 250 patent [P3] | Negative: no long loop at lags 3,500–16,000; long-line taps do not help (`model_s_results.md` §10.3, §12.3). K1 differs: the extension past Pᵢ is only about 1000 samples, on the loop's own line |
| One tap set shared by all loops | Ursa Major multi-tap buffer [P21]; Moorer's one FIR for all combs [P13] | Negative: 4–64 shared delays explain ≤ 13 dB (§12.1) |
| Free sparse integer taps per loop | Moorer's FIR, the 250's early path | Negative **in the short window only** (§8.3, §11.2); re-tested by K1 |

---

## 5. What this rules out

- **Chained or series long loops, and loop-sum feedback.** If loop j fed loop i, loop i's residues would contain 1/(1 − Gⱼ z^−Pⱼ) at loop i's poles. That changes strongly between RT steps (Gⱼ goes from 0.70 to 0.87 between RT 1.0 and 2.5). F1 shows xᵢ is the same at every RT, and the modes stay exactly harmonic. So the EMT 245's series loops can survive only as short stages *in front of* the long loops (K2). **S** (with F1's narrow held-out set). Figure-eight tanks, FDN coupling (Gerzon I, Stautner–Puckette, Smith) and in-loop all-passes were already excluded by the harmonic modes.
- **More than one feedback tap per loop.** That would put extra lags into the loop's characteristic polynomial and break the harmonic series; Griesinger's BBC history gives the same rule [P10]. **S**
- **A free dense FIR per loop.** At about 4 cycles per weighted tap (§2.3), 2 × 4 × ~2,500 taps is two orders of magnitude over the whole machine. On-chip FIRs (2 cycles per tap) are capped by the 144-word RAM and the cycle count. **C** (arithmetic)
- **A dense tapped read-out (K1 with tens of taps per loop).** The same arithmetic: K1 is viable only if the long-window excitations are sparse (≤ about 4–8 taps per loop and channel). **C** (arithmetic) / **H** (the split of the budget)
- **One long line shared by every loop, read at fixed offsets.** Already negative twice (F5); the overhang is per loop and unrelated across loops. **S**
- **Noise-like ±1 tap sequences.** Sign-only taps are cheap (no multiply), but a dense sequence still needs one access per tap. **S**
- **RT-driven injection tables.** F1: nothing in the injection moves with RT on flat settings. Tables, if any, act on HF/LF (K7) or on AP230 (already found). **S**
- **The overhang as a loop-gain effect.** F3: with the long window the G branch adds nothing; Q1's "x0 + G·x1" was the short window's approximation of fixed taps past Pᵢ. **S**
- **Component tolerances and other imperfections as the source of the excitation.** These have been looked at separately (unit facts §2, §2b; `analog_results.md`). None of them can produce a fixed, dense excitation that differs per loop:
  - **Analogue parts** (input and output filters, S&H stages, coupling capacitors, the clock crystal) sit before the digital section (one mono path) or after it (one path per channel). They act the same on every loop, so they cannot make loop-specific structure. They are already in the model per unit:
    - the S&H pole p (0.301 Tomkin, 0.325 CLA);
    - the filter-skirt scale (+0.25 % / −0.44 %);
    - the sub-20 Hz high-pass;
    - a clock ratio per capture (a few ppm in the library, about 98 ppm on CLA).
  - **Inside the digital section** nothing drifts or ages: the same input gives the same bits on every pass. The only departures from ideal arithmetic are deterministic fixed-point effects, all open under Q6/Q10:
    - 16-bit storage of every delay-line word, rounded or truncated at each store (SACH keeps the high word; EMT's 245 manual says rounding at every multiply [P4]);
    - coefficient grids (Q12 for MPYK, Q15 from RAM);
    - possible saturation and limit cycles.

    These errors are signal-dependent noise. The spec puts the reverb's own noise about 15 dB above the DELAY program's (73 against 88 dB S/N, unweighted, re nominal), so they lie near or below the −50…−60 dB null targets, and they matter for bit-exactness rather than for structure.

    **A designer's view of the same point** (Sean Costello, Valhalla DSP, forum post [P33], S):
    - old 16-bit fixed-point reverbs "spend a fair amount of time outside of the linear domain", so the order of EQ and reverb matters on them, unlike an LTI model;
    - post-filtering hides artifacts the hardware itself generates;
    - small coefficient word lengths make very low frequencies produce a little high-frequency noise ("put in a 20 Hz sine" on a 480L), and the noise falls a lot on later boxes.

    For the 246 this suggests (**H**):
    - The downstream LP (Studio Magazin 96) is the period-typical place.
    - The strong-setting misfits (LF 0.1 at 125–500 Hz, HF 0.2 at 8 kHz; `model_s_results.md` §14.6, §15.3–15.4) may include fixed-point error from the loop filters with heavy damping, besides law and form errors. Example: the LF shelf's pole 0.9167 with 16-bit state.
    - Such error is deterministic for a given input, so repeat captures of one setting would share it: the tail residuals of different captures do correlate at +0.89…+0.996.
    - It depends on the stimulus, and the sweep is unknown. So it can only be **bounded** (by the spec's 73 dB S/N and a fixed-point simulation of the loop filters), not nulled.
    - The K4 and K7 tests should be judged against that bound.
  - **Recirculation ("repeats")** is exact apart from that rounding: a loop pass does not degrade the signal the way a tape or BBD generation would. The loops' measured behaviour (harmonic modes, decay laws holding to < 1 % in RT) shows no sign of anything else.
  - **Time variance** is a designed feature in M and L ("stochastic control") and absent in S at the library's −57 dB repeat-capture level. The CLA unit's L-only ~1 % read-out jitter is unresolved (unit or capture chain).

  **S** (analogue and clock parts, measured) / **H** (the fixed-point size: the noise model is not built yet).

This note does **not** rule out K1–K7.

---

## 6. Order of work on the bench

*Steps 1–3 and the offline half of 5 are done (§7); the bench steps are listed in §7.8.*

1. **K1 sparsity, offline**, on `branches_flat_0_ext1700.npy`, with both controls. Minutes. This is the decisive first test: sparse means tapped read-out (and maybe K3); dense means a recursive generator (K2/K5).
2. **If sparse:** K3's sumset check on the selected taps, then the K1 bench (fixed positions, gains fitted across settings, held out by RT).
3. **If dense:** K2's envelope/autocorrelation screen and dictionary deconvolution, then K5. Then the K2 bench with a per-loop Eᵢ.
4. **K4** on the `all` subset with the long window and a filter-only branch, one run per split.
5. **K6** cross-ratio (offline), then **K7** in the `regress` stage if HF/LF dependence remains.
6. **Fixed-point side checks** once gains are resolved:
   - whether fixed gains sit on the MPYK Q12 grid (1/4096) rather than on Q15 (§2.1);
   - whether read-out taps share one magnitude (no multiply).
7. **Cross-program check (F8).** Whatever structure wins on S, with its delays scaled ×2 / ×4, should predict the first recirculation of M and L, where they are still time-invariant (`ml_structure.md` §4). The time variance does not touch that hold-out.

Report every bench result as the repo requires: held out, compensated, with the excess over the repeat-capture null and the parameter count.

---

## 7. First results: offline tests on the frozen flat injection (in-sample, pending the bench)

**What was run.**
- Tool: `.venv/bin/python tools/q2_candidates.py` → `analysis/q2/k1_sparsity.{json,log}`, `k2_stages.{json,log}`, `structure.{json,log}`. Runtime about 1 minute.
- Input: the frozen flat injection `analysis/q2/branches_flat_0_ext1700.npy`, 600 iterations, held-out flat nulls −54.3 / −54.0 / −48.8 dB (`model_s_results.md` §15.7).
- Scoring: every score is in band, using the bench's fit band |LFBAND|² (20 Hz–11 kHz) and band-limited atoms. Nothing is fitted to captures here.
- **Status: all of this is in-sample.** It describes the 27.9 k free coefficients. It is not evidence of the algorithm until the bench reproduces the held-out nulls with the compact model (§7.8).

### 7.1 K1: the excitations are sparse (S, in-sample)

In-band energy explained at K atoms (range over the 8 excitations: 4 loops × L, R):

| K | 1 | 2 | 4 | 8 | 16 | 64 |
|---|---|---|---|---|---|---|
| Excitations, integer taps | 16.5–19.9 % | 33–38 % | 66–69 % | 92.0–93.6 % | **98.3–99.4 %** | 98.6–99.5 % |
| Excitations, with adjacent-pair atoms | same, except 1587: 22–23 / 39–44 / 71–75 / 96–99 % | | | | 98.3–99.4 % | 98.6–99.5 % |
| Control: phase-randomised (same spectrum) | 0.6 % | 1.0–1.2 % | 2.0–2.2 % | 3.7–4.0 % | 6.6–7.2 % | 19.7–21.9 % |
| Control: 4 known taps + 1 % noise | – | – | 99.0 %, all 4 positions found exactly (fractional ones included), every run | | | |

- **The sixteen atoms are ten taps:**
  - **5 main taps of one magnitude** (0.28–0.30 in the bench's units), signs only;
  - **5 partners, each 460 samples after a main tap**, at about −0.44 of it.

  The partners are a factor (1 − g z⁻⁴⁶⁰) on every loop path, which cancels the bench's comb460 (§7.3). Adjacent-pair atoms help only in loop 1587, where a partner lands one sample from a main tap (567 / 568).
- **Divided by (1 − g z⁻⁴⁶⁰)**, five integer taps explain **97.96–99.25 %** per excitation. The best g is 0.435–0.44: the scan is flat to ±0.02, and the bench's comb gain is 0.445. The next atom is 25–35 dB below the main taps.
- **R is L shifted, tap by tap:** every R tap is an L tap delayed by **26 or 332** samples, with the same sign. The split between the two lags differs per loop (26 / 332 counts):

  | Loop | 1321 | 1389 | 1461 | 1587 |
  |---|---|---|---|---|
  | 26 / 332 | 2 / 3 | 3 / 2 | 2 / 3 | 3 / 2 |

  That is the earlier "lags 26 and 332 with per-loop weights" (`model_s_results.md` §5), resolved.
- **Tied model.** One magnitude a for all 40 taps, L taps and signs from OMP, R = L + 26 or 332, and the comb factor:

  > xᵢᶜ = a · Σₖ sₖ δ[tₖ + c·uₖ] · (1 − g z⁻⁴⁶⁰)

  It explains **98.66 %** in band overall (residual −18.7 dB re the excitations; 98.16–99.34 % per excitation). Parameters: 20 positions, 20 signs, 20 shift labels, a, g460.
- **The L taps (lag after cluster tap 0, sign) and R shifts:**

  | Loop | L taps | R shifts | Span |
  |---|---|---|---|
  | 1321 | 276 −, 913 −, 971 +, 1256 −, 1514 + | 332, 26, 332, 332, 26 | 1238 |
  | 1389 | 624 −, 768 −, 826 −, 1129 −, 1449 − | 26, 26, 332, 26, 332 | 825 |
  | 1461 | 244 +, 677 −, 1069 +, 1368 +, 1521 + | 332, 26, 26, 332, 332 | 1277 |
  | 1587 | 338 +, 797 −, 984 +, 1319 +, 1780 − | 332, 26, 332, 26, 26 | 1442 |

- **The "overhang" is explained:** every loop's five taps span **less than one loop length**. What ran past Pᵢ in the free excitations was the +460 partners and R's +332 shifts. Reads past the feedback point (§1.2a) are not needed.
- **Reading in hardware terms (H):**
  - per loop, 5 reads of the loop line, added or subtracted with no multiply (3 cycles each);
  - the sum goes to L;
  - each tap is also routed to one of **two output delay lines of 26 and 332 samples, shared by all loops**, whose outputs form R.

  That is about 20 reads + 4 accesses ≈ 24 accesses, inside the 20–35 left (F9). It is the EMT 250's tapped read-out plus its "output discorrelator and distributor" (delays t17–t19) [P3], and matches the 1977 drawing's "taps per comb into short delay lines" [P31]. The equal magnitude matches the 246's cluster (equal taps, signs only) and Dattorro's ±0.6 output taps [P9].

### 7.2 K2 / K5: no recursive stage per loop (S, in-sample)

Scan over every delay M = 2…1700 and g ∈ ±{0.3, 0.45, 0.6, 0.7, 0.8}, in four forms (extra comb, fed before a comb, extra all-pass, fed before an all-pass). Each candidate was scored with a fast 5-tap proxy, then by exact OMP on the best few.
- **On the raw excitations, the scan finds the comb factor blind.** In all 8, the best re-expression is "fed before a comb", M = 460, g = +0.45 (the grid point nearest 0.435). It explains 97.96–99.24 % with 5 taps, against 82.1–84.6 % for the raw excitations. Every runner-up is the same M = 460 at another gain.
- **On the comb-divided excitations, nothing helps.** The best candidate is 90.6–93.6 %, below the 97.96–99.25 % of the taps alone.
- **K5 (nested all-pass)** was therefore not scanned: no single stage improves on the taps, and the residual is small and unstructured (§7.5).
- **Reading:** no per-loop recursive generator exists. The one "stage" the data ask for is K2b in its simplest form: the loops are fed from a point that the bench's comb460 does not cover (§7.3).

### 7.3 First-pass FIR and comb460 (in-sample; H for the reading)

- **The first-pass FIR F is almost empty.** Its in-band energy is 0.3 % of a typical excitation's. It is not sparse (34–43 % at 32 taps, raw or comb-divided, against 17 % for its phase control). With the long window, every early arrival is carried by the loop excitations.
- **So no 460 recursion is left in the model's output.** The cluster does not pass `post`, F is negligible, and every loop excitation carries (1 − g z⁻⁴⁶⁰) with g within ±0.02 of the bench's 0.445.
- **Implication to test (H):** the "fifth recursion at 460 samples" (`EMT246_unit_facts.md` §5; `model_s_results.md` §3.3, §8.2, §11) may be a stand-in from the short-window era, not part of the unit. On the bench, removing comb460 from `post` should then change nothing in the nulls (with the long window and the excitations refitted). That is a cheap, decisive check. Until it runs, the unit-facts row stays as it is.

### 7.4 K3: no structure among the tap positions (in-sample)

- No pairwise difference repeats within a loop (0 of 40).
- Two differences are shared across loops (58: 1321 and 1389; 825: 1389 and 1461), against 0.65 expected by chance (95th percentile 2).
- So the five taps are not a sumset of an injection set and a read set. Each loop is written at one point and read at five. Multi-point injection (K3, Stautner–Puckette writes) is not needed.
- **H, unexplained:** loop 1587 has two tap pairs 459 and 461 samples apart, near the 460 lag. It could be coincidence, or a sign that the tap offsets are tied to the 460 line in the memory layout.

### 7.5 K6 and the residual (in-sample)

- **K6, one per-channel distributor for all loops** (R = (c0 + c26 z⁻²⁶ + c332 z⁻³³²) · L), explains only **48–52 %** of R. With c refitted per loop it explains 51–54 %. The per-tap 26 / 332 labels explain **98.2–99.3 %**. So R is not a filtered copy of L: the routing is per tap, onto two shared delay lines (§7.1). This also explains why §5 of `model_s_results.md` rejected a strict common distributor.
- **Where the tied model's residual sits:**
  - 10–15 % of the residual lies within ±8 samples of a tap or its partner, a window of 4–6 % of the samples. So small gain or position errors contribute somewhat.
  - The rest is spread in time and roughly flat across octaves (−15…−27 dB re each octave of the excitation).
  - That looks like estimation noise in the free fit rather than a missing stage (H). The bench decides: if the tied model's held-out nulls equal the free injection's, the residual was noise.

### 7.6 K4 and K7, from the other Q2 runs

These need the bench. The other Q2 session ran the relevant variants (`model_s_results.md` §15.8–15.9):
- **K4 (read taps on both sides of a mid-line filter): not supported (S).** A branch through the filters only gains about 1 dB on the full library.
- **K7 (per-setting tables for the injection): not supported (S).** Trained on 6 flat settings and frozen, the injection predicts all 128 settings with the plain laws better than any injection trained on the whole library. What remains at LF 0.1, HF 0.2 and short RT is in the loop filters, not the injection.

### 7.7 Scorecard

| Candidate | Result (in-sample unless noted) |
|---|---|
| **K1** tapped read-out | **Holds:** 5 equal signed taps per loop, 98.7 % with one tied model; no reads past the feedback point needed |
| K2 recursive pre-stage | Not needed; the blind scan only finds the comb460 cancellation (K2b: loops fed before comb460) |
| K3 tapped injection × read-out | Not needed: no sumset or shared-difference structure |
| K4 taps around a mid-line filter | Not supported (bench, §15.8: about 1 dB) |
| K5 nested all-pass | Not needed (no single stage helps; not scanned) |
| K6 common per-channel decorrelator | Rejected as a filter (≈ 50 %); replaced by per-tap routing to two shared 26 / 332 delay lines |
| K7 per-setting injection tables | Not supported (bench, §15.9) |

### 7.8 Next on the bench (to turn this into evidence)

1. **Sparse-tap mode in `q2_injection.py --stage branches`:** fixed tap positions, signs and 26 / 332 labels, with a and g460 (or free per-tap gains as a check) fitted on the flat training settings. Report held-out compensated nulls against the free injection's −54.3 / −54.0 / −48.8 dB, with the excess over the repeat-capture null and the parameter count (about 60 discrete + 2 continuous).
2. **comb460 out of `post`** with the long window: nothing should change (§7.3).
3. **All 128 settings** with the tied injection frozen and the laws, as in §15.9.
4. **Cross-program check (F8):** the same taps, delays scaled ×2 / ×4 (with the R lags 52 / 664 and 104 / 1328), against M and L's first recirculation.
5. **Then update** `EMT246_unit_facts.md` §5 (injection/read-out, the 460 row) and the hardware budget (§6) if the bench confirms.

---

## 8. Sources


**Read for this note** = primary text read by the research (scan, OCR or text layer). **Secondary** = description only. Patent texts were read from Google Patents' OCR, which has errors (for example "300 ns" for "300 ms"); check figures against the official PDFs.

| # | Source | URL | Status |
|---|---|---|---|
| P1 | Texas Instruments, *TMS32010 User's Guide* (1985); Rev B 1983 (SPRU001B) in the same directory | https://www.mirrorservice.org/sites/www.bitsavers.org/components/ti/TMS320xx/3201x/TMS32010_Users_Guide_1985.pdf | read |
| P2 | Texas Instruments, *Digital Signal Processing Applications with the TMS320 Family* (1986): ch. 3 Lovrich & Simar, "Implementation of FIR/IIR filters with the TMS32010/TMS32020" (= appnote SPRA003); ch. 11 J. Bradley, "Interfacing external memory to the TMS32010" | https://www.mirrorservice.org/sites/www.bitsavers.org/components/ti/TMS320xx/3202x/1986_TI_Digital_Signal_Processing_Applications_with_the_TMS320_Family_Theory_Algorithms_and_Implementations_198605.pdf ; SPRA003A scan: https://ardent-tool.com/unsorted/TI_TMS320C2x/spra003a.pdf | read |
| P3 | Blesser & Bäder, *Electric reverberation apparatus*, US 4,181,820 (priority DE 2719276, 1977), EMT 250 | https://patents.google.com/patent/US4181820A/en (local PDF in `sources/`) | read |
| P4 | EMT 245 service manual (German), introduction: 19 loops in mixed series/parallel, several processors, 16-bit words with rounding | local: `sources/emt-documents/EMT-245-service-manual-DE.pdf` | read (local) |
| P5 | Nichols (Quad Eight), US 4,219,880 (filed 1978): address counter + per-instruction offset | https://patents.google.com/patent/US4219880A/en | read |
| P6 | Nippon Gakki (Yamaha), US 4,472,993 (priority 1981): counter + ΔADR, early taps feeding combs, per-program tables | https://patents.google.com/patent/US4472993A/en | read |
| P7 | Matsushita, US 4,584,701 (priority 1982): one RAM, two 3-tap groups each feeding 2 of 4 loops | https://patents.google.com/patent/US4584701A/en | read |
| P8 | Yamaha, JPH0683515B2 (1985): tapped registers, per-hall tables in ROM | https://patents.google.com/patent/JPH0683515B2/en | read (machine translation) |
| P9 | Dattorro, J., "Effect Design Part 1: Reverberator and other filters", JAES 45(9):660–684, 1997, incl. Blesser's letter (App. 1) and Table 2 output taps | https://ccrma.stanford.edu/~dattorro/EffectDesignPart1.pdf | read |
| P10 | Griesinger, D., "Practical processors and programs for digital reverberation", AES 7th Int. Conf., Toronto, 1989 | https://freeverb3-vst.sourceforge.io/doc/10.1.1.72.5367.pdf | read |
| P11 | Anonymous Lexicon program-step and "room" notes, collected by freeverb3 | https://freeverb3-vst.sourceforge.io/tips/reverb.shtml | lore |
| P12 | Gardner, W. G., "A realtime multichannel room simulator" / MIT report on nested all-pass reverbs, 1992 | https://freeverb3-vst.sourceforge.io/doc/10.1.1.52.7008.pdf | read |
| P13 | Moorer, J. A., "About this reverberation business", Computer Music Journal 3(2):13–28, 1979 | https://freeverb3-vst.sourceforge.io/doc/aboutThisReverberationBusiness.pdf | read |
| P14 | Stautner, J. & Puckette, M., "Designing multi-channel reverberators", CMJ 6(1):52–65, 1982 | https://www.ee.columbia.edu/~dpwe/e4896/papers/StautP82-reverb.pdf | read |
| P15 | Schroeder, M. R. & Logan, B. F., "Colorless artificial reverberation", IRE Trans. Audio AU-9:209–214, 1961 | https://hajim.rochester.edu/ece/sites/zduan/teaching/ece472/reading/Schroeder_1961.pdf | read |
| P16 | Schroeder, M. R., "Natural sounding artificial reverberation", JAES 10(3):219–223, 1962 | https://hajim.rochester.edu/ece/sites/zduan/teaching/ece472/reading/Schroeder_1962.pdf | read |
| P17 | Gerzon, M. A., "Synthetic stereo reverberation", Studio Sound, Part 1 Dec 1971 pp. 632–635, Part 2 Jan 1972 pp. 24–28 | https://worldradiohistory.com/Archive-Studio-Sound/70s/Studio-Sound-1971-12.pdf ; https://worldradiohistory.com/Archive-Studio-Sound/70s/Studio-Sound-1972-01.pdf | read (OCR) |
| P18 | Kendall & Martens, US 4,731,848 (filed 1984) | https://patents.google.com/patent/US4731848 | read |
| P19 | Smith, J. O., "A new approach to digital reverberation using closed waveguide networks", ICMC 1985, pp. 47–53 | https://ccrma.stanford.edu/files/papers/stanm31.pdf | read (OCR) |
| P20 | Norlin (Gross), US 4,215,242 (1978) | https://patents.google.com/patent/US4215242A/en | read |
| P21 | Valhalla DSP on the Ursa Major Space Station SST-282 | https://valhalladsp.com/2010/05/14/stability-through-time-variation-ursa-major-space-station/ | read (Wayback); Costello's analysis of Moore's patent, with a comment by Moore |
| P22 | Blesser, Baeder & Zaorski, real-time audio computer (EMT 250 precursor), JAES 23, Nov 1975, pp. 698–707 | https://aes2.org/publications/elibrary-page/?id=2659 | abstract only |
| P23 | Kuttruff, H., "Künstlicher Nachhall", Frequenz 16(3):91–96, 1962 | https://www.degruyterbrill.com/journal/key/freq/16/3/html | secondary (via P3) |
| P24 | DEPATISnet (DPMA) searches: inventor Bäder / Blesser, applicant EMT / Franz Vertriebs, `Nachhall?`, G10K 15/12, 1978–95. Record pattern `https://depatisnet.dpma.de/DepatisNet/depatisnet?action=bibdat&docid=<DOCID>` (DE 3020247, DE 2842324) | https://depatisnet.dpma.de | searched (absence result) |
| P25 | EMT 251 brochure, "Digitales Hallgerät" (German) | https://bibliothek.echoschall.de/emt-251-digitales-hallgeraet-prospekt | read (scan) |
| P26 | EMT 245 brochure, "Digitales Hallgerät" (German) | https://bibliothek.echoschall.de/emt-245-digitales-hallgeraet-prospekt | read (scan) |
| P27 | Dittmar, Franz & Scheidegger, DE 3231925 / EP 0101895 (1982) | DEPATISnet, search the number (record URL pattern as P24); EP 0101895 on Espacenet | read |
| P28 | Frese (AKG), AT 378289 / US 4,475,229 / WO 81/03566 (1980–81) | https://patents.google.com/patent/US4475229 | read (OCR) |
| P29 | Dijkmans (Philips), DE 3015324, DE 3015357 (1979) | DEPATISnet, search the numbers (pattern as P24) | read (claims) |
| P30 | Keller & Adam (Deutsche ITT), DE 2360983, DE 2360984 (1973) | DEPATISnet, search the numbers (pattern as P24) | read (claims) |
| P31 | Gearspace, "Reverb Subculture" thread, page 9 (posts #242–#255, 2009; S. Costello and others on the EMT 250 drawing; text supplied by the user, the page itself was not reachable by script) | https://gearspace.com/threads/reverb-subculture.380233/page-9 | read (pasted copy); forum |
| P32 | Universal Audio, EMT 250 plug-in video with Barry Blesser (2009); the algorithm drawing is visible at about 2:32–2:36 | https://www.youtube.com/watch?v=slXPcE0ioX4 | not viewed (described in P31) |
| P33 | S. Costello (Valhalla DSP), posts in the Gearspace "Reverb Subculture" thread (page not recorded) on pre- vs post-EQ and fixed-point artifacts in older reverbs (480L, 960L, PCM96); text supplied by the user | https://gearspace.com/threads/reverb-subculture.380233/ | read (pasted copy); expert opinion |
| P34 | Valhalla DSP, "RIP Keith Barr" (2010), quoting Barr | https://valhalladsp.com/2010/08/25/rip-keith-barr/ | read (Wayback) |
| P35 | KVR forum, Costello post on the EMT 250's output-tap modulation (2012) | https://www.kvraudio.com/forum/viewtopic.php?p=5033275 | read |
| P36 | Valhalla DSP, "Reverbs: diffusion, allpass delays, and metallic artifacts" (2011) | https://valhalladsp.com/2011/01/21/reverbs-diffusion-allpass-delays-and-metallic-artifacts/ | read (Wayback) |
| P37 | Valhalla DSP, "Diffusion, vocals and real rooms" (2010) | https://valhalladsp.com/2010/06/07/diffusion-vocals-and-real-rooms/ | read (Wayback) |
| P38 | S. Costello, AES Seattle reverb presentation slides (2015) | http://valhalladsp.com/wp-content/uploads/2015/06/aes2015reverbpresentation.pdf | read (Wayback) |
| P39 | Valhalla DSP, "Getting started with reverb design, part 2: the foundations" (2021) | https://valhalladsp.com/2021/09/22/getting-started-with-reverb-design-part-2-the-foundations/ | read (Wayback) |
| P40 | Valhalla DSP, "ValhallaVintageVerb: the modes" (2023) and the 3.0.0 "Palace" notes | https://valhalladsp.com/2023/02/10/valhallavintageverb-the-modes/ ; https://valhalladsp.com/2023/01/25/valhallavintageverb-updated-to-3-0-0-new-reverb-mode-palace/ | read (Wayback) |
| P41 | Valhalla DSP, "Chamber1979 and Hall1984" (2023) | https://valhalladsp.com/2023/12/13/valhallavintageverb-4-0-0-new-reverb-modes-chamber1979-and-hall1984/ | read (Wayback) |

Corrections found on the way:
- US 4,955,057 is **not** a Griesinger/Lexicon patent: it is Tominari (Dynavector), an analogue all-pass dispersion delay.
- No Lexicon patent on reverb structure was found; Griesinger's patents are on LARES.
- The EMT 250 tap rule "Tₙ = ⁿ√k" depends on patent OCR; `algorithm_research.md` §1.1 checked it against the PDF image.
- No TI document or 1983–88 paper on reverb code for the TMS320 was found.
- "Vercoe & Puckette 1985, *Synthetic spaciousness*" could not be found as cited; the 1985 ICMC paper found is "Synthetic rehearsal".
