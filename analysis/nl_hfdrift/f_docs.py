"""Apply the docs edit for the 11-12.8 kHz chain term to the working file AND to HEAD's version (the staged blob), so the commit
carries only this edit, not another session's uncommitted hunks."""
import subprocess, sys
NEW = ("- **11–12.8 kHz after the warp fix: a term of the capture chain, per capture run (S; opt-in `NL_HF=1`, 27 October 2026, "
"`analysis/nl_hfdrift/REPORT.md`).** Headroom: the repeat ceiling (576 LATE-only pairs, alignment fitted on 8–9 kHz so every sub-band is held out) "
"is −56.5 / −56.0 / −52.4 / −40.4 dB at 9–11 / 11–11.6 / 11.6–12.3 / 12.3–12.8 kHz against the model's −52.6 / −42.2 / −36.0 / −28.4: 15–19 dB above 11 kHz, "
"model-side, not noise (62–76 dB down) and not repeat-limited; the 10–12 kHz ceiling is −55.8 dB. What is left is the same in every time window after 0.05 s, "
"SPAN, LATE, channel and trimming (the chain and its regularised inverse, not the model), and is measured from the front taps alone (R: −1.5 / −5.7 / −12.1 / −20.8 / −30.3 / −40.0 mrad "
"at 11.5 / 11.75 / 12 / 12.25 / 12.5 / 12.8 kHz, +0.03–0.09 dB). It follows the **capture run**, not SPAN: the WAV `bext` chunks carry the origination times (1,038 of 1,046 files, 9–12 April 2016), "
"which split the library into 12 runs; one scalar per run on a rank-1 shape of 9 complex knots explains 98 % of the run-to-run variance (a band-edge filter shift: not a delay, the best per-run delay leaves 28 % of the phase variance; "
"not the clock, ppm correlates +0.5 across runs and with changing sign within them). Re-fitting o, the front shift and ppm on a high band fails: it costs 5–10 kHz 4–9 dB for at most 2 dB at 11–12.8 kHz; the bench's alignment is optimal to 1·10⁻⁴ samples there, "
"and the +0.004-sample high-band delay is the same on the front taps and the branch (a chain phase). **`NL_HF=1`** (`analysis/nl/hf_corr.json`, built by `analysis/nl_hfdrift/f_table.py`) multiplies the model by exp(z) above 10 kHz, z = shared mean + scalar × shape at 9 knots (10.5–12.8 kHz; shape fitted on the other runs, scalar from the R front taps of the 15 captures nearest in time, the capture's own front taps excluded; 36 shared numbers + one scalar per capture run; "
"a measurement-side term like the clock and the gain, not the unit). Full 1,044-capture bench against the current default (`bench_warp2_bk_kf_a0.7_lb2_hf` vs `…_lb2`): **10–12 kHz −40.8 / −39.9 → −48.1 / −48.4 dB (median)**; every window and every band 20 Hz–10 kHz unchanged (worst channel-window +0.49 dB, worst band +0.20 dB). "
"**But 320 of 2,088 channel-bands at 10–12 kHz get worse by more than 0.5 dB (worst +13.9 dB)**, mostly captures that already nulled well (−47 → −42 dB), clustered in SPAN 0.3 (115 / 332), SPAN 0.1 (86 / 350), the trimmed captures (123 / 404, gain only 1.2 dB) and the (0.2, −0.2, *) family; so it is left **opt-in, default off**. "
"A pooled single correction (22 numbers from the reverb, or 20 from the front taps) gains only 1–4 dB with 22–30 % of captures worse.\n")
UF_ROW = ("| **11–12.8 kHz chain term:** the leftover above 11 kHz after the warp fix is a band-edge filter shift of the capture chain that follows the capture run (12 runs from the WAV timestamps; one scalar per run on a rank-1 shape), measurable from the front taps alone; "
"opt-in `NL_HF=1` takes 10–12 kHz from −40.8 / −39.9 to −48.1 / −48.4 dB (median) but 15 % of channel-bands get worse, so it is off | S (measured) / H (cause: the capture chain's analog edge, run to run) | NonLin §11.6 |\n")
def edit_main(s):
    a = s.index("### 11.7 Second capture")
    return s[:a].rstrip("\n") + "\n" + NEW + "\n" + s[a:]
def edit_uf(s):
    L = s.split("\n"); out = []
    for l in L:
        out.append(l)
        if l.startswith("| **Top end:**"): out.append(UF_ROW.rstrip("\n"))
    return "\n".join(out)
for path, fn in (("docs/NonLin_program.md", edit_main), ("docs/EMT246_unit_facts.md", edit_uf)):
    work = open(path).read(); head = subprocess.run(["git", "show", "HEAD:" + path], capture_output=True, text=True, check=True).stdout
    open(path, "w").write(fn(work))
    blob = subprocess.run(["git", "hash-object", "-w", "--stdin"], input=fn(head), capture_output=True, text=True, check=True).stdout.strip()
    subprocess.run(["git", "update-index", "--cacheinfo", f"100644,{blob},{path}"], check=True)
    print(path, "staged blob", blob[:10])
