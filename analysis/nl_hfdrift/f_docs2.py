"""Docs edit 2 (HF chain term, best-of-both estimator): applied to the working files and, separately, to HEAD's versions (the staged blobs)."""
import subprocess
OLD_N = ("**But 320 of 2,088 channel-bands at 10–12 kHz get worse by more than 0.5 dB (worst +13.9 dB)**, mostly captures that already nulled well (−47 → −42 dB), clustered in SPAN 0.3 (115 / 332), SPAN 0.1 (86 / 350), the trimmed captures (123 / 404, gain only 1.2 dB) and the (0.2, −0.2, *) family; so it is left **opt-in, default off**. ")
NEW_N = ("**With that first estimator 320 of 2,088 channel-bands at 10–12 kHz get worse by more than 0.5 dB (worst +13.9 dB)**, mostly captures that already nulled well, clustered in SPAN 0.3 (73 / 230), SPAN 0.1 (65 / 306) and the trimmed captures (123 / 404). "
"Why (`g_summary.txt`, `h2_est.log`): a perfect per-capture scalar (fitted on the capture's own reverb, an oracle) gains 6.1 dB and harms none, so the shape is right and the harm is the estimate of the scalar: |s_used − s_best| is 0.035 in the failing captures against 0.009 in the others, 19 % of the failing ones have neighbours from another capture run (8 % of the rest), and where they fail the front-tap estimates of the 15 nearest scatter widely; s changes between runs and little within ten minutes (median |Δs| 0.001). "
"A free complex gain on 3 or 9 knots per capture gains a further 2–3.5 dB beyond the scalar, so the term has more freedom than one scalar (open). "
"**Best of both (`NL_HF_FILE=analysis/nl/hf_corr_g5.json`, `analysis/nl_hfdrift/f_table2.py`): the scalar from the front taps weighted by a 5-minute Gaussian in time, and no correction where the 15 nearest front-tap estimates scatter above the library's 90th percentile (0.0164; front-tap data only, 205 of 1,038 captures skipped).** "
"Full bench `bench_warp2_bk_kf_a0.7_lb2_hf_g5`: 10–12 kHz −48.3 / −49.0 dB (median), worse than the default by more than 0.5 dB in 101 of 2,088 channel-bands (4.8 %; worst +8.2 dB; SPAN 0.3 none, SPAN 0.1 9.5 %, trimmed 4.0 %, the rest 4.9 %), the skipped captures unchanged; every window 0–2 s and every band 20 Hz–10 kHz unchanged (worst channel-window +0.49 dB, worst band +0.20 dB). "
"It stays **opt-in** (`NL_HF=1 NL_HF_FILE=…`): making it the default is a user decision. ")
OLD_U = "opt-in `NL_HF=1` takes 10–12 kHz from −40.8 / −39.9 to −48.1 / −48.4 dB (median) but 15 % of channel-bands get worse, so it is off"
NEW_U = "opt-in `NL_HF=1` with the Gaussian-weighted scalar and a front-tap skip rule (`hf_corr_g5.json`) takes 10–12 kHz from −40.8 / −39.9 to −48.3 / −49.0 dB (median), 4.8 % of channel-bands worse by more than 0.5 dB (15 % with the first estimator), nothing below 10 kHz moves; off by default"
def fix(old, new):
    def f(s):
        assert old in s, old[:40]; return s.replace(old, new, 1)
    return f
for path, fn in (("docs/NonLin_program.md", fix(OLD_N, NEW_N)), ("docs/EMT246_unit_facts.md", fix(OLD_U, NEW_U))):
    work = open(path).read(); head = subprocess.run(["git", "show", "HEAD:" + path], capture_output=True, text=True, check=True).stdout
    open(path, "w").write(fn(work))
    blob = subprocess.run(["git", "hash-object", "-w", "--stdin"], input=fn(head), capture_output=True, text=True, check=True).stdout.strip()
    subprocess.run(["git", "update-index", "--cacheinfo", f"100644,{blob},{path}"], check=True); print(path, blob[:10])
