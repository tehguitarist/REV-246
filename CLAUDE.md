# Agent guide: REV-246 (EMT 246 reverb reverse engineering)

**Read first:**
1. `docs/NEXT_STEPS.md`: current status, working rules and the open questions.
2. `docs/library_notes.md`: the Tomkin IR library, its naming and harness results.
3. `docs/EMT246_unit_facts.md`: the single source of truth for findings.

The full file map is in `README.md`.

## Hard constraints

- **No hardware access, ever.** There is no unit, no ROM/PROM dump and no new captures; repeat captures of a setting only exist where libraries already contain them. Evidence comes from existing IR libraries, documents and photos. Don't propose capture sessions or ROM dumps.
- **Accuracy bar:** the user wants a null of −50 to −60 dB against hardware IRs at minimum, ideally near bit-identical. Always report nulls:
  - as **held-out** predictions (time windows or settings not used in the fit);
  - with **excess over the noise floor**;
  - with the parameter count.

  An in-sample null from many free parameters is not evidence of the true algorithm.
- **Tomkin library: work on the 96 kHz set only.** The 48 kHz files are the same capture rendered at another rate (`docs/library_notes.md` §5), so they carry no extra information and must not be used as a held-out set or as repeat captures. Leave `--separate-rates` off. Files are trimmed and faded per setting: take RT from the label, not the file.
- **Efficiency matters.** The original ran on 2 × TMS32010 (about 50–70 weighted delay taps per sample). Prefer compact structural models over free FIRs.

## Conventions

- **Label claims:**
  - **C** = confirmed (primary source or direct measurement);
  - **S** = strong (several measurements agree);
  - **H** = hypothesis.

  When a finding changes, update `docs/EMT246_unit_facts.md`.
- **Device rate:** 19.2 MHz/728 = 26,373.6 Hz. All delays are integer device samples.
- **Mono in, stereo out.** One mono input feeds the whole algorithm, with separate L and R outputs (one multiplexed DAC; L is half a sample late). A stereo IR file is therefore the complete 1-in/2-out response for that setting. Model the structure as shared up to the read-out, with per-channel output taps.
- **Python:** use `.venv/bin/python`. Third-party downloads go in their own folders and are run with `-I`.
- **Don't commit** `sources/` or third-party IR libraries unless their licence allows it.
- **Git:** end commit messages with the attribution line given in the session's system reminder.
