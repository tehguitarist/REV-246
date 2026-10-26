"""Data defects in the Tomkin S library, kept out of fits, repeat references and tail checks (library_notes.md §2).

    BAD      four AMPL-0 captures that null against their own setting's other captures at only -10...-26 dB (model_s_results §26.1)
    SLIPPED  13 files with a sample slip at 99.9 ms file time, a capture-processing artefact (model_s_results §28)
    GLITCH   four more files with a waveform departure at the same 99.9 ms file time but no lag step (model_s_results §35.2,
             tools/loop_read/glitch_scan.py): 1.0_0.3_0.4 AMPL 0 and 90, 2.0_0.3_0.4 AMPL 30, 2.0_0.8_0.2 AMPL 30
    DEFECTS  BAD | SLIPPED | GLITCH (21 files), as (setting key, AMPL) pairs; in DEFECTS since model_s_results §35.4
             (1.0_0.3_0.4 is left with no clean capture: AMPL 30 / 60 SLIPPED, 0 / 90 GLITCH; the tools skip it)
    ML_SLIPPED  five M files with a one-sample lag step (-1.00 file samples, both channels) at 126-139 ms file time,
             before M's first recirculation (ml_structure §23, tools/ml_slip_scan.py), as (program, setting key, AMPL);
             none in L. Left out of the M / L pulse report and its medians"""
BAD = {("0.1_0.3_0.2", 0), ("1.0_0.6_0.6", 0), ("1.0_0.6_1.0", 0), ("0.4_1.0_0.4", 0)}
SLIPPED = {("0.1_0.3_0.2", 60), ("0.1_0.3_0.4", 60), ("0.4_0.3_0.2", 60), ("0.4_0.3_0.4", 30), ("0.4_0.3_0.4", 60),
           ("1.0_0.3_0.2", 30), ("1.0_0.3_0.2", 60), ("1.0_0.3_0.4", 30), ("1.0_0.3_0.4", 60), ("1.0_0.6_0.2", 30),
           ("1.0_0.6_0.2", 60), ("2.0_0.3_0.2", 30), ("2.0_0.6_0.2", 30)}
GLITCH = {("1.0_0.3_0.4", 0), ("1.0_0.3_0.4", 90), ("2.0_0.3_0.4", 30), ("2.0_0.8_0.2", 30)}
DEFECTS = BAD | SLIPPED | GLITCH
ML_SLIPPED = {("M", "0.7_3.0_0.4", 10), ("M", "1.2_1.8_0.4", 60), ("M", "2.0_0.6_0.7", 60), ("M", "2.0_0.8_0.4", 60),
              ("M", "2.0_1.2_0.4", 60)}
NL_DEFECTS = {(0.0, -0.5, 0.0),    # Non-Lin (SPAN, EARLY, LATE): a 0.27-sample slip between front taps 410 and 463 (NonLin_program 9.7)
              (0.1, -1.0, 0.0)}    # a foreign tail: a -52 dB plateau (re the early section) decaying ~17 dB over 5 s, ACF 537,
                                   # in a gated (LATE 0) capture whose neighbours end at 0.55 s (NonLin_program 11.1)
