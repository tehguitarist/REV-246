"""Data defects in the Tomkin S library, kept out of fits, repeat references and tail checks (library_notes.md §2).

    BAD      four AMPL-0 captures that null against their own setting's other captures at only -10...-26 dB (model_s_results §26.1)
    SLIPPED  13 files with a sample slip at 99.9 ms file time, a capture-processing artefact (model_s_results §28)
    DEFECTS  both, as (setting key, AMPL) pairs"""
BAD = {("0.1_0.3_0.2", 0), ("1.0_0.6_0.6", 0), ("1.0_0.6_1.0", 0), ("0.4_1.0_0.4", 0)}
SLIPPED = {("0.1_0.3_0.2", 60), ("0.1_0.3_0.4", 60), ("0.4_0.3_0.2", 60), ("0.4_0.3_0.4", 30), ("0.4_0.3_0.4", 60),
           ("1.0_0.3_0.2", 30), ("1.0_0.3_0.2", 60), ("1.0_0.3_0.4", 30), ("1.0_0.3_0.4", 60), ("1.0_0.6_0.2", 30),
           ("1.0_0.6_0.2", 60), ("2.0_0.3_0.2", 30), ("2.0_0.6_0.2", 30)}
DEFECTS = BAD | SLIPPED
