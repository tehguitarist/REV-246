# Plugin character settings against each IR

Null (dB, re the IR wet; lower = closer) per band after the best gain, against the IR wet, whose 2-number clock fit is redone for each variant. See `tools/ab_variants.py`.

| example | 16-bit | XFMR | 32-100 Hz | 100 Hz-4 kHz | 4-12 kHz |
|---|---|---|---|---|---|
| room-mic | off | off | -15.0 | -25.0 | -20.4 |
| room-mic | off | on | -11.4 | -21.0 | -20.2 |
| room-mic | on | off | -14.8 | -24.2 | -20.7 |
| room-mic | on | on | -11.4 | -20.6 | -20.4 |
| piano | off | off | -11.7 | -15.6 | -8.7 |
| piano | off | on | -8.8 | -15.2 | -8.7 |
| piano | on | off | -11.8 | -15.5 | -8.6 |
| piano | on | on | -8.8 | -15.2 | -8.6 |
| acoustic-guitar | off | off | -42.2 | -43.0 | -40.0 |
| acoustic-guitar | off | on | -17.6 | -25.2 | -23.3 |
| acoustic-guitar | on | off | -27.0 | -33.6 | -23.7 |
| acoustic-guitar | on | on | -16.9 | -24.5 | -24.2 |
| snare-ir2 | off | off | -8.3 | -20.4 | -12.3 |
| snare-ir2 | off | on | -8.1 | -20.5 | -12.1 |
| snare-ir2 | on | off | -7.8 | -19.8 | -11.7 |
| snare-ir2 | on | on | -7.7 | -20.0 | -11.5 |
| snare-ir3 | off | off | -36.7 | -34.1 | -17.8 |
| snare-ir3 | off | on | -15.4 | -26.6 | -17.7 |
| snare-ir3 | on | off | -24.7 | -27.7 | -15.1 |
| snare-ir3 | on | on | -16.1 | -24.9 | -15.2 |
