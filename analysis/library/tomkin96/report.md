# Harness report: tomkin96

2437 files from `/Users/leighpierce/Git/246/irs/EMT 246 captures/96kHz`. Per-file detail: `files.csv`. Pairs: `pairs.csv`. Label audit: `label_disagreements.csv`. Machine-readable summary: `library.json`.

## A. Ingest and QC

Sample rates: {96000: 2437}. Errors: 0.

| flag | files | meaning |
|---|---|---|
| `short_file` | 2003 | shorter than onset + RT |
| `fade_out` | 1751 | tail falls faster than the reverb decay: a fade was applied (biases RT) |
| `zero_tail` | 1622 | ends in digital silence (gated or trimmed tail) |
| `no_pre_roll` | 1323 | less than 2 ms before the onset (pre-onset check impossible; possible trimmed start) |
| `label_mismatch_rt` | 1190 | label audit, see section C4 |
| `pre_onset_energy` | 1091 | energy before the onset > noise + 6 dB (sweep-deconvolution artefacts?) |
| `label_mismatch_lf` | 907 | label audit, see section C4 |
| `label_mismatch_hf` | 843 | label audit, see section C4 |
| `truncated_before_noise_floor` | 632 | file ends while the reverb is still decaying |
| `duplicate_with_different_label` | 285 | same content as a file with different labels |
| `label_mismatch_delay` | 195 | label audit, see section C4 |
| `tail_too_short_for_loops` | 73 | not enough decay above the noise floor to identify loops |
| `clock_outlier` | 15 | capture clock differs from the library median by > 30 ppm (different interface/session, or resampled) |

**A4 cluster fingerprint:** {'no_cluster': 2038, 'ok': 384, 'novel_cluster': 15}. `swapped`/`inverted_*` files must be fixed before any stereo analysis; `novel_cluster` files have a cluster pattern not in the templates (new program or firmware).

**A5/A7 pair comparison:** 12094 pairs compared. Counts by (same label, class): {(False, 'different'): 10044, (False, 'near'): 1431, (False, 'repeat'): 619}. `copy` = same capture; `repeat` = separate capture of the same setting (residual at noise level); `repeat_or_resampled` = as `repeat` but across sample rates; `near` = same setting with extra differences (time variance? drift?); `different` = different settings.
Repeat/near pairs, median residual L -31.8 dB, R -32.0 dB. If L is consistently worse than R, that supports eugtone's L-only modulation report (spec E1).

**A6 normalisation:** 7% of files peak at -26.4 dBFS (peak-normalised likely: False). Cluster-gain spread within same-label groups: 0.0511 dB. Amplitude will be measured as **cluster_gain_db**.

## B. Calibration

- **B1 clock @ 96000 Hz:** median -13.2 ppm (IQR 4.8), n=399, off-grid files: 0.
- **B2 kernel @ 96000 Hz:** from 8 files; gain vs reference -22.93 dB, max magnitude deviation 0.64 dB (100 Hz–11 kHz), excess phase 0.202 rad, within-library spread 0.00847 dB.
- **B3 latency:** not_enough_delay_values (n=399, distinct delays=1).
- **B4 L/R offset:** median 19.19 us (MAD 0.00), outliers 0.
- **B5 unit fingerprint:** loop sets {'?': 1439, 'U:1108-1254-1386': 310, 'S': 101, 'U:2641': 37, 'U:231-628-803-1108-1292-2032': 37, 'U:628': 35, 'U:628-1108-1386': 28, 'U:5279': 26, 'U:1321-1389-1461': 26, 'U:1108': 23, 'U:1321-1389': 23, 'U:628-1108': 18, 'U:1108-1386': 16, 'U:5279-10826-12176': 15, 'U:5279-10826-12176-16105-16663-17455-17957': 15, 'U:924-1108-1254-1386': 15, 'U:231-1108-1292-2032': 14, 'U:231-803-1108-1292-2032': 14, 'U:5280': 11, 'U:231-1108-1292': 11, 'U:231-248-628-803-1108': 11, 'U:5279-10826-11620-12176-16105-16663-17455-17957': 10, 'U:1108-1254': 10, 'U:1321': 10, 'U:5279-10826': 8, 'U:2641-5416': 7, 'U:2641-5416-6090': 7, 'U:231-628-803-1108-1292': 7, 'U:231-628-1108': 6, 'U:2641-5416-11512': 5, 'U:2641-5416-5814-6090-8057-8335-8731-8983-11508-11624-14677-17040-17318-17714': 5, 'U:1108-1292': 5, 'U:231-248-310-318-438-446-524-628': 5, 'U:231-248-310-480-628-803-1108': 5, 'U:924-1108-1254-1292-1386': 4, 'U:231-248-340-628-803-1108': 4, 'U:231-628-803-1108': 4, 'U:248-628': 4, 'U:231-1108-1254-1292-2032': 4, 'U:5279-10826-12176-16105-17455': 3, 'U:5279-10827': 3, 'U:5279-10826-12176-16105-16663-17455': 3, 'U:5279-10827-12177': 3, 'U:920-2641': 3, 'U:231-248-480-628-803-1108': 3, 'U:248-438-446-628': 3, 'U:1389': 3, 'U:5279-10826-11625-12176': 2, 'U:2642': 2, 'U:2641-5416-6090-11512': 2, 'U:2641-5416-5814-6090-8057-8335-8731-11508-11624-14677-17040-17318': 2, 'U:231-248-310-340-480-628-803': 2, 'U:231-803-1108-1292': 2, 'U:1321-1389-2922': 2, 'U:5279-10826-12176-16106': 1, 'U:5279-10826-12176-16106-17455': 1, 'U:1836-5280': 1, 'U:5279-10826-11623-12176-16105-17455': 1, 'U:5279-10826-12176-17455': 1, 'U:5279-10826-12176-16106-17456': 1, 'U:5279-10826-12176-16106-16663-17455': 1, 'U:5279-10826-11623-12176': 1, 'U:5279-6341-10826-11620-12176-16105-16663-17455-17957': 1, 'U:920': 1, 'U:2641-5416-6090-8057-8335-11510': 1, 'U:2641-5416-6090-8057-8335-11510-11624': 1, 'U:2641-5416-5814-6090-8057-8335-8731-8983-11508-11624-14677-17040-17318': 1, 'U:2641-5416-6090-8057-8335-11510-17318': 1, 'U:2641-5416-6090-8057-8335-8731-8983-11510-11624-14677-17040-17318-17714': 1, 'U:2641-5416-6090-14677-17040-17318': 1, 'U:2641-5416-6090-8057-11510': 1, 'U:2641-5416-6090-8057-8335-11508-17040-17318': 1, 'U:2641-5416-5814-6090-8057-8335-8731-11508-11624-14677-17040-17318-17714': 1, 'U:2641-5416-6090-8057-8335-11510-14677-17040-17318': 1, 'U:2641-5416-6090-8057-8335-11510-11624-14677-17040-17318': 1, 'U:2641-5416-6090-11510': 1, 'U:2641-11512': 1, 'U:2641-5416-6090-8057-8335-11510-17040-17318': 1, 'U:2641-5416-6090-8057-8335-8731-8983-11508-11624-14677-17040-17318-17714': 1, 'U:2641-5416-6090-8057-8335-8731-11510-11624-14677-17040-17318-17714': 1, 'U:231-628-1108-1292-2032': 1, 'U:231-248-340-480-628-803-1108': 1, 'U:1108-1254-1292-2032': 1, 'U:231-248-256-310-318-438-446-524-628': 1, 'U:231-340-628-803-1108': 1, 'U:231-248-310-316-318-438-446-524-628': 1, 'U:231-256-628': 1, 'U:231-583-1108-1292-2032': 1, 'U:231-1108-1254-1292': 1, 'U:1108-1254-1292': 1, 'U:1017-2103': 1, 'U:1108-1254-1292-1386': 1, 'U:231-248-310-480-628': 1, 'U:1017-1108-1254-1386': 1, 'U:1108-1292-2032': 1, 'U:231-803-1108': 1, 'U:803-1108': 1, 'U:248-628-1108': 1, 'U:231-248-310-340-480-628-803-1108': 1, 'U:810-1108-1386': 1, 'U:810-1108': 1, 'U:810': 1, 'U:810-1108-1254-1386': 1, 'U:231-248-340-480-628': 1, 'U:534-1137-1428': 1, 'U:534-603-1108-1386-1428': 1, 'U:534-1108': 1, 'U:534-1108-1137-1254-1386-1428': 1, 'U:534-1108-1254-1386': 1, 'U:231-628-1108-1292': 1, 'U:248-438-446-524-628': 1, 'U:231-1108': 1, 'U:231-248-310-628': 1, 'U:924-1108-1254-1292-1386-2032': 1, 'U:231-248-340-480-628-803': 1, 'U:1321-1389-3174': 1, 'U:1321-2922': 1, 'U:1321-1389-1461-3174': 1}; L-R noise 0.423 dB (reference 6.6); L/R offset vs reference 0.589 us.
  - 116 loop set(s) not in config.known_loopsets (expected for M/L until added; a set close to S would suggest different firmware)
  - some cluster patterns differ from the reference templates (other programs, or different firmware)

![calibration](calibration.png)

## C. Settings inference

Label coverage (fraction of files with each label): program 100%, rt 57%, lf 57%, hf 57%, delay 100%, amp 57%, lp 100%, span 43%, early 43%, late 43%.

### C3 step tables

**Program L**
- rt (`rt_mid`): 21 measured clusters (expected 16); 9 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 2.2 | 40 | 2.06 | 0.154 | 1.67 | 2.3 |
| 2.6 | 40 | 2.38 | 0.256 | 1.3 | 2.71 |
| 3.6 | 40 | 3.36 | 0.22 | 1.55 | 3.64 |
| 4.5 | 40 | 4.1 | 0.361 | 1.62 | 4.62 |
| 6 | 40 | 5.42 | 0.497 | 1.64 | 6.16 |
| 8 | 40 | 7.27 | 0.497 | 2.76 | 8.36 |
| 10 | 40 | 9.21 | 0.669 | 4.17 | 10.4 |
| 15 | 40 | 14 | 0.958 | 7.22 | 15.6 |
| 20 | 40 | 18.3 | 1.6 | 9.89 | 20.8 |

- lf (`lf_ratio`): 9 measured clusters (expected 16); 4 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.1 | 90 | 0.253 | 0.0956 | 0.15 | 0.722 |
| 0.8 | 90 | 0.874 | 0.0558 | 0.781 | 1.62 |
| 1.2 | 90 | 1.17 | 0.0524 | 1.06 | 2.77 |
| 2 | 90 | 1.59 | 0.193 | 1.14 | 3.73 |

- hf (`hf_ratio`): 6 measured clusters (expected 10); 5 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.1 | 72 | 0.366 | 0.131 | 0.135 | 0.563 |
| 0.3 | 72 | 0.385 | 0.0483 | 0.3 | 0.565 |
| 0.5 | 72 | 0.543 | 0.0406 | 0.475 | 0.698 |
| 0.8 | 72 | 0.813 | 0.0348 | 0.732 | 1.13 |
| 1 | 72 | 1 | 0.0441 | 0.916 | 1.51 |

- lp (`lp_corner_khz`): 1 measured clusters (expected 8); 1 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 12 | 360 | 11 | 0 | 11 | 11 |

**Program M**
- rt (`rt_mid`): 7 measured clusters (expected 16); 8 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.6 | 65 | 0.553 | 0.0666 | 0.435 | 0.709 |
| 0.8 | 65 | 0.676 | 0.129 | 0.47 | 0.879 |
| 1.2 | 65 | 1.05 | 0.183 | 0.561 | 1.27 |
| 1.8 | 65 | 1.56 | 0.229 | 0.724 | 1.86 |
| 2.6 | 65 | 2.25 | 0.317 | 0.875 | 2.68 |
| 3 | 65 | 2.6 | 0.354 | 0.929 | 3.09 |
| 4 | 64 | 3.55 | 0.367 | 1.46 | 4.13 |
| 5 | 65 | 4.28 | 0.605 | 1.7 | 5.16 |

- lf (`lf_ratio`): 5 measured clusters (expected 16); 5 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.1 | 119 | 0.354 | 0.181 | 0.142 | 0.93 |
| 0.4 | 40 | 0.595 | 0.146 | 0.412 | 1.24 |
| 0.7 | 120 | 0.835 | 0.091 | 0.615 | 1.85 |
| 1.2 | 120 | 1.19 | 0.082 | 0.947 | 3.24 |
| 2 | 120 | 1.84 | 0.258 | 1.35 | 5.31 |

- hf (`hf_ratio`): 1 measured clusters (expected 10); 5 labelled values, monotonic: False

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.1 | 103 | 0.446 | 0.0799 | 0.186 | 0.699 |
| 0.2 | 104 | 0.366 | 0.0991 | 0.226 | 0.702 |
| 0.4 | 104 | 0.494 | 0.0482 | 0.387 | 0.699 |
| 0.7 | 104 | 0.742 | 0.0493 | 0.616 | 0.906 |
| 0.9 | 104 | 0.918 | 0.0435 | 0.79 | 1.18 |

- lp (`lp_corner_khz`): 1 measured clusters (expected 8); 1 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 12 | 519 | 11 | 0 | 11 | 11 |

**Program N**
- rt (`rt_mid`): 15 measured clusters (expected 16)
- lf (`lf_ratio`): 16 measured clusters (expected 16)
- hf (`hf_ratio`): 3 measured clusters (expected 10)
- amp (`cluster_gain_db`): 2 measured clusters (expected 11)
- delay (`tap0_ms`): 10 measured clusters (expected 1000); 1 labelled values, monotonic: True
- lp (`lp_corner_khz`): 1 measured clusters (expected 8); 1 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 12 | 1046 | 11 | 0 | 11 | 11 |

**Program S**
- rt (`rt_mid`): 8 measured clusters (expected 16); 8 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.3 | 64 | 0.281 | 0.0351 | 0.227 | 0.346 |
| 0.6 | 64 | 0.563 | 0.0402 | 0.302 | 0.634 |
| 0.8 | 64 | 0.741 | 0.0487 | 0.431 | 0.822 |
| 1 | 64 | 0.924 | 0.0603 | 0.603 | 1.02 |
| 1.4 | 64 | 1.27 | 0.103 | 0.925 | 1.43 |
| 1.8 | 64 | 1.65 | 0.114 | 1.18 | 1.85 |
| 2.2 | 64 | 2 | 0.149 | 1.44 | 2.27 |
| 2.5 | 64 | 2.26 | 0.149 | 1.64 | 2.55 |

- lf (`lf_ratio`): 5 measured clusters (expected 16); 4 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.1 | 128 | 0.282 | 0.0914 | 0.172 | 1.06 |
| 0.4 | 128 | 0.529 | 0.0608 | 0.441 | 1.08 |
| 1 | 128 | 1.04 | 0.0455 | 0.982 | 1.81 |
| 2 | 128 | 1.49 | 0.142 | 1.21 | 3.39 |

- hf (`hf_ratio`): 2 measured clusters (expected 10); 4 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 0.2 | 128 | 0.325 | 0.0589 | 0.222 | 0.689 |
| 0.4 | 128 | 0.475 | 0.0301 | 0.39 | 0.69 |
| 0.6 | 128 | 0.643 | 0.0428 | 0.553 | 0.77 |
| 1 | 128 | 1.03 | 0.083 | 0.91 | 1.34 |

- amp (`cluster_gain_db`): 3 measured clusters (expected 11); 3 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 30 | 128 | -9.55 | 0.00688 | -9.57 | -9.54 |
| 60 | 128 | -3.51 | 0.0012 | -3.56 | -3.51 |
| 90 | 128 | 0.0536 | 0.00996 | 0.0149 | 0.066 |

- delay (`tap0_ms`): 3 measured clusters (expected 1000); 1 labelled values, monotonic: True
- lp (`lp_corner_khz`): 1 measured clusters (expected 8); 1 labelled values, monotonic: True

| label | count | median | mad | min | max |
|---|---|---|---|---|---|
| 12 | 512 | 11 | 0 | 11 | 11 |


### C4 label audit: 3135 disagreements

| param | files |
|---|---|
| delay | 195 |
| hf | 843 |
| lf | 907 |
| rt | 1190 |

| relpath | param | label | measured | deviation |
|---|---|---|---|---|
| Long RT/12_0_0_01_22_01_dc.wav | rt | 2.2 | 1.68 | -0.382 |
| Long RT/12_0_0_01_22_03_dc.wav | rt | 2.2 | 1.67 | -0.389 |
| Long RT/12_0_0_01_22_05_dc.wav | rt | 2.2 | 1.88 | -0.184 |
| Long RT/12_0_0_01_22_08_dc.wav | rt | 2.2 | 2.13 | 0.0729 |
| Long RT/12_0_0_01_22_10_dc.wav | rt | 2.2 | 2.26 | 0.203 |
| Long RT/12_0_0_08_22_01_dc.wav | rt | 2.2 | 1.8 | -0.26 |
| Long RT/12_0_0_08_22_03_dc.wav | rt | 2.2 | 1.75 | -0.309 |
| Long RT/12_0_0_08_22_05_dc.wav | rt | 2.2 | 1.89 | -0.167 |
| Long RT/12_0_0_08_22_08_dc.wav | rt | 2.2 | 2.18 | 0.114 |
| Long RT/12_0_0_08_22_10_dc.wav | rt | 2.2 | 2.17 | 0.114 |
| Long RT/12_0_0_12_22_01_dc.wav | rt | 2.2 | 1.83 | -0.235 |
| Long RT/12_0_0_12_22_03_dc.wav | rt | 2.2 | 1.9 | -0.157 |
| Long RT/12_0_0_12_22_08_dc.wav | rt | 2.2 | 2.17 | 0.109 |
| Long RT/12_0_0_12_22_10_dc.wav | rt | 2.2 | 2.2 | 0.138 |
| Long RT/12_0_0_20_22_05_dc.wav | rt | 2.2 | 2.18 | 0.117 |
| Long RT/12_0_0_20_22_08_dc.wav | rt | 2.2 | 2.29 | 0.229 |
| Long RT/12_0_0_20_22_10_dc.wav | rt | 2.2 | 2.3 | 0.24 |
| Long RT/12_0_80_01_22_01_dc.wav | rt | 2.2 | 1.9 | -0.163 |
| Long RT/12_0_80_01_22_03_dc.wav | rt | 2.2 | 1.93 | -0.131 |
| Long RT/12_0_80_01_22_05_dc.wav | rt | 2.2 | 2 | -0.061 |
| Long RT/12_0_80_01_22_08_dc.wav | rt | 2.2 | 2.22 | 0.155 |
| Long RT/12_0_80_01_22_10_dc.wav | rt | 2.2 | 2.26 | 0.202 |
| Long RT/12_0_80_08_22_01_dc.wav | rt | 2.2 | 1.91 | -0.153 |
| Long RT/12_0_80_08_22_03_dc.wav | rt | 2.2 | 1.87 | -0.192 |
| Long RT/12_0_80_08_22_05_dc.wav | rt | 2.2 | 1.97 | -0.0922 |

![labels vs measured](labels_vs_measured.png)

### C5 coverage

Files per program: {'N': 1046, 'M': 519, 'S': 512, 'L': 360}

Best available sweep for each planned set: the largest group of files where only the swept parameter changes.

| planned set | program | have | need | swept values | held fixed |
|---|---|---|---|---|---|
| core_rt | S | 8 ⚠ | 16 | 0.3, 0.6, 0.8, 1, 1.4, 1.8, 2.2, 2.5 | lf=0.1, hf=0.2, amp=0.0, delay=0.0, lp=12.0 |
| core_rt | M | 8 ⚠ | 16 | 0.6, 0.8, 1.2, 1.8, 2.6, 3, 4, 5 | lf=0.1, hf=0.1, amp=60.0, delay=0.0, lp=12.0 |
| core_rt | L | 9 ⚠ | 16 | 2.2, 2.6, 3.6, 4.5, 6, 8, 10, 15, 20 | lf=0.1, hf=0.1, amp=0.0, delay=0.0, lp=12.0 |
| lf_law | S | 4 ⚠ | 16 | 0.1, 0.4, 1, 2 | rt=0.3, hf=0.2, amp=0.0, delay=0.0, lp=12.0 |
| lf_law | M | 5 ⚠ | 16 | 0.1, 0.4, 0.7, 1.2, 2 | rt=0.6, hf=0.1, amp=10.0, delay=0.0, lp=12.0 |
| lf_law | L | 4 ⚠ | 16 | 0.1, 0.8, 1.2, 2 | rt=10.0, hf=0.1, amp=0.0, delay=0.0, lp=12.0 |
| hf_law | S | 4 ⚠ | 10 | 0.2, 0.4, 0.6, 1 | rt=0.3, lf=0.1, amp=0.0, delay=0.0, lp=12.0 |
| hf_law | M | 5 ⚠ | 10 | 0.1, 0.2, 0.4, 0.7, 0.9 | rt=0.6, lf=0.1, amp=10.0, delay=0.0, lp=12.0 |
| hf_law | L | 5 ⚠ | 10 | 0.1, 0.3, 0.5, 0.8, 1 | rt=10.0, lf=0.1, amp=0.0, delay=0.0, lp=12.0 |
| cluster | S | 4 | 3 | 0, 30, 60, 90 | rt=0.3, lf=0.1, hf=0.2, delay=0.0, lp=12.0 |
| cluster | M | 3 | 3 | 10, 60, 99 | rt=0.6, lf=0.1, hf=0.1, delay=0.0, lp=12.0 |
| cluster | L | 2 ⚠ | 3 | 0, 80 | rt=10.0, lf=0.1, hf=0.1, delay=0.0, lp=12.0 |
| lowpass | S | 1 ⚠ | 8 | 12 | rt=0.3, lf=0.1, hf=0.2, amp=0.0, delay=0.0 |
| predelay | S | 1 ⚠ | 3 | 0 | rt=0.3, lf=0.1, hf=0.2, amp=0.0, lp=12.0 |

RT × other-parameter count matrices: `coverage_<program>_rt_x_<param>.csv`.
