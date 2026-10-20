"""Markdown report and plots for a harness run (stages A–C)."""
import os
from collections import Counter

import numpy as np
import pandas as pd

FLAG_HELP = {
    "error": "analysis raised an exception (see files.csv 'error')",
    "mono": "single-channel file (analysed as dual mono)",
    "clipped": "samples at or above -0.01 dBFS",
    "dc_offset": "DC above -60 dBFS",
    "zero_tail": "ends in digital silence (gated or trimmed tail)",
    "fade_out": "tail falls faster than the reverb decay: a fade was applied (biases RT)",
    "truncated_before_noise_floor": "file ends while the reverb is still decaying",
    "short_file": "shorter than onset + RT",
    "pre_onset_energy": "energy before the onset > noise + 6 dB (sweep-deconvolution artefacts?)",
    "no_pre_roll": "less than 2 ms before the onset (pre-onset check impossible; possible trimmed start)",
    "off_device_grid": "cluster taps do not land on the 26,373.6 Hz grid (different clock or variant?)",
    "tail_too_short_for_loops": "not enough decay above the noise floor to identify loops",
    "label_group_inconsistent": "files with identical labels differ in content",
    "duplicate_with_different_label": "same content as a file with different labels",
    "resampled_copy": "same capture as a same-label file at another sample rate",
    "duplicate_file": "same capture as a same-label file at the same rate",
    "has_repeat_capture": "an independent capture of the same setting exists (useful for E1)",
    "unlabelled_duplicate_of_labelled_file": "no labels, but same content as a labelled file (labels can be inferred)",
    "other_rate_version": "same setting exists at another sample rate (resampled copy or separate capture; not separable by residual)",
    "clock_outlier": "capture clock differs from the library median by > 30 ppm (different interface/session, or resampled)",
    "lr_offset_outlier": "L/R output offset differs from the library median by > 2 us",
}


def _fmt(v, nd=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "–"
    return f"{v:.{nd}g}" if isinstance(v, float) else str(v)


def _table(rows, cols):
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        out.append("| " + " | ".join(_fmt(r.get(c)) for c in cols) + " |")
    return "\n".join(out)


def plots(out, df, b3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    files = []
    prog = df["label_program"].where(df["label_program"].notna(), df["program_measured"]).fillna("?")
    pairs = [("rt", "rt_mid", "RT (s)"), ("lf", "lf_ratio", "LF x"), ("hf", "hf_ratio", "HF x"),
             ("amp", "cluster_gain_db", "AMPL (%)")]
    if any(df[f"label_{p}"].notna().any() for p, _, _ in pairs):
        fig, ax = plt.subplots(2, 2, figsize=(11, 9))
        for a, (p, col, lab) in zip(ax.flat, pairs):
            for k, g in df.groupby(prog):
                g = g.dropna(subset=[f"label_{p}", col])
                if len(g):
                    a.plot(g[f"label_{p}"], g[col], "o", ms=3, alpha=.6, label=k)
            a.set(xlabel=f"label {lab}", ylabel=f"measured {col}", title=lab)
            if p == "rt":
                a.set(xscale="log", yscale="log")
                lim = [0.1, 25]; a.plot(lim, lim, "k:", lw=.8)
            if p in ("lf", "hf"):
                a.plot([0, 2], [0, 2], "k:", lw=.8)
            a.grid(alpha=.3); a.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(out, "labels_vs_measured.png"), dpi=100); plt.close(fig)
        files.append("labels_vs_measured.png")
    c = df[df["has_cluster"] == True]  # noqa: E712
    if len(c):
        fig, ax = plt.subplots(1, 3, figsize=(14, 3.6))
        for fs, g in c.groupby("fs"):
            ax[0].hist(g["clock_ppm"].dropna(), bins=30, alpha=.6, label=f"{int(fs)} Hz")
        ax[0].set(xlabel="capture vs device clock (ppm)", title="B1 clock"); ax[0].legend(fontsize=7)
        ax[1].hist(c["lr_offset_us"].dropna(), bins=30); ax[1].set(xlabel="L lag (us)", title="B4 L/R offset")
        if b3.get("status") == "ok":
            g = c.dropna(subset=["label_delay"])
            ax[2].plot(g["label_delay"], g["tap0_ms"], "o", ms=3)
            xx = np.array([g["label_delay"].min(), g["label_delay"].max()])
            ax[2].plot(xx, b3["slope"] * xx + b3["intercept_ms"], "r-", lw=.8)
        ax[2].set(xlabel="label delay (ms)", ylabel="first tap (ms)", title="B3 latency regression")
        fig.tight_layout(); fig.savefig(os.path.join(out, "calibration.png"), dpi=100); plt.close(fig)
        files.append("calibration.png")
    return files


def write(out, df, pairs, dis, cov, s, b3):
    pngs = plots(out, df, b3)
    L = [f"# Harness report: {s['name']}", "",
         f"{s['n_files']} files from `{s['library']}`. Per-file detail: `files.csv`. Pairs: `pairs.csv`. "
         f"Label audit: `label_disagreements.csv`. Machine-readable summary: `library.json`.", ""]

    # ---------------- A
    L += ["## A. Ingest and QC", ""]
    L.append(f"Sample rates: {dict(Counter(int(v) for v in df['fs'].dropna()))}. "
             f"Errors: {sum('error' in f for f in df['flags'])}.")
    cnt = Counter(f for fl in df["flags"] for f in fl)
    if cnt:
        L += ["", "| flag | files | meaning |", "|---|---|---|"]
        for f, n in cnt.most_common():
            key = f if f in FLAG_HELP else f.split(":")[0]
            L.append(f"| `{f}` | {n} | {FLAG_HELP.get(key, 'label audit, see section C4' if f.startswith('label_mismatch') else '')} |")
    else:
        L.append("\nNo per-file flags.")
    st = df["cluster_status"].value_counts().to_dict()
    L += ["", f"**A4 cluster fingerprint:** {st}. `swapped`/`inverted_*` files must be fixed before any stereo analysis; "
               "`novel_cluster` files have a cluster pattern not in the templates (new program or firmware)."]
    if len(pairs):
        pc = pairs.groupby(["same_label", "class"]).size().to_dict()
        L += ["", f"**A5/A7 pair comparison:** {len(pairs)} pairs compared. Counts by (same label, class): {pc}. "
                  "`copy` = same capture; `repeat` = separate capture of the same setting (residual at noise level); "
                  "`repeat_or_resampled` = as `repeat` but across sample rates; "
                  "`near` = same setting with extra differences (time variance? drift?); `different` = different settings."]
        rep = pairs[pairs["class"].isin(["repeat", "near"])]
        if len(rep):
            L.append(f"Repeat/near pairs, median residual L {rep['residual_L_db'].median():.1f} dB, "
                     f"R {rep['residual_R_db'].median():.1f} dB. If L is consistently worse than R, that supports eugtone's "
                     "L-only modulation report (spec E1).")
    n = s["normalisation"]
    L += ["", f"**A6 normalisation:** {n['frac_files_at_peak_mode']:.0%} of files peak at {_fmt(n['peak_mode_dbfs'])} dBFS "
              f"(peak-normalised likely: {n['peak_normalised_likely']}). Cluster-gain spread within same-label groups: "
              f"{_fmt(n['cluster_gain_spread_db_median'])} dB. Amplitude will be measured as **{n['amp_measure']}**."]

    # ---------------- B
    L += ["", "## B. Calibration", ""]
    for fs, c in s["clock"].items():
        L.append(f"- **B1 clock @ {fs} Hz:** median {c['median_ppm']:+.1f} ppm (IQR {c['iqr_ppm']:.1f}), n={c['n']}, "
                 f"off-grid files: {c['off_grid_files']}.")
    for fs, k in s["kernels"].items():
        L.append(f"- **B2 kernel @ {fs} Hz:** from {k['n_files']} files; gain vs reference {k['gain_vs_ref_db']:+.2f} dB, "
                 f"max magnitude deviation {k['mag_dev_vs_ref_db_max']:.2f} dB (100 Hz–11 kHz), excess phase "
                 f"{k['excess_phase_vs_ref_rad_max']:.3f} rad, within-library spread {_fmt(k['within_library_spread_db_max'])} dB.")
    lat = s["latency"]
    if lat.get("status") == "ok":
        L.append(f"- **B3 latency:** first tap = {lat['slope']:.4f} x delay + {lat['intercept_ms']:.3f} ms "
                 f"(rms residual {lat['rms_residual_ms']:.3f} ms, n={lat['n']}, outliers {len(lat['outliers'])}).")
    else:
        L.append(f"- **B3 latency:** {lat.get('status')} (n={lat.get('n')}, distinct delays={lat.get('distinct_delays')}).")
    lr = s["lr_offset"]
    if lr.get("n"):
        L.append(f"- **B4 L/R offset:** median {lr['median_us']:.2f} us (MAD {lr['mad_us']:.2f}), outliers {len(lr['outliers'])}.")
    fp = s["unit_fingerprint"]
    L.append(f"- **B5 unit fingerprint:** loop sets {fp['loopset_counts']}; L-R noise {_fmt(fp['noise_L_minus_R_db'])} dB "
             f"(reference {fp['reference'].get('noise_L_minus_R_db')}); L/R offset vs reference {_fmt(fp['lr_offset_diff_us'])} us.")
    for note in fp["notes"]:
        L.append(f"  - {note}")
    if "calibration.png" in pngs:
        L += ["", "![calibration](calibration.png)"]

    # ---------------- C
    L += ["", "## C. Settings inference", ""]
    lbl = df[[c for c in df.columns if c.startswith("label_") and c != "label_patterns"]].notna().mean()
    L.append("Label coverage (fraction of files with each label): " + ", ".join(f"{k[6:]} {v:.0%}" for k, v in lbl.items()) + ".")
    unparsed = df["label_patterns"].eq("").sum()
    if unparsed:
        L.append(f"**{unparsed} files matched no naming pattern.** Run with `--templates` and add a pattern to config.json.")
    L += ["", "### C3 step tables", ""]
    for p, tables in s["step_tables"].items():
        L.append(f"**Program {p}**")
        for param, t in tables.items():
            line = f"- {param} (`{t['measured_col']}`): {t['n_clusters']} measured clusters (expected {t['expected_steps']})"
            if "labelled" in t:
                line += f"; {len(t['labelled'])} labelled values, monotonic: {t['monotonic']}"
            L.append(line)
            if "labelled" in t and param in ("rt", "lf", "hf", "amp", "lp"):
                L += ["", _table(t["labelled"], ["label", "count", "median", "mad", "min", "max"]), ""]
    L += ["", f"### C4 label audit: {len(dis)} disagreements", ""]
    if len(dis):
        L.append(_table([{"param": k, "files": v} for k, v in dis.groupby("param").size().items()], ["param", "files"]))
        L += ["", _table(dis.head(25).to_dict("records"), ["relpath", "param", "label", "measured", "deviation"])]
    if "labels_vs_measured.png" in pngs:
        L += ["", "![labels vs measured](labels_vs_measured.png)"]
    L += ["", "### C5 coverage", "", f"Files per program: {cov['files_per_program']}", "",
          "Best available sweep for each planned set: the largest group of files where only the swept parameter changes.", "",
          "| planned set | program | have | need | swept values | held fixed |", "|---|---|---|---|---|---|"]
    for name, res in cov["planned"].items():
        for p, r in res.items():
            mark = "" if r["have"] >= r["need"] else " ⚠"
            held = ", ".join(f"{k}={v}" for k, v in r.get("held", {}).items())
            vals = ", ".join(f"{v:g}" for v in r["values"][:16]) + (" …" if len(r["values"]) > 16 else "")
            L.append(f"| {name} | {p} | {r['have']}{mark} | {r['need']} | {vals} | {held} |")
    L += ["", "RT × other-parameter count matrices: `coverage_<program>_rt_x_<param>.csv`."]
    with open(os.path.join(out, "report.md"), "w") as f:
        f.write("\n".join(L) + "\n")
