"""EMT 246 IR-library harness, stages A–C (ingest/QC, calibration, settings inference, coverage).

Usage:
  .venv/bin/python tools/harness/run.py <library_dir> --name tomkin            # full run
  .venv/bin/python tools/harness/run.py <library_dir> --templates               # show naming scheme
  .venv/bin/python tools/harness/run.py <library_dir> --name t --limit 50       # quick trial
  .venv/bin/python tools/harness/run.py "irs/EMT 246 captures/96kHz" --name tomkin96   # the Tomkin working set
  .venv/bin/python tools/harness/run.py "irs/EMT 246 captures" --name trial --sample 60   # 60 settings x 2 rates

Outputs go to analysis/library/<name>/ (see report.md there). Per-file results are cached, so a
re-run after editing the naming patterns or library-stage code only redoes the cheap parts.
"""
import argparse
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import library as lib  # noqa: E402
import report  # noqa: E402
from core import ROOT, Kernel, analyze_file, file_key, load_config  # noqa: E402
from naming import parse_name, name_templates  # noqa: E402

EXTS = (".wav", ".aif", ".aiff", ".flac")


def discover(root):
    out = []
    for d, _, files in os.walk(root):
        for f in files:
            if f.lower().endswith(EXTS) and not f.startswith("._"):
                p = os.path.join(d, f)
                out.append((p, os.path.relpath(p, root)))
    return sorted(out, key=lambda t: t[1])


def sample_settings(files, n):
    """Keep every file of n settings spread evenly over the sorted settings (a setting = the path with its
    sample-rate folder removed, so both rates of a setting are kept together for the cross-rate tests)."""
    key = lambda rel: re.sub(r"(^|[/\\\\])\d+\s*k(hz)?(?=[/\\\\])", r"\1", rel, flags=re.IGNORECASE).lstrip("/\\")
    keys = sorted({key(r) for _, r in files})
    if n >= len(keys):
        return files
    pick = {keys[int(i)] for i in np.linspace(0, len(keys) - 1, n)}
    return [f for f in files if key(f[1]) in pick]


def _work(args):
    path, cfg, kernel, lr, cache = args
    key = file_key(path, kernel, lr)
    cj, cn = os.path.join(cache, key + ".json"), os.path.join(cache, key + ".npy")
    if os.path.exists(cj):
        with open(cj) as f:
            row = json.load(f)
        return row, (np.load(cn) if os.path.exists(cn) else None)
    try:
        row, sig = analyze_file(path, cfg, kernel, lr)
    except Exception as e:  # keep going; the file is reported as an error
        row, sig = {"path": path, "flags": ["error"], "error": repr(e)}, None
    with open(cj, "w") as f:
        json.dump(row, f, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    if sig is not None:
        np.save(cn, sig)
    return row, sig


def run_pass(jobs, workers, label):
    t = time.time()
    rows, sigs = [], []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (r, s) in enumerate(ex.map(_work, jobs, chunksize=4)):
            rows.append(r); sigs.append(s)
            if (i + 1) % 100 == 0 or i + 1 == len(jobs):
                print(f"  {label}: {i + 1}/{len(jobs)} files ({time.time() - t:.0f} s)", flush=True)
    return rows, sigs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("library")
    ap.add_argument("--name", default=None)
    ap.add_argument("--config", default=None)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-refine", action="store_true", help="skip pass 2 with library-derived kernels")
    ap.add_argument("--separate-rates", action="store_true",
                    help="same-label files at different sample rates are separate captures (sets config separate_rate_captures); only for libraries whose rates really are separate recordings, NOT the Tomkin library")
    ap.add_argument("--sample", type=int, default=None, metavar="N",
                    help="only N settings, spread evenly, with every sample rate of each (use for quick cross-rate trials)")
    ap.add_argument("--templates", action="store_true", help="print filename templates and exit")
    a = ap.parse_args()

    cfg = load_config(a.config)
    if a.separate_rates:
        cfg["separate_rate_captures"] = True
    files = discover(a.library)
    if a.sample:
        files = sample_settings(files, a.sample)
    files = files[: a.limit]
    if a.templates:
        for t, n, ex in name_templates([r for _, r in files]):
            print(f"{n:6d}  {t}\n        e.g. {ex}")
        return
    if not files:
        sys.exit(f"no audio files under {a.library}")

    name = a.name or os.path.basename(os.path.normpath(a.library))
    out = os.path.join(ROOT, "analysis", "library", name)
    cache = os.path.join(out, "cache")
    os.makedirs(cache, exist_ok=True)
    os.makedirs(os.path.join(out, "kernels"), exist_ok=True)
    print(f"{len(files)} files -> {out}")

    ref = Kernel.load(cfg["reference_kernel"])
    lr0 = cfg["default_lr_offset_us"]
    rows, sigs = run_pass([(p, cfg, ref, lr0, cache) for p, _ in files], a.workers, "pass 1 (reference kernel)")
    df = _frame(rows, files, cfg)

    # B2: library kernels per sample rate, then pass 2 with them
    kernels, kcmp = {}, {}
    for fs in sorted(df["fs"].dropna().unique()):
        k, cmp = lib.build_kernel(df, fs, ref)
        if k is not None:
            kernels[fs] = k
            kcmp[int(fs)] = cmp
            k.save(os.path.join(out, "kernels", f"kernel_{int(fs)}.npz"))
    lr_lib = lib.lr_offsets(df)
    lr1 = round(lr_lib["median_us"], 2) if lr_lib.get("n") else lr0
    if kernels and not a.no_refine:
        idx = [i for i, (p, _) in enumerate(files) if df.at[i, "fs"] in kernels]
        r2, s2 = run_pass([(files[i][0], cfg, kernels[df.at[i, "fs"]], lr1, cache) for i in idx], a.workers,
                          "pass 2 (library kernel)")
        for j, i in enumerate(idx):
            rows[i], sigs[i] = r2[j], s2[j]
        df = _frame(rows, files, cfg)
    sig = dict(enumerate(sigs))

    # A5–A7
    pairs = lib.compare_pairs(df, sig, cfg.get("separate_rate_captures", False))
    pf = lib.pair_flags(df, pairs) if len(pairs) else {}
    for i, extra in pf.items():
        df.at[i, "flags"] = sorted(set(df.at[i, "flags"]) | set(extra))
    norm = lib.normalisation(df)
    # B1, B3–B5
    clk = lib.clock(df)
    for i, r in df[df["clock_ppm"].notna()].iterrows():
        med = clk.get(int(r["fs"]), {}).get("median_ppm")
        if med is not None and abs(r["clock_ppm"] - med) > 30:
            df.at[i, "flags"] = df.at[i, "flags"] + ["clock_outlier"]
    lr = lib.lr_offsets(df)
    for i in lr.get("outliers", []):
        df.at[i, "flags"] = df.at[i, "flags"] + ["lr_offset_outlier"]
    b3 = lib.latency_regression(df)
    fp = lib.unit_fingerprint(df, cfg, lr, clk)
    # C3–C5
    steps = lib.step_tables(df, cfg, norm["amp_measure"])
    dis = lib.disagreements(df, cfg, norm["amp_measure"], b3)
    for _, d in dis.iterrows():
        i = df.index[df["relpath"] == d["relpath"]][0]
        df.at[i, "flags"] = df.at[i, "flags"] + [f"label_mismatch_{d['param']}"]
    cov = lib.coverage(df, cfg)

    summary = {"name": name, "library": os.path.abspath(a.library), "n_files": len(df), "kernels": kcmp,
               "lr_offset_used_in_pass2_us": lr1, "normalisation": norm, "clock": clk, "lr_offset": lr,
               "latency": {k: v for k, v in b3.items() if k != "residual_ms"}, "unit_fingerprint": fp,
               "step_tables": steps, "coverage": {k: v for k, v in cov.items() if k != "matrices"}}
    _write(out, df, pairs, dis, cov, summary)
    report.write(out, df, pairs, dis, cov, summary, b3)
    print(f"done: {os.path.join(out, 'report.md')}")


def _frame(rows, files, cfg):
    df = pd.DataFrame(rows)
    df["relpath"] = [r for _, r in files]
    lab = pd.DataFrame([parse_name(r, cfg["naming"]) for _, r in files])
    df = pd.concat([df.drop(columns=[c for c in lab.columns if c in df.columns]), lab], axis=1)
    for c in ("has_cluster", "fs", "cluster_status", "program_measured", "tap0_ms", "lr_offset_us", "clock_ppm",
              "grid_max_frac", "cluster_gain_db", "cluster_to_reverb_db", "cluster_sparsity", "rt_mid",
              "lf_ratio", "hf_ratio", "lp_corner_khz", "noise_rel_db", "noise_dbfs_L", "noise_dbfs_R",
              "peak_dbfs_L", "peak_dbfs_R"):
        if c not in df:
            df[c] = None
    df["flags"] = df["flags"].apply(lambda f: list(f) if isinstance(f, (list, tuple)) else [])
    return df


def _write(out, df, pairs, dis, cov, summary):
    flat = df.copy()
    for c in flat.columns:
        if flat[c].apply(lambda v: isinstance(v, (list, tuple))).any():
            flat[c] = flat[c].apply(lambda v: ";".join(map(str, v)) if isinstance(v, (list, tuple)) else v)
    flat.to_csv(os.path.join(out, "files.csv"), index=False)
    pairs.to_csv(os.path.join(out, "pairs.csv"), index=False)
    dis.to_csv(os.path.join(out, "label_disagreements.csv"), index=False)
    for p, mats in cov["matrices"].items():
        for other, m in mats.items():
            m.to_csv(os.path.join(out, f"coverage_{p}_rt_x_{other}.csv"))
    with open(os.path.join(out, "library.json"), "w") as f:
        json.dump(summary, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


if __name__ == "__main__":
    main()
