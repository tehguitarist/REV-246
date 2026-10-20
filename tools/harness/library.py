"""Library-level stages of the harness: A5–A7 (duplicates, normalisation, multi-rate copies),
B2–B5 (kernels, latency, L/R offset, unit fingerprint) and C3–C5 (step tables, label audit, coverage).

All functions take the per-file table (pandas DataFrame built from core.analyze_file rows).
"""
import itertools

import numpy as np
import pandas as pd
import soundfile as sf

from core import Kernel, fs_dev

LABELS = ["label_program", "label_rt", "label_lf", "label_hf", "label_delay", "label_amp", "label_lp",
          "label_span", "label_early", "label_late"]
MEASURED = {"rt": "rt_mid", "lf": "lf_ratio", "hf": "hf_ratio", "amp": "cluster_gain_db",
            "delay": "tap0_ms", "lp": "lp_corner_khz"}


def label_key(df):
    """One string per file joining all label fields (empty for missing)."""
    return df[LABELS].apply(lambda col: col.map(lambda v: "" if pd.isna(v) else str(v))).agg("|".join, axis=1)


def program_col(df):
    """Labelled program, falling back to the program identified from the loop set."""
    return df["label_program"].where(df["label_program"].notna(), df.get("program_measured"))


# ------------------------------------------------------------------ A5 / A7
def align_residual(a, b, up=16, max_lag=40, f_max=0.38):
    """Align b to a (integer + fractional lag), then residual energy relative to a, with and without a gain fit.

    Both are limited to f < f_max x fs_dev (about 10 kHz): above that the deconvolution inverts the steep
    converter roll-off and amplifies noise differently for each kernel, which would mask true copies.
    """
    n = min(len(a), len(b))
    N = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfftfreq(N)
    w = np.where(f < 0.9 * f_max, 1.0, np.where(f < f_max, 0.5 * (1 + np.cos(np.pi * (f - 0.9 * f_max) / (0.1 * f_max))), 0.0))
    A = np.fft.rfft(a[:n].astype(float), N, axis=0) * w[:, None]
    B = np.fft.rfft(b[:n].astype(float), N, axis=0) * w[:, None]
    a = np.fft.irfft(A, N, axis=0)[:n]
    X = np.sum(A * np.conj(B), axis=1)
    xc = np.fft.irfft(X, N * up) * up
    lags = np.arange(-max_lag * up, max_lag * up + 1)
    vals = np.r_[xc[-max_lag * up:], xc[:max_lag * up + 1]]
    j = int(np.argmax(vals))
    tau = lags[j] / up
    if 0 < j < len(vals) - 1:
        y0, y1, y2 = vals[j - 1], vals[j], vals[j + 1]
        den = y0 - 2 * y1 + y2
        tau += (0.5 * (y0 - y2) / den if den else 0.0) / up
    bs = np.fft.irfft(B * np.exp(-2j * np.pi * f * tau)[:, None], N, axis=0)[:n]
    ea = np.sum(a ** 2, axis=0)
    g = np.sum(a * bs) / (np.sum(bs ** 2) + 1e-30)
    raw = np.sum((a - bs) ** 2, axis=0)
    fit = np.sum((a - g * bs) ** 2, axis=0)
    db = lambda num, den: float(10 * np.log10(np.sum(num) / (np.sum(den) + 1e-30) + 1e-30))
    return {"lag": float(tau), "gain_db": float(20 * np.log10(abs(g) + 1e-30)),
            "residual_db": db(fit, ea), "residual_raw_db": db(raw, ea),
            "residual_L_db": db(fit[0], ea[0]), "residual_R_db": db(fit[1], ea[1])}


def classify_pair(r, noise_ref_db, same_rate=True):
    if r["residual_db"] < noise_ref_db - 20:
        return "copy"
    if r["residual_db"] < noise_ref_db + 6:
        # across sample rates, a vendor's SRC passband ripple can leave a residual near the noise
        # floor, so a resampled copy and an independent capture can't be told apart this way
        # (unless each rate is known to be a separate recording: config separate_rate_captures; not the case for the Tomkin library, where crossrate.py shows one capture at both rates)
        return "repeat" if same_rate else "repeat_or_resampled"
    if r["residual_db"] < -20:
        return "near"
    return "different"


def candidate_pairs(df):
    """Same-label pairs, plus pairs whose *measured* settings are indistinguishable."""
    pairs = set()
    lab = label_key(df)
    has_label = df[LABELS].notna().any(axis=1)
    for _, idx in df[has_label].groupby(lab[has_label]).groups.items():
        pairs.update(itertools.combinations(sorted(idx), 2))
    prog = program_col(df).fillna("?")
    for _, g in df.groupby(prog):
        g = g.dropna(subset=["rt_mid"]).sort_values("rt_mid")
        ix, rt = g.index.to_list(), g["rt_mid"].to_numpy()
        for a in range(len(ix)):
            b = a + 1
            while b < len(ix) and rt[b] <= rt[a] * 1.02:
                ra, rb = df.loc[ix[a]], df.loc[ix[b]]
                close = all(abs((ra[c] or 0) - (rb[c] or 0)) <= 0.05 * max(abs(ra[c] or 1), 1e-3)
                            for c in ("lf_ratio", "hf_ratio") if pd.notna(ra[c]) and pd.notna(rb[c]))
                if ra["has_cluster"] != rb["has_cluster"]:
                    close = False
                elif ra["has_cluster"] and abs(ra["cluster_gain_db"] - rb["cluster_gain_db"]) > 1.0:
                    close = False
                if close:
                    pairs.add(tuple(sorted((ix[a], ix[b]))))
                b += 1
    return sorted(pairs)


def compare_pairs(df, sigs, separate_rate_captures=False):
    rows = []
    lab = label_key(df)
    for i, j in candidate_pairs(df):
        if sigs.get(i) is None or sigs.get(j) is None:
            continue
        r = align_residual(sigs[i], sigs[j])
        noise_ref = max(df.at[i, "noise_rel_db"], df.at[j, "noise_rel_db"])
        same_label = lab[i] == lab[j] and df.loc[i, LABELS].notna().any()
        rows.append({"a": df.at[i, "relpath"], "b": df.at[j, "relpath"], "ia": i, "ib": j,
                     "fs_a": df.at[i, "fs"], "fs_b": df.at[j, "fs"], "same_label": same_label,
                     "noise_ref_db": noise_ref,
                     "class": classify_pair(r, noise_ref, df.at[i, "fs"] == df.at[j, "fs"] or
                                               (separate_rate_captures and same_label)), **r})
    return pd.DataFrame(rows)


def pair_flags(df, pairs):
    """A5 / A7 per-file flags from the pair table."""
    add = {i: [] for i in df.index}
    no_delay = [c for c in LABELS if c != "label_delay"]
    key_nd = df[no_delay].apply(lambda col: col.map(lambda v: "" if pd.isna(v) else str(v))).agg("|".join, axis=1)
    unlabelled = ~df[LABELS].notna().any(axis=1)
    # same-label groups: blame only files that differ from most of their group
    same = pairs[pairs["same_label"]]
    for i in df.index:
        mine = same[(same["ia"] == i) | (same["ib"] == i)]
        if len(mine) and (mine["class"] == "different").sum() > len(mine) / 2:
            add[i].append("label_group_inconsistent")
    # different-label duplicates: blame the file whose duplicates mostly carry other labels
    # (anchoring removes pre-delay, so delay-only label differences are expected and ignored)
    dup = pairs[pairs["class"].isin(["copy", "repeat", "repeat_or_resampled"])]
    for i in df.index:
        if unlabelled[i]:
            if ((dup["ia"] == i) | (dup["ib"] == i)).any():
                add[i].append("unlabelled_duplicate_of_labelled_file")
            continue
        partners = [b if a == i else a for a, b in zip(dup["ia"], dup["ib"]) if i in (a, b)]
        partners = [j for j in partners if not unlabelled[j]]
        other = sum(key_nd[j] != key_nd[i] for j in partners)
        if other and other >= len(partners) - other:
            add[i].append("duplicate_with_different_label")
    for _, p in pairs.iterrows():
        a, b = p["ia"], p["ib"]
        if p["same_label"] and p["class"] == "copy":
            tag = "resampled_copy" if p["fs_a"] != p["fs_b"] else "duplicate_file"
            add[b].append(tag)
        if p["same_label"] and p["class"] == "repeat":
            add[a].append("has_repeat_capture"); add[b].append("has_repeat_capture")
        if p["same_label"] and p["class"] == "repeat_or_resampled":
            add[a].append("other_rate_version"); add[b].append("other_rate_version")
    return {i: sorted(set(v)) for i, v in add.items()}


def normalisation(df):
    """A6: were files level-normalised by the vendor?"""
    peak = df[["peak_dbfs_L", "peak_dbfs_R"]].max(axis=1).round(1)
    mode = peak.mode().iloc[0] if len(peak) else None
    frac_at_mode = float((abs(peak - mode) <= 0.05).mean()) if mode is not None else 0.0
    out = {"peak_mode_dbfs": None if mode is None else float(mode), "frac_files_at_peak_mode": frac_at_mode,
           "peak_normalised_likely": frac_at_mode >= 0.5 and len(df) >= 4}
    c = df[df["has_cluster"] == True].copy()  # noqa: E712
    c["prog"] = program_col(c)
    if "label_amp" in c and c["label_amp"].notna().any():
        g = c.dropna(subset=["label_amp"]).groupby(["prog", "label_amp"])["cluster_gain_db"]
        spread = g.agg(lambda s: float(s.max() - s.min()) if len(s) > 1 else np.nan).dropna()
        out["cluster_gain_spread_db_median"] = float(spread.median()) if len(spread) else None
        out["cluster_gain_inconsistent"] = bool(len(spread) and spread.median() > 0.5)
    else:
        out["cluster_gain_spread_db_median"] = None
        out["cluster_gain_inconsistent"] = None
    out["amp_measure"] = "cluster_to_reverb_db" if (out["peak_normalised_likely"] or out["cluster_gain_inconsistent"]) else "cluster_gain_db"
    return out


# ------------------------------------------------------------------ B2 – B5
def build_kernel(df, fs, ref, n_max=8, kernel_ms=5.2):
    """B2: average the R-channel tap-0 waveform over the cleanest cluster IRs at this sample rate."""
    bad = df["flags"].apply(lambda f: bool({"clipped", "mono", "error", "fade_out"} & set(f)))
    c = df[(df["fs"] == fs) & (df["cluster_status"] == "ok") & ~bad]
    c = c.sort_values(["cluster_sparsity", "cluster_to_reverb_db"], ascending=False).head(n_max)
    if not len(c):
        return None, None
    L = int(kernel_ms / 1000 * fs)
    gs = []
    for _, r in c.iterrows():
        x, _ = sf.read(r["path"], always_2d=True)
        s = x[:, 1] if x.shape[1] > 1 else x[:, 0]
        t0 = r["tap0_ms"] / 1000 * fs
        i0 = int(np.floor(t0))
        frac = t0 - i0
        s = np.r_[np.zeros(64), s]  # guard so files trimmed right at the onset keep the 64-sample lead
        seg = s[i0:i0 + L + 128]
        N = 1 << int(np.ceil(np.log2(2 * len(seg))))
        f = np.fft.rfftfreq(N)
        sh = np.fft.irfft(np.fft.rfft(seg, N) * np.exp(2j * np.pi * f * frac), N)[:len(seg)]
        gs.append(sh[64:64 + L])
    nt = int(0.0015 * fs)
    taper = np.ones(L); taper[-nt:] = np.hanning(2 * nt)[nt:]
    gs = [q * taper for q in gs]
    gains = np.array([np.sqrt(np.sum(q ** 2)) for q in gs])
    gs = [q / gn * np.median(gains) for q, gn in zip(gs, gains)]  # average shapes at the median level
    g = np.mean(gs, axis=0)
    k = Kernel(g, fs, source=f"library mean of {len(gs)} R tap-0 waveforms")
    # comparison with the reference kernel (in-band magnitude and excess group delay)
    gr = ref.at(fs)
    N = 1 << 16
    f = np.fft.rfftfreq(N, 1 / fs)
    m = (f > 100) & (f < 11000)
    H = np.fft.rfft(g, N) / np.fft.rfft(gr, N)
    mag = 20 * np.log10(np.abs(H[m]))
    ph = np.unwrap(np.angle(H[m]))
    p = np.polyfit(f[m], ph, 1)
    exc = ph - np.polyval(p, f[m])
    spread = [20 * np.log10(np.abs(np.fft.rfft(q, N)[m]) / np.abs(np.fft.rfft(g, N)[m])) for q in gs]
    cmp = {"n_files": len(gs), "gain_vs_ref_db": float(np.median(mag)),
           "candidate_gain_spread_db": float(20 * np.log10(gains.max() / gains.min())),
           "mag_dev_vs_ref_db_max": float(np.max(np.abs(mag - np.median(mag)))),
           "excess_phase_vs_ref_rad_max": float(np.max(np.abs(exc))),
           "delay_vs_ref_us": float(-p[0] / (2 * np.pi) * 1e6),
           "within_library_spread_db_max": float(np.max(np.abs(spread))) if spread else None}
    return k, cmp


def latency_regression(df):
    """B3: first cluster tap time against the labelled pre-delay."""
    c = df[(df["has_cluster"] == True) & df["label_delay"].notna() & df["tap0_ms"].notna()]  # noqa: E712
    out = {"n": int(len(c)), "distinct_delays": int(c["label_delay"].nunique()) if len(c) else 0}
    if out["distinct_delays"] < 3:
        if len(df[df["has_cluster"] == True]) >= 3 and df.loc[df["has_cluster"] == True, "tap0_ms"].std() < 0.05:  # noqa: E712
            out["status"] = "files_trimmed_to_onset_or_single_delay"
        else:
            out["status"] = "not_enough_delay_values"
        return out
    x, y = c["label_delay"].to_numpy(float), c["tap0_ms"].to_numpy(float)
    # robust: fit through the per-delay medians, then flag files far from the line
    ux = np.unique(x)
    p = np.polyfit(ux, [np.median(y[x == v]) for v in ux], 1)
    res = y - np.polyval(p, x)
    mad = np.median(np.abs(res - np.median(res))) * 1.4826
    keep = np.abs(res) < max(4 * mad, 0.2)
    out.update({"status": "ok", "slope": float(p[0]), "intercept_ms": float(p[1]),
                "rms_residual_ms": float(np.sqrt(np.mean(res[keep] ** 2))),
                "outliers": c.index[~keep].to_list(), "residual_ms": dict(zip(c.index, res))})
    if y.std() < 0.05:
        out["status"] = "files_trimmed_to_onset"
    elif abs(p[0] - 1) > 0.01:
        out["status"] = "slope_not_unity"
    return out


def lr_offsets(df):
    mono = df["flags"].apply(lambda f: "mono" in f)
    c = df[(df["cluster_status"] == "ok") & ~mono].dropna(subset=["lr_offset_us"])
    if not len(c):
        return {"n": 0}
    med = float(c["lr_offset_us"].median())
    return {"n": int(len(c)), "median_us": med, "mad_us": float((c["lr_offset_us"] - med).abs().median()),
            "outliers": c.index[(c["lr_offset_us"] - med).abs() > 2.0].to_list()}


def clock(df):
    c = df[df["has_cluster"] == True].dropna(subset=["clock_ppm"])  # noqa: E712
    out = {}
    for fs, g in c.groupby("fs"):
        out[int(fs)] = {"n": int(len(g)), "median_ppm": float(g["clock_ppm"].median()),
                        "iqr_ppm": float(g["clock_ppm"].quantile(0.75) - g["clock_ppm"].quantile(0.25)),
                        "off_grid_files": int((g["grid_max_frac"] > 0.25).sum())}
    return out


def unit_fingerprint(df, cfg, lr, clk):
    ref = cfg.get("reference_unit", {})
    out = {"reference": ref}
    out["lr_offset_diff_us"] = (lr.get("median_us") - ref["lr_offset_us"]) if lr.get("n") and "lr_offset_us" in ref else None
    nd = (df["noise_dbfs_L"] - df["noise_dbfs_R"]).median()
    out["noise_L_minus_R_db"] = float(nd) if pd.notna(nd) else None
    progs = df.get("program_measured", pd.Series(dtype=str)).fillna("?")
    out["loopset_counts"] = progs.value_counts().to_dict()
    statuses = df["cluster_status"].value_counts().to_dict()
    out["cluster_status_counts"] = statuses
    notes = []
    if out["lr_offset_diff_us"] is not None and abs(out["lr_offset_diff_us"]) > 2:
        notes.append("L/R output offset differs from the reference unit")
    unknown = [p for p in out["loopset_counts"] if str(p).startswith("U:")]
    if unknown:
        notes.append(f"{len(unknown)} loop set(s) not in config.known_loopsets (expected for M/L until added; a set close to S would suggest different firmware)")
    if statuses.get("novel_cluster"):
        notes.append("some cluster patterns differ from the reference templates (other programs, or different firmware)")
    out["notes"] = notes
    return out


# ------------------------------------------------------------------ C3 – C5
def gap_clusters(values, scale="lin", rel_tol=0.05, abs_tol=None):
    v = np.sort(np.asarray([q for q in values if q is not None and np.isfinite(q)], float))
    if not len(v):
        return []
    t = np.log(v) if scale == "log" else v
    groups, cur = [], [0]
    for i in range(1, len(v)):
        gap = t[i] - t[i - 1]
        lim = abs_tol if abs_tol is not None else (rel_tol if scale == "log" else rel_tol * max(abs(v[i - 1]), 1e-9))
        if gap > lim:
            groups.append(cur); cur = []
        cur.append(i)
    groups.append(cur)
    return [{"center": float(np.median(v[g])), "n": len(g), "min": float(v[g].min()), "max": float(v[g].max())} for g in groups]


def step_tables(df, cfg, amp_measure):
    out = {}
    prog = program_col(df).fillna("?")
    meas = dict(MEASURED, amp=amp_measure)
    for p, gdf in df.groupby(prog):
        tables = {}
        for param, col in meas.items():
            if col not in gdf or gdf[col].notna().sum() == 0:
                continue
            spec = cfg["parameters"][param]
            lab = f"label_{param}"
            entry = {"measured_col": col}
            if lab in gdf and gdf[lab].notna().any():
                g = gdf.dropna(subset=[lab, col]).groupby(lab)[col]
                t = g.agg(["count", "median", "min", "max"])
                t["mad"] = g.agg(lambda s: float((s - s.median()).abs().median()))
                entry["labelled"] = t.reset_index().rename(columns={lab: "label"}).to_dict("records")
                med = t["median"].to_numpy()
                entry["monotonic"] = bool(len(med) < 2 or np.all(np.diff(med) > 0) or np.all(np.diff(med) < 0))
            if param == "amp":
                tol = dict(abs_tol=spec["abs_tol_db"])
            elif param == "delay":
                tol = dict(abs_tol=spec["abs_tol"])
            else:
                tol = dict(rel_tol=spec["rel_tol"], scale=spec["scale"])
            cl = gap_clusters(gdf[col].to_list(), **tol)
            exp = spec.get("steps_per_program", spec.get("steps"))
            entry["clusters"] = cl
            entry["n_clusters"] = len(cl)
            entry["expected_steps"] = exp
            tables[param] = entry
        out[p] = tables
    return out


def disagreements(df, cfg, amp_measure, b3):
    rows = []
    prog = program_col(df).fillna("?")
    meas = dict(MEASURED, amp=amp_measure)
    for param, col in meas.items():
        lab = f"label_{param}"
        if lab not in df or col not in df:
            continue
        spec = cfg["parameters"][param]
        sub = df.dropna(subset=[lab, col]).assign(_prog=lambda d: prog[d.index])
        extra = [k for k in ("label_rt", "label_lf", "label_hf") if param == "amp" and sub[k].notna().any()]
        sub = sub.assign(**{f"_{k}": sub[k].map(lambda v: "" if pd.isna(v) else str(v)) for k in extra})
        for keys, g in sub.groupby(["_prog", lab] + [f"_{k}" for k in extra]):
            p, v = keys[0], keys[1]
            if param == "delay" and b3.get("status") == "ok":
                dev = pd.Series({i: b3["residual_ms"].get(i, 0.0) for i in g.index})
                bad = dev.abs() > spec["abs_tol"]
            elif len(g) >= 3:
                med = g[col].median()
                dev = g[col] - med
                if param == "amp":
                    bad = dev.abs() > spec["abs_tol_db"]
                elif param == "delay":
                    bad = dev.abs() > spec["abs_tol"]
                else:
                    bad = (dev / med).abs() > spec["rel_tol"]
            elif param == "rt":  # RT is in seconds, so a lone file can be checked against its own label
                dev = g[col] - v
                bad = (dev / v).abs() > 0.15 if v else dev.abs() > 0.15
            else:
                continue
            for i in g.index[bad.to_numpy()]:
                rows.append({"relpath": df.at[i, "relpath"], "param": param, "label": v,
                             "measured": float(df.at[i, col]), "deviation": float(dev[i]), "program": p})
    # program label against loop-set identity
    if "program_measured" in df:
        known = set(cfg["known_loopsets"])
        for i, r in df.iterrows():
            lp_, mp = r["label_program"], r.get("program_measured")
            if pd.notna(lp_) and isinstance(mp, str) and mp in known and lp_ != mp:
                rows.append({"relpath": r["relpath"], "param": "program", "label": lp_, "measured": mp,
                             "deviation": None, "program": lp_})
    return pd.DataFrame(rows)


def _value(df, param):
    return df[f"label_{param}"]


def coverage(df, cfg):
    prog = program_col(df).fillna("?")
    out = {"files_per_program": prog.value_counts().to_dict(), "planned": {}, "matrices": {}}
    neutral = {k: v["neutral"] for k, v in cfg["parameters"].items() if not k.startswith("_")}
    for name, s in cfg["planned_sets"].items():
        if name.startswith("_"):
            continue
        res = {}
        vary = f"label_{s['vary']}"
        fixed = [f"label_{k}" for k in neutral if k != s["vary"]]
        for p in s["programs"]:
            g = df[prog == p].dropna(subset=[vary])
            best = {"have": 0, "need": s["need"], "values": [], "held": {}}
            if len(g):
                keys = g[fixed].apply(lambda col: col.map(lambda v: "?" if pd.isna(v) else str(v)))
                for k, grp in g.groupby(keys.agg("|".join, axis=1)):
                    vals = sorted(float(v) for v in grp[vary].unique())
                    if len(vals) > best["have"]:
                        best = {"have": len(vals), "need": s["need"], "values": vals,
                                "held": dict(zip([f[6:] for f in fixed], k.split("|")))}
            res[p] = best
        out["planned"][name] = res
    for p, g in df.groupby(prog):
        mats = {}
        for other in ("lf", "hf", "amp", "delay", "lp"):
            sub = g.dropna(subset=["label_rt", f"label_{other}"])
            if len(sub):
                mats[other] = pd.crosstab(sub["label_rt"], sub[f"label_{other}"])
        out["matrices"][p] = mats
    return out
