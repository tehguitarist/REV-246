"""Render the plugin at every Altiverb factory preset's snapped settings and compare with the preset's own Altiverb IR.

    .venv/bin/python tools/memguard.py --limit-gb 2 -- .venv/bin/python tools/altiverb_preset_check.py [--stage measure|report|all]
        [--only <substring> ...]

Inputs: analysis/altiverb/presets.json (snapped values), analysis/altiverb/cache/*.npy (deconvolved Altiverb IRs, 48 kHz,
zero latency), plugin/build-render/REV246Render_artefacts/Release/REV246Render (render at rate 48000, mod=1, sixteen=0:
the plugin default 1 adds a -96 dBFS floor that wrecks Schroeder T60 of an 0 dBFS impulse render, kept as render16on; the engine latency printed per render is removed). Output: analysis/altiverb/preset_check/
(renders/*.wav, measure.json, report.md, summary.json).

Metrics (the IR and the render pass through the same function; time-variant programs M / L / E: octave T60 and EDT are the
mean of three renders at different input phases, so they are statistics, never nulls):
  T60 / EDT   octave bands 125 Hz-8 kHz (altiverb_compare.band_t60's fit, -5..-25 dB; EDT 0..-10 dB), per-band ln(render / IR)
  AMPL        cluster-tap / reverb ratio (S / M / L: altiverb_presets.tap_ratio_at; E: e_ratio), render over IR in dB
  DELAY       onset (altiverb_presets.classify) of the render minus the IR's, ms
  LOW PASS    tap-0 spectrum over the family reference's (1-12 kHz), IR minus render
  N           5 ms band envelopes (altiverb_presets.n_envelope_data), rms dB residual of the render against the IR
Family-relative: per program the median over its family presets of each log ratio (T60 per band, AMPL, EDT) is the
"chain bias" (Altiverb's 8 kHz anomaly, low-band chain, clusters' level convention) and is removed before judging.
Flags use the snapping step: RT miss > half the RT step's ratio, LF / HF > half their step's ratio (band residual after the
RT part), DELAY > 0.5 ms, AMPL > 5 (half the 10 % step), LOW PASS > 1.5 dB at 10 kHz relative to the family.
Flagged presets are re-rendered at the suggested step to verify the suggestion.
"""
import argparse, json, os, subprocess, sys
import numpy as np
import soundfile as sf
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import altiverb_compare as ac  # noqa: E402
import altiverb_presets as ap  # noqa: E402

FS = ac.FS
OUT = os.path.join(ROOT, "analysis", "altiverb", "preset_check")
RENDERS = os.path.join(OUT, "renders")
EXE = os.path.join(ROOT, "plugin", "build-render", "REV246Render_artefacts", "Release", "REV246Render")
PRESETS = json.load(open(os.path.join(ROOT, "analysis", "altiverb", "presets.json")))["presets"]
LEN12 = 12 * FS
PRES = [0.0, 0.37, 0.9]          # input phases for the time-variant programs (seconds, integer samples at 48 kHz)
BANDS = [butter(4, [fc / 2 ** 0.5, min(fc * 2 ** 0.5, 0.45 * FS)], btype="band", fs=FS, output="sos") for fc in ac.CF]


# ---------------------------------------------------------------- rendering
def render(prog, st, tag, pre=0.0, length=12.5, sixteen=0):
    """Render the plugin at settings st ({delay, lowpass, rt, lf, hf, ampl}), remove latency (and pre): (2, n) float64."""
    path = os.path.join(RENDERS, tag + ".wav")
    a = [EXE, f"out={path}", "rate=48000", f"program={prog}", f"lf={st['lf']}", f"hf={st['hf']}", f"ampl={st['ampl']}",
         f"delay={st['delay']}", f"lp={st['lowpass']}", f"sixteen={sixteen}", "mod=1", f"len={length + pre}", f"pre={pre}"]
    if st.get("rt") is not None:
        a.append(f"rt={st['rt']}")
    o = subprocess.run(a, capture_output=True, text=True, check=True).stdout
    lat = int(o.split("latency=")[1].split()[0])
    x, fs = sf.read(path, always_2d=True)
    assert fs == FS
    k = lat + int(round(pre * FS))
    return x[k:].T.astype(np.float64), lat


# ---------------------------------------------------------------- metrics
def decay(ir):
    """Octave T60 and EDT per band (channel mean): lists of 7."""
    t60, edt = [], []
    for sos in BANDS:
        a, b = [], []
        for c in (0, 1):
            x = sosfiltfilt(sos, ir[c][:LEN12])
            e = np.cumsum((x ** 2)[::-1])[::-1]
            e = 10 * np.log10(e / e[0] + 1e-30)
            tt = np.nonzero((e < -5) & (e > -25))[0]
            a.append(np.nan if len(tt) < 200 or tt[-1] / FS > len(x) / FS - 0.02 else -60 / np.polyfit(tt / FS, e[tt], 1)[0])
            te = np.nonzero((e <= 0) & (e > -10))[0]
            b.append(-60 / np.polyfit(te / FS, e[te], 1)[0] if len(te) > 20 else np.nan)
        t60.append(float(np.nanmean(a)) if not np.all(np.isnan(a)) else np.nan)
        edt.append(float(np.nanmean(b)) if not np.all(np.isnan(b)) else np.nan)
    return t60, edt


def tap_spec(ir, prog, onset_ms):
    ch = 1 if prog == "E" else 0
    k = int(round(onset_ms * 1e-3 * FS + (41 / ap.FSD * FS if prog == "N" else 0)))
    k = k + int(np.argmax(np.abs(ir[ch][k - 3: k + 4]))) - 3
    return ap.tap_spectrum(ir, ch, k, ap.HALF[prog])[: int(12500 / (FS / (1 << 15)))]


def features(ir, want_decay=True):
    prog, onset, how = ap.classify(ir)
    r = {"program": prog, "onset_ms": onset, "how": how}
    if want_decay and prog != "N":
        r["t60"], r["edt"] = decay(ir)
    if prog in ("S", "M", "L"):
        has = not how.startswith("no cluster")
        r["tap_ratio"] = ap.tap_ratio_at(ir, prog, int(round(onset * 1e-3 * FS))) if has else 0.0
    if prog == "E":
        r["e_ratio"] = ap.e_ratio(ir, int(round(onset * 1e-3 * FS)))
    if prog is not None and (prog != "L" or "no cluster" not in how):
        r["tap_spec"] = tap_spec(ir, prog, onset).tolist()
    if prog == "N":
        r["n_env"] = ap.n_envelope_data(ir, onset).tolist()
    return r


def render_features(prog, st, tag, tv):
    """Features of the render; decay averaged (log) over the input phases for the time-variant programs."""
    pres = PRES if tv else PRES[:1]
    ir0, lat = render(prog, st, tag, 0.0)
    f = features(ir0, want_decay=False)
    f["latency"] = lat
    if prog != "N":
        ts, es = [], []
        for i, p in enumerate(pres):
            ir = ir0 if p == 0.0 else render(prog, st, f"{tag}~pre{p}", p)[0]
            t, e = decay(ir)
            ts.append(t)
            es.append(e)
        ts, es = np.array(ts), np.array(es)
        f["t60"] = np.exp(np.nanmean(np.log(ts), 0)).tolist()
        f["edt"] = np.exp(np.nanmean(np.log(es), 0)).tolist()
        f["t60_sd"] = (np.nanstd(np.log(ts), 0) if len(ts) > 1 else np.zeros(7)).tolist()
        f["edt_sd"] = (np.nanstd(np.log(es), 0) if len(es) > 1 else np.zeros(7)).tolist()
    return f


def safe(n):
    return n.replace(" ", "_").replace("/", "_")


def settings_of(name):
    s = PRESETS[name]["settings"]
    return {k: (s[k]["value"] if s[k].get("value") is not None else None) for k in ("delay", "lowpass", "rt", "lf", "hf", "ampl")}


# ---------------------------------------------------------------- N: local search with the plugin's own render
def n_cost(env_d, env_m):
    cs, ns = 0.0, 0
    for ch in (0, 1):
        d, m = np.array(env_d[ch]), np.array(env_m[ch])
        sel = d > d.max() - 50
        r = d[sel] - np.maximum(m[sel], m.max() - 90)
        cs += np.sum((r - r.mean()) ** 2)
        ns += sel.sum()
    return float(np.sqrt(cs / ns))


def n_neighbours(st):
    early = ap.N_EARLY
    ie = min(range(len(early)), key=lambda i: abs(early[i] - st["lf"]))
    out = []
    for dsp in (-0.1, -0.05, 0, 0.05, 0.1):
        sp = round(st["ampl"] + dsp, 2)
        if sp < 0 or sp > 0.5:
            continue
        for de in (-1, 0, 1):
            if not 0 <= ie + de < len(early):
                continue
            for dl in (-0.2, 0, 0.2):
                la = round(st["hf"] + dl, 2)
                if la < 0 or la > 2.0:
                    continue
                out.append({**st, "ampl": sp, "lf": early[ie + de], "hf": la})
    return out


def fa_of(ir):
    """First arrival (48 kHz samples): first sample above 1 % of the largest peak in the first 0.6 s, both channels."""
    x = np.abs(np.array(ir[:, : int(0.6 * FS)])).max(0)
    return int(np.argmax(x > 0.01 * x.max()))


def stage_fa():
    """Add the first-arrival onset of IR and render to measure.json (classification-independent DELAY check)."""
    path = os.path.join(OUT, "measure.json")
    res = json.load(open(path))
    for n, r in res.items():
        ir = np.load(os.path.join(ap.CACHE, n + ".npy"), mmap_mode="r")
        rd, _ = render(r["program"], r["settings"], safe(n) + "~fa", 0.0, length=1.0)
        r["fa_ir"], r["fa_render"] = fa_of(ir), fa_of(rd)
        r["render_program"] = ap.classify(rd)[0]
    json.dump(res, open(path, "w"))


def stage_measure(only):
    path = os.path.join(OUT, "measure.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    names = [n for n, v in PRESETS.items() if v["included"] and (not only or any(o in n for o in only))]
    for n in names:
        prog = PRESETS[n]["program"]
        st = settings_of(n)
        ir = np.load(os.path.join(ap.CACHE, n + ".npy"))
        r = {"program": prog, "settings": st, "ir": features(ir)}
        r["render"] = render_features(prog, st, safe(n), prog in ("M", "L", "E"))
        r["render16on"] = None
        if prog != "N":
            ir0, _ = render(prog, st, safe(n) + "~s16", 0.0, sixteen=1)
            r["render16on"] = {"t60": decay(ir0)[0], "onset_ms": ap.classify(ir0)[1]}
        # LOW PASS candidates against the family reference (spectrum only, short render)
        r["lp_cand"] = {}
        for lp in (6, 8, 10, 12):
            ir0, _ = render(prog, {**st, "lowpass": lp}, safe(n) + f"~lp{lp}", 0.0, length=0.4)
            o = ap.classify(ir0)[1]
            r["lp_cand"][str(lp)] = tap_spec(ir0, prog, o).tolist()
        if prog == "N":
            cand = n_neighbours(st)
            best = []
            for i, c in enumerate(cand):
                ir0, _ = render("N", c, safe(n) + "~ng", 0.0, length=2.2)
                o = ap.classify(ir0)[1]
                best.append((n_cost(ap.n_envelope_data(ir, ap.classify(ir)[1]), ap.n_envelope_data(ir0, o)), c))
            best.sort(key=lambda t: t[0])
            r["n_search"] = [{"cost": b[0], "st": b[1]} for b in best[:8]]
            r["n_search_chosen"] = next(b[0] for b in best if b[1] == st)
        res[n] = r
        print(f"{n:24s} {prog}  onset IR {r['ir']['onset_ms']:.2f} render {r['render']['onset_ms']:.2f}  latency {r['render']['latency']}", flush=True)
        json.dump(res, open(path, "w"))
    return res


# ---------------------------------------------------------------- judging
def lnr(a, b):
    return np.log(np.array(a, float) / np.array(b, float))


def family_key(n, prog):
    f = ap.family_of(n)
    return f if f else None


def step_half(table, val, direction):
    """Half the log ratio of the step from val towards direction (+1 up, -1 down); falls back to the other side at an end."""
    i = min(range(len(table)), key=lambda j: abs(np.log(table[j] / val)))
    j = i + (1 if direction > 0 else -1)
    if not 0 <= j < len(table):
        j = i - (1 if direction > 0 else -1)
    return abs(np.log(table[j] / table[i])) / 2


def nearest(table, val):
    return min(table, key=lambda t: abs(np.log(t / val)))


def metrics(m):
    """Per-preset raw render/IR comparison (no family removal yet)."""
    ir, rd, prog = m["ir"], m["render"], m["program"]
    o = {"onset_diff_ms": rd["onset_ms"] - ir["onset_ms"]}
    if prog != "N":
        o["t60_ln"] = lnr(rd["t60"], ir["t60"]).tolist()
        o["edt_ln"] = lnr(rd["edt"], ir["edt"]).tolist()
    if prog in ("S", "M", "L") and ir["tap_ratio"] > 0:
        o["ampl_db"] = 20 * np.log10(rd["tap_ratio"] / ir["tap_ratio"])
    if prog == "E":
        o["ampl_db"] = 20 * np.log10(rd["e_ratio"] / ir["e_ratio"])
    return o


def lp_eval(m, ref, name):
    """LP: IR relative to family ref IR vs render relative to family ref render; per candidate. ref = measure of the reference."""
    f = ap.F_SPEC[: int(12500 / (FS / (1 << 15)))]
    ir_rel = 20 * np.log10(np.array(m["ir"]["tap_spec"]) / np.array(ref["ir"]["tap_spec"]))
    sel = (f >= 1000) & (f <= 12000)
    out = {}
    cur = str(int(m["settings"]["lowpass"]))
    base = np.array(ref["render"]["tap_spec"])
    for lp, sp in m["lp_cand"].items():
        rr = 20 * np.log10(np.array(sp) / base)
        d = ir_rel - rr
        lvl = float(np.mean(d[(f >= 1000) & (f <= 2000)]))
        d = d - lvl                                   # level offset (tap gain) removed: shape only, referred to 1-2 kHz
        out[lp] = {"level_db": lvl, "rms": float(np.sqrt(np.mean(d[sel] ** 2))),
                   "d10": float(np.interp(10000, f, d)), "d8": float(np.interp(8000, f, d)), "d12": float(np.interp(12000, f, d))}
    return out


def judge(n, m, bias, res):
    prog = m["program"]
    st = m["settings"]
    o = metrics(m)
    row = {"program": prog, "settings": st, "label": PRESETS[n]["label_string"], "flags": [], "suggest": {}, "notes": []}
    row["onset_diff_ms"] = o["onset_diff_ms"]
    row["latency"] = m["render"]["latency"]
    row["onset_fa_diff_ms"] = (m["fa_render"] - m["fa_ir"]) / FS * 1e3 if "fa_ir" in m else None
    od = o["onset_diff_ms"]
    if m.get("render_program", prog) != prog and row["onset_fa_diff_ms"] is not None:
        od = row["onset_fa_diff_ms"]          # the render no longer matches the program's tap template (steep LP): first arrival
        row["notes"].append("render's tap template lost (LP): onset from the first arrival")
    row["onset_used_diff_ms"] = od
    if abs(od) > 0.5:
        row["flags"].append("DELAY")
        row["suggest"]["delay"] = int(round(st["delay"] - od))
    if prog != "N":
        b = bias[prog]
        rel = np.array(o["t60_ln"]) - b["t60"]          # render/IR, family chain removed; + = plugin too long
        erel = np.array(o["edt_ln"]) - b["edt"]
        row["t60_rel_pct"] = (np.expm1(rel) * 100).round(1).tolist()
        row["edt_rel_pct"] = (np.expm1(erel) * 100).round(1).tolist()
        row["t60_raw_pct"] = (np.expm1(o["t60_ln"]) * 100).round(1).tolist()
        row["t60_ir"] = np.round(m["ir"]["t60"], 3).tolist()
        row["t60_render"] = np.round(m["render"]["t60"], 3).tolist()
        row["t60_sd_pct"] = (np.array(m["render"]["t60_sd"]) * 100).round(1).tolist()
        idx = {125: 0, 250: 1, 500: 2, 1000: 3, 2000: 4, 4000: 5, 8000: 6}
        mid = [idx[k] for k in (500, 1000, 2000)]
        g = float(np.nanmean(rel[mid]))
        rt = st["rt"]
        tab = ap.RT_STEPS[prog]
        need = rt * np.exp(-g)
        hs = step_half(tab, rt, -1 if g > 0 else 1)
        row["rt_miss_pct"] = float(np.expm1(g) * 100)
        row["rt_half_step_pct"] = float(np.expm1(hs) * 100)
        if abs(g) > hs and not np.isnan(g):
            row["flags"].append("RT")
            row["suggest"]["rt"] = nearest(tab, need)
        lfd = float(np.nanmean(rel[[0, 1]]) - g)
        hfd = float(np.nanmean(rel[[4, 5]]) - g)     # 2 and 4 kHz (8 kHz excluded: Altiverb anomaly)
        row["lf_miss_pct"], row["hf_miss_pct"] = float(np.expm1(lfd) * 100), float(np.expm1(hfd) * 100)
        lft, hft = ap.LF_STEPS, ap.HF_STEPS
        lhs = step_half(lft, st["lf"], -1 if lfd > 0 else 1)
        hhs = step_half(hft, st["hf"], -1 if hfd > 0 else 1)
        row["lf_half_step_pct"], row["hf_half_step_pct"] = float(np.expm1(lhs) * 100), float(np.expm1(hhs) * 100)
        if abs(lfd) > lhs and not np.isnan(lfd):
            row["flags"].append("LF")
            row["suggest"]["lf"] = nearest(lft, st["lf"] * np.exp(-lfd))
        if abs(hfd) > hhs and not np.isnan(hfd):
            row["flags"].append("HF")
            row["suggest"]["hf"] = nearest(hft, st["hf"] * np.exp(-hfd))
        row["edt_mid_pct"] = float(np.expm1(np.nanmean(erel[mid])) * 100)
        row["rel8k_pct"] = float(np.expm1(rel[6]) * 100)
        row["sixteen_on_t60_max_pct"] = float(np.nanmax(np.abs(np.expm1(lnr(m["render16on"]["t60"], m["render"]["t60"]))) * 100))
        if "ampl_db" in o:
            rd = o["ampl_db"] - b["ampl"]
            row["ampl_rel_db"] = float(rd)
            row["ampl_raw_db"] = float(o["ampl_db"])
            imp = st["ampl"] * 10 ** (-rd / 20)
            row["ampl_implied"] = float(imp)
            if abs(imp - st["ampl"]) > 5:
                row["flags"].append("AMPL")
                row["suggest"]["ampl"] = int(np.clip(round(imp / 10) * 10, 0, 100))
        # LOW PASS vs the family reference
        refn = ap.FAMILY_REF[prog]
        if n != refn and "tap_spec" in m["ir"]:
            lpe = lp_eval(m, res[refn], n)
            cur = str(int(st["lowpass"]))
            best = min(lpe, key=lambda k: lpe[k]["rms"])
            row["lp"] = {"current": lpe[cur], "best": int(best), "best_rms": lpe[best]["rms"]}
            if int(best) != int(cur) and lpe[cur]["rms"] - lpe[best]["rms"] > 1.0:
                row["flags"].append("LOWPASS")
                row["suggest"]["lowpass"] = int(best)
        elif n != refn:
            row["notes"].append("no cluster / tap 0: LOW PASS not measurable")
    else:
        row["n_cost_db"] = m["n_search_chosen"]
        bs = m["n_search"][0]
        row["n_best"] = {"cost": bs["cost"], **{k: bs["st"][k] for k in ("ampl", "lf", "hf")}}
        if m["n_search_chosen"] - bs["cost"] > 0.5 and bs["cost"] < 0.85 * m["n_search_chosen"]:
            row["flags"].append("N")
            row["suggest"].update({"ampl": bs["st"]["ampl"], "lf": bs["st"]["lf"], "hf": bs["st"]["hf"]})
        refn = "NON-LIN"
        if n != refn:
            lpe = lp_eval(m, res[refn], n)
            cur = str(int(st["lowpass"]))
            best = min(lpe, key=lambda k: lpe[k]["rms"])
            row["lp"] = {"current": lpe[cur], "best": int(best), "best_rms": lpe[best]["rms"]}
            if int(best) != int(cur) and lpe[cur]["rms"] - lpe[best]["rms"] > 1.0:
                row["flags"].append("LOWPASS")
                row["suggest"]["lowpass"] = int(best)
    return row


def stage_report():
    res = json.load(open(os.path.join(OUT, "measure.json")))
    rows, flagged = {}, []
    bias = {}
    for prog in ("S", "M", "L", "E"):
        fam = [n for n in res if ap.family_of(n) == prog]
        ms = {n: metrics(res[n]) for n in fam}
        b = {"t60": np.nanmedian([ms[n]["t60_ln"] for n in fam], 0), "edt": np.nanmedian([ms[n]["edt_ln"] for n in fam], 0),
             "ampl": float(np.median([ms[n]["ampl_db"] for n in fam])), "n": len(fam)}
        bias[prog] = b
    for n, m in res.items():
        row = judge(n, m, bias, res)
        rows[n] = row
        if row["flags"]:
            flagged.append(n)
    json.dump({"bias": {k: {kk: (vv.tolist() if hasattr(vv, "tolist") else vv) for kk, vv in v.items()} for k, v in bias.items()},
               "rows": rows}, open(os.path.join(OUT, "summary_raw.json"), "w"), indent=1, default=float)
    return res, rows, bias, flagged


def verify(res, rows, bias):
    """Re-render each flagged non-N preset at the suggested step(s) and judge again (same family bias)."""
    out = {}
    for n, row in rows.items():
        if not row["flags"] or not row["suggest"] or row["program"] == "N":
            continue
        m = res[n]
        st2 = {**m["settings"], **{k: v for k, v in row["suggest"].items() if k in m["settings"]}}
        f2 = render_features(m["program"], st2, safe(n) + "~verify", m["program"] in ("M", "L", "E"))
        m2 = {**m, "settings": st2, "render": f2}
        if st2["delay"] != m["settings"]["delay"]:
            m2["fa_render"] = m["fa_render"] + int(round((st2["delay"] - m["settings"]["delay"]) * FS / 1e3))
        r2 = judge(n, m2, bias, res)
        out[n] = {"settings": st2, "flags": r2["flags"],
                  "rt_miss_pct": r2.get("rt_miss_pct"), "lf_miss_pct": r2.get("lf_miss_pct"), "hf_miss_pct": r2.get("hf_miss_pct"),
                  "ampl_rel_db": r2.get("ampl_rel_db"), "ampl_implied": r2.get("ampl_implied"),
                  "onset_diff_ms": r2["onset_used_diff_ms"]}
    return out


def fmt(v, f="{:+.1f}"):
    return "-" if v is None or (isinstance(v, float) and np.isnan(v)) else f.format(v)


def write_outputs():
    res, rows, bias, flagged = stage_report()
    ver = verify(res, rows, bias)
    # numbered-preset systematic AMPL offset per program (family-relative dB)
    num = {}
    for prog in ("S", "M", "L", "E"):
        v = [r["ampl_rel_db"] for n, r in rows.items() if r["program"] == prog and ap.family_of(n) is None and "ampl_rel_db" in r]
        num[prog] = {"n": len(v), "median_db": float(np.median(v)), "min": float(min(v)), "max": float(max(v))}
    for n, r in rows.items():
        if "ampl_rel_db" in r and ap.family_of(n) is None:
            r["ampl_rel_numbered_db"] = r["ampl_rel_db"] - num[r["program"]]["median_db"]
    a2 = {}
    for n, r in rows.items():
        if "ampl_rel_numbered_db" in r and r["settings"]["ampl"] > 0:
            imp = r["settings"]["ampl"] * 10 ** (-r["ampl_rel_numbered_db"] / 20)
            a2[n] = {"implied": float(imp), "set": r["settings"]["ampl"], "off": bool(abs(imp - r["settings"]["ampl"]) > 5)}
    passing = [n for n, r in rows.items() if not r["flags"]]
    cat = {}
    for n, r in rows.items():
        for f in r["flags"]:
            cat.setdefault(f, []).append(n)
    fam_pass = [n for n in rows if ap.family_of(n) and not rows[n]["flags"]]
    summary = {"n_presets": len(rows), "n_pass": len(passing), "pass": passing, "flag_counts": {k: len(v) for k, v in cat.items()},
               "flagged": {n: {"flags": r["flags"], "suggest": r["suggest"], "verify": ver.get(n)} for n, r in rows.items() if r["flags"]},
               "family_chain_bias": {k: {"t60_pct": (np.expm1(v["t60"]) * 100).round(2).tolist(), "ampl_db": v["ampl"]} for k, v in bias.items()},
               "numbered_ampl_offset_db": num, "ampl_after_numbered_offset": a2, "setup": {"rate": 48000, "mod": 1, "sixteen": 0, "latency_removed_samples": 50,
                                                          "t60_renders": "3 input phases (0, 0.37, 0.9 s) for M / L / E"},
               "rows": rows}
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)

    L = ["# Altiverb factory presets rendered through the plugin, against their own Altiverb IR", "",
         "Generated by `tools/altiverb_preset_check.py` (`--stage all`). Every preset rendered with `REV246Render` at its snapped "
         "panel values (`presets.json`), 48 kHz (the IRs' rate), `mod=1`, **`sixteen=0`**, impulse at 0 dBFS, 50 samples of engine latency "
         "removed (the same for every render). `sixteen=1` (the plugin default) puts a -96 dBFS quantisation floor under the tail that "
         "wrecks a Schroeder T60 of an impulse render (T60 reads 1-2 s for a 5 s decay); the 16-bit render of every non-N preset is kept in "
         "`measure.json` (`render16on`) and the check is therefore not valid at `sixteen=1`. M / L / E are time-variant: T60 and EDT are the "
         "mean of 3 renders at different input phases, so they are statistics (jitter in the table), never nulls.", "",
         f"**{len(passing)} of {len(rows)} presets pass** every check (no flag); {len(rows) - len(passing)} flagged. "
         f"Family presets (the calibration): {len(fam_pass)} of {sum(1 for n in rows if ap.family_of(n))} pass.", "",
         "Metrics are judged **family-relative**: per program the median over its family presets of each log ratio (render / IR) is the "
         "'chain bias' and is removed (Altiverb 8 kHz anomaly, its low-band chain, level conventions). Chain bias of T60, 125 Hz-8 kHz, % "
         "(S / M / L / E): " + "; ".join(f"{k} " + " ".join(f"{x:+.1f}" for x in np.expm1(v["t60"]) * 100) for k, v in bias.items()) + ".", "",
         "Flag thresholds (the snapping step): **RT** = mean of 500 Hz-2 kHz T60 miss > half the RT table step's ratio in the direction of the "
         "miss; **LF** = 125 / 250 Hz miss after the RT part > half the LF step's ratio; **HF** = 2 / 4 kHz likewise (8 kHz never counted); "
         "**AMPL** = implied AMPL (set x 10^(-dB/20), family-relative cluster-tap / reverb ratio) differs from the set value by > 5; "
         "**DELAY** = onset difference > 0.5 ms; **LOWPASS** = a candidate in {6, 8, 10, 12} kHz fits the tap-0 shape (1-12 kHz, relative "
         "to the family reference, level removed) better than the set one by > 1 dB rms; **N** = a neighbouring (SPAN, EARLY, LATE) rendered "
         "by the plugin lowers the 5 ms band-envelope residual by > 0.5 dB and 15 %.", ""]
    L += ["## Systematic findings", ""]
    L += [f"- **M and L: the plugin's cluster / reverb ratio is {-bias['M']['ampl']:.1f} dB (M) / {-bias['L']['ampl']:.1f} dB (L) below the "
          f"Altiverb's family presets at the same AMPL**, S and E are within {abs(bias['S']['ampl']):.1f} / {abs(bias['E']['ampl']):.1f} dB. "
          "Not Altiverb-specific: the same measure on Tomkin M files gives the plugin -9 to -12 dB against the hardware capture at the same "
          "settings (M: all four cluster taps 3.75x smaller relative to the 40-120 ms reverb RMS, constant over AMPL 10-99 and RT / LF / HF), "
          "while on Tomkin S it is a constant +2.3 dB. The family-relative method hides it; absolute AMPL readings in M / L are therefore not "
          "comparable to Tomkin's. Not a preset error, a plugin/engine level question (reverb level against cluster in M and L).",
          f"- **Numbered presets carry a common AMPL offset against the family** (render minus IR, dB, median): "
          + ", ".join(f"{k} {v['median_db']:+.1f} (range {v['min']:+.1f}..{v['max']:+.1f}, n={v['n']})" for k, v in num.items()) +
          ". Every numbered S / M / L preset has the plugin's cluster louder than the IR's, relative to the family; E is within about 1 dB. "
          "Either the numbered AMPLs are all high relative to the family anchors (80 / 80 / 60), or the anchors are low. The column "
          "'AMPL vs numbered median' removes it to show the scatter.", ""]
    L += ["## Table", "",
          "Onset = render minus IR (ms). T60 / LF / HF = family-relative miss, plugin vs IR, % (half-step in brackets; + = plugin longer). "
          "AMPL = family-relative cluster / reverb ratio, render minus IR (dB) and the implied AMPL. EDT = mid-band, family-relative. "
          "jit = render-to-render sd of T60 (%), 8k = 8 kHz miss (not counted).", "",
          "| Preset | Prog | Set RT / LF / HF / AMPL / DELAY / LP | onset | RT % (hs) | LF % (hs) | HF % (hs) | 8k % | EDT % | AMPL dB (impl.) | vs numbered | LP | jit | flags |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, r in rows.items():
        st = r["settings"]
        if r["program"] == "N":
            sets = f"SPAN {st['ampl']} / EARLY {st['lf']} / LATE {st['hf']} / {st['delay']} ms / LP {st['lowpass']}"
            L.append(f"| {n} | N | {sets} | {fmt(r['onset_used_diff_ms'], '{:+.2f}')} | - | - | - | - | - | env rms {r['n_cost_db']:.2f} dB (best nbr {r['n_best']['cost']:.2f}) | - | "
                     f"{fmt(r.get('lp', {}).get('current', {}).get('rms'), '{:.1f}')} dB rms | - | {', '.join(r['flags']) or 'ok'} |")
            continue
        sets = f"{st['rt']} / {st['lf']} / {st['hf']} / {st['ampl']} / {st['delay']} / {st['lowpass']}"
        ampl = f"{fmt(r.get('ampl_rel_db'), '{:+.2f}')} ({fmt(r.get('ampl_implied'), '{:.0f}')})" if "ampl_rel_db" in r else "-"
        vn = fmt(r.get("ampl_rel_numbered_db"), "{:+.2f}") if "ampl_rel_numbered_db" in r else "-"
        lp = f"{r['lp']['current']['rms']:.1f}" if "lp" in r else "-"
        L.append(f"| {n} | {r['program']} | {sets} | {fmt(r['onset_used_diff_ms'], '{:+.2f}')} | {fmt(r['rt_miss_pct'])} ({r['rt_half_step_pct']:.1f}) | "
                 f"{fmt(r['lf_miss_pct'])} ({r['lf_half_step_pct']:.1f}) | {fmt(r['hf_miss_pct'])} ({r['hf_half_step_pct']:.1f}) | {fmt(r['rel8k_pct'])} | "
                 f"{fmt(r['edt_mid_pct'])} | {ampl} | {vn} | {lp} | {np.mean(r['t60_sd_pct'][2:5]):.1f} | {', '.join(r['flags']) or 'ok'} |")
    L += ["", "## Flagged presets: miss, suggested step, and the suggestion re-rendered", "",
          "| Preset | Flags | Miss | Suggested | After the suggestion (re-render) |", "|---|---|---|---|---|"]
    for n, r in rows.items():
        if not r["flags"]:
            continue
        miss = []
        if "RT" in r["flags"]: miss.append(f"RT {r['rt_miss_pct']:+.1f} % (half step {r['rt_half_step_pct']:.1f})")
        if "LF" in r["flags"]: miss.append(f"LF {r['lf_miss_pct']:+.1f} % ({r['lf_half_step_pct']:.1f})")
        if "HF" in r["flags"]: miss.append(f"HF {r['hf_miss_pct']:+.1f} % ({r['hf_half_step_pct']:.1f})")
        if "AMPL" in r["flags"]: miss.append(f"AMPL {r['ampl_rel_db']:+.2f} dB, implied {r['ampl_implied']:.0f} vs {r['settings']['ampl']:.0f}")
        if "DELAY" in r["flags"]: miss.append(f"onset {r['onset_used_diff_ms']:+.2f} ms")
        if "LOWPASS" in r["flags"]: miss.append(f"LP rms {r['lp']['current']['rms']:.1f} dB vs {r['lp']['best_rms']:.1f} at {r['lp']['best']}")
        if "N" in r["flags"]: miss.append(f"env rms {r['n_cost_db']:.2f} vs {r['n_best']['cost']:.2f}")
        v = ver.get(n)
        after = "-"
        if v:
            after = (f"flags: {', '.join(v['flags']) or 'none'}; RT {fmt(v['rt_miss_pct'])} %, LF {fmt(v['lf_miss_pct'])} %, HF {fmt(v['hf_miss_pct'])} %, "
                     f"AMPL impl. {fmt(v['ampl_implied'], '{:.0f}')}, onset {fmt(v['onset_diff_ms'], '{:+.2f}')} ms")
        L.append(f"| {n} | {', '.join(r['flags'])} | {'; '.join(miss)} | {r['suggest']} | {after} |")
    L += ["", "## Notes", "",
          "- PRESET23 (LP 8 kHz): the plugin's 2nd-order Butterworth does not match the Altiverb's steep cut-off (tap-0 shape residual "
          f"{rows['PRESET23 SM SECT1']['lp']['current']['rms']:.1f} dB rms at 8, {rows['PRESET23 SM SECT1']['lp']['best_rms']:.1f} at the best step): a plugin filter-shape limit, not a step choice. "
          "Its render also no longer matches the N front-tap template, so its onset comes from the first arrival (+0.3 ms: the Butterworth's group delay).",
          "- AMPL with the numbered-preset offset removed (per program median): " + str(sum(v["off"] for v in a2.values())) + " of " + str(len(a2)) +
          " numbered presets still differ from the set AMPL by more than 5: " + ", ".join(f"{n} ({v['set']:.0f} -> {v['implied']:.0f})" for n, v in a2.items() if v["off"]) + ".",
          "- PRESET23 (N): the broadband envelope residual is 2.6 dB rms (the others 0.16-0.17): the plugin's Butterworth LP 8 against the Altiverb's steeper cut-off, no neighbouring SPAN / EARLY / LATE improves it.",
          "- PRESET90 has no cluster (AMPL 0): AMPL and LOW PASS cannot be checked from tap 0.",
          "- `sixteen=1` (default): see the header; the T60 of the 16-bit render is in `measure.json`.",
          "- Time variance: M / L / E T60 are 3-render means; `jit` is the sd across renders, to compare a miss with."]
    open(os.path.join(OUT, "report.md"), "w").write("\n".join(L) + "\n")
    return summary


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--stage", default="all")
    ap_.add_argument("--only", nargs="*")
    a = ap_.parse_args()
    os.makedirs(RENDERS, exist_ok=True)
    if a.stage in ("measure", "all"):
        stage_measure(a.only)
    if a.stage in ("fa", "all"):
        stage_fa()
    if a.stage in ("report", "all"):
        write_outputs()


if __name__ == "__main__":
    main()
