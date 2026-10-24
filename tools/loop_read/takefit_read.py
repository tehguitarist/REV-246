"""What §31's per-setting refits (lowband_shelf.py --stage takefit) are made of (model_s_results 32).

    EMT_KERNEL=chain3 EMT_PMAX=0.8475 EMT_HFRT_MIN=0.12 EMT_CV_FOLDS=$PWD/analysis/loop_cv/cap0.8475_folds \
        SHELF_CELLS=$PWD/analysis/loop_cv/hfcells_chain3_cap0.8475_hfrt/cells.json \
        .venv/bin/python tools/loop_read/takefit_read.py            # analysis/lowband/chain3_takefit_read.{log,json}

1. Each loop's correction (broadband gain, DC section, HF pole: 12 numbers per setting, fitted per take) as a per-pass
   magnitude change in dB at 125 Hz ... 8 kHz, and the unit's absolute per-pass DC and Nyquist gains (baseline + correction).
2. Transfer: take B's values predicted from take A's at the other settings of a cell ((LF, RT), (HF, RT), with or
   without the loop), and by an LF effect + HF effect per (RT step, loop) (additive), each setting left out.
3. The DC floor per loop: law DC gain (global law ratios) -> unit DC gain, LF 0.1 / 0.4.
4. The additive effects at RT >= 1.0 per loop (1 and 4 kHz).
No library data are loaded (laws, fold terms and the takefit files only); seconds, ~0.2 GB.
"""
import glob, json, os, sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lowband_shelf as ls, q2_injection as q2, q1_generalise as q1, model_s as m  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(HERE))
TF = ls.TAKEFIT_DIR
OUT = os.path.join(ROOT, "analysis", "lowband", "chain3_takefit_read")
F = np.array([125, 250, 500, 1000, 2000, 4000, 8000.])
LFS, HFS = [0.1, 0.4, 1.0, 2.0], [0.2, 0.4, 0.6, 1.0]
LAW_DC = {0.1: 0.1317, 0.4: 0.414}                # the global law's DC ratios (LF 0.4: §25.4)


def corrections():
    """Per setting, loop and take: dB change per pass at F, and the unit's DC / Nyquist per-pass gains (HF 1.0 form)."""
    _, law, _ = q2.final_injection(compact=ls.lb.lf_.COMPACT)
    w = np.r_[2 * np.pi * F / q1.FS, 0.0, np.pi]
    a, z1 = law["a"], np.exp(-1j * w)
    op = lambda p: (1 - p) / (1 - p * z1)
    rows = []
    for f in sorted(glob.glob(os.path.join(TF, "*.json"))):
        d = json.load(open(f))
        if "A" not in d["fit"] or "B" not in d["fit"]:
            continue
        t, lfdc, kw = ls.base_kw(d["key"], ls.HF_CELLS_JSON)
        L = json.loads(json.dumps(law))
        L["lf"] = {float(k): v for k, v in L["lf"].items()}
        for k_, x in lfdc.items():
            L["lf"][float(k_)]["dc"] = float(x)
        lab = dict(lf=d["lf"], rt=d["rt"], hf=d["hf"])
        args = dict(lfrt=t["lfrt"], gamma=t["gamma"], hfdp=t["hfdp"], pramp=t["pramp"], b230=0.0, **kw)
        G1 = q2.final_Gs(dict(lab, hf=1.0), L, w, **args)
        pb = None
        if d["hf"] < 1:
            G = q2.final_Gs(lab, L, w, **args)
            pb = [(1 - h) / (1 + h) for h in [float(np.real(g[-1] / g1[-1])) for g, g1 in zip(G, G1)]]
        for i in range(4):
            r = dict(key=d["key"], **lab, loop=i)
            for tk in "AB":
                u = np.r_[d["fit"][tk], np.zeros(12 - len(d["fit"][tk]))]
                dB = u[i] + 20 * np.log10(np.abs(m.shelf(w, 10 ** (u[4 + i] / 20), 1.0, a)))
                if pb is not None:
                    dB = dB + 20 * np.log10(np.abs(op(pb[i] + u[8 + i]) / op(pb[i])))
                r[tk] = dB[:len(F)].tolist()
                r["dc_" + tk] = float(abs(G1[i][-2]) * 10 ** (dB[-2] / 20))
                r["ny_" + tk] = float(abs(G1[i][-1]) * 10 ** (u[i] / 20))     # broadband part only (the pole is HF)
            rows.append(r)
    return rows


def explained(B, P):
    return [1 - np.var(B[:, j] - P[:, j]) / np.var(B[:, j]) for j in range(B.shape[1])]


def transfer(rows):
    out = []
    A, B = np.array([r["A"] for r in rows]), np.array([r["B"] for r in rows])
    setting = lambda r: (r["lf"], r["rt"], r["hf"])
    out.append("## 2. Transfer of the per-pass corrections (take B predicted from take A at other settings; explained variance)")
    out.append("   freq (Hz)          " + " ".join(f"{f:>7.0f}" for f in F))
    out.append("   rms B (dB/pass)    " + " ".join(f"{x:7.4f}" for x in np.sqrt((B ** 2).mean(0))))
    out.append("   corr(A, B)         " + " ".join(f"{np.corrcoef(A[:, j], B[:, j])[0, 1]:7.2f}" for j in range(len(F))))
    for cell in [("lf", "rt"), ("hf", "rt"), ("lf", "rt", "loop"), ("hf", "rt", "loop")]:
        P = np.zeros_like(B)
        for i, r in enumerate(rows):
            o = [j for j, s in enumerate(rows) if all(s[c] == r[c] for c in cell) and setting(s) != setting(r)]
            P[i] = A[o].mean(0) if o else 0
        out.append(f"   cell {'+'.join(cell):14s}" + " ".join(f"{x:7.2f}" for x in explained(B, P)))
    P = additive(rows, A)
    ok = ~np.isnan(P[:, 0])
    out.append(f"   additive LF + HF   " + " ".join(f"{x:7.2f}" for x in explained(B[ok], P[ok])) + "   (per (RT, loop))")
    for rt in sorted({r["rt"] for r in rows}):
        k = ok & np.array([r["rt"] == rt for r in rows])
        out.append(f"     RT {rt:<4}          " + " ".join(f"{x:7.2f}" for x in explained(B[k], P[k])) +
                   f"   rms 1k / 4k {np.sqrt(np.mean(B[k, 3] ** 2)):.4f} / {np.sqrt(np.mean(B[k, 5] ** 2)):.4f}")
    return out, P


def additive(rows, A):
    """Each row predicted from the other settings of its (RT, loop) block as an LF effect + an HF effect."""
    P = np.full_like(A, np.nan)
    for i, r in enumerate(rows):
        blk = [j for j, s in enumerate(rows) if s["rt"] == r["rt"] and s["loop"] == r["loop"] and (s["lf"], s["hf"]) != (r["lf"], r["hf"])]
        if not any(rows[j]["lf"] == r["lf"] for j in blk) or not any(rows[j]["hf"] == r["hf"] for j in blk):
            continue
        oh = lambda s: np.r_[[s["lf"] == x for x in LFS], [s["hf"] == x for x in HFS]].astype(float)
        X = np.array([oh(rows[j]) for j in blk])
        P[i] = oh(r) @ np.linalg.lstsq(X, A[blk], rcond=None)[0]
    return P


def knee(rows):
    out = ["## 3. The DC floor per loop (LF 0.1 / 0.4, sorted by the law's DC loss): law DC gain -> unit (A, B mean); dB per pass"]
    c = defaultdict(list)
    for r in rows:
        if r["lf"] in LAW_DC:
            c[(r["lf"], r["rt"], r["loop"])].append(r)
    lines = []
    for (lf, rt, l), u in c.items():
        x = 10 ** (-3 * q1.P[l] / (q1.FS * LAW_DC[lf] * rt))
        y = np.mean([(r["dc_A"] + r["dc_B"]) / 2 for r in u])
        s = np.std([(r["dc_A"] + r["dc_B"]) / 2 for r in u])
        lines.append((x, f"   LF {lf} RT {rt} loop {l + 1}: law {x:.4f} ({-20 * np.log10(x):6.2f} dB) -> unit {y:.4f} ({-20 * np.log10(y):6.3f} dB), "
                         f"sd over HF {s:.4f}"))
    out += [s for x, s in sorted(lines) if x < 0.36]
    return out


def effects(rows):
    out = ["## 4. Additive effects at RT >= 1.0 (both takes, all settings of the block), dB per pass, loops 1-4, at 1 kHz | 4 kHz"]
    Y = np.array([(np.array(r["A"]) + np.array(r["B"])) / 2 for r in rows])
    for rt in (1.0, 1.4, 1.8, 2.2, 2.5):
        eff = np.zeros((4, 7, len(F)))
        for l in range(4):
            k = [i for i, r in enumerate(rows) if r["rt"] == rt and r["loop"] == l]
            X = np.array([np.r_[[rows[i]["lf"] == x for x in LFS], [rows[i]["hf"] == x for x in HFS[:3]]] for i in k], float)
            eff[l] = np.linalg.lstsq(X, Y[k], rcond=None)[0]
        out.append(f"   RT {rt}")
        for e, nm in enumerate([f"LF {x}" for x in LFS] + [f"HF {x}" for x in HFS[:3]]):
            out.append(f"     {nm:6s} " + " ".join(f"{x:+.4f}" for x in eff[:, e, 3]) + "  | " + " ".join(f"{x:+.4f}" for x in eff[:, e, 5]))
    return out


def main():
    rows = corrections()
    tr, _ = transfer(rows)
    out = [f"takefit files: {TF} ({len(rows) // 4} settings)", ""] + tr + [""] + knee(rows) + [""] + effects(rows)
    json.dump(rows, open(OUT + ".json", "w"))
    open(OUT + ".log", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
