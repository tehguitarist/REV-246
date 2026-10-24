"""loop_fir.py --current fits, held out by capture (AMPL 0 + 60): R excess over the repeat null per (HF, RT) cell,
baseline -> + free correction; and for HF 0.2 short RT, where in frequency the 150-300 ms residual sits (dense_mid prep,
next AMPL capture as the repeat). model_s_results §25.14.

    EMT_KERNEL=chain3 .venv/bin/python tools/loop_read/heldout.py"""
import collections, glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import butter, sosfiltfilt
import dense_mid as dm
rep = {(r["key"], r["amp"]): r for r in json.load(open(os.path.join(ROOT, "analysis/q2/nulls_free_chain2.json")))}
D = [json.load(open(p)) for p in glob.glob(os.path.join(ROOT, "analysis/loop_fir/chain3_current_kc8/*.json"))]
ws = ("50-150", "150-300", "300-600")
def ex(rr, w):
    e = [r["null_comp"][w][1] - rep[(r["key"], r["amp"])]["repeat_comp"][w][1] for r in rr
         if r["null_comp"].get(w) and rep.get((r["key"], r["amp"])) and rep[(r["key"], r["amp"])]["repeat_comp"].get(w)]
    return np.median(e) if e else np.nan
G = collections.defaultdict(list)
for d in D:
    G[(d["hf"], d["rt"])].append(d)
print("held out by capture, R excess over repeat (chain2 refs): baseline -> + free correction (+-8 lags per loop)\nHF  RT  | b    | 50-150        150-300        300-600")
for k in sorted(G):
    dd = G[k]
    L = [r for d in dd for r in d["rows"] if r["loops"] == "law" and r["held_out"]]
    F = [r for d in dd for r in d["rows"] if r["loops"] == "fir" and r["held_out"]]
    print(f"{k[0]:<3} {k[1]:<3} | {dd[0]['b230']:.2f} | " + "  ".join(f"{ex(L, w):+5.1f}->{ex(F, w):+5.1f}" for w in ws))
FS, PRE = dm.FS, dm.PRE
bands = [(88, 177), (177, 354), (354, 707), (707, 1414), (1414, 2828), (2828, 5657), (5657, 8000), (8000, 10000), (10000, 11300)]
sos = {b: butter(3, b, btype="band", fs=FS, output="sos") for b in bands}
for keys, nm in (([f"{lf}_{rt}_0.2" for lf in ("1.0", "2.0") for rt in ("0.6", "0.8", "1.0")], "HF 0.2, RT 0.6-1.0, LF 1.0 / 2.0"),
                 ([f"{lf}_{rt}_0.2" for lf in ("0.1", "0.4") for rt in ("0.6", "0.8", "1.0")], "HF 0.2, RT 0.6-1.0, LF 0.1 / 0.4"),
                 ([f"1.0_{rt}_1.0" for rt in ("0.6", "0.8", "1.0")], "HF 1.0 (control), RT 0.6-1.0, LF 1.0")):
    acc = {b: [[], []] for b in bands}
    for k in keys:
        z = np.load(os.path.join(dm.OUT, "prep", k + ".npz"))
        caps = sorted(int(n[3:]) for n in z.files if n.startswith("res"))
        for j in caps:
            y, m = z[f"y{j}"][:, 1].astype(float), z[f"m{j}"][:, 1].astype(float)
            yb = z[f"y{caps[(caps.index(j) + 1) % len(caps)]}"][:, 1].astype(float)
            n = min(len(y), len(yb)); i0, i1 = PRE + int(.15 * FS), min(PRE + int(.3 * FS), n)
            for b in bands:
                Y, M, YB = (sosfiltfilt(sos[b], v[:n])[i0:i1] for v in (y, m, yb))
                acc[b][0].append(np.sum((Y - M) ** 2)); acc[b][1].append(np.sum((Y - YB) ** 2))
    tot = sum(np.sum(acc[b][0]) for b in bands)
    print(f"\n{nm}: 150-300 ms, R; per band: share of the model's residual energy | model / repeat residual (dB)")
    print("   " + "  ".join(f"{b[0]}-{b[1]}: {100 * np.sum(acc[b][0]) / tot:.0f}% {10 * np.log10(np.median(np.array(acc[b][0]) / np.array(acc[b][1]))):+.1f}" for b in bands))
