# theta depth re-check: the Python sweep sim (== the plugin's core, checked to float level static and -48...-62 dB modulated)
# with N evenly spaced start phases (no random-draw luck), at several depths, against Tomkin's AMPL pairs per group.
import sys, os, json, numpy as np
os.environ.setdefault("EMT_KERNEL", "chain3")
sys.path.insert(0, "tools")
import ml_variance as mv
from concurrent.futures import ProcessPoolExecutor
NPH = 12
DEPTHS = [0.11, 0.12, 0.13, 0.14, 0.15]
WINS = ("recirc1", "recirc2", "0.15-0.3")
def job(a):
    prog, key, rms = a
    M = mv.mll.Model(prog)
    A = rms * np.sqrt(2)
    irs = [mv.sweep_measure(M, key, lambda t, ph=2 * np.pi * k / NPH: A * np.sin(2 * np.pi * mv.THETA_FM * t + ph)) for k in range(NPH)]
    pairs = [mv.pair_nulls(prog, irs[i], irs[j], len(irs[i]) - 1000) for i in range(NPH) for j in range(i + 1, NPH)]
    return prog, key, rms, {w: float(np.median([p[w] for p in pairs if p[w] is not None])) for w in WINS}
def tomkin(prog, key):
    fit = mv.mll.fitted(prog, key); amps = list(fit); out = []
    for i in range(len(amps)):
        for j in range(i + 1, len(amps)):
            e = min(min(fit[amps[i]][2]), min(fit[amps[j]][2]))
            out.append(mv.pair_nulls(prog, fit[amps[i]][0], fit[amps[j]][0], e))
    return out
if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    groups = [(p, k) for p, ks in mv.SIM_KEYS.items() for k in ks]
    if only: groups = groups[:1]; DEPTHS[:] = [0.15]
    jobs = [(p, k, r) for p, k in groups for r in DEPTHS]
    res = {}
    with ProcessPoolExecutor(int(os.environ.get("W", "4")), max_tasks_per_child=1) as ex:
        for p, k, r, v in ex.map(job, jobs):
            res.setdefault(f"{p}_{k}", {})[str(r)] = v; print(p, k, r, v, flush=True)
    out = {}
    for p, k in groups:
        T = tomkin(p, k)
        t = {w: float(np.median([x[w] for x in T if x[w] is not None])) for w in WINS}
        out[f"{p}_{k}"] = dict(tomkin=t, n_tomkin_pairs=len(T), sim=res[f"{p}_{k}"])
    gaps = {r: [] for r in map(str, DEPTHS)}
    for g, o in out.items():
        for r in gaps:
            for w in WINS:
                gaps[r].append(o["tomkin"][w] - o["sim"][r][w])   # + : Tomkin less decorrelated than the sim
    out["gap_median"] = {r: float(np.median(v)) for r, v in gaps.items()}
    print(json.dumps(out["gap_median"]))
    if not only:
        json.dump(out, open("analysis/plugin_check/theta_depth.json", "w"), indent=1)
