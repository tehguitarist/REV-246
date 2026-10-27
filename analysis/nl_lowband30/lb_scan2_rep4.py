import json, sys, numpy as np, collections
sys.path.insert(0, "/Users/leighpierce/Git/246/tools")
import os; os.environ["NL_BK"] = "0"
import nl_model as M
d = json.load(open("scan2.json"))
rows = []
for k, v in d.items():
    s, e, l = (float(x) for x in k.split("_"))
    ld = M.late_db(s, e) if l > 0 else -999
    rows.append((s, e, l, ld, v))
def b_late(ld): return "none" if ld == -999 else ("quiet<-40" if ld < -40 else ("mid -40..-15" if ld < -15 else "loud>=-15"))
def b_L(l): return "LATE<=0.6" if l <= 0.6 else ("0.8-1.2" if l <= 1.2 else ">=1.4")
g = collections.defaultdict(list)
for s, e, l, ld, v in rows:
    sess = "fs0" if abs(v["fs"]) < 0.5 else "fs>0"
    g[(b_late(ld), b_L(l) if ld != -999 else "-", sess)].append(0.5 * (v["w30_100-1k_0"] + v["w30_100-1k_1"]))
print("100Hz-1kHz first-30ms null (median dB, L/R mean), by late level (model law), LATE, session")
for key in sorted(g):
    print(f"  late {key[0]:13s} {key[1]:10s} {key[2]:5s} n {len(g[key]):3d} median {np.median(g[key]):6.1f}")
