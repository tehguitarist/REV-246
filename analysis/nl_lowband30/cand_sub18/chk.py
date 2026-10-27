from common import *
import pickle, json
D = pickle.load(open(os.path.join(HERE, "smoke2.pkl"), "rb"))
for k, d in D.items():
    b = json.load(open(os.path.join(BENCH, "_".join(map(str, k)) + ".json")))
    g = d["g"]["band"]
    mine = [10*np.log10(g[("20-100",c)]["ee_raw"]/g[("20-100",c)]["yy"]) for c in (0,1)]
    cor = [10*np.log10(g[("20-100",c)]["ee_cor"]/g[("20-100",c)]["yy"]) for c in (0,1)]
    print(k, d["sess"], "mine raw", np.round(mine,1), "cor", np.round(cor,1), "bench", np.round(b["bands"]["20-100"],1), "c 1e-4", np.round(np.array(d["g"]["c"])*1e4,1))
