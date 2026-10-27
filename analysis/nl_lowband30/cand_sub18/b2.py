from b1 import *
import b1
def longwin(k, v, lo): 
    w = v["win"].get("w45f"); return w is not None and w[3] >= lo * 26373.6
b1.DES["+LATE>0.2, windows >= 0.30 s only"] = lambda k, v: "w45f" if k[0] >= 0.2 and k[1] >= 0.3 and longwin(k, v, 0.30) else None
b1.DES["+LATE>0.2, windows >= 0.30 s, window to 0.6"] = lambda k, v: "w60f" if k[0] >= 0.2 and k[1] >= 0.3 and (k[2] <= 0.2 or longwin(k, v, 0.30)) else None
b1.DES["base, windows >= 0.30 s only"] = lambda k, v: "w45f" if k[0] >= 0.2 and k[2] <= 0.2 and k[1] >= 0.3 and longwin(k, v, 0.30) else None
b1.DES["base + EARLY<0 + SPAN0.1 (all LATE<=0.2)"] = lambda k, v: "w45f" if k[0] >= 0.1 and k[2] <= 0.2 and (k[1] >= 0.3 or enough(v)) else None
base_res, nb = evaluate(list(b1.DES)[0]); ref = row(base_res); pr("base", ref, nb, None)
for des in list(b1.DES)[-4:]:
    r, n = evaluate(des); pr(des, row(r), n, ref)
# also report window-length mix of the base training set
import numpy as np
print("window lengths (s) of base design:", np.round(sorted(NN[k]["win"]["w45f"][3] / 26373.6 for k in NN if DES[list(DES)[0]](k, NN[k])), 2)[::5])
