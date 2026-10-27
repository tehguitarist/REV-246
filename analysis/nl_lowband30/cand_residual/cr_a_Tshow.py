from cr_common import *
exec(open(os.path.join(HERE, "cr_a_lti2.py")).read().split("print(\"EARLY |")[0].replace("from cr_common import *", ""))
EAR2 = (1.0, 1.5, 2.0, 3.0, 10.0)
Ts = {e_: estimate([k for k in sel if k[1] == e_]) for e_ in EAR2}
cf = np.sqrt(edges[:-1] * edges[1:])
for c in (0, 1):
    print(f"\nchannel {'LR'[c]}: T-1 as amplitude %/phase mrad (the model * T = capture), per EARLY (1.0 1.5 2.0 | 3.0 10.0 are the other session)")
    for i in range(ng):
        if cf[i] < 19 or cf[i] > 1000: continue
        print(f"{cf[i]:7.1f}: " + "  ".join(f"{100*(abs(Ts[e][c][i])-1):+6.2f}%/{1000*np.angle(Ts[e][c][i]):+6.1f}" for e in EAR2))
