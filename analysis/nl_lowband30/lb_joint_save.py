import sys, json
import os; sys.argv = [sys.argv[0], os.environ.get("LB_SM", "2"), "1"]
exec(open("lb_joint.py").read().split("one = np.ones")[0])
Sall = tot([S[q] for q in S_keys()]); Nall = tot(list(N.values()))
T = {"joint": T_from([(*Sall, 1.0), (*Nall, 1.0)]), "N": T_from([(*Nall, 1.0)]), "S": T_from([(*Sall, 1.0)])}
m = (fc >= 12) & (fc <= 400)
out = {"f_hz": fc[m].tolist(), "note": "complex data/model transfer ratio on 2 Hz bins (centres), smoothed +-2 bins; joint = N and S with equal program weight (each normalised by its 20-1600 Hz |M|^2); N = SPAN>=0.2 LATE<=0.2 EARLY>=0.3 current model (warp2, branch kernel 0.7); S = LF>=1.0 all RT, _apml baseline"}
for nm, t in T.items():
    for c in (0, 1):
        out[f"{nm}_db{c}"] = (20 * np.log10(abs(t[c][m]))).tolist(); out[f"{nm}_deg{c}"] = np.degrees(np.angle(t[c][m])).tolist()
json.dump(out, open("lowband_joint_T.json", "w"))
print("saved", len(out["f_hz"]), "bins")
