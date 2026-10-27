"""Does a per-EARLY-class transfer beat the pooled one?  N leave-one-SPAN-out from joint_N_sums.pkl (SM 0): T_joint(N other SPANs all EARLY + S) vs T estimated with the SAME-class N captures only (class A: EARLY <= 1.0, class B: EARLY >= 1.5; plus S in both), scored on the held-out SPAN's captures of that class."""
import sys
sys.argv = [sys.argv[0], "0", "1"]
exec(open("lb_joint.py").read().split("one = np.ones")[0])
one = np.ones((2, TNB), complex); Sall = tot([S[q] for q in S_keys()])
cls = lambda k: "A" if k[1] <= 1.0 else "B"
print("N captures by class:", {c: sum(1 for k in N if cls(k) == c) for c in "AB"})
res = {}
for sp in (0.2, 0.3, 0.4, 0.5):
    tr = {k: d for k, d in N.items() if d["span"] != sp}
    T_pool = T_from([(*Sall, 1.0), (*tot(list(tr.values())), 1.0)])
    T_cls = {c: T_from([(*Sall, 1.0), (*tot([d for k, d in tr.items() if cls(k) == c]), 1.0)]) for c in "AB"}
    T_clsN = {c: T_from([(*tot([d for k, d in tr.items() if cls(k) == c]), 1.0)]) for c in "AB"}
    for k, d in N.items():
        if d["span"] != sp: continue
        for c in (0, 1):
            for bi, (lo, hi) in enumerate(BD[:4]):
                m = (fc >= lo) & (fc < hi); den = d["syy"][c][m].sum()
                for nm, T in (("none", one), ("pooled", T_pool), ("class+S", T_cls[cls(k)]), ("classN", T_clsN[cls(k)])):
                    res.setdefault((cls(k), nm, bi), []).append(10 * np.log10(max(err(d, T, m, c), 1e-30) / den))
for c in "AB":
    print(f"class {c} ({'EARLY<=1' if c=='A' else 'EARLY>=1.5'}): median held-out band error (L+R pooled), dB:   bands {[b for b in BD[:4]]}")
    for nm in ("none", "pooled", "class+S", "classN"):
        print(f"   {nm:8s} " + "  ".join(f"{np.median(res[(c, nm, bi)]):6.1f}" for bi in range(4)))
