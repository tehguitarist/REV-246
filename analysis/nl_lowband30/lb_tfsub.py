exec(open("lb_tfcorr.py").read().split("def report")[0])
fs0 = [k for k, d in S.items() if abs(d["fs"]) < .5]; fs23 = [k for k in S if k not in fs0]
for lab, tr, te in (("fs23->fs0", fs23, fs0), ("fs0->fs23", fs0, fs23)):
    print(lab, "sub-band null before -> after (held-out T, 40 log bands x 2 numbers)")
    for lo, hi in ((20, 25), (25, 32), (32, 50), (50, 100), (100, 300), (300, 1000)):
        a = []; b = []
        for c in (0, 1):
            Tf = est(tr, c)
            for k in te: a.append(null(k, c, None, lo, hi)); b.append(null(k, c, Tf, lo, hi))
        print(f"   {lo:4d}-{hi:<4d} {np.median(a):6.1f} -> {np.median(b):6.1f}")
