import json, sys, numpy as np
for tag in sys.argv[1:] or ("lb0", "lb2"):
    d = json.load(open(f"/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline/p6_{tag}.json"))
    print("=====", tag, "captures", len(d))
    chk = max(abs(v["base"]["0.0-0.05"][0] - v["bench_json_0.0-0.05"][0]) for v in d.values())
    chk2 = max(abs(v["base"]["20-100@0.05-0.6"][c] - v["bench_json_20-100"][c]) for v in d.values() for c in (0, 1))
    print(f"baseline reproduces bench json: max |diff| 0.0-0.05 {chk:.3f} dB, 20-100 band {chk2:.3f} dB")
    for w in ("0.0-0.05", "0.05-0.15", "0.15-0.3", "20-100@0.0-0.05", "20-100@0.05-0.6", "100-1000@0.0-0.05", "100-1000@0.05-0.6"):
        b = np.array([v["base"][w] for v in d.values() if w in v["base"]]).ravel()
        f = np.array([v["dcfix"][w] for v in d.values() if w in v["dcfix"]]).ravel()
        dl = f - b
        if "dctaper" in next(iter(d.values())):
            g = np.array([v["dctaper"][w] for v in d.values() if w in v["dctaper"]]).ravel(); dg = g - b
            print(f"{'':18s}   taper: median {np.median(g):7.2f} median delta {np.median(dg):6.2f} better>0.5 {int((dg<-0.5).sum()):2d} worse>0.5 {int((dg>0.5).sum()):2d} worst {dg.max():5.2f}")
        print(f"{w:18s} n{len(b):3d} median base {np.median(b):7.2f} -> fix {np.median(f):7.2f}  median delta {np.median(dl):6.2f}  mean {dl.mean():6.2f}  better>0.5: {int((dl<-0.5).sum()):2d}  worse>0.5: {int((dl>0.5).sum()):2d} worst {dl.max():5.2f} best {dl.min():6.2f}")
