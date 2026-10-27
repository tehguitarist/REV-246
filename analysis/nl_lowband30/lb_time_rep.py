import pickle, numpy as np, collections
D = pickle.load(open("time.pkl", "rb"))
G = collections.defaultdict(list)
for (sp, fs, k), rec in D.items(): G[(sp, fs)].append(rec)
bins = [(-40, 0), (0, 40), (40, 80), (80, 160), (160, 320), (320, 450), (450, 520), (520, 680), (680, 800)]  # samples re z (front frame zero)
print("rms of 100Hz-1kHz residual, dB re model peak, per time bin (samples re front zero; taps at 41, 410, 463; first main 685)")
for c in (0, 1):
    print(f"\n channel {'LR'[c]}   bins: " + " ".join(f"{a:>4d}..{b:<4d}" for a, b in bins))
    for key in sorted(G):
        row = []
        for a, b in bins:
            v = []
            for rec in G[key]:
                z = rec[c]["z"]; eb = rec[c]["eb"][z + a:z + b]
                v.append(np.mean(eb ** 2))
            row.append(10 * np.log10(np.median(v)))
        print(f" {key[0]}@{key[1]:5.1f} n={len(G[key]):2d} " + " ".join(f"{x:9.1f}" for x in row))
