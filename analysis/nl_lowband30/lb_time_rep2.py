import pickle, numpy as np, collections
D = pickle.load(open("time.pkl", "rb"))
G = collections.defaultdict(list)
for (sp, fs, k), rec in D.items(): G[(sp, fs)].append(rec)
edges = list(range(0, 130, 10)) + [150, 200, 300]
for name, field in (("100Hz-1kHz residual", "eb"), ("full-band residual", "e"), ("model (100-1k)", "mb"), ("capture y full", "y")):
    for c in (0, 1):
        print(f"\n{name}, ch {'LR'[c]}: rms dB re model peak by ABSOLUTE index (file sample 0 = index 50; first front tap at 50+fs+41)  bins: {edges}")
        for key in sorted(G):
            row = []
            for a, b in zip(edges[:-1], edges[1:]):
                v = [np.mean(rec[c][field][a:b] ** 2) for rec in G[key]]
                row.append(10 * np.log10(np.median(v) + 1e-30))
            print(f" {key[0]}@{key[1]:5.1f} " + " ".join(f"{x:6.1f}" for x in row))
