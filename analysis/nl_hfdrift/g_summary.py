import sys, pickle, numpy as np, collections
sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_hfdrift")
rs = pickle.load(open("/Users/leighpierce/Git/246/analysis/nl_hfdrift/g_oracle_30.pkl", "rb"))
def cl(r):
    k = r["key"]
    if k[0] == 0.2 and k[1] == -0.2: return "fam(0.2,-0.2,*)"
    if r["trimmed"]: return "trimmed"
    if k[0] == 0.3: return "SPAN0.3"
    if k[0] == 0.1: return "SPAN0.1"
    return "rest"
m = lambda r, k: float(np.mean(r[k]))
for r in rs: r["fail"] = m(r, "used") - m(r, "none") > 0.5
print(len(rs), "captures; null = mean of L,R, 10-12 kHz, dB")
def table(sel, title):
    print("\n" + title)
    print("%-17s %-5s %3s | %-24s | oracle  free3  free9 | med|su-sb| med|sown-sb| med sb  med su | oracle<none-1 | noise-null  HF share dB" % ("cluster", "fail", "n", "none   used"))
    g = collections.defaultdict(list)
    for r in rs:
        g[(cl(r), r["fail"])].append(r)
    for (c, f), v in sorted(g.items()):
        so = [r for r in v if "s_own" in r]
        print("%-17s %-5s %3d | %6.1f %6.1f         | %6.1f %6.1f %6.1f | %6.3f    %6s    %6.3f %6.3f | %3d%% | %6.1f   %6.1f" % (
            c, f, len(v), np.median([m(r, "none") for r in v]), np.median([m(r, "used") for r in v]),
            np.median([m(r, "oracle") for r in v]), np.median([m(r, "free3") for r in v]), np.median([m(r, "free9") for r in v]),
            np.median([abs(r["s_used"] - r["s_best"]) for r in v]), ("%.3f" % np.median([abs(r["s_own"] - r["s_best"]) for r in so])) if so else "  -",
            np.median([r["s_best"] for r in v]), np.median([r["s_used"] for r in v]),
            100 * np.mean([m(r, "oracle") < m(r, "none") - 1 for r in v]), np.median([np.mean(r["noise_null"]) for r in v]), np.median([r["frac_hf_energy_db"] for r in v])))
table(rs, "A. by cluster x (used worse than none by > 0.5 dB)")
# overall
print("\nB. by cluster (all captures of the cluster in the sample): none / used / oracle / free3 / free9, fail rate with s_used, fail rate with oracle (worse than none)")
g = collections.defaultdict(list)
for r in rs: g[cl(r)].append(r)
for c, v in sorted(g.items()):
    print("%-17s n %3d  %6.1f %6.1f %6.1f %6.1f %6.1f | used worse>0.5: %3d%% oracle worse>0.5: %3d%% | s_best sd %.3f  s_used sd %.3f  corr %.2f" % (c, len(v),
        *[np.median([m(r, k) for r in v]) for k in ("none", "used", "oracle", "free3", "free9")],
        100 * np.mean([m(r, "used") > m(r, "none") + 0.5 for r in v]), 100 * np.mean([m(r, "oracle") > m(r, "none") + 0.5 for r in v]),
        np.std([r["s_best"] for r in v]), np.std([r["s_used"] for r in v]), np.corrcoef([r["s_best"] for r in v], [r["s_used"] for r in v])[0, 1]))
# Q1: failing captures: is the s_best near zero / is the oracle better than none?
f = [r for r in rs if r["fail"]]
print("\nC. the failing captures (n=%d): oracle - none median %.1f dB; free3 - none %.1f; free9 - none %.1f; median none %.1f; |s_used - s_best| median %.3f vs non-failing %.3f; s_best median %.3f, sd %.3f"
      % (len(f), np.median([m(r, "oracle") - m(r, "none") for r in f]), np.median([m(r, "free3") - m(r, "none") for r in f]), np.median([m(r, "free9") - m(r, "none") for r in f]),
         np.median([m(r, "none") for r in f]), np.median([abs(r["s_used"] - r["s_best"]) for r in f]), np.median([abs(r["s_used"] - r["s_best"]) for r in rs if not r["fail"]]),
         np.median([r["s_best"] for r in f]), np.std([r["s_best"] for r in f])))
# Q1b: failing captures split by what the oracle achieves
print("   failing, by 'none' null: already good (< -45 dB):", sum(m(r, "none") < -45 for r in f), "of", len(f), "; the oracle can still help (>1 dB) in", sum(m(r, "oracle") < m(r, "none") - 1 for r in f))
bad = [r for r in f if abs(r["s_used"] - r["s_best"]) > 0.03]
print("   failing with |s_used - s_best| > 0.03: %d of %d ; with s_best within 0.03 of s_used: %d (-> shape mismatch / level, not estimator)" % (len(bad), len(f), len(f) - len(bad)))
for r in f[:0]: pass
# Q3 smoothness
print("\nD. smoothness of s_best in time: sampled captures with timestamps, nearest sampled neighbour in time")
ts = [(r["t"], r["s_best"], r["s_used"], r["key"]) for r in rs if r["t"] is not None]
ts.sort()
T = np.array([x[0] for x in ts]); S = np.array([x[1] for x in ts]); SU = np.array([x[2] for x in ts])
dt, ds, dsu = [], [], []
for i in range(len(ts)):
    j = [k for k in (i - 1, i + 1) if 0 <= k < len(ts)]
    k = min(j, key=lambda q: abs(T[q] - T[i]))
    dt.append(abs(T[k] - T[i])); ds.append(abs(S[k] - S[i])); dsu.append(abs(SU[k] - SU[i]))
dt = np.array(dt); ds = np.array(ds); dsu = np.array(dsu)
for lo, hi in ((0, 120), (120, 600), (600, 3600), (3600, 4e4), (4e4, 1e9)):
    s = (dt >= lo) & (dt < hi)
    if s.sum(): print("   dt %6.0f-%7.0f s  n %3d  median |ds_best| %.3f   median |ds_used| %.3f" % (lo, hi, s.sum(), np.median(ds[s]), np.median(dsu[s])))
print("   overall sd of s_best %.3f, of s_used %.3f" % (S.std(), SU.std()))
# Q4: trimmed
print("\nE. trimmed vs not (all clusters)")
for tag, sel in (("trimmed", [r for r in rs if r["trimmed"]]), ("not trimmed", [r for r in rs if not r["trimmed"]])):
    print("   %-12s n %3d  none %.1f used %.1f oracle %.1f free3 %.1f free9 %.1f | s_best-s_used median %+.3f sd %.3f | HF share %.1f dB | noise-null %.1f" % (tag, len(sel),
          *[np.median([m(r, k) for r in sel]) for k in ("none", "used", "oracle", "free3", "free9")], np.median([r["s_best"] - r["s_used"] for r in sel]), np.std([r["s_best"] - r["s_used"] for r in sel]),
          np.median([r["frac_hf_energy_db"] for r in sel]), np.median([np.mean(r["noise_null"]) for r in sel])))
tr = [r for r in rs if r["trimmed"]]
x = np.array([r["t_resp0"] for r in tr]); y_ = np.array([r["s_best"] - r["s_used"] for r in tr]); y2 = np.array([m(r, "none") - m(r, "oracle") for r in tr])
print("   trimmed: window start in response time %.2f..%.2f s ; corr(t_resp0, s_best - s_used) = %.2f ; corr(t_resp0, oracle gain) = %.2f ; corr(t_resp0, null none) = %.2f" % (x.min(), x.max(), np.corrcoef(x, y_)[0, 1], np.corrcoef(x, y2)[0, 1], np.corrcoef(x, [m(r, 'none') for r in tr])[0, 1]))
allx = np.array([r["t_resp0"] for r in rs]); ally = np.array([r["s_best"] - r["s_used"] for r in rs])
print("   all: corr(t_resp0, s_best - s_used) = %.2f ; corr(t_resp0, s_best) = %.2f" % (np.corrcoef(allx, ally)[0, 1], np.corrcoef(allx, [r["s_best"] for r in rs])[0, 1]))
# free9 shape check: average complex d on the failing cluster vs rank-1 direction
print("\nF. what the free fit adds beyond the scalar (median free3 - oracle, free9 - oracle):", np.round(np.median([m(r, "free3") - m(r, "oracle") for r in rs]), 1), np.round(np.median([m(r, "free9") - m(r, "oracle") for r in rs]), 1))

# G. estimator variants (second pass): null with the capture's own run mean, own front scalar, the line's closest point to z = 0
print("\nG. scalar rules (median 10-12 kHz null, dB; 'worse%' = worse than none by > 0.5 dB), by cluster x fail")
g = collections.defaultdict(list)
for r in rs: g[(cl(r), r["fail"])].append(r)
print("%-17s %-5s %3s | none   used   run-mean own    zline oracle | worse%%: used run own zline | med |s_run-sb| neigh_out_frac neigh_span_min n_run_front" % ("cluster", "fail", "n"))
for (c, f), v in sorted(g.items()):
    rr = [r for r in v if "run_null" in r]; oo = [r for r in v if "own_null" in r]
    w = lambda lst, k: 100 * np.mean([m(r, k) > m(r, "none") + 0.5 for r in lst]) if lst else float("nan")
    print("%-17s %-5s %3d | %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f |  %3.0f %3.0f %3.0f %3.0f | %6.3f   %.2f   %6.1f   %4d" % (c, f, len(v),
        np.median([m(r, "none") for r in v]), np.median([m(r, "used") for r in v]), np.median([m(r, "run_null") for r in rr]) if rr else np.nan,
        np.median([m(r, "own_null") for r in oo]) if oo else np.nan, np.median([m(r, "zline_null") for r in v]), np.median([m(r, "oracle") for r in v]),
        w(v, "used"), w(rr, "run_null"), w(oo, "own_null"), w(v, "zline_null"),
        np.median([abs(r["s_run"] - r["s_best"]) for r in rr]) if rr else np.nan,
        np.mean([r.get("neigh_out_frac", np.nan) for r in v]), np.nanmedian([r.get("neigh_span_s", np.nan) for r in v]) / 60, int(np.median([r["n_run_front"] for r in v]))))
allr = [r for r in rs if "run_null" in r]
print("\nH. all sampled captures with >= 3 front-tap captures in their run (n=%d): median none %.1f used %.1f run-mean %.1f oracle %.1f ; worse>0.5: used %.0f%% run-mean %.0f%%" % (
    len(allr), *[np.median([m(r, k) for r in allr]) for k in ("none", "used", "run_null", "oracle")],
    100 * np.mean([m(r, "used") > m(r, "none") + 0.5 for r in allr]), 100 * np.mean([m(r, "run_null") > m(r, "none") + 0.5 for r in allr])))
fo = [r for r in rs if "neigh_out_frac" in r]
print("   used-scalar neighbours from another run: %.0f%% of captures have >= 1 of 15 neighbours outside their run; failing: %.0f%%, non-failing: %.0f%%" % (
    100 * np.mean([r["neigh_out_frac"] > 0 for r in fo]), 100 * np.mean([r["neigh_out_frac"] > 0 for r in fo if r["fail"]]), 100 * np.mean([r["neigh_out_frac"] > 0 for r in fo if not r["fail"]])))
print("   |s_used - s_best| vs neighbours outside run: fraction out = 0: %.3f ; > 0: %.3f" % (
    np.median([abs(r["s_used"] - r["s_best"]) for r in fo if r["neigh_out_frac"] == 0]), np.median([abs(r["s_used"] - r["s_best"]) for r in fo if r["neigh_out_frac"] > 0])))
print("   distance of z = 0 from the shape line (split-norm units): %.3f ; s at the closest point %.3f (s_best median %.3f)" % (rs[0]["zline_dist"], rs[0]["s_zero"], np.median([r["s_best"] for r in rs])))
