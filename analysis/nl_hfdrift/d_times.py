"""Capture timestamps (bext OriginationDate/Time, Deconvolver) of the 96 kHz Non-Lin files -> d_times.json {key: epoch_s}"""
import os, struct, json, datetime, collections
D = "/Users/leighpierce/Git/246/irs/EMT 246 captures/96kHz/Non-Lin"
out = {}; nob = 0
for f in sorted(os.listdir(D)):
    if not f.endswith(".wav"): continue
    b = open(os.path.join(D, f), "rb").read(); i = 12; bext = None
    while i < len(b) - 8:
        cid = b[i:i + 4]; sz = struct.unpack("<I", b[i + 4:i + 8])[0]
        if cid == b"bext": bext = b[i + 8:i + 8 + sz]
        i += 8 + sz + (sz & 1)
    if bext is None: nob += 1; continue
    ds = bext[320:330].decode("latin1"); ts = bext[330:338].decode("latin1")
    t = datetime.datetime.strptime(ds + " " + ts, "%Y-%m-%d %H:%M:%S").timestamp()
    lp, delay, span, early, late = (int(v) for v in f.replace("_dc.wav", "").split("_"))
    out["_".join(str(v) for v in (span / 10, early / 10, late / 10))] = t
json.dump(out, open("/Users/leighpierce/Git/246/analysis/nl_hfdrift/d_times.json", "w"))
ts = sorted(out.values()); print(len(out), "files; no bext", nob, "; first", datetime.datetime.fromtimestamp(ts[0]), "last", datetime.datetime.fromtimestamp(ts[-1]))
by = collections.defaultdict(list)
for k, t in out.items(): by[k.split("_")[0]].append(t)
for s in sorted(by): print("SPAN", s, len(by[s]), datetime.datetime.fromtimestamp(min(by[s])).strftime("%m-%d %H:%M"), "..", datetime.datetime.fromtimestamp(max(by[s])).strftime("%m-%d %H:%M"))
