"""Offline bundle of the UAD EMT 250 captures, so fits run with no plugin and no UAD hardware.

    .venv/bin/python tools/uad250/offline_bundle.py manifest   # write analysis/uad250/manifest.json (every capture: params, sha1, size)
    .venv/bin/python tools/uad250/offline_bundle.py pack       # tar.zst of captures/ (inputs excluded) + manifest, in analysis/uad250/bundle/
    .venv/bin/python tools/uad250/offline_bundle.py check [--full]   # every manifest entry present (sizes; --full also sha1)
    .venv/bin/python tools/uad250/offline_bundle.py inputs     # regenerate the (seeded) input WAVs and verify their sha1 against the manifest

Analysis scripts in tools/uad250/analyze_*.py, comb_coherence.py and compare_altiverb.py read only captures/ and the Altiverb caches; none
imports or needs the plugin. Only `uad250_suite.py run` and au_render need the plugin. The captures are renders of a licensed plugin: they are
gitignored and must not be committed or redistributed; the manifest (hashes and settings only) is committed.
"""
import argparse, hashlib, json, os, subprocess, sys, time
import numpy as np, soundfile as sf
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import uad250_suite as su

MAN = os.path.join(su.OUT, "manifest.json")

def sha1(path, n=1 << 22):
    h = hashlib.sha1()
    with open(path, "rb") as f:
        while b := f.read(n): h.update(b)
    return h.hexdigest()

def cmd_manifest():
    ent = []
    for g in sorted(os.listdir(su.CAP)):
        d = os.path.join(su.CAP, g)
        if g == "inputs" or not os.path.isdir(d): continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".json") or fn.endswith(".events.json"): continue
            meta = json.load(open(os.path.join(d, fn))); w = os.path.join(d, fn[:-5] + ".wav")
            if not os.path.exists(w): continue
            ent.append(dict(group=g, name=meta["name"], rate=meta["rate"], params=meta["params"], events=meta["events"], tail=meta["tail"],
                            input_sha1=meta["input_sha1"], plugin_version=meta["plugin_version"], latency=meta["latency"], stereo_input=meta["stereo_input"],
                            note=meta.get("note", ""), bytes=os.path.getsize(w), sha1=sha1(w)))
    json.dump(dict(created=time.strftime("%Y-%m-%d"), plugin="UAD EMT 250 AU " + su.plugin_version(), captures=ent), open(MAN, "w"), indent=1)
    print(f"{len(ent)} captures, {sum(e['bytes'] for e in ent) / 1e9:.2f} GB -> {MAN}")

def cmd_pack():
    out = os.path.join(su.OUT, "bundle"); os.makedirs(out, exist_ok=True)
    tarf = os.path.join(out, f"uad250_captures_{time.strftime('%Y-%m-%d')}.tar.zst")
    cmd = f"tar -C '{su.OUT}' --exclude='captures/inputs' -cf - captures manifest.json | zstd -T0 -3 -o '{tarf}' -f"
    subprocess.check_call(cmd, shell=True); print(f"{tarf}: {os.path.getsize(tarf) / 1e9:.2f} GB")
    open(os.path.join(out, "README.txt"), "w").write(
        "UAD EMT 250 captures, rendered headlessly (tools/uad250/au_render) from the UAD EMT 250 AU. Licensed-plugin renders: do not redistribute.\n"
        "Restore: cd analysis/uad250 && zstd -dc bundle/<file>.tar.zst | tar -xf -\nVerify: .venv/bin/python tools/uad250/offline_bundle.py check --full\n"
        "Regenerate inputs: .venv/bin/python tools/uad250/offline_bundle.py inputs\n")

def cmd_check(full):
    m = json.load(open(MAN)); bad = 0
    for e in m["captures"]:
        w = os.path.join(su.CAP, e["group"], e["name"] + ".wav")
        if not os.path.exists(w) or os.path.getsize(w) != e["bytes"] or (full and sha1(w) != e["sha1"]): bad += 1; print("BAD", e["group"], e["name"])
    print(f"{len(m['captures'])} manifest entries, {bad} bad")

def cmd_inputs():
    su.build_cases(); m = {(e["group"], e["name"]): e for e in json.load(open(MAN))["captures"]}; bad = 0; n = 0
    for c in su.CASES:
        e = m.get((c["group"], c["name"]))
        if not e: continue
        n += 1; ok = hashlib.sha1(c["x"].tobytes()).hexdigest() == e["input_sha1"]
        if not ok: bad += 1; print("input differs:", c["group"], c["name"])
        else:
            os.makedirs(os.path.dirname(su.path_in(c)), exist_ok=True)
            if not os.path.exists(su.path_in(c)): sf.write(su.path_in(c), c["x"], c["r"], subtype="FLOAT")
    print(f"{n} inputs regenerated/verified, {bad} differ")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["manifest", "pack", "check", "inputs"]); ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    {"manifest": cmd_manifest, "pack": cmd_pack, "check": lambda: cmd_check(a.full), "inputs": cmd_inputs}[a.cmd]()
