"""A/B audio examples for plugin/README.md: dry, impulse response (IR) and REV-246 versions of four recordings.

    .venv/bin/python tools/ab_samples.py [--stage render|encode] [--only name ...]

Per example: the recording is summed to mono (the unit is mono in), rendered through REV246Render (the plugin's own
engine, latency removed as a DAW would) and convolved with the IR (both channels, resampled to 48 kHz). The IR wet is
scaled to the plugin wet's RMS, so any difference heard is the reverb's character, not its level; both wets are then
set to the dry's RMS and mixed at MIX (0.8 dry + 0.2 wet). The three files of an example share one gain, so no clip
and no level change between them. Outputs: docs/listen/audio/<name>_{dry,ir,plugin}.mp3 plus `levels.json`.

IR sources (kept unnamed in the README): 1 = CLA capture, 2 = Altiverb, 3 = Tomkin library.
"""
import argparse, json, subprocess, sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import fftconvolve, resample_poly

ROOT = Path(__file__).resolve().parent.parent
RENDER = ROOT / "plugin/build-render/REV246Render_artefacts/Release/REV246Render"
OUT = ROOT / "docs/listen/audio"
WORK = ROOT / "analysis/ab_samples"
RATE = 48000
MIX = 0.2
TOMKIN = ROOT / "irs/EMT 246 captures/96kHz"
ALT = ROOT / "irs/Altiverb EMT256 IRs"

# name, source recording, IR source number, IR file, plugin settings (RenderIr keys), tail seconds, description
EXAMPLES = [
    dict(name="room-mic", wav="Room Mic", src=1, ir=ROOT / "irs/eugtone_EMT246_CLA_shortRT.wav",
         prog="S", args=dict(lf=1.0, rt=1.0, hf=0.9, ampl=80, delay=10, lp=12), tail=3.0,
         text="Program S, RT 1.0 s, LF 1.0, HF 0.9, AMPL 80 %, DELAY 10 ms, LP 12 kHz"),
    dict(name="piano", wav="Piano", src=2, preset="PRESET1 PIANO", ir=ROOT / "analysis/altiverb/cache/PRESET1 PIANO.npy",
         prog="M", args=dict(lf=1.2, rt=1.5, hf=0.7, ampl=60, delay=40, lp=12), tail=5.0,
         text="Program M, RT 1.5 s, LF 1.2, HF 0.7, AMPL 60 %, DELAY 40 ms, LP 12 kHz"),
    dict(name="acoustic-guitar", wav="Acoustic Guitar", src=3, ir=TOMKIN / "Short RT/12_0_60_10_14_06_dc.wav",
         prog="S", args=dict(lf=1.0, rt=1.4, hf=0.6, ampl=60, delay=0, lp=12), tail=3.0,
         text="Program S, RT 1.4 s, LF 1.0, HF 0.6, AMPL 60 %, DELAY 0, LP 12 kHz"),
    dict(name="snare-ir2", wav="Snare", src=2, preset="PRESET54 GATED", ir=ROOT / "analysis/altiverb/cache/PRESET54 GATED.npy",
         prog="N", args=dict(lf=-3.0, hf=0.8, ampl=0.2, delay=39, lp=12), tail=2.0,
         text="Program N, SPAN 0.2 s, EARLY RT -3.0 s, LATE RT 0.8 s, DELAY 39 ms, LP 12 kHz"),
    dict(name="snare-ir3", wav="Snare", src=3, ir=TOMKIN / "Non-Lin/12_0_02_-30_08_dc.wav",
         prog="N", args=dict(lf=-3.0, hf=0.8, ampl=0.2, delay=0, lp=12), tail=2.0,
         text="Program N, SPAN 0.2 s, EARLY RT -3.0 s, LATE RT 0.8 s, DELAY 0, LP 12 kHz"),
]


def rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


def read_ir(path):
    if str(path).endswith(".npy"):  # Altiverb: the deconvolved sweep capture (tools/altiverb_compare.py), 2 x N at 48 kHz
        return np.load(path).T
    h, r = sf.read(path, always_2d=True)
    if r != RATE:
        g = np.gcd(int(r), RATE)
        h = resample_poly(h, RATE // g, int(r) // g, axis=0)
    return h[:, :2] if h.shape[1] >= 2 else np.repeat(h, 2, axis=1)


def render_plugin(ex, mono_path, n_out, tag):
    out = WORK / f"{tag}_plugin_raw.wav"
    cmd = [str(RENDER), f"out={out}", f"rate={RATE}", f"program={ex['prog']}", f"len={n_out / RATE + 0.5:.4f}",
           f"input={mono_path}", "mod=1"] + [f"{k}={v}" for k, v in ex["args"].items()]
    line = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip().splitlines()[-1]
    lat = int(line.split("latency=")[1].split()[0])
    y, r = sf.read(out, always_2d=True)
    assert r == RATE, r
    return y[lat:lat + n_out, :2], lat


def onset(h, rel=0.02):
    """First sample at which any channel reaches rel x the file's peak."""
    a = np.abs(h).max(axis=1)
    return int(np.argmax(a >= rel * a.max()))


def plugin_onset(ex, tag):
    """The plugin's own first tap (latency removed): the IR is shifted so both start their reverb together."""
    out = WORK / f"{tag}_imp_raw.wav"
    cmd = [str(RENDER), f"out={out}", f"rate={RATE}", f"program={ex['prog']}", "len=1.5", "mod=1"] \
        + [f"{k}={v}" for k, v in ex["args"].items()]
    line = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip().splitlines()[-1]
    lat = int(line.split("latency=")[1].split()[0])
    y, _ = sf.read(out, always_2d=True)
    return onset(y[lat:]), lat


def build(ex, dry, tag):
    mono = dry.mean(axis=1)
    n_out = len(mono) + int(ex["tail"] * RATE)
    mono_path = WORK / f"{tag}_mono.wav"
    sf.write(mono_path, mono.astype(np.float32), RATE, subtype="FLOAT")

    wet_p, lat = render_plugin(ex, mono_path, n_out, tag)
    h = read_ir(ex["ir"])
    op, _ = plugin_onset(ex, tag)
    if ex["src"] == 2:  # Altiverb: use the cluster onset found by tools/altiverb_presets.py
        on = json.loads((ROOT / "analysis/altiverb/presets.json").read_text())
        on = (on.get("presets") or on)[ex["preset"]]["onset_ms"]
        shift = int(round(on * 1e-3 * RATE)) - op
    else:
        shift = onset(h) - op  # >0: the IR file starts its reverb later than the plugin (leading silence / pre-roll), trim it
    h = h[shift:] if shift >= 0 else np.pad(h, ((-shift, 0), (0, 0)))
    wet_i = np.stack([fftconvolve(mono, h[:, c]) for c in range(2)], axis=1)
    wet_i = np.pad(wet_i, ((0, max(0, n_out - len(wet_i))), (0, 0)))[:n_out]

    # refine the alignment on the wets themselves (the onset criterion is only good to a few samples; N's first taps
    # are small): the lag of the cross-correlation peak within +-20 ms, applied to the IR wet
    m = min(len(wet_p), 20 * RATE)
    c = fftconvolve(wet_i[:m, 0], wet_p[:m, 0][::-1])[m - 1 - 960:m + 960]
    lag = int(np.argmax(np.abs(c))) - 960
    if lag > 0:
        wet_i = np.pad(wet_i[lag:], ((0, lag), (0, 0)))
    elif lag < 0:
        wet_i = np.pad(wet_i, ((-lag, 0), (0, 0)))[:n_out]
    shift_total = shift + lag
    kernel = h[lag:] if lag >= 0 else np.pad(h, ((-lag, 0), (0, 0)))  # the IR as aligned, for ab_analysis's clock fit

    # equal RMS over the whole render, then the dry's RMS; the tails end with the same short fade
    fade = np.ones(n_out)
    nf = int(0.5 * RATE)
    fade[-nf:] = np.linspace(1, 0, nf)
    wet_p, wet_i = wet_p * fade[:, None], wet_i * fade[:, None]
    dry_p = np.pad(dry, ((0, n_out - len(dry)), (0, 0)))
    target = rms(dry) / np.sqrt(2) * np.sqrt(2)  # stereo RMS of the dry
    gp, gi = target / rms(wet_p), target / rms(wet_i)
    mix = lambda w: (1 - MIX) * dry_p + MIX * w
    outs = dict(dry=dry_p, ir=mix(wet_i * gi), plugin=mix(wet_p * gp))
    peak = max(np.abs(v).max() for v in outs.values())
    k = min(1.0, 0.89 / peak)  # one gain for the three, no clipping
    info = dict(latency=lat, ir_shift_ms=1000 * shift_total / RATE, plugin_wet_gain_db=20 * np.log10(gp), ir_wet_gain_db=20 * np.log10(gi),
                natural_plugin_vs_ir_db=20 * np.log10(gi / gp), common_gain_db=20 * np.log10(k))
    info['kernel'] = kernel
    return {a: b * k for a, b in outs.items()}, wet_p * gp, wet_i * gi, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    levels = {}
    for ex in EXAMPLES:
        if a.only and ex["name"] not in a.only:
            continue
        dry, r = sf.read(ROOT / "ref" / f"{ex['wav']}.wav", always_2d=True)
        assert r == RATE
        if dry.shape[1] == 1:
            dry = np.repeat(dry, 2, axis=1)
        outs, wp, wi, info = build(ex, dry, ex["name"])
        np.savez_compressed(WORK / f"{ex['name']}_wets.npz", plugin=wp.astype(np.float32), ir=wi.astype(np.float32),
                            kernel=info.pop('kernel').astype(np.float32), gain=np.float32(10 ** (info['ir_wet_gain_db'] / 20)))
        for kind, y in outs.items():
            if kind == "dry" and ex["name"] == "snare-ir3":
                continue  # shared with snare-ir2
            base = "snare" if ex["name"].startswith("snare") and kind == "dry" else ex["name"]
            wav = WORK / f"{base}_{kind}.wav"
            sf.write(wav, y.astype(np.float32), RATE, subtype="FLOAT")
            subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(wav), "-codec:a", "libmp3lame", "-b:a", "320k",
                            str(OUT / f"{base}_{kind}.mp3")], check=True)
        levels[ex["name"]] = dict(info, settings=ex["text"], ir_source=ex["src"])
        print(ex["name"], {k: round(v, 2) if isinstance(v, float) else v for k, v in info.items()})
    prev = json.loads((OUT / "levels.json").read_text()) if (OUT / "levels.json").exists() else {}
    prev.update(levels)
    (OUT / "levels.json").write_text(json.dumps(prev, indent=1))


if __name__ == "__main__":
    sys.exit(main())
