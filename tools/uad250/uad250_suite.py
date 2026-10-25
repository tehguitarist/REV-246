"""UAD EMT 250 headless capture suite (docs/UAD250_test_plan.md).

Generates every test input, then renders each case through the UAD EMT 250 AU with tools/uad250/au_render (one
process at a time, under tools/memguard.py, checkpointed: finished cases are skipped).

    .venv/bin/python tools/uad250/uad250_suite.py build                  # compile au_render
    .venv/bin/python tools/uad250/uad250_suite.py list [--only PAT]      # cases, audio seconds, render estimate
    .venv/bin/python tools/uad250/uad250_suite.py run [--only PAT] [--smoke]

Outputs: analysis/uad250/captures/<group>/<case>.wav (+ .json: params, events, latency, plugin version, render time).
Inputs:  analysis/uad250/captures/inputs/<case>.wav. Both are gitignored (large); the manifest is analysis/uad250/manifest.json.
All inputs are mono float32 and are fed identically to L and R (au_render --mono-in) unless the case says "stereo".
"""
import argparse, fnmatch, hashlib, json, os, plistlib, subprocess, sys, time
import numpy as np
import soundfile as sf

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "uad250")
CAP = os.path.join(OUT, "captures")
BIN = os.path.join(OUT, "au_render")
PLUGIN = "/Library/Audio/Plug-Ins/Components/UAD EMT 250.component/Contents/Info.plist"

TIMES = ["0.4 s", "0.6 s", "0.8 s", "1.0 s", "1.2 s", "1.5 s", "1.8 s", "2.0 s", "2.2 s", "2.5 s", "2.8 s", "3.0 s",
         "3.3 s", "3.6 s", "4.0 s", "4.5 s"]
LFS = ["x 0.5", "x 1.0", "x 1.5", "x 2.0"]
HFS = ["x 0.25", "x 0.33", "x 0.50", "max"]
RT = {t: float(t.split()[0]) for t in TIMES}
# Base state: Reverb program, 100 % wet, no pre-delay, LF x1.0, HF x0.50, modulation ("Noise") as given per case.
BASE = {"Program": "Reverb", "Mix": 1, "WetSolo": 1,  # WetSolo=1 is 100 % wet regardless of Mix; cases that need dry set WetSolo=0
         "Predelay": "0 ms", "LFDecay": "x 1.0", "HFDecay": "x 0.50",
        "Outputs": "Front", "Power": 1, "Noise": 1, "Time": "2.5 s"}


# ---------------------------------------------------------------- signals (mono, float32)
def db(x): return 10 ** (x / 20)
def silence(sec, r): return np.zeros(int(sec * r), np.float32)
def impulse(amp_db, r, at=100, n=None, sec=0.5):
    x = np.zeros(n if n else int(sec * r), np.float32); x[at] = db(amp_db); return x
def train(r, n, spacing, amp_db=-6, jitter=0, lead=0.1):
    x = np.zeros(int((lead + n * spacing) * r), np.float32)
    for k in range(n): x[int(lead * r) + int(round(k * spacing * r)) + k * jitter] = db(amp_db)
    return x
def tone(f, sec, r, level=-20, fade=1.0):
    t = np.arange(int(sec * r)) / r; x = np.sin(2 * np.pi * f * t) * db(level)
    n = int(fade * r); x[:n] *= np.linspace(0, 1, n); return x.astype(np.float32)
def two_tone(f1, f2, sec, r, level=-26):
    return (tone(f1, sec, r, level) + tone(f2, sec, r, level)).astype(np.float32)
def noise_rms(sec, r, rms_db=-18, seed=2):
    x = np.random.default_rng(seed).standard_normal(int(sec * r)); return (x / x.std() * db(rms_db)).astype(np.float32)
def noise(sec, r, level=-20, seed=1):
    return (np.random.default_rng(seed).standard_normal(int(sec * r)) * db(level) / 3).astype(np.float32)
def ess(r, sec=10.0, f0=20.0, f1=20000.0, level=-12):
    f1 = min(f1, 0.48 * r); t = np.arange(int(sec * r)) / r; K = sec / np.log(f1 / f0)
    x = np.sin(2 * np.pi * f0 * K * (np.exp(t / K) - 1)) * db(level)
    n = int(0.01 * r); x[:n] *= np.linspace(0, 1, n); x[-n:] *= np.linspace(1, 0, n); return x.astype(np.float32)
def sweeps(r, n, gap, level=-12):
    s = ess(r, level=level); x = np.zeros(int(0.5 * r) + n * (len(s) + int(gap * r)), np.float32)
    for k in range(n): a = int(0.5 * r) + k * (len(s) + int(gap * r)); x[a:a + len(s)] = s
    return x
def stereo_one(x, ch): y = np.zeros((len(x), 2), np.float32); y[:, ch] = x; return y


# ---------------------------------------------------------------- cases
CASES = []
def case(group, name, x, r=48000, tail=8.0, params=None, events=None, note=""):
    p = dict(BASE); p.update(params or {})
    CASES.append(dict(group=group, name=name, x=x, r=r, tail=tail, params=p, events=events or [], note=note,
                      stereo=(x.ndim == 2)))
def rt_tail(t, extra=2.0): return min(3.0 * RT[t] + extra, 16.0)   # tail to the -98 dBFS floor, capped

def build_cases():
    r = 48000
    # A. static structure, modulation off (LTI): all 16 Time x 4 LF x 4 HF as an impulse, plus sweep cross-checks
    for t in TIMES:
        for lf in LFS:
            for hf in HFS:
                n = f"imp_T{t.split()[0]}_LF{lf[2:]}_HF{hf.split()[-1]}"
                case("A_static", n, impulse(0, r), tail=rt_tail(t), params=dict(Time=t, LFDecay=lf, HFDecay=hf, Noise=0))
    for t in ["0.4 s", "1.0 s", "2.5 s", "4.5 s"]:
        case("A_static", f"ess_T{t.split()[0]}", sweeps(r, 1, 0), tail=rt_tail(t, 3), params=dict(Time=t, Noise=0))
    # B. non-linear / level: impulse and sine ladders, noise floor, dead zone, overload
    for lv in [12, 6, 0, -6, -12, -18, -24, -30, -36, -42, -48, -54, -60, -66, -72, -78, -84, -90]:
        for t in ["0.4 s", "4.5 s"]:
            case("B_level", f"imp_{lv:+d}dB_T{t.split()[0]}", impulse(lv, r), tail=rt_tail(t, 4), params=dict(Time=t, Noise=0))
    for f in [250, 1000, 4000]:
        for lv in [-60, -40, -20, -6, 0, 6]:
            case("B_level", f"sine{f}_{lv:+d}dB", tone(f, 8, r, lv), tail=1.5, params=dict(Time="0.4 s", Noise=0))
    case("B_level", "twotone_1000_1100_-26dB", two_tone(1000, 1100, 8, r), tail=1.5, params=dict(Time="0.4 s", Noise=0))
    for nz in (0, 1):
        for pw in (1, 0):
            case("B_level", f"silence_noise{nz}_power{pw}", silence(20, r), tail=0, params=dict(Noise=nz, Power=pw, Time="4.5 s"))
    for pw, mx in [(0, 1), (1, 0), (1, 1)]:
        case("B_level", f"whitenoise_-20dB_power{pw}_mix{mx}", noise(20, r), tail=2, params=dict(Time="4.5 s", Power=pw, Mix=mx, WetSolo=0, Noise=0))
    # Ghost: below about -66 dBFS an impulse gives a level-independent response (~ -100 dBFS). Polarity, depth, density, channel.
    for lv in (-60, -66, -72, -90, -110, -130, -150):
        for sg, sn in ((1, "pos"), (-1, "neg")):
            case("B_ghost", f"imp_{lv}dB_{sn}", sg * impulse(lv, r), tail=3, params=dict(Time="0.4 s", Noise=0))
    case("B_ghost", "imp_zero", impulse(-300, r), tail=3, params=dict(Time="0.4 s", Noise=0))
    case("B_ghost", "imp_-90dB_Lonly", stereo_one(impulse(-90, r), 0), tail=3, params=dict(Time="0.4 s", Noise=0))
    case("B_ghost", "imp_-90dB_Ronly", stereo_one(impulse(-90, r), 1), tail=3, params=dict(Time="0.4 s", Noise=0))
    case("B_ghost", "imp_-90dB_x2", (impulse(-90, r) + impulse(-90, r, at=2000)).astype(np.float32), tail=3, params=dict(Time="0.4 s", Noise=0))
    for lv in (-60, -80, -100):
        case("B_ghost", f"noise_{lv}rms_10s", noise_rms(10, r, lv), tail=2, params=dict(Time="0.4 s", Noise=0))
    for t in ("0.4 s", "4.5 s"):
        case("B_ghost", f"zero_input_T{t.split()[0]}", impulse(-300, r, n=int(0.5 * r)), tail=rt_tail(t, 6), params=dict(Time=t, Noise=0), note="start-up ghost, subtract from other captures at the same setting")
    for lv in (-17, -16, -15, -14, -13):
        case("B_ghost", f"imp_{lv}dB_T4.5", impulse(lv, r), tail=14, params=dict(Time="4.5 s", Noise=0))
    case("B_ghost", "zero_input_T2.5", impulse(-300, r, n=int(0.5 * r)), tail=10, params=dict(Time="2.5 s", Noise=0))
    for lv in (-36, -30, -24, -18, -12, -6):
        case("B_ghost", f"noise_{lv}rms_30s_T2.5", noise_rms(30, r, lv), tail=10, params=dict(Time="2.5 s", Noise=0))
    case("B_ghost", "imp_-90dB_T4.5", impulse(-90, r), tail=14, params=dict(Time="4.5 s", Noise=0))
    # Calibration: 0 VU = -18 dBFS RMS (user, unit calibration). Sine peak -15 dBFS = -18 dBFS RMS; noise below is exactly -18 RMS.
    for lv, nm in [(-15, "0VU"), (-9, "+6VU"), (-3, "+12VU")]:
        for f in (250, 1000, 4000):
            case("B_level", f"cal_sine{f}_{nm}", tone(f, 8, r, lv), tail=1.5, params=dict(Time="0.4 s", Noise=0))
    for t in ("0.4 s", "2.5 s", "4.5 s"):
        case("B_level", f"cal_noise_-18rms_T{t.split()[0]}", noise_rms(30, r), tail=rt_tail(t, 2), params=dict(Time=t, Noise=0))
        case("B_level", f"cal_noise_-18rms_T{t.split()[0]}_noise1", noise_rms(30, r), tail=rt_tail(t, 2), params=dict(Time=t, Noise=1))
    # C. modulation / time variance (Noise on unless noted)
    for t in ["1.0 s", "2.5 s", "4.5 s"]:
        case("C_mod", f"train12s_T{t.split()[0]}_noise1", train(r, 24, 12.0), tail=12, params=dict(Time=t, Noise=1))
        case("C_mod", f"train12s_T{t.split()[0]}_noise0", train(r, 4, 12.0), tail=12, params=dict(Time=t, Noise=0))
    case("C_mod", "train12s_plus1sample_T2.5_noise1", train(r, 24, 12.0, jitter=1), tail=12, params=dict(Time="2.5 s"))
    case("C_mod", "train12s_plus1sample_T2.5_noise0", train(r, 4, 12.0, jitter=1), tail=12, params=dict(Time="2.5 s", Noise=0))
    for f in [250, 1000, 4000]:
        for nz in (1, 0):
            case("C_mod", f"tone{f}_60s_T2.5_noise{nz}", tone(f, 60, r), tail=4, params=dict(Time="2.5 s", Noise=nz))
    for lv, nm in [(-15, "0VU"), (-3, "+12VU")]:
        case("C_mod", f"tone1000_60s_T2.5_{nm}_noise1", tone(1000, 60, r, lv), tail=4, params=dict(Time="2.5 s"))
    case("C_mod", "tone1000_60s_T4.5_noise1", tone(1000, 60, r), tail=4, params=dict(Time="4.5 s"))
    case("C_mod", "tone1000_60s_T0.4_noise1", tone(1000, 60, r), tail=4, params=dict(Time="0.4 s"))
    for k, at in enumerate([0, 1, 2, 3, 100, 511, 512, 513, 1000, 24000, 48000, 96000, 177600]):
        case("C_mod", f"onset_{at:07d}smp_T2.5_noise1", impulse(0, r, at=at, n=at + 100), tail=8, params=dict(Time="2.5 s"))
    case("C_mod", "onset_repeat_control_T2.5_noise1", impulse(0, r, at=100, n=200), tail=8, params=dict(Time="2.5 s"))
    case("C_mod", "reset_reseed_T2.5_noise1", train(r, 2, 5.0, lead=0.1), tail=8, params=dict(Time="2.5 s"),
         events=[dict(t=5.0, action="reset")], note="impulses at 0.1 / 5.1 s, plugin reset at 5.0 s: equal responses = reset re-seeds")
    case("C_mod", "noise_toggle_off_on_tone1000", tone(1000, 20, r), tail=3, params=dict(Time="2.5 s", Noise=1),
         events=[dict(t=8.0, param="Noise", value=0), dict(t=14.0, param="Noise", value=1)])
    case("C_mod", "noise_toggle_midtail_imp", impulse(0, r, at=4800, n=r * 2), tail=10, params=dict(Time="4.5 s", Noise=1),
         events=[dict(t=2.0, param="Noise", value=0)])
    for n_ in (1, 0):
        case("C_mod", f"sweep_x6_gap8s_T2.5_noise{n_}", sweeps(r, 6, 8.0), tail=6, params=dict(Time="2.5 s", Noise=n_))
    # C_dense: the time-varying response h(t_onset, tau) sampled densely. The state follows absolute stream time, tails must not
    # overlap, so Time 0.4 s with a 1 s train at four 0.25 s offsets = 0.25 s sampling over 240 s; Time 2.5 s at 3 s spacing.
    for k in range(4):
        case("C_dense", f"dense_T0.4_off{k * 250:03d}ms", train(r, 240, 1.0, lead=0.1 + 0.25 * k), tail=2, params=dict(Time="0.4 s", Noise=0))
    for k in range(2):
        case("C_dense", f"dense_T2.5_sp3_off{k * 1500:04d}ms", train(r, 80, 3.0, lead=0.1 + 1.5 * k), tail=4, params=dict(Time="2.5 s", Noise=0))
    # D. other programs and time/rate coverage
    for i, nm in [(0, "I"), (1, "II"), (2, "III"), (3, "IV")]:
        for nz in (1, 0):
            case("D_programs", f"chorus{nm}_train2s_noise{nz}", train(r, 8, 2.0), tail=3, params=dict(Program="Chorus", ChorusMode=nm, Noise=nz))
            case("D_programs", f"chorus{nm}_tone1000_30s_noise{nz}", tone(1000, 30, r), tail=2, params=dict(Program="Chorus", ChorusMode=nm, Noise=nz))
    for nz in (1, 0):
        case("D_programs", f"space_train2s_noise{nz}", train(r, 8, 2.0), tail=4, params=dict(Program="Space", Noise=nz))
        case("D_programs", f"space_tone1000_30s_noise{nz}", tone(1000, 30, r), tail=3, params=dict(Program="Space", Noise=nz))
    for ld, rd in [(0, 0), (1, 0), (0, 1), (2, 2), (31, 31), (63, 63), (63, 0), (0, 63)]:
        case("D_programs", f"delay_L{ld}_R{rd}_96k", impulse(0, 96000, sec=1), r=96000, tail=1,
             params=dict(Program="Delay", LDelayTime=f"{ld * 5} ms", RDelayTime=f"{rd * 5} ms", Noise=0))
    for e in [1, 2, 31, 63]:
        case("D_programs", f"echo_E{e}_96k", impulse(0, 96000, sec=1), r=96000, tail=3,
             params=dict(Program="Echo", Echo=f"{e * 5} ms", LDelayTime="0 ms", RDelayTime="0 ms", Noise=0))
    for i in range(16):
        case("D_programs", f"phase_{i:02d}_96k", impulse(0, 96000, sec=0.5), r=96000, tail=0.5, params=dict(Program="Phase", Phase=None, Noise=0))
        CASES[-1]["params"]["Phase"] = ["0.0 ms", "0.1 ms", "0.2 ms", "0.4 ms", "0.8 ms", "1.6 ms", "3.2 ms", "6.0 ms", "9.0 ms",
                                        "11.8 ms", "13.4 ms", "14.2 ms", "14.6 ms", "14.8 ms", "14.9 ms", "15.0 ms"][i]
    for pd_ in ["0 ms", "20 ms", "40 ms", "60 ms"]:
        case("D_programs", f"predelay_{pd_.split()[0]}ms_96k", impulse(0, 96000, sec=0.5), r=96000, tail=4, params=dict(Predelay=pd_, Time="1.0 s", Noise=0))
    for rr in (44100, 48000, 88200, 96000, 192000):
        case("D_rate", f"imp_{rr}_T1.0", impulse(0, rr, sec=0.5), r=rr, tail=4, params=dict(Time="1.0 s", Noise=0))
        case("D_rate", f"imp_{rr}_T1.0_noise1", impulse(0, rr, sec=0.5), r=rr, tail=4, params=dict(Time="1.0 s", Noise=1))
        case("D_rate", f"whitenoise_{rr}_mix0", noise(10, rr), r=rr, tail=0.5, params=dict(Time="1.0 s", Mix=0, Noise=0))
    for mx in (0, 0.25, 0.5, 0.75, 1):
        case("D_route", f"imp_mix{mx}", impulse(0, r), tail=6, params=dict(Time="1.0 s", Mix=mx, WetSolo=0, Noise=0))
    case("D_route", "imp_mix1_wetsolo1_check", impulse(0, r), tail=6, params=dict(Time="1.0 s", Mix=0, WetSolo=1, Noise=0), note="WetSolo=1 forces 100 % wet whatever Mix says")
    case("D_route", "imp_outputs_rear", impulse(0, r), tail=6, params=dict(Time="1.0 s", Outputs="Rear", Noise=0))
    for ch, nm in [(0, "Lonly"), (1, "Ronly")]:
        case("D_route", f"imp_in_{nm}", stereo_one(impulse(0, r), ch), tail=6, params=dict(Time="1.0 s", Noise=0))
    case("D_route", "imp_in_both", impulse(0, r), tail=6, params=dict(Time="1.0 s", Noise=0))
    # F. parameter changes while ringing (automation)
    imp = impulse(0, r, at=4800, n=r)
    for a, b, at in [("4.5 s", "0.4 s", 3.0), ("0.4 s", "4.5 s", 0.3), ("4.5 s", "2.5 s", 3.0)]:
        case("F_auto", f"time_{a.split()[0]}_to_{b.split()[0]}_at{at}", imp, tail=16, params=dict(Time=a, Noise=0), events=[dict(t=at, param="Time", value=b)])
    for prm, a, b in [("LFDecay", "x 0.5", "x 2.0"), ("HFDecay", "x 0.25", "max")]:
        case("F_auto", f"{prm}_{a.split()[-1]}_to_{b.split()[-1]}_at3", imp, tail=14, params=dict(Time="4.5 s", Noise=0, **{prm: a}), events=[dict(t=3.0, param=prm, value=b)])
    for a, b in [("0 ms", "60 ms"), ("60 ms", "0 ms")]:
        case("F_auto", f"predelay_{a.split()[0]}_to_{b.split()[0]}_tone1000_mix0", tone(1000, 10, r), tail=2, params=dict(Mix=0, WetSolo=0, Predelay=a, Noise=0), events=[dict(t=5.0, param="Predelay", value=b)])
        case("F_auto", f"predelay_{a.split()[0]}_to_{b.split()[0]}_tone1000_wet", tone(1000, 10, r), tail=3, params=dict(Time="1.0 s", Predelay=a, Noise=0), events=[dict(t=5.0, param="Predelay", value=b)])
    case("F_auto", "delay_L0_to_63_tone1000", tone(1000, 10, r), tail=2, params=dict(Program="Delay", LDelayTime="0 ms", RDelayTime="0 ms", Noise=0), events=[dict(t=5.0, param="LDelayTime", value="315 ms")])
    case("F_auto", "mix_ramp_noise_5s", noise(8, r), tail=1, params=dict(Time="1.0 s", Mix=0, WetSolo=0, Noise=0),
         events=[dict(t=1.0 + i * 0.04, param="Mix", value=round(i / 100, 4)) for i in range(101)])
    case("F_auto", "power_toggle_midtail", imp, tail=10, params=dict(Time="4.5 s", Noise=0),
         events=[dict(t=1.0, param="Power", value=0), dict(t=2.0, param="Power", value=1)])
    case("F_auto", "program_switch_midtail", imp, tail=10, params=dict(Time="4.5 s", Noise=0),
         events=[dict(t=1.0, param="Program", value="Delay"), dict(t=2.0, param="Program", value="Reverb")])


def audio_seconds(c): return len(c["x"]) / c["r"] + c["tail"]
def path_in(c): return os.path.join(CAP, "inputs", f"{c['group']}__{c['name']}.wav")
def path_out(c): return os.path.join(CAP, c["group"], f"{c['name']}.wav")

def select(only, smoke):
    cs = [c for c in CASES if not only or any(fnmatch.fnmatch(f"{c['group']}/{c['name']}", p) for p in only.split(","))]
    if smoke:  # one case per group, shortest first: proves each pipeline end to end
        seen = {}
        for c in sorted(cs, key=audio_seconds): seen.setdefault(c["group"], c)
        cs = list(seen.values())
    return cs

def plugin_version():
    with open(PLUGIN, "rb") as f: return plistlib.load(f).get("CFBundleShortVersionString")

def run_case(c, limit_gb):
    os.makedirs(os.path.dirname(path_out(c)), exist_ok=True); os.makedirs(os.path.dirname(path_in(c)), exist_ok=True)
    sf.write(path_in(c), c["x"], c["r"], subtype="FLOAT")
    ev = path_out(c)[:-4] + ".events.json"
    json.dump(c["events"], open(ev, "w"))
    cmd = [sys.executable, os.path.join(ROOT, "tools", "memguard.py"), "--limit-gb", str(limit_gb), "--", BIN,
           path_in(c), path_out(c), "--rate", str(c["r"]), "--tail", str(c["tail"])]
    if not c["stereo"]: cmd += ["--mono-in"]
    for k, v in c["params"].items(): cmd += ["--set", f"{k}={v}"]
    if c["events"]: cmd += ["--auto", ev]
    t0 = time.time(); p = subprocess.run(cmd, capture_output=True, text=True); dt = time.time() - t0
    if p.returncode != 0 or not os.path.exists(path_out(c)):
        return False, (p.stdout + p.stderr)[-400:]
    lat = next((l for l in p.stdout.splitlines() if l.startswith("reported latency")), "")
    json.dump(dict(group=c["group"], name=c["name"], rate=c["r"], tail=c["tail"], params=c["params"], events=c["events"],
                   input_sha1=hashlib.sha1(c["x"].tobytes()).hexdigest(), plugin_version=plugin_version(),
                   latency=lat, render_s=round(dt, 1), note=c["note"], stereo_input=c["stereo"]),
              open(path_out(c)[:-4] + ".json", "w"), indent=1)
    return True, f"{dt:.1f}s"

def main():
    ap = argparse.ArgumentParser(usage=__doc__)
    ap.add_argument("cmd", choices=["build", "list", "run"]); ap.add_argument("--only", default="")
    ap.add_argument("--smoke", action="store_true"); ap.add_argument("--limit-gb", type=float, default=6.0)
    a = ap.parse_args()
    if a.cmd == "build":
        os.makedirs(OUT, exist_ok=True)
        subprocess.check_call(["swiftc", "-O", os.path.join(ROOT, "tools/uad250/au_render.swift"), "-o", BIN]); return
    build_cases(); cs = select(a.only, a.smoke)
    if a.cmd == "list":
        by = {}
        for c in cs: by.setdefault(c["group"], []).append(c)
        for g, l in by.items(): print(f"{g:12s} {len(l):4d} cases  {sum(audio_seconds(c) for c in l) / 60:7.1f} min audio")
        tot = sum(audio_seconds(c) for c in cs) / 60
        print(f"{'total':12s} {len(cs):4d} cases  {tot:7.1f} min audio (renders at roughly 4-5x real time)"); return
    done = fail = 0
    for i, c in enumerate(cs):
        if os.path.exists(path_out(c)[:-4] + ".json"): done += 1; continue   # checkpoint
        ok, msg = run_case(c, a.limit_gb)
        print(f"[{i + 1}/{len(cs)}] {c['group']}/{c['name']}: {'ok ' + msg if ok else 'FAILED ' + msg}", flush=True)
        fail += (not ok)
    print(f"done; {done} already present, {fail} failed")

if __name__ == "__main__":
    main()
