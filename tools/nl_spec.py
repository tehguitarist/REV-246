"""Program N (Non-Lin) for the plugin: a sample-by-sample reference implementation, the specification as data, and
reference IRs to check a port against (docs/NonLin_program.md §12).

    .venv/bin/python tools/nl_spec.py

simulate(x, span, early, late) runs the static model the way a plugin would, at the device rate (26,373.6 Hz), mono in,
(n, 2) out, with integer delay lines, all-passes and loops (branch_delay = 1 sample, the plugin's choice). main() checks it
against the frequency-domain bench model (tools/nl_model.render, band=False) and writes
analysis/nl/nl_plugin_spec.json and analysis/nl/ref/nl_<SPAN>_<EARLY>_<LATE>.wav (float32, device rate, (L, R), the
input impulse at sample 0, main-tap units: a main tap = 1).
"""
import json, os, sys
import numpy as np
import soundfile as sf
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ["NL_KERNEL"] = "0"          # references are the kernel-free device network
os.environ["NL_BK"] = "0"              # no branch read kernel in the references
os.environ["NL_WARP2"] = "1"           # the adopted long-kernel warp's kernel metadata
os.environ["NL_BK_ALPHA"] = "0.7"
os.environ["NL_BK_LP"] = "11000"
os.environ.pop("NL_BK_FILE", None)
import nl_early as ne  # noqa: E402
import nl_late as nlt  # noqa: E402
import nl_model as nm  # noqa: E402
import nl_structure as nl  # noqa: E402

FS = nl.FSD
REF_SETTINGS = [(0.0, 10.0, 1.0), (0.2, 1.0, 1.0), (0.5, -1.0, 2.0), (0.3, 3.0, 0.0), (0.4, -3.0, 0.6), (0.1, 0.3, 1.4)]
SPEC = os.path.join(nl.OUT, "nl_plugin_spec.json")
REF = os.path.join(nl.OUT, "ref")


def delay(x, d):
    if d <= 0:
        return x.copy()
    return np.r_[np.zeros(min(d, len(x))), x[:max(len(x) - d, 0)]]


def ap(x, M, g=ne.G_AP):
    """Schroeder all-pass (-g + z^-M) / (1 - g z^-M), g = -1/sqrt(2)."""
    return lfilter(np.r_[-g, np.zeros(M - 1), 1.0], np.r_[1.0, np.zeros(M - 1), -g], x)


BRANCH = 1   # the AP(124) branch one sample after the table positions (the captures' modal offset; R's first-tap direct term then coincides with its tap)


def simulate(x, span, early, late, branch=BRANCH):
    """Program N's static network on a mono input x (device rate); returns (len(x), 2), main-tap units. branch: the
    AP(124) branch's delay re the front taps beyond the tables (1 for the plugin; 0 reproduces the bench's frame).
    Kernel-free: the adopted branch_read_kernel is not applied here, so the reference IRs are kernel-free (see the
    spec field's reference_irs)."""
    ns = ne.N_STAGES[span]
    out = np.zeros((len(x), 2))
    u = delay(ap(x, ne.AP_IN), branch)                            # the branch: AP(124) at the input
    for ch in (0, 1):
        fdb, sdb = ne.law_db(span, early, ch)
        # front taps (no diffuser)
        for i, f in enumerate(ne.FRONT[ch]):
            out[:, ch] += ne.FRONT_AMP[ch][i] * 10 ** (fdb[i] / 20) * delay(x, f)
        # main taps on the AP(124) line, summed, through AP(180) L / AP(188) R
        acc = np.zeros(len(x))
        for t, s, i in ne.taps(ch, ns):
            g = s * 10 ** (sdb[[j for j, st in enumerate(ne.STAGES) if i in st][0]] / 20)
            acc += g * delay(u, t)
            if ch == 1 and i == 0:
                # R's first tap: the diffuser's direct term (0.5 x) one sample later, in the front taps' frame
                d0 = 0.5
                out[:, ch] += g * d0 * (delay(x, ne.R_FIRST_DIRECT) - delay(x, t + branch))   # 0 at branch = 1
        out[:, ch] += ap(acc, ne.AP_OUT[ch])
        # late: the AP(124) line at the feed, AP(314), three flat loops, four read taps each
        if late > 0:
            v = 10 ** (nm.late_db(span, early) / 20) * ap(delay(u, nlt.FEED[span] + (nlt.R_LATE if ch else 0)),
                                                           nlt.AP_LATE)
            for P, qs in nlt.READS.items():
                line = lfilter([1.0], np.r_[1.0, np.zeros(P - 1), -nlt.loop_gain(P, late)], v)
                for s, q in zip(nlt.READ_SIGNS[ch][P], qs):
                    out[:, ch] += nm.LATE_AMP * s * delay(line, q)
    return out


def spec():
    lab = lambda v, c, src: dict(value=v, label=c, source=src)
    return {
        "program": "N (Non-Lin), docs/NonLin_program.md §9-12",
        "device_rate_hz": FS,
        "units": "a main tap = 1; delays in integer device samples from the input; mono in, (L, R) out",
        "front_taps": lab(dict(L=dict(delay=list(ne.FRONT[0]), amp=list(ne.FRONT_AMP[0]), law_tick=list(ne.FRONT_TICKS[0])),
                               R=dict(delay=list(ne.FRONT[1]), amp=list(ne.FRONT_AMP[1]), law_tick=list(ne.FRONT_TICKS[1])),
                               note="plain taps, no diffuser, gains on the EARLY law at law_tick x 4 ms"), "S", "§9.2-9.3"),
        "branch_input": lab(dict(allpass=dict(D=ne.AP_IN, g=ne.G_AP), form="(-g + z^-D) / (1 - g z^-D)"), "S", "§9.2"),
        "main_taps": lab(dict(L_delay=list(ne.TAPS), R_delay=[t + ne.R_OFF for t in ne.TAPS], R_offset=ne.R_OFF,
                              L_sign=ne.SIGNS[0], R_sign=ne.SIGNS[1],
                              stages=ne.STAGES, out_allpass=dict(L=ne.AP_OUT[0], R=ne.AP_OUT[1], g=ne.G_AP),
                              branch_delay=dict(value=BRANCH, note="the whole AP(124) branch (main taps and late) sits this many samples after "
                                                "the table positions re the front taps: the captures' modal offset; with it R's first tap (464 + 1) "
                                                "is where its direct term was measured (465), so no special case is needed"),
                              R_first_direct=ne.R_FIRST_DIRECT,
                              R_first_tap="in captures whose branch is jittered its direct term stays at 465 (§9.4); none at branch = 1"),
                         "S", "§9.2"),
        "span": lab(dict(stages_running={str(k): v for k, v in ne.N_STAGES.items()},
                         between_steps="H: a stage runs if its last tap is at least 15-33 ms before SPAN"),
                    "C (library steps) / H (between them)", "§9.2"),
        "early_law": lab(dict(tick_s=ne.TICK, stage_s=round(9 * ne.TICK, 6),
                              positive="stage j at tau = 16 + 36 j ms, front taps at their tick; gain -60 tau / EARLY dB",
                              negative="last stage at tau_end = SPAN + 12 ms, others 36 ms apart before it; gain +60 (tau - tau_end) / |EARLY| dB",
                              flat="EARLY 10: 0 dB everywhere"), "S", "§9.3"),
        "late": lab(dict(feed=dict(delay_on_branch_line={str(k): v for k, v in nlt.FEED.items()}, R="220 samples earlier"),
                         gain="EARLY > 0: the early law at stage min(ns, 13) (tau = 16 + 36 min(ns, 13) ms); EARLY < 0 or 10: 0 dB",
                         allpass=dict(D=nlt.AP_LATE, g=ne.G_AP),
                         loops=dict(P=nlt.LOOPS, form="line[n] = v[n] + G line[n - P], G = 10^(-3 P / (fs LATE)), no filter"),
                         reads={str(P): dict(delay=list(q), L_sign=list(nlt.READ_SIGNS[0][P]), R_sign=list(nlt.READ_SIGNS[1][P]))
                                for P, q in nlt.READS.items()},
                         read_amp=nm.LATE_AMP, late_0="LATE 0: no late section at all (gated)"), "S", "§10"),
        "output_eq": lab(dict(table=None, form="none: the plugin's N output is flat to the 12.8 kHz analysis edge",
                               note="the -0.24 / -1.03 / -1.85 / -3.0 dB roll-off at 11.5 / 12 / 12.3 / 12.6 kHz that 11.3 first adopted as an EQ is the harness's own "
                                   "96 kHz -> device-rate resampler (resample_poly 25/91, Kaiser 5; -0.22 / -0.90 / -1.50 / -2.33 / -3.65 dB at "
                                   "11.5 / 12 / 12.25 / 12.5 / 12.8 kHz, nl_structure.resampler_gain), not a property of the unit: refitted after "
                                   "the clock warp was fixed it matches that filter to 0.05 dB up to 12 kHz. The LOW PASS therefore needs no "
                                   "special 12 kHz shape on N. (The plugin's low-band CaptureEq is the separate model_s_results 27 transfer.) The band-edge term above 10 kHz that "
                                   "the bench applies per capture run (NonLin_program 11.6, hf_corr_g5.json) belongs to the capture chain, not to this spec."),
                         "S", "§11.3, §11.6"),
        "branch_read_kernel": lab(dict(taps={("L", "R")[int(c)]: v for c, v in json.load(open(nm.BK_FILE)).items()},
                                       centre=nm.BK_PAD, band_limit_hz=nm.BK_LP, strength=float(nm.BK_ALPHA),
                                       use="kernel = delta + strength * (taps - delta), per channel (the stored taps are the fitted, strength 1 kernel)",
                                       applies_to="the whole AP(124) branch (main taps, R's first-tap direct term, late section); not the front taps",
                                        reference_irs="kernel-free: analysis/nl/ref/*.wav come from simulate(), which omits this kernel. Null a "
                                                      "kernel-on port with the kernel bypassed, or expect about -59 dB worst on these six settings "
                                                      "(-62 dB median) against these files; the port's own reference should be rendered with the kernel "
                                                      "and BK_PAD samples of look-ahead (nl_model.render applies it zero-latency via a -BK_PAD roll)",
                                       latency_samples=nm.BK_PAD,
                                       status="adopted (27 Oct 2026; refitted and set to strength 0.7 under the long-kernel clock warp, NonLin_program 11.5-11.6): "
                                              "0.05-0.15 / 0.15-0.3 / 0.3-0.6 s medians -58.7 / -57.6 / -52.2 dB (1,044 captures held out; the 0.65 kernel on the old warp "
                                               "gave -58.4 / -57.3 / -52.0), no channel-window worse than the previous adopted bench by more than +0.07 dB (0.0-0.6 s); against "
                                               "the same bench without the kernel the worst loss is +0.81 dB (22 channel-windows over +0.5 at 0.05-0.15 s, plus 2+2 at "
                                               "0.15-0.3 / 0.3-0.6 s); the strength-0.5 screen used a 151-capture subset and was not benched in full; the medians mix L and R",
                                        note="deviations from a unit tap are <= 0.00073 (-62.8 dB): inaudible, only for null-level fidelity; centred, "
                                             "so a plugin needs BK_PAD samples of look-ahead or the same delay on the front taps"),
                                  "S (measured, held out) / H (cause: the main taps' read)", "§11.5"),
        "controls": {
            "DELAY": lab("a pure delay ahead of everything (card: DELAY precedes the early section); the library is all DELAY 0",
                         "H", "card, §0"),
            "SPAN (AMPL slot)": lab(dict(values_s=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5]), "C (these exist) / H (finer steps)", "§1, §9.2"),
            "EARLY RT (LF slot)": lab(dict(values_s=[-3.0, -2.0, -1.0, -0.7, -0.5, -0.2, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 10.0],
                                           note="16 values = the LF column's 16 steps; 10 means flat"), "S", "§1, §9.3"),
            "RT (RT slot)": lab("blank in N (front-panel photo)", "C", "unit facts §3"),
            "LATE RT (HF slot)": lab(dict(values_s=[round(0.2 * i, 1) for i in range(11)],
                                          note="RT = LATE exactly, flat in frequency; 0 = no late; 11 values (the HF column has 10 steps in S)"),
                                     "C (these values) / H (beyond 2.0)", "§1, §10"),
            "LOW PASS": lab("applies to N as to the other programs (card: 2-12 kHz)", "C (card) / H (response)", "§0"),
        },
        "per_capture_only": "branch offset o = 1 in the plugin/reference IRs (0 = the bench frame), front-tap trim shift, clock, gain",
        "cost_per_sample": "per channel: 3 front + up to 27 main taps, 4 all-passes total (AP124 shared; AP180 L / AP188 R; AP314), 3 loops, 12 reads",
        "references": dict(files={f"{s}_{e}_{l}": f"analysis/nl/ref/nl_{s}_{e}_{l}.wav" for s, e, l in REF_SETTINGS},
                           wav_sample_rate_hz=int(round(FS)),
                           note="device rate (the exact device_rate_hz; the WAV header is the integer 26,374 Hz), float32, (L, R), "
                                "input impulse at sample 0, main-tap units; kernel-free and branch = 1. Null them with the adopted "
                                "branch read kernel bypassed, or render the port's own reference with the kernel and BK_PAD samples of "
                                "look-ahead"),
        "validation": {},
        "not_resolved": ["the remaining early-section residual after the deconvolution response and branch read kernel (§9.5, §11)",
                         "late gains below about -55 dB run 6-12 dB above the law (§10.1; §11.6 reports 3-5 dB below)",
                         "SPAN steps between 0.1, LATE beyond 2.0, DELAY, LOW PASS", "the per-capture branch offset's cause (§9.4)"],
    }


def main():
    os.makedirs(REF, exist_ok=True)
    n = int(2.5 * FS)
    x = np.zeros(n)
    x[0] = 1.0
    worst = -999.0
    for key in REF_SETTINGS:
        y = simulate(x, *key)
        m = nm.render(key, n + nm.K0, float(BRANCH), 0.0, band=False, proc=False)[nm.K0:]
        null = 10 * np.log10(np.sum((y - m) ** 2) / np.sum(y ** 2))
        worst = max(worst, null)
        print(key, f"time-domain simulate vs frequency-domain model: {null:.1f} dB")
        sf.write(os.path.join(REF, f"nl_{key[0]}_{key[1]}_{key[2]}.wav"), y.astype(np.float32), int(round(FS)), subtype="FLOAT")
    s = spec()
    s["validation"] = dict(simulate_vs_model_worst_null_db=float(worst),
                           note=f"simulate() = the bench model to {worst:.0f} dB on these settings "
                                "(nl_model.render, band=False, proc=False, branch = 1)")
    json.dump(s, open(SPEC, "w"), indent=1)
    print("wrote", SPEC)


if __name__ == "__main__":
    main()
