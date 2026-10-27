#pragma once

#include "params/Parameters.h"

#include <juce_core/juce_core.h>

// The text of the panel's tooltips (shown while the "?" is on). Each function returns what a control does in the state
// it is in, so the same lever explains itself differently in S, M, L, E and N, and the stereo controls say why they
// are unavailable on a mono bus. The numbers that belong to the parameters (ranges, steps) come from
// params::joyInputHint, so they cannot drift from what the lever does.
namespace rv::ui::help
{
// The bus layout and routing the stereo controls depend on.
struct Routing
{
    bool stereoToStereo = true; // a stereo input and a stereo output: the stereo button is available
    bool stereoOut = true;      // a stereo output: width and pan are available
    bool trueStereo = false;    // the TRUE ST choice is active
};

inline juce::String program(params::Program p)
{
    using Program = params::Program;
    switch (p)
    {
    case Program::s:
        return "S, short reverb: small rooms and booths. The densest echo build-up of the three, RT 0.2 to 2.5 s "
               "(best between 0.2 and 1 s).";
    case Program::m:
        return "M, medium reverb: the all-round concert hall. RT 0.5 to 5 s (best between 0.8 and 4 s).";
    case Program::l:
        return "L, long reverb: the cathedral. Widely spaced early reflections and a lower echo density, RT 2 to 20 s "
               "(best between 3 and 10 s).";
    case Program::e:
        return "E, the EMT 250 program: the older EMT 250 algorithm (recirculating delay loops with slowly moving "
               "taps) with extended pre-delay. The controls are those of S, M and L.";
    case Program::n:
        return "N, non-linear: DELAY, then an early section lasting SPAN with its own decay (EARLY RT; a negative "
               "value makes it rise, for reverse reverb), then the late reverb (LATE RT). A short SPAN cuts the early "
               "section off abruptly, which gives gated reverb. RT, LF and HF are not used.";
    case Program::count:
    default: return {};
    }
}

// The lever, the funnel above it and its display all describe the same parameter, in the current program.
inline juce::String joy(params::Program p, params::Joy j)
{
    using Joy = params::Joy;
    const bool n = p == params::Program::n;
    juce::String what;
    switch (j)
    {
    case Joy::lowpass:
        what = "Lowpass: a low-pass filter on the whole reverb output, after the reflections. Lower settings give a "
               "darker, softer reverb.";
        break;
    case Joy::delay:
        what = n ? "Delay: the pre-delay before the early section starts."
                 : "Reflections delay (pre-delay): the time between the dry sound and the first reflections.";
        break;
    case Joy::ampl:
        what = n ? "Span: how long the early section lasts. When it ends the late reverb takes over, abruptly if the "
                   "early decay is still loud."
                 : "Reflections amplitude: the level of the burst of early reflections against the start of the "
                   "reverb. 0 leaves the reverb alone.";
        break;
    case Joy::lf:
        what = n ? "Early RT: the decay time of the early section, which can be negative: the section then rises "
                   "instead of falling (reverse reverb)."
                 : "Decay LF: a multiplier on the decay time at low frequencies (RT at 50 Hz = RT x LF). 1.0 is "
                   "flat; above 1 the bass rings longer, below 1 it dies away sooner.";
        break;
    case Joy::rt:
        what = n ? "RT is not used in the non-linear program: set EARLY RT and LATE RT instead."
                 : "Decay RT: the reverberation time at 1 kHz, the time the tail takes to die away.";
        break;
    case Joy::hf:
        what = n ? "Late RT: the decay time of the late reverb that follows the early section."
                 : "Decay HF: a multiplier on the decay time at high frequencies (RT at 7 kHz = RT x HF). 1.0 adds "
                   "no damping; lower is darker, as air absorbs the top end of a real room.";
        break;
    case Joy::count:
    default: return {};
    }
    if (j == params::Joy::rt && n)
        return what;
    return what + " " + params::joyInputHint(p, j)
           + " Lever: press up or down, hold to repeat (faster after a second). Funnel: click or drag. "
                 "Option/Alt-click resets.";
}

inline juce::String power(bool on)
{
    return on ? "Power: on. Press to bypass: the input passes through unchanged (delayed by the plugin's latency, so "
                "bypassing does not shift the timing)."
              : "Power: bypassed. Press to switch the reverb on.";
}

inline juce::String link()
{
    return "Link (default on): INPUT and OUTPUT move together, the other by the same number of dB in the opposite "
           "direction, so the overall level stays put while you change how hard the reverb is driven.";
}

inline juce::String input()
{
    return "Input level, +12 dB to -infinity. It sets how hard the reverb is driven (and what the VU meter reads, and "
           "where the 16-bit A/D clips); the dry signal in the MIX is not affected. Double-click or Option/Alt-click "
           "for 0 dB; Shift for fine.";
}

inline juce::String output()
{
    return "Output level, +12 dB to -infinity, applied to the dry and reverb together after MIX. Double-click or "
           "Option/Alt-click for 0 dB; Shift for fine.";
}

inline juce::String mix()
{
    return "Mix: the balance of dry signal to reverb, 0 % (all dry) to 100 % (all reverb, the default for use on a "
           "send or aux). Double-click or Option/Alt-click for 100 %; Shift for fine.";
}

inline juce::String pan(const Routing& r)
{
    if (!r.stereoOut)
        return "Pan: has no effect with a mono output, which carries the dry signal in the centre.";
    return juce::String("Pan: where the dry signal sits in the stereo field; the reverb is not moved (use WIDTH for "
                        "that). ")
           + (r.trueStereo ? "In TRUE ST it balances the dry left and right channels. " : "")
           + "Double-click or Option/Alt-click to centre; Shift for fine.";
}

inline juce::String width(const Routing& r)
{
    if (!r.stereoOut)
        return "Width: needs a stereo output. With a mono output the reverb is a single channel.";
    return "Width of the reverb (the dry signal is not affected): 0 = mono, 100 = as the unit makes it, 200 = extra "
           "wide. Above 100 the bass stays mono. Option/Alt-click for 100; double-click to type; Shift for fine.";
}

inline juce::String stereoMode(const Routing& r)
{
    if (!r.stereoToStereo)
        return "ORIGINAL / TRUE ST: needs a stereo input and a stereo output. Here the plugin runs as the unit "
               "itself does: one mono input into the reverb.";
    if (r.trueStereo)
        return "TRUE ST (current): the left and right inputs each run through their own reverb, so the input's stereo "
               "image is kept (twice the CPU). The dry signal keeps its channels. Click for Original.";
    return "Original (current): the stereo input is summed to mono at -6 dB and fed to one reverb, then the reverb's "
           "left and right outputs are heard, as in the unit (mono in, stereo out). Click for TRUE ST, which runs "
           "a separate reverb for each input channel.";
}
} // namespace rv::ui::help
