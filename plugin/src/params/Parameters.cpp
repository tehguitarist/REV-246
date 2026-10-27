#include "params/Parameters.h"

#include <array>
#include <cmath>
#include <limits>

namespace rv::params
{
namespace
{
// Tick positions measured on ui/main.png (design pixels): the cap centre at each printed value.
struct Tick
{
    float db, y;
};
constexpr std::array<Tick, 11> gainTicks{{{12.f, 218.f},
                                          {6.f, 350.f},
                                          {0.f, 480.f},
                                          {-5.f, 593.f},
                                          {-10.f, 688.f},
                                          {-15.f, 782.f},
                                          {-20.f, 872.f},
                                          {-30.f, 982.f},
                                          {-40.f, 1076.f},
                                          {-50.f, 1160.f},
                                          {-100.f, 1238.f}}}; // the bottom is -infinity
constexpr float faderTopY = 218.f, faderBottomY = 1238.f; // the travel's ends in the same design pixels

// Step tables (H until the unit's exact tables are known; the counts and ranges are from the data sheet).
constexpr std::array<double, 8> lowpassStepsKHz{2, 3, 4, 5, 6, 8, 10, 12};
constexpr std::array<double, 16> lfStepsX{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0};
constexpr int rtStepCount = 16;       // RT, all programs
constexpr int hfStepCount = 10;       // HF x 0.1 ... 1.0; also N's late RT
constexpr int delayMaxStepMs = 999;   // DELAY: 0 ... 999 ms in 1 ms steps
constexpr int amplMaxStep = 10;       // AMPL: 0 ... 100 % (N: span 0 ... 0.5 s) in 11 steps
constexpr int nEarlyRtZeroStep = 8;   // N's EARLY RT runs -4.0 ... +3.5 s in 0.5 s steps, zero at this step

struct RtRange
{
    double lo, hi;
};
// S, M and L from the data sheet (16 steps); E and N are guesses.
constexpr RtRange rtRange(Program p)
{
    switch (p)
    {
    case Program::s: return {0.2, 2.5};
    case Program::m: return {0.5, 5.0};
    case Program::l: return {2.0, 20.0};
    case Program::e: return {0.4, 5.0};
    case Program::n:
    case Program::count:
    default: return {0.5, 10.0}; // N: LATE RT, on the HF column
    }
}

// Step `step` of `steps` geometrically spaced values from r.lo to r.hi.
double geometricStep(RtRange r, int step, int steps)
{
    return r.lo * std::pow(r.hi / r.lo, (double)step / (double)(steps - 1));
}

// The number as the display prints it: one decimal below 10, whole numbers from 10 up. The display right-aligns it
// in its digits, so no padding here.
juce::String displayNumber(double v)
{
    const auto rounded = std::round(v * 10.0) / 10.0;
    return rounded < 9.95 ? juce::String(rounded, 1) : juce::String((int)std::round(rounded));
}
} // namespace

const char* joyId(Joy j)
{
    static constexpr const char* ids[] = {id::lowpass, id::delay, id::ampl, id::lf, id::rt, id::hf};
    return ids[(int)j];
}

int joyMaxStep(Joy j)
{
    static constexpr int maxima[] = {(int)lowpassStepsKHz.size() - 1, delayMaxStepMs, amplMaxStep,
                                     (int)lfStepsX.size() - 1, rtStepCount - 1, hfStepCount - 1};
    return maxima[(int)j];
}

int joyDefaultStep(Joy j)
{
    // The CLA reference setting (S): LP 12, DELAY 10 ms, AMPL 80 %, LF x1.0, RT about 1 s, HF x0.9. RT step 8 reads
    // 0.8 s on the placeholder geometric table (1.0 s falls between steps 9 and 10).
    static constexpr int defaults[] = {7, 10, 8, 9, 8, 8};
    return defaults[(int)j];
}

double joyNumber(Program p, Joy j, int step) noexcept
{
    step = juce::jlimit(0, joyMaxStep(j), step);
    const bool n = p == Program::n;
    switch (j)
    {
    case Joy::lowpass: return lowpassStepsKHz[(size_t)step];
    case Joy::delay: return step;
    case Joy::ampl: return n ? step * 0.05 : step * 10.0; // N: SPAN 0 to 0.5 s
    case Joy::lf: return n ? (step - nEarlyRtZeroStep) * 0.5 : lfStepsX[(size_t)step]; // N: EARLY RT, negative allowed
    case Joy::rt: return geometricStep(rtRange(p), step, rtStepCount);
    case Joy::hf: return n ? geometricStep(rtRange(p), step, hfStepCount) : (step + 1) * 0.1; // N: LATE RT
    case Joy::count:
    default: return 0.0;
    }
}

JoyValue joyValue(Program p, Joy j, int step)
{
    JoyValue v;
    step = juce::jlimit(0, joyMaxStep(j), step);
    v.value = joyNumber(p, j, step);
    const bool n = p == Program::n;
    switch (j)
    {
    case Joy::lowpass:
        v.text = displayNumber(v.value);
        v.unit = "kHz";
        break;
    case Joy::delay:
        v.text = juce::String(step);
        v.unit = "ms";
        break;
    case Joy::ampl:
        v.text = n ? displayNumber(v.value) : juce::String((int)v.value);
        v.unit = n ? "s" : "%";
        break;
    case Joy::lf:
        v.text = v.value < 0.0 ? "-" + displayNumber(-v.value) : displayNumber(v.value);
        v.unit = n ? "s" : "x";
        break;
    case Joy::rt:
        v.used = !n;
        v.text = v.used ? displayNumber(v.value) : juce::String();
        v.unit = "s";
        break;
    case Joy::hf:
        v.text = displayNumber(v.value);
        v.unit = n ? "s" : "x";
        break;
    case Joy::count:
    default: break;
    }
    return v;
}

std::optional<int> joyStepFromValue(Program p, Joy j, double value)
{
    if (!std::isfinite(value))
        return std::nullopt;
    if (j == Joy::rt && p == Program::n)
        return std::nullopt; // RT is blank in the non-linear program

    // The geometric decays are not a table the user can recite, so accept any value in range and use the nearest
    // step. Everything else has a visible, exact step (whole ms, 10 %, 0.1 x, the lowpass table, ...).
    const bool nearest = j == Joy::rt || (j == Joy::hf && p == Program::n);
    if (nearest)
    {
        const auto range = rtRange(p);
        if (value < range.lo - 1e-9 || value > range.hi + 1e-9)
            return std::nullopt;
        int best = 0;
        auto bestDistance = std::numeric_limits<double>::max();
        for (int step = 0; step <= joyMaxStep(j); ++step)
        {
            const auto distance = std::abs(joyNumber(p, j, step) - value);
            if (distance < bestDistance)
            {
                bestDistance = distance;
                best = step;
            }
        }
        return best;
    }

    constexpr double tolerance = 1e-4; // display precision is at most two decimals
    for (int step = 0; step <= joyMaxStep(j); ++step)
        if (std::abs(joyNumber(p, j, step) - value) <= tolerance)
            return step;
    return std::nullopt;
}

juce::String joyInputName(Program p, Joy j)
{
    const bool n = p == Program::n;
    switch (j)
    {
    case Joy::lowpass: return "lowpass frequency in kHz";
    case Joy::delay: return "reflections delay in ms";
    case Joy::ampl: return n ? "span in s" : "reflections amplitude in %";
    case Joy::lf: return n ? "early RT in s" : "decay LF multiplier (x)";
    case Joy::rt: return n ? "RT" : "decay RT in s";
    case Joy::hf: return n ? "late RT in s" : "decay HF multiplier (x)";
    case Joy::count:
    default: return {};
    }
}

juce::String joyInputHint(Program p, Joy j)
{
    const bool n = p == Program::n;
    switch (j)
    {
    case Joy::lowpass:
    {
        juce::String valid;
        for (int step = 0; step <= joyMaxStep(j); ++step)
        {
            if (step > 0)
                valid << (step == joyMaxStep(j) ? " or " : ", ");
            valid << juce::String((int)lowpassStepsKHz[(size_t)step]);
        }
        return "Valid: " + valid + " kHz.";
    }
    case Joy::delay: return "Valid: whole numbers from 0 to 999 ms.";
    case Joy::ampl: return n ? "Valid: 0 to 0.5 s in steps of 0.05 s."
                             : "Valid: 0 to 100 % in steps of 10 %.";
    case Joy::lf: return n ? "Valid: -4.0 to 3.5 s in steps of 0.5 s."
                           : "Valid: 0.1 to 1.2 x in 0.1 steps, then 1.4, 1.6, 1.8 and 2.0 x.";
    case Joy::rt:
        if (n)
            return "RT is not used in the non-linear program.";
        return "Valid: " + displayNumber(rtRange(p).lo) + " to " + displayNumber(rtRange(p).hi)
               + " s; any value in that range snaps to the nearest of " + juce::String(joyMaxStep(j) + 1)
               + " steps.";
    case Joy::hf:
        if (n)
            return "Valid: " + displayNumber(rtRange(p).lo) + " to " + displayNumber(rtRange(p).hi)
                   + " s; any value in that range snaps to the nearest of " + juce::String(joyMaxStep(j) + 1)
                   + " steps.";
        return "Valid: 0.1 to 1.0 in steps of 0.1 (x).";
    case Joy::count:
    default: return {};
    }
}

juce::String programName(Program p)
{
    static const char* names[] = {"SHT.RT", "MED.RT", "LNG.RT", "E.250", "N.LIN"};
    return names[(int)p];
}

juce::String programLongName(Program p)
{
    static const char* names[] = {"S: short reverb", "M: medium reverb", "L: long reverb", "E: 250 mode",
                                  "N: non-linear"};
    return names[(int)p];
}

//==============================================================================
float gainDbFromPosition(float position)
{
    const float y = faderBottomY - juce::jlimit(0.f, 1.f, position) * (faderBottomY - faderTopY);
    for (size_t i = 1; i < gainTicks.size(); ++i)
        if (y <= gainTicks[i].y)
        {
            const auto& a = gainTicks[i - 1];
            const auto& b = gainTicks[i];
            return a.db + (b.db - a.db) * (y - a.y) / (b.y - a.y);
        }
    return gainTicks.back().db;
}

float gainPositionFromDb(float db)
{
    db = juce::jlimit(gainTicks.back().db, gainTicks.front().db, db);
    for (size_t i = 1; i < gainTicks.size(); ++i)
        if (db >= gainTicks[i].db)
        {
            const auto& a = gainTicks[i - 1];
            const auto& b = gainTicks[i];
            const float y = a.y + (b.y - a.y) * (db - a.db) / (b.db - a.db);
            return (faderBottomY - y) / (faderBottomY - faderTopY);
        }
    return 0.f;
}

float gainLinearFromPosition(float position)
{
    const auto db = gainDbFromPosition(position);
    return db <= gainOffDb ? 0.f : juce::Decibels::decibelsToGain(db);
}

juce::String gainText(float position)
{
    const auto db = gainDbFromPosition(position);
    return db <= gainOffDb ? juce::String(juce::CharPointer_UTF8("-\xe2\x88\x9e dB")) : juce::String(db, 2) + " dB";
}

juce::AudioProcessorValueTreeState::ParameterLayout createLayout()
{
    using namespace juce;
    AudioProcessorValueTreeState::ParameterLayout layout;
    auto pid = [](const char* s) { return ParameterID(s, parameterVersion); };

    StringArray programs;
    for (int i = 0; i < numPrograms; ++i)
        programs.add(programLongName((Program)i));
    layout.add(std::make_unique<AudioParameterChoice>(pid(id::program), "Program", programs, 0));

    // Indexed by Joy. The text is the value in program S's reading; the host cannot see the current program here.
    static const char* const joyNames[] = {"Lowpass", "Reflections delay", "Reflections amplitude", "Decay LF",
                                           "Decay RT", "Decay HF"};
    for (int i = 0; i < numJoys; ++i)
    {
        const auto joy = (Joy)i;
        layout.add(std::make_unique<AudioParameterInt>(
            pid(joyId(joy)), joyNames[i], 0, joyMaxStep(joy), joyDefaultStep(joy),
            AudioParameterIntAttributes().withStringFromValueFunction([joy](int v, int) {
                const auto jv = joyValue(Program::s, joy, v);
                return jv.text.trim() + " " + jv.unit;
            })));
    }

    layout.add(std::make_unique<AudioParameterFloat>(
        pid(id::mix), "Mix", NormalisableRange<float>(0.f, 1.f), 1.f,
        AudioParameterFloatAttributes().withStringFromValueFunction([](float v, int) { return String(v * 100.f, 2) + " %"; })));
    layout.add(std::make_unique<AudioParameterFloat>(
        pid(id::pan), "Pan", NormalisableRange<float>(-1.f, 1.f), 0.f,
        AudioParameterFloatAttributes().withStringFromValueFunction([](float v, int) {
            if (std::abs(v) < 0.005f)
                return String("C");
            return String(v < 0.f ? "L " : "R ") + String(std::abs(v) * 100.f, 2) + " %";
        })));
    layout.add(std::make_unique<AudioParameterFloat>(
        pid(id::stereoWidth), "Stereo width", NormalisableRange<float>(0.f, 200.f), 100.f,
        AudioParameterFloatAttributes().withStringFromValueFunction(
            [](float v, int) { return String(juce::roundToInt(v)) + " %"; })));
    // The faders' value is the position along the travel (0 ... 1); the text is the printed dB scale.
    struct GainFader
    {
        const char* paramId;
        const char* name;
    };
    for (const auto& g : {GainFader{id::inputGain, "Input"}, GainFader{id::outputGain, "Output"}})
        layout.add(std::make_unique<AudioParameterFloat>(
            pid(g.paramId), g.name, NormalisableRange<float>(0.f, 1.f), defaultGainPosition,
            AudioParameterFloatAttributes().withStringFromValueFunction([](float v, int) { return gainText(v); })));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::power), "Power", true));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::link), "Link input/output", true));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::transformer), "Transformer", true));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::trueStereo), "True stereo", false));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::sixteenBit), "16-bit processing", true));
    return layout;
}
} // namespace rv::params
