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

// Step tables. Counts and ranges are the data sheet's (C, docs/EMT246_unit_facts.md section 3); the values are the
// Tomkin library's where it has them (C: real panel steps) and guesses in between (H).
constexpr std::array<double, 8> lowpassStepsKHz{2, 3, 4, 5, 6, 8, 10, 12}; // H between the observed 2 / 8 / 10 / 12
// LF x (data sheet 0.1-2.0, 16 steps; 0.1-1.2 as before and every library LF of S / M / L: 0.1, 0.4, 1.0, 2.0 (H).
// Values above 1.2 are H except 1.5, a step on the EMT 250's front panel and the Altiverb family's fitted 1.49;
// 1.3 and 1.8 fill the gaps. Decided 30 Oct, docs/EMT250_program.md section 19).
constexpr std::array<double, 16> lfStepsX{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.8, 2.0};
constexpr int hfStepCount = 10;       // HF x 0.1 ... 1.0 in S / M / L / E (holds every library HF)
constexpr int delayMaxStepMs = 999;   // DELAY: 0 ... 999 ms in 1 ms steps
constexpr int amplMaxStep = 10;       // AMPL: 0 ... 100 % (N: span 0 ... 0.5 s) in 11 steps

// RT: 16 steps per program. The library's RT values are real panel steps (C); the range ends are the data sheet's
// (S, M, L: C; E: H); the steps between them (H) split the widest gaps in log RT, each at the geometric mean of its
// neighbours rounded to the display's precision (0.1 s below 10 s, whole seconds above), until there are 16. Library
// steps: S 0.3 0.6 0.8 1.0 1.4 1.8 2.2 2.5; M 0.6 0.8 1.2 1.8 2.6 3.0 4.0 5.0; L 2.2 2.6 3.6 4.5 6 8 10 15 20
// (analysis/ml/ml_plugin_spec.json); E the EMT 250's 16 front-panel steps 0.4-4.5 (C for the 250's panel, H that the 246's E uses
// them; docs/EMT250_program.md section 19; every Altiverb-family E step is among them). The default step is S's 1.0 s (index 8).
using RtTable = std::array<double, 16>;
constexpr RtTable rtStepsS{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5};
constexpr RtTable rtStepsM{0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.3, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0};
constexpr RtTable rtStepsL{2.0, 2.2, 2.6, 3.1, 3.6, 4.0, 4.5, 5.2, 6.0, 6.9, 8.0, 10.0, 12.0, 15.0, 17.0, 20.0};
constexpr RtTable rtStepsE{0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.3, 3.6, 4.0, 4.5};
constexpr int rtDefaultStep = 8; // 1.0 s in S

// N (analysis/nl/nl_plugin_spec.json controls). EARLY RT in the LF slot: its 16 values (S), 10 = flat. LATE RT in the
// HF slot: 0 ... 2.0 s by 0.2 (C), 0 = no late section; 11 values, one more than S / M / L / E's 10 HF steps, so the
// HF parameter runs 0 ... 10 and the other programs stop at step 9 (joyMaxStep(Program, Joy); a higher step reads as
// their step 9). SPAN in the AMPL slot: 0 ... 0.5 s by 0.05 (0.1 s steps C, the 0.05 s ones H).
constexpr std::array<double, 16> nEarlySteps{-3.0, -2.0, -1.0, -0.7, -0.5, -0.2, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 10.0};
constexpr int nLateSteps = 11;
constexpr double nLateStepS = 0.2;

const RtTable& rtTable(Program p)
{
    switch (p)
    {
    case Program::m: return rtStepsM;
    case Program::l: return rtStepsL;
    case Program::e: return rtStepsE;
    case Program::s:
    case Program::n:
    case Program::count:
    default: return rtStepsS;
    }
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
                                     (int)lfStepsX.size() - 1, (int)RtTable{}.size() - 1, nLateSteps - 1};
    return maxima[(int)j];
}

int joyMaxStep(Program p, Joy j)
{
    if (j == Joy::hf && p != Program::n)
        return hfStepCount - 1;
    return joyMaxStep(j);
}

int joyDefaultStep(Joy j)
{
    // The CLA reference setting (S): LP 12, DELAY 10 ms, AMPL 80 %, LF x1.0, RT 1.0 s, HF x0.9.
    static constexpr int defaults[] = {7, 10, 8, 9, rtDefaultStep, 8};
    return defaults[(int)j];
}

double joyNumber(Program p, Joy j, int step) noexcept
{
    step = juce::jlimit(0, joyMaxStep(p, j), step);
    const bool n = p == Program::n;
    switch (j)
    {
    case Joy::lowpass: return lowpassStepsKHz[(size_t)step];
    case Joy::delay: return step;
    case Joy::ampl: return n ? step * 0.05 : step * 10.0; // N: SPAN 0 to 0.5 s
    case Joy::lf: return n ? nEarlySteps[(size_t)step] : lfStepsX[(size_t)step]; // N: EARLY RT
    case Joy::rt: return n ? 0.0 : rtTable(p)[(size_t)step];                    // N has no RT
    case Joy::hf: return n ? step * nLateStepS : (step + 1) * 0.1;               // N: LATE RT
    case Joy::count:
    default: return 0.0;
    }
}

JoyValue joyValue(Program p, Joy j, int step)
{
    JoyValue v;
    step = juce::jlimit(0, joyMaxStep(p, j), step);
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

int joyNearestStep(Program p, Joy j, double value)
{
    // RT is spaced roughly geometrically, so it is compared in log; every other column linearly.
    const bool logScale = j == Joy::rt && p != Program::n && value > 0.0;
    int best = 0;
    auto bestDistance = std::numeric_limits<double>::max();
    for (int step = 0; step <= joyMaxStep(p, j); ++step)
    {
        const auto v = joyNumber(p, j, step);
        const auto distance = logScale ? std::abs(std::log(v / value)) : std::abs(v - value);
        if (distance < bestDistance)
        {
            bestDistance = distance;
            best = step;
        }
    }
    return best;
}

std::optional<int> joyStepFromValue(Program p, Joy j, double value)
{
    if (!std::isfinite(value))
        return std::nullopt;
    if (j == Joy::rt && p == Program::n)
        return std::nullopt; // RT is blank in the non-linear program

    // RT's steps are not a table the user can recite, so accept any value in range and use the nearest step.
    // Everything else has a visible, exact step (whole ms, 10 %, 0.1 x, the lowpass table, N's EARLY / LATE ...).
    if (j == Joy::rt)
    {
        const auto& t = rtTable(p);
        if (value < t.front() - 1e-9 || value > t.back() + 1e-9)
            return std::nullopt;
        return joyNearestStep(p, j, value);
    }

    constexpr double tolerance = 1e-4; // display precision is at most two decimals
    for (int step = 0; step <= joyMaxStep(p, j); ++step)
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
    case Joy::lf:
        if (n)
        {
            juce::String valid;
            for (size_t i = 0; i < nEarlySteps.size(); ++i)
            {
                if (i > 0)
                    valid << (i + 1 == nEarlySteps.size() ? " or " : ", ");
                valid << displayNumber(nEarlySteps[i]);
            }
            return "Valid: " + valid + " s (10 = flat).";
        }
        return "Valid: 0.1 to 1.2 x in 0.1 steps, then 1.3, 1.5, 1.8 and 2.0 x.";
    case Joy::rt:
        if (n)
            return "RT is not used in the non-linear program.";
        return "Valid: " + displayNumber(rtTable(p).front()) + " to " + displayNumber(rtTable(p).back())
               + " s; any value in that range snaps to the nearest of " + juce::String(joyMaxStep(p, j) + 1)
               + " steps.";
    case Joy::hf:
        if (n)
            return "Valid: 0 to 2.0 s in steps of 0.2 s (0 = no late reverb).";
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

void migrateState(juce::ValueTree state)
{
    const int from = state.getProperty(stateVersionProperty, 0);
    if (from > stateVersion)
        return;
    // from 0: nothing to change (0 and 1 share the layout). A future step goes here, e.g.
    //   if (from < 2) { ... remap the old steps of the parameter that changed ... }
    state.setProperty(stateVersionProperty, stateVersion, nullptr);
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
