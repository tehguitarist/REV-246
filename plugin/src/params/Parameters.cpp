#include "params/Parameters.h"

#include <array>
#include <cmath>

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
constexpr float yTop = 218.f, yBottom = 1238.f;

constexpr std::array<double, 8> lowpassKHz{2, 3, 4, 5, 6, 8, 10, 12};
constexpr std::array<double, 16> lfX{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0};

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

double geometricStep(RtRange r, int step, int steps)
{
    return r.lo * std::pow(r.hi / r.lo, (double)step / (double)(steps - 1));
}

// The number as the display prints it: one decimal below 10, whole numbers from 10 up. The display right-aligns it
// in its digits, so no padding here.
juce::String threeChars(double v)
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
    static constexpr int maxima[] = {7, 999, 10, 15, 15, 9};
    return maxima[(int)j];
}

int joyDefaultStep(Joy j)
{
    // The CLA reference setting (S): LP 12, DELAY 10 ms, AMPL 80 %, LF x1.0, RT 1.0 s, HF x0.9.
    static constexpr int defaults[] = {7, 10, 8, 9, 8, 8};
    return defaults[(int)j];
}

JoyValue joyValue(Program p, Joy j, int step)
{
    JoyValue v;
    step = juce::jlimit(0, joyMaxStep(j), step);
    const bool n = p == Program::n;
    switch (j)
    {
    case Joy::lowpass:
        v.value = lowpassKHz[(size_t)step];
        v.text = threeChars(v.value);
        v.unit = "kHz";
        break;
    case Joy::delay:
        v.value = step;
        v.text = juce::String(step);
        v.unit = "ms";
        break;
    case Joy::ampl:
        v.value = n ? step * 0.05 : step * 10.0; // N: SPAN 0 to 0.5 s
        v.text = n ? threeChars(v.value) : juce::String((int)v.value);
        v.unit = n ? "s" : "%";
        break;
    case Joy::lf:
        v.value = n ? (step - 8) * 0.5 : lfX[(size_t)step]; // N: EARLY RT, negative allowed
        v.text = n ? (v.value < 0.0 ? "-" + threeChars(-v.value) : threeChars(v.value)) : threeChars(v.value);
        v.unit = n ? "s" : "x";
        break;
    case Joy::rt:
        v.used = !n;
        v.value = geometricStep(rtRange(p), step, 16);
        v.text = v.used ? threeChars(v.value) : juce::String();
        v.unit = "s";
        break;
    case Joy::hf:
        v.value = n ? geometricStep(rtRange(p), step, 10) : (step + 1) * 0.1; // N: LATE RT
        v.text = threeChars(v.value);
        v.unit = n ? "s" : "x";
        break;
    case Joy::count:
    default: break;
    }
    return v;
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
    const float y = yBottom - juce::jlimit(0.f, 1.f, position) * (yBottom - yTop);
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
            return (yBottom - y) / (yBottom - yTop);
        }
    return 0.f;
}

float gainLinear(float position)
{
    const auto db = gainDbFromPosition(position);
    return db <= -99.f ? 0.f : juce::Decibels::decibelsToGain(db);
}

juce::String gainText(float position)
{
    const auto db = gainDbFromPosition(position);
    return db <= -99.f ? juce::String(juce::CharPointer_UTF8("-\xe2\x88\x9e dB")) : juce::String(db, 2) + " dB";
}

juce::AudioProcessorValueTreeState::ParameterLayout createLayout()
{
    using namespace juce;
    AudioProcessorValueTreeState::ParameterLayout layout;
    auto pid = [](const char* s) { return ParameterID(s, version); };

    StringArray programs;
    for (int i = 0; i < numPrograms; ++i)
        programs.add(programLongName((Program)i));
    layout.add(std::make_unique<AudioParameterChoice>(pid(id::program), "Program", programs, 0));

    for (int i = 0; i < numJoys; ++i)
    {
        const auto joy = (Joy)i;
        static const char* names[] = {"Lowpass", "Reflections delay", "Reflections amplitude", "Decay LF",
                                      "Decay RT", "Decay HF"};
        layout.add(std::make_unique<AudioParameterInt>(
            pid(joyId(joy)), names[i], 0, joyMaxStep(joy), joyDefaultStep(joy), AudioParameterIntAttributes().withStringFromValueFunction([joy](int v, int) {
                const auto jv = joyValue(Program::s, joy, v);
                return jv.text.trim() + " " + jv.unit;
            })));
    }

    layout.add(std::make_unique<AudioParameterFloat>(
        pid(id::mix), "Mix", NormalisableRange<float>(0.f, 1.f), 0.25f,
        AudioParameterFloatAttributes().withStringFromValueFunction([](float v, int) { return String(v * 100.f, 2) + " %"; })));
    layout.add(std::make_unique<AudioParameterFloat>(
        pid(id::pan), "Pan", NormalisableRange<float>(-1.f, 1.f), 0.f,
        AudioParameterFloatAttributes().withStringFromValueFunction([](float v, int) {
            if (std::abs(v) < 0.005f)
                return String("C");
            return String(v < 0.f ? "L " : "R ") + String(std::abs(v) * 100.f, 2) + " %";
        })));
    for (const auto* g : {id::inputGain, id::outputGain})
        layout.add(std::make_unique<AudioParameterFloat>(
            pid(g), g == id::inputGain ? "Input" : "Output", NormalisableRange<float>(0.f, 1.f), defaultGainPosition,
            AudioParameterFloatAttributes().withStringFromValueFunction([](float v, int) { return gainText(v); })));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::power), "Power", true));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::link), "Link input/output", false));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::transformer), "Transformer", false));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::trueStereo), "True stereo", false));
    layout.add(std::make_unique<AudioParameterBool>(pid(id::sixteenBit), "16-bit processing", false));
    return layout;
}
} // namespace rv::params
