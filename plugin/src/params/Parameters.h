#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

// Parameter IDs, ranges and text conversions. The IDs are permanent: never rename one, and add new parameters with a
// new version number instead of changing an old one.
//
// The six joystick parameters are step indices, like the hardware's remote levers; what a step means depends on the
// program (`valueText`, `stepValue`). The tables below are PLACEHOLDERS (H) until the unit's exact step tables are
// known (docs/EMT246_unit_facts.md, open question 3); only the step counts and ranges are from the data sheet.
namespace rv::params
{
inline constexpr int version = 1; // juce::ParameterID version hint for every parameter below

namespace id
{
inline constexpr const char* program = "program";
inline constexpr const char* lowpass = "lowpass";
inline constexpr const char* delay = "delay";
inline constexpr const char* ampl = "ampl"; // N: SPAN (confirmed by the front-panel photo: SPAN in the AMPL slot)
inline constexpr const char* lf = "lf";     // N: EARLY RT (in the LF slot, RT blank, LATE in the HF slot)
inline constexpr const char* rt = "rt";
inline constexpr const char* hf = "hf";     // N: LATE RT
inline constexpr const char* mix = "mix";
inline constexpr const char* pan = "pan";
inline constexpr const char* inputGain = "inputGain";
inline constexpr const char* outputGain = "outputGain";
inline constexpr const char* power = "power"; // true = processing, false = bypassed
inline constexpr const char* link = "link";   // input and output faders move together
// The transformer-fitted unit (dsp/Transformer.h): the commercial captures' low band and guessed transformer
// distortion on the wet path; off = faithful to Tomkin and the stock unit.
inline constexpr const char* transformer = "transformer";
} // namespace id

enum class Program
{
    s,
    m,
    l,
    e,
    n,
    count
};

// The joystick columns, left to right.
enum class Joy
{
    lowpass,
    delay,
    ampl,
    lf,
    rt,
    hf,
    count
};

inline constexpr int numPrograms = (int)Program::count;
inline constexpr int numJoys = (int)Joy::count;

const char* joyId(Joy);
int joyMaxStep(Joy); // the highest step index (all start at 0)
int joyDefaultStep(Joy);

// What a joystick's value reads as: the display text (3 characters, a '.' after a digit lights its point) and the
// number behind it. `used` is false for a slot the program doesn't have (N has no RT).
struct JoyValue
{
    bool used = true;
    double value = 0.0;
    juce::String text;
    const char* unit = "";
};
JoyValue joyValue(Program, Joy, int step);

// 5-character program name for the display ("SHT.RT", "N.LIN", "E.250" with the points lit separately).
juce::String programName(Program);
juce::String programLongName(Program);

//==============================================================================
// Fader mapping: parameter value 0..1 is the position along the fader's travel, so the artwork and the host agree.
// The gain faders' scale is not linear, so the position <-> dB table is the printed scale, measured on the artwork.
float gainDbFromPosition(float position);  // below -99 dB reads as -infinity
float gainPositionFromDb(float db);
float gainLinear(float position);          // 0 at the bottom of the travel
juce::String gainText(float position);     // "-3.25 dB", "-inf dB"
inline constexpr float defaultGainPosition = 0.7431f; // 0 dB

juce::AudioProcessorValueTreeState::ParameterLayout createLayout();
} // namespace rv::params
