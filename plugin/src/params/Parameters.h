#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

#include <optional>

// Parameter IDs, ranges and text conversions. The IDs are permanent: never rename one, and add new parameters with a
// new version number instead of changing an old one.
//
// The six joystick parameters are step indices, like the hardware's remote levers; what a step means depends on the
// program (`joyValue`, `joyNumber`). The step counts and ranges are the data sheet's; the values are the Tomkin
// library's where it has them (real panel steps, C) and guesses between them (H): see the tables in Parameters.cpp.
namespace rv::params
{
inline constexpr int parameterVersion = 1; // juce::ParameterID version hint for every parameter below
// Parameters added after the first release use a new version hint, so a host never confuses one with an old ID.
inline constexpr int addedParameterVersion = 2;

// The layout of the saved state (sessions and preset files), written as the root's `stateVersion`. Bump it whenever a
// stored value changes meaning (a step table's range or order, a parameter's unit) and teach migrateState() the step.
// 0 = written before versioning (the pre-release builds, the same layout as 1).
inline constexpr int stateVersion = 2; // 2: pre-delay and the ducker were added
inline const juce::Identifier stateVersionProperty{"stateVersion"};
// Brings a state written by an older build up to stateVersion, in place (a ValueTree is a shared reference: pass a
// copy you own). A newer state is left as it is: every parameter is matched by id, unknown ones are ignored.
void migrateState(juce::ValueTree state);

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
// Stereo in / stereo out only: false (default) sums the input to the unit's mono input ("ORIGINAL"), true processes
// left and right through separate reverbs. A routing preference, so presets leave it alone like link and transformer.
inline constexpr const char* trueStereo = "trueStereo";
// True (default) = the original's 16-bit storage of the reverb loops; false = full-precision (float) processing. A
// sound setting, so unlike the transformer it is part of a preset.
inline constexpr const char* sixteenBit = "sixteenBit";
// Wet stereo width: 0 = mono, 100 = unchanged, 200 = extra wide. A sound setting; presets store it.
inline constexpr const char* stereoWidth = "stereoWidth";
// An extra delay before the reverb's wet signal (0-200 ms), after the input-keyed duck.
inline constexpr const char* preDelayMs = "preDelayMs";
// Input-keyed ducker on the wet path: lowers the reverb while the dry input is over the threshold.
inline constexpr const char* duckEnabled = "duckEnabled";
inline constexpr const char* duckRangeDb = "duckRangeDb";         // attenuation while ducked, -60 ... 0 dB
inline constexpr const char* duckThresholdDb = "duckThresholdDb"; // -80 ... 0 dB
inline constexpr const char* duckAttackMs = "duckAttackMs";       // 0 ... 1000 ms
inline constexpr const char* duckReleaseMs = "duckReleaseMs";     // 0 ... 1000 ms
} // namespace id

inline constexpr float preDelayMaxMs = 200.0f;
inline constexpr float duckMinRangeDb = -60.0f, duckMaxRangeDb = 0.0f, duckDefaultRangeDb = -10.0f;
inline constexpr float duckMinThresholdDb = -80.0f, duckMaxThresholdDb = 0.0f, duckDefaultThresholdDb = -20.0f;
inline constexpr float duckMaxTimeMs = 1000.0f;
inline constexpr float duckDefaultAttackMs = 10.0f, duckDefaultReleaseMs = 200.0f;

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
int joyMaxStep(Joy);          // the parameter's highest step index (all start at 0), over every program
int joyMaxStep(Program, Joy); // this program's highest step (HF: 9 in S / M / L / E, 10 in N); higher reads as it
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
// Just the number of joyValue(...).value: allocation-free, so the audio thread can use it (joyValue builds a String).
double joyNumber(Program, Joy, int step) noexcept;

// Manual entry (double-click a display field). The value text is in the field's unit. Discrete controls accept only an
// exact step; RT (S / M / L / E) accepts any value in range and picks the nearest step. Nullopt means the
// value is not usable (including RT in the non-linear program, which has no RT).
std::optional<int> joyStepFromValue(Program, Joy, double value);
// The program's step nearest to `value` (RT in log, the rest linearly), for any value.
int joyNearestStep(Program, Joy, double value);
juce::String joyInputName(Program, Joy); // "lowpass frequency in kHz", etc.
juce::String joyInputHint(Program, Joy); // the valid values / range for the dialog

// 5-character program name for the display ("SHT.RT", "N.LIN", "E.250" with the points lit separately).
juce::String programName(Program);
juce::String programLongName(Program);

//==============================================================================
// Fader mapping: parameter value 0..1 is the position along the fader's travel, so the artwork and the host agree.
// The gain faders' scale is not linear, so the position <-> dB table is the printed scale, measured on the artwork.
inline constexpr float gainOffDb = -99.f;  // at or below this the fader reads -infinity (linear gain 0)
float gainDbFromPosition(float position);
float gainPositionFromDb(float db);
float gainLinearFromPosition(float position); // 0 at the bottom of the travel
juce::String gainText(float position);        // "-3.25 dB", "-inf dB"
inline constexpr float defaultGainPosition = 0.74313725f; // 0 dB

juce::AudioProcessorValueTreeState::ParameterLayout createLayout();
} // namespace rv::params
