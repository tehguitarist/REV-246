#pragma once

#include "params/Parameters.h"

#include <juce_audio_basics/juce_audio_basics.h>

namespace rv::dsp
{
// The reverb's settings as the panel shows them, already converted from step indices to values.
struct Settings
{
    params::Program program = params::Program::s;
    double lowpassKHz = 12.0;
    double delayMs = 10.0;
    double amplPercent = 80.0; // N: span in seconds in `ampl`'s place, see params::joyValue
    double lfX = 1.0;
    double rt = 1.0;
    double hfX = 0.9;
    bool sixteenBit = false; // quantise the loop storage to 16 bits (the real algorithm); no-op in the placeholder
    std::array<int, params::numJoys> steps{}; // the raw step indices, for engines that want the table index
};

// The EMT 246 model. This is a PLACEHOLDER with the final interface: mono in, stereo out, at the host rate. The real
// algorithm (docs/model_s_results.md, docs/ml_structure.md) is being worked out in tools/ and will replace
// `process`; until then it passes the mono input to both outputs so the panel, meters and mix path can be tried.
// `process` must stay allocation- and lock-free.
class ReverbEngine
{
  public:
    void prepare(double sampleRate, int maxBlockSize)
    {
        rate = sampleRate;
        (void)maxBlockSize;
        reset();
    }
    void reset() {}
    void setSettings(const Settings& s) { settings = s; }

    // `in` is the mono input; writes `n` samples of the wet signal.
    void process(const float* in, float* outL, float* outR, int n)
    {
        for (int i = 0; i < n; ++i)
            outL[i] = outR[i] = in[i];
    }

  private:
    double rate = 48000.0;
    Settings settings;
};
} // namespace rv::dsp
