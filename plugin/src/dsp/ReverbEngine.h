#pragma once

#include "dsp/engine/ProgramCore.h"
#include "params/Parameters.h"

#include <juce_audio_basics/juce_audio_basics.h>

#include <array>

namespace rv::dsp
{
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
