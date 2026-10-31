#pragma once

#include <juce_audio_basics/juce_audio_basics.h>

#include <cmath>
#include <vector>

namespace rv::dsp
{
// A stereo delay for the wet path. The delay time is smoothed so a change glides instead of clicking; at 0 ms the
// output is the input exactly (no extra latency).
class PreDelay
{
  public:
    static constexpr float maxDelayMs = 200.0f;

    void prepare(double newSampleRate)
    {
        sampleRate = newSampleRate;
        const auto maxSamples = juce::jmax(2, (int)std::ceil(sampleRate * maxDelayMs * 0.001) + 2);
        left.assign((size_t)maxSamples, 0.0f);
        right.assign((size_t)maxSamples, 0.0f);
        delaySamples.reset(sampleRate, 0.03);
        delaySamples.setCurrentAndTargetValue(0.0f);
        writePos = 0;
    }

    void reset()
    {
        std::fill(left.begin(), left.end(), 0.0f);
        std::fill(right.begin(), right.end(), 0.0f);
        writePos = 0;
        delaySamples.setCurrentAndTargetValue(delaySamples.getTargetValue());
    }

    void setDelayMs(float ms)
    {
        const auto maxSamples = (float)juce::jmax(0, (int)left.size() - 2);
        delaySamples.setTargetValue(juce::jlimit(0.0f, maxSamples, ms * 0.001f * (float)sampleRate));
    }

    void process(float* wetLeft, float* wetRight, int numSamples)
    {
        const auto size = (int)left.size();
        if (size < 2)
            return;
        for (int i = 0; i < numSamples; ++i)
        {
            left[(size_t)writePos] = wetLeft[i];
            right[(size_t)writePos] = wetRight[i];

            const auto delay = delaySamples.getNextValue();
            auto readPos = (float)writePos - delay;
            while (readPos < 0.0f)
                readPos += (float)size;
            const auto index = (int)readPos;
            const auto fraction = readPos - (float)index;
            const auto next = (index + 1) % size;

            wetLeft[i] = left[(size_t)index] + fraction * (left[(size_t)next] - left[(size_t)index]);
            wetRight[i] = right[(size_t)index] + fraction * (right[(size_t)next] - right[(size_t)index]);

            writePos = (writePos + 1) % size;
        }
    }

  private:
    double sampleRate = 48000.0;
    std::vector<float> left, right;
    juce::SmoothedValue<float> delaySamples;
    int writePos = 0;
};
} // namespace rv::dsp
