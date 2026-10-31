#pragma once

#include <juce_audio_basics/juce_audio_basics.h>

#include <cmath>
#include <vector>

namespace rv::dsp
{
// A stereo delay for the wet path. The delay time is slew-limited so a change glides instead of clicking; at 0 ms the
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
        currentDelay = targetDelay = 0.0f;
        writePos = 0;
    }

    void reset()
    {
        std::fill(left.begin(), left.end(), 0.0f);
        std::fill(right.begin(), right.end(), 0.0f);
        writePos = 0;
        currentDelay = targetDelay;
    }

    void setDelayMs(float ms)
    {
        const auto maxSamples = (float)juce::jmax(0, (int)left.size() - 2);
        targetDelay = (juce::jlimit(0.0f, maxSamples, ms * 0.001f * (float)sampleRate));
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

            // Slew-limited to half a sample per sample, so the read head always moves forward (a change glides as a
            // pitch bend of at most an octave, instead of running backwards).
            currentDelay += juce::jlimit(-maxSlew, maxSlew, targetDelay - currentDelay);
            const auto delay = currentDelay;
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
    static constexpr float maxSlew = 0.5f;
    float currentDelay = 0.0f, targetDelay = 0.0f;
    int writePos = 0;
};
} // namespace rv::dsp
