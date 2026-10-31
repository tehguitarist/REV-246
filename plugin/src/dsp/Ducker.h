#pragma once

#include <juce_audio_basics/juce_audio_basics.h>

#include <cmath>

namespace rv::dsp
{
// A simple input-keyed ducker for the wet path: while the input is over the threshold the wet signal is lowered by
// `rangeDb`, with one-pole attack/release smoothing on the gain. 0 ms means no smoothing (an instant step).
class Ducker
{
  public:
    void prepare(double newSampleRate)
    {
        sampleRate = newSampleRate;
        updateCoefficients();
        reset();
    }

    void reset() { gain = 1.0f; }

    void setEnabled(bool shouldDuck) { enabled = shouldDuck; }
    // A negative dB value; 0 dB is no attenuation.
    void setRangeDb(float db) { rangeGain = db >= 0.0f ? 1.0f : std::pow(10.0f, db / 20.0f); }
    void setThresholdDb(float db) { thresholdGain = std::pow(10.0f, db / 20.0f); }
    void setAttackMs(float ms)
    {
        attackMs = ms;
        attackCoefficient = coefficientFor(ms);
    }
    void setReleaseMs(float ms)
    {
        releaseMs = ms;
        releaseCoefficient = coefficientFor(ms);
    }

    // `key` is one sample per frame (the detector input); `left` and `right` are the wet channels to duck.
    void process(const float* key, float* left, float* right, int numSamples)
    {
        for (int i = 0; i < numSamples; ++i)
        {
            const auto target = enabled && std::abs(key[i]) > thresholdGain ? rangeGain : 1.0f;
            const auto coefficient = target < gain ? attackCoefficient : releaseCoefficient;
            gain += coefficient * (target - gain);
            left[i] *= gain;
            right[i] *= gain;
        }
    }

  private:
    float coefficientFor(float ms) const
    {
        return ms <= 0.0f || sampleRate <= 0.0 ? 1.0f
                                               : 1.0f - std::exp(-1.0f / ((float)sampleRate * ms * 0.001f));
    }

    void updateCoefficients()
    {
        attackCoefficient = coefficientFor(attackMs);
        releaseCoefficient = coefficientFor(releaseMs);
    }

    double sampleRate = 48000.0;
    float attackMs = 10.0f, releaseMs = 200.0f;
    float attackCoefficient = 1.0f, releaseCoefficient = 1.0f;
    float rangeGain = 0.31622777f; // -10 dB
    float thresholdGain = 0.1f;    // -20 dB
    float gain = 1.0f;
    bool enabled = false;
};
} // namespace rv::dsp
