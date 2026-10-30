#pragma once

#include <juce_audio_basics/juce_audio_basics.h>
#include <juce_dsp/juce_dsp.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace rv::dsp
{
// Stereo width, 0..200 (100 = untouched, zero CPU).
//
// No late-tail decorrelator (31 October 2026): the unit's late L / R are already uncorrelated (|r| <= 0.04 over
// 0.1-1 s in every program), so one adds no width, and the design's Schroeder all-pass version delayed the tail ~8 ms
// against the early reflections when engaged. The M/S gains below move nothing in time, so the early reflections are
// never smeared (docs/DEVELOPMENT.md, Parameters).
//
// StereoWidener is level-matched for equal Mid/Side energy:
//   gm = sqrt2 * cos(theta), gs = sqrt2 * sin(theta)
//   width 0..100: tan(theta) = width/100 (0 = mono, 100 = 1:1)
//   width 100..200: side/mid ratio rises to 6 dB (mid -4 dB, side +2 dB at 200)
// Above 100 the low band's angle is pulled towards mono as the width rises.
class StereoWidener
{
  public:
    void prepare(double sampleRate, float crossoverHz = 150.0f)
    {
        const auto fc = juce::jmin(crossoverHz, (float)(sampleRate * 0.45));
        const auto lp = juce::dsp::IIR::Coefficients<float>::makeLowPass(sampleRate, fc, kButterworthQ);
        const auto hp = juce::dsp::IIR::Coefficients<float>::makeHighPass(sampleRate, fc, kButterworthQ);

        mLo.set(lp);
        sLo.set(lp);
        mHi.set(hp);
        sHi.set(hp);

        width.reset(sampleRate, 0.03);
        width.setCurrentAndTargetValue(100.0f);
        engage.reset(sampleRate, 0.02);
        engage.setCurrentAndTargetValue(0.0f);
        resetFilters();
    }

    /** 0..200, 100 = bypass. `snap` sets the current value too (for prepareToPlay, not for live changes). */
    void setWidth(float w, bool snap = false) noexcept
    {
        w = juce::jlimit(0.0f, 200.0f, w);
        if (snap)
        {
            width.setCurrentAndTargetValue(w);
            engage.setCurrentAndTargetValue(w > 100.0f ? 1.0f : 0.0f);
        }
        else
            width.setTargetValue(w);
    }

    void process(float* l, float* r, int n) noexcept
    {
        if (!width.isSmoothing() && juce::exactlyEqual(width.getTargetValue(), 100.0f) && !engage.isSmoothing()
            && engage.getCurrentValue() <= 0.0f)
        {
            releaseFilters();
            return;
        }

        int pos = 0;
        while (pos < n)
        {
            const auto c = std::min(kChunk, n - pos);
            const auto w = width.skip(c);
            engage.setTargetValue(w > 100.0f || width.getTargetValue() > 100.0f ? 1.0f : 0.0f);
            processChunk(l + pos, r + pos, c, w, engage.skip(c));
            pos += c;
        }
    }

  private:
    static constexpr int kChunk = 32;          // control-rate update interval (samples)
    static constexpr float kSqrt2 = juce::MathConstants<float>::sqrt2;
    static constexpr float kButterworthQ = 0.70710678f; // 1/sqrt2: two of these in cascade make an LR4

    struct LR4 // Linkwitz-Riley 4th order = two cascaded Butterworth biquads
    {
        juce::dsp::IIR::Filter<float> a, b;

        void set(juce::dsp::IIR::Coefficients<float>::Ptr c)
        {
            a.coefficients = c;
            b.coefficients = c;
        }

        void reset()
        {
            a.reset();
            b.reset();
        }

        float process(float x) noexcept { return b.processSample(a.processSample(x)); }
    };

    // The Mid/Side gains for one width: the high band (all of it on the broadband path) and the low band.
    struct Gains
    {
        float gmH = 1.0f, gsH = 1.0f, gmL = 1.0f, gsL = 1.0f;
    };

    static Gains gainsFor(float w) noexcept
    {
        Gains g;
        const auto u = std::max(0.0f, (w - 100.0f) / 100.0f);

        // Side/Mid ratio: linear below 100, 0 -> 6 dB above (0.3 = 6 dB / 20).
        const auto t = (w <= 100.0f) ? w * 0.01f : std::pow(10.0f, 0.3f * u);
        const auto thH = std::atan(t);
        g.gmH = kSqrt2 * std::cos(thH);
        g.gsH = kSqrt2 * std::sin(thH);

        // Low band: the angle is pulled towards 0 (mono) as the width rises above 100.
        const auto thL = thH * (1.0f - u);
        g.gmL = kSqrt2 * std::cos(thL);
        g.gsL = kSqrt2 * std::sin(thL);
        return g;
    }

    // xb: the weight of the crossover path against the broadband one. It fades in time (20 ms) when the width leaves
    // 100 upwards and back once it returns, never with the width itself: the two LR4 halves sum to an all-pass that is
    // -1 at the crossover, so a blend held at xb = 0.5 would notch the wet there.
    void processChunk(float* l, float* r, int n, float w, float xb) noexcept
    {
        if (! juce::exactlyEqual(w, cachedWidth)) // the trig only changes while the width moves
        {
            cachedGains = gainsFor(w);
            cachedWidth = w;
        }
        const auto gmH = cachedGains.gmH, gsH = cachedGains.gsH, gmL = cachedGains.gmL, gsL = cachedGains.gsL;

        if (xb <= 0.0f)
            releaseFilters();
        else
            filtersDirty = true;

        for (int i = 0; i < n; ++i)
        {
            const auto M = 0.5f * (l[i] + r[i]);
            const auto S = 0.5f * (l[i] - r[i]);

            auto oM = gmH * M;
            auto oS = gsH * S;

            if (xb > 0.0f)
            {
                const auto cM = gmH * mHi.process(M) + gmL * mLo.process(M);
                const auto cS = gsH * sHi.process(S) + gsL * sLo.process(S);
                oM += xb * (cM - oM);
                oS += xb * (cS - oS);
            }

            l[i] = oM + oS;
            r[i] = oM - oS;
        }
    }

    void releaseFilters() noexcept
    {
        if (filtersDirty)
            resetFilters();
    }

    void resetFilters() noexcept
    {
        mLo.reset();
        mHi.reset();
        sLo.reset();
        sHi.reset();
        filtersDirty = false;
    }

    LR4 mLo, mHi, sLo, sHi;
    juce::SmoothedValue<float, juce::ValueSmoothingTypes::Linear> width, engage;
    Gains cachedGains;
    float cachedWidth = -1.0f; // outside 0..200: nothing cached
    bool filtersDirty = false;
};
} // namespace rv::dsp
