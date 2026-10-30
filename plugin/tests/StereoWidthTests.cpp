#include "dsp/StereoWidth.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

namespace
{
float rms(const std::vector<float>& v, size_t from)
{
    double sum = 0.0;
    for (size_t i = from; i < v.size(); ++i)
        sum += (double)v[i] * v[i];
    return (float)std::sqrt(sum / (double)(v.size() - from));
}

void fillSine(std::vector<float>& l, std::vector<float>& r, double hz, double rate, float amplitude)
{
    for (size_t i = 0; i < l.size(); ++i)
    {
        const auto x = amplitude * (float)std::sin(2.0 * juce::MathConstants<double>::pi * hz * (double)i / rate);
        l[i] = x;
        r[i] = -x;
    }
}
} // namespace

TEST_CASE("stereo width 100 is a bit-exact bypass")
{
    rv::dsp::StereoWidener w;
    w.prepare(48000.0);
    w.setWidth(100.0f);

    std::vector<float> l(256), r(256);
    juce::Random random(1234);
    for (size_t i = 0; i < l.size(); ++i)
    {
        l[i] = random.nextFloat() * 2.0f - 1.0f;
        r[i] = random.nextFloat() * 2.0f - 1.0f;
    }
    const auto originalL = l, originalR = r;
    w.process(l.data(), r.data(), (int)l.size());
    CHECK(l == originalL);
    CHECK(r == originalR);
}

TEST_CASE("stereo width 0 sums to mono at equal power")
{
    rv::dsp::StereoWidener w;
    w.prepare(48000.0);
    w.setWidth(0.0f);

    std::vector<float> l(8192, 0.25f), r(8192, -0.25f);
    w.process(l.data(), r.data(), (int)l.size());
    CHECK(std::abs(l.back()) < 1e-6f);
    CHECK(std::abs(r.back()) < 1e-6f);

    std::fill(l.begin(), l.end(), 0.25f);
    std::fill(r.begin(), r.end(), 0.25f);
    w.process(l.data(), r.data(), (int)l.size());
    CHECK(std::abs(l.back() - 0.35355339f) < 1e-5f);
    CHECK(std::abs(r.back() - l.back()) < 1e-6f);
}

TEST_CASE("stereo width 200 boosts a high-frequency side signal")
{
    rv::dsp::StereoWidener w;
    w.prepare(48000.0);
    w.setWidth(200.0f);
    std::vector<float> l(16384), r(16384);
    fillSine(l, r, 1000.0, 48000.0, 0.5f);
    const auto inputRms = rms(l, 8192);
    w.process(l.data(), r.data(), (int)l.size());
    const auto ratio = rms(l, 8192) / inputRms;
    CHECK(std::abs(ratio - 1.2649f) < 0.05f); // +2 dB side gain at 200
}

TEST_CASE("stereo width 200 keeps bass mono")
{
    rv::dsp::StereoWidener w;
    w.prepare(48000.0);
    w.setWidth(200.0f);
    std::vector<float> l(16384), r(16384);
    fillSine(l, r, 50.0, 48000.0, 0.5f);
    w.process(l.data(), r.data(), (int)l.size());
    CHECK(rms(l, 8192) < 0.05f); // the side bass is pulled to mono
    CHECK(rms(r, 8192) < 0.05f);
}

TEST_CASE("stereo width keeps the level across 0 to 200")
{
    juce::Random random(5678);
    std::vector<float> l(16384), r(16384);
    for (size_t i = 0; i < l.size(); ++i)
    {
        l[i] = random.nextFloat() * 2.0f - 1.0f;
        r[i] = random.nextFloat() * 2.0f - 1.0f;
    }
    const auto inputRms = std::sqrt((rms(l, 0) * rms(l, 0) + rms(r, 0) * rms(r, 0)) * 0.5f);

    for (const auto width : {0.0f, 50.0f, 100.0f, 150.0f, 200.0f})
    {
        rv::dsp::StereoWidener w;
        w.prepare(48000.0);
        w.setWidth(width);
        auto left = l, right = r;
        w.process(left.data(), right.data(), (int)left.size());
        const auto outputRms = std::sqrt((rms(left, 8192) * rms(left, 8192) + rms(right, 8192) * rms(right, 8192)) * 0.5f);
        CHECK(std::abs(20.0f * std::log10(outputRms / inputRms)) < 1.0f);
    }
}

TEST_CASE("stereo width keeps the level of pink noise across 0 to 200")
{
    juce::Random random(24680);
    std::vector<float> l(24000), r(24000);
    float bl0 = 0, bl1 = 0, bl2 = 0, br0 = 0, br1 = 0, br2 = 0;
    for (size_t i = 0; i < l.size(); ++i)
    {
        const auto nl = random.nextFloat() * 2.0f - 1.0f;
        const auto nr = random.nextFloat() * 2.0f - 1.0f;
        bl0 = 0.99765f * bl0 + nl * 0.0990460f;
        bl1 = 0.96300f * bl1 + nl * 0.2965164f;
        bl2 = 0.57000f * bl2 + nl * 1.0526913f;
        br0 = 0.99765f * br0 + nr * 0.0990460f;
        br1 = 0.96300f * br1 + nr * 0.2965164f;
        br2 = 0.57000f * br2 + nr * 1.0526913f;
        l[i] = 0.1f * (bl0 + bl1 + bl2 + nl * 0.1848f);
        r[i] = 0.1f * (br0 + br1 + br2 + nr * 0.1848f);
    }
    const auto inputRms = std::sqrt((rms(l, 0) * rms(l, 0) + rms(r, 0) * rms(r, 0)) * 0.5f);

    for (const auto width : {0.0f, 50.0f, 100.0f, 150.0f, 200.0f})
    {
        rv::dsp::StereoWidener w;
        w.prepare(48000.0);
        w.setWidth(width);
        auto left = l, right = r;
        w.process(left.data(), right.data(), (int)left.size());
        const auto outputRms = std::sqrt((rms(left, 12000) * rms(left, 12000) + rms(right, 12000) * rms(right, 12000)) * 0.5f);
        CHECK(std::abs(20.0f * std::log10(outputRms / inputRms)) < 1.0f);
    }
}

TEST_CASE("stereo width just above 100 keeps the crossover frequency")
{
    // The crossover path fades in over time, not with the width: the LR4 halves sum to an all-pass that is -1 at 150 Hz,
    // so a width-driven blend held half way (the old 102.5) notched the wet there.
    const double rate = 48000.0;
    for (const auto width : {101.0f, 102.5f, 104.0f})
    {
        rv::dsp::StereoWidener w;
        w.prepare(rate);
        w.setWidth(width);
        std::vector<float> l(48000), r(48000);
        for (size_t i = 0; i < l.size(); ++i) // a 150 Hz mid (L = R) signal
            l[i] = r[i] = 0.25f * (float)std::sin(2.0 * juce::MathConstants<double>::pi * 150.0 * (double)i / rate);
        const auto in = rms(l, 24000);
        w.process(l.data(), r.data(), (int)l.size());
        INFO("width " << width);
        CHECK(std::abs(20.0f * std::log10(rms(l, 24000) / in)) < 0.5f);
    }
}
