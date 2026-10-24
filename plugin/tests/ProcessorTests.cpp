#include "PluginProcessor.h"

#include <catch2/catch_test_macros.hpp>

namespace
{
void setParam(Rev246Processor& p, const char* id, float real)
{
    auto* param = p.apvts.getParameter(id);
    param->setValueNotifyingHost(param->convertTo0to1(real));
}

juce::AudioBuffer<float> ramp(int n)
{
    juce::AudioBuffer<float> b(2, n);
    for (int i = 0; i < n; ++i)
    {
        b.setSample(0, i, 0.001f * (float)i);
        b.setSample(1, i, -0.0005f * (float)i);
    }
    return b;
}
} // namespace

TEST_CASE("power off passes the input through")
{
    Rev246Processor p;
    p.prepareToPlay(48000.0, 512);
    setParam(p, rv::params::id::power, 0.0f);
    juce::MidiBuffer midi;
    // Let the 20 ms crossfade finish, then compare.
    for (int blocks = 0; blocks < 4; ++blocks)
    {
        auto buffer = ramp(512);
        const auto original = ramp(512);
        p.processBlock(buffer, midi);
        if (blocks == 3)
            for (int ch = 0; ch < 2; ++ch)
                for (int i = 0; i < 512; ++i)
                    CHECK(buffer.getSample(ch, i) == original.getSample(ch, i));
    }
}

TEST_CASE("the VU meter reads the input after the INPUT fader")
{
    Rev246Processor p;
    p.prepareToPlay(48000.0, 256);
    juce::MidiBuffer midi;

    auto run = [&](float level, float seconds) {
        for (int blocks = 0; blocks < (int)(seconds * 48000.0f / 256.0f); ++blocks)
        {
            juce::AudioBuffer<float> b(2, 256);
            for (int ch = 0; ch < 2; ++ch)
                for (int i = 0; i < 256; ++i)
                    b.setSample(ch, i, level);
            p.processBlock(b, midi);
        }
    };

    run(0.5f, 0.6f); // let the VU ballistic settle
    CHECK(p.getInputVu() > 0.45f);
    CHECK(p.getInputVu() < 0.55f);

    setParam(p, rv::params::id::inputGain, rv::params::gainPositionFromDb(-12.0f));
    run(0.5f, 0.6f);
    const auto scaled = p.getInputVu();
    CHECK(scaled > 0.10f);
    CHECK(scaled < 0.15f); // -12 dB ≈ 0.126 of the 0.5 mono sum
}

TEST_CASE("the clip flag catches a full-scale input")
{
    Rev246Processor p;
    p.prepareToPlay(48000.0, 256);
    juce::MidiBuffer midi;
    juce::AudioBuffer<float> b(2, 256);
    for (int ch = 0; ch < 2; ++ch)
        for (int i = 0; i < 256; ++i)
            b.setSample(ch, i, 1.0f);
    p.processBlock(b, midi);
    CHECK(p.takeInputClip());
    CHECK_FALSE(p.takeInputClip()); // reset after the read
}

TEST_CASE("the pan keeps the centre at unity")
{
    Rev246Processor p;
    p.prepareToPlay(48000.0, 256);
    setParam(p, rv::params::id::mix, 0.0f);
    juce::MidiBuffer midi;
    for (int i = 0; i < 6; ++i)
    {
        juce::AudioBuffer<float> b(2, 256);
        for (int ch = 0; ch < 2; ++ch)
            for (int k = 0; k < 256; ++k)
                b.setSample(ch, k, 0.25f);
        p.processBlock(b, midi);
        if (i == 5)
        {
            CHECK(std::abs(b.getSample(0, 255) - 0.25f) < 1e-3f);
            CHECK(std::abs(b.getSample(1, 255) - 0.25f) < 1e-3f);
        }
    }
}
