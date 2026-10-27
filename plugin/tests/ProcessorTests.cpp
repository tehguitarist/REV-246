#include "PluginProcessor.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

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

TEST_CASE("power off passes the input through, delayed by the reported latency")
{
    Rev246Processor p;
    p.prepareToPlay(48000.0, 512);
    const auto latency = p.getLatencySamples();
    REQUIRE(latency > 0);
    REQUIRE(latency < 512);
    setParam(p, rv::params::id::power, 0.0f);
    juce::MidiBuffer midi;
    // Let the 20 ms crossfade finish, then compare. Every block is the same ramp, so the delayed signal is the ramp
    // shifted by the latency, wrapping into the previous block's copy.
    for (int blocks = 0; blocks < 4; ++blocks)
    {
        auto buffer = ramp(512);
        const auto original = ramp(512);
        p.processBlock(buffer, midi);
        if (blocks == 3)
            for (int ch = 0; ch < 2; ++ch)
                for (int i = 0; i < 512; ++i)
                    CHECK(buffer.getSample(ch, i) == original.getSample(ch, (i - latency + 512) % 512));
    }
}

TEST_CASE("the factory defaults are period-accurate, fully wet and linked")
{
    Rev246Processor p;
    auto value = [&](const char* id) { return p.apvts.getRawParameterValue(id)->load(); };
    CHECK(value(rv::params::id::power) >= 0.5f);
    CHECK(value(rv::params::id::mix) == 1.0f);
    CHECK(value(rv::params::id::link) >= 0.5f);
    CHECK(value(rv::params::id::transformer) >= 0.5f);
    CHECK(value(rv::params::id::sixteenBit) >= 0.5f);
    CHECK(value(rv::params::id::stereoWidth) == 100.0f);
    CHECK(std::abs(rv::params::gainDbFromPosition(value(rv::params::id::inputGain))) < 0.01f);
    CHECK(std::abs(rv::params::gainDbFromPosition(value(rv::params::id::outputGain))) < 0.01f);
    CHECK(rv::params::joyDefaultStep(rv::params::Joy::lowpass) == rv::params::joyMaxStep(rv::params::Joy::lowpass));
    CHECK(rv::params::joyValue(rv::params::Program::s, rv::params::Joy::lowpass,
                               rv::params::joyDefaultStep(rv::params::Joy::lowpass))
              .text == "12");
}

TEST_CASE("the latency is the same in float and 16-bit, and dry and wet stay aligned")
{
    for (const auto rate : {44100.0, 96000.0})
        for (const bool sixteen : {false, true})
        {
            Rev246Processor p;
            p.prepareToPlay(rate, 256);
            setParam(p, rv::params::id::transformer, 0.0f); // isolate the latency from the accuracy stages
            setParam(p, rv::params::id::sixteenBit, sixteen ? 1.0f : 0.0f);
            const auto latency = p.getLatencySamples();
            rv::dsp::AdcClipper adc;
            adc.prepare(rate);
            CHECK(latency == adc.latency()); // the A/D stage's, about 1.4 ms
            juce::MidiBuffer midi;
            // One small impulse, after the ramps have settled; the dry and the (pass-through) wet must coincide.
            for (const auto mix : {0.0f, 1.0f})
            {
                setParam(p, rv::params::id::mix, mix);
                std::vector<float> out;
                for (int block = 0; block < 20; ++block)
                {
                    juce::AudioBuffer<float> b(2, 256);
                    b.clear();
                    if (block == 10)
                        b.setSample(0, 0, 0.1f), b.setSample(1, 0, 0.1f);
                    p.processBlock(b, midi);
                    for (int i = 0; i < 256; ++i)
                        out.push_back(b.getSample(0, i));
                }
                const auto peak = std::max_element(out.begin(), out.end(), [](float a, float b) { return std::abs(a) < std::abs(b); });
                INFO("rate " << rate << ", 16-bit " << sixteen << ", mix " << mix);
                CHECK((int)(peak - out.begin()) == 10 * 256 + latency);
            }
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

TEST_CASE("the buses allow mono/mono, mono/stereo and stereo/stereo")
{
    Rev246Processor p;
    auto supports = [&](juce::AudioChannelSet in, juce::AudioChannelSet out) {
        juce::AudioProcessor::BusesLayout layout;
        layout.inputBuses.add(in);
        layout.outputBuses.add(out);
        return p.isBusesLayoutSupported(layout);
    };
    CHECK(supports(juce::AudioChannelSet::mono(), juce::AudioChannelSet::mono()));
    CHECK(supports(juce::AudioChannelSet::mono(), juce::AudioChannelSet::stereo()));
    CHECK(supports(juce::AudioChannelSet::stereo(), juce::AudioChannelSet::stereo()));
    CHECK_FALSE(supports(juce::AudioChannelSet::stereo(), juce::AudioChannelSet::mono()));
}

namespace
{
void setLayout(Rev246Processor& p, int in, int out)
{
    juce::AudioProcessor::BusesLayout layout;
    layout.inputBuses.add(in == 1 ? juce::AudioChannelSet::mono() : juce::AudioChannelSet::stereo());
    layout.outputBuses.add(out == 1 ? juce::AudioChannelSet::mono() : juce::AudioChannelSet::stereo());
    REQUIRE(p.setBusesLayout(layout));
}

// Runs a constant block long enough for the 20 ms ramps, then returns the last sample of each output.
std::vector<float> runConstant(Rev246Processor& p, int channels, float l, float r)
{
    setParam(p, rv::params::id::mix, 1.0f);
    setParam(p, rv::params::id::transformer, 0.0f);
    setParam(p, rv::params::id::sixteenBit, 0.0f);
    juce::MidiBuffer midi;
    for (int block = 0; block < 6; ++block)
    {
        juce::AudioBuffer<float> b(channels, 256);
        for (int i = 0; i < 256; ++i)
        {
            b.setSample(0, i, l);
            if (channels > 1)
                b.setSample(1, i, r);
        }
        p.processBlock(b, midi);
        if (block == 5)
            return {b.getSample(0, 255), channels > 1 ? b.getSample(1, 255) : 0.0f};
    }
    return {};
}
} // namespace

TEST_CASE("mono in / mono out uses the reverb's left side")
{
    Rev246Processor p;
    setLayout(p, 1, 1);
    p.prepareToPlay(48000.0, 256);
    const auto out = runConstant(p, 1, 0.25f, 0.0f);
    CHECK(std::abs(out[0] - 0.25f) < 1e-3f);
}

TEST_CASE("mono in / stereo out feeds both outputs")
{
    Rev246Processor p;
    setLayout(p, 1, 2);
    p.prepareToPlay(48000.0, 256);
    const auto out = runConstant(p, 2, 0.25f, -0.25f); // the second input is ignored in a mono layout
    CHECK(std::abs(out[0] - 0.25f) < 1e-3f);
    CHECK(std::abs(out[1] - 0.25f) < 1e-3f);
}

TEST_CASE("stereo in / stereo out defaults to the original mono sum with -6 dB summing gain")
{
    Rev246Processor p;
    CHECK(p.apvts.getRawParameterValue(rv::params::id::trueStereo)->load() < 0.5f);
    setLayout(p, 2, 2);
    p.prepareToPlay(48000.0, 256);
    CHECK(p.isStereoToStereo());
    // 0.5 * (0.2 + -0.1) = 0.05, not 0.1: summing two correlated channels does not add 6 dB.
    const auto out = runConstant(p, 2, 0.2f, -0.1f);
    CHECK(std::abs(out[0] - 0.05f) < 1e-3f);
    CHECK(std::abs(out[1] - 0.05f) < 1e-3f);
}

TEST_CASE("the stereo width defaults to 100 and can make the wet mono")
{
    Rev246Processor p;
    auto* width = p.apvts.getParameter(rv::params::id::stereoWidth);
    CHECK(std::abs(width->getDefaultValue() - 0.5f) < 1e-6f); // 100 of 0..200
    CHECK(std::abs(width->convertFrom0to1(0.0f) - 0.0f) < 1e-6f);
    CHECK(std::abs(width->convertFrom0to1(1.0f) - 200.0f) < 1e-6f);

    setLayout(p, 2, 2);
    setParam(p, rv::params::id::trueStereo, 1.0f);
    setParam(p, rv::params::id::stereoWidth, 0.0f);
    p.prepareToPlay(48000.0, 256);
    const auto out = runConstant(p, 2, 0.2f, -0.1f);
    CHECK(std::abs(out[0] - out[1]) < 1e-3f);
}

TEST_CASE("true stereo runs one reverb per input channel")
{
    Rev246Processor p;
    setLayout(p, 2, 2);
    setParam(p, rv::params::id::trueStereo, 1.0f);
    p.prepareToPlay(48000.0, 256);
    const auto out = runConstant(p, 2, 0.2f, -0.1f);
    CHECK(std::abs(out[0] - 0.2f) < 1e-3f);
    CHECK(std::abs(out[1] + 0.1f) < 1e-3f);
}
