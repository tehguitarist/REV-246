#include "PluginProcessor.h"
#include "dsp/CaptureEq.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;

namespace
{
constexpr double deviceRate = 19.2e6 / 728.0;

// 20 log10 |bq_resp(CLA, f)| from tools/unit_transfer.py at the device rate.
struct Point
{
    double hz, db;
};
constexpr Point reference[] = {{50, 2.7839},  {60, 2.1453},   {80, 1.0834},   {100, 0.3708},
                               {160, -0.6366}, {200, -0.4105}, {276, 0.4065},  {400, -0.0046},
                               {700, -0.0120}, {1000, -0.0062}, {2000, -0.0016}};

double db(std::complex<double> h)
{
    return 20.0 * std::log10(std::abs(h));
}
} // namespace

TEST_CASE("the capture EQ's coefficients match bq_resp at every host rate")
{
    for (double rate : {22050.0, deviceRate, 44100.0, 48000.0, 88200.0, 96000.0, 192000.0})
    {
        rv::dsp::CaptureEq eq;
        eq.prepare(rate);
        for (const auto& p : reference)
        {
            INFO("rate " << rate << " Hz, " << p.hz << " Hz");
            CHECK_THAT(db(eq.response(p.hz, rate)), WithinAbs(p.db, 0.005)); // 0.05 dB required; the table has 4 dp
        }
    }
}

TEST_CASE("the capture EQ's processing matches bq_resp")
{
    // The response measured from the processed impulse response, so the filter code is checked, not just its design.
    for (double rate : {48000.0, 96000.0})
    {
        rv::dsp::CaptureEq eq;
        eq.prepare(rate);
        std::vector<float> ir((size_t)rate); // 1 s: every section has decayed below -150 dB
        for (size_t i = 0; i < ir.size(); ++i)
            ir[i] = eq.process(0, i == 0 ? 1.0f : 0.0f);
        for (const auto& p : reference)
        {
            std::complex<double> h;
            const auto w = -2.0 * juce::MathConstants<double>::pi * p.hz / rate;
            for (size_t i = 0; i < ir.size(); ++i)
                h += (double)ir[i] * std::polar(1.0, w * (double)i);
            INFO("rate " << rate << " Hz, " << p.hz << " Hz");
            CHECK_THAT(db(h), WithinAbs(p.db, 0.01));
        }
    }
}

TEST_CASE("the capture EQ is off by default and leaves the output alone")
{
    Rev246Processor p;
    CHECK(p.apvts.getRawParameterValue(rv::params::id::captureEq)->load() < 0.5f);

    // Mix 100 %: the output is the wet path. The placeholder engine passes the mono input, so off the output must be
    // exactly mono x INPUT x OUTPUT (the default fader positions are within 0.05 dB of unity, not exactly 1); on must
    // change it.
    auto run = [](bool eqOn) {
        Rev246Processor proc;
        proc.prepareToPlay(48000.0, 256);
        auto* mix = proc.apvts.getParameter(rv::params::id::mix);
        mix->setValueNotifyingHost(1.0f);
        proc.apvts.getParameter(rv::params::id::captureEq)->setValueNotifyingHost(eqOn ? 1.0f : 0.0f);
        juce::MidiBuffer midi;
        float maxDiff = 0.0f;
        for (int block = 0; block < 8; ++block)
        {
            juce::AudioBuffer<float> b(2, 256);
            for (int ch = 0; ch < 2; ++ch)
                for (int i = 0; i < 256; ++i)
                    b.setSample(ch, i, 0.3f * std::sin(0.02f * (float)(block * 256 + i))); // ~150 Hz
            const auto in = b;
            proc.processBlock(b, midi);
            const auto g = rv::params::gainLinear(rv::params::defaultGainPosition);
            if (block >= 4) // after the 20 ms ramps
                for (int ch = 0; ch < 2; ++ch)
                    for (int i = 0; i < 256; ++i)
                    {
                        const auto expected = 0.5f * (in.getSample(0, i) + in.getSample(1, i)) * g * g;
                        maxDiff = std::max(maxDiff, std::abs(b.getSample(ch, i) - expected));
                    }
        }
        return maxDiff;
    };
    CHECK_FALSE(run(false) > 0.0f); // exactly the input
    CHECK(run(true) > 1e-3f);
}

TEST_CASE("presets leave the capture EQ alone")
{
    Rev246Processor p;
    auto* eq = p.apvts.getParameter(rv::params::id::captureEq);
    CHECK(p.presets.load("Init S"));
    eq->setValueNotifyingHost(1.0f);
    CHECK_FALSE(p.presets.isModified()); // a preference, not a setting
    CHECK(p.presets.load("Init M"));
    CHECK(eq->getValue() >= 0.5f);
}
