#include "PluginProcessor.h"
#include "dsp/Transformer.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <complex>
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

// Harmonics 1..7 of a sine at `hz`, amplitude `amp`, through one transformer core (the second of two seconds).
std::array<double, 8> harmonics(double rate, double hz, float amp)
{
    rv::dsp::TransformerCore core;
    core.prepare(rate, 1.0);
    const auto n = (int)(2.0 * rate);
    std::array<std::complex<double>, 8> h{};
    for (int i = 0; i < n; ++i)
    {
        const auto t = (double)i / rate;
        const auto y = (double)core.process(amp * (float)std::sin(juce::MathConstants<double>::twoPi * hz * t));
        if (i >= n / 2)
            for (int m = 1; m < 8; ++m)
                h[(size_t)m] += y * std::polar(1.0, -juce::MathConstants<double>::twoPi * m * hz * t);
    }
    std::array<double, 8> mag{};
    for (int m = 1; m < 8; ++m)
        mag[(size_t)m] = 2.0 * std::abs(h[(size_t)m]) / (double)(n - n / 2);
    return mag;
}

double thdPercent(const std::array<double, 8>& h)
{
    double sum = 0.0;
    for (size_t m = 2; m < h.size(); ++m)
        sum += h[m] * h[m];
    return 100.0 * std::sqrt(sum) / h[1];
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

TEST_CASE("the transformer is off by default and leaves the output alone")
{
    Rev246Processor p;
    CHECK(p.apvts.getRawParameterValue(rv::params::id::transformer)->load() < 0.5f);

    // Mix 100 %: the output is the wet path. The placeholder engine passes the mono input, so off the output must be
    // exactly mono x INPUT x OUTPUT (the default fader positions are within 0.05 dB of unity, not exactly 1); on must
    // change it.
    auto run = [](bool eqOn) {
        Rev246Processor proc;
        proc.prepareToPlay(48000.0, 256);
        auto* mix = proc.apvts.getParameter(rv::params::id::mix);
        mix->setValueNotifyingHost(1.0f);
        proc.apvts.getParameter(rv::params::id::transformer)->setValueNotifyingHost(eqOn ? 1.0f : 0.0f);
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

TEST_CASE("presets leave the transformer alone")
{
    Rev246Processor p;
    auto* eq = p.apvts.getParameter(rv::params::id::transformer);
    CHECK(p.presets.load("Init S"));
    eq->setValueNotifyingHost(1.0f);
    CHECK_FALSE(p.presets.isModified()); // a preference, not a setting
    CHECK(p.presets.load("Init M"));
    CHECK(eq->getValue() >= 0.5f);
}

TEST_CASE("the transformer core's distortion follows its calibration")
{
    // 0 dBFS = +12 dBu. Calibrated to 0.08 % at 40 Hz (Transformer.h); 0.078 % with the 5th and 7th included.
    const auto at40 = harmonics(48000.0, 40.0, 1.0f);
    CHECK_THAT(thdPercent(at40), WithinAbs(0.078, 0.004));
    CHECK(at40[5] < 0.15 * at40[3]); // mostly 3rd: the 5th is about 20 dB under it
    for (size_t even : {2u, 4u, 6u})
        CHECK(at40[even] < 1e-6 * at40[1]); // odd only

    // Falls steeply with frequency and with level; the same at another rate.
    CHECK(thdPercent(harmonics(48000.0, 20.0, 1.0f)) > 4.0 * thdPercent(at40));
    CHECK(thdPercent(harmonics(48000.0, 1000.0, 1.0f)) < 0.001);
    CHECK(thdPercent(harmonics(48000.0, 40.0, 0.25f)) < 0.3 * thdPercent(at40));
    CHECK_THAT(thdPercent(harmonics(96000.0, 40.0, 1.0f)), WithinAbs(thdPercent(at40), 0.002));
}
