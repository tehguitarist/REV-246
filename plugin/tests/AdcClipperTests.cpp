#include "PluginProcessor.h"
#include "dsp/AdcClipper.h"
#include "ui/Design.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <vector>

namespace
{
constexpr double twoPi = juce::MathConstants<double>::twoPi;
constexpr double deviceRate = rv::dsp::AdcClipper::deviceRate;

std::vector<float> sine(double rate, double hz, double amp, int n)
{
    std::vector<float> x((size_t)n);
    for (int i = 0; i < n; ++i)
        x[(size_t)i] = (float)(amp * std::sin(twoPi * hz * i / rate));
    return x;
}

std::vector<float> run(rv::dsp::AdcClipper& c, const std::vector<float>& x)
{
    std::vector<float> y(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        y[i] = c.process(x[i]);
    return y;
}

// Amplitude of the component at `hz` over the second half of `y` (Hann window).
double amplitude(const std::vector<float>& y, double rate, double hz)
{
    const auto start = y.size() / 2, n = y.size() - start;
    std::complex<double> sum;
    double wsum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const auto w = 0.5 - 0.5 * std::cos(twoPi * (double)i / (double)n);
        sum += w * (double)y[start + i] * std::polar(1.0, -twoPi * hz * (double)(start + i) / rate);
        wsum += w;
    }
    return 2.0 * std::abs(sum) / wsum;
}

// Harmonic k of a sine of amplitude `amp` hard-clipped at +-fs (continuous time).
double clippedHarmonic(double amp, double fs, int k)
{
    constexpr int n = 200000;
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const auto ph = twoPi * (i + 0.5) / n;
        sum += std::clamp(amp * std::sin(ph), -fs, fs) * std::sin(k * ph);
    }
    return 2.0 * std::abs(sum) / n;
}

double db(double ratio)
{
    return 20.0 * std::log10(ratio);
}
} // namespace

TEST_CASE("below full scale the A/D stage is an exact delay")
{
    const auto rate = GENERATE(44100.0, 48000.0, 96000.0);
    rv::dsp::AdcClipper c;
    c.prepare(rate);
    c.setFullScale(0.5f);
    auto x = sine(rate, 1000.0, 0.3, (int)rate / 4);
    const auto hf = sine(rate, 10000.0, 0.05, (int)x.size()); // +5 dB through the input filter: still under
    for (size_t i = 0; i < x.size(); ++i)
        x[i] += hf[i];
    const auto y = run(c, x);
    REQUIRE(c.latency() > 0);
    CHECK(c.latency() < (int)(1.6e-3 * rate));
    for (size_t i = (size_t)c.latency(); i < x.size(); ++i)
        REQUIRE(y[i] == x[i - (size_t)c.latency()]);
    CHECK_FALSE(c.takeOverload());
}

TEST_CASE("a sine 3 dB over full scale gives the hard-clipped harmonics, back through the input filter")
{
    const auto rate = GENERATE(44100.0, 48000.0, 96000.0);
    rv::dsp::AdcClipper c;
    c.prepare(rate);
    const double fs = 0.3, amp = fs * std::pow(10.0, 3.0 / 20.0);
    c.setFullScale((float)fs);
    const auto y = run(c, sine(rate, 1000.0, amp, (int)rate));
    CHECK(c.takeOverload());
    for (const int k : {1, 3, 5, 7})
    {
        const auto expected = clippedHarmonic(amp, fs, k) / rv::dsp::AdcClipper::inputFilterGain(1000.0 * k);
        INFO("rate " << rate << ", harmonic " << k);
        CHECK(std::abs(db(amplitude(y, rate, 1000.0 * k) / expected)) < 0.2);
    }
    CHECK(amplitude(y, rate, 2000.0) < 1e-4 * amp); // symmetric clip: no even harmonics
}

TEST_CASE("the input filter's pre-emphasis makes high frequencies clip under full scale")
{
    rv::dsp::AdcClipper c;
    c.prepare(48000.0);
    c.setFullScale(0.3f);
    // 9 kHz at -0.9 dB re full scale is +3.8 dB at the A/D (+4.7 dB pre-emphasis): clips.
    const auto y = run(c, sine(48000.0, 9000.0, 0.3 * std::pow(10.0, -0.9 / 20.0), 24000));
    CHECK(c.takeOverload());
    CHECK(amplitude(y, 48000.0, 9000.0) < 0.95 * 0.3 * std::pow(10.0, -0.9 / 20.0));
    // 1 kHz at the same level does not, once the onset has passed (switching a sine on abruptly rings the
    // filter's 12 kHz edge, which can clip: a real effect).
    c.reset();
    const auto low = sine(48000.0, 1000.0, 0.3 * std::pow(10.0, -0.9 / 20.0), 24000);
    run(c, std::vector<float>(low.begin(), low.begin() + 12000));
    c.takeOverload();
    for (size_t i = 12000; i < low.size(); ++i)
        c.process(low[i]);
    CHECK_FALSE(c.takeOverload());
}

TEST_CASE("the clip is sampled at the device rate, so its harmonics fold around 13.2 kHz")
{
    rv::dsp::AdcClipper c;
    c.prepare(48000.0);
    const double fs = 0.2, amp = 2.0 * fs / rv::dsp::AdcClipper::inputFilterGain(5000.0); // +6 dB at the A/D
    c.setFullScale((float)fs);
    const auto y = run(c, sine(48000.0, 5000.0, amp, 48000));
    const auto fundamental = amplitude(y, 48000.0, 5000.0);
    // The 3rd harmonic (15 kHz) lands at 26373.6 - 15000 = 11373.6 Hz, at its hard-clip level re the A/D, and is
    // 40 dB weaker at 15 kHz itself, where a clip at the host rate would put it (what is there is the leakage of
    // the truncated reinjection kernel, -38 dB of its energy).
    const auto alias = deviceRate - 15000.0;
    const auto expected = clippedHarmonic(2.0 * fs, fs, 3) / rv::dsp::AdcClipper::inputFilterGain(alias);
    CHECK(std::abs(db(amplitude(y, 48000.0, alias) / expected)) < 0.5);
    CHECK(db(amplitude(y, 48000.0, 15000.0) / amplitude(y, 48000.0, alias)) < -40.0);
    CHECK(db(amplitude(y, 48000.0, 15000.0) / fundamental) < -50.0);
}

TEST_CASE("the meter calibration places the A/D full scale 6 dB over 0 VU")
{
    CHECK(std::abs(rv::dsp::AdcClipper::fullScaleForZeroVu(-18.0f) - std::sqrt(2.0f) * std::pow(10.0f, -12.0f / 20.0f)) < 1e-6f);

    Rev246Processor p;
    p.prepareToPlay(48000.0, 480);
    auto set = [&](const char* id, float real) {
        auto* param = p.apvts.getParameter(id);
        param->setValueNotifyingHost(param->convertTo0to1(real));
    };
    set(rv::params::id::mix, 1.0f);
    set(rv::params::id::transformer, 0.0f);
    set(rv::params::id::sixteenBit, 1.0f);
    juce::MidiBuffer midi;
    // A 1 kHz sine whose mono sum peaks at 0.5: 3 dB over the default A/D full scale (CAL -18: 0.355 peak).
    auto play = [&](bool& clipped) {
        std::vector<float> out;
        clipped = false;
        for (int block = 0; block < 100; ++block)
        {
            juce::AudioBuffer<float> b(2, 480);
            for (int i = 0; i < 480; ++i)
            {
                const auto v = 0.5f * (float)std::sin(twoPi * 1000.0 * (block * 480 + i) / 48000.0);
                b.setSample(0, i, v);
                b.setSample(1, i, v);
            }
            p.processBlock(b, midi);
            clipped |= p.takeInputClip();
            for (int i = 0; i < 480; ++i)
                out.push_back(b.getSample(0, i));
        }
        return amplitude(out, 48000.0, 1000.0);
    };
    bool clipped = false;
    const auto fs = rv::dsp::AdcClipper::fullScaleForZeroVu(-18.0f);
    CHECK(std::abs(db(play(clipped) / clippedHarmonic(0.5, fs, 1))) < 0.2);
    CHECK(clipped);

    // CAL -12: full scale 0.71, no clipping.
    p.setMeterCal(5);
    CHECK(std::abs(db(play(clipped) / 0.5)) < 0.05);
    CHECK_FALSE(clipped);

    // Float mode still lights OVL where the A/D would clip, but leaves the signal alone.
    p.setMeterCal(rv::design::defaultMeterCal);
    set(rv::params::id::sixteenBit, 0.0f);
    CHECK(std::abs(db(play(clipped) / 0.5)) < 0.05);
    CHECK(clipped);

    // The calibration is kept with the session.
    p.setMeterCal(4);
    juce::MemoryBlock state;
    p.getStateInformation(state);
    Rev246Processor q;
    q.setStateInformation(state.getData(), (int)state.getSize());
    CHECK(q.getMeterCal() == 4);
}

TEST_CASE("at low host rates the A/D stage still passes a sub-clip sine at its level")
{
    // Zeros and poles of the input filter above Nyquist used to fold into the band and amplify a 1 kHz sine ~80x.
    const auto rate = GENERATE(22050.0, 24000.0, 32000.0);
    rv::dsp::AdcClipper c;
    c.prepare(rate);
    c.setFullScale(0.3f);
    const auto y = run(c, sine(rate, 1000.0, 0.1, (int)rate / 2));
    REQUIRE(amplitude(y, rate, 1000.0) == Catch::Approx(0.1).margin(0.005));
}
