#include "PluginProcessor.h"
#include "TestEngine.h"
#include "dsp/AdcCalibration.h"
#include "ui/Design.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <complex>
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
                    CHECK(std::equal_to<float>{}(buffer.getSample(ch, i), original.getSample(ch, (i - latency + 512) % 512)));
    }
}

TEST_CASE("the factory defaults are period-accurate, fully wet and linked")
{
    Rev246Processor p;
    auto value = [&](const char* id) { return p.apvts.getRawParameterValue(id)->load(); };
    CHECK(value(rv::params::id::power) >= 0.5f);
    CHECK(std::equal_to<float>{}(value(rv::params::id::mix), 1.0f));
    CHECK(value(rv::params::id::link) >= 0.5f);
    CHECK(value(rv::params::id::transformer) >= 0.5f);
    CHECK(value(rv::params::id::sixteenBit) >= 0.5f);
    CHECK(std::equal_to<float>{}(value(rv::params::id::stereoWidth), 100.0f));
    CHECK(std::abs(rv::params::gainDbFromPosition(value(rv::params::id::inputGain))) < 0.01f);
    CHECK(std::abs(rv::params::gainDbFromPosition(value(rv::params::id::outputGain))) < 0.01f);
    CHECK(rv::params::joyDefaultStep(rv::params::Joy::lowpass) == rv::params::joyMaxStep(rv::params::Joy::lowpass));
    CHECK(rv::params::joyValue(rv::params::Program::s, rv::params::Joy::lowpass,
                               rv::params::joyDefaultStep(rv::params::Joy::lowpass))
              .text == "12");
    CHECK(std::equal_to<double>{}(rv::params::joyNumber(rv::params::Program::s, rv::params::Joy::rt,
                                                        rv::params::joyDefaultStep(rv::params::Joy::rt)),
                                  1.0)); // exactly 1.0 s
}

TEST_CASE("the latency is the engine's, the same in every program and mode, and dry and wet stay aligned")
{
    for (const auto rate : {44100.0, 96000.0})
    {
        rv::dsp::ReverbEngine engine;
        engine.prepare(rate, 256);
        int latency = -1;
        for (int program = 0; program < rv::params::numPrograms; ++program)
            for (const bool sixteen : {false, true})
            {
                Rev246Processor p;
                setParam(p, rv::params::id::program, (float)program);
                setParam(p, rv::params::id::sixteenBit, sixteen ? 1.0f : 0.0f);
                p.prepareToPlay(rate, 256);
                INFO("rate " << rate << ", program " << program << ", 16-bit " << sixteen);
                CHECK(p.getLatencySamples() == engine.latency());
                if (latency < 0)
                    latency = p.getLatencySamples();
                CHECK(p.getLatencySamples() == latency);
            }

        // The dry path (mix 0) is the input delayed by the latency; the wet's first arrival (S, DELAY 0, AMPL 100: the
        // cluster's tap 0, a unit tap) lands there too, plus the modelled analogue chain's own peak delay (a few
        // samples, not part of the latency: the unit has it as well).
        Rev246Processor p;
        p.prepareToPlay(rate, 256);
        setParam(p, rv::params::id::transformer, 0.0f);
        setParam(p, rv::params::id::sixteenBit, 0.0f);
        setParam(p, rv::params::id::delay, 0.0f);
        setParam(p, rv::params::id::ampl, 10.0f);
        juce::MidiBuffer midi;
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
                    out.push_back(b.getSample(1, i));
            }
            // the first arrival only: S's next R tap is 143 device samples (5.4 ms) later
            const auto from = out.begin() + 10 * 256, to = from + latency + (int)(0.003 * rate);
            const auto peak = std::max_element(from, to, [](float a, float b) { return std::abs(a) < std::abs(b); });
            const auto at = (int)(peak - from);
            INFO("rate " << rate << ", mix " << mix << ", peak at " << at << ", latency " << latency);
            if (mix == 0.0f)
                CHECK(at == latency);
            else
            {
                CHECK(at >= latency);
                CHECK(at <= latency + (int)std::ceil(0.00025 * rate)); // the chain's own peak delay: about 0.17 ms
            }
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

constexpr int routingBlock = 256, routingBlocks = 24;

// A short noise burst per channel (different on each), then silence.
std::vector<float> burst(unsigned seed)
{
    std::vector<float> x((size_t)(routingBlock * routingBlocks), 0.0f);
    juce::Random random((juce::int64)seed);
    for (size_t i = 0; i < 1024; ++i)
        x[i] = 0.2f * (random.nextFloat() - 0.5f);
    return x;
}

// The processor at its defaults (S) with the wet path isolated: mix 100 %, transformer off, float; returns each
// output channel.
std::vector<std::vector<float>> runProcessor(Rev246Processor& p, int channels, const std::vector<float>& l,
                                             const std::vector<float>& r)
{
    setParam(p, rv::params::id::mix, 1.0f);
    setParam(p, rv::params::id::transformer, 0.0f);
    setParam(p, rv::params::id::sixteenBit, 0.0f);
    p.prepareToPlay(48000.0, routingBlock);
    std::vector<std::vector<float>> out((size_t)p.getTotalNumOutputChannels());
    juce::MidiBuffer midi;
    for (int block = 0; block < routingBlocks; ++block)
    {
        juce::AudioBuffer<float> b(std::max(channels, p.getTotalNumOutputChannels()), routingBlock);
        b.clear();
        for (int i = 0; i < routingBlock; ++i)
        {
            const auto k = (size_t)(block * routingBlock + i);
            b.setSample(0, i, l[k]);
            if (channels > 1)
                b.setSample(1, i, r[k]);
        }
        p.processBlock(b, midi);
        for (size_t ch = 0; ch < out.size(); ++ch)
            for (int i = 0; i < routingBlock; ++i)
                out[ch].push_back(b.getSample((int)ch, i));
    }
    return out;
}

std::array<std::vector<float>, 2> runEngine(const std::vector<float>& in)
{
    return rv::test::runEngine(in, 48000.0, routingBlock);
}

// `got` matches `want` to within the INPUT / OUTPUT faders' 0 dB (a few 1e-6 off unity): relative to the peak.
void checkSame(const std::vector<float>& got, const std::vector<float>& want)
{
    REQUIRE(got.size() == want.size());
    float peak = 0.0f, err = 0.0f;
    for (size_t i = 0; i < got.size(); ++i)
    {
        peak = std::max(peak, std::abs(want[i]));
        err = std::max(err, std::abs(got[i] - want[i]));
    }
    REQUIRE(peak > 1e-4f);
    CHECK(err <= 1e-4f * peak);
}

std::vector<float> scaledSum(const std::vector<float>& a, const std::vector<float>& b, float g)
{
    std::vector<float> y(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        y[i] = g * (a[i] + b[i]);
    return y;
}
} // namespace

TEST_CASE("mono in / mono out folds the two read-outs down at equal power")
{
    Rev246Processor p;
    setLayout(p, 1, 1);
    const auto x = burst(1);
    const auto out = runProcessor(p, 1, x, x);
    const auto e = runEngine(x);
    checkSame(out[0], scaledSum(e[0], e[1], 1.0f / std::sqrt(2.0f)));
}

TEST_CASE("mono in / stereo out is the unit's own path")
{
    Rev246Processor p;
    setLayout(p, 1, 2);
    const auto x = burst(2);
    const auto out = runProcessor(p, 1, x, x);
    const auto e = runEngine(x);
    checkSame(out[0], e[0]);
    checkSame(out[1], e[1]);
}

TEST_CASE("stereo in / stereo out defaults to the original mono sum with -6 dB summing gain")
{
    Rev246Processor p;
    CHECK(p.apvts.getRawParameterValue(rv::params::id::trueStereo)->load() < 0.5f);
    setLayout(p, 2, 2);
    CHECK(p.isStereoToStereo());
    const auto l = burst(3), r = burst(4);
    const auto out = runProcessor(p, 2, l, r);
    // 0.5 * (L + R): summing two correlated channels does not add 6 dB.
    const auto e = runEngine(scaledSum(l, r, 0.5f));
    checkSame(out[0], e[0]);
    checkSame(out[1], e[1]);
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
    const auto out = runProcessor(p, 2, burst(5), burst(6));
    float diff = 0.0f, peak = 0.0f;
    for (size_t i = 0; i < out[0].size(); ++i)
    {
        diff = std::max(diff, std::abs(out[0][i] - out[1][i]));
        peak = std::max(peak, std::abs(out[0][i]));
    }
    REQUIRE(peak > 1e-4f);
    CHECK(diff < 1e-5f * peak);
}

TEST_CASE("true stereo runs one reverb per input channel")
{
    Rev246Processor p;
    setLayout(p, 2, 2);
    setParam(p, rv::params::id::trueStereo, 1.0f);
    const auto l = burst(7), r = burst(8);
    const auto out = runProcessor(p, 2, l, r);
    checkSame(out[0], runEngine(l)[0]); // engine 0's left
    checkSame(out[1], runEngine(r)[1]); // engine 1's right
}

namespace
{
constexpr double twoPi = juce::MathConstants<double>::twoPi;

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

double db(double x)
{
    return 20.0 * std::log10(x);
}
} // namespace

TEST_CASE("the meter calibration places the A/D full scale 6 dB over 0 VU, and OVL follows it")
{
    CHECK(std::abs(rv::dsp::adcFullScaleForZeroVu(-18.0f) - std::sqrt(2.0f) * std::pow(10.0f, -12.0f / 20.0f)) < 1e-6f);

    Rev246Processor p;
    p.prepareToPlay(48000.0, 480);
    setParam(p, rv::params::id::mix, 1.0f);
    setParam(p, rv::params::id::transformer, 0.0f);
    juce::MidiBuffer midi;
    // A 1 kHz sine with a mono-sum peak `amp` through S (linear and time-invariant, so the wet's 1 kHz amplitude is
    // proportional to the A/D's fundamental); returns that amplitude.
    auto play = [&](float amp, bool sixteen, bool& clipped) {
        setParam(p, rv::params::id::sixteenBit, sixteen ? 1.0f : 0.0f);
        std::vector<float> out;
        clipped = false;
        for (int block = 0; block < 200; ++block)
        {
            juce::AudioBuffer<float> b(2, 480);
            for (int i = 0; i < 480; ++i)
            {
                const auto v = amp * (float)std::sin(twoPi * 1000.0 * (block * 480 + i) / 48000.0);
                b.setSample(0, i, v);
                b.setSample(1, i, v);
            }
            p.processBlock(b, midi);
            if (block >= 20) // past any previous run's tail
                clipped |= p.takeInputClip();
            else
                p.takeInputClip();
            for (int i = 0; i < 480; ++i)
                out.push_back(b.getSample(1, i));
        }
        return amplitude(out, 48000.0, 1000.0);
    };
    bool clipped = false;
    const auto fs = (double)rv::dsp::adcFullScaleForZeroVu(-18.0f); // 0.355 peak at the default CAL
    // 3 dB over full scale: 16-BIT clips (the fundamental drops as the clip's), FLOAT does not, and OVL lights in both.
    const auto floatHot = play(0.5f, false, clipped);
    CHECK(clipped);
    const auto clippedHot = play(0.5f, true, clipped);
    CHECK(clipped);
    CHECK(std::abs(db(clippedHot / floatHot) - db(clippedHarmonic(0.5, fs, 1) / 0.5)) < 0.3);

    // CAL -12: full scale 0.71, no clipping and no OVL.
    p.setMeterCal(5);
    const auto clean = play(0.5f, true, clipped);
    CHECK_FALSE(clipped);
    CHECK(std::abs(db(clean / floatHot)) < 0.1);

    // The calibration is kept with the session.
    p.setMeterCal(4);
    juce::MemoryBlock state;
    p.getStateInformation(state);
    Rev246Processor q;
    q.setStateInformation(state.getData(), (int)state.getSize());
    CHECK(q.getMeterCal() == 4);
}
