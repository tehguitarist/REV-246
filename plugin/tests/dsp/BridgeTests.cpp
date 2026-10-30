#include "dsp/engine/DeviceBridge.h"
#include "dsp/generated/bridge_spec.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

// DeviceBridge against the modelled analogue chain (tools/plugin_codegen/bridge.py: the MNA solver's input filter x
// output filter x S&H lag with the unit's coupling corners, unity at 1 kHz, times the output-side correction), through an
// identity program; and the corrected chain against its target, the unit part of the measured chain3
// (analysis/plugin_check/REPORT_chainfix.md).

namespace
{
using rv::dsp::DeviceBridge;
namespace gen = rv::dsp::gen::bridge;
using C = std::complex<double>;

constexpr std::array<double, 4> rates{44100.0, 48000.0, 96000.0, 192000.0};
constexpr double twoPi = 6.283185307179586;

void identity(const float* in, float* l, float* r, int n)
{
    std::copy(in, in + n, l);
    std::copy(in, in + n, r);
}

struct Out
{
    std::vector<float> l, r;
};

template <typename F>
Out render(DeviceBridge& b, const std::vector<float>& x, int block, F&& device)
{
    Out o{std::vector<float>(x.size()), std::vector<float>(x.size())};
    for (size_t i = 0; i < x.size(); i += (size_t)block)
    {
        const auto n = (int)std::min((size_t)block, x.size() - i);
        b.process(x.data() + i, o.l.data() + i, o.r.data() + i, n, device);
    }
    return o;
}

std::vector<float> cosine(double fs, double f, double amp, size_t n)
{
    std::vector<float> x(n);
    for (size_t m = 0; m < n; ++m)
        x[m] = (float)(amp * std::cos(twoPi * f * (double)m / fs));
    return x;
}

// Complex amplitude of y at f over [start, start + len) with a 4-term Blackman-Harris window; phase re time `origin`.
C toneAt(const std::vector<float>& y, double fs, double f, size_t start, size_t len, double origin)
{
    C acc(0.0, 0.0);
    double wsum = 0.0;
    for (size_t k = 0; k < len; ++k)
    {
        const auto x = twoPi * (double)k / (double)len;
        const auto w = 0.35875 - 0.48829 * std::cos(x) + 0.14128 * std::cos(2 * x) - 0.01168 * std::cos(3 * x);
        const auto m = (double)(start + k);
        acc += w * (double)y[start + k] * std::polar(1.0, -twoPi * f * (m - origin) / fs);
        wsum += w;
    }
    return 2.0 * acc / wsum;
}

// Bit-for-bit equality (no float == in the checks).
bool sameBits(const std::vector<float>& a, const std::vector<float>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0;
}

double db(double x)
{
    return 20.0 * std::log10(std::max(x, 1e-300));
}
} // namespace

TEST_CASE("Bridge: chain response through an identity program, every host rate", "[bridge]")
{
    REQUIRE(std::abs(gen::refShPole - DeviceBridge::shPole) < 1e-15);
    struct Band
    {
        double lo, hi, worstAbs = 0.0, worstRel = 0.0;
    };
    for (const auto fs : rates)
    {
        DeviceBridge b;
        b.prepare(fs, 512);
        const auto d = (size_t)b.latency();
        std::array<Band, 5> bands{{{20, 100}, {100, 1000}, {1000, 5000}, {5000, 10000}, {10000, 12800.1}}};
        double worstLr = 0.0, phase1k = 0.0, gain1k = 0.0;
        for (size_t i = 0; i < gen::refHz.size(); ++i)
        {
            const auto f = gen::refHz[i];
            const C ref(gen::refRe[i], gen::refIm[i]);
            const auto pre = (size_t)((f < 200.0 ? 0.6 : 0.1) * fs);
            const auto len = (size_t)((f < 200.0 ? 1.0 : 0.25) * fs);
            constexpr double amp = 0.5;
            b.reset();
            const auto o = render(b, cosine(fs, f, amp, pre + len + d), 512, identity);
            const auto hr = toneAt(o.r, fs, f, pre + d, len, (double)d) / amp;
            const auto hl = toneAt(o.l, fs, f, pre + d, len, (double)d) / amp;
            const auto absErr = std::abs(hr - ref);
            const auto relErr = std::abs(hr / ref - 1.0);
            for (auto& band : bands)
                if (f >= band.lo && f < band.hi)
                {
                    band.worstAbs = std::max(band.worstAbs, absErr);
                    band.worstRel = std::max(band.worstRel, relErr);
                }
            // L: its own correction, lLateDeviceSamples later.
            worstLr = std::max(worstLr, std::abs(hl - C(gen::refLRe[i], gen::refLIm[i])));
            if (std::abs(f - 1000.0) < 1e-9)
            {
                phase1k = std::arg(hr / ref);
                gain1k = std::abs(hr);
            }
            INFO("fs " << fs << " f " << f << " |err| " << db(absErr) << " dB, rel " << db(relErr) << " dB");
            CHECK(absErr < std::pow(10.0, -90.0 / 20.0));
        }
        std::printf("Bridge %6.0f Hz: latency %d host samples (%.3f ms)\n", fs, b.latency(), 1000.0 * b.latency() / fs);
        for (const auto& band : bands)
            std::printf("  %6.0f-%6.0f Hz: |error| %7.1f dB re unity, %7.1f dB re the chain\n", band.lo, band.hi,
                        db(band.worstAbs), db(band.worstRel));
        std::printf("  L against its reference: %.1f dB; 1 kHz gain %.8f, phase error %.2e rad\n", db(worstLr),
                    gain1k, phase1k);
        CHECK(worstLr < std::pow(10.0, -90.0 / 20.0));
        CHECK(std::abs(gain1k - 1.0) < 1e-5);
        CHECK(std::abs(phase1k) < 1e-5);
    }
}

TEST_CASE("Bridge: a tone above device Nyquist folds at the level the input filter allows", "[bridge]")
{
    for (const auto fs : rates)
    {
        DeviceBridge b;
        b.prepare(fs, 512);
        const auto d = (size_t)b.latency();
        for (size_t i = 0; i < gen::aliasHz.size(); ++i)
        {
            const auto f = gen::aliasHz[i];
            const auto folded = rv::dsp::deviceRate - f;
            const C expected(gen::aliasRe[i], gen::aliasIm[i]);
            const auto pre = (size_t)(0.1 * fs), len = (size_t)(0.25 * fs);
            b.reset();
            const auto o = render(b, cosine(fs, f, 1.0, pre + len + d), 512, identity);
            const auto got = toneAt(o.r, fs, folded, pre + d, len, (double)d);
            const auto rel = std::abs(got / expected - 1.0);
            std::printf("Bridge alias %6.0f Hz: %.0f Hz in -> %.1f Hz out at %.2f dB (input filter %.2f dB), "
                        "error %.4f dB / %.2e\n",
                        fs, f, folded, db(std::abs(got)), gen::aliasInputGainDb[i], db(std::abs(got)) - db(std::abs(expected)),
                        rel);
            // The host-rate interpretation (phi) is flat to 12.8 kHz and rolls off towards host Nyquist: 16 kHz at
            // 44.1 kHz is 1 % low (unit: an input above 13.2 kHz at a 44.1 kHz host is -60 dB and lower after folding).
            CHECK(rel < (f < 15000.0 ? 1e-3 : 0.02));
        }
    }
}

TEST_CASE("Bridge: L lags R by lLateDeviceSamples, with its own correction", "[bridge]")
{
    // From the ratio L / R on steady tones against the references' ratio (the response test checks each to -90 dB). The
    // two corrections differ by about 0.1 dB and well under a degree below 8 kHz, so the lateness reads directly at 1 kHz.
    REQUIRE(std::abs(DeviceBridge::lLateDeviceSamples - 0.5037) < 1e-12);
    for (const auto fs : {44100.0, 192000.0})
    {
        DeviceBridge b;
        b.prepare(fs, 512);
        const auto d = (size_t)b.latency();
        for (size_t i = 0; i < gen::refHz.size(); ++i)
        {
            const auto f = gen::refHz[i];
            if (std::abs(f - 1000.0) > 1e-9 && std::abs(f - 10000.0) > 1e-9)
                continue;
            b.reset();
            const auto pre = (size_t)(0.1 * fs), len = (size_t)(0.25 * fs);
            const auto o = render(b, cosine(fs, f, 0.5, pre + len + d), 512, identity);
            const auto ratio = toneAt(o.l, fs, f, pre + d, len, 0.0) / toneAt(o.r, fs, f, pre + d, len, 0.0);
            const auto expected = C(gen::refLRe[i], gen::refLIm[i]) / C(gen::refRe[i], gen::refIm[i]);
            const auto lateUs = -std::arg(ratio) / (twoPi * f) * 1e6;
            std::printf("Bridge %6.0f Hz, %5.0f Hz tone: L %.6f us after R (lLateDeviceSamples %.6f us), level %.4f dB, "
                        "error re the references %.1f dB\n",
                        fs, f, lateUs, 1e6 * DeviceBridge::lLateDeviceSamples / rv::dsp::deviceRate, db(std::abs(ratio)),
                        db(std::abs(ratio / expected - 1.0)));
            CHECK(std::abs(ratio / expected - 1.0) < 1e-4);
        }
    }
}

TEST_CASE("Bridge: corrected chain against the unit part of the measured chain", "[bridge]")
{
    // targetHz holds 6 points in each scoring band of plugin_check (20 Hz-12.3 kHz). Error energy over target energy per
    // band, per channel. Below 80 Hz the target keeps chain3's own shape, which the captures do not support (its kernel is
    // too short for the 1.8 Hz corner; with the analytic corner the plugin nulls the captures far better there, REPORT_
    // chainfix.md), so those two bands sit near -36 / -45 dB; from 80 Hz up the correction holds -50 dB.
    constexpr size_t perBand = 6;
    const auto numBands = gen::targetEdgesHz.size() - 1;
    REQUIRE(gen::targetHz.size() == perBand * numBands);
    for (const auto fs : {48000.0, 96000.0})
    {
        DeviceBridge b;
        b.prepare(fs, 512);
        const auto d = (size_t)b.latency();
        std::vector<double> errE[2], refE[2];
        for (auto c = 0; c < 2; ++c)
        {
            errE[c].assign(numBands, 0.0);
            refE[c].assign(numBands, 0.0);
        }
        for (size_t i = 0; i < gen::targetHz.size(); ++i)
        {
            const auto f = gen::targetHz[i];
            const auto pre = (size_t)((f < 200.0 ? 0.6 : 0.1) * fs);
            const auto len = (size_t)((f < 200.0 ? 1.0 : 0.25) * fs);
            b.reset();
            const auto o = render(b, cosine(fs, f, 0.5, pre + len + d), 512, identity);
            const C got[2]{toneAt(o.l, fs, f, pre + d, len, (double)d) / 0.5, toneAt(o.r, fs, f, pre + d, len, (double)d) / 0.5};
            const C tgt[2]{C(gen::targetLRe[i], gen::targetLIm[i]), C(gen::targetRRe[i], gen::targetRIm[i])};
            for (auto c = 0; c < 2; ++c)
            {
                errE[c][i / perBand] += std::norm(got[c] - tgt[c]);
                refE[c][i / perBand] += std::norm(tgt[c]);
            }
        }
        std::printf("Bridge %6.0f Hz against the unit part of chain3 (error re target, dB, L / R):\n", fs);
        for (size_t k = 0; k < numBands; ++k)
        {
            const double e[2]{10.0 * std::log10(errE[0][k] / refE[0][k]), 10.0 * std::log10(errE[1][k] / refE[1][k])};
            std::printf("  %6.0f-%6.0f Hz: %6.1f / %6.1f\n", gen::targetEdgesHz[k], gen::targetEdgesHz[k + 1], e[0], e[1]);
            const auto bar = gen::targetEdgesHz[k] < 39.0 ? -32.0 : gen::targetEdgesHz[k] < 79.0 ? -42.0 : -50.0;
            CHECK(e[0] < bar);
            CHECK(e[1] < bar);
        }
    }
}

TEST_CASE("Bridge: latency, causality of the reported delay", "[bridge]")
{
    for (const auto fs : rates)
    {
        DeviceBridge b;
        b.prepare(fs, 512);
        const auto d = b.latency();
        std::vector<float> x((size_t)(0.05 * fs), 0.0f);
        x[0] = 1.0f;
        const auto o = render(b, x, 512, identity);
        double before = 0.0, total = 0.0;
        for (size_t m = 0; m < x.size(); ++m)
        {
            const auto e = (double)o.r[m] * o.r[m];
            total += e;
            if ((int)m < d)
                before += e;
        }
        // The chain is causal; what precedes the reported delay is only the band-limiting sinc's pre-ringing.
        std::printf("Bridge %6.0f Hz: energy before the latency %.1f dB\n", fs, 10.0 * std::log10(before / total + 1e-300));
        CHECK(before < 1e-8 * total);
    }
}

TEST_CASE("Bridge: A/D clip hook, overload flag and float mode", "[bridge]")
{
    DeviceBridge b;
    b.prepare(48000.0, 512);
    b.setFullScale(0.5f);
    std::vector<float> seen;
    auto capture = [&](const float* in, float* l, float* r, int n) {
        seen.insert(seen.end(), in, in + n);
        identity(in, l, r, n);
    };
    const auto lo = -0.5f, hi = (float)(0.5 * 32767.0 / 32768.0);

    // 16-bit mode: a 1 kHz sine at 1.0 (twice full scale) clips at the device rate, at exactly the limits.
    b.setClipping(true);
    render(b, cosine(48000.0, 1000.0, 1.0, 9600), 512, capture);
    const auto [mn, mx] = std::minmax_element(seen.begin() + 1000, seen.end());
    CHECK(std::abs(*mn - lo) < 1e-9f);
    CHECK(std::abs(*mx - hi) < 1e-9f);
    CHECK(b.takeOverload());
    CHECK_FALSE(b.takeOverload()); // cleared by the read

    // Float mode: no clip, but the flag still says where the A/D would have clipped.
    b.reset();
    seen.clear();
    b.setClipping(false);
    render(b, cosine(48000.0, 1000.0, 1.0, 9600), 512, capture);
    CHECK(*std::max_element(seen.begin() + 1000, seen.end()) > 0.99f);
    CHECK(b.takeOverload());

    // Below full scale at 1 kHz: no overload. The same level at 10 kHz overloads: the input filter's +5 dB lift.
    for (const auto clip : {false, true})
    {
        b.setClipping(clip);
        b.reset();
        seen.clear();
        // (after 50 ms: the cosine's onset is a step, and the filter's step response overshoots)
        const auto x = cosine(48000.0, 1000.0, 0.45, 9600);
        const std::vector<float> onset(x.begin(), x.begin() + 2400), rest(x.begin() + 2400, x.end());
        render(b, onset, 512, capture);
        b.takeOverload();
        render(b, rest, 512, capture);
        CHECK_FALSE(b.takeOverload());
        CHECK(*std::max_element(seen.begin() + 1000, seen.end()) > 0.44f);
        b.reset();
        render(b, cosine(48000.0, 10000.0, 0.45, 9600), 512, capture);
        CHECK(b.takeOverload());
    }
}

TEST_CASE("Bridge: output is independent of the block split", "[bridge]")
{
    for (const auto fs : rates)
    {
        std::mt19937 rng(7);
        std::uniform_real_distribution<float> noise(-0.5f, 0.5f);
        std::vector<float> x((size_t)(0.3 * fs));
        for (auto& v : x)
            v = noise(rng);
        DeviceBridge b;
        b.prepare(fs, 512);
        const auto ref = render(b, x, 512, identity);
        for (const auto block : {1, 7, 333, 4096})
        {
            b.reset();
            const auto o = render(b, x, block, identity);
            CHECK(sameBits(o.l, ref.l));
            CHECK(sameBits(o.r, ref.r));
        }
        // Random splits.
        b.reset();
        Out o{std::vector<float>(x.size()), std::vector<float>(x.size())};
        std::uniform_int_distribution<int> len(0, 700);
        for (size_t i = 0; i < x.size();)
        {
            const auto n = (int)std::min((size_t)len(rng), x.size() - i);
            b.process(x.data() + i, o.l.data() + i, o.r.data() + i, n, identity);
            i += (size_t)n;
        }
        CHECK(sameBits(o.l, ref.l));
        CHECK(sameBits(o.r, ref.r));
    }
}

TEST_CASE("Bridge: CPU", "[bridge]")
{
    for (const auto fs : rates)
    {
        DeviceBridge b;
        b.prepare(fs, 512);
        std::mt19937 rng(3);
        std::uniform_real_distribution<float> noise(-0.5f, 0.5f);
        std::vector<float> x(512), l(512), r(512);
        for (auto& v : x)
            v = noise(rng);
        const auto blocks = (int)(10.0 * fs / 512);
        const auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < blocks; ++k)
            b.process(x.data(), l.data(), r.data(), 512, identity);
        const auto sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("Bridge %6.0f Hz: %.3f %% of one core (10 s of audio in %.1f ms)\n", fs, 100.0 * sec / 10.0, 1000.0 * sec);
        CHECK(sec < 1.0); // under 10 % of a core
    }
}
