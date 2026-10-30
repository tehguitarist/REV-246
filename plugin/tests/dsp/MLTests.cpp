#include "dsp/RefIr.h"
#include "dsp/programs/MLCore.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <random>
#include <string>
#include <vector>

// Programs M and L against analysis/ml/ml_plugin_spec.json and its reference IRs (docs/ENGINE.md, Tests).

namespace
{
using rv::dsp::MLCore;
using rv::params::Program;

juce::var loadSpec()
{
    return juce::JSON::parse(rv::test::repoFile("analysis/ml/ml_plugin_spec.json"));
}

rv::dsp::Settings settingsFor(Program p, double lf, double rt, double hf)
{
    rv::dsp::Settings s;
    s.program = p;
    s.lfX = lf;
    s.rt = rt;
    s.hfX = hf;
    s.amplPercent = 90.0; // the reference IRs' unit cluster
    s.sixteenBit = false;
    return s;
}

// A fresh core at one setting, prepared and reset.
void setUp(MLCore& core, Program p, double lf, double rt, double hf, bool modulation)
{
    core.setProgram(p);
    core.setModulation(modulation);
    core.prepare();
    core.setSettings(settingsFor(p, lf, rt, hf));
    core.reset();
}

double relErr(double got, double want)
{
    return std::abs(got - want) / std::max(std::abs(want), 1e-300);
}

double energy(const std::vector<double>& x, size_t a, size_t b)
{
    double e = 0.0;
    for (size_t i = a; i < std::min(b, x.size()); ++i)
        e += x[i] * x[i];
    return e;
}

// Bit-identical (the renders must repeat exactly, not approximately).
bool identical(double a, double b)
{
    return std::equal_to<double>{}(a, b);
}

size_t at(double seconds)
{
    return (size_t)(seconds * rv::dsp::deviceRate);
}
} // namespace

TEST_CASE("ML: the coefficient law reproduces every resolved setting", "[ml]")
{
    const auto spec = loadSpec();
    REQUIRE(spec.isObject());
    double worstGain = 0.0, worstDc = 0.0, worstPole = 0.0, worstB = 0.0;
    int count = 0;
    for (const auto* name : {"M", "L"})
    {
        const auto prog = std::string(name) == "L" ? Program::l : Program::m;
        const auto* settings = spec[name]["settings"].getDynamicObject();
        REQUIRE(settings != nullptr);
        for (const auto& kv : settings->getProperties())
        {
            const auto& v = kv.value;
            const auto c = MLCore::law(prog, (double)v["lf"], (double)v["rt"], (double)v["hf"]);
            for (int j = 0; j < 4; ++j)
            {
                const auto& line = v["lines"][j];
                const auto& got = c.lines[(size_t)j];
                worstGain = std::max(worstGain, relErr(got.gain, (double)line["gain"]));
                worstDc = std::max(worstDc, relErr(got.shelfDc, (double)line["shelf_dc"]));
                const double p = line["hf_pole"];
                worstPole = std::max(worstPole, juce::exactlyEqual(p, 0.0) ? std::abs(got.hfPole) : relErr(got.hfPole, p));
            }
            const double b = v["b"];
            worstB = std::max(worstB, juce::exactlyEqual(b, 0.0) ? std::abs(c.damping) : relErr(c.damping, b));
            ++count;
        }
    }
    std::printf("ML law vs %d resolved settings, max relative error: gain %.2e, shelf DC %.2e, HF pole %.2e, b %.2e\n", count,
                worstGain, worstDc, worstPole, worstB);
    CHECK(count == 340);
    CHECK(worstGain < 1e-9);
    CHECK(worstDc < 1e-9);
    CHECK(worstPole < 1e-9);
    CHECK(identical(worstB, 0.0));
}

TEST_CASE("ML: the law is continuous between the library's steps", "[ml]")
{
    // LF between S's steps (log-interpolated) and RT / HF off the steps give finite, stable coefficients.
    for (auto prog : {Program::m, Program::l})
        for (double lf : {0.1, 0.25, 0.7, 0.8, 1.5, 2.0})
            for (double rt : {0.6, 0.9, 2.2, 2.3, 7.0, 20.0})
                for (double hf : {0.1, 0.15, 0.55, 1.0})
                {
                    const auto c = MLCore::law(prog, lf, rt, hf);
                    for (const auto& l : c.lines)
                    {
                        CHECK(std::isfinite(l.gain));
                        CHECK(l.gain * l.shelfDc < 1.0); // the per-pass DC gain decays
                        CHECK(l.gain < 1.0);
                        CHECK(l.hfPole >= 0.0);
                        CHECK(l.hfPole < 0.86);
                    }
                }
}

TEST_CASE("ML: the static network (theta = 0) nulls every reference IR", "[ml]")
{
    const auto spec = loadSpec();
    MLCore core;
    double worst = -400.0;
    for (const auto* name : {"M", "L"})
    {
        const auto prog = std::string(name) == "L" ? Program::l : Program::m;
        const auto* files = spec[name]["references"]["files"].getDynamicObject();
        REQUIRE(files != nullptr);
        for (const auto& kv : files->getProperties())
        {
            const auto key = kv.name.toString();
            const auto parts = juce::StringArray::fromTokens(key, "_", "");
            const auto ref = rv::test::loadWav(kv.value["file"].toString().toStdString());
            REQUIRE(ref.size() > 1000);
            setUp(core, prog, parts[0].getDoubleValue(), parts[1].getDoubleValue(), parts[2].getDoubleValue(), false);
            const auto got = rv::test::renderImpulse(core, ref.size());
            const auto w = rv::test::worstNull(std::string(name) + " " + key.toStdString(), ref, got);
            // Beyond the spec windows: 1 s to the end of the file.
            const double tl = rv::test::nullDb(ref.l, got.l, 1000.0, 1e9), tr = rv::test::nullDb(ref.r, got.r, 1000.0, 1e9);
            std::printf("    1000 ms-end %.1f/%.1f\n", tl, tr);
            worst = std::max({worst, w, tl, tr});
        }
    }
    std::printf("ML static null, worst over every file, window and channel: %.1f dB\n", worst);
    CHECK(worst <= -100.0);
}

TEST_CASE("ML: theta has the spec's rate and depth, and F(theta) stays orthogonal", "[ml]")
{
    MLCore core;
    setUp(core, Program::m, 1.2, 1.8, 0.9, true);
    const auto n = at(24.1); // 20 periods of 0.83 Hz and a little
    std::vector<double> th(n);
    float in = 0.0f, l = 0.0f, r = 0.0f;
    for (size_t i = 0; i < n; ++i)
    {
        core.process(&in, &l, &r, 1);
        th[i] = core.theta();
    }
    // Rate from the upward zero crossings (linear between samples), rms over whole periods.
    std::vector<double> up;
    for (size_t i = 1; i < n; ++i)
        if (th[i - 1] < 0.0 && th[i] >= 0.0)
            up.push_back((double)(i - 1) + th[i - 1] / (th[i - 1] - th[i]));
    REQUIRE(up.size() >= 10);
    const double period = (up.back() - up.front()) / (double)(up.size() - 1);
    const double rate = rv::dsp::deviceRate / period;
    double ss = 0.0, peak = 0.0;
    for (auto i = (size_t)up.front(); i < (size_t)up.back(); ++i)
    {
        ss += th[i] * th[i];
        peak = std::max(peak, std::abs(th[i]));
    }
    const double rms = std::sqrt(ss / (double)((size_t)up.back() - (size_t)up.front()));
    std::printf("ML theta: rate %.5f Hz (spec 0.83), rms %.5f rad (spec 0.15), peak %.5f rad (spec 0.21)\n", rate, rms, peak);
    CHECK(std::abs(rate - rv::dsp::gen::ml::thetaRateHz) < 1e-4);
    CHECK(std::abs(rms - rv::dsp::gen::ml::thetaRms) < 1e-4);
    CHECK(std::abs(peak - rv::dsp::gen::ml::thetaAmplitude) < 1e-4);

    double worstOrth = 0.0;
    for (double t = -0.25; t <= 0.25; t += 0.01)
    {
        const auto F = MLCore::couplingAt(t);
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = 0; j < 4; ++j)
            {
                double d = 0.0;
                for (size_t k = 0; k < 4; ++k)
                    d += F[i * 4 + k] * F[j * 4 + k];
                worstOrth = std::max(worstOrth, std::abs(d - (i == j ? 1.0 : 0.0)));
            }
    }
    std::printf("ML F(theta) F(theta)^T - I, worst entry over |theta| <= 0.25: %.2e\n", worstOrth);
    CHECK(worstOrth < 1e-14);
    // theta = 0 is the spec's F exactly.
    const auto F0 = MLCore::couplingAt(0.0);
    for (size_t i = 0; i < 16; ++i)
        CHECK(identical(F0[i], rv::dsp::gen::ml::coupling[i]));
}

TEST_CASE("ML: modulated renders repeat, start as the static network, decorrelate later, keep the energy", "[ml]")
{
    struct Case
    {
        Program prog;
        double lf, rt, hf, seconds;
    };
    for (const auto& c : {Case{Program::m, 1.2, 5.0, 0.9, 8.0}, Case{Program::l, 2.0, 20.0, 1.0, 12.0}})
    {
        const auto n = at(c.seconds);
        MLCore core(c.prog);
        setUp(core, c.prog, c.lf, c.rt, c.hf, false);
        const auto fixed = rv::test::renderImpulse(core, n);
        setUp(core, c.prog, c.lf, c.rt, c.hf, true);
        const auto moving = rv::test::renderImpulse(core, n, 256);
        core.reset();
        const auto again = rv::test::renderImpulse(core, n, 97); // another block size: the same samples
        CHECK(again.l == moving.l);
        CHECK(again.r == moving.r);

        // Before the first recirculation reaches a read tap the outputs are identical.
        const auto& t = c.prog == Program::l ? rv::dsp::gen::ml::l::lineDelay : rv::dsp::gen::ml::m::lineDelay;
        const auto& tapL = c.prog == Program::l ? rv::dsp::gen::ml::l::tapLDelay : rv::dsp::gen::ml::m::tapLDelay;
        const auto& tapR = c.prog == Program::l ? rv::dsp::gen::ml::l::tapRDelay : rv::dsp::gen::ml::m::tapRDelay;
        const auto first = (size_t)(*std::min_element(t.begin(), t.end())
                                    + std::min(*std::min_element(tapL.begin(), tapL.end()), *std::min_element(tapR.begin(), tapR.end())));
        bool same = true;
        for (size_t i = 0; i < first; ++i)
            same = same && identical(fixed.l[i], moving.l[i]) && identical(fixed.r[i], moving.r[i]);
        CHECK(same);

        std::printf("ML %s %g_%g_%g, theta on vs off: window, null L/R dB, energy ratio L/R dB\n",
                    c.prog == Program::l ? "L" : "M", c.lf, c.rt, c.hf);
        const std::vector<std::array<double, 2>> wins{{0.15, 0.3}, {0.3, 0.6}, {0.6, 1.0}, {1.0, 2.0}, {2.0, 4.0}, {4.0, c.seconds}};
        for (const auto& w : wins)
        {
            const auto a = at(w[0]), b = at(w[1]);
            const double nl = rv::test::nullDb(fixed.l, moving.l, w[0] * 1000.0, w[1] * 1000.0);
            const double nr = rv::test::nullDb(fixed.r, moving.r, w[0] * 1000.0, w[1] * 1000.0);
            const double el = 10.0 * std::log10(energy(moving.l, a, b) / energy(fixed.l, a, b));
            const double er = 10.0 * std::log10(energy(moving.r, a, b) / energy(fixed.r, a, b));
            std::printf("  %4.2f-%5.2f s  %6.1f/%6.1f  %+5.2f/%+5.2f\n", w[0], w[1], nl, nr, el, er);
            if (w[0] >= 1.0)
                CHECK(std::max(nl, nr) > -20.0); // the moving coupling has decorrelated the tail
            CHECK(std::abs(el) < 1.5);           // and kept its energy (F(theta) is orthogonal)
            CHECK(std::abs(er) < 1.5);
        }
    }
}

TEST_CASE("ML: noise through L at RT 20 s with theta on stays bounded and keeps its level", "[ml]")
{
    const auto n = at(30.0);
    std::mt19937 rng(246);
    std::uniform_real_distribution<float> u(-0.5f, 0.5f);
    std::vector<float> in(n);
    for (auto& v : in)
        v = u(rng);
    std::array<double, 2> rmsTail{};
    double peak = 0.0;
    bool finite = true;
    for (int mod = 0; mod < 2; ++mod)
    {
        MLCore core(Program::l);
        setUp(core, Program::l, 2.0, 20.0, 1.0, mod == 1);
        std::vector<float> l(n), r(n);
        for (size_t i = 0; i < n; i += 512)
        {
            const auto m = (int)std::min<size_t>(512, n - i);
            core.process(in.data() + i, l.data() + i, r.data() + i, m);
        }
        double ss = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            finite = finite && std::isfinite(l[i]) && std::isfinite(r[i]);
            if (mod == 1)
                peak = std::max({peak, (double)std::abs(l[i]), (double)std::abs(r[i])});
            if (i >= at(20.0))
                ss += (double)l[i] * l[i] + (double)r[i] * r[i];
        }
        rmsTail[(size_t)mod] = std::sqrt(ss / (double)(2 * (n - at(20.0))));
    }
    const double ratio = 20.0 * std::log10(rmsTail[1] / rmsTail[0]);
    std::printf("ML L 2.0_20.0_1.0, 30 s noise: steady rms theta off %.4f, on %.4f (%+.2f dB), peak on %.3f\n", rmsTail[0],
                rmsTail[1], ratio, peak);
    CHECK(finite);
    CHECK(std::abs(ratio) < 0.5);
    CHECK(peak < 50.0);
}

TEST_CASE("ML: switching the program matches a fresh core", "[ml]")
{
    MLCore a(Program::m), b(Program::l);
    setUp(a, Program::m, 1.2, 1.8, 0.9, false);
    rv::test::renderImpulse(a, 5000);
    setUp(a, Program::l, 1.2, 6.0, 0.8, false);
    setUp(b, Program::l, 1.2, 6.0, 0.8, false);
    const auto ya = rv::test::renderImpulse(a, at(1.0)), yb = rv::test::renderImpulse(b, at(1.0));
    CHECK(ya.l == yb.l);
    CHECK(ya.r == yb.r);
}

TEST_CASE("ML: cost per device sample", "[ml][.bench]")
{
    const auto n = at(20.0);
    std::vector<float> in(n), l(n), r(n);
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> u(-0.5f, 0.5f);
    for (auto& v : in)
        v = u(rng);
    for (auto prog : {Program::m, Program::l})
    {
        MLCore core(prog);
        setUp(core, prog, 1.2, prog == Program::l ? 6.0 : 1.8, 0.9, true);
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < n; i += 256)
            core.process(in.data() + i, l.data() + i, r.data() + i, (int)std::min<size_t>(256, n - i));
        const auto t1 = std::chrono::steady_clock::now();
        const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / (double)n;
        std::printf("ML %s: %.1f ns per device sample (%.2f %% of one core in real time)\n", prog == Program::l ? "L" : "M", ns,
                    ns * rv::dsp::deviceRate * 1e-7);
        CHECK(ns > 0.0);
    }
}
