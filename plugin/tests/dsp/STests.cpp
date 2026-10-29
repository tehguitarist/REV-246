#include "dsp/RefIr.h"
#include "dsp/programs/SCore.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

// Program S against its spec (analysis/s/s_plugin_spec.json): the coefficient law against the 128 resolved settings,
// the reference IRs (tools/s_spec.py simulate()) nulled in the spec windows, and the law between the library's steps.

namespace
{
juce::var loadSpec()
{
    return juce::JSON::parse(rv::test::repoFile("analysis/s/s_plugin_spec.json"));
}

std::vector<double> numbers(const juce::var& a)
{
    std::vector<double> v;
    for (const auto& x : *a.getArray())
        v.push_back((double)x);
    return v;
}

std::vector<double> convolve(const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<double> out(x.size() + y.size() - 1, 0.0);
    for (size_t i = 0; i < x.size(); ++i)
        for (size_t j = 0; j < y.size(); ++j)
            out[i + j] += x[i] * y[j];
    return out;
}

// G_i as direct-form polynomials (b with the gain, a with a[0] = 1), as s_spec.loop_poly.
void polynomial(const rv::dsp::SCore::LoopFilter& f, std::vector<double>& b, std::vector<double>& a)
{
    b = {f.gain};
    a = {1.0};
    for (int k = 0; k < f.count; ++k)
    {
        const auto& s = f.sections[(size_t)k];
        b = convolve(b, {s.b0, s.b1});
        a = convolve(a, {1.0, s.a1});
    }
}

// max |got - ref| / max |ref|, the shorter one padded with zeros.
double relError(std::vector<double> got, std::vector<double> ref)
{
    const auto n = std::max(got.size(), ref.size());
    got.resize(n, 0.0);
    ref.resize(n, 0.0);
    double e = 0.0, m = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        e = std::max(e, std::abs(got[i] - ref[i]));
        m = std::max(m, std::abs(ref[i]));
    }
    return e / m;
}

rv::dsp::Settings settingsFor(double lf, double rt, double hf, double ampl = 90.0)
{
    rv::dsp::Settings s;
    s.program = rv::params::Program::s;
    s.lfX = lf;
    s.rt = rt;
    s.hfX = hf;
    s.amplPercent = ampl;
    return s;
}
} // namespace

TEST_CASE("S: the coefficient law reproduces every resolved setting of the spec", "[s]")
{
    const auto spec = loadSpec();
    auto* settings = spec["settings"].getDynamicObject();
    REQUIRE(settings != nullptr);
    REQUIRE(settings->getProperties().size() == 128);
    double worstB = 0.0, worstA = 0.0, worstAp = 0.0;
    for (const auto& kv : settings->getProperties())
    {
        const auto& v = kv.value;
        const auto c = rv::dsp::SCore::law((double)v["lf"], (double)v["rt"], (double)v["hf"]);
        const auto& loops = *v["loops"].getArray();
        REQUIRE(loops.size() == 4);
        for (size_t i = 0; i < 4; ++i)
        {
            std::vector<double> b, a;
            polynomial(c.loops[i], b, a);
            const auto eb = relError(b, numbers(loops[(int)i]["b"])), ea = relError(a, numbers(loops[(int)i]["a"]));
            if (eb > 1e-9 || ea > 1e-9)
                std::printf("S law: %s loop %zu b %.2e a %.2e\n", kv.name.toString().toRawUTF8(), i, eb, ea);
            worstB = std::max(worstB, eb);
            worstA = std::max(worstA, ea);
        }
        const double b230 = v["ap230_b"];
        worstAp = std::max(worstAp, std::abs(c.ap230b - b230) / std::max(b230, 1e-300));
    }
    std::printf("S law vs 128 resolved settings: max relative error b %.2e, a %.2e, ap230_b %.2e\n", worstB, worstA, worstAp);
    CHECK(worstB < 1e-9);
    CHECK(worstA < 1e-9);
    CHECK(worstAp < 1e-12);
}

TEST_CASE("S: the reference IRs null at -100 dB or lower in every window", "[s]")
{
    const auto spec = loadSpec();
    auto* files = spec["references"]["files"].getDynamicObject();
    REQUIRE(files != nullptr);
    REQUIRE(files->getProperties().size() == 9);
    double worst = -400.0;
    for (const auto& kv : files->getProperties())
    {
        const auto& r = kv.value;
        const auto ref = rv::test::loadWav(r["file"].toString().toStdString());
        REQUIRE(ref.size() > 10000);
        const auto& v = spec["settings"][juce::Identifier(r["setting"].toString())];
        rv::dsp::SCore core;
        core.prepare();
        core.setSettings(settingsFor(v["lf"], v["rt"], v["hf"], r["ampl"]));
        core.reset();
        const auto got = rv::test::renderImpulse(core, ref.size());
        worst = std::max(worst, rv::test::worstNull(kv.name.toString().toStdString(), ref, got));
    }
    std::printf("S worst null over all references, windows and channels: %.1f dB\n", worst);
    CHECK(worst <= -100.0);
}

TEST_CASE("S: the law between the library's steps gives stable loops at every panel step", "[s]")
{
    using rv::params::Joy;
    using rv::params::Program;
    constexpr double twoPi = juce::MathConstants<double>::twoPi;
    double worstGain = 0.0;
    int count = 0;
    for (int il = 0; il <= rv::params::joyMaxStep(Joy::lf); ++il)
        for (int ir = 0; ir <= rv::params::joyMaxStep(Joy::rt); ++ir)
            for (int ih = 0; ih <= rv::params::joyMaxStep(Joy::hf); ++ih)
            {
                const auto c = rv::dsp::SCore::law(rv::params::joyNumber(Program::s, Joy::lf, il),
                                                   rv::params::joyNumber(Program::s, Joy::rt, ir),
                                                   rv::params::joyNumber(Program::s, Joy::hf, ih));
                ++count;
                for (const auto& f : c.loops)
                {
                    for (int k = 0; k < f.count; ++k)
                        REQUIRE(std::abs(f.sections[(size_t)k].a1) < 1.0);
                    // |G_i| < 1 at every frequency: the loop decays.
                    for (int j = 0; j <= 64; ++j)
                    {
                        const auto z1 = std::polar(1.0, -twoPi * 0.5 * j / 64.0);
                        std::complex<double> h = f.gain;
                        for (int k = 0; k < f.count; ++k)
                        {
                            const auto& s = f.sections[(size_t)k];
                            h *= (s.b0 + s.b1 * z1) / (1.0 + s.a1 * z1);
                        }
                        worstGain = std::max(worstGain, std::abs(h));
                    }
                }
                REQUIRE(c.ap230b >= 0.0);
                REQUIRE(c.ap230b < 0.5);
            }
    std::printf("S law at %d panel settings: largest loop gain %.4f\n", count, worstGain);
    CHECK(worstGain < 1.0);
}

TEST_CASE("S: settings change while running without blowing up, and reset repeats the render", "[s]")
{
    rv::dsp::SCore core;
    core.prepare();
    core.setSettings(settingsFor(0.1, 2.5, 1.0));
    core.reset();
    const auto a = rv::test::renderImpulse(core, 20000);
    juce::Random rng(246);
    std::vector<float> in(256), l(256), r(256);
    for (int block = 0; block < 400; ++block)
    {
        core.setSettings(settingsFor(0.1 + 1.9 * rng.nextDouble(), 0.2 + 2.3 * rng.nextDouble(), 0.1 + 0.9 * rng.nextDouble(),
                                     100.0 * rng.nextDouble()));
        for (auto& x : in)
            x = (float)(rng.nextDouble() * 2.0 - 1.0);
        core.process(in.data(), l.data(), r.data(), (int)in.size());
        for (size_t i = 0; i < in.size(); ++i)
            REQUIRE((std::isfinite(l[i]) && std::isfinite(r[i]) && std::abs(l[i]) < 100.0f && std::abs(r[i]) < 100.0f));
    }
    core.setSettings(settingsFor(0.1, 2.5, 1.0));
    core.reset();
    const auto b = rv::test::renderImpulse(core, 20000);
    CHECK(a.l == b.l);
    CHECK(a.r == b.r);
}
