#include "dsp/RefIr.h"
#include "dsp/generated/n_spec.h"
#include "dsp/programs/NCore.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
using rv::dsp::NCore;
using rv::dsp::Settings;
namespace gen = rv::dsp::gen::n;

Settings nSettings(double span, double early, double late)
{
    Settings s;
    s.program = rv::params::Program::n;
    s.amplPercent = span;
    s.lfX = early;
    s.hfX = late;
    return s;
}

struct RefSetting
{
    double span, early, late;
    std::string file() const
    {
        char b[96];
        std::snprintf(b, sizeof b, "analysis/nl/ref/nl_%.1f_%.1f_%.1f.wav", span, early, late);
        return b;
    }
};

// analysis/nl/nl_plugin_spec.json references.files
const std::array<RefSetting, 6> refs{{{0.0, 10.0, 1.0}, {0.2, 1.0, 1.0}, {0.5, -1.0, 2.0}, {0.3, 3.0, 0.0}, {0.4, -3.0, 0.6}, {0.1, 0.3, 1.4}}};

void configure(NCore& core, const Settings& s, bool kernel)
{
    core.setBranchKernel(kernel);
    core.prepare();
    core.setSettings(s);
    core.reset();
}
} // namespace

TEST_CASE("N: the generated tables agree with the engine", "[n]")
{
    CHECK(std::abs(gen::deviceRateHz - rv::dsp::deviceRate) < 1e-9);
    CHECK(gen::bkTapsL.size() == 2 * (size_t)gen::bkCentre + 1);
    CHECK(NCore().latency() == gen::bkCentre); // the adopted kernel is on by default
}

TEST_CASE("N: the SPAN rule gives the library's stage counts, and is monotonic in between", "[n]")
{
    for (size_t i = 0; i < gen::libSpan.size(); ++i)
    {
        CHECK(NCore::stagesForSpan(gen::libSpan[i]) == gen::libStages[i]);
        CHECK(NCore::feedFor(gen::libStages[i]) == gen::libFeed[i]);
    }
    int prev = 0;
    for (int i = 0; i <= 500; ++i)
    {
        const auto ns = NCore::stagesForSpan(i / 1000.0);
        CHECK(ns >= prev);
        CHECK(ns - prev <= 2);
        prev = ns;
        if (ns > 0)
            CHECK(NCore::feedFor(ns) > 0);
    }
    CHECK(NCore::stagesForSpan(0.25) == 5);
    CHECK(NCore::stagesForSpan(7.0) == NCore::numStages); // out of range is clamped
}

TEST_CASE("N: the law reproduces the Python bench's gains at all 96 library (SPAN, EARLY) pairs", "[n]")
{
    double worst = 0.0;
    size_t row = 0;
    const auto w = (size_t)gen::lawVectorRow;
    for (auto span : gen::spanValues)
        for (auto early : gen::earlyValues)
        {
            const auto* v = gen::lawVectors.data() + row++ * w;
            const auto law = NCore::earlyLaw(span, early);
            CHECK(law.stages == (int)v[0]);
            auto cmp = [&](double got, double want) {
                worst = std::max(worst, std::abs(got - want) / std::max(1.0, std::abs(want)));
            };
            for (size_t i = 0; i < 3; ++i)
            {
                cmp(law.frontDbL[i], v[1 + i]);
                cmp(law.frontDbR[i], v[4 + i]);
            }
            for (int j = 0; j < law.stages; ++j)
                cmp(law.stageDb[(size_t)j], v[7 + (size_t)j]);
            cmp(law.lateDb, v[7 + (size_t)NCore::numStages]);
        }
    CHECK(row == gen::spanValues.size() * gen::earlyValues.size());
    std::printf("N law vs bench, worst relative error over 96 pairs: %.3g\n", worst);
    CHECK(worst < 1e-12);
}

TEST_CASE("N: any value in range gives a finite law and a stable render", "[n]")
{
    for (double early : {-3.0, -1.7, -0.35, -0.01, 0.0, 0.05, 0.8, 2.6, 9.99, 10.0})
        for (double span : {0.0, 0.07, 0.25, 0.37, 0.5})
        {
            const auto law = NCore::earlyLaw(span, early);
            for (auto d : law.stageDb)
                CHECK(std::isfinite(d));
            CHECK(std::isfinite(law.lateDb));
        }
    NCore core;
    configure(core, nSettings(0.25, 0.8, 0.9), true);
    const auto ir = rv::test::renderImpulse(core, 20000);
    double e = 0.0;
    for (size_t i = 0; i < ir.size(); ++i)
        e += ir.l[i] * ir.l[i] + ir.r[i] * ir.r[i];
    CHECK(std::isfinite(e));
    CHECK(e > 0.0);
}

TEST_CASE("N: kernel bypassed, every reference IR nulls at -100 dB or lower, both channels, every window", "[n]")
{
    double worst = -400.0;
    for (const auto& r : refs)
    {
        const auto ref = rv::test::loadWav(r.file());
        REQUIRE(ref.size() > 60000);
        NCore core;
        configure(core, nSettings(r.span, r.early, r.late), false);
        REQUIRE(core.latency() == 0);
        const auto got = rv::test::renderImpulse(core, ref.size());
        const auto w = rv::test::worstNull(r.file().substr(16), ref, got);
        // the tail the spec windows stop short of (1.0 s to the end of the file)
        const auto tl = rv::test::nullDb(ref.l, got.l, 1000, 2500), tr = rv::test::nullDb(ref.r, got.r, 1000, 2500);
        std::printf("    1000-2500 ms %.1f/%.1f\n", tl, tr);
        CHECK(w <= -100.0);
        CHECK(std::max(tl, tr) <= -100.0);
        worst = std::max({worst, w, tl, tr});
    }
    std::printf("N worst null, kernel bypassed: %.1f dB\n", worst);
}

TEST_CASE("N: a SPAN a hair off a library step renders the same IR", "[n]")
{
    const auto ref = rv::test::loadWav(refs[5].file());
    NCore core;
    configure(core, nSettings(0.1 + 1e-9, 0.3 + 1e-12, 1.4), false);
    CHECK(rv::test::worstNull("0.1+eps", ref, rv::test::renderImpulse(core, ref.size())) <= -100.0);
}

TEST_CASE("N: kernel on, the look-ahead is reported and the level against the kernel-free references holds", "[n]")
{
    double worst = -400.0;
    for (const auto& r : refs)
    {
        const auto ref = rv::test::loadWav(r.file());
        NCore core;
        configure(core, nSettings(r.span, r.early, r.late), true);
        REQUIRE(core.latency() == 40);
        const auto got = rv::test::renderImpulse(core, ref.size()); // latency removed by the helper
        worst = std::max(worst, rv::test::worstNull(std::string("kernel on ") + r.file().substr(16), ref, got));
    }
    std::printf("N worst null, kernel on vs kernel-free refs: %.1f dB (spec: about -59 dB worst)\n", worst);
    CHECK(worst <= -57.0);
    CHECK(worst >= -75.0); // the kernel is doing something
}

TEST_CASE("N: the output does not depend on the block size, and reset repeats the render", "[n]")
{
    const auto s = nSettings(0.5, 0.7, 1.2);
    NCore core;
    configure(core, s, true);
    const auto a = rv::test::renderImpulse(core, 12000, 256);
    core.reset();
    const auto b = rv::test::renderImpulse(core, 12000, 256);
    core.reset();
    const auto c = rv::test::renderImpulse(core, 12000, 37);
    core.reset();
    const auto d = rv::test::renderImpulse(core, 12000, 1);
    double diff = 0.0; // bit-identical: the largest difference is 0
    for (size_t i = 0; i < a.size(); ++i)
        diff = std::max({diff, std::abs(a.l[i] - b.l[i]), std::abs(a.r[i] - b.r[i]), std::abs(a.l[i] - c.l[i]), std::abs(a.r[i] - d.r[i]),
                         std::abs(a.l[i] - d.l[i]), std::abs(a.r[i] - c.r[i])});
    CHECK(diff <= 0.0);
}

TEST_CASE("N: gating the late section off clears it, and settings can change while it runs", "[n]")
{
    NCore core;
    configure(core, nSettings(0.3, 1.0, 1.5), false);
    std::vector<float> in(4096, 0.0f), l(4096), r(4096);
    in[0] = 1.0f;
    core.process(in.data(), l.data(), r.data(), 4096);
    core.setSettings(nSettings(0.3, 1.0, 0.0)); // late off
    in[0] = 0.0f;
    for (int k = 0; k < 20; ++k)
        core.process(in.data(), l.data(), r.data(), 4096);
    core.setSettings(nSettings(0.3, -2.0, 1.5)); // late back on, other EARLY
    double e = 0.0;
    for (int k = 0; k < 4; ++k)
    {
        core.process(in.data(), l.data(), r.data(), 4096);
        for (size_t i = 0; i < 4096; ++i)
            e += (double)l[i] * l[i] + (double)r[i] * r[i];
    }
    CHECK(e <= 0.0); // nothing left in the loops from before the gate (and the main taps have run out)
}

// Run with "[.n-cost]": the cost of one device-rate sample at the heaviest setting (all 27 main taps, the late section).
TEST_CASE("N: cost per sample", "[.n-cost]")
{
    for (const bool kernel : {false, true})
    {
        NCore core;
        configure(core, nSettings(0.5, -1.0, 2.0), kernel);
        std::vector<float> in(1 << 16, 0.0f), l(in.size()), r(in.size());
        in[0] = 1.0f;
        const int reps = 40;
        const auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < reps; ++k)
        {
            core.process(in.data(), l.data(), r.data(), (int)in.size());
            in[0] = 0.0f;
        }
        const auto dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const auto ns = dt / (double)reps / (double)in.size() * 1e9;
        std::printf("N kernel %s: %.1f ns/sample = %.3f %% of one core at the device rate\n", kernel ? "on " : "off", ns,
                    ns * 1e-9 * rv::dsp::deviceRate * 100.0);
        CHECK(std::isfinite(l.back()));
    }
}

TEST_CASE("N: at LATE 0.2 the late loops run at the unit's measured gains; elsewhere at the RT law", "[n]")
{
    // NonLin_program 10.4: at LATE 0.2 the unit's three loop gains sit 2 / 5 / 12 % above 10^(-3 P / (fs LATE))
    // (measured pass by pass, held out on the bench); from LATE 0.4 the law holds.
    for (const double late : {0.2, 0.4, 1.0, 2.0})
    {
        NCore core;
        configure(core, nSettings(0.2, 10.0, late), false);
        for (size_t k = 0; k < gen::loopDelay.size(); ++k)
        {
            const auto law = std::pow(10.0, -3.0 * gen::loopDelay[k] / (rv::dsp::deviceRate * late));
            const auto want = juce::exactlyEqual(late, 0.2) ? gen::lateTableGain[k] : law;
            INFO("LATE " << late << ", loop " << k);
            CHECK(std::abs(core.lateLoopGain(k) - want) < 1e-12);
        }
    }
}
