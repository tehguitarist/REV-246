#include "dsp/RefIr.h"
#include "dsp/programs/ECore.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace
{
namespace g = rv::dsp::gen::e;

rv::dsp::Settings eSettings(double rt, double lf = g::familyLf, double hf = g::familyHf, double ampl = 80.0)
{
    rv::dsp::Settings s;
    s.program = rv::params::Program::e;
    s.rt = rt;
    s.lfX = lf;
    s.hfX = hf;
    s.amplPercent = ampl;
    return s;
}

void runCore(rv::dsp::ECore& core, const std::vector<float>& in, std::vector<float>& l, std::vector<float>& r)
{
    l.assign(in.size(), 0.0f);
    r.assign(in.size(), 0.0f);
    for (size_t i = 0; i < in.size(); i += 256)
    {
        const auto m = (int)std::min<size_t>(256, in.size() - i);
        core.process(in.data() + i, l.data() + i, r.data() + i, m);
    }
}

// The impulse response after `preroll` samples of silence (a different phase of the moving reads), as doubles.
rv::test::StereoIr impulseAfter(rv::dsp::ECore& core, size_t preroll, size_t n)
{
    std::vector<float> in(preroll + n, 0.0f), l, r;
    in[preroll] = 1.0f;
    runCore(core, in, l, r);
    rv::test::StereoIr ir;
    ir.l.assign(l.begin() + (long)preroll, l.end());
    ir.r.assign(r.begin() + (long)preroll, r.end());
    return ir;
}

double fullNull(const std::vector<double>& ref, const std::vector<double>& got)
{
    return rv::test::nullDb(ref, got, 0.0, (double)ref.size() / rv::dsp::deviceRate * 1000.0 + 1.0);
}
} // namespace

TEST_CASE("E: the loop law reproduces ml_lines.line_coefs", "[e]")
{
    constexpr size_t nv = g::lawVectorIn.size() / 3;
    double worst = 0.0;
    for (size_t v = 0; v < nv; ++v)
        for (int k = 0; k < rv::dsp::ECore::numLoops; ++k)
        {
            const auto c = rv::dsp::ECore::loopLaw(k, g::lawVectorIn[3 * v], g::lawVectorIn[3 * v + 1], g::lawVectorIn[3 * v + 2]);
            const auto* want = &g::lawVectorOut[(v * rv::dsp::ECore::numLoops + (size_t)k) * 3];
            const double got[3] = {c.g, c.r0, c.p};
            for (int j = 0; j < 3; ++j)
                worst = std::max(worst, std::abs(got[j] - want[j]) / std::max(std::abs(want[j]), 1e-3));
        }
    std::printf("E law: %zu settings x 3 loops, worst relative error %.2e\n", nv, worst);
    CHECK(worst < 1e-9);
}

TEST_CASE("E: the static network nulls the reference IRs", "[e]")
{
    for (const double rt : {0.4, 1.0, 2.5, 4.5})
    {
        char name[64];
        std::snprintf(name, sizeof name, "analysis/e250/ref/e_static_rt%.1f.wav", rt);
        const auto ref = rv::test::loadWav(name);
        REQUIRE(ref.size() > 10000);
        rv::dsp::ECore core;
        core.setModulation(false);
        core.prepare();
        core.setSettings(eSettings(rt));
        core.reset();
        const auto got = rv::test::renderImpulse(core, ref.size());
        const double worst = rv::test::worstNull(name, ref, got);
        std::printf("  whole file (%.2f s): L %.1f / R %.1f dB\n", (double)ref.size() / rv::dsp::deviceRate,
                    fullNull(ref.l, got.l), fullNull(ref.r, got.r));
        CHECK(worst <= -100.0);
    }
}

TEST_CASE("E: AMPL scales the direct taps only", "[e]")
{
    rv::dsp::ECore a, b;
    for (auto* c : {&a, &b})
    {
        c->setModulation(false);
        c->prepare();
    }
    a.setSettings(eSettings(1.0, g::familyLf, g::familyHf, 80.0));
    b.setSettings(eSettings(1.0, g::familyLf, g::familyHf, 40.0));
    a.reset();
    b.reset();
    const auto ya = rv::test::renderImpulse(a, 4000), yb = rv::test::renderImpulse(b, 4000);
    CHECK(std::abs(yb.r[(size_t)g::k0] - 0.5 * ya.r[(size_t)g::k0]) < 1e-7);
    // the first static read (loop B at 805) is the reverb's first arrival: unchanged
    const auto j = (size_t)(g::k0 + 805 + 146 + 46 + 392);
    CHECK(std::abs(ya.r[j]) > 1e-4);
    CHECK(std::abs(yb.r[j] - ya.r[j]) < 1e-9);
}

TEST_CASE("E: moving reads follow e_sim's law (depth, rate, integer addresses)", "[e]")
{
    rv::dsp::ECore core;
    core.prepare();
    core.setSettings(eSettings(2.5));
    core.reset();
    const double secs = 120.0;
    const auto n = (size_t)(secs * rv::dsp::deviceRate);
    const auto nr = (size_t)g::numReads;
    std::vector<int> lo(nr, 1 << 30), hi(nr, -(1 << 30)), ups(nr, 0), prev(nr);
    for (int r = 0; r < g::numReads; ++r)
        prev[(size_t)r] = core.readOffset(r);
    const float zero = 0.0f;
    float l = 0.0f, rr = 0.0f;
    bool stepOk = true;
    for (size_t i = 0; i < n; ++i)
    {
        core.process(&zero, &l, &rr, 1);
        for (int r = 0; r < g::numReads; ++r)
        {
            const auto o = core.readOffset(r), c = g::readCentre[(size_t)r];
            lo[(size_t)r] = std::min(lo[(size_t)r], o);
            hi[(size_t)r] = std::max(hi[(size_t)r], o);
            ups[(size_t)r] += (prev[(size_t)r] < c && o >= c) ? 1 : 0;
            stepOk = stepOk && std::abs(o - prev[(size_t)r]) <= 1; // one sample per step at most (triangle: 4 A f / fs = 0.004 at most)
            prev[(size_t)r] = o;
        }
    }
    CHECK(stepOk);
    double meanRate = 0.0;
    for (size_t r = 0; r < (size_t)g::numReads; ++r)
    {
        const double rate = ups[r] / secs;
        meanRate += g::readRateHz[r] / g::numReads;
        std::printf("E read %zu: centre %d, range %+d..%+d, rate %.3f Hz measured / %.4f Hz set\n", r, g::readCentre[r],
                    lo[r] - g::readCentre[r], hi[r] - g::readCentre[r], rate, g::readRateHz[r]);
        CHECK(lo[r] == g::readCentre[r] - (int)g::readDepth);
        CHECK(hi[r] == g::readCentre[r] + (int)g::readDepth);
        CHECK(std::abs(rate - g::readRateHz[r]) < 1.5 / secs);
        CHECK(g::readRateHz[r] >= g::readRate * (1.0 - g::readRateSpread));
        CHECK(g::readRateHz[r] <= g::readRate * (1.0 + g::readRateSpread));
    }
    std::printf("E reads: mean rate %.3f Hz (law: 0.553 x U(0.5, 1.5))\n", meanRate);
}

TEST_CASE("E: moving reads match e_sim.render (triangle, seed 0)", "[e]")
{
    // analysis/e250/ref/e_moving_rt*.wav (tools/plugin_codegen/e_moving_ref.py): e_sim.render with ADOPTED, seed 0 (the
    // realisation in e_spec.h), the impulse at network time k0 so the reads' clock is the core's since reset()
    for (const double rt : {1.0, 4.5})
    {
        char name[64];
        std::snprintf(name, sizeof name, "analysis/e250/ref/e_moving_rt%.1f.wav", rt);
        REQUIRE(rv::test::repoFile(name).existsAsFile());
        const auto ref = rv::test::loadWav(name);
        REQUIRE(ref.size() > 10000);
        rv::dsp::ECore core;
        core.prepare();
        core.setSettings(eSettings(rt));
        core.reset();
        const auto got = rv::test::renderImpulse(core, ref.size());
        const double worst = rv::test::worstNull(name, ref, got);
        std::printf("  whole file: L %.1f / R %.1f dB\n", fullNull(ref.l, got.l), fullNull(ref.r, got.r));
        CHECK(worst <= -100.0);
    }
}

TEST_CASE("E: moving / static energy ratios against the spec's targets", "[e]")
{
    // e_sim.stats' e_ratio on R: (moving IR - static IR) energy over the static IR's, windows in samples after R's
    // tap. The targets are sweep captures' (spec references.targets); impulse renders at eight phases of the reads.
    const std::array<std::array<int, 2>, 4> win{{{300, 805}, {805, 2215}, {2215, 4400}, {4400, 9000}}};
    const std::array<double, 8> rts{{1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5}};
    std::array<double, 4> mean{};
    for (size_t i = 0; i < rts.size(); ++i)
    {
        rv::dsp::ECore mov, stat;
        stat.setModulation(false);
        for (auto* c : {&mov, &stat})
        {
            c->prepare();
            c->setSettings(eSettings(rts[i]));
            c->reset();
        }
        const auto pre = (size_t)((0.37 + 1.13 * (double)i) * rv::dsp::deviceRate);
        const auto ym = impulseAfter(mov, pre, 9200 + (size_t)g::k0), ys = impulseAfter(stat, pre, 9200 + (size_t)g::k0);
        for (size_t w = 0; w < win.size(); ++w)
        {
            double e = 0.0, s = 0.0;
            for (auto t = (size_t)(g::k0 + win[w][0]); t < (size_t)(g::k0 + win[w][1]); ++t)
            {
                e += (ym.r[t] - ys.r[t]) * (ym.r[t] - ys.r[t]);
                s += ys.r[t] * ys.r[t];
            }
            mean[w] += 10.0 * std::log10(e / s) / (double)rts.size();
        }
    }
    const std::array<double, 4> target{{-15.715453576478803, 4.3253274602586895, 5.439992679370786, 4.783781565105903}};
    const std::array<double, 4> sd{{0.30994565407781627, 1.8651687091555038, 1.322350727495855, 1.3491773424704672}};
    for (size_t w = 0; w < win.size(); ++w)
    {
        std::printf("E e_ratio %d-%d: %+.2f dB (target %+.2f, sd %.2f)\n", win[w][0], win[w][1], mean[w], target[w], sd[w]);
        CHECK(std::abs(mean[w] - target[w]) < 2.0 * sd[w] + 1.0);
    }
}

TEST_CASE("E: bounded and decaying at the extremes, moving reads on", "[e]")
{
    std::mt19937 rng(7);
    std::normal_distribution<float> noise(0.0f, 0.1f);
    const auto n = (size_t)(20.0 * rv::dsp::deviceRate);
    struct Case
    {
        double rt, lf, hf;
    };
    for (const auto& c : {Case{5.0, 2.0, 1.0}, Case{5.0, 100.0, 100.0}, Case{0.4, 0.1, 0.05}, Case{5.0, 0.1, 0.1}})
    {
        rv::dsp::ECore core;
        core.prepare();
        core.setSettings(eSettings(c.rt, c.lf, c.hf, 100.0));
        core.reset();
        std::vector<float> in(2 * n, 0.0f), l, r;
        for (size_t i = 0; i < n; ++i)
            in[i] = noise(rng);
        runCore(core, in, l, r);
        double peak = 0.0, tail = 0.0;
        bool finite = true;
        for (size_t i = 0; i < 2 * n; ++i)
        {
            finite = finite && std::isfinite(l[i]) && std::isfinite(r[i]);
            peak = std::max({peak, (double)std::abs(l[i]), (double)std::abs(r[i])});
            if (i >= 2 * n - (size_t)rv::dsp::deviceRate)
                tail = std::max({tail, (double)std::abs(l[i]), (double)std::abs(r[i])});
        }
        std::printf("E RT %g LF %g HF %g: peak %.3f, last second %.2e\n", c.rt, c.lf, c.hf, peak, tail);
        CHECK(finite);
        CHECK(peak < 10.0);
        CHECK(tail < 1e-6);
    }
}

TEST_CASE("E: cost per sample", "[e][cost]")
{
    rv::dsp::ECore core;
    core.prepare();
    core.setSettings(eSettings(2.5));
    core.reset();
    std::mt19937 rng(1);
    std::normal_distribution<float> noise(0.0f, 0.1f);
    const auto n = (size_t)(60.0 * rv::dsp::deviceRate);
    std::vector<float> in(n), l, r;
    for (auto& v : in)
        v = noise(rng);
    const auto t0 = std::chrono::steady_clock::now();
    runCore(core, in, l, r);
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double sink = 0.0;
    for (size_t i = 0; i < n; i += 997)
        sink += l[i] + r[i];
    std::printf("E cost: %.1f ns per device sample, %.3f %% of one core in real time (checksum %.3g)\n", s / (double)n * 1e9,
                s / 60.0 * 100.0, sink);
    CHECK(s < 60.0);
}
