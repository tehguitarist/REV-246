#include "dsp/engine/AnalogueHiss.h"
#include "dsp/engine/ProgramCore.h"
#include "dsp/generated/sixteen_ref.h"
#include "dsp/programs/ECore.h"
#include "dsp/programs/MLCore.h"
#include "dsp/programs/NCore.h"
#include "dsp/programs/SCore.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

// The 16-bit storage mode (docs/ENGINE.md, 16-bit; model_s_results §36.5) and the analogue hiss: ProgramCore's
// scheme against tools/loop_read/c/s_loops_fx.c bit for bit, float mode untouched, and per program the 16-bit
// error against float, the idle tail (dead band, limit cycle) and the idle noise with the hiss against the spec's
// REVERB S/N (73 dB unweighted re nominal; reported, not a pass bar).

namespace
{
using rv::dsp::deviceRate;

constexpr double lsb = 1.0 / 32768.0; // one A/D LSB at full scale 1.0 (the tests' setStorageFullScale)
const double nominalRms = std::pow(10.0, rv::dsp::AnalogueHiss::nominalDb / 20.0) / std::sqrt(2.0);

double dbOf(double x) { return 20.0 * std::log10(std::max(x, 1e-300)); }

uint64_t fnv1a(uint64_t h, int32_t w)
{
    auto u = (uint32_t)w;
    for (int b = 0; b < 4; ++b)
    {
        h ^= (u >> (8 * b)) & 0xFFu;
        h *= 0x100000001B3ull;
    }
    return h;
}

// s_loops_fx.c's network, built from ProgramCore's 16-bit primitives (the same store / cut / q15 the cores use).
class LoopsFx final : public rv::dsp::ProgramCore
{
  public:
    explicit LoopsFx(int c) : kase(c) {}
    void prepare() override
    {
        for (auto& d : lines)
            d.assign(2048, 0.0f);
    }
    void reset() override {}
    void setSettings(const rv::dsp::Settings&) override {}
    void process(const float*, float*, float*, int) override {}

    // Runs the case's input; returns the output's 16-bit words.
    std::vector<int32_t> run()
    {
        namespace g = rv::dsp::gen::sixteen;
        quantise = true;
        setStorageCut(g::mode[(size_t)kase] == 2 ? StorageCut::nearest : StorageCut::floor);
        setStorageFullScale(1.0);
        setGridGain(0, g::inputGain[(size_t)kase]); // the loops' input gain b: their stored word = b x the A/D word
        std::array<double, 4> c0{}, c1{}, a{}, p{}, p1{}, s{}, y{}, rp{};
        for (size_t i = 0; i < 4; ++i)
        {
            const auto* k = &g::coef[((size_t)kase * 4 + i) * 4];
            c0[i] = q15(k[0]), c1[i] = q15(k[1]), a[i] = q15(k[2]), p[i] = q15(k[3]), p1[i] = q15(1.0 - k[3]);
        }
        const auto n = (size_t)g::length[(size_t)kase];
        auto seed = g::inputSeed[(size_t)kase];
        const int amp = g::inputAmp[(size_t)kase];
        std::vector<int32_t> words(n);
        for (size_t t = 0; t < n; ++t)
        {
            int xi = 0; // the A/D word (the inputs are whole LSBs, so s_loops_fx.c's llround(x 2^30) >> 15 is exact)
            if (g::impulse[(size_t)kase] != 0)
                xi = t == 0 ? amp : 0;
            else
            {
                seed ^= seed << 13;
                seed ^= seed >> 17;
                seed ^= seed << 5;
                xi = (int)(seed % (uint32_t)(2 * amp + 1)) - amp;
            }
            const double x = xi * lsb;
            double o = 0.0;
            for (size_t i = 0; i < 4; ++i)
            {
                auto& d = lines[i];
                const double r = d[(t - (size_t)g::loopDelay[i]) & 2047];
                s[i] = cut16(c0[i] * r + c1[i] * rp[i] + a[i] * s[i], 0);
                rp[i] = r;
                y[i] = cut16(p1[i] * s[i] + p[i] * y[i], 0);
                d[t & 2047] = store(x + y[i], 0);
                o += g::loopSign[i] * r;
            }
            words[t] = (int32_t)std::llround(cut16(o, 0) / gridStep(0)); // the output word, saturated as s_loops_fx.c
        }
        return words;
    }

  private:
    int kase;
    std::array<std::vector<float>, 4> lines;
};

// ---- the program cores ----
struct Program
{
    std::string name;
    rv::params::Program prog;
    rv::dsp::Settings set;
};

using Cut = rv::dsp::ProgramCore::StorageCut;

std::unique_ptr<rv::dsp::ProgramCore> makeCore(const Program& p, bool sixteen, Cut cut = Cut::floor)
{
    using P = rv::params::Program;
    std::unique_ptr<rv::dsp::ProgramCore> c;
    if (p.prog == P::s)
        c = std::make_unique<rv::dsp::SCore>();
    else if (p.prog == P::m || p.prog == P::l)
        c = std::make_unique<rv::dsp::MLCore>(p.prog);
    else if (p.prog == P::e)
        c = std::make_unique<rv::dsp::ECore>();
    else
        c = std::make_unique<rv::dsp::NCore>();
    c->prepare();
    c->setStorageFullScale(1.0);
    c->setStorageCut(cut);
    auto s = p.set;
    s.program = p.prog;
    s.sixteenBit = sixteen;
    c->setSettings(s);
    c->reset();
    return c;
}

struct Stereo
{
    std::vector<double> l, r;
};

Stereo render(rv::dsp::ProgramCore& core, const std::vector<float>& in)
{
    const auto lat = (size_t)core.latency();
    const size_t total = in.size() + lat;
    std::vector<float> x(total, 0.0f), l(total), r(total);
    std::copy(in.begin(), in.end(), x.begin());
    for (size_t i = 0; i < total; i += 256)
    {
        const auto m = (int)std::min<size_t>(256, total - i);
        core.process(x.data() + i, l.data() + i, r.data() + i, m);
    }
    Stereo out;
    out.l.assign(l.begin() + (long)lat, l.end());
    out.r.assign(r.begin() + (long)lat, r.end());
    return out;
}

// The analogue chain's DC blocking for the measurements: a second-order Butterworth high-pass at 20 Hz (§36.5).
std::vector<double> highPass20(const std::vector<double>& x)
{
    const double pi = std::acos(-1.0), k = std::tan(pi * 20.0 / deviceRate), q = std::sqrt(0.5);
    const double norm = 1.0 / (1.0 + k / q + k * k);
    const double b0 = norm, b1 = -2.0 * norm, b2 = norm, a1 = 2.0 * (k * k - 1.0) * norm, a2 = (1.0 - k / q + k * k) * norm;
    std::vector<double> y(x.size());
    double z1 = 0.0, z2 = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        const double o = b0 * x[i] + z1;
        z1 = b1 * x[i] - a1 * o + z2;
        z2 = b2 * x[i] - a2 * o;
        y[i] = o;
    }
    return y;
}

size_t at(double s) { return (size_t)(s * deviceRate); }

double rms(const std::vector<double>& x, double t0, double t1)
{
    const auto a = at(t0), b = std::min(at(t1), x.size());
    double e = 0.0;
    for (size_t i = a; i < b; ++i)
        e += x[i] * x[i];
    return b > a ? std::sqrt(e / (double)(b - a)) : 0.0;
}

double mean(const std::vector<double>& x, double t0, double t1)
{
    const auto a = at(t0), b = std::min(at(t1), x.size());
    double e = 0.0;
    for (size_t i = a; i < b; ++i)
        e += x[i];
    return b > a ? e / (double)(b - a) : 0.0;
}

// Dominant frequency from the zero crossings of a (high-passed) signal.
double crossingHz(const std::vector<double>& x, double t0, double t1)
{
    const auto a = at(t0), b = std::min(at(t1), x.size());
    int n = 0;
    for (size_t i = a + 1; i < b; ++i)
        if ((x[i - 1] < 0.0) != (x[i] < 0.0))
            ++n;
    return n / 2.0 / ((double)(b - a) / deviceRate);
}

std::vector<float> uniformNoise(size_t n, double rmsLevel, uint32_t seed)
{
    std::vector<float> x(n);
    const double a = rmsLevel * std::sqrt(3.0);
    for (auto& v : x)
    {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        v = (float)(a * (double)(int32_t)seed / 2147483648.0);
    }
    return x;
}

std::vector<Program> programs()
{
    auto st = [](double lf, double rt, double hf) {
        rv::dsp::Settings s;
        s.lfX = lf, s.rt = rt, s.hfX = hf, s.amplPercent = 80.0;
        return s;
    };
    rv::dsp::Settings n1 = st(10.0, 1.0, 1.0), n2 = st(10.0, 1.0, 2.0); // N: EARLY 10 (flat), LATE 1.0 / 2.0 s
    n1.amplPercent = n2.amplPercent = 0.3;                               // SPAN 0.3 s
    using P = rv::params::Program;
    return {{"S 1.0_0.3_1.0", P::s, st(1.0, 0.3, 1.0)}, {"S 1.0_1.0_1.0", P::s, st(1.0, 1.0, 1.0)},
            {"S 1.0_2.5_1.0", P::s, st(1.0, 2.5, 1.0)}, {"S 2.0_0.8_0.2", P::s, st(2.0, 0.8, 0.2)},
            {"S 1.0_1.0_0.9", P::s, st(1.0, 1.0, 0.9)}, {"M 1.0_1.0_1.0", P::m, st(1.0, 1.0, 1.0)},
            {"M 1.0_5.0_0.9", P::m, st(1.0, 5.0, 0.9)}, {"L 1.0_2.0_1.0", P::l, st(1.0, 2.0, 1.0)},
            {"L 1.0_20_0.9", P::l, st(1.0, 20.0, 0.9)}, {"E 1.49_1.0_0.69", P::e, st(1.49, 1.0, 0.69)},
            {"E 1.49_5.0_0.69", P::e, st(1.49, 5.0, 0.69)}, {"N late 1.0", P::n, n1}, {"N late 2.0", P::n, n2}};
}
} // namespace

TEST_CASE("16-bit storage reproduces s_loops_fx.c bit for bit", "[sixteen]")
{
    namespace g = rv::dsp::gen::sixteen;
    for (int c = 0; c < g::numCases; ++c)
    {
        LoopsFx fx(c);
        fx.prepare();
        const auto words = fx.run();
        uint64_t h = 0xCBF29CE484222325ull;
        int64_t ss = 0;
        for (auto w : words)
        {
            h = fnv1a(h, w);
            ss += (int64_t)w * w;
        }
        std::printf("s_loops_fx case %d (%s, b %g, %s, %ld saturated in C): hash %016llx vs %016llx, sum of squares %lld vs %lld\n", c,
                    g::mode[(size_t)c] == 2 ? "rounding" : "truncation", g::inputGain[(size_t)c], g::impulse[(size_t)c] ? "impulse" : "noise", g::saturated[(size_t)c],
                    (unsigned long long)h, (unsigned long long)g::outputHash[(size_t)c], (long long)ss,
                    (long long)g::outputSumSq[(size_t)c]);
        CHECK(h == g::outputHash[(size_t)c]);
        CHECK(ss == g::outputSumSq[(size_t)c]);
    }
}

TEST_CASE("16-bit: switching back to float leaves float renders bit-identical", "[sixteen]")
{
    std::vector<float> imp(at(0.5), 0.0f);
    imp[0] = 1.0f;
    for (const auto& p : programs())
    {
        auto fresh = makeCore(p, false);
        const auto ref = render(*fresh, imp);
        auto used = makeCore(p, true);
        render(*used, imp); // run in 16-bit, then back to float
        used->setStorageFullScale(0.37);
        auto s = p.set;
        s.program = p.prog;
        s.sixteenBit = false;
        used->setSettings(s);
        used->reset();
        const auto got = render(*used, imp);
        auto q = makeCore(p, true);
        const auto q16 = render(*q, imp);
        double diff16 = 0.0;
        for (size_t i = 0; i < ref.r.size(); ++i)
            diff16 = std::max(diff16, std::abs(q16.r[i] - ref.r[i]));
        INFO(p.name);
        CHECK(got.l == ref.l);
        CHECK(got.r == ref.r);
        CHECK(diff16 > 0.0); // 16-bit mode does change the render
    }
}

TEST_CASE("16-bit: per-program error, idle tail and noise (report)", "[sixteen][report]")
{
    const double hissRms = rv::dsp::AnalogueHiss::unitRms(); // the hiss in the band, full scale 1.0
    std::printf("\nfull scale 1.0, LSB = 2^-15; nominal rms %.4f (a sine 6 dB under full scale); hiss %.1f dB re nominal\n",
                nominalRms, dbOf(nominalRms / hissRms));
    std::printf("impulse 0.5 (-6 dBFS): 16-bit minus float after a 20 Hz high-pass, dB re float's 50-300 ms rms; "
                "creep = the error's mean 1-2 s before the high-pass (LSB); tail 5-6 s: mean (LSB), rms after the "
                "high-pass (LSB and dB re nominal), its zero-crossing frequency\n");
    std::printf("idle: white input at the hiss level (H: the unit's A/D noise), output 20 Hz high-passed, 1-4 s; S/N re nominal "
                "of the core alone (float / 16-bit) and of 16-bit + hiss (spec REVERB 73 dB unweighted)\n");
    for (const auto cut : {Cut::floor, Cut::toZero, Cut::nearest})
    {
    std::printf("\n-- cut rule: %s --\n", cut == Cut::floor ? "floor (s_loops_fx.c truncation, the default)"
                                          : cut == Cut::toZero ? "toward zero (magnitude truncation)" : "nearest (s_loops_fx.c rounding)");
    for (const auto& p : programs())
    {
        const double T = 6.0;
        std::vector<float> imp(at(T), 0.0f);
        imp[0] = 0.5f;
        auto cf = makeCore(p, false);
        auto cq = makeCore(p, true, cut);
        const auto f = render(*cf, imp), q = render(*cq, imp);
        std::string line = p.name + ":";
        char buf[200];
        double worstErr = -400.0;
        for (int ch = 0; ch < 2; ++ch)
        {
            const auto& fr = ch == 0 ? f.l : f.r;
            const auto& qr = ch == 0 ? q.l : q.r;
            std::vector<double> e(fr.size());
            for (size_t i = 0; i < e.size(); ++i)
                e[i] = qr[i] - fr[i];
            const auto fh = highPass20(fr), eh = highPass20(e), qh = highPass20(qr);
            const double ref = rms(fh, 0.05, 0.3);
            const double e1 = dbOf(rms(eh, 0.05, 0.3) / ref), e2 = dbOf(rms(eh, 0.3, 0.6) / ref), e3 = dbOf(rms(eh, 0.6, 1.0) / ref),
                         e4 = dbOf(rms(eh, 1.0, 2.0) / ref);
            worstErr = std::max({worstErr, e1, e2});
            bool finite = true;
            for (auto v : qr)
                finite = finite && std::isfinite(v);
            CHECK(finite);
            std::snprintf(buf, sizeof buf, " %s err %.1f/%.1f/%.1f/%.1f creep %+.1f tail %+.1f LSB, %.2f LSB %.1f dB %.0f Hz;",
                          ch == 0 ? "L" : "R", e1, e2, e3, e4, mean(e, 1.0, 2.0) / lsb, mean(qr, 5.0, 6.0) / lsb,
                          rms(qh, 5.0, 6.0) / lsb, dbOf(rms(qh, 5.0, 6.0) / nominalRms), crossingHz(qh, 5.0, 6.0));
            line += buf;
        }
        std::printf("%s\n", line.c_str());

        // idle noise
        const auto noise = uniformNoise(at(4.0), hissRms, 0x1234567u);
        auto nf = makeCore(p, false);
        auto nq = makeCore(p, true, cut);
        const auto yf = render(*nf, noise), yq = render(*nq, noise);
        double pf = 0.0, pq = 0.0;
        for (int ch = 0; ch < 2; ++ch)
        {
            const auto a = highPass20(ch == 0 ? yf.l : yf.r), b = highPass20(ch == 0 ? yq.l : yq.r);
            pf += 0.5 * std::pow(rms(a, 1.0, 4.0), 2.0);
            pq += 0.5 * std::pow(rms(b, 1.0, 4.0), 2.0);
        }
        std::printf("    idle S/N re nominal: float core %.1f dB, 16-bit core %.1f dB, 16-bit + hiss %.1f dB (float + hiss %.1f)\n",
                    dbOf(nominalRms / std::sqrt(pf)), dbOf(nominalRms / std::sqrt(pq)),
                    dbOf(nominalRms / std::sqrt(pq + hissRms * hissRms)), dbOf(nominalRms / std::sqrt(pf + hissRms * hissRms)));
    }
    }
}

TEST_CASE("16-bit: full-scale noise and impulses stay bounded", "[sixteen]")
{
    for (const auto& p : programs())
    {
        auto c = makeCore(p, true);
        auto x = uniformNoise(at(4.0), 1.0 / std::sqrt(3.0), 0xBEEFu); // uniform in [-1, 1): the A/D's whole range
        std::fill(x.begin() + (long)at(2.0), x.end(), 0.0f);
        x[at(2.5)] = 1.0f;
        const auto y = render(*c, x);
        double peak = 0.0;
        bool finite = true;
        for (size_t i = 0; i < y.l.size(); ++i)
        {
            finite = finite && std::isfinite(y.l[i]) && std::isfinite(y.r[i]);
            peak = std::max({peak, std::abs(y.l[i]), std::abs(y.r[i])});
        }
        std::printf("%s: full-scale noise 2 s + impulse, 16-bit: peak output %.2f FS, last 0.5 s rms %.2g LSB\n", p.name.c_str(),
                    peak, rms(highPass20(y.r), 3.5, 4.0) / lsb);
        INFO(p.name);
        CHECK(finite);
        CHECK(peak < 100.0);

        // switching the mode while running (no reset): the filter state carries over, nothing jumps or runs away
        auto t = makeCore(p, false);
        const auto x1 = uniformNoise(at(1.0), 0.1, 0xF00Du);
        auto set = p.set;
        set.program = p.prog;
        double before = 0.0, during = 0.0;
        for (int k = 0; k < 4; ++k)
        {
            set.sixteenBit = (k % 2) == 1;
            t->setSettings(set);
            const auto yk = render(*t, x1);
            for (size_t i = 0; i < yk.l.size(); ++i)
            {
                const double m = std::max(std::abs(yk.l[i]), std::abs(yk.r[i]));
                (k == 0 ? before : during) = std::max(k == 0 ? before : during, std::isfinite(m) ? m : 1e30);
            }
        }
        std::printf("    mode switched every second (float, 16, float, 16): peak %.3f against %.3f in float alone\n", during, before);
        CHECK(during < 2.0 * before);
    }
}

TEST_CASE("16-bit: cost per device sample (report)", "[sixteen]")
{
    const auto x = uniformNoise(at(10.0), 0.1, 0x5EEDu);
    for (const auto& p : programs())
    {
        double ns[2] = {0.0, 0.0};
        for (int m = 0; m < 2; ++m)
        {
            auto c = makeCore(p, m == 1);
            const auto t0 = std::chrono::steady_clock::now();
            render(*c, x);
            const auto t1 = std::chrono::steady_clock::now();
            ns[m] = std::chrono::duration<double, std::nano>(t1 - t0).count() / (double)x.size();
        }
        std::printf("%s: %.1f ns per device sample float, %.1f 16-bit (render + copy)\n", p.name.c_str(), ns[0], ns[1]);
    }
}

TEST_CASE("Analogue hiss: level, L/R independence, determinism", "[sixteen][hiss]")
{
    for (double rate : {44100.0, 48000.0, 96000.0})
    {
        rv::dsp::AnalogueHiss h;
        h.prepare(rate);
        h.setFullScale(2.0);
        const int n = (int)(10.0 * rate);
        std::vector<float> l((size_t)n, 0.0f), r((size_t)n, 0.0f);
        for (int i = 0; i < n; i += 512)
            h.add(l.data() + i, r.data() + i, std::min(512, n - i), true);
        double el = 0.0, er = 0.0, lr = 0.0;
        for (int i = 0; i < n; ++i)
        {
            el += (double)l[(size_t)i] * l[(size_t)i];
            er += (double)r[(size_t)i] * r[(size_t)i];
            lr += (double)l[(size_t)i] * r[(size_t)i];
        }
        const double rl = std::sqrt(el / n), rr = std::sqrt(er / n), rho = lr / std::sqrt(el * er);
        const double target = h.targetRms();
        std::printf("hiss at %g Hz, fs 2.0: rms L %.4g R %.4g (target in 20 Hz-20 kHz %.4g: %+.2f / %+.2f dB; %.1f dB re a full-scale sine), "
                    "L/R correlation %+.4f\n",
                    rate, rl, rr, target, dbOf(rl / target), dbOf(rr / target), dbOf(rr / (2.0 / std::sqrt(2.0))), rho);
        // the whole band holds a little more than 20 Hz-20 kHz (the low-pass's skirt above 20 kHz)
        CHECK(dbOf(rl / target) > -0.1);
        CHECK(dbOf(rl / target) < 0.6);
        CHECK(dbOf(rr / target) > -0.1);
        CHECK(dbOf(rr / target) < 0.6);
        CHECK(std::abs(rho) < 0.01);

        // deterministic after reset(); `on = false` adds nothing
        h.reset();
        std::vector<float> l2(512, 0.0f), r2(512, 0.0f);
        h.add(l2.data(), r2.data(), 512, true);
        CHECK(std::equal(l2.begin(), l2.end(), l.begin()));
        std::vector<float> z(512, 0.0f), z2(512, 0.0f);
        h.add(z.data(), z2.data(), 512, false);
        CHECK(std::all_of(z.begin(), z.end(), [](float v) { return v == 0.0f; }));
    }
}
