#include "dsp/AdcCalibration.h"
#include "dsp/ReverbEngine.h"
#include "params/Parameters.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <random>
#include <vector>

// ReverbEngine end to end at host rates (docs/ENGINE.md, Tests): every program through the bridge, the common time
// grid and latency, the program crossfade, same-program changes, DELAY, LOW PASS, block-size independence, no
// allocation in process(), and the cost per program.

//==============================================================================
// Allocation counter: this runner replaces the global operator new / delete (standard C++; nothing else in the DSP
// runner does), and a test counts the calls made while `countAllocations` is set.
namespace
{
std::atomic<bool> countAllocations{false};
std::atomic<long> allocations{0};

void* allocate(std::size_t n)
{
    if (countAllocations.load(std::memory_order_relaxed))
        allocations.fetch_add(1, std::memory_order_relaxed);
    if (auto* p = std::malloc(n == 0 ? 1 : n))
        return p;
    throw std::bad_alloc();
}
} // namespace

void* operator new(std::size_t n)
{
    return allocate(n);
}
void* operator new[](std::size_t n)
{
    return allocate(n);
}
void operator delete(void* p) noexcept
{
    std::free(p);
}
void operator delete[](void* p) noexcept
{
    std::free(p);
}
void operator delete(void* p, std::size_t) noexcept
{
    std::free(p);
}
void operator delete[](void* p, std::size_t) noexcept
{
    std::free(p);
}

//==============================================================================
namespace
{
using rv::dsp::ReverbEngine;
using rv::dsp::Settings;
using rv::params::Joy;
using rv::params::Program;

constexpr std::array<double, 3> rates{44100.0, 48000.0, 96000.0};
constexpr std::array<Program, 5> programs{Program::s, Program::m, Program::l, Program::e, Program::n};
constexpr double twoPi = 6.283185307179586;

const char* name(Program p)
{
    static const char* names[] = {"S", "M", "L", "E", "N"};
    return names[(int)p];
}

// The program's default panel (as the processor builds it), with DELAY / AMPL / LP overridable (in their units).
Settings settingsFor(Program p, double delayMs = 0.0, double amplStep = 10, double lpKHz = 12.0, bool sixteen = false)
{
    using namespace rv::params;
    Settings s;
    s.program = p;
    for (int j = 0; j < numJoys; ++j)
        s.steps[(size_t)j] = std::min(joyDefaultStep((Joy)j), joyMaxStep(p, (Joy)j));
    s.steps[(size_t)Joy::ampl] = (int)amplStep;
    s.steps[(size_t)Joy::delay] = (int)delayMs;
    s.amplPercent = joyNumber(p, Joy::ampl, s.steps[(size_t)Joy::ampl]);
    s.lfX = joyNumber(p, Joy::lf, s.steps[(size_t)Joy::lf]);
    s.rt = joyNumber(p, Joy::rt, s.steps[(size_t)Joy::rt]);
    s.hfX = joyNumber(p, Joy::hf, s.steps[(size_t)Joy::hf]);
    s.delayMs = delayMs;
    s.lowpassKHz = lpKHz;
    s.sixteenBit = sixteen;
    return s;
}

void start(ReverbEngine& e, double rate, int block, const Settings& s)
{
    e.prepare(rate, block);
    e.setFullScale(rv::dsp::adcFullScaleForZeroVu(-18.0f));
    e.setSettings(s);
    e.reset(); // start in s.program rather than crossfading into it
}

struct Stereo
{
    std::vector<float> l, r;
};

Stereo render(ReverbEngine& e, const std::vector<float>& x, int block)
{
    Stereo o{std::vector<float>(x.size()), std::vector<float>(x.size())};
    for (size_t i = 0; i < x.size(); i += (size_t)block)
    {
        const auto n = (int)std::min((size_t)block, x.size() - i);
        e.process(x.data() + i, o.l.data() + i, o.r.data() + i, n);
    }
    return o;
}

std::vector<float> impulse(size_t n, float amp = 0.1f)
{
    std::vector<float> x(n, 0.0f);
    x[0] = amp;
    return x;
}

std::vector<float> noise(size_t n, float amp, unsigned seed)
{
    std::mt19937 g(seed);
    std::uniform_real_distribution<float> u(-amp, amp);
    std::vector<float> x(n);
    for (auto& v : x)
        v = u(g);
    return x;
}

std::vector<float> sine(double rate, double hz, float amp, size_t n)
{
    std::vector<float> x(n);
    for (size_t i = 0; i < n; ++i)
        x[i] = amp * (float)std::sin(twoPi * hz * (double)i / rate);
    return x;
}

// A sine faded in over 50 ms (raised cosine), so its onset excites no broadband reverb to mask a later click.
std::vector<float> softSine(double rate, double hz, float amp, size_t n)
{
    auto x = sine(rate, hz, amp, n);
    const auto ramp = (size_t)(0.05 * rate);
    for (size_t i = 0; i < std::min(ramp, n); ++i)
        x[i] *= (float)(0.5 - 0.5 * std::cos(0.5 * twoPi * (double)i / (double)ramp));
    return x;
}

double energy(const std::vector<float>& y, size_t a, size_t b)
{
    double e = 0.0;
    for (size_t i = a; i < std::min(b, y.size()); ++i)
        e += (double)y[i] * y[i];
    return e;
}

double db(double ratio)
{
    return 10.0 * std::log10(std::max(ratio, 1e-300));
}

bool allFinite(const Stereo& y)
{
    for (size_t i = 0; i < y.l.size(); ++i)
        if (!std::isfinite(y.l[i]) || !std::isfinite(y.r[i]))
            return false;
    return true;
}

// Index of the largest |y| in [a, b).
size_t peakIn(const std::vector<float>& y, size_t a, size_t b)
{
    size_t k = a;
    for (size_t i = a; i < std::min(b, y.size()); ++i)
        if (std::abs(y[i]) > std::abs(y[k]))
            k = i;
    return k;
}

// Amplitude of the component at `hz` over [a, b) (Hann window).
double amplitude(const std::vector<float>& y, double rate, double hz, size_t a, size_t b)
{
    std::complex<double> sum;
    double wsum = 0.0;
    const auto n = b - a;
    for (size_t i = 0; i < n; ++i)
    {
        const auto w = 0.5 - 0.5 * std::cos(twoPi * (double)i / (double)n);
        sum += w * (double)y[a + i] * std::polar(1.0, -twoPi * hz * (double)(a + i) / rate);
        wsum += w;
    }
    return 2.0 * std::abs(sum) / wsum;
}

// Click detector: each channel through a 2-8 kHz band-pass (two 2nd-order high-passes at 2 kHz, two low-passes at
// 8 kHz, RBJ), which takes a 200 Hz sine down 80 dB and the analogue chain's images of it (21-22 kHz, about -57 dB)
// down 35 dB, but passes most of a click; returns the band's largest |sample| in [a, b), both channels.
struct Biquad
{
    double b0, b1, b2, a1, a2, z1 = 0.0, z2 = 0.0;
    static Biquad make(double rate, double hz, bool highPass)
    {
        const double w = twoPi * hz / rate, c = std::cos(w), al = std::sin(w) / (2.0 * std::sqrt(0.5)), a0 = 1.0 + al;
        const double k = highPass ? (1.0 + c) / 2.0 : (1.0 - c) / 2.0;
        return {k / a0, (highPass ? -2.0 : 2.0) * k / a0, k / a0, -2.0 * c / a0, (1.0 - al) / a0};
    }
    double run(double x)
    {
        const double y = b0 * x + z1;
        z1 = b1 * x - a1 * y + z2;
        z2 = b2 * x - a2 * y;
        return y;
    }
};

float clickBand(const Stereo& y, double rate, size_t a, size_t b)
{
    float m = 0.0f;
    for (const auto* ch : {&y.l, &y.r})
    {
        std::array<Biquad, 4> f{Biquad::make(rate, 2000.0, true), Biquad::make(rate, 2000.0, true),
                                Biquad::make(rate, 8000.0, false), Biquad::make(rate, 8000.0, false)};
        for (size_t i = 0; i < std::min(b, ch->size()); ++i)
        {
            double v = (*ch)[i];
            for (auto& q : f)
                v = q.run(v);
            if (i >= a)
                m = std::max(m, (float)std::abs(v));
        }
    }
    return m;
}
} // namespace

//==============================================================================
TEST_CASE("every program's impulse response is finite at 44.1 / 48 / 96 kHz, with one latency for all")
{
    for (const auto rate : rates)
    {
        int latency = -1;
        for (const auto p : programs)
        {
            ReverbEngine e;
            start(e, rate, 512, settingsFor(p));
            if (latency < 0)
            {
                latency = e.latency();
                std::printf("engine latency at %g Hz: %d host samples (%.3f ms; alignment %d device samples)\n", rate,
                            latency, 1000.0 * latency / rate, e.alignmentDelay());
            }
            INFO(name(p) << " at " << rate);
            CHECK(e.latency() == latency);
            CHECK(e.alignmentDelay() == 54);
            const auto y = render(e, impulse((size_t)(1.5 * rate)), 512);
            CHECK(allFinite(y));
            CHECK(energy(y.l, 0, y.l.size()) > 1e-6);
            CHECK(energy(y.r, 0, y.r.size()) > 1e-6);
            // nothing before the latency but the analogue chain's own pre-ringing (its sinc), far below the signal
            const auto pre = energy(y.r, 0, (size_t)latency - 8) + energy(y.l, 0, (size_t)latency - 8);
            CHECK(db(pre / (energy(y.r, 0, y.r.size()) + energy(y.l, 0, y.l.size()))) < -60.0);
        }
    }
}

TEST_CASE("the programs share one time grid: E's direct tap arrives with M's first cluster tap")
{
    // Device rate first: M's R cluster tap 0 comes out at its kernel centre, E's R direct tap at k0.
    {
        rv::dsp::MLCore m(Program::m);
        rv::dsp::ECore ec;
        for (rv::dsp::ProgramCore* c : {(rv::dsp::ProgramCore*)&m, (rv::dsp::ProgramCore*)&ec})
        {
            c->prepare();
            c->setSettings(settingsFor(c == &m ? Program::m : Program::e));
            c->reset();
        }
        std::vector<float> in(200, 0.0f), l(200), r(200);
        in[0] = 1.0f;
        m.process(in.data(), l.data(), r.data(), 200);
        CHECK(peakIn(r, 0, 200) == 10);
        ec.process(in.data(), l.data(), r.data(), 200);
        CHECK(peakIn(r, 0, 200) == 64);
    }
    // Host rate: the first arrivals (cluster / direct tap 0 at AMPL 100, DELAY 0) of M and E coincide; S's is 10 device
    // samples earlier (its cluster has no kernel delay) and L's 20 later (its kernel's centre is 30); and S's lands at
    // the latency plus what the analogue chain itself takes to peak (a bridge with an identity program).
    for (const auto rate : rates)
    {
        const auto r = rate / rv::dsp::deviceRate;
        std::array<double, 5> at{};
        int latency = 0;
        for (const auto p : {Program::s, Program::m, Program::l, Program::e})
        {
            ReverbEngine e;
            start(e, rate, 512, settingsFor(p));
            latency = e.latency();
            const auto y = render(e, impulse((size_t)(0.1 * rate)), 512);
            at[(size_t)p] = (double)peakIn(y.r, 0, (size_t)latency + (size_t)(0.003 * rate));
        }
        rv::dsp::DeviceBridge b;
        b.prepare(rate, 512);
        std::vector<float> x = impulse((size_t)(0.05 * rate)), l(x.size()), rr(x.size());
        b.process(x.data(), l.data(), rr.data(), (int)x.size(), [](const float* in, float* dl, float* dr, int n) {
            std::copy(in, in + n, dl);
            std::copy(in, in + n, dr);
        });
        const auto chainPeak = (double)peakIn(rr, 0, rr.size()) + 54.0 * r; // the engine's alignment, unrounded
        INFO("rate " << rate << ": S " << at[0] << " M " << at[1] << " L " << at[2] << " E " << at[3] << ", latency "
                     << latency << ", chain peak " << chainPeak);
        CHECK(std::abs(at[(size_t)Program::e] - at[(size_t)Program::m]) <= 1.0);
        CHECK(std::abs(at[(size_t)Program::m] - at[(size_t)Program::s] - 10.0 * r) <= 1.0);
        CHECK(std::abs(at[(size_t)Program::l] - at[(size_t)Program::s] - 30.0 * r) <= 1.0);
        CHECK(std::abs(at[(size_t)Program::s] - chainPeak) <= 1.0);
    }
}

TEST_CASE("the output is the same for any block size")
{
    for (const auto p : programs)
    {
        INFO(name(p));
        auto x = noise(48000, 0.1f, 3);
        std::fill(x.begin() + 4800, x.end(), 0.0f);
        Stereo ref;
        for (const int block : {512, 1, 37, 300})
        {
            ReverbEngine e;
            start(e, 48000.0, 512, settingsFor(p, 10.0, 8, 6.0));
            const auto y = render(e, x, block);
            if (block == 512)
                ref = y;
            else
            {
                bool same = true;
                for (size_t i = 0; i < x.size() && same; ++i)
                    same = std::equal_to<float>{}(y.l[i], ref.l[i]) && std::equal_to<float>{}(y.r[i], ref.r[i]);
                CHECK(same);
            }
        }
    }
}

namespace
{
struct SwitchResult
{
    float before, across, after;
};

// A soft 200 Hz sine through `from`, switched to `to` (and then to `then`, queued, if given) 4 s in; the 2-8 kHz
// click band before, across and well after the switch.
SwitchResult runSwitch(Program from, Program to, Program then = Program::count)
{
    const double rate = 48000.0;
    ReverbEngine e;
    start(e, rate, 256, settingsFor(from, 10.0, 8));
    const auto x = softSine(rate, 200.0, 0.1f, (size_t)(6.0 * rate));
    Stereo y{std::vector<float>(x.size()), std::vector<float>(x.size())};
    const size_t switchAt = 750 * 256; // 4 s in (the onset's tail is long gone at RT 1 s), on a block boundary
    for (size_t i = 0; i < x.size(); i += 256)
    {
        if (i == switchAt)
        {
            e.setSettings(settingsFor(to, 10.0, 8));
            CHECK(e.isCrossfading());
            if (then != Program::count)
                e.setSettings(settingsFor(then, 10.0, 8)); // queued behind the first
        }
        e.process(x.data() + i, y.l.data() + i, y.r.data() + i, 256);
        if (then != Program::count && i == switchAt + 10 * 256) // 53 ms: the first fade is over, the queued one on
            CHECK(e.activeProgram() == then);
    }
    CHECK(allFinite(y));
    CHECK(!e.isCrossfading());
    CHECK(e.activeProgram() == (then != Program::count ? then : to));
    return {clickBand(y, rate, switchAt - (size_t)(0.5 * rate), switchAt),
            clickBand(y, rate, switchAt, switchAt + (size_t)(0.2 * rate)),
            clickBand(y, rate, x.size() - (size_t)(0.5 * rate), x.size())};
}
} // namespace

TEST_CASE("a program switch crossfades without a click, and a second switch waits its turn")
{
    struct Case
    {
        Program a, b, c;
    };
    for (const auto& k : {Case{Program::s, Program::m, Program::count}, Case{Program::m, Program::s, Program::count},
                          Case{Program::s, Program::e, Program::count}, Case{Program::e, Program::n, Program::count},
                          Case{Program::n, Program::l, Program::count}, Case{Program::s, Program::m, Program::l}})
    {
        const auto r = runSwitch(k.a, k.b, k.c);
        std::printf("program switch %s -> %s%s%s (200 Hz sine): 2-8 kHz peak %.3g across, %.3g before, %.3g after\n",
                    name(k.a), name(k.b), k.c != Program::count ? " -> " : "", k.c != Program::count ? name(k.c) : "",
                    r.across, r.before, r.after);
        CHECK(r.across < 3.0f * std::max(r.before, r.after));
    }
}

TEST_CASE("same-program setting changes apply at once, with no blow-up")
{
    // Every step of LF, RT, HF and AMPL in turn, one per 256-sample block, with noise running: the output stays
    // finite and bounded (no state is lost, so nothing rings up), and each core sees the change at once.
    const double rate = 48000.0;
    for (const auto p : programs)
    {
        ReverbEngine e;
        start(e, rate, 256, settingsFor(p, 10.0, 8));
        const auto x = noise((size_t)(rate * 4.0), 0.05f, 11);
        Stereo y{std::vector<float>(x.size()), std::vector<float>(x.size())};
        size_t block = 0;
        float steadyPeak = 0.0f;
        for (size_t i = 0; i < x.size(); i += 256, ++block)
        {
            if (i >= (size_t)rate) // 1 s of the default first
            {
                auto s = settingsFor(p, 10.0, 8);
                const auto k = (int)(block % 64);
                const auto j = (Joy)((int)Joy::ampl + (k / 16) % 4);
                s.steps[(size_t)j] = std::min(k % 16, rv::params::joyMaxStep(p, j));
                s.amplPercent = rv::params::joyNumber(p, Joy::ampl, s.steps[(size_t)Joy::ampl]);
                s.lfX = rv::params::joyNumber(p, Joy::lf, s.steps[(size_t)Joy::lf]);
                s.rt = rv::params::joyNumber(p, Joy::rt, s.steps[(size_t)Joy::rt]);
                s.hfX = rv::params::joyNumber(p, Joy::hf, s.steps[(size_t)Joy::hf]);
                e.setSettings(s);
                CHECK_FALSE(e.isCrossfading());
            }
            e.process(x.data() + i, y.l.data() + i, y.r.data() + i, 256);
            if (i < (size_t)rate)
                for (size_t k = i; k < i + 256; ++k)
                    steadyPeak = std::max({steadyPeak, std::abs(y.l[k]), std::abs(y.r[k])});
        }
        float changingPeak = 0.0f;
        for (size_t k = (size_t)rate; k < x.size(); ++k)
            changingPeak = std::max({changingPeak, std::abs(y.l[k]), std::abs(y.r[k])});
        std::printf("%s: settings stepping every 5.3 ms: peak %.3g against %.3g at the default (x%.2f)\n", name(p),
                    changingPeak, steadyPeak, changingPeak / steadyPeak);
        INFO(name(p));
        CHECK(allFinite(y));
        CHECK(changingPeak < 8.0f * steadyPeak);
    }
}

TEST_CASE("DELAY is a pure delay of the whole output in device samples, and a change crossfades")
{
    const double rate = 48000.0;
    const auto x = impulse((size_t)(1.6 * rate));
    Stereo y0, y1;
    {
        ReverbEngine e;
        start(e, rate, 512, settingsFor(Program::s, 0.0));
        y0 = render(e, x, 512);
    }
    {
        ReverbEngine e;
        start(e, rate, 512, settingsFor(Program::s, 100.0));
        y1 = render(e, x, 512);
    }
    const auto shift = (double)std::lround(0.1 * rv::dsp::deviceRate) * rate / rv::dsp::deviceRate; // 2637 device samples
    // tap 0 of the cluster (S's next R tap is 143 device samples, 5.4 ms, later)
    const auto window = (size_t)(0.003 * rate);
    const auto p0 = (double)peakIn(y0.r, 0, (size_t)(0.003 * rate) + 200);
    const auto p1 = (double)peakIn(y1.r, (size_t)(p0 + shift) - window, (size_t)(p0 + shift) + window);
    CHECK(std::abs(p1 - p0 - shift) <= 1.0);
    // S is time-invariant, so the shifted response has the same energy (the bridge's sub-sample grid only).
    const auto n = (size_t)(1.0 * rate);
    CHECK(std::abs(db(energy(y1.r, (size_t)p1 - 50, (size_t)p1 - 50 + n) / energy(y0.r, (size_t)p0 - 50, (size_t)p0 - 50 + n))) < 0.01);
    CHECK(std::equal_to<float>{}(y1.r[(size_t)p1 - 200], 0.0f)); // and nothing ahead of it but the bridge's own pre-ringing tail
    // A held lever: DELAY stepping 1 ms per block under a sine stays smooth.
    ReverbEngine e;
    start(e, rate, 256, settingsFor(Program::s, 20.0, 8));
    const auto s = softSine(rate, 200.0, 0.1f, (size_t)(5.0 * rate));
    Stereo z{std::vector<float>(s.size()), std::vector<float>(s.size())};
    for (size_t i = 0, b = 0; i < s.size(); i += 256, ++b)
    {
        if (i >= (size_t)(4.0 * rate))
            e.setSettings(settingsFor(Program::s, 20.0 + (double)(b % 40), 8));
        e.process(s.data() + i, z.l.data() + i, z.r.data() + i, 256);
    }
    const auto steady = clickBand(z, rate, (size_t)(3.5 * rate), (size_t)(4.0 * rate));
    const auto moving = clickBand(z, rate, (size_t)(4.0 * rate), s.size());
    std::printf("DELAY stepping 1 ms per block (200 Hz sine): 2-8 kHz peak %.3g, steady %.3g (x%.2f)\n", moving, steady,
                moving / steady);
    CHECK(moving < 3.0f * steady);
}

TEST_CASE("LOW PASS: 12 kHz is an exact bypass, below it a 2nd-order Butterworth at the device rate")
{
    // The filter on its own.
    rv::dsp::DeviceLowPass lp;
    lp.setCutoff(12.0);
    lp.reset();
    auto a = noise(4096, 0.5f, 5), b = noise(4096, 0.5f, 6);
    const auto a0 = a, b0 = b;
    lp.process(a.data(), b.data(), (int)a.size());
    CHECK(a == a0);
    CHECK(b == b0);
    for (const double fc : {2.0, 5.0, 10.0})
    {
        rv::dsp::DeviceLowPass f;
        f.setCutoff(fc);
        f.reset();
        for (const double hz : {100.0, fc * 1000.0, 2.0 * fc * 1000.0})
        {
            if (hz >= 0.5 * rv::dsp::deviceRate)
                continue;
            auto l = sine(rv::dsp::deviceRate, hz, 1.0f, 26374), r = l;
            f.reset();
            f.process(l.data(), r.data(), (int)l.size());
            const auto got = 20.0 * std::log10(amplitude(l, rv::dsp::deviceRate, hz, 13187, 26374));
            // bilinear: the Butterworth magnitude at the warped frequency
            const auto w = std::tan(0.5 * twoPi * hz / rv::dsp::deviceRate) / std::tan(0.5 * twoPi * fc * 1000.0 / rv::dsp::deviceRate);
            const auto want = -10.0 * std::log10(1.0 + std::pow(w, 4.0));
            INFO("fc " << fc << " kHz at " << hz << " Hz: " << got << " dB, want " << want);
            CHECK(std::abs(got - want) < 0.05);
        }
    }
    // Through the engine (S, time-invariant): an 8 kHz sine loses what the filter takes at 8 kHz, a 300 Hz one nothing.
    const double rate = 48000.0;
    for (const double hz : {300.0, 8000.0})
    {
        std::array<double, 2> amp{};
        for (int k = 0; k < 2; ++k)
        {
            ReverbEngine e;
            start(e, rate, 512, settingsFor(Program::s, 0.0, 8, k == 0 ? 12.0 : 2.0));
            const auto y = render(e, sine(rate, hz, 0.05f, (size_t)(3.0 * rate)), 512);
            amp[(size_t)k] = amplitude(y.r, rate, hz, (size_t)(1.5 * rate), (size_t)(3.0 * rate));
        }
        const auto w = std::tan(0.5 * twoPi * hz / rv::dsp::deviceRate) / std::tan(0.5 * twoPi * 2000.0 / rv::dsp::deviceRate);
        const auto want = -10.0 * std::log10(1.0 + std::pow(w, 4.0));
        const auto got = 20.0 * std::log10(amp[1] / amp[0]);
        INFO(hz << " Hz: LP 2 kHz takes " << got << " dB, want " << want);
        CHECK(std::abs(got - want) < 0.5);
    }
}

TEST_CASE("the mono fold-down (L + R) / sqrt(2) keeps each side's level")
{
    // L and R are uncorrelated read-outs of one network, so their equal-power sum is as loud as either side.
    const double rate = 48000.0;
    for (const auto p : programs)
    {
        ReverbEngine e;
        start(e, rate, 512, settingsFor(p, 0.0, 8));
        const auto y = render(e, noise((size_t)(3.0 * rate), 0.1f, 9), 512);
        const auto a = (size_t)(1.0 * rate), b = y.l.size();
        const auto el = energy(y.l, a, b), er = energy(y.r, a, b);
        double em = 0.0;
        for (size_t i = a; i < b; ++i)
        {
            const auto m = (y.l[i] + y.r[i]) / std::sqrt(2.0);
            em += m * m;
        }
        const auto rel = db(em / (0.5 * (el + er)));
        std::printf("%s: L / R %.2f dB, fold-down re their mean %+.2f dB\n", name(p), db(el / er), rel);
        INFO(name(p));
        CHECK(std::abs(rel) < 1.0);
    }
}

TEST_CASE("16-bit mode runs clean and flags the A/D")
{
    const double rate = 48000.0;
    for (const auto p : programs)
    {
        ReverbEngine e;
        start(e, rate, 512, settingsFor(p, 10.0, 8, 12.0, true));
        const auto y = render(e, sine(rate, 1000.0, 0.2f, (size_t)rate), 512);
        INFO(name(p));
        CHECK(allFinite(y));
        CHECK_FALSE(e.takeOverload()); // 0.2 peak is under the default full scale (0.355)
        const auto z = render(e, sine(rate, 1000.0, 0.8f, (size_t)(0.2 * rate)), 512);
        CHECK(allFinite(z));
        CHECK(e.takeOverload());
        CHECK_FALSE(e.takeOverload());
    }
}

TEST_CASE("process() does not allocate, in any program or across a switch")
{
    const double rate = 48000.0;
    ReverbEngine e;
    start(e, rate, 512, settingsFor(Program::s, 10.0, 8, 6.0, true));
    std::vector<float> x = noise(512, 0.1f, 1), l(512), r(512);
    allocations = 0;
    countAllocations = true;
    for (int round = 0; round < 3; ++round)
        for (const auto p : programs)
        {
            e.setSettings(settingsFor(p, 10.0 + round, 8, round == 1 ? 12.0 : 4.0, round != 2));
            for (int b = 0; b < 8; ++b)
                e.process(x.data(), l.data(), r.data(), 512);
        }
    e.reset();
    e.process(x.data(), l.data(), r.data(), 512);
    countAllocations = false;
    CHECK(allocations.load() == 0);
}

TEST_CASE("cost per engine and program")
{
    for (const auto rate : {48000.0, 96000.0})
        for (const auto p : programs)
        {
            ReverbEngine e;
            start(e, rate, 512, settingsFor(p, 10.0, 8, 12.0, true));
            const auto seconds = 5.0;
            const auto x = noise((size_t)(rate * seconds), 0.1f, 2);
            std::vector<float> l(512), r(512);
            const auto t0 = std::chrono::steady_clock::now();
            for (size_t i = 0; i + 512 <= x.size(); i += 512)
                e.process(x.data() + i, l.data(), r.data(), 512);
            const auto t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::printf("cost %s at %g Hz (16-bit, block 512): %.2f %% of one core\n", name(p), rate, 100.0 * t / seconds);
            CHECK(t < seconds); // faster than real time, whatever the machine
        }
}
