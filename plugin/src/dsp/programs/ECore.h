#pragma once

// Program E (the EMT 250 algorithm as the 246 runs it) at the device rate, from analysis/e250/e_plugin_spec.json
// (docs/EMT250_program.md 11-17) and the sample-by-sample reference tools/e250/e_sim.py render(). Mono in, (L, R) out:
//
//   direct     taps 0 / 429 / 718 (+1 / -0.709 / +0.589) x AMPL / 80; L: the same 314 later, inverted, x gL
//   all-passes 46 / 146 / 392 in series, (-g + z^-D) / (1 - g z^-D), before the loops
//   loops      three, P = 2215 / 2507 / 2776, negative feedback: v[t] = u[t] - line[t - P], line = G{v}, G one pass
//              of the 246's loop law at RT / LF / HF (ml_lines.line_coefs: shelf + HF pole). The delay lines hold v;
//              G runs on each read of the line (the feedback at P, the static read at d), which is the same
//              signal as filtering once and storing `line` (G is linear and time invariant)
//   static     one read per loop (2223 / 805 / 1802) of line, signed gains; L the sum 315 later, x gL
//   moving     nine reads of v (before G, 17.2) at centre + 32 tri(f_j t + phase_j), tri(x) = 4 |frac(x) - 0.5| - 1, rounded to an integer
//              address (no interpolation, as e_sim.e_reads); L reads the same address with its own final amplitude
//              (gL included), one sample further for the reads whose signs agree, and the L sum arrives 314 later;
//              a one-pole 0.07 on both sums (17.7). setModulation(false) removes them: the static network, which is
//              what the reference IRs are (e_sim.render_ir moving=False)
//   kernel     [0.77, 0.23] on both outputs (14.1)
//
// R's direct tap comes out at gen::e::k0 = 64 samples, as in the reference IRs (16 of them are the output kernel's
// offset); latency() is 0. The panel's LF x / HF x go straight into the loop law (the Altiverb family is LF 1.49 /
// HF 0.69); RT is clamped to 0.4-5 s and the law is evaluated at any value in it.
//
// Precision: the recirculating storage (the all-passes' and the loops' delay lines, written through
// ProgramCore::store) and the input line are float; every coefficient, filter state and sum is double. The reference
// IRs are float32 renders of a double model; tests/dsp/ETests.cpp reports where float storage puts the null.
// Cost per sample: 6 direct taps, 3 all-passes (2 MACs), 3 loops x 2 runs of G (4 MACs) + the static read, 9 moving
// reads (2 MACs, a phase step and the triangle, one rounding), the one-poles and the kernel: about 90
// multiply-adds; measured 23 ns per device sample (0.06 % of one core) on an Apple-silicon Mac, Release.

#include "dsp/engine/ProgramCore.h"
#include "dsp/generated/e_spec.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <juce_core/juce_core.h>

namespace rv::dsp
{
class ECore final : public ProgramCore
{
  public:
    static constexpr int numLoops = (int)gen::e::loopLength.size();
    static constexpr int numReads = gen::e::numReads;
    static constexpr double rtMin = 0.4, rtMax = 5.0;

    // One pass of a loop's filter as the law gives it (ml_lines.line_coefs): gain g, LF shelf r0, HF pole p.
    struct LoopLaw
    {
        double g = 1.0, r0 = 1.0, p = 0.0;
    };

    // The 246's loop law for loop k at (RT, LF, HF) (ml_lines.line_coefs with its caps; model_s_results 14.2,
    // 16.15-16.16, 29, 32.6). Any values; the LF law is interpolated in log LF between S's steps and held at the ends.
    static LoopLaw loopLaw(int k, double rt, double lf, double hf)
    {
        namespace g = gen::e;
        double hi = 0.0, dc = 0.0, kl = 0.0;
        lfLaw(lf, hi, dc, kl);
        hf = std::max(hf, g::hfFloor[(size_t)k] / rt); // the HF control's floor (model_s_results 29)
        const double rtHi = hi * rt, rtDc = dc * rt;
        const double kk = -3.0 * g::loopLength[(size_t)k] / deviceRate;
        LoopLaw c;
        c.g = std::pow(10.0, kk / rtHi);
        c.r0 = std::pow(10.0, kk * (1.0 / rtDc - 1.0 / rtHi));
        if (c.g * c.r0 < g::dcMin) // the per-pass DC gain limit (model_s_results 16.16)
        {
            const double lift = 20.0 * std::log10(g::dcMin / (c.g * c.r0));
            c.g *= std::pow(10.0, -g::beta * lift / 20.0);
            c.r0 = g::dcMin / c.g;
        }
        const double t = hf < 1.0 ? std::pow(10.0, kk / rt * kl * (1.0 / hf - 1.0)) : 1.0;
        double p = t < 1.0 ? onePoleFor(t, g::lawFRef) : 0.0;
        if (p > 0.0) // the HF pole curve with its drop-back past the fold (model_s_results 32.6)
            p = std::max(p < g::poleFold ? interp(p, g::poleKnots, g::poleCurve) : g::poleDropBack, 0.0);
        c.p = p;
        return c;
    }

    void prepare() override
    {
        inLine.assign(pow2Above(gen::e::k0 + gen::e::lDirectDelay + gen::e::directTaps.back() + 1), 0.0f);
        for (size_t i = 0; i < apLine.size(); ++i)
            apLine[i].assign(pow2Above(gen::e::apDelay[i] + 1), 0.0f);
        int longest = 0; // feedback read P + 1, static read d + 1, moving reads centre + depth + 1 extra
        for (int k = 0; k < numLoops; ++k)
            longest = std::max({longest, gen::e::loopLength[(size_t)k] + 1, gen::e::staticReadDelay[(size_t)k] + 1});
        for (int i = 0; i < numReads; ++i)
            longest = std::max(longest, gen::e::readCentre[(size_t)i] + (int)std::ceil(gen::e::readDepth) + 2);
        for (auto& v : vLine)
            v.assign(pow2Above(longest + 1), 0.0f);
        lateLine.assign(pow2Above(gen::e::movingLDelay + 1), 0.0);
        setSettings(Settings{});
        reset();
    }

    void reset() override
    {
        std::fill(inLine.begin(), inLine.end(), 0.0f);
        for (auto& a : apLine)
            std::fill(a.begin(), a.end(), 0.0f);
        for (auto& v : vLine)
            std::fill(v.begin(), v.end(), 0.0f);
        std::fill(lateLine.begin(), lateLine.end(), 0.0);
        fbState = {};
        statState = {};
        pos = 0;
        movR = movL = 0.0;
        prevStat = 0.0;
        yR1 = yL1 = 0.0;
        moving = modulationOn;
        for (size_t i = 0; i < (size_t)numReads; ++i)
        {
            // fixed start phase (H: e_sim.e_reads' seed-0 draw), advanced one device sample at a time
            phase[i] = gen::e::readPhase[i];
            phaseStep[i] = gen::e::readRateHz[i] / deviceRate;
            offset[i] = readAddress(i);
        }
    }

    void setSettings(const Settings& s) override
    {
        quantise = s.sixteenBit;
        const double rt = std::clamp(s.rt, rtMin, rtMax);
        const double lf = std::clamp(s.lfX, 0.01, 100.0), hf = std::clamp(s.hfX, 0.01, 100.0);
        const double a = gen::e::lawShelfPole;
        for (size_t k = 0; k < (size_t)numLoops; ++k)
        {
            // e_sim.loop_iir: G(z) = g (1 - p) (b0 + b1 z^-1) / ((1 - a z^-1) (1 - p z^-1))
            const auto c = loopLaw((int)k, rt, lf, hf);
            const double b0 = (c.r0 * (1.0 - a) + (1.0 + a)) / 2.0, b1 = (c.r0 * (1.0 - a) - (1.0 + a)) / 2.0;
            const double gp = c.g * (1.0 - c.p);
            loop[k] = {gp * b0, gp * b1, a + c.p, -a * c.p};
        }
        const double ampl = std::max(s.amplPercent, 0.0) / gen::e::amplRef; // direct taps x AMPL / 80 (spec controls.AMPL, H)
        for (size_t j = 0; j < directR.size(); ++j)
        {
            directR[j] = gen::e::directGains[j] * ampl;
            directL[j] = gen::e::lDirectSign * gen::e::lGain * directR[j];
        }
    }

    void process(const float* in, float* outL, float* outR, int n) override
    {
        namespace g = gen::e;
        const auto inMask = inLine.size() - 1, lateMask = lateLine.size() - 1;
        const double gLate = g::lReverbSign * g::lGain, mlp = g::movingOnePole;
        static_assert(g::lReverbDelay == g::movingLDelay + 1, "L's static reverb rides one sample behind the moving sum");
        for (int i = 0; i < n; ++i)
        {
            inLine[pos & inMask] = in[i];
            const auto x = [&](int d) { return (double)inLine[(pos - (uint32_t)d) & inMask]; };

            // direct path, R and L (314 later, inverted, x gL)
            double yR = 0.0, yL = 0.0;
            for (size_t j = 0; j < directR.size(); ++j)
            {
                yR += directR[j] * x(g::k0 + g::directTaps[j]);
                yL += directL[j] * x(g::k0 + g::lDirectDelay + g::directTaps[j]);
            }

            // all-pass chain: w = u + g w[t - D], u' = -g w + w[t - D]
            double u = x(g::k0);
            for (size_t j = 0; j < apLine.size(); ++j)
            {
                auto& line = apLine[j];
                const auto mask = line.size() - 1;
                const double wd = line[(pos - (uint32_t)g::apDelay[j]) & mask];
                const double w = u + g::apGain[j] * wd;
                line[pos & mask] = store((float)w);
                u = wd - g::apGain[j] * w;
            }

            // loops: the feedback and the static read each run G on their read of v
            double stat = 0.0;
            for (size_t k = 0; k < (size_t)numLoops; ++k)
            {
                auto& v = vLine[k];
                const auto mask = v.size() - 1;
                const auto& c = loop[k];
                const auto at = [&](int d) { return (double)v[(pos - (uint32_t)d) & mask]; };
                auto& f = fbState[k];
                const double fb = c.b0 * at(g::loopLength[k]) + c.b1 * at(g::loopLength[k] + 1) + c.a1 * f[0] + c.a2 * f[1];
                f[1] = f[0];
                f[0] = fb;
                auto& s = statState[k];
                const double line = c.b0 * at(g::staticReadDelay[k]) + c.b1 * at(g::staticReadDelay[k] + 1) + c.a1 * s[0] + c.a2 * s[1];
                s[1] = s[0];
                s[0] = line;
                stat += g::staticReadGain[k] * line;
                v[pos & mask] = store((float)(u - fb));
            }

            // moving reads of v, R now and L (its final amplitude) 314 later; one-pole on both sums
            double mR = 0.0, mL = 0.0;
            if (moving)
            {
                for (size_t r = 0; r < (size_t)numReads; ++r)
                {
                    const auto& v = vLine[(size_t)g::readLoop[r]];
                    const auto mask = v.size() - 1;
                    const auto j = pos - (uint32_t)offset[r];
                    mR += g::readR[r] * v[j & mask];
                    mL += g::readL[r] * v[(j - (uint32_t)g::readLExtra[r]) & mask];
                    phase[r] += phaseStep[r];
                    if (phase[r] >= 1.0)
                        phase[r] -= 1.0;
                    offset[r] = readAddress(r); // the address for the next sample
                }
                movR = (1.0 - mlp) * mR + mlp * movR;
                movL = (1.0 - mlp) * mL + mlp * movL;
            }

            // L's late part: the moving sum 314 later and the static reverb 315 later (the previous sample's, x gL)
            lateLine[pos & lateMask] = movL + gLate * prevStat;
            prevStat = stat;
            yL += lateLine[(pos - (uint32_t)g::movingLDelay) & lateMask];
            yR += stat + movR;

            // output kernel on both channels
            outR[i] = (float)(g::outKernel[0] * yR + g::outKernel[1] * yR1);
            outL[i] = (float)(g::outKernel[0] * yL + g::outKernel[1] * yL1);
            yR1 = yR;
            yL1 = yL;
            ++pos;
        }
    }

    // Moving read r's current address (output delay after R's direct tap, device samples), for the tests.
    int readOffset(int r) const { return offset[(size_t)r]; }

  private:
    struct Biquad
    {
        double b0 = 0.0, b1 = 0.0, a1 = 0.0, a2 = 0.0; // y = b0 x + b1 x[-1] + a1 y[-1] + a2 y[-2]
    };

    static size_t pow2Above(int n)
    {
        size_t s = 1;
        while (s < (size_t)n)
            s <<= 1;
        return s;
    }

    // S's LF law (hi, dc, k) at an LF value: log-log interpolation between S's steps, held at the ends
    // (ml_lines.lf_law)
    static void lfLaw(double lf, double& hi, double& dc, double& k)
    {
        namespace g = gen::e;
        constexpr size_t n = g::lawLf.size();
        const double x = std::log(lf); // at a step this returns the table's value to a few ulps
        const auto logInterp = [&](const std::array<double, n>& y) {
            if (x <= std::log(g::lawLf[0]))
                return y[0];
            if (x >= std::log(g::lawLf[n - 1]))
                return y[n - 1];
            size_t i = 1;
            while (x > std::log(g::lawLf[i]))
                ++i;
            const double x0 = std::log(g::lawLf[i - 1]), x1 = std::log(g::lawLf[i]);
            const double y0 = std::log(y[i - 1]), y1 = std::log(y[i]);
            return std::exp(y0 + (x - x0) * (y1 - y0) / (x1 - x0));
        };
        hi = logInterp(g::lawLf_hi);
        dc = logInterp(g::lawLf_dc);
        k = logInterp(g::lawLf_k);
    }

    // The one-pole low-pass whose gain at f is T (ml_lines.onepole_for)
    static double onePoleFor(double t, double f)
    {
        const double c = std::cos(2.0 * juce::MathConstants<double>::pi * f / deviceRate);
        const double t2 = std::min(t * t, gen::e::onepoleT2Max);
        return ((1.0 - t2 * c) - std::sqrt((1.0 - t2 * c) * (1.0 - t2 * c) - (1.0 - t2) * (1.0 - t2))) / (1.0 - t2);
    }

    // np.interp: piecewise linear, held at the ends
    template <size_t N>
    static double interp(double x, const std::array<double, N>& xs, const std::array<double, N>& ys)
    {
        if (x <= xs[0])
            return ys[0];
        if (x >= xs[N - 1])
            return ys[N - 1];
        size_t i = 1;
        while (x > xs[i])
            ++i;
        return ys[i - 1] + (x - xs[i - 1]) * (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
    }

    // round(centre + A tri(phase)), tri = 4 |x - 0.5| - 1 as e_sim.e_reads; half to even like np.round (the default
    // rounding mode)
    int readAddress(size_t r) const
    {
        return (int)std::lrint((double)gen::e::readCentre[r] + gen::e::readDepth * (4.0 * std::abs(phase[r] - 0.5) - 1.0));
    }

    std::vector<float> inLine;
    std::array<std::vector<float>, gen::e::apDelay.size()> apLine;
    std::array<std::vector<float>, (size_t)numLoops> vLine;
    std::vector<double> lateLine;
    std::array<Biquad, (size_t)numLoops> loop{};
    std::array<std::array<double, 2>, (size_t)numLoops> fbState{}, statState{};
    std::array<double, gen::e::directTaps.size()> directR{}, directL{};
    std::array<double, (size_t)numReads> phase{}, phaseStep{};
    std::array<int, (size_t)numReads> offset{};
    uint32_t pos = 0;
    double movR = 0.0, movL = 0.0, prevStat = 0.0, yR1 = 0.0, yL1 = 0.0;
    bool moving = true;
};
} // namespace rv::dsp
