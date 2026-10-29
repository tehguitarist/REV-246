#pragma once

#include "dsp/engine/ProgramCore.h"
#include "dsp/generated/s_spec.h"

#include <juce_core/juce_core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace rv::dsp
{
// Program S at the device rate (analysis/s/s_plugin_spec.json, tools/s_spec.py; docs/model_s_results.md §17, §22,
// §32, §38.1). The spec's signal path is: input -> tap injection -> four loops (G_i in the feedback, delay P_i) per
// channel -> sum -> AP20 -> AP72 -> AP230 (damping b) -> output, plus the direct taps x AMPL / 90. Every part is
// linear and time-invariant, so the core runs the same network in the unit's order (§14.3: AP230 is on the input
// side; tools/loop_read/c/s_loops_fx.c): input -> AP20 -> AP72 -> AP230 -> four shared loops, and each channel reads
// its five signed taps per loop out of the loop lines (the spec's injection, as read-out). One diffuser and four
// loops instead of two and eight; the impulse response is the spec's.
//
// Precision: the delay lines are float (every write goes through ProgramCore::store(), which is float); filter
// states, coefficients and sums are double. Against the float32 reference IRs the null is -140 dB or lower in every
// window (STests.cpp), set by float storage on both sides, far below the -100 dB bar.
//
// Cost per device sample: 3 all-passes (8 multiplies with AP230's damping), 4 loops (gain + 3 first-order sections,
// 10 multiplies each), 40 signed read-out taps (adds, then one multiply per channel for the injection magnitude), 8
// direct taps (10 multiplies with AMPL): about 60 multiplies and 56 delay-line reads, the 2 x TMS32010's 50-70 taps.
class SCore final : public ProgramCore
{
  public:
    // One first-order section (b0 + b1 z^-1) / (1 + a1 z^-1).
    struct Section
    {
        double b0 = 1.0, b1 = 0.0, a1 = 0.0;
    };
    // G_i: scalar gain times up to three sections (the law's shelf, the HF one-pole, the cell's DC shelf), in that
    // order; unused sections stay identity.
    struct LoopFilter
    {
        double gain = 1.0;
        std::array<Section, 3> sections{};
        int count = 0;
    };
    struct Coefs
    {
        std::array<LoopFilter, 4> loops{};
        double ap230b = 0.0;
    };

    static constexpr int numLoops = 4;

    void prepare() override
    {
        for (auto& l : lines)
            l.assign(lineSize, 0.0f);
        for (size_t k = 0; k < aps.size(); ++k)
            aps[k].assign(apSize[k], 0.0f);
        history.assign(historySize, 0.0f);
        reset();
    }

    void reset() override
    {
        for (auto& l : lines)
            std::fill(l.begin(), l.end(), 0.0f);
        for (auto& a : aps)
            std::fill(a.begin(), a.end(), 0.0f);
        std::fill(history.begin(), history.end(), 0.0f);
        for (auto& s : sectionState)
            s.fill(0.0);
        pos = 0;
    }

    void setSettings(const Settings& s) override
    {
        quantise = s.sixteenBit;
        coefs = law(s.lfX, s.rt, s.hfX);
        directScale = s.amplPercent / gen::s::amplReference;
    }

    void process(const float* in, float* outL, float* outR, int n) override
    {
        namespace g = gen::s;
        const double ag = g::allpassGain, b = coefs.ap230b, f0 = 1.0 - b;
        for (int i = 0; i < n; ++i)
        {
            const double x = in[i];
            history[pos & historyMask] = in[i];

            // Diffuser, H = (-g + z^-M F) / (1 - g z^-M F): q[n] = x + g F(q[n - M]), y = -g q[n] + F(q[n - M]).
            double v = x;
            for (size_t k = 0; k < 2; ++k)
            {
                auto& a = aps[k];
                const auto mask = apSize[k] - 1;
                const double d = a[(pos - (size_t)g::allpassDelay[k]) & mask];
                const float q = store((float)(v + ag * d));
                a[pos & mask] = q;
                v = d - ag * q;
            }
            {
                auto& a = aps[2];
                const auto mask = apSize[2] - 1;
                const auto m = (size_t)g::allpassDelay[2];
                const double d = f0 * a[(pos - m) & mask] + b * a[(pos - m - 1) & mask];
                const float q = store((float)(v + ag * d));
                a[pos & mask] = q;
                v = d - ag * q;
            }

            // Loops: w_i[n] = v[n] + G_i(w_i[n - P_i]).
            for (size_t l = 0; l < (size_t)numLoops; ++l)
            {
                auto& line = lines[l];
                const auto& lf = coefs.loops[l];
                double y = lf.gain * line[(pos - (size_t)g::loopDelay[l]) & lineMask];
                for (int k = 0; k < lf.count; ++k)
                {
                    const auto& sec = lf.sections[(size_t)k];
                    auto& z = sectionState[l][(size_t)k];
                    const double o = sec.b0 * y + z; // transposed direct form II
                    z = sec.b1 * y - sec.a1 * o;
                    y = o;
                }
                line[pos & lineMask] = store((float)(v + y));
            }

            // Read-out taps (the spec's injection) and the direct path.
            double acc[2] = {0.0, 0.0};
            for (size_t t = 0; t < g::tapDelay.size(); ++t)
            {
                const double r = lines[(size_t)g::tapLoop[t]][(pos - (size_t)g::tapDelay[t]) & lineMask];
                acc[g::tapChannel[t]] += g::tapSign[t] > 0 ? r : -r;
            }
            double dl = 0.0, dr = 0.0;
            for (size_t t = 0; t < g::directLDelay.size(); ++t)
            {
                dl += g::directLGain[t] * history[(pos - (size_t)g::directLDelay[t]) & historyMask];
                dr += g::directRGain[t] * history[(pos - (size_t)g::directRDelay[t]) & historyMask];
            }
            outL[i] = (float)(g::tapMagnitude * acc[0] + directScale * dl);
            outR[i] = (float)(g::tapMagnitude * acc[1] + directScale * dr);
            ++pos;
        }
    }

    const Coefs& currentCoefs() const { return coefs; }

    // The loop-filter law and AP230's damping step for any LF / RT / HF (tools/s_spec.py loop_filter; spec
    // loop_filter_law). At the library's 128 steps it reproduces the spec's resolved settings. Between them (H, spec
    // "nearest RT fold, nearest sampled cell", §32.8): the law runs at the panel values, LF's law numbers (hi, dc, k)
    // log-interpolated between the LF steps as M / L's lf_law does; the RT fold, the cell's shelf terms and LF's
    // step-keyed corrections come from the nearest library step (LF and RT nearest in log, HF in linear; ties go up);
    // the HF pole is the cell's own only at a library setting, elsewhere §32.6's pole curve at the law's pole; AP230's
    // step is the nearest (RT, HF) cell's.
    static Coefs law(double lfIn, double rtIn, double hfIn)
    {
        namespace g = gen::s;
        constexpr double twoPi = juce::MathConstants<double>::twoPi;
        const double fs = g::deviceRate;

        double lf = std::clamp(lfIn, g::lfSteps.front(), g::lfSteps.back());
        double rt = std::clamp(rtIn, 0.05, 20.0);
        double hf = std::clamp(hfIn, 0.01, 1.0);
        bool lfOn = true, on = true;
        const auto iLf = nearest(g::lfSteps, lf, true, lfOn);
        const auto iRt = nearest(g::rtSteps, rt, true, on);
        const auto iHf = nearest(g::hfSteps, hf, false, on);
        on = on && lfOn; // a library setting: all three on their steps
        const auto cell = (iLf * g::rtSteps.size() + iRt) * g::hfSteps.size() + iHf;
        const auto foldLf = iRt * g::lfSteps.size(); // + LF step

        // The law's LF nodes with the fold's and the cell's overrides [s_spec.loop_filter].
        auto hi = g::lawHi, dcn = g::lawDc, kn = g::lawK;
        for (size_t j = 0; j < g::lfSteps.size(); ++j)
            if (g::foldLfDcOn[foldLf + j] != 0)
                dcn[j] = g::foldLfDc[foldLf + j];
        double gamma = g::foldGamma[iRt], dcMin = g::dcMin, beta = g::dcBeta;
        bool gammaDc = true;
        const int kind = g::cellKind[cell];
        if (kind == 1)
        {
            dcn[0] = g::cellDc[cell]; // the cell's own DC ratio for LF 0.1's law
            gamma = g::cellGamma[cell];
            dcMin = g::cellDcMin[cell];
            beta = g::cellBeta[cell];
            gammaDc = g::cellGammaDc[cell] != 0;
        }
        double Lhi, Ldc, Lk;
        if (lfOn)
            Lhi = hi[iLf], Ldc = dcn[iLf], Lk = kn[iLf];
        else
            Lhi = logInterp(hi, lf), Ldc = logInterp(dcn, lf), Lk = logInterp(kn, lf);

        Coefs c;
        for (size_t i = 0; i < (size_t)numLoops; ++i)
        {
            const double P = g::loopDelay[i];
            const double kk = -3.0 * P / fs;
            double gain = std::pow(10.0, kk / (Lhi * rt));
            LoopFilter f;
            if (kind == 2 && Ldc < Lhi)
            {
                // The floor_form branch: a soft DC floor in dB per pass.
                const double r0 = std::pow(10.0, kk * (1.0 / (Ldc * rt) - 1.0 / (Lhi * rt)));
                const double x = 20.0 * std::log10(gain * r0), fd = 20.0 * std::log10(g::cellFloorF[cell]);
                const double cc = g::cellFloorC[cell];
                const double dd = std::max(x, fd - cc) + cc * std::exp(-std::max(x - fd + cc, 0.0) / g::cellFloorTau[cell]);
                const double u = 20.0 * std::log10(g::cellFloorF0[cell]) - x;
                if (u > 0.0)
                    gain *= std::pow(10.0, (-g::cellFloorBeta[cell] * u + g::cellFloorGamma[cell] * 60.0 / (Lhi * rt) * P / fs) / 20.0);
                f.sections[(size_t)f.count++] = shelf(std::pow(10.0, dd / 20.0) / gain);
            }
            else if (Ldc < Lhi || Ldc > Lhi) // the spec's dc != hi
            {
                // The shelf with the DC floor g r0 >= dc_min [§22.3].
                const double r0 = std::pow(10.0, kk * (1.0 / (Ldc * rt) - 1.0 / (Lhi * rt)));
                bool floored = false;
                if (gain * r0 < dcMin)
                {
                    floored = true;
                    const double lift = 20.0 * std::log10(dcMin / (gain * r0));
                    gain *= std::pow(10.0, -beta * lift / 20.0);
                }
                const double gammaDb = gamma * 60.0 / (Lhi * rt) * P / fs;
                if (floored && gamma != 0.0 && !gammaDc)
                    gain *= std::pow(10.0, gammaDb / 20.0);
                f.sections[(size_t)f.count++] = shelf(std::max(r0, dcMin / gain));
                if (floored && gamma != 0.0 && gammaDc)
                    gain *= std::pow(10.0, gammaDb / 20.0); // the spec's "gain" section, folded into the scalar
            }
            if (hf < 1.0)
            {
                // The HF one-pole at f_ref, the control floored at HFRT_MIN / RT [§29].
                double p;
                if (on && g::cellPolesOn[cell] != 0)
                    p = g::cellPoles[cell * 4 + i]; // the cell's held-out curve value [§32.6, §32.10]
                else
                {
                    const double hfE = std::max(hf, g::hfRtMin / rt);
                    const double T = std::pow(10.0, kk / rt * Lk * (1.0 / hfE - 1.0));
                    p = poleCurve(onePoleFor(T, std::cos(twoPi * g::hfRefHz / fs)));
                }
                f.sections[(size_t)f.count++] = Section{1.0 - p, 0.0, -p};
            }
            const auto fl = foldLf + iLf;
            if (g::foldLfRtOn[fl] != 0)
            {
                // LF 0.4 / 2.0's correction in RT [§22.2].
                const double a = g::foldLfRtA[fl];
                double rel = a + g::foldLfRtB[fl] * rt;
                if (g::foldLfRtClip[fl] != 0 && rel * a <= 0.0)
                    rel = 0.0;
                gain *= std::pow(10.0, rel * (60.0 / (Lhi * rt)) * P / (20.0 * fs));
            }
            // The cell's per-loop gain and DC shelf [§32.2-32.3].
            gain *= std::pow(10.0, g::cellLoopGainDb[cell * 4 + i] / 20.0);
            f.sections[(size_t)f.count++] = shelf(std::pow(10.0, g::cellLoopDcDb[cell * 4 + i] / 20.0));
            f.gain = gain;
            c.loops[i] = f;
        }
        c.ap230b = g::ap230BSteps[(size_t)g::ap230Step[iRt * g::hfSteps.size() + iHf]];
        return c;
    }

  private:
    // Index of the nearest step (log or linear distance, ties to the upper step); a value within 1e-6 relative of a
    // step snaps to it exactly, and `onStep` stays true only while every value has snapped.
    template <size_t N>
    static size_t nearest(const std::array<double, N>& steps, double& v, bool logScale, bool& onStep)
    {
        size_t j = 0;
        while (j + 1 < N && (logScale ? v >= std::sqrt(steps[j] * steps[j + 1]) : v >= 0.5 * (steps[j] + steps[j + 1])))
            ++j;
        if (std::abs(v - steps[j]) <= 1e-6 * steps[j])
            v = steps[j];
        else
            onStep = false;
        return j;
    }

    // Piecewise linear in log LF between the LF steps (clamped at the ends), in log of the value.
    static double logInterp(const std::array<double, 4>& y, double lf)
    {
        const auto& x = gen::s::lfSteps;
        size_t j = 0;
        while (j + 2 < x.size() && lf > x[j + 1])
            ++j;
        const double t = (std::log(lf) - std::log(x[j])) / (std::log(x[j + 1]) - std::log(x[j]));
        return std::exp(std::log(y[j]) + t * (std::log(y[j + 1]) - std::log(y[j])));
    }

    // Shelf (b0 + b1 z^-1) / (1 - a z^-1) with DC gain `dc`, Nyquist gain 1 and the law's pole a (model_s.shelf).
    static Section shelf(double dc)
    {
        const double a = gen::s::shelfPole;
        return {(dc * (1.0 - a) + (1.0 + a)) / 2.0, (dc * (1.0 - a) - (1.0 + a)) / 2.0, -a};
    }

    // Pole p of (1 - p) / (1 - p z^-1) with |H| = T at the frequency whose cosine is c (q1.onepole_for).
    static double onePoleFor(double T, double c)
    {
        const double t2 = std::min(T * T, 0.999999);
        const double u = 1.0 - t2 * c;
        return (u - std::sqrt(u * u - (1.0 - t2) * (1.0 - t2))) / (1.0 - t2);
    }

    // §32.6's all-data pole curve of the law's pole, with the drop-back past its fold point [§32.4] (ml_lines).
    static double poleCurve(double p)
    {
        namespace g = gen::s;
        if (p >= g::poleCurveFold)
            return g::poleCurveDropBack;
        const auto& x = g::poleCurveKnots;
        const auto& y = g::poleCurveValues;
        double out;
        if (p <= x.front())
            out = y.front();
        else if (p >= x.back())
            out = y.back();
        else
        {
            size_t j = 0;
            while (p > x[j + 1])
                ++j;
            out = y[j] + (p - x[j]) / (x[j + 1] - x[j]) * (y[j + 1] - y[j]);
        }
        return std::max(out, 0.0);
    }

    static constexpr size_t lineSize = 2048; // > every P_i and read-out tap (max 1853)
    static constexpr size_t lineMask = lineSize - 1;
    static constexpr std::array<size_t, 3> apSize{32, 128, 256};
    static constexpr size_t historySize = 512; // > the longest direct tap (353)
    static constexpr size_t historyMask = historySize - 1;

    std::array<std::vector<float>, numLoops> lines;
    std::array<std::vector<float>, 3> aps;
    std::vector<float> history;
    std::array<std::array<double, 3>, numLoops> sectionState{};
    Coefs coefs = law(1.0, 1.0, 1.0);
    double directScale = 1.0;
    size_t pos = 0;
};
} // namespace rv::dsp
