#pragma once

// Programs M and L: one class, a program flag (docs/ENGINE.md; analysis/ml/ml_plugin_spec.json, docs/ml_structure.md
// §7-9, §14-26). At the device rate, in the spec's read-out form (modulation.application):
//
//   input -> diffuser (three all-passes, the last damped by b) -> x
//   u_i = x + sum_j F(theta)_ji v_j,  v_j = G_j applied to u_j delayed P_j     (four lines, coupling F(theta))
//   reverb_ch = sum of the channel's 20 sign-only read taps on u
//   out_ch = h * (cluster gain x direct taps_ch + reverb_ch)                    (causal output kernel)
//
// With theta = 0 this is the static network of the reference IRs (tools/ml_spec.py simulate(), the transposed form:
// 1^T (I - Lambda F)^-1 T = T^T (I - F^T Lambda)^-1 1, and the diffuser commutes with the lines). theta(t) is one
// angle of both 2 x 2 blocks of F, F(theta) = F_b R(theta) per block, so the coupling stays orthogonal at every sample.
//
// Precision: the delay lines (the four line buffers and the three all-pass buffers) are float, written through
// ProgramCore::store(); every filter state, the coupling sums, the LFO and the output kernel's history are double.
// That nulls every reference IR at -141 dB or better in every window up to L's RT 20 s (tests/dsp/MLTests.cpp), at
// the float32 references' own floor and far under the -100 dB bar; the float line storage (store()'s type) is the
// only rounding inside the loops, and doubles elsewhere keep the shelf (pole 0.917) and HF pole (up to 0.85) from
// adding their own rounding to it. Cost: about 160 (M) / 240 (L) multiply-adds per device sample, the output kernel
// (2 x 21 / 2 x 61 taps) the largest part.
//
// 16-bit mode (Settings::sixteenBit; S's scheme, model_s_results §36.5, ProgramCore): the all-pass and line writes
// and each line filter's two stored outputs (shelf, then HF one-pole: ProgramCore::ShelfPole16) are truncated,
// saturated 16-bit words; coefficients, coupling and LFO stay double (H: their format is unknown, §36.4 shows the
// loops need no more than Q12). Scaling (H, the no-overflow rule of §36.5's scaling A): the lines at
// k = (1 - max_j max|G_j|) / 2 (the coupling is orthogonal, so for a steady full-scale input |u_j| <= ||u|| <= 2 /
// (1 - max|G|)), the all-passes at k = 1 - |g|. The float filter and the 16-bit cascade keep one input / output
// history, so switching the mode carries the filter state over.

#include "dsp/engine/ProgramCore.h"
#include "dsp/generated/ml_spec.h"

#include <juce_core/juce_core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace rv::dsp
{
namespace mlsize
{
template <size_t N>
constexpr int maxOf(const std::array<int, N>& a)
{
    int m = 0;
    for (auto v : a)
        m = std::max(m, v);
    return m;
}
constexpr uint32_t pow2Above(int n)
{
    uint32_t p = 1;
    while ((int)p <= n)
        p <<= 1;
    return p;
}
// The longest read of a line buffer (a line's delay or a read tap) in either program.
constexpr int longestLineRead()
{
    namespace g = gen::ml;
    return std::max({maxOf(g::m::lineDelay), maxOf(g::m::tapLDelay), maxOf(g::m::tapRDelay), maxOf(g::l::lineDelay),
                     maxOf(g::l::tapLDelay), maxOf(g::l::tapRDelay)});
}
} // namespace mlsize

class MLCore final : public ProgramCore
{
  public:
    // L's next library RT step above 2.2: where its 0.40 HF floor has blended into the general one (law(), H).
    static constexpr double l22BlendEnd = 2.6;

    // One line's per-pass filter: gain, shelf DC gain r0, HF pole p (spec lines.filter_form; p = 0 bypasses).
    struct LineFilter
    {
        double gain = 0.0, shelfDc = 1.0, hfPole = 0.0;
    };
    // Everything the panel's LF / RT / HF set: the four line filters and the last all-pass's damping b.
    struct Coefs
    {
        std::array<LineFilter, 4> lines{};
        double damping = 0.0;
    };

    explicit MLCore(params::Program p = params::Program::m) : isL(p == params::Program::l)
    {
        applyCoefs(law(program(), lfX, rtS, hfX));
    }

    // M or L (anything but L means M). Clears the state, as the delays change. Real-time safe after prepare().
    void setProgram(params::Program p)
    {
        const bool l = p == params::Program::l;
        if (l == isL)
            return;
        isL = l;
        applyCoefs(law(p, lfX, rtS, hfX));
        reset();
    }
    params::Program program() const { return isL ? params::Program::l : params::Program::m; }

    void prepare() override
    {
        // The worst case of both programs, so switching never allocates.
        lineBuf.assign((size_t)(4 * lineSize), 0.0f);
        for (auto& b : apBuf)
            b.assign((size_t)apSize, 0.0f);
        inHist.assign((size_t)inSize, 0.0f);
        for (auto& h : outHist)
            h.assign((size_t)kernelSize, 0.0);
        reset();
    }

    void reset() override
    {
        std::fill(lineBuf.begin(), lineBuf.end(), 0.0f);
        for (auto& b : apBuf)
            std::fill(b.begin(), b.end(), 0.0f);
        std::fill(inHist.begin(), inHist.end(), 0.0f);
        for (auto& h : outHist)
            std::fill(h.begin(), h.end(), 0.0);
        for (auto& f : filt)
            f.s1 = f.s2 = 0.0;
        for (auto& f : filt16)
            f.s1 = 0.0;
        hist = {};
        pos = 0;
        modActive = modulationOn;
        // theta's start phase is fixed (0: theta starts at 0 and rising), so renders repeat (H: the unit's is random).
        lfoSin = 0.0;
        lfoCos = 1.0;
        thetaNow = 0.0;
    }

    void setSettings(const Settings& s) override
    {
        lfX = s.lfX;
        rtS = s.rt;
        hfX = s.hfX;
        // AMPL scales the cluster only; AMPL 90 = the reference IRs' unit cluster (H: the absolute scale is not
        // resolved, spec controls.AMPL; the same convention as S's AMPL / 90).
        clusterGain = s.amplPercent / 90.0;
        const bool was16 = quantise;
        quantise = s.sixteenBit;
        applyCoefs(law(program(), lfX, rtS, hfX));
        if (quantise != was16)
            carryFilterState();
    }

    // M or L comes from the constructor / setProgram(), not from Settings::program (one core per program).
    void process(const float* in, float* outL, float* outR, int n) override
    {
        if (lineBuf.empty()) // not prepared: silence rather than reading unallocated buffers
        {
            std::fill(outL, outL + n, 0.0f);
            std::fill(outR, outR + n, 0.0f);
            return;
        }
        const auto& t = tables();
        const double g = gen::ml::allpassG;
        const double b = coefs.damping, b0 = 1.0 - b;
        for (int i = 0; i < n; ++i, ++pos)
        {
            // Direct path (the cluster), on the raw input.
            inHist[pos & inMask] = in[i];
            std::array<double, 2> cl{};
            for (int ch = 0; ch < 2; ++ch)
                for (size_t k = 0; k < 4; ++k)
                    cl[(size_t)ch] += t.clusterSign[(size_t)ch][k] * (double)inHist[(pos - (uint32_t)t.clusterDelay[(size_t)ch][k]) & inMask];

            // Diffuser: (-g + z^-M F) / (1 - g z^-M F), F = 1 for the first two, (1 - b) + b z^-1 for the last.
            double x = in[i];
            for (size_t k = 0; k < 3; ++k)
            {
                auto& buf = apBuf[k];
                const auto d = (uint32_t)t.diffuserDelay[k];
                double f = buf[(pos - d) & apMask];
                if (k == 2)
                    f = b0 * f + b * (double)buf[(pos - d - 1) & apMask];
                const double w = x + g * f;
                buf[pos & apMask] = store(w, gridAp);
                x = f - g * w;
            }

            // The lines: v_j = G_j (u_j delayed P_j); G_j = (B0 + B1 z^-1) / (1 + A1 z^-1 + A2 z^-2), transposed DF II.
            std::array<double, 4> v{};
            for (size_t j = 0; j < 4; ++j)
            {
                const double d = lineBuf[j * lineSize + ((pos - (uint32_t)t.lineDelay[j]) & lineMask)];
                auto& h = hist[j];
                double y;
                if (quantise)
                    y = runShelfPole16(filt16[j], d, h.x1, h.y1, gridLine);
                else
                {
                    auto& f = filt[j];
                    y = f.b0 * d + f.s1;
                    f.s1 = f.b1 * d - f.a1 * y + f.s2;
                    f.s2 = -f.a2 * y;
                }
                h = {d, y, h.y1};
                v[j] = y;
            }

            // theta: a quadrature LFO at 0.83 Hz (renormalised every sample), cos / sin of theta by series (|theta|
            // <= 0.22 rad: error below 1e-14).
            double c = 1.0, s = 0.0;
            if (modActive)
            {
                const double ns = lfoSin * lfoStepCos + lfoCos * lfoStepSin;
                const double nc = lfoCos * lfoStepCos - lfoSin * lfoStepSin;
                const double k = 1.5 - 0.5 * (ns * ns + nc * nc);
                lfoSin = ns * k;
                lfoCos = nc * k;
                thetaNow = gen::ml::thetaAmplitude * lfoSin;
                rotation(thetaNow, c, s);
            }

            // Coupling, block by block: u_cols = R(theta)^T F_b^T v_rows [spec modulation.coupling].
            std::array<double, 4> u{};
            for (const auto& blk : blocks)
            {
                const double w0 = F(blk.rows[0], blk.cols[0]) * v[blk.rows[0]] + F(blk.rows[1], blk.cols[0]) * v[blk.rows[1]];
                const double w1 = F(blk.rows[0], blk.cols[1]) * v[blk.rows[0]] + F(blk.rows[1], blk.cols[1]) * v[blk.rows[1]];
                u[blk.cols[0]] = x + c * w0 + s * w1;
                u[blk.cols[1]] = x - s * w0 + c * w1;
            }
            for (size_t j = 0; j < 4; ++j)
                lineBuf[j * lineSize + (pos & lineMask)] = store(u[j], gridLine);

            // Read taps, the cluster, then the output kernel (causal: its centre lands on the cluster's tap 0).
            for (int ch = 0; ch < 2; ++ch)
            {
                double rev = 0.0;
                for (size_t k = 0; k < numTaps; ++k)
                    rev += t.tapSign[(size_t)ch][k]
                           * (double)lineBuf[(size_t)t.tapLoop[(size_t)ch][k] * lineSize
                                             + ((pos - (uint32_t)t.tapDelay[(size_t)ch][k]) & lineMask)];
                auto& h = outHist[(size_t)ch];
                h[pos & kernelMask] = clusterGain * cl[(size_t)ch] + rev;
                double y = 0.0;
                for (int k = 0; k < t.kernelTaps; ++k)
                    y += t.kernel[k] * h[(pos - (uint32_t)k) & kernelMask];
                (ch == 0 ? outL : outR)[i] = (float)y;
            }
        }
    }

    // The coefficient law [spec lines.law, diffuser.damping] for any LF / RT / HF in range. At the library's steps it
    // is the spec's (its resolved settings are the test vectors); between them H: LF log-interpolated between S's LF
    // steps, RT and HF continuous, b from the nearest (RT, HF) cell (RT nearest in log, HF in value).
    static Coefs law(params::Program prog, double lf, double rt, double hf)
    {
        namespace g = gen::ml;
        const auto& t = tables(prog == params::Program::l);
        rt = std::max(rt, 1e-3);
        lf = std::max(lf, 1e-3);
        hf = std::max(hf, 1e-3);
        Coefs out;
        out.damping = g::dampingTable[(size_t)dampingStep(t, rt, hf)];

        double hi = 0.0, dc = 0.0, k = 0.0;
        lfLaw(lf, hi, dc, k);
        // The HF control saturates: RT x HF >= 0.12 s x the program scale (x2 M, x4 L); L at its RT 2.2 step acts as
        // max(HF, 0.40) [§8.6, model_s_results §29]. Off the library's steps (H): L keeps the 0.40 floor below 2.2
        // (its lowest library RT; the panel's 2.0 would otherwise be darker than 2.2) and blends it in log RT to the
        // general floor at the next library step, 2.6, so neither the panel's steps nor typed values jump.
        hf = std::max(hf, g::hfrtMin * t.scale / rt);
        if (t.scale == 4 && rt < l22BlendEnd)
        {
            const double w = std::clamp(std::log(rt / 2.2) / std::log(l22BlendEnd / 2.2), 0.0, 1.0);
            const double general = g::hfrtMin * t.scale / l22BlendEnd;
            hf = std::max(hf, std::exp((1.0 - w) * std::log(g::l22Hf) + w * std::log(general)));
        }
        const double rtHi = hi * rt, rtDc = dc * rt;
        for (size_t j = 0; j < 4; ++j)
        {
            const double kk = -3.0 * t.lineDelay[j] / deviceRate;
            double gain = std::pow(10.0, kk / rtHi);
            double r0 = std::pow(10.0, kk * (1.0 / rtDc - 1.0 / rtHi));
            if (gain * r0 < g::dcMin) // the DC floor lifts g a little and sets the shelf [§16.15-16.16]
            {
                const double lift = 20.0 * std::log10(g::dcMin / (gain * r0));
                gain *= std::pow(10.0, -g::beta * lift / 20.0);
                r0 = g::dcMin / gain;
            }
            const double T = hf < 1.0 ? std::pow(10.0, kk / rt * k * (1.0 / hf - 1.0)) : 1.0;
            double p = T < 1.0 ? onePoleFor(T) : 0.0;
            if (p > 0.0) // the pole curve with its drop-back past the fold point [§25.20, model_s_results §32.6]
                p = std::max(p < g::poleFold ? interp(p, g::poleKnots.data(), g::poleValues.data(), (int)g::poleKnots.size())
                                             : g::poleDropBack,
                             0.0);
            out.lines[j] = {gain, r0, p};
        }
        return out;
    }

    // F(theta): both 2 x 2 blocks of F rotated by theta (F_b R(theta)), rows: lines, columns: loops.
    static std::array<double, 16> couplingAt(double theta)
    {
        double c = 1.0, s = 0.0;
        rotation(theta, c, s);
        std::array<double, 16> out{};
        for (const auto& blk : blocks)
            for (size_t r = 0; r < 2; ++r)
            {
                const double f0 = F(blk.rows[r], blk.cols[0]), f1 = F(blk.rows[r], blk.cols[1]);
                out[blk.rows[r] * 4 + blk.cols[0]] = f0 * c + f1 * s;
                out[blk.rows[r] * 4 + blk.cols[1]] = -f0 * s + f1 * c;
            }
        return out;
    }

    const Coefs& coefficients() const { return coefs; }
    // theta at the last processed sample (rad); 0 with the modulation off.
    double theta() const { return thetaNow; }

  private:
    // One program's tables from the generated header, as views.
    struct Tables
    {
        const int* lineDelay;
        int scale;
        const int* diffuserDelay;
        std::array<const int*, 2> clusterDelay, clusterSign, tapDelay, tapSign, tapLoop;
        const double* kernel;
        int kernelTaps;
        const double* rtSteps;
        int numRt;
        const double* hfSteps;
        int numHf;
        const int* dampingStep;
    };

    static constexpr size_t numTaps = gen::ml::m::tapLDelay.size();
    static_assert(gen::ml::l::tapLDelay.size() == numTaps && gen::ml::m::tapRDelay.size() == numTaps
                      && gen::ml::l::tapRDelay.size() == numTaps,
                  "both programs read 20 taps per channel");
    static_assert(gen::ml::thetaAmplitude < 0.5, "rotation()'s series is accurate for |theta| < 0.5 only");

    static const Tables& tables(bool l)
    {
        namespace m = gen::ml::m;
        namespace L = gen::ml::l;
        static constexpr Tables tm{m::lineDelay.data(), m::programScale, m::diffuserDelay.data(),
                               {m::clusterLDelay.data(), m::clusterRDelay.data()}, {m::clusterLSign.data(), m::clusterRSign.data()},
                               {m::tapLDelay.data(), m::tapRDelay.data()}, {m::tapLSign.data(), m::tapRSign.data()},
                               {m::tapLLoop.data(), m::tapRLoop.data()}, m::kernel.data(), (int)m::kernel.size(),
                               m::rtSteps.data(), (int)m::rtSteps.size(), m::hfSteps.data(), (int)m::hfSteps.size(),
                               m::dampingStep.data()};
        static constexpr Tables tl{L::lineDelay.data(), L::programScale, L::diffuserDelay.data(),
                               {L::clusterLDelay.data(), L::clusterRDelay.data()}, {L::clusterLSign.data(), L::clusterRSign.data()},
                               {L::tapLDelay.data(), L::tapRDelay.data()}, {L::tapLSign.data(), L::tapRSign.data()},
                               {L::tapLLoop.data(), L::tapRLoop.data()}, L::kernel.data(), (int)L::kernel.size(),
                               L::rtSteps.data(), (int)L::rtSteps.size(), L::hfSteps.data(), (int)L::hfSteps.size(),
                               L::dampingStep.data()};
        return l ? tl : tm;
    }
    const Tables& tables() const { return tables(isL); }

    // Buffer sizes: powers of two over the longest read of either program.
    static constexpr uint32_t lineSize = mlsize::pow2Above(mlsize::longestLineRead());
    static constexpr uint32_t apSize = mlsize::pow2Above(std::max(mlsize::maxOf(gen::ml::m::diffuserDelay),
                                                                  mlsize::maxOf(gen::ml::l::diffuserDelay)) + 1);
    static constexpr uint32_t inSize = mlsize::pow2Above(std::max({mlsize::maxOf(gen::ml::m::clusterLDelay),
                                                                   mlsize::maxOf(gen::ml::m::clusterRDelay),
                                                                   mlsize::maxOf(gen::ml::l::clusterLDelay),
                                                                   mlsize::maxOf(gen::ml::l::clusterRDelay)}));
    static constexpr uint32_t kernelSize = mlsize::pow2Above((int)std::max(gen::ml::m::kernel.size(), gen::ml::l::kernel.size()));
    static constexpr uint32_t lineMask = lineSize - 1, apMask = apSize - 1, inMask = inSize - 1, kernelMask = kernelSize - 1;

    // The coupling's two blocks: (lines, loops) [ml_variance.BLOCKS, ml_structure §7.3].
    struct Block
    {
        std::array<size_t, 2> rows, cols;
    };
    static constexpr std::array<Block, 2> blocks{{{{0, 2}, {1, 3}}, {{1, 3}, {0, 2}}}};
    static constexpr double F(size_t row, size_t col) { return gen::ml::coupling[row * 4 + col]; }

    static void rotation(double th, double& c, double& s)
    {
        const double t2 = th * th;
        s = th * (1.0 - t2 / 6.0 * (1.0 - t2 / 20.0 * (1.0 - t2 / 42.0 * (1.0 - t2 / 72.0 * (1.0 - t2 / 110.0)))));
        c = 1.0 - t2 / 2.0 * (1.0 - t2 / 12.0 * (1.0 - t2 / 30.0 * (1.0 - t2 / 56.0 * (1.0 - t2 / 90.0 * (1.0 - t2 / 132.0)))));
    }

    // numpy.interp: linear between the knots, the end values outside.
    static double interp(double x, const double* xp, const double* fp, int n)
    {
        if (x <= xp[0])
            return fp[0];
        if (x >= xp[n - 1])
            return fp[n - 1];
        int j = 0;
        while (x >= xp[j + 1])
            ++j;
        const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
        return slope * (x - xp[j]) + fp[j];
    }

    // S's LF law numbers at an LF value, log-interpolated between S's LF steps (ml_lines.lf_law).
    static void lfLaw(double lf, double& hi, double& dc, double& k)
    {
        namespace g = gen::ml;
        constexpr int n = (int)g::lawLf.size();
        for (size_t i = 0; i < (size_t)n; ++i)
            if (std::equal_to<double>{}(lf, g::lawLf[i])) // a step's own numbers, exactly (as ml_lines)
            {
                hi = g::lawHi[i], dc = g::lawDc[i], k = g::lawK[i];
                return;
            }
        std::array<double, n> x{}, yh{}, yd{}, yk{};
        for (size_t i = 0; i < (size_t)n; ++i)
        {
            x[i] = std::log(g::lawLf[i]);
            yh[i] = std::log(g::lawHi[i]);
            yd[i] = std::log(g::lawDc[i]);
            yk[i] = std::log(g::lawK[i]);
        }
        const double l = std::log(lf);
        hi = std::exp(interp(l, x.data(), yh.data(), n));
        dc = std::exp(interp(l, x.data(), yd.data(), n));
        k = std::exp(interp(l, x.data(), yk.data(), n));
    }

    // The HF one-pole with magnitude T at f_ref (ml_lines.onepole_for).
    static double onePoleFor(double T)
    {
        const double c = std::cos(2.0 * juce::MathConstants<double>::pi * gen::ml::hfRefHz / deviceRate);
        const double t2 = std::min(T * T, 0.999999);
        const double q = 1.0 - t2 * c;
        return (q - std::sqrt(q * q - (1.0 - t2) * (1.0 - t2))) / (1.0 - t2);
    }

    // The damping step of the nearest (RT, HF) cell.
    static int dampingStep(const Tables& t, double rt, double hf)
    {
        int ir = 0, ih = 0;
        for (int i = 1; i < t.numRt; ++i)
            if (std::abs(std::log(rt / t.rtSteps[i])) < std::abs(std::log(rt / t.rtSteps[ir])))
                ir = i;
        for (int i = 1; i < t.numHf; ++i)
            if (std::abs(hf - t.hfSteps[i]) < std::abs(hf - t.hfSteps[ih]))
                ih = i;
        return t.dampingStep[ir * t.numHf + ih];
    }

    void applyCoefs(const Coefs& c)
    {
        coefs = c;
        const double a = gen::ml::shelfPole;
        for (size_t j = 0; j < 4; ++j)
        {
            const auto& l = c.lines[j];
            const double gp = l.gain * (1.0 - l.hfPole);
            filt[j].b0 = gp * (l.shelfDc * (1.0 - a) + (1.0 + a)) / 2.0;
            filt[j].b1 = gp * (l.shelfDc * (1.0 - a) - (1.0 + a)) / 2.0;
            filt[j].a1 = -(a + l.hfPole);
            filt[j].a2 = a * l.hfPole;
        }
        // 16-bit mode: the same G_j as shelf (gain folded in) then one-pole, and the scaling (see the header comment)
        double peak = 0.0;
        for (size_t j = 0; j < 4; ++j)
        {
            const auto& l = c.lines[j];
            auto& f = filt16[j];
            f.c0 = l.gain * (l.shelfDc * (1.0 - a) + (1.0 + a)) / 2.0;
            f.c1 = l.gain * (l.shelfDc * (1.0 - a) - (1.0 + a)) / 2.0;
            f.a = a;
            f.p = l.hfPole;
            for (int k = 0; k < gridPoints; ++k)
            {
                const double cw = gridCos(k);
                peak = std::max(peak, sectionPower(f.c0, f.c1, -a, cw) * sectionPower(1.0 - f.p, 0.0, -f.p, cw));
            }
        }
        setGridGain(gridLine, std::max(1.0 - std::sqrt(peak), 1e-6) / 2.0);
        setGridGain(gridAp, 1.0 - std::abs(gen::ml::allpassG));
    }

    // On a mode switch: the new form's filter state from the shared history (x[-1], y[-1], y[-2]).
    void carryFilterState()
    {
        for (size_t j = 0; j < 4; ++j)
        {
            const auto& h = hist[j];
            auto& f = filt[j];
            f.s1 = f.b1 * h.x1 - f.a1 * h.y1 - f.a2 * h.y2; // transposed DF II: s1[n-1], s2[n-1]
            f.s2 = -f.a2 * h.y1;
            auto& q = filt16[j];
            q.s1 = q.p < 1.0 ? (h.y1 - q.p * h.y2) / (1.0 - q.p) : 0.0; // the shelf's output one sample back
        }
    }

    struct Biquad
    {
        double b0 = 0.0, b1 = 0.0, a1 = 0.0, a2 = 0.0, s1 = 0.0, s2 = 0.0;
    };

    bool isL = false;
    double lfX = 1.0, rtS = 1.0, hfX = 0.9, clusterGain = 80.0 / 90.0;
    Coefs coefs;
    std::array<Biquad, 4> filt{};
    std::array<ShelfPole16, 4> filt16{}; // 16-bit mode's form of the same filters
    struct History
    {
        double x1 = 0.0, y1 = 0.0, y2 = 0.0; // a line filter's last input and last two outputs, in either mode
    };
    std::array<History, 4> hist{};
    static constexpr int gridAp = 0, gridLine = 1;

    std::vector<float> lineBuf;            // four line buffers of lineSize, u_j (recirculating: store())
    std::array<std::vector<float>, 3> apBuf; // the all-passes' internal signal (recirculating: store())
    std::vector<float> inHist;             // the raw input, for the direct taps
    std::array<std::vector<double>, 2> outHist; // the output kernel's input per channel
    uint32_t pos = 0;                      // write index; the buffers are powers of two, so it wraps with them

    bool modActive = true;
    double lfoSin = 0.0, lfoCos = 1.0, thetaNow = 0.0;
    const double lfoStepSin = std::sin(2.0 * juce::MathConstants<double>::pi * gen::ml::thetaRateHz / deviceRate);
    const double lfoStepCos = std::cos(2.0 * juce::MathConstants<double>::pi * gen::ml::thetaRateHz / deviceRate);
};
} // namespace rv::dsp
