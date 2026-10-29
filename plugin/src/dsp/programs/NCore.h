#pragma once

// Program N (Non-Lin) at the device rate, from analysis/nl/nl_plugin_spec.json (docs/NonLin_program.md 9-12) and the
// sample-by-sample reference tools/nl_spec.py simulate(). Mono in, (L, R) out:
//
//   front taps   three plain taps per channel on the input, gains on the EARLY law
//   branch       input -> AP(124) -> one sample -> the "branch line"; main taps on it (up to 27 per channel, SPAN
//                switches whole stages, EARLY sets one gain per stage) -> AP(180) L / AP(188) R
//   late         the branch line at the feed -> AP(314) -> three flat loops (1108 / 1254 / 1386) at RT = LATE ->
//                four signed reads per loop and channel. One late section for both channels: R's feed is 220
//                samples earlier, which is the same loops read 220 samples nearer
//   kernel       the branch read kernel (adopted, spec branch_read_kernel): a centred 77-tap FIR per channel on the
//                whole branch (main taps and late, not the front taps); its look-ahead is latency(), the front taps
//                carry the same delay. setBranchKernel(false) is for the null tests (the reference IRs omit it)
//
// Panel slots (Settings): amplPercent = SPAN in seconds (0-0.5), lfX = EARLY RT in seconds (-3...10, 10 = flat),
// hfX = LATE RT in seconds (0-2.0, 0 = no late section). rt is unused. All three laws accept any value in range.
//
// Precision: the recirculating storage (the four all-passes' and the three loops' delay lines, written through
// ProgramCore::store) and the branch line are float; every sum and coefficient is double. The reference IRs are
// float32 renders of a double model, so float storage sets the null floor (tests/dsp/NTests.cpp reports it).
// Cost per sample: 6 front + up to 54 main taps, 4 all-passes (2 MACs each), 3 loops, 24 reads, and with the kernel
// 154 MACs more: about 100 MACs without it, about 250 with it.

#include "dsp/engine/ProgramCore.h"
#include "dsp/generated/n_spec.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace rv::dsp
{
class NCore final : public ProgramCore
{
  public:
    static constexpr int numStages = (int)gen::n::stageEnd.size();
    static constexpr int numMainTaps = (int)gen::n::mainDelayL.size();
    static constexpr int numLoops = (int)gen::n::loopDelay.size();
    static constexpr int numReads = (int)gen::n::readDelay.size();
    static constexpr int readsPerLoop = numReads / numLoops;

    // The resolved law of one (SPAN, EARLY): gains in dB of the three front taps per channel and of each running
    // stage, re a main tap at EARLY 10, and the late section's gain (docs/NonLin_program.md 9.3, 10).
    struct EarlyLaw
    {
        int stages = 0;
        std::array<double, 3> frontDbL{}, frontDbR{};
        std::array<double, numStages> stageDb{};
        double lateDb = 0.0;
    };

    // Stages that run at SPAN (s). The six library steps give 0 / 2 / 4 / 7 / 10 / 14. Between them (H): a stage
    // runs while its last tap ends at least `spanMarginS` before SPAN (the spec's 15-33 ms; the middle of it).
    static int stagesForSpan(double span)
    {
        const double limit = (std::clamp(span, 0.0, 0.5) - spanMarginS) * deviceRate;
        int ns = 0;
        for (int j = 0; j < numStages; ++j)
        {
            if ((double)gen::n::mainDelayL[(size_t)gen::n::stageEnd[(size_t)j]] > limit)
                break;
            ns = j + 1;
        }
        return ns;
    }

    // The late feed's delay on the branch line for `ns` running stages. The library steps are the spec's table; between
    // them (H) it is one sample before the last running tap, as the table is one before a main tap.
    static int feedFor(int ns)
    {
        for (size_t i = 0; i < gen::n::libStages.size(); ++i)
            if (gen::n::libStages[i] == ns)
                return gen::n::libFeed[i];
        return gen::n::mainDelayL[(size_t)gen::n::stageEnd[(size_t)(ns - 1)]] - 1;
    }

    // The EARLY law (spec early_law). EARLY > 0: stage j at tau = 16 + 36 j ms, gain -60 tau / EARLY dB, front taps at
    // their tick; EARLY < 0: the last stage at SPAN + 12 ms, the others 36 ms apart before it, gain +60 (tau - tau_end)
    // / |EARLY|; EARLY 10: flat. |EARLY| is held at 1 ms so EARLY = 0 stays finite (the unit has no such step).
    static EarlyLaw earlyLaw(double span, double early)
    {
        using namespace gen::n;
        EarlyLaw l;
        span = std::clamp(span, 0.0, 0.5);
        early = std::clamp(early, -3.0, 10.0);
        l.stages = stagesForSpan(span);
        if (early >= 10.0)
            return l;
        const double mag = std::max(std::abs(early), 1e-3);
        const double stageS = stageTicks * tickS;
        if (early >= 0.0)
        {
            for (size_t i = 0; i < 3; ++i)
            {
                l.frontDbL[i] = -60.0 * frontTickL[i] * tickS / mag;
                l.frontDbR[i] = -60.0 * frontTickR[i] * tickS / mag;
            }
            for (int j = 0; j < l.stages; ++j)
                l.stageDb[(size_t)j] = -60.0 * (stageTick0 * tickS + j * stageS) / mag;
            l.lateDb = -60.0 * (stageTick0 * tickS + std::min(l.stages, numStages - 1) * stageS) / mag;
            return l;
        }
        const double tEnd = span + negativeEndTicks * tickS;
        for (size_t i = 0; i < 3; ++i)
        {
            l.frontDbL[i] = 60.0 * (frontTickL[i] * tickS - tEnd) / mag;
            l.frontDbR[i] = 60.0 * (frontTickR[i] * tickS - tEnd) / mag;
        }
        for (int j = 0; j < l.stages; ++j)
            l.stageDb[(size_t)j] = 60.0 * (-stageS * (l.stages - 1 - j)) / mag;
        return l;
    }

    // Jot gain of a flat loop of `p` samples for RT = `late` s (spec late.loops.form).
    static double loopGain(int p, double late) { return std::pow(10.0, -3.0 * p / (deviceRate * late)); }

    // Test-only: bypass the branch read kernel (the reference IRs are kernel-free). Off, latency() is 0. Call before
    // reset(): it moves the front taps' read position and so the output's delay.
    void setBranchKernel(bool on) { kernelOn = on; }
    bool branchKernel() const { return kernelOn; }

    void prepare() override
    {
        using namespace gen::n;
        xHist.allocate((size_t)(frontDelayL.back() + bkCentre));  // front taps + the kernel's look-ahead
        branch.allocate((size_t)(mainDelayL.back() + branchDelay)); // the longest read: the last main tap
        apIn.allocate((size_t)apInDelay);
        apOut[0].allocate((size_t)apOutDelayL);
        apOut[1].allocate((size_t)apOutDelayR);
        apLate.allocate((size_t)apLateDelay);
        for (size_t l = 0; l < (size_t)numLoops; ++l)
            loops[l].allocate((size_t)loopDelay[l]);
        for (auto& k : kernelHist)
            k.assign(kernelRing, 0.0);
        for (size_t ch = 0; ch < 2; ++ch)           // delta + strength (taps - delta)
        {
            const auto& t = ch == 0 ? bkTapsL : bkTapsR;
            for (size_t i = 0; i < t.size(); ++i)
                kernel[ch][i] = (i == (size_t)bkCentre ? 1.0 : 0.0) * (1.0 - bkStrength) + bkStrength * t[i];
        }
        apIn.d = apInDelay;
        apOut[0].d = apOutDelayL;
        apOut[1].d = apOutDelayR;
        apLate.d = apLateDelay;
        prepared = true;
        reset();
    }

    void reset() override
    {
        xHist.clear();
        branch.clear();
        apIn.clear();
        apOut[0].clear();
        apOut[1].clear();
        apLate.clear();
        for (auto& l : loops)
            l.clear();
        for (auto& k : kernelHist)
            std::fill(k.begin(), k.end(), 0.0);
        kernelPos = 0;
    }

    void setSettings(const Settings& s) override
    {
        using namespace gen::n;
        quantise = s.sixteenBit;
        const double span = std::clamp(s.amplPercent, 0.0, 0.5);
        const double late = std::clamp(s.hfX, 0.0, 2.0);
        const auto law = earlyLaw(span, s.lfX);
        nTaps = law.stages == 0 ? 0 : stageEnd[(size_t)(law.stages - 1)] + 1;

        for (size_t i = 0; i < 3; ++i)
        {
            frontDelay[0][i] = frontDelayL[i];
            frontDelay[1][i] = frontDelayR[i];
            frontGain[0][i] = frontAmpL[i] * std::pow(10.0, law.frontDbL[i] / 20.0);
            frontGain[1][i] = frontAmpR[i] * std::pow(10.0, law.frontDbR[i] / 20.0);
        }
        int stage = 0;
        for (int i = 0; i < nTaps; ++i)
        {
            if (i > stageEnd[(size_t)stage])
                ++stage;
            const double g = std::pow(10.0, law.stageDb[(size_t)stage] / 20.0);
            // the main taps read the branch line one sample after the table's position (branch_delay)
            tapDelay[0][(size_t)i] = mainDelayL[(size_t)i] + branchDelay;
            tapDelay[1][(size_t)i] = mainDelayR[(size_t)i] + branchDelay;
            tapGain[0][(size_t)i] = mainSignL[(size_t)i] * g;
            tapGain[1][(size_t)i] = mainSignR[(size_t)i] * g;
        }

        const bool lateNow = late > 0.0;
        if (lateWasOn && !lateNow)                  // gated: nothing may ring when it comes back
        {
            apLate.clear();
            for (auto& l : loops)
                l.clear();
        }
        lateWasOn = lateNow;
        lateOn = lateNow;
        if (lateNow)
        {
            lateFeed = feedFor(law.stages) + branchDelay;
            lateGain = std::pow(10.0, law.lateDb / 20.0);
            for (size_t l = 0; l < (size_t)numLoops; ++l)
                loopG[l] = loopGain(loopDelay[l], late);
            for (size_t r = 0; r < (size_t)numReads; ++r)
            {
                readPos[0][r] = readDelay[r];
                readPos[1][r] = readDelay[r] - lateRShift;    // R's feed is lateRShift earlier: the same loop, nearer
                readGain[0][r] = readAmp * readSignL[r];
                readGain[1][r] = readAmp * readSignR[r];
            }
        }
    }

    void process(const float* in, float* outL, float* outR, int n) override
    {
        using namespace gen::n;
        if (!prepared)
        {
            std::fill(outL, outL + n, 0.0f);
            std::fill(outR, outR + n, 0.0f);
            return;
        }
        const unsigned pad = kernelOn ? (unsigned)bkCentre : 0u;
        for (int i = 0; i < n; ++i)
        {
            const double x = in[i];
            xHist.push(in[i]);
            const double y124 = allpass(apIn, x);
            branch.push((float)y124);

            std::array<double, 2> y{};
            for (size_t ch = 0; ch < 2; ++ch)
            {
                double acc = 0.0;
                for (int k = 0; k < nTaps; ++k)
                    acc += tapGain[ch][(size_t)k] * branch.at((unsigned)tapDelay[ch][(size_t)k]);
                y[ch] = allpass(apOut[ch], acc);
            }
            if (lateOn)
            {
                const double v = lateGain * allpass(apLate, branch.at((unsigned)lateFeed));
                for (size_t l = 0; l < (size_t)numLoops; ++l)
                {
                    const double line = v + loopG[l] * loops[l].at((unsigned)loopDelay[l] - 1);
                    loops[l].push(store((float)line));
                }
                // the reads, loop-major as generated
                for (size_t ch = 0; ch < 2; ++ch)
                {
                    double acc = 0.0;
                    for (size_t r = 0; r < (size_t)numReads; ++r)
                        acc += readGain[ch][r] * loops[r / (size_t)readsPerLoop].at((unsigned)readPos[ch][r]);
                    y[ch] += acc;
                }
            }
            if (kernelOn)
            {
                kernelPos = (kernelPos + 1) & (kernelRing - 1);
                for (size_t ch = 0; ch < 2; ++ch)
                {
                    auto& h = kernelHist[ch];
                    h[kernelPos] = y[ch];
                    double acc = 0.0;
                    for (size_t k = 0; k < kernelLen; ++k)
                        acc += kernel[ch][k] * h[(kernelPos - k) & (kernelRing - 1)];
                    y[ch] = acc;
                }
            }
            double l = y[0], r = y[1];
            for (size_t k = 0; k < 3; ++k)
            {
                l += frontGain[0][k] * xHist.at((unsigned)frontDelay[0][k] + pad);
                r += frontGain[1][k] * xHist.at((unsigned)frontDelay[1][k] + pad);
            }
            outL[i] = (float)l;
            outR[i] = (float)r;
        }
    }

    int latency() const override { return kernelOn ? gen::n::bkCentre : 0; }

  private:
    static constexpr double spanMarginS = 0.024;    // H: between the spec's 15 and 33 ms (docs/NonLin_program.md 9.2)
    static constexpr size_t kernelLen = gen::n::bkTapsL.size();
    static constexpr size_t kernelRing = 128;       // power of two >= kernelLen

    // A power-of-two ring of float: at(0) is the newest sample, at(d) the one d writes earlier.
    struct Ring
    {
        std::vector<float> b;
        unsigned mask = 0, pos = 0;
        // Room for reads up to `maxDelay` samples back.
        void allocate(size_t maxDelay)
        {
            size_t n = 1;
            while (n < maxDelay + 1)
                n <<= 1;
            b.assign(n, 0.0f);
            mask = (unsigned)(n - 1);
            pos = 0;
        }
        void clear()
        {
            std::fill(b.begin(), b.end(), 0.0f);
            pos = 0;
        }
        float at(unsigned d) const { return b[(pos - d) & mask]; }
        void push(float x)
        {
            pos = (pos + 1) & mask;
            b[pos] = x;
        }
    };

    // The state of one Schroeder all-pass (-g + z^-D) / (1 - g z^-D) in its single-delay form: w[n] = x[n] + g w[n-D],
    // y[n] = -g w[n] + w[n-D]. The delay line is recirculating storage.
    struct AllPass : Ring
    {
        int d = 0;
    };

    double allpass(AllPass& a, double x)
    {
        const double wd = a.at((unsigned)a.d - 1);
        const float w = store((float)(x + gen::n::apG * wd));
        a.push(w);
        return -gen::n::apG * w + wd;
    }

    bool prepared = false, kernelOn = true, lateOn = false, lateWasOn = false;
    Ring xHist, branch;
    AllPass apIn, apLate;
    std::array<AllPass, 2> apOut;
    std::array<Ring, numLoops> loops;
    std::array<std::vector<double>, 2> kernelHist;
    unsigned kernelPos = 0;
    std::array<std::array<double, kernelLen>, 2> kernel{};

    // The coefficients setSettings() resolves (fixed-size: no allocation there)
    int nTaps = 0, lateFeed = 0;
    double lateGain = 0.0;
    std::array<std::array<int, 3>, 2> frontDelay{};
    std::array<std::array<double, 3>, 2> frontGain{};
    std::array<std::array<int, numMainTaps>, 2> tapDelay{};
    std::array<std::array<double, numMainTaps>, 2> tapGain{};
    std::array<double, numLoops> loopG{};
    std::array<std::array<int, numReads>, 2> readPos{};
    std::array<std::array<double, numReads>, 2> readGain{};
};
} // namespace rv::dsp
