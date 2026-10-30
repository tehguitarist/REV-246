#pragma once

#include "dsp/engine/BridgeFilter.h"
#include "dsp/engine/ProgramCore.h"
#include "dsp/generated/bridge_spec.h"

#include <juce_core/juce_core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

namespace rv::dsp
{
// Host rate <-> device rate, carrying the unit's analogue chain (docs/ENGINE.md, Signal path; docs/NEXT_STEPS.md section 4;
// unit facts section 2b). One per engine:
//
//   host mono -> input filter 7 246 131 x coupling high-passes 3.35, 0.23, 0.72 Hz -> sampled at the device instants
//             -> [A/D clip] -> the program (device rate, the caller's) -> (L, R)
//             -> per channel: correction (1.8 Hz high-pass, FIR) -> output S&H lag (device-rate one-pole)
//             -> output filter 7 246 130 x 0.16 Hz high-pass -> host rate; L 0.5037 device samples late.
//
// Design: each analogue filter is discretised exactly for band-limited signals (BridgeFilter): the host signal is the
// band-limited one its samples stand for (a Kaiser sinc flat to 1e-6 up to 12.8 kHz), and the input filter's response to
// it is read at the device instants; the device samples drive the output filter as an impulse train (the analogue model
// has no zero-order hold: p absorbs it, NEXT_STEPS section 4), whose response is read at the host instants after the same
// sinc as the capture's anti-alias low-pass. No anti-alias filter is added at the device rate: whatever passes the input
// filter above 13.19 kHz folds, as in the unit. Short tables for the last few samples, one recursion per pole for the rest,
// so nothing is truncated; double throughout. The sampling grid is exact rational arithmetic (host rate in mHz against
// 19.2 MHz / 728), so the device instants never drift and any block split gives bit-identical output.
//
// Correction (analysis/plugin_check/REPORT_chainfix.md): the netlists miss the measured chain (chain3) by up to 70 degrees
// in the output filter's skirt above 8 kHz and by about a degree in the mid band. Their ratio's unit part runs per channel
// on the device samples before the S&H lag: a first-order high-pass at 1.8 Hz (H: Tomkin's 6.25 Hz total over the
// schematic's 4.46 Hz, a DC blocker in the firmware or the capture rig; analog_results.md section 3) and a 65-tap FIR with
// a bulk delay of correctionDelay() device samples, which the engine adds to its latency. The capture rig's own residue
// below 80 Hz is left out.
//
// Gain: unity at 1 kHz on each side (the input filter's inversion is not modelled); the correction is unity at 1 kHz on
// R, L keeps its measured level (+0.1 dB) and shape re R. Full scale for the clip is a peak in that input domain (the
// A/D's, AdcCalibration.h). Latency: latency() host samples on R plus correctionDelay() device samples (L another
// lLateDeviceSamples).
class DeviceBridge
{
  public:
    // The output S&H's charging lag y[n] = y[n-1] + (1 - p)(x[n] - y[n-1]) at the device rate (unit facts 2b): p = 0.301
    // fits Tomkin's unit (the library this plugin is checked against); CLA's unit fits 0.325. A per-unit trim (P202).
    static constexpr double shPole = 0.301;
    // L's S&H pulse comes half a sample period after R's (the timing PROM; measured 0.49-0.51, unit facts 2b); against
    // the library, 0.5037 (plugin_check REPORT section 4: a fixed L / R chain-delay difference folded in).
    static constexpr double lLateDeviceSamples = gen::bridge::lLateDeviceSamples;
    // The correction's FIR length and bulk delay (device samples; the engine compensates it).
    static constexpr int correctionTaps = (int)gen::bridge::fixFirR.size();
    static constexpr int correctionDelay() { return gen::bridge::fixDelayDeviceSamples; }

    // Not real-time: builds the tables (a few ms) and allocates for blocks of up to maxBlock host samples (longer blocks
    // are split).
    void prepare(double hostRate, int maxBlock)
    {
        namespace spec = gen::bridge;
        hostRateHz = hostRate;
        chunk = std::max(1, maxBlock);
        // The sampling grid: host samples per device sample = num / den exactly.
        num = (std::int64_t)std::llround(hostRate * 1000.0) * 728;
        den = (std::int64_t)19200000 * 1000;
        const auto r = (double)num / (double)den;
        halfNum = (std::int64_t)std::llround(lLateDeviceSamples * (double)num);

        // phi: a Kaiser sinc at host Nyquist whose transition is symmetric about it, from the 12.8 kHz band edge to the
        // image of 12.8 kHz, 120 dB design attenuation (pass ripple about 1e-6 at every rate).
        const auto hostT = 1.0 / hostRate;
        const auto edge = std::min(passEdgeHz, 0.4 * hostRate);
        const auto dw = juce::MathConstants<double>::twoPi * (1.0 - 2.0 * edge / hostRate);
        const auto halfWidth = std::max(minHalfWidth, (designDb - 8.0) / (2.285 * dw) / 2.0); // host samples
        const auto beta = 0.1102 * (designDb - 8.7);
        const auto phiHalf = halfWidth * hostT;

        std::vector<BridgeFilter::Mode> in, out;
        for (size_t j = 0; j < spec::inPoleRe.size(); ++j)
            in.push_back({{spec::inPoleRe[j], spec::inPoleIm[j]}, {spec::inResRe[j], spec::inResIm[j]}, spec::inWeight[j]});
        for (size_t j = 0; j < spec::outPoleRe.size(); ++j)
            out.push_back({{spec::outPoleRe[j], spec::outPoleIm[j]}, {spec::outResRe[j], spec::outResIm[j]}, spec::outWeight[j]});

        const auto devT = 1.0 / deviceRate;
        // The S&H lag's gain at 1 kHz is taken out on the output side, so the chain is unity there.
        const auto w1k = juce::MathConstants<double>::twoPi * 1000.0 * devT;
        const auto lag1k = (1.0 - shPole) / std::abs(std::complex<double>(1.0 - shPole * std::cos(w1k), shPole * std::sin(w1k)));
        // The correction's high-pass: bilinear, corner prewarped (as chain_fix.json's fit).
        const auto k = std::tan(juce::MathConstants<double>::pi * spec::fixHighPassHz * devT);
        hpGain = 1.0 / (1.0 + k);
        hpPole = (1.0 - k) / (1.0 + k);
        inFilter.build(in, hostT, hostT, phiHalf, beta, std::max(64, (int)std::ceil(hostT / tableStep)), 1.0);
        outFilter.build(out, devT, hostT, phiHalf, beta, std::max(64, (int)std::ceil(devT / tableStep)), 1.0 / lag1k);
        wi = inFilter.halfTaps();
        wo = outFilter.halfTaps();
        // The output at host m reads device samples up to floor((m - delay) / r) + wo; the input side has made those once
        // host sample floor(n r) + wi has arrived. So delay >= wi + wo r.
        delay = wi + (int)std::ceil((double)wo * r - 1e-9);

        xRing.size = juce::nextPowerOfTwo(2 * wi + (int)std::ceil(r) + 8);
        xRing.data.assign((size_t)xRing.size * 2, 0.0);
        maxDevice = (int)std::ceil(chunk / r) + 2;
        const auto ySize = juce::nextPowerOfTwo((int)std::ceil((chunk + delay) / r) + 2 * wo + 16);
        for (auto* ring : {&yRing[0], &yRing[1]})
        {
            ring->size = ySize;
            ring->data.assign((size_t)ySize * 2, 0.0);
        }
        devIn.assign((size_t)maxDevice, 0.0f);
        devOut[0].assign((size_t)maxDevice, 0.0f);
        devOut[1].assign((size_t)maxDevice, 0.0f);
        reset();
    }

    // Clears all state; the device grid restarts with device instant 0 at host sample 0. Real-time safe.
    void reset()
    {
        for (auto* ring : {&xRing, &yRing[0], &yRing[1]})
            std::fill(ring->data.begin(), ring->data.end(), 0.0);
        inState.clear();
        outState[0].clear();
        outState[1].clear();
        lag[0] = lag[1] = 0.0;
        for (auto& c : fix)
            c = {};
        hostIn = 0;
        devHost = 0;
        devRem = 0;
        tailIn = -wi; // inputs at m <= -wi are all zero
        devPushed = 0;
        // The first output (host 0) is analogue time -delay: device position -delay den / num.
        const auto p0 = -(std::int64_t)delay * den;
        outJ = floorDiv(p0, num);
        outRem = p0 - outJ * num;
        tailOut[0] = outJ - wo - 2; // R; zeros before device 0
        tailOut[1] = outJ - wo - 2; // L
        overload = false;
    }

    // Host samples by which R lags the modelled analogue chain (L: plus lLateDeviceSamples device samples).
    int latency() const { return delay; }
    // The most device samples one call of the device callback can get.
    int maxDeviceBlock() const { return maxDevice; }
    double hostRate() const { return hostRateHz; }

    // The A/D (unit facts section 2): full scale is a peak in the input domain (unity at 1 kHz); two's complement, so one
    // code short at the top. clip = the 16-bit mode's hard clip; off (float) only reports where it would clip.
    void setFullScale(float peak)
    {
        fsLow = -(double)peak;
        fsHigh = (double)peak * 32767.0 / 32768.0;
    }
    void setClipping(bool on) { clip = on; }
    // True if an A/D sample was outside full scale (clipped or, in float mode, would have been) since the last call.
    bool takeOverload()
    {
        const auto o = overload;
        overload = false;
        return o;
    }

    // Processes n host samples. `device(const float* in, float* outL, float* outR, int numDevice)` runs the program on the
    // device samples this block produced (numDevice <= maxDeviceBlock(), possibly 0). Real-time safe.
    template <typename DeviceProcess>
    void process(const float* in, float* outL, float* outR, int n, DeviceProcess&& device)
    {
        for (int start = 0; start < n; start += chunk)
        {
            const auto len = std::min(chunk, n - start);
            const auto numDevice = inputPass(in + start, len);
            device(devIn.data(), devOut[0].data(), devOut[1].data(), numDevice);
            outputPass(numDevice, outL + start, outR + start, len);
        }
    }

  private:
    struct Ring // doubled, so any window shorter than size is contiguous
    {
        std::vector<double> data;
        int size = 1;
        void write(std::int64_t i, double v)
        {
            const auto k = (size_t)(i & (size - 1));
            data[k] = data[k + (size_t)size] = v;
        }
        double at(std::int64_t i) const { return data[(size_t)(i & (size - 1))]; }
        const double* from(std::int64_t i) const { return data.data() + (i & (size - 1)); }
    };

    static constexpr double passEdgeHz = 12.8e3, designDb = 120.0, minHalfWidth = 6.0;
    static constexpr double tableStep = 5e-8; // head table phase step (s): interpolation error below -115 dB

    static std::int64_t floorDiv(std::int64_t a, std::int64_t b)
    {
        const auto q = a / b;
        return (a % b != 0 && (a < 0) != (b < 0)) ? q - 1 : q;
    }

    int inputPass(const float* in, int len)
    {
        int count = 0;
        for (int i = 0; i < len; ++i)
        {
            const auto x = std::isfinite(in[i]) ? (double)in[i] : 0.0; // a NaN would stay in the recursions for good
            xRing.write(hostIn, x);
            // Every device instant whose head window is now complete (device n at host position devHost + u).
            while (devHost + wi <= hostIn)
            {
                while (tailIn < devHost - wi)
                    inFilter.push(inState, xRing.at(++tailIn));
                const auto u = (double)devRem / (double)den;
                auto v = inFilter.head(xRing.from(devHost - wi + 1), u) + inFilter.tail(inState, u);
                if (v < fsLow || v > fsHigh)
                {
                    overload = true;
                    if (clip)
                        v = std::clamp(v, fsLow, fsHigh);
                }
                devIn[(size_t)count++] = (float)v;
                devRem += num;
                devHost += devRem / den;
                devRem %= den;
            }
            ++hostIn;
        }
        jassert(count <= maxDevice);
        return count;
    }

    void outputPass(int numDevice, float* outL, float* outR, int len)
    {
        for (int k = 0; k < numDevice; ++k)
        {
            for (int c = 0; c < 2; ++c)
            {
                const auto x = devOut[(size_t)c][(size_t)k];
                const auto v = correct(c, std::isfinite(x) ? (double)x : 0.0);
                lag[c] += (1.0 - shPole) * (v - lag[c]);
                yRing[c].write(devPushed, lag[c]);
            }
            ++devPushed;
        }
        for (int i = 0; i < len; ++i)
        {
            // R at device position outJ + outRem / num; L half a device sample earlier in the sequence (= later output).
            outR[i] = (float)readOut(1, outJ, (double)outRem / (double)num);
            auto remL = outRem - halfNum;
            auto jL = outJ;
            if (remL < 0)
            {
                remL += num;
                --jL;
            }
            outL[i] = (float)readOut(0, jL, (double)remL / (double)num);
            outRem += den;
            outJ += outRem / num;
            outRem %= num;
        }
    }

    // The correction on channel c's next device sample: high-pass, then the FIR (newest sample first in the ring).
    double correct(int c, double x)
    {
        auto& s = fix[(size_t)c];
        const auto y = hpGain * (x - s.hpIn) + hpPole * s.hpOut;
        s.hpIn = x;
        s.hpOut = y;
        s.pos = (s.pos + fixRing - 1) & (fixRing - 1);
        s.ring[(size_t)s.pos] = s.ring[(size_t)(s.pos + fixRing)] = y;
        const auto& h = c == 0 ? gen::bridge::fixFirL : gen::bridge::fixFirR;
        const auto* r = s.ring.data() + s.pos;
        double acc = 0.0;
        for (int i = 0; i < correctionTaps; ++i)
            acc += h[(size_t)i] * r[i];
        return acc;
    }

    // Channel c (0 = L, 1 = R) at device position j + v.
    double readOut(int c, std::int64_t j, double v)
    {
        jassert(j + wo < devPushed); // latency covers the head's look-ahead
        auto& ring = yRing[c];
        auto& st = outState[c];
        while (tailOut[c] < j - wo)
            outFilter.push(st, ring.at(++tailOut[c]));
        return outFilter.head(ring.from(j - wo + 1), v) + outFilter.tail(st, v);
    }

    static constexpr int fixRing = 128; // a power of two >= correctionTaps; doubled so the FIR reads contiguously
    static_assert(fixRing >= correctionTaps && gen::bridge::fixFirL.size() == gen::bridge::fixFirR.size());
    struct Fix
    {
        double hpIn = 0.0, hpOut = 0.0;
        std::array<double, 2 * fixRing> ring{};
        int pos = 0;
    };

    BridgeFilter inFilter, outFilter;
    std::array<Fix, 2> fix{}; // [0] = L, [1] = R
    double hpGain = 1.0, hpPole = 0.0;
    BridgeFilter::State inState, outState[2];
    Ring xRing, yRing[2]; // yRing[0] = L, [1] = R
    std::vector<float> devIn, devOut[2];
    double hostRateHz = 48000.0, lag[2]{};
    double fsLow = -1.0, fsHigh = 1.0;
    std::int64_t num = 1, den = 1, halfNum = 0;
    std::int64_t hostIn = 0, devHost = 0, devRem = 0, tailIn = 0; // input side: next device instant at devHost + devRem / den
    std::int64_t devPushed = 0, outJ = 0, outRem = 0, tailOut[2]{};
    int chunk = 512, maxDevice = 1, wi = 1, wo = 1, delay = 0;
    bool clip = false, overload = false;
};
} // namespace rv::dsp
