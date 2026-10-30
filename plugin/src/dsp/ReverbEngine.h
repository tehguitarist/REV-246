#pragma once

#include "dsp/engine/AnalogueHiss.h"
#include "dsp/engine/DeviceBridge.h"
#include "dsp/engine/ProgramCore.h"
#include "dsp/generated/e_spec.h"
#include "dsp/generated/ml_spec.h"
#include "dsp/programs/ECore.h"
#include "dsp/programs/MLCore.h"
#include "dsp/programs/NCore.h"
#include "dsp/programs/SCore.h"
#include "params/Parameters.h"

#include <juce_audio_basics/juce_audio_basics.h>
#include <juce_core/juce_core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <cstdint>
#include <vector>

namespace rv::dsp
{
// LOW PASS (H, e_plugin_spec controls.LOW_PASS): a 2nd-order Butterworth low-pass at the device rate on the whole
// reverb output (the operation card draws it after the reflections), bilinear with the cutoff prewarped. 12 kHz is
// "open": an exact bypass. Switching between open and a cutoff fades the filter in or out over `fadeSamples` (the
// filter's state is cleared while it is bypassed, so it starts silent behind the fade); a cutoff change while it is
// in swaps the coefficients at once, as the unit's coefficient loads would.
class DeviceLowPass
{
  public:
    static constexpr double openKHz = 12.0;
    static constexpr int fadeSamples = 264; // 10 ms at the device rate

    void reset()
    {
        for (auto& s : state)
            s = {};
        weight = target * fadeSamples;
    }

    // Real-time safe (one tan()).
    void setCutoff(double kHz)
    {
        if (kHz >= openKHz - 1e-9)
        {
            target = 0;
            return;
        }
        target = 1;
        if (std::equal_to<double>{}(kHz, cutoffKHz))
            return;
        cutoffKHz = kHz;
        const double k = std::tan(juce::MathConstants<double>::pi * kHz * 1000.0 / deviceRate);
        const double sq2 = juce::MathConstants<double>::sqrt2;
        const double norm = 1.0 / (1.0 + sq2 * k + k * k);
        b0 = k * k * norm;
        a1 = 2.0 * (k * k - 1.0) * norm;
        a2 = (1.0 - sq2 * k + k * k) * norm;
    }

    bool isOpen() const { return target == 0 && weight == 0; }

    void process(float* l, float* r, int n)
    {
        if (isOpen())
            return; // exact bypass
        float* ch[2] = {l, r};
        for (int i = 0; i < n; ++i)
        {
            if (weight != target * fadeSamples)
                weight += target == 1 ? 1 : -1;
            const double w = (double)weight / fadeSamples;
            for (size_t c = 0; c < 2; ++c)
            {
                const double x = ch[c][i];
                // transposed direct form II; b1 = 2 b0, b2 = b0
                auto& s = state[c];
                const double y = b0 * x + s.z1;
                s.z1 = 2.0 * b0 * x - a1 * y + s.z2;
                s.z2 = b0 * x - a2 * y;
                ch[c][i] = (float)(x + w * (y - x));
            }
            if (weight == 0 && target == 0)
            {
                for (auto& s : state)
                    s = {};
                return; // faded out: the rest of the block is bypassed
            }
        }
    }

  private:
    struct State
    {
        double z1 = 0.0, z2 = 0.0;
    };
    std::array<State, 2> state{};
    double cutoffKHz = -1.0, b0 = 1.0, a1 = 0.0, a2 = 0.0;
    int target = 0, weight = 0; // target 0 = open, 1 = filtering; weight counts 0 ... fadeSamples
};

// DELAY: a pure delay of the whole reverb output in integer device samples (0-999 ms). A change crossfades from the
// old read position to the new one over `fadeSamples` (sin^2 / cos^2, which sums to 1 for the correlated low band):
// the unit's read pointer would jump, which in a running tail is a click per step, and a held DELAY lever steps every
// few milliseconds. Changes arriving during a fade are taken one at a time, the latest last.
class OutputDelay
{
  public:
    static constexpr int fadeSamples = 264; // 10 ms at the device rate

    void prepare(int maxDelay)
    {
        size = (int)juce::nextPowerOfTwo(maxDelay + 2);
        for (auto& b : buf)
            b.assign((size_t)size, 0.0f);
        reset();
    }

    void reset()
    {
        for (auto& b : buf)
            std::fill(b.begin(), b.end(), 0.0f);
        pos = 0;
        current = from = target;
        fade = 0;
    }

    void setDelay(int samples) { target = std::clamp(samples, 0, size - 2); }
    int delay() const { return target; }

    void process(float* l, float* r, int n)
    {
        const auto mask = size - 1;
        float* ch[2] = {l, r};
        for (int i = 0; i < n; ++i)
        {
            if (fade == 0 && current != target)
            {
                from = current;
                current = target;
                fade = fadeSamples;
            }
            for (size_t c = 0; c < 2; ++c)
                buf[c][(size_t)pos] = ch[c][i];
            for (size_t c = 0; c < 2; ++c)
            {
                const auto now = buf[c][(size_t)((pos - current) & mask)];
                if (fade == 0)
                    ch[c][i] = now;
                else
                {
                    const auto old = buf[c][(size_t)((pos - from) & mask)];
                    const auto g = gainIn(fadeSamples - fade + 1);
                    ch[c][i] = old + g * (now - old);
                }
            }
            if (fade > 0)
                --fade;
            pos = (pos + 1) & mask;
        }
    }

  private:
    static float gainIn(int k) // sin^2 over the fade, k = 1 ... fadeSamples
    {
        const auto s = std::sin(0.5 * juce::MathConstants<double>::pi * (double)k / (double)fadeSamples);
        return (float)(s * s);
    }

    std::array<std::vector<float>, 2> buf;
    int size = 1, pos = 0, current = 0, from = 0, target = 0, fade = 0;
};

// The EMT 246 model, one per mono input: host-rate mono in, host-rate (L, R) wet out. docs/ENGINE.md has the contract.
//
//   host in -> DeviceBridge input side (analogue input chain, A/D: clip in 16-bit mode, overload flag)
//           -> device rate: alignment delay -> program core (S / M / L / E / N; two during a program crossfade)
//              -> DELAY -> LOW PASS
//           -> DeviceBridge output side (correction, S&H lag, output filter, L 0.5037 device samples late) -> [16-bit: hiss]
//
// Timing: every program comes out on one common grid, so the reported latency is the same for all of them and a
// program switch never moves the timing. N's centred branch kernel looks 40 device samples ahead (NCore::latency), and
// E's R direct tap sits at k0 = 64 device samples where M's first cluster tap sits at 10 (its output kernel's centre);
// the spec puts them together at DELAY 0 (e_plugin_spec controls.onset, docs/EMT250_program.md 17.8, S). The core
// inputs are delayed so that each program's own time origin lands `alignDelay` device samples late: S / M / L by
// alignDelay, N by alignDelay - 40, E by alignDelay - 54, with alignDelay = 54. latency() is the bridge's plus
// alignDelay and the bridge's correction delay (17 device samples) in host samples, rounded (the fraction left, under
// half a host sample, is not compensated).
//
// Program switch: the new core is reset and its input faded in over `programFadeSamples` (30 ms) while the old one
// keeps running on the input and its output fades out (equal power: the two reverbs are uncorrelated). Fading the new
// core's input rather than its output keeps each of its taps' onsets smooth. A switch requested during a fade waits
// for it to end. Settings within a program reach the running core at once, as the unit's coefficient loads
// did (docs/ENGINE.md, Decisions). Real-time safe after prepare(): no allocation, no locks.
class ReverbEngine
{
  public:
    static constexpr int programFadeSamples = 791; // 30 ms at the device rate
    static constexpr double maxDelayMs = 999.0;

    ReverbEngine()
    {
        cores = {&sCore, &mCore, &lCore, &eCore, &nCore};
    }
    ReverbEngine(const ReverbEngine&) = delete; // `cores` points into this object
    ReverbEngine& operator=(const ReverbEngine&) = delete;

    // Not real-time: allocates every core for its worst case and builds the bridge's tables.
    void prepare(double sampleRate, int maxBlockSize)
    {
        bridge.prepare(sampleRate, std::max(1, maxBlockSize));
        for (auto* c : cores)
            c->prepare();
        alignDelay = 0;
        for (int p = 0; p < params::numPrograms; ++p)
            alignDelay = std::max(alignDelay, lateness((params::Program)p));
        const auto maxDev = bridge.maxDeviceBlock();
        inRing.size = (int)juce::nextPowerOfTwo(maxDev + alignDelay + 2);
        inRing.data.assign((size_t)inRing.size * 2, 0.0f);
        for (auto* b : {&fadeBuf[0], &fadeBuf[1], &fadeGain[0], &fadeGain[1], &fadeIn})
            b->assign((size_t)maxDev, 0.0f);
        outputDelay.prepare((int)std::ceil(maxDelayMs / 1000.0 * deviceRate) + 1);
        hostAlign = (int)std::lround((double)(alignDelay + DeviceBridge::correctionDelay()) * sampleRate / deviceRate);
        hiss.prepare(sampleRate);
        prepared = true;
        applyFullScale();
        reset();
    }

    // Clears all state (bridge, cores, delay, filters) and ends any program crossfade. Keeps the settings. Real-time
    // safe.
    void reset()
    {
        if (!prepared)
            return;
        bridge.reset();
        for (auto* c : cores)
            c->reset();
        std::fill(inRing.data.begin(), inRing.data.end(), 0.0f);
        inCount = 0;
        if (pendingSwitch) // a queued switch takes effect now: nothing is left to fade
        {
            settings = pending;
            pendingSwitch = false;
        }
        active = settings.program;
        fadeFrom = active;
        fadeLeft = 0;
        cores[(size_t)active]->setSettings(settings);
        applyEngineSettings(settings);
        outputDelay.reset();
        lowPass.reset();
        hiss.reset();
    }

    // New panel values. Within the running program they apply at once; a new program starts a crossfade (or waits
    // for the running one). DELAY, LOW PASS, the A/D clip and the hiss apply at once in every case. Real-time safe.
    void setSettings(const Settings& s)
    {
        applyEngineSettings(s);
        if (!prepared)
        {
            settings = s;
            active = fadeFrom = s.program;
            return;
        }
        if (fadeLeft > 0 && s.program != active)
        {
            pending = s; // the switch after this one
            pendingSwitch = true;
            return;
        }
        pendingSwitch = false;
        settings = s;
        if (s.program == active)
            cores[(size_t)active]->setSettings(s);
        else
            startSwitch(s);
    }

    // The A/D's full scale (a peak in the input domain, from CAL): the bridge's clip level, the cores' 16-bit storage
    // scale and the hiss level all follow it.
    void setFullScale(float peak)
    {
        if (std::equal_to<float>{}(peak, fullScale))
            return;
        fullScale = peak;
        applyFullScale();
    }

    // True if the A/D was out of range (clipped, or in float mode would have) since the last call.
    bool takeOverload() { return bridge.takeOverload(); }

    // Host samples by which the wet output (R; L is 0.5037 device samples later inside the bridge) lags the input,
    // the same for every program.
    int latency() const { return bridge.latency() + hostAlign; }
    // Device samples every program is delayed to its common grid (the largest program look-ahead).
    int alignmentDelay() const { return alignDelay; }
    params::Program activeProgram() const { return active; }
    bool isCrossfading() const { return fadeLeft > 0 || pendingSwitch; }

    // For offline checks against the reference and hardware IRs (tools/render/RenderIr.cpp), not for the plugin:
    // the cores' time variance switch (ProgramCore::setModulation; takes effect at the next reset()), and the
    // device-rate path alone (alignment, core, DELAY, LOW PASS; no bridge, no analogue chain, no hiss), in chunks of
    // at most maxDeviceBlock() samples.
    void setModulation(bool on)
    {
        for (auto* c : cores)
            c->setModulation(on);
    }
    int maxDeviceBlock() const { return bridge.maxDeviceBlock(); }
    void processDeviceRate(const float* devIn, float* devL, float* devR, int nDev) { processDevice(devIn, devL, devR, nDev); }

    // `in` is the mono input; writes `n` samples of the wet signal per channel. Real-time safe.
    void process(const float* in, float* outL, float* outR, int n)
    {
        juce::ScopedNoDenormals noDenormals;
        if (!prepared)
        {
            std::fill(outL, outL + n, 0.0f);
            std::fill(outR, outR + n, 0.0f);
            return;
        }
        bridge.process(in, outL, outR, n,
                       [this](const float* devIn, float* devL, float* devR, int nDev) { processDevice(devIn, devL, devR, nDev); });
        hiss.add(outL, outR, n, sixteenBit);
    }

    // The device samples by which program p's core output is late against the common time origin (see the class
    // comment): N's look-ahead, E's direct tap against M's first cluster tap.
    int lateness(params::Program p) const
    {
        switch (p)
        {
        case params::Program::e:
            return gen::e::k0 - (gen::ml::m::clusterRDelay[0] + gen::ml::m::kernelCentre);
        case params::Program::n:
            return nCore.latency();
        case params::Program::s:
        case params::Program::m:
        case params::Program::l:
        case params::Program::count:
        default:
            return 0;
        }
    }

  private:
    struct Ring // doubled, so any window up to `size` is contiguous
    {
        std::vector<float> data;
        int size = 1;
        void write(std::int64_t i, float v)
        {
            const auto k = (size_t)(i & (size - 1));
            data[k] = data[k + (size_t)size] = v;
        }
        const float* from(std::int64_t i) const { return data.data() + (i & (size - 1)); }
    };

    void applyEngineSettings(const Settings& s)
    {
        const auto ms = std::clamp(s.delayMs, 0.0, maxDelayMs);
        outputDelay.setDelay((int)std::lround(ms / 1000.0 * deviceRate));
        lowPass.setCutoff(s.lowpassKHz);
        bridge.setClipping(s.sixteenBit);
        sixteenBit = s.sixteenBit;
    }

    void applyFullScale()
    {
        bridge.setFullScale(fullScale);
        for (auto* c : cores)
            c->setStorageFullScale((double)fullScale);
        hiss.setFullScale((double)fullScale);
    }

    void startSwitch(const Settings& s)
    {
        fadeFrom = active;
        active = s.program;
        auto& core = *cores[(size_t)active];
        core.reset();
        core.setSettings(s);
        fadeLeft = programFadeSamples;
    }

    // The last nDev device inputs, delayed by program p's alignment (contiguous).
    const float* coreInput(params::Program p, int nDev) const
    {
        return inRing.from(inCount - nDev - (alignDelay - lateness(p)));
    }

    void processDevice(const float* devIn, float* devL, float* devR, int nDev)
    {
        if (nDev <= 0)
            return;
        if (fadeLeft == 0 && pendingSwitch)
        {
            pendingSwitch = false;
            settings = pending;
            if (pending.program == active)
                cores[(size_t)active]->setSettings(pending);
            else
                startSwitch(pending);
        }
        for (int k = 0; k < nDev; ++k)
            inRing.write(inCount + k, devIn[k]);
        inCount += nDev;

        if (fadeLeft == 0)
            cores[(size_t)active]->process(coreInput(active, nDev), devL, devR, nDev);
        else
        {
            // The new core's INPUT fades in (so every delayed copy inside it starts smoothly: a reset core fed an
            // abrupt signal would click at each of its taps), and its output builds up as it would from silence; the
            // old core keeps its input and its OUTPUT fades out. Equal power (sin / cos of a sin^2 ramp).
            auto* gIn = fadeGain[0].data();
            auto* gOut = fadeGain[1].data();
            auto* newIn = fadeIn.data();
            const auto* x = coreInput(active, nDev);
            const auto halfPi = 0.5 * juce::MathConstants<double>::pi;
            for (int k = 0; k < nDev; ++k)
            {
                const auto t = std::min(1.0, (double)(programFadeSamples - fadeLeft + 1 + k) / (double)programFadeSamples);
                const auto u = std::sin(halfPi * t);
                const auto v = halfPi * u * u; // sin^2: flat at both ends, so neither gain has a kink there
                gIn[k] = (float)std::sin(v);
                gOut[k] = (float)std::cos(v);
                newIn[k] = gIn[k] * x[k];
            }
            cores[(size_t)active]->process(newIn, devL, devR, nDev);
            auto* oldL = fadeBuf[0].data();
            auto* oldR = fadeBuf[1].data();
            cores[(size_t)fadeFrom]->process(coreInput(fadeFrom, nDev), oldL, oldR, nDev);
            for (int k = 0; k < nDev; ++k)
            {
                devL[k] += gOut[k] * oldL[k];
                devR[k] += gOut[k] * oldR[k];
            }
            fadeLeft = std::max(0, fadeLeft - nDev);
        }
        outputDelay.process(devL, devR, nDev);
        lowPass.process(devL, devR, nDev);
    }

    DeviceBridge bridge;
    SCore sCore;
    MLCore mCore{params::Program::m}, lCore{params::Program::l};
    ECore eCore;
    NCore nCore;
    std::array<ProgramCore*, params::numPrograms> cores{};
    AnalogueHiss hiss;

    Settings settings, pending;
    params::Program active = params::Program::s, fadeFrom = params::Program::s;
    bool pendingSwitch = false, prepared = false, sixteenBit = false; // sixteenBit: the A/D clip and hiss, at once
    int fadeLeft = 0, alignDelay = 0, hostAlign = 0;
    float fullScale = 1.0f;

    Ring inRing;
    std::int64_t inCount = 0;
    std::array<std::vector<float>, 2> fadeBuf, fadeGain; // the old core's output; the fade's gains (in, out)
    std::vector<float> fadeIn;                           // the new core's faded-in input
    OutputDelay outputDelay;
    DeviceLowPass lowPass;
};
} // namespace rv::dsp
