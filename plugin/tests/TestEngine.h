#pragma once

// Shared by the processor-level tests: one ReverbEngine on its own, set up as the processor sets up its engines at
// the default panel, so a test can compare the processor's wet path with it.

#include "dsp/AdcCalibration.h"
#include "dsp/ReverbEngine.h"
#include "params/Parameters.h"
#include "ui/Design.h"

#include <algorithm>
#include <array>
#include <vector>

namespace rv::test
{
// The engine's settings at the parameters' defaults (program S, the CLA reference levers), as the processor builds them.
inline dsp::Settings defaultSettings(bool sixteenBit = false)
{
    using namespace rv::params;
    dsp::Settings s;
    for (int j = 0; j < numJoys; ++j)
        s.steps[(size_t)j] = joyDefaultStep((Joy)j);
    s.lowpassKHz = joyNumber(Program::s, Joy::lowpass, s.steps[(size_t)Joy::lowpass]);
    s.delayMs = joyNumber(Program::s, Joy::delay, s.steps[(size_t)Joy::delay]);
    s.amplPercent = joyNumber(Program::s, Joy::ampl, s.steps[(size_t)Joy::ampl]);
    s.lfX = joyNumber(Program::s, Joy::lf, s.steps[(size_t)Joy::lf]);
    s.rt = joyNumber(Program::s, Joy::rt, s.steps[(size_t)Joy::rt]);
    s.hfX = joyNumber(Program::s, Joy::hf, s.steps[(size_t)Joy::hf]);
    s.sixteenBit = sixteenBit;
    return s;
}

// The engine's (L, R) for the mono input `in` at the default panel and CAL, in blocks of `block`.
inline std::array<std::vector<float>, 2> runEngine(const std::vector<float>& in, double rate = 48000.0, int block = 256)
{
    dsp::ReverbEngine e;
    e.prepare(rate, block);
    e.setFullScale(dsp::adcFullScaleForZeroVu(design::meterCalDbFs[design::defaultMeterCal]));
    e.setSettings(defaultSettings());
    e.reset();
    std::array<std::vector<float>, 2> out{std::vector<float>(in.size()), std::vector<float>(in.size())};
    for (size_t i = 0; i < in.size(); i += (size_t)block)
    {
        const auto n = (int)std::min((size_t)block, in.size() - i);
        e.process(in.data() + i, out[0].data() + i, out[1].data() + i, n);
    }
    return out;
}
} // namespace rv::test
