#pragma once

#include <juce_core/juce_core.h>

#include <cmath>

namespace rv::dsp
{
// The A/D's full scale against the meter calibration (docs/EMT246_unit_facts.md section 2). The unit's nominal level
// (its 0 dB LED) is 6 dB below the A/D's overload point for a sine below 2 kHz (C). CAL says what 0 VU is in dBFS, so
// the full scale (a peak in the input domain, unity gain at 1 kHz) follows it: DeviceBridge::setFullScale, the cores'
// 16-bit storage and the transformer cores all take this number.
inline constexpr double adcHeadroomDb = 6.0;

// The full-scale peak for a meter whose 0 VU reads `zeroVuDbFs` (RMS): a sine 6 dB over 0 VU just fits.
inline float adcFullScaleForZeroVu(float zeroVuDbFs)
{
    return juce::MathConstants<float>::sqrt2 * std::pow(10.0f, (zeroVuDbFs + (float)adcHeadroomDb) / 20.0f);
}
} // namespace rv::dsp
