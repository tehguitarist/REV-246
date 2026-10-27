#pragma once

#include "dsp/CaptureEq.h"

#include <juce_core/juce_core.h>

#include <array>
#include <cmath>

namespace rv::dsp
{
// The core of one audio transformer, as a flux model. The voltage across the winding is integrated to the core flux;
// the magnetising current's nonlinear part, dropped across the source and winding resistance, is the distortion
// voltage, which is subtracted from the signal. Odd only (no DC reaches the core: the drivers are capacitor coupled),
// mostly 3rd harmonic, and it grows with level and falls steeply with frequency, as audio transformers do.
//
// GUESSED (H): the service manual (pp. 37, 46) offers input and output transformers only as "slightly" raising the
// distortion, with no type or values. The two terms are a Rayleigh-like a·φ|φ| (low and mid level) and a saturation
// onset b·φ³; they are set so that at a peak of 1 (= +12 dBu: the unit's +6 dBu nominal plus its 6 dB headroom) and 40 Hz
// the THD is 0.08 %, three quarters from the Rayleigh term. Against datasheets: 0.11 % at 20 Hz / +4 dBu (Jensen
// JT-11-FLCF, a modern 50 % nickel part: 0.035 %), 0.064 % at 50 Hz / +14 dBu (Lundahl LL1935: < 0.1 %), 0.31 % at
// 50 Hz / +24 dBu (LL1935: < 1 % at +26). 0.011 % at 100 Hz, 0.0001 % at 1 kHz; 5th harmonic 20 dB under the 3rd.
class TransformerCore
{
  public:
    // `drive` scales the distortion (the source-plus-winding resistance relative to the output transformer's).
    void prepare(double sampleRate, double drive)
    {
        const auto twoPi = juce::MathConstants<double>::twoPi;
        fluxGain = twoPi * refHz / sampleRate;
        pole = std::exp(-twoPi * cornerHz / sampleRate);
        blockGain = 0.5 * (1.0 + pole);
        rayleighGain = drive * rayleigh;
        cubicGain = drive * saturation;
        reset();
    }

    void reset() { blockerIn = blockerOut = flux = 0.0; }

    float process(float in)
    {
        const double v = in;
        blockerOut = blockGain * (v - blockerIn) + pole * blockerOut; // DC blocker at cornerHz (the coupling capacitor)
        blockerIn = v;
        flux = pole * flux + fluxGain * blockerOut; // a leaky integral, |flux| = |v| for a sine at refHz
        const auto p = juce::jlimit(-maxFlux, maxFlux, flux);
        return (float)(v - (rayleighGain * p * std::abs(p) + cubicGain * p * p * p));
    }

  private:
    static constexpr double refHz = 40.0, cornerHz = 5.0;
    static constexpr double maxFlux = 8.0; // bounds the cubic term for a huge input (flux 1 = the full-scale sine at refHz)
    static constexpr double rayleigh = 0.00353429, saturation = 0.0008; // THD 0.08 % at 40 Hz, peak 1 (see above)
    double fluxGain = 0, pole = 0, blockGain = 0, rayleighGain = 0, cubicGain = 0;
    double blockerIn = 0, blockerOut = 0, flux = 0;
};

// The transformer-fitted unit (the manual's option, H): an input transformer on the mono input, and on each output
// the measured low-band transfer of the commercial captures (CaptureEq) followed by an output transformer.
//
// Level: the cores are calibrated with a peak of 1 = +12 dBu, the A/D's full scale (nominal +6 dBu plus its 6 dB
// headroom; the unit is unity gain in and out, so the outputs share the reference). The A/D's full scale follows the
// meter calibration (AdcClipper::fullScaleForZeroVu), so setFullScale moves the cores' drive with it.
class Transformer
{
  public:
    void prepare(double sampleRate)
    {
        eq.prepare(sampleRate);
        input.prepare(sampleRate, 0.5); // bridging input driven from a low-impedance studio source: half the drive
        for (auto& o : output)
            o.prepare(sampleRate, 1.0);
    }

    void reset()
    {
        eq.reset();
        input.reset();
        for (auto& o : output)
            o.reset();
    }

    // The A/D's full-scale peak in host units: the level the cores treat as +12 dBu.
    void setFullScale(float peak)
    {
        fullScale = peak;
        invFullScale = 1.0f / peak;
    }

    float processInput(float in) { return input.process(finiteOrZero(in) * invFullScale) * fullScale; }
    float processOutput(int ch, float in)
    {
        return output[(size_t)ch].process(eq.process(ch, finiteOrZero(in)) * invFullScale) * fullScale;
    }

    const CaptureEq& captureEq() const { return eq; }

  private:
    // A NaN/inf would stay in the recursive filters' state for good.
    static float finiteOrZero(float x) { return std::isfinite(x) ? x : 0.0f; }

    CaptureEq eq;
    TransformerCore input;
    std::array<TransformerCore, 2> output;
    float fullScale = 1.0f, invFullScale = 1.0f;
};
} // namespace rv::dsp
