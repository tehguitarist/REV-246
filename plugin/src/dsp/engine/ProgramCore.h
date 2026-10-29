#pragma once

#include "params/Parameters.h"

#include <array>

namespace rv::dsp
{
// The unit's sample rate: 19.2 MHz / 728. Every program core runs here; all its delays are integer device samples.
inline constexpr double deviceRate = 19.2e6 / 728.0;

// The reverb's settings as the panel shows them, already converted from step indices to values.
struct Settings
{
    params::Program program = params::Program::s;
    double lowpassKHz = 12.0; // applied by the engine on the reverb path, not by the cores
    double delayMs = 10.0;    // a pure delay of the whole output, applied by the engine, not by the cores
    double amplPercent = 80.0; // N: span in seconds in `ampl`'s place, see params::joyValue
    double lfX = 1.0;          // N: EARLY RT in seconds
    double rt = 1.0;           // unused in N
    double hfX = 0.9;          // N: LATE RT in seconds
    bool sixteenBit = false;   // quantise the loop storage to 16 bits (ProgramCore::store)
    std::array<int, params::numJoys> steps{}; // the raw step indices, for cores that want the table index
};

// One program's algorithm at the device rate: mono in, (L, R) out, both outputs at the same device instants (L's
// half-sample lag is the engine's output stage, not the core's). The contract is in docs/ENGINE.md.
class ProgramCore
{
  public:
    virtual ~ProgramCore() = default;

    // Allocates every buffer the program can need at any setting. Not real-time.
    virtual void prepare() = 0;
    // Clears all state (delay lines, filters, modulation phase back to its start). Real-time safe.
    virtual void reset() = 0;
    // Recomputes the coefficients for new panel values. Real-time safe: no allocation, no locks; a few hundred
    // transcendental calls are fine. Takes effect from the next sample; the engine handles any crossfade.
    virtual void setSettings(const Settings& s) = 0;
    // Processes `n` device-rate samples. Real-time safe.
    virtual void process(const float* in, float* outL, float* outR, int n) = 0;

    // Samples by which the core's output is late against the reference IRs (N's branch read kernel needs a
    // centred look-ahead; every other core returns 0). The engine adds this to its latency.
    virtual int latency() const { return 0; }
    // Time variance (M / L's theta, E's moving reads). On by default; the null tests switch it off to compare
    // the static network with the reference IRs. Takes effect at the next reset().
    void setModulation(bool on) { modulationOn = on; }

  protected:
    // Every write into recirculating storage (a loop's or all-pass's delay line) goes through here, so the 16-bit
    // mode can quantise it in one place. Float mode passes the value through.
    float store(float x) const { return quantise ? quantiseTo16(x) : x; }
    // Placeholder for the 16-bit mode (docs/ENGINE.md, 16-bit): set from Settings::sixteenBit by each core.
    bool quantise = false;
    bool modulationOn = true;

  private:
    static float quantiseTo16(float x) { return x; } // the fixed-point scheme lands in phase 3
};
} // namespace rv::dsp
