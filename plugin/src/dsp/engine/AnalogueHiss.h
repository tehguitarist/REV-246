#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace rv::dsp
{
// The unit's analogue / converter hiss on the wet outputs, at the host rate (docs/ENGINE.md, 16-bit mode). The engine
// adds it with on = Settings::sixteenBit; FLOAT stays clean.
//
// Level: the data sheet's pass-through (DELAY program) S/N is 88 dB unweighted rms re nominal (C, unit facts §2); the
// plugin uses 91 dB (user choice, 30 Oct 2026: a better-than-spec unit). Nominal is 6 dB below the A/D's full scale
// (C, unit facts §2: the 0 dB LED sits 6 dB under the REG overload) and is read as a sine (H: the spec's reference
// signal is not stated), so the hiss is 91 + 6 = 97 dB below a full-scale sine's rms, fs / sqrt(2) x 10^(-97/20),
// with fs the A/D full scale (a peak) in host units from CAL. The REVERB program's 73 dB is not calibrated here: its
// extra 15 dB is the 16-bit storage's own noise, which the cores' 16-bit mode makes by itself
// (tests/dsp/SixteenBitTests.cpp reports the total against it).
//
// Spectrum (H): white, through a second-order Butterworth low-pass at 12 kHz (the output low-pass's band, unit facts
// §2: 12 kHz in -> -2.5 dB), calibrated so its power in 20 Hz-20 kHz (up to Nyquist at low host rates) is the level
// above (an unweighted measurement band). L and R are independent (H: the output stages' noise; the A/D's share would
// be common to both). Deterministic: xorshift32 per channel, seeds restored by reset(). Real-time safe, no allocation;
// cost per sample and channel: one xorshift, one biquad.
class AnalogueHiss
{
  public:
    static constexpr double snrDb = 91.0;     // re nominal (data sheet 88 dB, C; plugin 91 dB, user choice 30 Oct)
    static constexpr double nominalDb = -6.0; // nominal re full scale (C)
    static constexpr double cutoffHz = 12000.0;

    void prepare(double hostRate)
    {
        rate = hostRate > 0.0 ? hostRate : 48000.0;
        // Butterworth low-pass by the bilinear transform with prewarping, cut-off held under 0.45 x the host rate.
        const double pi = std::acos(-1.0);
        const double fc = std::min(cutoffHz, 0.45 * rate);
        const double k = std::tan(pi * fc / rate), q = std::sqrt(0.5);
        const double norm = 1.0 / (1.0 + k / q + k * k);
        b0 = k * k * norm;
        b1 = 2.0 * b0;
        b2 = b0;
        a1 = 2.0 * (k * k - 1.0) * norm;
        a2 = (1.0 - k / q + k * k) * norm;
        // Power of the filtered unit-variance white noise in the measurement band: (2 / rate) x integral |H|^2 df.
        const double f0 = 20.0, f1 = std::min(20000.0, 0.5 * rate);
        constexpr int steps = 2000;
        double band = 0.0;
        for (int i = 0; i <= steps; ++i)
        {
            const double f = f0 + (f1 - f0) * i / steps;
            const double w = 2.0 * pi * f / rate;
            const double c1 = std::cos(w), c2 = std::cos(2.0 * w), s1 = std::sin(w), s2 = std::sin(2.0 * w);
            const double nr = b0 + b1 * c1 + b2 * c2, ni = -(b1 * s1 + b2 * s2);
            const double dr = 1.0 + a1 * c1 + a2 * c2, di = -(a1 * s1 + a2 * s2);
            const double h2 = (nr * nr + ni * ni) / (dr * dr + di * di);
            band += (i == 0 || i == steps ? 0.5 : 1.0) * h2;
        }
        band *= (f1 - f0) / steps * 2.0 / rate;
        // uniform in [-1, 1) has variance 1/3: this gain makes the band rms equal to unitRms (per unit full scale)
        unitGain = unitRms() / std::sqrt(band / 3.0);
        setFullScale(fullScale);
        reset();
    }

    void reset()
    {
        seed = {seedL, seedR};
        for (auto& s : state)
            s = {};
    }

    // The A/D full scale (a peak) in host units, from CAL. Real-time safe.
    void setFullScale(double fs)
    {
        fullScale = fs > 0.0 ? fs : 1.0;
        gain = unitGain * fullScale;
    }

    // Adds the hiss to the wet outputs when `on`; otherwise leaves them (and the generator) untouched.
    void add(float* l, float* r, int n, bool on)
    {
        if (!on)
            return;
        float* out[2] = {l, r};
        for (int ch = 0; ch < 2; ++ch)
        {
            if (out[ch] == nullptr)
                continue;
            auto& s = state[(size_t)ch];
            auto x = seed[(size_t)ch];
            for (int i = 0; i < n; ++i)
            {
                x ^= x << 13;
                x ^= x >> 17;
                x ^= x << 5;
                const double u = (double)((std::int64_t)x - (x >> 31 ? 4294967296LL : 0)) * (1.0 / 2147483648.0); // x as two's complement: [-1, 1)
                const double y = b0 * u + s.z1;                            // transposed direct form II
                s.z1 = b1 * u - a1 * y + s.z2;
                s.z2 = b2 * u - a2 * y;
                out[ch][i] += (float)(gain * y);
            }
            seed[(size_t)ch] = x;
        }
    }

    // The target rms in the measurement band for a full scale of 1.0.
    static double unitRms() { return std::pow(10.0, (nominalDb - snrDb) / 20.0) / std::sqrt(2.0); }
    // The current target rms in host units.
    double targetRms() const { return unitRms() * fullScale; }

  private:
    struct State
    {
        double z1 = 0.0, z2 = 0.0;
    };
    static constexpr uint32_t seedL = 0x9E3779B9u, seedR = 0x85EBCA6Bu;

    double rate = 48000.0, fullScale = 1.0, unitGain = 0.0, gain = 0.0;
    double b0 = 0.0, b1 = 0.0, b2 = 0.0, a1 = 0.0, a2 = 0.0;
    std::array<uint32_t, 2> seed{seedL, seedR};
    std::array<State, 2> state{};
};
} // namespace rv::dsp
