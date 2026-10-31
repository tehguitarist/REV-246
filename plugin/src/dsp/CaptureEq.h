#pragma once

#include <juce_core/juce_core.h>

#include <array>
#include <cmath>
#include <complex>

namespace rv::dsp
{
// The optional "commercial capture" low-band stage (docs/model_s_results.md §27): CLA's and the Altiverb's IRs share
// a minimum-phase transfer below ~400 Hz (+3 dB at 40-60 Hz, -0.6 dB at 160 Hz, ripple to 400 Hz) that is in neither
// Tomkin's library nor the stock unit; possibly the manual's transformer option (EMT246_unit_facts.md §2b). Part
// of the TRANSFORMER stage (Transformer.h), on by default.
//
// The 13 numbers are CLA's capture transfer fitted at the device rate (analysis/unit_chain/cla13.json, "CLA13T";
// docs/model_s_results.md §37.1-37.3, adopted by the user's decision on 31 October 2026 over bq.json's 10-number
// "CLA" entry: 6-12 dB deeper nulls at 160-400 Hz and 2 dB at 400-800 Hz against about 1 dB worse at 80-160 Hz and
// 0.8-1.6 kHz, a footprint of at most 0.17 dB): a first-order high-pass, an RBJ low shelf and three RBJ peaking
// sections, the formulas of tools/unit_transfer.py (hp1_sos, lowshelf_sos, peak_sos). The sections are designed at the host rate from the same parameters; that matches
// bq_resp within 0.001 dB at 50 Hz-2 kHz from 22.05 to 192 kHz (tests/TransformerTests.cpp). Below ~50 Hz the fit
// is not defined by the data. Double precision: the corners sit at 40-380 Hz, far below the host rate.
class CaptureEq
{
  public:
    // cla13.json "CLA13T": hp fc, shelf (f0, gain dB, Q), peak (f0, gain dB, Q) x 3.
    static constexpr std::array<double, 13> cla{39.99999999999999,  52.824460722524165, 9.326356298633552,
                                                0.5990898451122239, 163.42759767700667, -1.0610196325227694,
                                                1.1137864203744967, 295.39939134402744, 0.8637236387197627,
                                                1.39871337943611,   377.78388410094027, -0.5656165432870692,
                                                2.801949733068588};

    void prepare(double sampleRate)
    {
        const auto& x = cla;
        sections[0] = highPass1(x[0], sampleRate);
        sections[1] = lowShelf(x[1], x[2], x[3], sampleRate);
        sections[2] = peak(x[4], x[5], x[6], sampleRate);
        sections[3] = peak(x[7], x[8], x[9], sampleRate);
        sections[4] = peak(x[10], x[11], x[12], sampleRate);
        reset();
    }

    void reset()
    {
        for (auto& ch : state)
            for (auto& s : ch)
                s = {};
    }

    // One sample of channel `ch` (0 or 1).
    float process(int ch, float in)
    {
        double v = in;
        auto& chState = state[(size_t)ch];
        for (size_t k = 0; k < sections.size(); ++k)
        {
            const auto& c = sections[k];
            auto& s = chState[k];
            const auto y = c.b0 * v + s.z1; // transposed direct form II
            s.z1 = c.b1 * v - c.a1 * y + s.z2;
            s.z2 = c.b2 * v - c.a2 * y;
            v = y;
        }
        return (float)v;
    }

    // The cascade's response at `hz`, from the coefficients (as bq_resp computes it).
    std::complex<double> response(double hz, double sampleRate) const
    {
        const auto z1 = std::polar(1.0, -2.0 * juce::MathConstants<double>::pi * hz / sampleRate);
        const auto z2 = z1 * z1;
        std::complex<double> h{1.0, 0.0};
        for (const auto& c : sections)
            h *= (c.b0 + c.b1 * z1 + c.b2 * z2) / (1.0 + c.a1 * z1 + c.a2 * z2);
        return h;
    }

  private:
    struct Coeffs
    {
        double b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0; // normalised, a0 = 1
    };
    struct State
    {
        double z1 = 0, z2 = 0;
    };

    static Coeffs normalise(double b0, double b1, double b2, double a0, double a1, double a2)
    {
        return {b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0};
    }

    // First-order high-pass (bilinear).
    static Coeffs highPass1(double fc, double fs)
    {
        const auto k = std::tan(juce::MathConstants<double>::pi * fc / fs);
        return normalise(1.0, -1.0, 0.0, 1.0 + k, k - 1.0, 0.0);
    }

    // RBJ low shelf.
    static Coeffs lowShelf(double f0, double gainDb, double q, double fs)
    {
        const auto A = std::pow(10.0, gainDb / 40.0), w = juce::MathConstants<double>::twoPi * f0 / fs;
        const auto al = std::sin(w) / (2.0 * q), c = std::cos(w), sA = 2.0 * std::sqrt(A) * al;
        return normalise(A * ((A + 1) - (A - 1) * c + sA), 2 * A * ((A - 1) - (A + 1) * c),
                         A * ((A + 1) - (A - 1) * c - sA), (A + 1) + (A - 1) * c + sA, -2 * ((A - 1) + (A + 1) * c),
                         (A + 1) + (A - 1) * c - sA);
    }

    // RBJ peaking EQ.
    static Coeffs peak(double f0, double gainDb, double q, double fs)
    {
        const auto A = std::pow(10.0, gainDb / 40.0), w = juce::MathConstants<double>::twoPi * f0 / fs;
        const auto al = std::sin(w) / (2.0 * q), c = std::cos(w);
        return normalise(1 + al * A, -2 * c, 1 - al * A, 1 + al / A, -2 * c, 1 - al / A);
    }

    std::array<Coeffs, 5> sections;
    std::array<std::array<State, 5>, 2> state;
};
} // namespace rv::dsp
