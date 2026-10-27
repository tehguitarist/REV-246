#pragma once

#include <juce_core/juce_core.h>

#include <array>
#include <cmath>
#include <complex>

namespace rv::dsp
{
// The optional "commercial capture" low-band stage (docs/model_s_results.md §27): CLA's and the Altiverb's IRs share
// a minimum-phase transfer below ~400 Hz (+3 dB at 40-60 Hz, -0.7 dB at 160 Hz, ripple to 400 Hz) that is in neither
// Tomkin's library nor the stock unit; possibly the manual's transformer option (EMT246_unit_facts.md §2b). Part
// of the TRANSFORMER stage (Transformer.h), on by default.
//
// The 10 numbers are the "CLA" entry of analysis/unit_transfer/bq.json, fitted at the device rate: a first-order
// high-pass, an RBJ low shelf and two RBJ peaking sections, the formulas of tools/unit_transfer.py (hp1_sos,
// lowshelf_sos, peak_sos). The sections are designed at the host rate from the same parameters; that matches
// bq_resp within 0.001 dB at 50 Hz-2 kHz from 22.05 to 192 kHz (tests/TransformerTests.cpp). Below ~50 Hz the fit
// is not defined by the data. Double precision: the corners sit at 40-280 Hz, far below the host rate.
class CaptureEq
{
  public:
    // bq.json "CLA": hp fc, shelf (f0, gain dB, Q), peak (f0, gain dB, Q) x 2.
    static constexpr std::array<double, 10> cla{39.999999999999986, 52.81516142596459,  9.30246847774469,
                                                0.601569759607263,  154.6646171955935,  -0.8993904376426464,
                                                1.2483159003121012, 276.0498874880384,  0.6201888324477692,
                                                3.5043409976570743};

    void prepare(double sampleRate)
    {
        const auto& x = cla;
        sections[0] = highPass1(x[0], sampleRate);
        sections[1] = lowShelf(x[1], x[2], x[3], sampleRate);
        sections[2] = peak(x[4], x[5], x[6], sampleRate);
        sections[3] = peak(x[7], x[8], x[9], sampleRate);
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

    std::array<Coeffs, 4> sections;
    std::array<std::array<State, 4>, 2> state;
};
} // namespace rv::dsp
