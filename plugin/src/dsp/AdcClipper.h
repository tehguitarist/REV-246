#pragma once

#include <juce_core/juce_core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

namespace rv::dsp
{
// The unit's A/D converter running out of range (16-bit mode). Design and checks: tools/analog/adc_clip.py,
// analysis/analog/adc_clip.json; the facts: docs/EMT246_unit_facts.md section 2.
//
// The A/D is a discrete 16-bit successive-approximation converter (PCM-type DAC, comparator, two 8-bit registers)
// behind the input filter 7 246 131 and an integrate/hold S&H (C). A SAR cannot wrap: past full scale every one of
// its 16 comparisons gives the same answer, so the code sticks at the end of the range, a hard clip of each held
// sample at the device rate (H: the stages before it are guessed to have more swing than the converter). Two things
// follow, and both are modelled:
//  - the A/D sees the input filter's output, which lifts the HF by up to +5.1 dB (10-11.5 kHz) before an elliptic
//    skirt at 12.2 kHz. So high frequencies clip up to 5 dB earlier than a level reading suggests (the manual gives
//    the 6 dB headroom for f < 2 kHz only), and content above about 13 kHz never reaches it;
//  - the clip happens on samples at 26.4 kHz, so its harmonics fold back around 13.2 kHz into inharmonic tones.
//
// The engine has no input filter yet, so this works beside the signal: it runs the input filter (matched-z, every
// pole and zero at its own frequency, plus one zero tuned for the flattest magnitude), samples its output at the
// device instants with a windowed sinc, clips there, and turns each clip error back into the input domain with
// g0 = IFT{T(w) / H(jw)} (H the input filter, T a taper to 0 by 12.6 kHz that is also the reconstruction low-pass).
// The output is the input delayed by latency() samples plus `amount` times those corrections (0 = the input delayed:
// float mode; it runs anyway, for warm state and the overload flag), so with no clipping it is the input delayed exactly, and passing the output through the input filter gives the clipped A/D samples to the taper's
// floor (-35 dB of the clip error for a sine). The future engine must apply the input filter AFTER this stage.
//
// Full scale is a peak level in the input domain at 1 kHz (where the input filter's gain is taken as 1).
class AdcClipper
{
  public:
    static constexpr double deviceRate = 19.2e6 / 728.0;
    // The unit's nominal level (its 0 dB LED) is 6 dB below the A/D's overload point, for a sine below 2 kHz (C).
    static constexpr double headroomDb = 6.0;

    // The full-scale peak for a meter whose 0 VU reads `zeroVuDbFs` (RMS): a sine 6 dB over 0 VU just fits.
    static float fullScaleForZeroVu(float zeroVuDbFs)
    {
        return juce::MathConstants<float>::sqrt2 * std::pow(10.0f, (zeroVuDbFs + (float)headroomDb) / 20.0f);
    }

    void prepare(double sampleRate)
    {
        using C = std::complex<double>;
        rate = sampleRate;
        const auto T = 1.0 / rate;
        const auto twoPi = juce::MathConstants<double>::twoPi;

        // Matched-z sections: (pole pair, jw zero pair) x 3, (pole pair, real zero and -alpha), pole pair, real pole.
        const auto zp = [T](C r) { return std::exp(r * T); };
        auto pairDen = [&](C p, Section& s) {
            const auto z = zp(p);
            s.a1 = -2.0 * z.real();
            s.a2 = std::norm(z);
        };
        constexpr std::array<size_t, 3> withJwZero{3, 4, 2}; // the 12.13 a, 12.13 b and 10.54 kHz pole pairs
        for (int i = 0; i < 3; ++i)
        {
            auto& s = sections[(size_t)i];
            pairDen(polePairs[withJwZero[(size_t)i]], s);
            const auto z = zp(C(0.0, jwZeros[(size_t)i]));
            s.b0 = 1.0;
            s.b1 = -2.0 * z.real();
            s.b2 = std::norm(z);
        }
        pairDen(polePairs[1], sections[3]);
        pairDen(polePairs[0], sections[4]);
        sections[4].b0 = 1.0;
        sections[4].b1 = sections[4].b2 = 0.0;
        sections[5].b0 = 1.0;
        sections[5].b1 = sections[5].b2 = sections[5].a2 = 0.0;
        sections[5].a1 = -zp(C(realPole)).real();

        // The extra zero at -alpha: the smallest peak error against the analogue magnitude, 20 Hz to 12.6 kHz.
        const auto z0 = zp(C(realZero)).real();
        auto fixedPart = [&](double hz) {
            const auto zi = std::polar(1.0, -twoPi * hz * T); // z^-1
            C h(1.0, 0.0);
            for (size_t i = 0; i < sections.size(); ++i)
            {
                const auto& s = sections[i];
                const auto num = i == 3 ? 1.0 - z0 * zi : s.b0 + zi * (s.b1 + zi * s.b2);
                h *= num / (1.0 + zi * (s.a1 + zi * s.a2));
            }
            return h;
        };
        constexpr int numFit = 64;
        std::array<C, numFit> fixed{};
        std::array<double, numFit> target{};
        const auto a1k = analog(1000.0);
        for (int k = 0; k < numFit; ++k)
        {
            const auto hz = 20.0 + (12.6e3 - 20.0) * k / (numFit - 1);
            fixed[(size_t)k] = fixedPart(hz);
            target[(size_t)k] = std::abs(analog(hz) / a1k);
        }
        const auto f1k = fixedPart(1000.0);
        const auto zi1k = std::polar(1.0, -twoPi * 1000.0 * T);
        double bestErr = 1e9;
        for (int a = 0; a <= 200; ++a)
        {
            const auto alpha = a / 200.0;
            const auto norm = std::abs(f1k * (1.0 + alpha * zi1k));
            double err = 0.0;
            for (int k = 0; k < numFit; ++k)
            {
                const auto zi = std::polar(1.0, -twoPi * (20.0 + (12.6e3 - 20.0) * k / (numFit - 1)) * T);
                const auto mag = std::abs(fixed[(size_t)k] * (1.0 + alpha * zi)) / norm;
                err = std::max(err, std::abs(std::log(mag / target[(size_t)k])));
            }
            if (err < bestErr)
                bestErr = err, extraZero = alpha;
        }
        sections[3].b0 = 1.0;
        sections[3].b1 = extraZero - z0;
        sections[3].b2 = -extraZero * z0;
        const auto h1k = f1k * (1.0 + extraZero * zi1k);
        gain = 1.0 / std::abs(h1k);
        // Matched-z is off from the analogue filter by a near-constant delay (its phase is linear in f): the A/D
        // samples belong to analogue time t + lead.
        leadSamples = std::arg(h1k / a1k) / (twoPi * 1000.0) * rate;

        hpPole = std::exp(-twoPi * couplingHz * T);

        // The device-instant sampler: a Kaiser-windowed sinc, tabulated over |tau| = 0 .. interpHalf, times 1 / rate.
        halfW = interpHalf * rate;
        const auto fc = std::min(interpFc, 0.5 * rate);
        interpTable.resize(tableSize + 2);
        for (int i = 0; i < tableSize + 2; ++i)
        {
            const auto u = std::min(1.0, (double)i / tableSize);
            const auto tau = u * interpHalf;
            const auto x = 2.0 * fc * tau;
            const auto sinc = x == 0.0 ? 1.0 : std::sin(juce::MathConstants<double>::pi * x) / (juce::MathConstants<double>::pi * x);
            interpTable[(size_t)i] = 2.0 * fc * sinc * besselI0(interpBeta * std::sqrt(1.0 - u * u)) / besselI0(interpBeta) / rate;
        }
        step = rate / deviceRate;

        (void)g0Table(); // built once, here rather than at the first clip on the audio thread
        preS = pre * rate;
        postS = post * rate;
        delay = (int)std::ceil(halfW + preS + std::max(0.0, -leadSamples)) + 2;
        vBuf.assign((size_t)juce::nextPowerOfTwo((int)(2.0 * halfW + step) + 8), 0.0);
        dryBuf.assign((size_t)juce::nextPowerOfTwo(delay + 1), 0.0f);
        acc.assign((size_t)juce::nextPowerOfTwo(delay + (int)postS + 8), 0.0);
        reset();
    }

    void reset()
    {
        for (auto& s : sections)
            s.s1 = s.s2 = 0.0;
        hpX = hpY = 0.0;
        std::fill(vBuf.begin(), vBuf.end(), 0.0);
        std::fill(dryBuf.begin(), dryBuf.end(), 0.0f);
        std::fill(acc.begin(), acc.end(), 0.0);
        vPos = dryPos = accPos = 0;
        nextDevice = 0.0;
        overload = false;
    }

    void setFullScale(float peak)
    {
        fsLow = -(double)peak;
        fsHigh = (double)peak * 32767.0 / 32768.0; // two's complement: one code short at the top
    }

    int latency() const { return delay; }

    // The input filter's gain re 1 kHz (the analogue model the A/D side chain follows).
    static double inputFilterGain(double hz) { return std::abs(analog(hz) / analog(1000.0)); }

    // True if the A/D clipped since the last call.
    bool takeOverload()
    {
        const auto o = overload;
        overload = false;
        return o;
    }

    float process(float in, float amount = 1.0f)
    {
        const auto dryMask = dryBuf.size() - 1;
        dryBuf[dryPos] = in;
        const auto delayed = dryBuf[(dryPos - (size_t)delay) & dryMask];
        dryPos = (dryPos + 1) & dryMask;

        // What the A/D sees: the input filter, then the coupling into its summing point (C106 / R117, 3.35 Hz).
        double v = in * gain;
        for (auto& s : sections)
        {
            const auto y = s.b0 * v + s.s1;
            s.s1 = s.b1 * v - s.a1 * y + s.s2;
            s.s2 = s.b2 * v - s.a2 * y;
            v = y;
        }
        const auto hp = 0.5 * (1.0 + hpPole) * (v - hpX) + hpPole * hpY;
        hpX = v;
        hpY = hp;
        const auto vMask = vBuf.size() - 1;
        vPos = (vPos + 1) & vMask;
        vBuf[vPos] = hp;

        // Every device instant whose window is complete: sample, clip, and place the correction.
        nextDevice -= 1.0;
        while (nextDevice <= -halfW)
        {
            const auto p = nextDevice;
            double sample = 0.0;
            for (auto j = (int)std::ceil(p - halfW); j <= (int)std::floor(p + halfW); ++j)
                sample += vBuf[(size_t)((long long)vPos + j) & vMask] * interp(std::abs(p - j));
            if (sample < fsLow || sample > fsHigh)
            {
                overload = true;
                place(std::clamp(sample, fsLow, fsHigh) - sample, p + delay + leadSamples);
            }
            nextDevice += step;
        }

        const auto accMask = acc.size() - 1;
        const auto out = delayed + amount * (float)acc[accPos];
        acc[accPos] = 0.0;
        accPos = (accPos + 1) & accMask;
        return out;
    }

  private:
    struct Section
    {
        double b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0, s1 = 0, s2 = 0;
    };

    // The input filter's poles and zeros (rad/s), from the transcribed netlist by tools/analog/adc_clip.py; its own
    // sub-Hz coupling (a pole at 0.62 Hz, a zero at DC) is left out. Pole pairs by frequency: 5.23, 8.16, 10.54,
    // 12.13, 12.13 kHz; the zeros on the jw axis are at 14.31, 15.76 and 23.21 kHz.
    static constexpr double realPole = -12461.369753754747, realZero = -17006.802721088432;
    static inline const std::array<std::complex<double>, 5> polePairs{{{-19405.14045501846, 26494.395519821675},
                                                                       {-15216.10370999179, 48963.1628948609},
                                                                       {-10017.828922011184, 65476.47699572555},
                                                                       {-4451.567515479799, 76083.9747094202},
                                                                       {-4460.209234086727, 76083.46860203434}}};
    static constexpr std::array<double, 3> jwZeros{89917.94340602914, 99035.025119406, 145852.83765769782};

    static constexpr double taperLo = 11.5e3, taperHi = 12.6e3; // T(w): 1 below, raised cosine, 0 above
    static constexpr double pre = 1.2e-3, post = 0.3e-3;        // g0's support around its centre (mostly anticausal)
    static constexpr double interpHalf = 0.2e-3, interpFc = 22e3, interpBeta = 7.0;
    static constexpr double couplingHz = 3.35;
    static constexpr int tableSize = 4096;
    static constexpr double g0Step = 0.5e-6;

    static std::complex<double> analog(double hz)
    {
        const std::complex<double> s(0.0, juce::MathConstants<double>::twoPi * hz);
        auto h = (s - realZero) / (s - realPole);
        for (const auto w : jwZeros)
            h *= (s * s + w * w);
        for (const auto& p : polePairs)
            h /= (s - p) * (s - std::conj(p));
        return h;
    }

    static double besselI0(double x)
    {
        double sum = 1.0, term = 1.0;
        for (int k = 1; k < 40; ++k)
        {
            term *= (x / (2.0 * k)) * (x / (2.0 * k));
            sum += term;
        }
        return sum;
    }

    // g0(tau) / deviceRate on a 0.5 us grid from -pre to +post: the same for every host rate, so built once.
    static const std::vector<double>& g0Table()
    {
        static const std::vector<double> table = [] {
            const auto twoPi = juce::MathConstants<double>::twoPi;
            constexpr double df = 5.0;
            const auto numF = (int)(taperHi / df) + 1;
            const auto norm = std::abs(analog(1000.0));
            std::vector<std::complex<double>> G((size_t)numF);
            for (int i = 0; i < numF; ++i)
            {
                const auto f = i * df;
                const auto t = f <= taperLo ? 1.0
                                            : 0.5 * (1.0 + std::cos(juce::MathConstants<double>::pi * (f - taperLo) / (taperHi - taperLo)));
                const auto w = (i == 0 || i == numF - 1) ? 0.5 * df : df; // trapezoid
                G[(size_t)i] = t * norm / analog(f) * w;
            }
            const auto n = (int)std::lround((pre + post) / g0Step) + 1;
            std::vector<double> out((size_t)n);
            for (int k = 0; k < n; ++k)
            {
                const auto tau = -pre + k * g0Step;
                const auto rot = std::polar(1.0, twoPi * df * tau);
                std::complex<double> ph(1.0, 0.0);
                double sum = 0.0;
                for (const auto& g : G)
                {
                    sum += (g * ph).real();
                    ph *= rot;
                }
                // A raised-cosine fade over the last 50 us at each end.
                const auto edge = std::clamp(std::min(tau + pre, post - tau) / 50e-6, 0.0, 1.0);
                out[(size_t)k] = 2.0 * sum * 0.5 * (1.0 - std::cos(juce::MathConstants<double>::pi * edge)) / deviceRate;
            }
            return out;
        }();
        return table;
    }

    double interp(double distance) const
    {
        const auto x = distance / halfW * tableSize;
        const auto i = (int)x;
        if (i >= tableSize)
            return 0.0;
        const auto f = x - i;
        return interpTable[(size_t)i] + f * (interpTable[(size_t)i + 1] - interpTable[(size_t)i]);
    }

    // Adds e * g0 centred at `centre` (host samples after the sample being output now) into the accumulator.
    void place(double e, double centre)
    {
        const auto& g0 = g0Table();
        const auto accMask = acc.size() - 1;
        const auto last = (int)std::floor(centre + postS);
        for (auto m = std::max(0, (int)std::ceil(centre - preS)); m <= last; ++m)
        {
            const auto x = ((m - centre) / rate + pre) / g0Step;
            const auto i = (int)x;
            if (i < 0 || i + 1 >= (int)g0.size())
                continue;
            const auto f = x - i;
            acc[(accPos + (size_t)m) & accMask] += e * (g0[(size_t)i] + f * (g0[(size_t)i + 1] - g0[(size_t)i]));
        }
    }

    double rate = 48000.0;
    std::array<Section, 6> sections{};
    double gain = 1.0, extraZero = 0.0, leadSamples = 0.0;
    double hpPole = 0.0, hpX = 0.0, hpY = 0.0;
    double halfW = 1.0, step = 1.0, preS = 0.0, postS = 0.0, nextDevice = 0.0;
    double fsLow = -1.0, fsHigh = 1.0;
    int delay = 0;
    std::vector<double> interpTable, vBuf, acc;
    std::vector<float> dryBuf;
    size_t vPos = 0, dryPos = 0, accPos = 0;
    bool overload = false;
};
} // namespace rv::dsp
