#pragma once

#include <juce_core/juce_core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

namespace rv::dsp
{
// One analogue filter of the unit's chain, discretised exactly for band-limited signals (DeviceBridge's building block).
//
// The filter is given as partial fractions h(t) = sum_j w_j Re{R_j e^(p_j t)} (t >= 0; one entry per real pole, w = 1,
// or conjugate pair, w = 2). It runs over samples x[m] at spacing `unit` and is read at arbitrary instants. The samples
// stand for a band-limited signal through phi, a Kaiser-windowed sinc at the HOST Nyquist frequency (half-width W, flat to
// about 1e-6 up to 12.8 kHz); phi is the host's continuous signal on the input side and the capture's anti-alias low-pass on
// the output side. The response to x at time t is then sum_m x[m] k(t - t_m), k = unit * (h * phi), and since phi has
// finite support, k(t) is a sum of pure exponentials for t >= W:
//  - head: the last few samples (|t - t_m| < W plus up to one sample) through a tabulated k, linear between phases;
//  - tail: every older sample through one recursion per mode, q_j = e^(p_j unit) q_j + x, read with e^(p_j u unit),
//    u the fractional position, from a table and a short Taylor series.
// Nothing is truncated: the result is the analogue filter's exact response to phi's signal, to the table's precision
// (about -115 dB). State and arithmetic are double: the two 12.1 kHz pole pairs of the input filter almost coincide and
// their modes cancel by about 50 dB.
class BridgeFilter
{
  public:
    struct Mode
    {
        std::complex<double> pole, residue; // rad/s, 1/s
        int weight = 1;                     // 2 for a conjugate pair
    };
    static constexpr int maxModes = 12;

    struct State
    {
        std::array<double, maxModes> re{}, im{};
        void clear()
        {
            re.fill(0.0);
            im.fill(0.0);
        }
    };

    // phiHalf: phi's half-width (s); hostT: the host sample period (phi's sinc spacing); phases: table rows per unit;
    // gain: a scale for the whole filter. Not real-time.
    void build(const std::vector<Mode>& modes, double unit, double hostT, double phiHalf, double beta, int phases, double gain)
    {
        jassert(modes.size() <= (size_t)maxModes);
        using C = std::complex<double>;
        numModes = (int)modes.size();
        numPhases = std::max(8, phases);
        half = std::max(1, (int)std::ceil(phiHalf / unit - 1e-9));
        const auto taps = 2 * half;
        const auto pi = juce::MathConstants<double>::pi;

        // phi (1/s) on a grid of step unit / (2 P) over [-half, half] units; Simpson's rule on pairs of steps.
        const auto dt = unit / (2.0 * numPhases);
        const auto numGrid = 4 * half * numPhases;
        const auto i0Beta = besselI0(beta);
        std::vector<double> phi((size_t)numGrid + 1);
        for (int g = 0; g <= numGrid; ++g)
        {
            const auto sigma = -half * unit + g * dt;
            const auto r = sigma / phiHalf;
            if (std::abs(r) >= 1.0)
                continue;
            const auto x = sigma / hostT;
            const auto sinc = x == 0.0 ? 1.0 : std::sin(pi * x) / (pi * x);
            phi[(size_t)g] = sinc * besselI0(beta * std::sqrt(1.0 - r * r)) / i0Beta / hostT;
        }

        // k(t) = sum_j w_j Re{R_j e^(p_j t) F_j(t)}, F_j(t) = int_{-W}^{t} e^(-p_j s) phi(s) ds, on the even grid points.
        const auto numK = 2 * half * numPhases;
        std::vector<double> kappa((size_t)numK + 1, 0.0);
        std::array<C, maxModes> phiHat{};
        for (int j = 0; j < numModes; ++j)
        {
            const auto& m = modes[(size_t)j];
            auto f = [&](int g) { return std::exp(-m.pole * (-half * unit + g * dt)) * phi[(size_t)g]; };
            C F(0.0, 0.0);
            auto fPrev = f(0);
            for (int i = 0; i <= numK; ++i)
            {
                if (i > 0)
                {
                    const auto fMid = f(2 * i - 1), fNext = f(2 * i);
                    F += dt / 3.0 * (fPrev + 4.0 * fMid + fNext);
                    fPrev = fNext;
                }
                const auto t = -half * unit + i * (unit / numPhases);
                kappa[(size_t)i] += m.weight * (m.residue * std::exp(m.pole * t) * F).real();
            }
            phiHat[(size_t)j] = F; // phi's transform at -p_j: the tail's weight
        }

        // Head table, rows ph = 0..P (u = ph / P), taps oldest first: k((u + half - 1 - tap) * unit).
        table.assign((size_t)(numPhases + 1) * (size_t)taps, 0.0);
        for (int ph = 0; ph <= numPhases; ++ph)
            for (int tap = 0; tap < taps; ++tap)
                table[(size_t)(ph * taps + tap)] = gain * unit * kappa[(size_t)(ph + (taps - 1 - tap) * numPhases)];

        // Tail: the sample `half` units before the read-out instant and older, through the modes.
        expTable.assign((size_t)numModes * (expSteps + 1) * 2, 0.0);
        for (int j = 0; j < numModes; ++j)
        {
            const auto& m = modes[(size_t)j];
            const auto c = gain * unit * m.weight * m.residue * phiHat[(size_t)j] * std::exp(m.pole * (half * unit));
            for (int k = 0; k <= expSteps; ++k)
            {
                const auto e = c * std::exp(m.pole * (unit * k / expSteps));
                expTable[(size_t)((j * (expSteps + 1) + k) * 2)] = e.real();
                expTable[(size_t)((j * (expSteps + 1) + k) * 2 + 1)] = e.imag();
            }
            const auto a = std::exp(m.pole * unit);
            aRe[(size_t)j] = a.real();
            aIm[(size_t)j] = a.imag();
            const auto z = m.pole * (unit / expSteps);
            zRe[(size_t)j] = z.real();
            zIm[(size_t)j] = z.imag();
        }
    }

    // Head taps either side of the read-out instant: the head reads 2 * halfTaps() samples.
    int halfTaps() const { return half; }

    // The head: x holds the 2 * halfTaps() samples around the instant, oldest first; the instant is u (0 <= u < 1) units
    // after the (halfTaps() - 1)-th of them.
    double head(const double* x, double u) const
    {
        const auto taps = 2 * half;
        const auto pos = u * numPhases;
        const auto ph = std::min((int)pos, numPhases - 1);
        const auto alpha = pos - ph;
        const auto* r0 = table.data() + (size_t)(ph * taps);
        const auto* r1 = r0 + taps;
        double s0 = 0.0, s1 = 0.0;
        for (int k = 0; k < taps; ++k)
        {
            s0 += x[k] * r0[k];
            s1 += x[k] * r1[k];
        }
        return s0 + alpha * (s1 - s0);
    }

    // Feeds the sample `halfTaps()` units older than the head's newest into the tail recursions.
    void push(State& s, double x) const
    {
        for (int j = 0; j < numModes; ++j)
        {
            const auto re = aRe[(size_t)j] * s.re[(size_t)j] - aIm[(size_t)j] * s.im[(size_t)j] + x;
            s.im[(size_t)j] = aRe[(size_t)j] * s.im[(size_t)j] + aIm[(size_t)j] * s.re[(size_t)j];
            s.re[(size_t)j] = re;
        }
    }

    // The tail's contribution at fractional position u (as for head()).
    double tail(const State& s, double u) const
    {
        const auto pos = u * expSteps;
        const auto k = std::min((int)pos, expSteps - 1);
        const auto d = pos - k;
        double acc = 0.0;
        for (int j = 0; j < numModes; ++j)
        {
            // e^(z d), z = p unit / steps, |z d| < 0.005: Taylor to the 3rd order (error < 3e-11 relative).
            const auto zr = zRe[(size_t)j] * d, zi = zIm[(size_t)j] * d;
            auto er = 1.0 + zr / 3.0, ei = zi / 3.0;  // 1 + z/3
            auto tr = 1.0 + 0.5 * (zr * er - zi * ei); // 1 + z/2 (1 + z/3)
            auto ti = 0.5 * (zr * ei + zi * er);
            er = 1.0 + zr * tr - zi * ti; // 1 + z (...)
            ei = zr * ti + zi * tr;
            const auto* c = expTable.data() + (size_t)((j * (expSteps + 1) + k) * 2);
            const auto cr = c[0] * er - c[1] * ei, ci = c[0] * ei + c[1] * er;
            acc += cr * s.re[(size_t)j] - ci * s.im[(size_t)j];
        }
        return acc;
    }

  private:
    static constexpr int expSteps = 1024;

    static double besselI0(double x)
    {
        double sum = 1.0, term = 1.0;
        for (int k = 1; k < 60; ++k)
        {
            term *= (x / (2.0 * k)) * (x / (2.0 * k));
            sum += term;
            if (term < 1e-18 * sum)
                break;
        }
        return sum;
    }

    int numModes = 0, numPhases = 8, half = 1;
    std::vector<double> table, expTable;
    std::array<double, maxModes> aRe{}, aIm{}, zRe{}, zIm{};
};
} // namespace rv::dsp
