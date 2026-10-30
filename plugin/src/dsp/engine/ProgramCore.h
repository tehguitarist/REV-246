#pragma once

#include "params/Parameters.h"

#include <array>
#include <cmath>

namespace rv::dsp
{
// The unit's sample rate: 19.2 MHz / 728. Every program core runs here; all its delays are integer device samples.
inline constexpr double deviceRate = 19.2e6 / 728.0;

// sum a[k] b[k], k < n, in four partial sums: one running sum is a chain of dependent adds the compiler may not
// reorder (without fast-math), about 4x slower for the FIR kernels this is used for.
inline double dotProduct(const double* a, const double* b, int n)
{
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    int k = 0;
    for (; k + 4 <= n; k += 4)
    {
        s0 += a[k] * b[k];
        s1 += a[k + 1] * b[k + 1];
        s2 += a[k + 2] * b[k + 2];
        s3 += a[k + 3] * b[k + 3];
    }
    for (; k < n; ++k)
        s0 += a[k] * b[k];
    return (s0 + s1) + (s2 + s3);
}

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
    // transcendental calls are fine. Takes effect from the next sample, with no smoothing (as the unit's coefficient
    // loads; docs/ENGINE.md, The engine).
    virtual void setSettings(const Settings& s) = 0;
    // Processes `n` device-rate samples. Real-time safe.
    virtual void process(const float* in, float* outL, float* outR, int n) = 0;

    // Samples by which the core's output is late against the reference IRs (N's branch read kernel needs a
    // centred look-ahead; every other core returns 0). The engine adds this to its latency.
    virtual int latency() const { return 0; }
    // Reads every input-side path `samples` device samples sooner: the whole output moves that much earlier, exactly,
    // because every path through the core is at least maxAdvance() samples deep (E: its front delay; N: its branch
    // and front taps, which absorbs the branch kernel's look-ahead). The engine uses it to put every program on one
    // grid without delaying the others (docs/ENGINE.md, One time grid); the reference-IR tests leave it at 0. Call
    // before reset().
    virtual int maxAdvance() const { return 0; }
    virtual void setAdvance(int samples) { juce::ignoreUnused(samples); }
    // Time variance (M / L's theta, E's moving reads). On by default; the null tests switch it off to compare
    // the static network with the reference IRs. Takes effect at the next reset().
    void setModulation(bool on) { modulationOn = on; }

    // The A/D's full scale (a peak) in the core's input units, from CAL: one A/D LSB is fs / 32768. The 16-bit
    // storage grids scale with it. Real-time safe; default 1.0.
    void setStorageFullScale(double fs)
    {
        fullScale = fs > 0.0 ? fs : 1.0;
        for (size_t i = 0; i < grids.size(); ++i)
            grids[i] = makeGrid(fullScale, gridGain[i]);
    }
    double storageFullScale() const { return fullScale; }

    // How a stored value is cut to 16 bits. `toZero` (magnitude truncation, not in s_loops_fx.c) is the default (user,
    // 30 October): floor and nearest leave frozen tails louder than the unit's noise spec. `floor` is s_loops_fx.c's
    // truncation (mode 1, the TMS32010's SACH); `nearest` its rounding (mode 2, + half an LSB first). floor and nearest are bit-exact against s_loops_fx.c (tests/dsp/SixteenBitTests.cpp). Measured
    // there: with floor or nearest a loop can end in a dead band that holds a frozen, loop-periodic pattern for ever
    // after the input stops (S at HF 1.0, M / L / E at long RT); toZero always decays to 0. Which one the unit used is
    // H. Real-time safe.
    enum class StorageCut
    {
        floor,
        toZero,
        nearest
    };
    void setStorageCut(StorageCut c) { cutRule = c; }
    StorageCut storageCut() const { return cutRule; }

  protected:
    // ---- 16-bit storage (docs/ENGINE.md, 16-bit; tools/loop_read/c/s_loops_fx.c; model_s_results §36.5) ----
    // The unit stores every intermediate result as a 16-bit word (unit facts §2, C: 32-bit accumulate, 16-bit
    // storage). As s_loops_fx.c: products are summed in the (wide) accumulator, and each stored value is the
    // accumulator cut to 16 bits (setStorageCut: toZero by default; floor is the TMS32010's SACH, and S, §36.5: only
    // truncation reproduces the spec's 15 dB REVERB noise penalty at a no-overflow scaling) and saturated at -32768...32767 (s_loops_fx.c's
    // store(); H: that the firmware runs the TMS32010 in overflow (saturation) mode rather than letting words wrap).
    //
    // A storage group's internal gain k (the stored word = k x the A/D word for the same signal) is the core's
    // scaling, set per group with setGridGain(); its grid step in core units is fs / (32768 k), so the core's float
    // arithmetic is unchanged and only the stored values are cut. The step is rounded to 8 significant bits (at most
    // 0.2 %, 0.017 dB): then every stored value n x step (|n| <= 2^15) is exact in a float delay line, a stored value
    // re-stored is unchanged, and with Q15 coefficients the floor is the integer floor of s_loops_fx.c bit for bit.
    struct Grid16
    {
        double step = 1.0 / 32768.0, inv = 32768.0;
    };
    static constexpr int maxGrids = 8;

    // Sets storage group `i`'s internal gain k (> 0). Real-time safe.
    void setGridGain(int i, double k)
    {
        gridGain[(size_t)i] = k > 0.0 ? k : 1.0;
        grids[(size_t)i] = makeGrid(fullScale, gridGain[(size_t)i]);
    }
    double gridStep(int i) const { return grids[(size_t)i].step; }

    // Every write into recirculating storage (a loop's or all-pass's delay line) goes through here. Float mode passes
    // the value through, as float; 16-bit mode truncates it to storage group `grid`'s 16-bit word.
    float store(float x) const { return quantise ? (float)cut(x, grids[0]) : x; }
    float store(double x, int grid) const { return quantise ? (float)cut(x, grids[(size_t)grid]) : (float)x; }
    // The same for a recirculating filter state held in double (a loop filter's stored output): the identity in
    // float mode, so float renders stay bit-identical.
    double storeAcc(double acc, int grid) const { return quantise ? cut(acc, grids[(size_t)grid]) : acc; }
    // Unconditional cut, for code paths that only run in 16-bit mode.
    double cut16(double acc, int grid) const { return cut(acc, grids[(size_t)grid]); }

    // A coefficient as a 16-bit Q15 word, rounded to nearest (s_loops_fx.c's q15(): llround(c 2^15)). Values past
    // +-1 keep their 2^-15 step (H: the unit's coefficient format is not known, model_s_results §36.1-36.4).
    static double q15(double c) { return (double)std::llround(c * 32768.0) / 32768.0; }

    // |H(e^jw)|^2 of the first-order section (b0 + b1 z^-1) / (1 + a1 z^-1) at cos w = c, for the cores' scalings.
    static double sectionPower(double b0, double b1, double a1, double c)
    {
        return (b0 * b0 + b1 * b1 + 2.0 * b0 * b1 * c) / (1.0 + a1 * a1 + 2.0 * a1 * c);
    }
    // cos(pi j / (gridPoints - 1)), j = 0..gridPoints-1: the frequencies a core's scaling checks for a loop filter's
    // peak gain (DC and Nyquist included; the filters are first-order cascades, smooth on this grid).
    static constexpr int gridPoints = 33;
    static double gridCos(int j) { return std::cos(std::acos(-1.0) * j / (gridPoints - 1)); }

    // The 16-bit loop filter of s_loops_fx.c for M / L / E (S runs its own sections through storeAcc): the shelf with
    // the loop gain folded in, s = cut(c0 x + c1 x[-1] + a s[-1]), then the HF one-pole y = cut((1 - p) s + p y[-1]);
    // both outputs are stored words (three requantisations per loop with the delay write, model_s_results §36.5).
    struct ShelfPole16
    {
        double c0 = 0.0, c1 = 0.0, a = 0.0, p = 0.0; // coefficients
        double s1 = 0.0;                              // the shelf's stored output, one sample back
    };
    // One sample: x the line read, x1 the previous read, y1 the filter's previous output. 16-bit mode only.
    double runShelfPole16(ShelfPole16& f, double x, double x1, double y1, int grid) const
    {
        const auto& g = grids[(size_t)grid];
        const double s = cut(f.c0 * x + f.c1 * x1 + f.a * f.s1, g);
        f.s1 = s;
        return cut((1.0 - f.p) * s + f.p * y1, g);
    }

    bool quantise = false; // set from Settings::sixteenBit by each core in setSettings()
    bool modulationOn = true;

  private:
    static Grid16 makeGrid(double fs, double k)
    {
        int e = 0;
        const double m = std::frexp(fs / (32768.0 * k), &e); // [0.5, 1)
        Grid16 g;
        g.step = std::ldexp(std::round(m * 256.0) / 256.0, e);
        g.inv = 1.0 / g.step;
        return g;
    }

    // x / step cut to an integer by the StorageCut rule, saturated to a 16-bit word, times step. The 1e-9 LSB guard
    // absorbs the rounding of x * inv for a value already on the grid (exact sums of Q15 x grid products are 2^-15 LSB
    // apart, far coarser). NaN -> -32768.
    double cut(double x, const Grid16& g) const
    {
        const double u = x * g.inv;
        double n;
        if (cutRule == StorageCut::floor)
            n = std::floor(u + 1e-9);
        else if (cutRule == StorageCut::nearest)
            n = std::floor(u + 0.5);
        else
            n = std::trunc(u + (u < 0.0 ? -1e-9 : 1e-9));
        if (!(n >= -32768.0))
            n = -32768.0;
        else if (n > 32767.0)
            n = 32767.0;
        return n * g.step;
    }

    double fullScale = 1.0;
    StorageCut cutRule = StorageCut::toZero; // the plugin's default (user, 30 October): tails decay to silence
    std::array<double, maxGrids> gridGain{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    std::array<Grid16, maxGrids> grids{};
};
} // namespace rv::dsp
