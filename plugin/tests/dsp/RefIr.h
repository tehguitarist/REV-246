#pragma once

// Shared helpers for the DSP tests: load a reference IR, render a core's impulse response, and null the two in the
// spec tools' windows (tools/s_spec.py, tools/ml_spec.py: 10 log10 of the error energy over the reference energy,
// per window and channel).

#include "dsp/engine/ProgramCore.h"

#include <juce_audio_formats/juce_audio_formats.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace rv::test
{
struct StereoIr
{
    std::vector<double> l, r;
    size_t size() const { return l.size(); }
};

// A file under the repository root (RV_REPO_ROOT is set by CMake).
inline juce::File repoFile(const std::string& relative)
{
    return juce::File(RV_REPO_ROOT).getChildFile(relative);
}

// Loads a stereo WAV (the reference IRs are float32 at the device rate, (L, R), impulse at sample 0).
inline StereoIr loadWav(const std::string& relative)
{
    StereoIr ir;
    juce::WavAudioFormat wav;
    auto stream = repoFile(relative).createInputStream(); // null for a missing file: return an empty IR, not a crash
    if (stream == nullptr)
        return ir;
    std::unique_ptr<juce::AudioFormatReader> reader(wav.createReaderFor(stream.release(), true));
    if (reader == nullptr)
        return ir;
    const auto n = (int)reader->lengthInSamples;
    juce::AudioBuffer<float> b((int)reader->numChannels, n);
    reader->read(&b, 0, n, 0, true, true);
    ir.l.resize((size_t)n);
    ir.r.resize((size_t)n);
    for (int i = 0; i < n; ++i)
    {
        ir.l[(size_t)i] = b.getSample(0, i);
        ir.r[(size_t)i] = b.getSample(b.getNumChannels() > 1 ? 1 : 0, i);
    }
    return ir;
}

// The core's response to a unit impulse at sample 0, `n` samples, with its latency() removed.
inline StereoIr renderImpulse(dsp::ProgramCore& core, size_t n, int block = 256)
{
    const auto lat = (size_t)core.latency();
    const auto total = n + lat;
    std::vector<float> in(total, 0.0f), l(total), r(total);
    in[0] = 1.0f;
    for (size_t i = 0; i < total; i += (size_t)block)
    {
        const auto m = (int)std::min((size_t)block, total - i);
        core.process(in.data() + i, l.data() + i, r.data() + i, m);
    }
    StereoIr ir;
    ir.l.assign(l.begin() + (long)lat, l.end());
    ir.r.assign(r.begin() + (long)lat, r.end());
    return ir;
}

// Null in dB of `got` against `ref` over [t0, t1) ms; window energy of the reference as the denominator.
inline double nullDb(const std::vector<double>& ref, const std::vector<double>& got, double t0Ms, double t1Ms)
{
    const auto a = (size_t)(t0Ms / 1000.0 * dsp::deviceRate);
    const auto b = std::min({(size_t)(t1Ms / 1000.0 * dsp::deviceRate), ref.size(), got.size()});
    double e = 0.0, s = 0.0;
    for (size_t i = a; i < b; ++i)
    {
        const auto d = got[i] - ref[i];
        e += d * d;
        s += ref[i] * ref[i];
    }
    if (s <= 0.0)
        return e <= 0.0 ? -400.0 : 400.0;
    return 10.0 * std::log10(std::max(e, 1e-300) / s);
}

inline constexpr std::array<std::array<double, 2>, 5> specWindows{{{0, 50}, {50, 150}, {150, 300}, {300, 600}, {600, 1000}}};

// The worst (largest) null over the spec windows and both channels, printing one line per file.
inline double worstNull(const std::string& name, const StereoIr& ref, const StereoIr& got)
{
    double worst = -400.0;
    std::string line = name + ":";
    for (const auto& w : specWindows)
    {
        if (w[0] / 1000.0 * dsp::deviceRate >= (double)ref.size())
            break;
        const auto nl = nullDb(ref.l, got.l, w[0], w[1]), nr = nullDb(ref.r, got.r, w[0], w[1]);
        worst = std::max({worst, nl, nr});
        char buf[96];
        std::snprintf(buf, sizeof buf, "  %g-%g ms %.1f/%.1f", w[0], w[1], nl, nr);
        line += buf;
    }
    std::printf("%s\n", line.c_str());
    return worst;
}
} // namespace rv::test
