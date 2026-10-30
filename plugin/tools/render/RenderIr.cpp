// REV246Render: offline renders of the plugin's reverb engine (src/dsp/ReverbEngine.h) for the phase 5 check
// against the hardware IRs (tools/plugin_check/, analysis/plugin_check/). No editor, no processor, no parameter tables:
// the Settings are given as values.
//
//   REV246Render <list.txt>        one render per non-empty line that does not start with '#'
//   REV246Render key=value ...     one render
//
// Keys (defaults in brackets):
//   out=<path.wav>       output, float32 stereo (L, R) WAV (required)
//   rate=<Hz>            host rate [96000]
//   program=S|M|L|E|N    [S]
//   lf= rt= hf= ampl=    panel values (N: lf = EARLY RT s, hf = LATE RT s, ampl = SPAN s) [1.0 1.0 0.9 80]
//   delay=<ms> [0]  lp=<kHz> [12, the exact bypass]  sixteen=0|1 [0]
//   mod=0|1              time variance (M / L theta, E's moving reads) [1]; 0 = the cores' null-test switch
//   len=<s>              output length after the input starts [3]
//   pre=<s>              silence before the input (moves the time variance's phase at the input) [0]
//   input=impulse|<path.wav>   a unit impulse at host sample round(pre x rate), or a mono WAV (first channel) [impulse]
//   gain=<x>             input gain [1]
//   domain=host|device   host: the whole engine at the host rate (bridge, analogue chain, latency left in the file);
//                        device: the device-rate core output (alignment, DELAY, LOW PASS; no bridge, no chain), input
//                        impulse at device sample round(pre x deviceRate), file rate the integer nearest the device rate
//                        (the true rate is 19.2 MHz / 728) [host]
//   block=<n>            host block size [512]
//
// Each render prints one line: out, the engine's latency (host samples, R) and alignment (device samples).
// The output file starts at the first input sample (the latency is in the file, not removed).

#include "dsp/ReverbEngine.h"

#include <juce_audio_basics/juce_audio_basics.h>
#include <juce_core/juce_core.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
using rv::dsp::ReverbEngine;
using rv::dsp::Settings;
using rv::params::Program;

struct Wav
{
    double rate = 0.0;
    int channels = 0;
    std::vector<float> data; // interleaved
};

bool readWav(const std::string& path, Wav& w)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    std::vector<char> b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    auto u16 = [&](size_t i) { return (uint32_t)(uint8_t)b[i] | ((uint32_t)(uint8_t)b[i + 1] << 8); };
    auto u32 = [&](size_t i) { return u16(i) | (u16(i + 2) << 16); };
    if (b.size() < 12 || std::memcmp(b.data(), "RIFF", 4) != 0 || std::memcmp(b.data() + 8, "WAVE", 4) != 0)
        return false;
    int fmt = 0, bits = 0;
    size_t p = 12;
    while (p + 8 <= b.size())
    {
        const auto sz = (size_t)u32(p + 4);
        if (std::memcmp(b.data() + p, "fmt ", 4) == 0)
        {
            fmt = (int)u16(p + 8);
            w.channels = (int)u16(p + 10);
            w.rate = (double)u32(p + 12);
            bits = (int)u16(p + 22);
            if (fmt == 0xFFFE && sz >= 40)
                fmt = (int)u16(p + 32); // WAVE_FORMAT_EXTENSIBLE: the sub-format's first two bytes
        }
        else if (std::memcmp(b.data() + p, "data", 4) == 0)
        {
            const auto bytes = std::min(sz, b.size() - p - 8);
            const auto* d = b.data() + p + 8;
            const int bps = bits / 8;
            const auto n = bytes / (size_t)bps;
            w.data.resize(n);
            for (size_t i = 0; i < n; ++i)
            {
                const auto* s = d + i * (size_t)bps;
                if (fmt == 3 && bits == 32)
                {
                    float v;
                    std::memcpy(&v, s, 4);
                    w.data[i] = v;
                }
                else if (fmt == 3 && bits == 64)
                {
                    double v;
                    std::memcpy(&v, s, 8);
                    w.data[i] = (float)v;
                }
                else if (fmt == 1 && bits == 16)
                    w.data[i] = (float)(int16_t)u16((size_t)(s - b.data())) / 32768.0f;
                else if (fmt == 1 && bits == 24)
                {
                    int32_t v = (int32_t)((uint32_t)(uint8_t)s[0] << 8 | (uint32_t)(uint8_t)s[1] << 16 | (uint32_t)(uint8_t)s[2] << 24);
                    w.data[i] = (float)((double)(v >> 8) / 8388608.0);
                }
                else if (fmt == 1 && bits == 32)
                {
                    int32_t v;
                    std::memcpy(&v, s, 4);
                    w.data[i] = (float)((double)v / 2147483648.0);
                }
                else
                    return false;
            }
            return w.channels > 0;
        }
        p += 8 + sz + (sz & 1);
    }
    return false;
}

bool writeWav(const std::string& path, int rate, const std::vector<float>& l, const std::vector<float>& r)
{
    std::ofstream f(path, std::ios::binary);
    if (!f)
        return false;
    const auto n = (uint32_t)l.size();
    const uint32_t dataBytes = n * 2 * 4;
    auto put32 = [&](uint32_t v) { f.write(reinterpret_cast<const char*>(&v), 4); };
    auto put16 = [&](uint16_t v) { f.write(reinterpret_cast<const char*>(&v), 2); };
    f.write("RIFF", 4);
    put32(36 + dataBytes);
    f.write("WAVEfmt ", 8);
    put32(16);
    put16(3); // IEEE float
    put16(2);
    put32((uint32_t)rate);
    put32((uint32_t)rate * 8);
    put16(8);
    put16(32);
    f.write("data", 4);
    put32(dataBytes);
    for (uint32_t i = 0; i < n; ++i)
    {
        f.write(reinterpret_cast<const char*>(&l[i]), 4);
        f.write(reinterpret_cast<const char*>(&r[i]), 4);
    }
    return (bool)f;
}

using Args = std::map<std::string, std::string>;

std::string get(const Args& a, const std::string& k, const std::string& def)
{
    const auto it = a.find(k);
    return it == a.end() ? def : it->second;
}
double num(const Args& a, const std::string& k, double def)
{
    const auto it = a.find(k);
    return it == a.end() ? def : std::stod(it->second);
}

Program programOf(const std::string& s)
{
    const auto c = s.empty() ? 'S' : (char)std::toupper((unsigned char)s[0]);
    switch (c)
    {
    case 'M': return Program::m;
    case 'L': return Program::l;
    case 'E': return Program::e;
    case 'N': return Program::n;
    default: return Program::s;
    }
}

int render(const Args& a)
{
    const auto out = get(a, "out", "");
    if (out.empty())
    {
        std::fprintf(stderr, "missing out=\n");
        return 1;
    }
    const double rate = num(a, "rate", 96000.0);
    const int block = (int)num(a, "block", 512);
    const bool device = get(a, "domain", "host") == "device";
    const double fileRate = device ? rv::dsp::deviceRate : rate;

    Settings s;
    s.program = programOf(get(a, "program", "S"));
    s.lfX = num(a, "lf", 1.0);
    s.rt = num(a, "rt", 1.0);
    s.hfX = num(a, "hf", 0.9);
    s.amplPercent = num(a, "ampl", 80.0);
    s.delayMs = num(a, "delay", 0.0);
    s.lowpassKHz = num(a, "lp", 12.0);
    s.sixteenBit = num(a, "sixteen", 0.0) > 0.5;
    const bool mod = num(a, "mod", 1.0) > 0.5;

    // Input at the file's rate.
    const auto pre = (size_t)std::llround(num(a, "pre", 0.0) * fileRate);
    const auto len = (size_t)std::llround(num(a, "len", 3.0) * fileRate);
    const double gain = num(a, "gain", 1.0);
    std::vector<float> in(pre + len, 0.0f);
    const auto inputSpec = get(a, "input", "impulse");
    if (inputSpec == "impulse")
        in[pre] = (float)gain;
    else
    {
        Wav w;
        if (!readWav(inputSpec, w))
        {
            std::fprintf(stderr, "cannot read %s\n", inputSpec.c_str());
            return 1;
        }
        if (std::abs(w.rate - fileRate) > 0.5 && !device)
        {
            std::fprintf(stderr, "%s: rate %.0f, render rate %.0f\n", inputSpec.c_str(), w.rate, fileRate);
            return 1;
        }
        const auto n = w.data.size() / (size_t)w.channels;
        for (size_t i = 0; i < n && pre + i < in.size(); ++i)
            in[pre + i] = (float)(gain * w.data[i * (size_t)w.channels]);
    }

    auto engine = std::make_unique<ReverbEngine>();
    engine->prepare(rate, block);
    engine->setSettings(s);
    engine->setModulation(mod);
    engine->reset(); // start in s.program (no crossfade) with the modulation switch applied

    std::vector<float> l(in.size()), r(in.size());
    if (device)
    {
        // The device-rate path alone: alignment ring, core, DELAY, LOW PASS.
        const int chunk = engine->maxDeviceBlock();
        for (size_t i = 0; i < in.size(); i += (size_t)chunk)
        {
            const int n = (int)std::min((size_t)chunk, in.size() - i);
            engine->processDeviceRate(in.data() + i, l.data() + i, r.data() + i, n);
        }
    }
    else
        for (size_t i = 0; i < in.size(); i += (size_t)block)
        {
            const int n = (int)std::min((size_t)block, in.size() - i);
            engine->process(in.data() + i, l.data() + i, r.data() + i, n);
        }

    // The file starts at the input's first sample (pre-roll dropped).
    std::vector<float> ol(l.begin() + (long)pre, l.end()), orr(r.begin() + (long)pre, r.end());
    if (!writeWav(out, (int)std::lround(fileRate), ol, orr))
    {
        std::fprintf(stderr, "cannot write %s\n", out.c_str());
        return 1;
    }
    std::printf("%s latency=%d align=%d rate=%.6f\n", out.c_str(), device ? 0 : engine->latency(), engine->alignmentDelay(),
                fileRate);
    return 0;
}

Args parse(const std::vector<std::string>& toks)
{
    Args a;
    for (const auto& t : toks)
    {
        const auto eq = t.find('=');
        if (eq != std::string::npos)
            a[t.substr(0, eq)] = t.substr(eq + 1);
    }
    return a;
}
} // namespace

int main(int argc, char** argv)
{
    std::vector<std::string> toks(argv + 1, argv + argc);
    if (toks.empty())
    {
        std::fprintf(stderr, "usage: REV246Render <list.txt> | key=value ...\n");
        return 2;
    }
    if (toks.size() == 1 && toks[0].find('=') == std::string::npos)
    {
        std::ifstream f(toks[0]);
        if (!f)
        {
            std::fprintf(stderr, "cannot open %s\n", toks[0].c_str());
            return 1;
        }
        std::string line;
        int bad = 0;
        while (std::getline(f, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream ss(line);
            std::vector<std::string> t;
            for (std::string w; ss >> w;)
                t.push_back(w);
            bad += render(parse(t)) != 0;
            std::fflush(stdout);
        }
        return bad == 0 ? 0 : 1;
    }
    return render(parse(toks));
}
