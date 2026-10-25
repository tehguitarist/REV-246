#include "PluginProcessor.h"

#include "PluginEditor.h"

#include <cmath>

namespace
{
constexpr double rampSeconds = 0.02;
constexpr float vuRiseSeconds = 0.3f; // VU: 99% of a step in ~300 ms
// Stereo-to-mono summation (original mode only): 0.5 = -6 dB, so two identical/correlated channels stay at unity
// instead of gaining 6 dB. The equal-power alternative, 1/sqrt(2), would keep uncorrelated stereo at its power but
// add 3 dB to centred material; input headroom matters more here, so use the conservative sum.
constexpr float monoSumGain = 0.5f;
}

Rev246Processor::Rev246Processor()
    : AudioProcessor(BusesProperties()
                         .withInput("Input", juce::AudioChannelSet::stereo(), true)
                         .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      apvts(*this, nullptr, "STATE", rv::params::createLayout())
{
    namespace id = rv::params::id;
    pProgram = apvts.getRawParameterValue(id::program);
    pMix = apvts.getRawParameterValue(id::mix);
    pPan = apvts.getRawParameterValue(id::pan);
    pIn = apvts.getRawParameterValue(id::inputGain);
    pOut = apvts.getRawParameterValue(id::outputGain);
    pPower = apvts.getRawParameterValue(id::power);
    pTransformer = apvts.getRawParameterValue(id::transformer);
    pTrueStereo = apvts.getRawParameterValue(id::trueStereo);
    for (int i = 0; i < rv::params::numJoys; ++i)
        pJoy[(size_t)i] = apvts.getRawParameterValue(rv::params::joyId((rv::params::Joy)i));
}

bool Rev246Processor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    // Mono in/mono out, mono in/stereo out (the unit's own 1-in/2-out path), or stereo in/stereo out (with the
    // original / true-stereo choice). Stereo in/mono out has no sensible meaning for the model.
    const auto& in = layouts.getMainInputChannelSet();
    const auto& out = layouts.getMainOutputChannelSet();
    return (in == juce::AudioChannelSet::mono() && out == juce::AudioChannelSet::mono()) ||
           (in == juce::AudioChannelSet::mono() && out == juce::AudioChannelSet::stereo()) ||
           (in == juce::AudioChannelSet::stereo() && out == juce::AudioChannelSet::stereo());
}

void Rev246Processor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    for (auto& e : engines)
        e.prepare(sampleRate, samplesPerBlock);
    for (auto& t : transformers)
        t.prepare(sampleRate);
    for (auto* s : {&inGain, &outGain, &mixAmount, &panLeft, &panRight, &wetOn, &transformerOn})
        s->reset(sampleRate, rampSeconds);
    inGain.setCurrentAndTargetValue(rv::params::gainLinear(*pIn));
    outGain.setCurrentAndTargetValue(rv::params::gainLinear(*pOut));
    mixAmount.setCurrentAndTargetValue(*pMix);
    wetOn.setCurrentAndTargetValue(*pPower >= 0.5f ? 1.0f : 0.0f);
    transformerOn.setCurrentAndTargetValue(*pTransformer >= 0.5f ? 1.0f : 0.0f);
    scratch.setSize(7, samplesPerBlock, false, false, true);
    // VU meter: a first-order RMS integrator. The classic VU reaches 99% in ~300 ms, so tau = 300 ms / ln 100.
    vuTauSamples = (vuRiseSeconds / std::log(100.0f)) * sampleRate;
    vu = 0.0f;
    inputVu.store(0.0f);
    inputClip.store(false);
    lastTrueStereo = isStereoToStereo() && *pTrueStereo >= 0.5f;
    setLatencySamples(0);
}

void Rev246Processor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    using namespace rv::params;
    const auto n = buffer.getNumSamples();
    const auto numIn = getTotalNumInputChannels();
    const auto numOut = getTotalNumOutputChannels();
    const bool stereoIn = numIn > 1, stereoOut = numOut > 1;
    // The original mode (the default, and the only option outside stereo in / stereo out) sums the input to the unit's mono
    // input. True stereo runs one engine per input channel.
    const bool trueStereo = stereoIn && stereoOut && *pTrueStereo >= 0.5f;
    if (n > scratch.getNumSamples())
        scratch.setSize(7, n, false, false, true); // only if a host exceeds its promised block size

    // A new routing has no usable state in the second engine: start it silent rather than with an old tail.
    if (trueStereo != lastTrueStereo)
    {
        engines[1].reset();
        transformers[1].reset();
        lastTrueStereo = trueStereo;
    }

    // Settings from the panel, converted to values for the engine.
    rv::dsp::Settings s;
    s.program = (Program)juce::jlimit(0, numPrograms - 1, juce::roundToInt(pProgram->load()));
    for (int i = 0; i < numJoys; ++i)
        s.steps[(size_t)i] = juce::roundToInt(pJoy[(size_t)i]->load());
    s.lowpassKHz = joyValue(s.program, Joy::lowpass, s.steps[0]).value;
    s.delayMs = joyValue(s.program, Joy::delay, s.steps[1]).value;
    s.amplPercent = joyValue(s.program, Joy::ampl, s.steps[2]).value;
    s.lfX = joyValue(s.program, Joy::lf, s.steps[3]).value;
    s.rt = joyValue(s.program, Joy::rt, s.steps[4]).value;
    s.hfX = joyValue(s.program, Joy::hf, s.steps[5]).value;
    for (auto& e : engines)
        e.setSettings(s);

    inGain.setTargetValue(gainLinear(*pIn));
    outGain.setTargetValue(gainLinear(*pOut));
    mixAmount.setTargetValue(*pMix);
    wetOn.setTargetValue(*pPower >= 0.5f ? 1.0f : 0.0f);
    transformerOn.setTargetValue(*pTransformer >= 0.5f ? 1.0f : 0.0f);
    // The pan places the dry signal; constant power with unity at the centre.
    const auto angle = (juce::jlimit(-1.0f, 1.0f, pPan->load()) + 1.0f) * juce::MathConstants<float>::pi * 0.25f;
    panLeft.setTargetValue(std::cos(angle) * juce::MathConstants<float>::sqrt2);
    panRight.setTargetValue(std::sin(angle) * juce::MathConstants<float>::sqrt2);

    auto* monoL = scratch.getWritePointer(0);
    auto* monoR = scratch.getWritePointer(1);
    auto* wetL0 = scratch.getWritePointer(2);
    auto* wetR0 = scratch.getWritePointer(3);
    auto* wetL1 = scratch.getWritePointer(4);
    auto* wetR1 = scratch.getWritePointer(5);
    auto* xfWeight = scratch.getWritePointer(6);
    const auto* inL = buffer.getReadPointer(0);
    const auto* inR = buffer.getReadPointer(stereoIn ? 1 : 0);

    // The unit has one mono input: sum (0.5 times, so correlated channels keep their level rather than gaining
    // 6 dB), apply INPUT, and read the level there. Then the optional input transformer. The transformer always
    // runs so its state is warm when switched on; off, its crossfade weight is exactly 0 and the signal is unchanged.
    float peak = 0.0f, sumSq = 0.0f;
    for (int i = 0; i < n; ++i)
    {
        const auto g = inGain.getNextValue();
        const auto w = transformerOn.getNextValue();
        const auto l = inL[i], r = inR[i];
        const auto meterInput = monoSumGain * (l + r) * g;
        peak = std::max(peak, std::abs(meterInput));
        sumSq += meterInput * meterInput;
        xfWeight[i] = w;
        if (trueStereo)
        {
            monoL[i] = l * g;
            monoR[i] = r * g;
            monoL[i] += w * (transformers[0].processInput(monoL[i]) - monoL[i]);
            monoR[i] += w * (transformers[1].processInput(monoR[i]) - monoR[i]);
        }
        else
        {
            monoL[i] = meterInput;
            monoL[i] += w * (transformers[0].processInput(monoL[i]) - monoL[i]);
        }
    }
    const auto rms = std::sqrt(sumSq / (float)n);
    const auto a = 1.0f - std::exp(-(float)n / (float)vuTauSamples);
    vu = a * rms + (1.0f - a) * vu;
    inputVu.store(vu, std::memory_order_relaxed);
    if (peak >= 0.999f)
        inputClip.store(true, std::memory_order_relaxed);

    engines[0].process(monoL, wetL0, wetR0, n);
    if (trueStereo)
        engines[1].process(monoR, wetL1, wetR1, n);

    // Mono out has no read-out of its own: use the left side of the reverb. (Revisit when the algorithm's output
    // taps are implemented, docs/model_s_results.md.)
    const float* wetL = wetL0;
    const float* wetR = trueStereo ? wetR1 : wetR0;

    // The transformer stage's outputs: the captures' low band and the output transformers. In true stereo each
    // channel has its own transformer, and engine 0's left / engine 1's right are the two sides used.
    if (trueStereo)
    {
        for (int i = 0; i < n; ++i)
        {
            const auto w = xfWeight[i];
            const auto xl = transformers[0].processOutput(0, wetL0[i]);
            const auto xr = transformers[1].processOutput(1, wetR1[i]);
            wetL0[i] += w * (xl - wetL0[i]);
            wetR1[i] += w * (xr - wetR1[i]);
        }
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            const auto w = xfWeight[i];
            const auto xl = transformers[0].processOutput(0, wetL0[i]);
            const auto xr = transformers[0].processOutput(1, wetR0[i]);
            wetL0[i] += w * (xl - wetL0[i]);
            wetR0[i] += w * (xr - wetR0[i]);
        }
    }

    auto* outL = buffer.getWritePointer(0);
    auto* outR = stereoOut ? buffer.getWritePointer(1) : nullptr;
    for (int i = 0; i < n; ++i)
    {
        const auto mixV = mixAmount.getNextValue();
        const auto og = outGain.getNextValue();
        const auto on = wetOn.getNextValue();
        const auto pl = panLeft.getNextValue(), pr = panRight.getNextValue();
        const auto dryMono = monoSumGain * (inL[i] + inR[i]);
        // True stereo keeps the dry signal's channels; the original mode pans the mono sum as before.
        const auto dryL = trueStereo ? inL[i] * pl : dryMono * pl;
        const auto dryR = trueStereo ? inR[i] * pr : dryMono * pr;
        if (stereoOut)
        {
            const auto processedL = ((1.0f - mixV) * dryL + mixV * wetL[i]) * og;
            const auto processedR = ((1.0f - mixV) * dryR + mixV * wetR[i]) * og;
            const auto bypassL = inL[i], bypassR = inR[i]; // power off: the host's signal untouched
            outL[i] = on * processedL + (1.0f - on) * bypassL;
            outR[i] = on * processedR + (1.0f - on) * bypassR;
        }
        else
        {
            // A mono output cannot carry the pan, so the dry is the mono sum.
            const auto processed = ((1.0f - mixV) * dryMono + mixV * wetL[i]) * og;
            outL[i] = on * processed + (1.0f - on) * inL[i];
        }
    }
}

juce::AudioProcessorEditor* Rev246Processor::createEditor()
{
    return new Rev246Editor(*this);
}

void Rev246Processor::getStateInformation(juce::MemoryBlock& dest)
{
    auto state = apvts.copyState();
    auto ui = state.getOrCreateChildWithName(UiProps::node, nullptr);
    ui.setProperty("preset", presets.currentName(), nullptr);
    if (auto xml = state.createXml())
        copyXmlToBinary(*xml, dest);
}

void Rev246Processor::setStateInformation(const void* data, int size)
{
    if (auto xml = getXmlFromBinary(data, size))
        if (xml->hasTagName(apvts.state.getType()))
        {
            const auto tree = juce::ValueTree::fromXml(*xml);
            apvts.replaceState(tree);
            presets.setCurrentName(uiState().getProperty("preset", "Init S").toString());
        }
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Rev246Processor();
}
