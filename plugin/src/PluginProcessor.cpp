#include "PluginProcessor.h"

#include "PluginEditor.h"

namespace
{
constexpr double rampSeconds = 0.02;
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
    for (int i = 0; i < rv::params::numJoys; ++i)
        pJoy[(size_t)i] = apvts.getRawParameterValue(rv::params::joyId((rv::params::Joy)i));
}

bool Rev246Processor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    // Mono or stereo in (the unit has one mono input), stereo out.
    const auto& in = layouts.getMainInputChannelSet();
    return layouts.getMainOutputChannelSet() == juce::AudioChannelSet::stereo() &&
           (in == juce::AudioChannelSet::stereo() || in == juce::AudioChannelSet::mono());
}

void Rev246Processor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    engine.prepare(sampleRate, samplesPerBlock);
    for (auto* s : {&inGain, &outGain, &mixAmount, &panLeft, &panRight, &wetOn})
        s->reset(sampleRate, rampSeconds);
    inGain.setCurrentAndTargetValue(rv::params::gainLinear(*pIn));
    outGain.setCurrentAndTargetValue(rv::params::gainLinear(*pOut));
    mixAmount.setCurrentAndTargetValue(*pMix);
    wetOn.setCurrentAndTargetValue(*pPower >= 0.5f ? 1.0f : 0.0f);
    scratch.setSize(3, samplesPerBlock, false, false, true);
    setLatencySamples(0);
}

void Rev246Processor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    using namespace rv::params;
    const auto n = buffer.getNumSamples();
    const auto inChannels = getTotalNumInputChannels();
    if (n > scratch.getNumSamples())
        scratch.setSize(3, n, false, false, true); // only if a host exceeds its promised block size

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
    engine.setSettings(s);

    inGain.setTargetValue(gainLinear(*pIn));
    outGain.setTargetValue(gainLinear(*pOut));
    mixAmount.setTargetValue(*pMix);
    wetOn.setTargetValue(*pPower >= 0.5f ? 1.0f : 0.0f);
    // The pan places the dry signal; constant power with unity at the centre.
    const auto angle = (juce::jlimit(-1.0f, 1.0f, pPan->load()) + 1.0f) * juce::MathConstants<float>::pi * 0.25f;
    panLeft.setTargetValue(std::cos(angle) * juce::MathConstants<float>::sqrt2);
    panRight.setTargetValue(std::sin(angle) * juce::MathConstants<float>::sqrt2);

    auto* mono = scratch.getWritePointer(0);
    auto* wetL = scratch.getWritePointer(1);
    auto* wetR = scratch.getWritePointer(2);
    const auto* inL = buffer.getReadPointer(0);
    const auto* inR = buffer.getReadPointer(inChannels > 1 ? 1 : 0);

    // The unit has one mono input: sum, apply INPUT, and read the level there.
    float peak = 0.0f;
    for (int i = 0; i < n; ++i)
    {
        mono[i] = 0.5f * (inL[i] + inR[i]) * inGain.getNextValue();
        peak = std::max(peak, std::abs(mono[i]));
    }
    float prior = inputPeak.load(std::memory_order_relaxed);
    while (peak > prior && !inputPeak.compare_exchange_weak(prior, peak)) {}

    engine.process(mono, wetL, wetR, n);

    auto* outL = buffer.getWritePointer(0);
    auto* outR = buffer.getWritePointer(1);
    for (int i = 0; i < n; ++i)
    {
        const auto dryIn = 0.5f * (inL[i] + inR[i]);
        const auto mixV = mixAmount.getNextValue();
        const auto og = outGain.getNextValue();
        const auto on = wetOn.getNextValue();
        const auto dryL = dryIn * panLeft.getNextValue(), dryR = dryIn * panRight.getNextValue();
        const auto processedL = ((1.0f - mixV) * dryL + mixV * wetL[i]) * og;
        const auto processedR = ((1.0f - mixV) * dryR + mixV * wetR[i]) * og;
        const auto bypassL = inL[i], bypassR = inR[i]; // power off: the host's signal untouched
        outL[i] = on * processedL + (1.0f - on) * bypassL;
        outR[i] = on * processedR + (1.0f - on) * bypassR;
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
