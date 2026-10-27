#include "PluginProcessor.h"

#include "PluginEditor.h"
#include "ui/Design.h"

#include <algorithm>
#include <cmath>

namespace
{
constexpr double rampSeconds = 0.02; // 20 ms: the gain, mix, pan and crossfade ramps
constexpr float vuRiseSeconds = 0.3f; // VU: 99% of a step in ~300 ms
// Stereo-to-mono summation (original mode only): 0.5 = -6 dB, so two identical/correlated channels stay at unity
// instead of gaining 6 dB. The equal-power alternative, 1/sqrt(2), would keep uncorrelated stereo at its power but
// add 3 dB to centred material; input headroom matters more here, so use the conservative sum.
constexpr float monoSumGain = 0.5f;
constexpr float ovlPeak = 0.999f; // linear: an input sample at or above this lights OVL (full scale to within a code)
constexpr float switchOnThreshold = 0.5f; // a bool parameter's raw value is 0 or 1

// The channels of `Rev246Processor::scratch`.
enum ScratchChannel
{
    chEngineIn0,         // the mono feed of engine 0 (the summed input, or the left channel in true stereo)
    chEngineIn1,         // engine 1's feed (the right channel, true stereo only)
    chWetL0,             // engine 0's outputs
    chWetR0,
    chWetL1,             // engine 1's outputs
    chWetR1,
    chTransformerWeight, // the transformer stage's crossfade weight, per sample, shared by input and output
    chDryInL,            // the input, delayed by the reported latency
    chDryInR,
    numScratchChannels
};

bool isOn(const std::atomic<float>& rawValue)
{
    return rawValue.load(std::memory_order_relaxed) >= switchOnThreshold;
}

// The pan places the dry signal: constant power, with unity gain (not -3 dB) at the centre.
struct PanGains
{
    float left, right;
};

PanGains panGains(float pan)
{
    const auto angle = (juce::jlimit(-1.0f, 1.0f, pan) + 1.0f) * juce::MathConstants<float>::pi * 0.25f;
    return {std::cos(angle) * juce::MathConstants<float>::sqrt2, std::sin(angle) * juce::MathConstants<float>::sqrt2};
}
} // namespace

// Typed views of the scratch buffer for one chunk.
struct Rev246Processor::BlockBuffers
{
    explicit BlockBuffers(juce::AudioBuffer<float>& s)
        : engineIn0(s.getWritePointer(chEngineIn0)),
          engineIn1(s.getWritePointer(chEngineIn1)),
          wetL0(s.getWritePointer(chWetL0)),
          wetR0(s.getWritePointer(chWetR0)),
          wetL1(s.getWritePointer(chWetL1)),
          wetR1(s.getWritePointer(chWetR1)),
          transformerWeight(s.getWritePointer(chTransformerWeight)),
          dryInL(s.getWritePointer(chDryInL)),
          dryInR(s.getWritePointer(chDryInR))
    {
    }

    float *engineIn0, *engineIn1, *wetL0, *wetR0, *wetL1, *wetR1, *transformerWeight, *dryInL, *dryInR;
};

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
    pSixteenBit = apvts.getRawParameterValue(id::sixteenBit);
    pStereoWidth = apvts.getRawParameterValue(id::stereoWidth);
    for (int i = 0; i < rv::params::numJoys; ++i)
        pJoy[(size_t)i] = apvts.getRawParameterValue(rv::params::joyId((rv::params::Joy)i));
    applyMeterCal(rv::design::defaultMeterCal);
}

void Rev246Processor::applyMeterCal(int index)
{
    index = juce::jlimit(0, rv::design::numMeterCal - 1, index);
    meterCal.store(index, std::memory_order_relaxed);
    adcFullScale.store(rv::dsp::AdcClipper::fullScaleForZeroVu(rv::design::meterCalDbFs[index]), std::memory_order_relaxed);
}

void Rev246Processor::setMeterCal(int index)
{
    applyMeterCal(index);
    uiState().setProperty(UiProps::meterCal, getMeterCal(), nullptr);
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
    const auto maxBlock = juce::jmax(1, samplesPerBlock);
    for (auto& e : engines)
        e.prepare(sampleRate, maxBlock);
    for (auto& t : transformers)
        t.prepare(sampleRate);
    stereoWidener.prepare(sampleRate);
    stereoWidener.setWidth(pStereoWidth->load(), true);
    for (auto& c : clippers)
        c.prepare(sampleRate);
    for (auto* s : {&inGain, &outGain, &mixAmount, &panLeft, &panRight, &powerOn, &transformerOn, &sixteenBitOn})
        s->reset(sampleRate, rampSeconds);
    inGain.setCurrentAndTargetValue(rv::params::gainLinearFromPosition(pIn->load()));
    outGain.setCurrentAndTargetValue(rv::params::gainLinearFromPosition(pOut->load()));
    mixAmount.setCurrentAndTargetValue(pMix->load());
    const auto pan = panGains(pPan->load());
    panLeft.setCurrentAndTargetValue(pan.left); // else the dry would fade in from 0 over the first 20 ms
    panRight.setCurrentAndTargetValue(pan.right);
    powerOn.setCurrentAndTargetValue(isOn(*pPower) ? 1.0f : 0.0f);
    transformerOn.setCurrentAndTargetValue(isOn(*pTransformer) ? 1.0f : 0.0f);
    sixteenBitOn.setCurrentAndTargetValue(isOn(*pSixteenBit) ? 1.0f : 0.0f);
    // processBlock works in chunks of at most this many samples, so it never has to grow the scratch buffer.
    scratch.setSize(numScratchChannels, maxBlock, false, false, true);
    // VU meter: a first-order RMS integrator. The classic VU reaches 99% in ~300 ms, so tau = 300 ms / ln 100.
    vuTauSamples = (double)(vuRiseSeconds / std::log(100.0f)) * sampleRate;
    vuLevel = 0.0f;
    inputVu.store(0.0f);
    inputClip.store(false);
    lastTrueStereo = isStereoToStereo() && isOn(*pTrueStereo);
    engineKey = {}; // the engines are prepared afresh: send them their settings again
    const auto latency = clippers[0].latency();
    for (auto& d : dryDelay)
        d.assign((size_t)latency, 0.0f);
    dryPos = 0;
    setLatencySamples(latency);
}

void Rev246Processor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    const auto numSamples = buffer.getNumSamples();
    const auto chunkCapacity = scratch.getNumSamples();
    if (numSamples <= 0 || buffer.getNumChannels() <= 0)
        return; // some hosts send empty blocks; there is nothing to process or to meter
    if (chunkCapacity <= 0)
    {
        buffer.clear(); // not prepared
        return;
    }

    // A host should never exceed the block size it promised in prepareToPlay; if one does, process the block in
    // pieces rather than growing the scratch buffer here. The pieces are views of `buffer` (no allocation).
    for (int start = 0; start < numSamples; start += chunkCapacity)
    {
        juce::AudioBuffer<float> chunk(buffer.getArrayOfWritePointers(), buffer.getNumChannels(), start,
                                       juce::jmin(chunkCapacity, numSamples - start));
        processChunk(chunk);
    }
}

void Rev246Processor::processChunk(juce::AudioBuffer<float>& buffer)
{
    const auto n = buffer.getNumSamples();
    const bool stereoIn = getTotalNumInputChannels() > 1, stereoOut = getTotalNumOutputChannels() > 1;
    // The original mode (the default, and the only option outside stereo in / stereo out) sums the input to the
    // unit's mono input. True stereo runs one engine per input channel.
    const bool trueStereo = stereoIn && stereoOut && isOn(*pTrueStereo);

    // A new routing has no usable state in the second engine: start it silent rather than with an old tail.
    if (trueStereo != lastTrueStereo)
    {
        engines[1].reset();
        transformers[1].reset();
        clippers[1].reset();
        engineKey = {}; // a reset may drop the settings: send them again
        lastTrueStereo = trueStereo;
    }

    applyParameters();

    const BlockBuffers b(scratch);
    processInputStage(buffer.getReadPointer(0), buffer.getReadPointer(stereoIn ? 1 : 0), b, n, trueStereo);

    engines[0].process(b.engineIn0, b.wetL0, b.wetR0, n);
    if (trueStereo)
        engines[1].process(b.engineIn1, b.wetL1, b.wetR1, n);

    // Mono out uses the left side of the reverb for now; with the real read-outs it becomes the (L + R) / sqrt(2)
    // fold-down (docs/DEVELOPMENT.md, Channel modes). In true stereo the two sides come from the two engines.
    processTransformerOutputs(b, n, trueStereo, stereoOut);

    // Stereo width on the wet path. The late-tail decorrelator from the design is not in the chain yet: the placeholder
    // engine has no early/late split, and its allpass delay (~8 ms) would also land on the early part. When the engine
    // exposes the late tail, run rv::dsp::LateDecorrelator on that tail before this widener (DEVELOPMENT.md, Open
    // items).
    stereoWidener.setWidth(pStereoWidth->load());
    if (stereoOut)
        stereoWidener.process(b.wetL0, trueStereo ? b.wetR1 : b.wetR0, n);

    mixOutput(buffer, b, trueStereo, stereoOut);
}

// Panel settings -> the engines' settings and the audio-rate ramps' targets.
void Rev246Processor::applyParameters()
{
    using namespace rv::params;

    EngineKey key;
    key.program = juce::jlimit(0, numPrograms - 1, juce::roundToInt(pProgram->load()));
    key.sixteenBit = isOn(*pSixteenBit);
    for (int i = 0; i < numJoys; ++i)
        key.steps[(size_t)i] = juce::jlimit(0, joyMaxStep((Joy)i), juce::roundToInt(pJoy[(size_t)i]->load()));
    if (!(key == engineKey)) // the settings follow from the key alone; the engines only hear about a change
    {
        engineKey = key;
        const auto program = (Program)key.program;
        auto stepOf = [&key](Joy j) { return key.steps[(size_t)j]; };
        rv::dsp::Settings s;
        s.program = program;
        s.steps = key.steps;
        s.lowpassKHz = joyNumber(program, Joy::lowpass, stepOf(Joy::lowpass));
        s.delayMs = joyNumber(program, Joy::delay, stepOf(Joy::delay));
        s.amplPercent = joyNumber(program, Joy::ampl, stepOf(Joy::ampl));
        s.lfX = joyNumber(program, Joy::lf, stepOf(Joy::lf));
        s.rt = joyNumber(program, Joy::rt, stepOf(Joy::rt));
        s.hfX = joyNumber(program, Joy::hf, stepOf(Joy::hf));
        s.sixteenBit = key.sixteenBit;
        for (auto& e : engines)
            e.setSettings(s);
    }

    inGain.setTargetValue(gainLinearFromPosition(pIn->load()));
    outGain.setTargetValue(gainLinearFromPosition(pOut->load()));
    mixAmount.setTargetValue(pMix->load());
    powerOn.setTargetValue(isOn(*pPower) ? 1.0f : 0.0f);
    transformerOn.setTargetValue(isOn(*pTransformer) ? 1.0f : 0.0f);
    sixteenBitOn.setTargetValue(key.sixteenBit ? 1.0f : 0.0f);
    const auto fullScale = adcFullScale.load(std::memory_order_relaxed);
    for (auto& c : clippers)
        c.setFullScale(fullScale);
    for (auto& t : transformers)
        t.setFullScale(fullScale);
    const auto pan = panGains(pPan->load());
    panLeft.setTargetValue(pan.left);
    panRight.setTargetValue(pan.right);
}

// The unit has one mono input: sum (monoSumGain times), apply INPUT, and read the level there. Then the optional input
// transformer, then the A/D: in 16-bit mode its clipping (AdcClipper). Both stages always run so their state is warm
// when switched on, and the A/D also feeds the OVL LED in float mode; off, a stage's crossfade weight is exactly 0.
// The A/D stage delays the wet path by the plugin's latency in either mode; the dry and bypassed signals get the same
// delay here (b.dryIn*).
void Rev246Processor::processInputStage(const float* inL, const float* inR, const BlockBuffers& b, int n,
                                        bool trueStereo)
{
    float peak = 0.0f, sumSq = 0.0f;
    const auto ringSize = dryDelay[0].size();
    auto ringPos = dryPos;
    for (int i = 0; i < n; ++i)
    {
        const auto g = inGain.getNextValue();
        const auto transformerMix = transformerOn.getNextValue();
        const auto q = sixteenBitOn.getNextValue();
        const auto l = inL[i], r = inR[i];
        b.dryInL[i] = dryDelay[0][ringPos];
        b.dryInR[i] = dryDelay[1][ringPos];
        dryDelay[0][ringPos] = l;
        dryDelay[1][ringPos] = r;
        ringPos = ringPos + 1 == ringSize ? 0 : ringPos + 1;
        const auto monoInput = monoSumGain * (l + r) * g;
        peak = std::max(peak, std::abs(monoInput));
        sumSq += monoInput * monoInput;
        b.transformerWeight[i] = transformerMix;
        if (trueStereo)
        {
            auto left = l * g, right = r * g;
            left += transformerMix * (transformers[0].processInput(left) - left);
            right += transformerMix * (transformers[1].processInput(right) - right);
            b.engineIn0[i] = clippers[0].process(left, q);
            b.engineIn1[i] = clippers[1].process(right, q);
        }
        else
        {
            auto mono = monoInput;
            mono += transformerMix * (transformers[0].processInput(mono) - mono);
            b.engineIn0[i] = clippers[0].process(mono, q);
        }
    }
    dryPos = ringPos;

    const auto rms = std::sqrt(sumSq / (float)n);
    const auto a = 1.0f - std::exp(-(float)n / (float)vuTauSamples);
    vuLevel = a * rms + (1.0f - a) * vuLevel;
    inputVu.store(vuLevel, std::memory_order_relaxed);
    // OVL: a full-scale input, or the A/D clipping (the unit's REG LED; in float mode, where it would clip).
    const bool over0 = clippers[0].takeOverload(), over1 = clippers[1].takeOverload();
    if (peak >= ovlPeak || over0 || (trueStereo && over1))
        inputClip.store(true, std::memory_order_relaxed);
}

// The transformer stage's outputs: the captures' low band and the output transformers, crossfaded with the dry-through
// wet by the stage's weight. In true stereo each channel has its own transformer, and engine 0's left / engine 1's
// right are the two sides used. With a mono output only the left side exists, and the right is left alone.
void Rev246Processor::processTransformerOutputs(const BlockBuffers& b, int n, bool trueStereo, bool stereoOut)
{
    auto& transformerL = transformers[0];
    auto& transformerR = trueStereo ? transformers[1] : transformers[0];
    float* wetL = b.wetL0;
    float* wetR = trueStereo ? b.wetR1 : b.wetR0;
    for (int i = 0; i < n; ++i)
    {
        const auto w = b.transformerWeight[i];
        const auto xl = transformerL.processOutput(0, wetL[i]);
        wetL[i] += w * (xl - wetL[i]);
        if (stereoOut)
        {
            const auto xr = transformerR.processOutput(1, wetR[i]);
            wetR[i] += w * (xr - wetR[i]);
        }
    }
}

// Dry/wet mix, OUTPUT, and the power crossfade into the bypassed (host) signal, written over `buffer`.
void Rev246Processor::mixOutput(juce::AudioBuffer<float>& buffer, const BlockBuffers& b, bool trueStereo,
                                bool stereoOut)
{
    const float* wetL = b.wetL0;
    const float* wetR = trueStereo ? b.wetR1 : b.wetR0;
    const auto n = buffer.getNumSamples();
    auto* outL = buffer.getWritePointer(0);
    auto* outR = stereoOut ? buffer.getWritePointer(1) : nullptr;
    for (int i = 0; i < n; ++i)
    {
        const auto mix = mixAmount.getNextValue();
        const auto outputGain = outGain.getNextValue();
        const auto power = powerOn.getNextValue();
        const auto panL = panLeft.getNextValue(), panR = panRight.getNextValue();
        const auto dryMono = monoSumGain * (b.dryInL[i] + b.dryInR[i]);
        if (stereoOut)
        {
            // True stereo keeps the dry signal's channels; the original mode pans the mono sum as before.
            const auto dryL = trueStereo ? b.dryInL[i] * panL : dryMono * panL;
            const auto dryR = trueStereo ? b.dryInR[i] * panR : dryMono * panR;
            const auto processedL = ((1.0f - mix) * dryL + mix * wetL[i]) * outputGain;
            const auto processedR = ((1.0f - mix) * dryR + mix * wetR[i]) * outputGain;
            // Power off: the host's signal, delayed only.
            outL[i] = power * processedL + (1.0f - power) * b.dryInL[i];
            outR[i] = power * processedR + (1.0f - power) * b.dryInR[i];
        }
        else
        {
            // A mono output cannot carry the pan, so the dry is the mono sum.
            const auto processed = ((1.0f - mix) * dryMono + mix * wetL[i]) * outputGain;
            outL[i] = power * processed + (1.0f - power) * b.dryInL[i];
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
    ui.setProperty(UiProps::preset, presets.currentName(), nullptr);
    if (auto xml = state.createXml())
        copyXmlToBinary(*xml, dest);
}

void Rev246Processor::setStateInformation(const void* data, int size)
{
    const auto xml = getXmlFromBinary(data, size);
    if (xml == nullptr || !xml->hasTagName(apvts.state.getType()))
        return; // not our state (empty, garbage or another plugin's): keep the current settings

    const auto tree = juce::ValueTree::fromXml(*xml);
    apvts.replaceState(tree);
    const auto ui = tree.getChildWithName(UiProps::node); // may be absent (older sessions): defaults then
    applyMeterCal((int)ui.getProperty(UiProps::meterCal, rv::design::defaultMeterCal));
    presets.setCurrentName(ui.getProperty(UiProps::preset, rv::PresetManager::factoryName(0)).toString());
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Rev246Processor();
}
