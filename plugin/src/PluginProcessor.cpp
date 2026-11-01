#include "PluginProcessor.h"

#include "PluginEditor.h"
#include "dsp/AdcCalibration.h"
#include "ui/Design.h"

#include <algorithm>
#include <cmath>

namespace
{
constexpr double rampSeconds = 0.02; // 20 ms: the gain, mix, pan and crossfade ramps
constexpr float vuRiseSeconds = 0.3f; // VU: 99% of a step in ~300 ms (ANSI C16.5 / IEC 60268-17)
// Two cascaded one-poles of equal tau are critically damped: a step reaches 1 - (1 + x) e^-x of its value at x = t / tau,
// which is 99% at x = 6.64. That is the slow start and the unhurried fall of a VU needle, which a single pole lacks.
constexpr float vuCriticalDampingX99 = 6.64f;
// Stereo-to-mono summation (original mode only): 0.5 = -6 dB, so two identical/correlated channels stay at unity
// instead of gaining 6 dB. The equal-power alternative, 1/sqrt(2), would keep uncorrelated stereo at its power but
// add 3 dB to centred material; input headroom matters more here, so use the conservative sum.
constexpr float monoSumGain = 0.5f;
constexpr float monoFoldGain = 0.70710678f; // 1 / sqrt(2): the mono output's (L + R) fold-down of the wet
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
    chDuckKey,           // the dry input level that keys the ducker (one value per sample)
    chDryInL,            // the input, delayed by the reported latency
    chDryInR,
    chDryGain,           // the INPUT gain that went with each delayed dry sample
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
          duckKey(s.getWritePointer(chDuckKey)),
          dryInL(s.getWritePointer(chDryInL)),
          dryInR(s.getWritePointer(chDryInR)),
          dryGain(s.getWritePointer(chDryGain))
    {
    }

    float *engineIn0, *engineIn1, *wetL0, *wetR0, *wetL1, *wetR1, *transformerWeight, *duckKey, *dryInL, *dryInR, *dryGain;
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
    pDuckEnabled = apvts.getRawParameterValue(id::duckEnabled);
    pDuckRangeDb = apvts.getRawParameterValue(id::duckRangeDb);
    pDuckThresholdDb = apvts.getRawParameterValue(id::duckThresholdDb);
    pDuckAttackMs = apvts.getRawParameterValue(id::duckAttackMs);
    pDuckReleaseMs = apvts.getRawParameterValue(id::duckReleaseMs);
    for (int i = 0; i < rv::params::numJoys; ++i)
        pJoy[(size_t)i] = apvts.getRawParameterValue(rv::params::joyId((rv::params::Joy)i));
    applyMeterCal(rv::design::defaultMeterCal);
}

void Rev246Processor::applyMeterCal(int index)
{
    index = juce::jlimit(0, rv::design::numMeterCal - 1, index);
    meterCal.store(index, std::memory_order_relaxed);
    adcFullScale.store(rv::dsp::adcFullScaleForZeroVu(rv::design::meterCalDbFs[index]), std::memory_order_relaxed);
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
    ducker.prepare(sampleRate);
    for (auto* s : {&inGain, &outGain, &mixAmount, &panLeft, &panRight, &powerOn, &transformerOn})
        s->reset(sampleRate, rampSeconds);
    inGain.setCurrentAndTargetValue(rv::params::gainLinearFromPosition(pIn->load()));
    outGain.setCurrentAndTargetValue(rv::params::gainLinearFromPosition(pOut->load()));
    mixAmount.setCurrentAndTargetValue(pMix->load());
    const auto pan = panGains(pPan->load());
    panLeft.setCurrentAndTargetValue(pan.left); // else the dry would fade in from 0 over the first 20 ms
    panRight.setCurrentAndTargetValue(pan.right);
    powerOn.setCurrentAndTargetValue(isOn(*pPower) ? 1.0f : 0.0f);
    transformerOn.setCurrentAndTargetValue(isOn(*pTransformer) ? 1.0f : 0.0f);
    routeGain.reset(sampleRate, 0.5 * rampSeconds);
    routeGain.setCurrentAndTargetValue(1.0f);
    dryStereo.reset(sampleRate, rampSeconds);
    // processBlock works in chunks of at most this many samples, so it never has to grow the scratch buffer.
    scratch.setSize(numScratchChannels, maxBlock, false, false, true);
    // VU meter: a critically damped (two-pole) RMS integrator that reaches 99% of a step in 300 ms, like the standard.
    vuTauSamples = (double)(vuRiseSeconds / vuCriticalDampingX99) * sampleRate;
    vuLevel = vuStage1 = 0.0f;
    inputVu.store(0.0f);
    inputClip.store(false);
    lastTrueStereo = isStereoToStereo() && isOn(*pTrueStereo);
    dryStereo.setCurrentAndTargetValue(lastTrueStereo ? 1.0f : 0.0f);
    engineKey = {}; // the engines are prepared afresh: send them their settings, and start in the current program
    applyParameters(); // rather than crossfading into it from the default one
    ducker.reset();
    for (auto& e : engines)
        e.reset();
    const auto latency = engines[0].latency(); // the same for every program, both modes and both engines
    for (size_t k = 0; k < dryDelay.size(); ++k)
        dryDelay[k].assign((size_t)latency, k == 2 ? 1.0f : 0.0f);
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
    const bool wantTrueStereo = stereoIn && stereoOut && isOn(*pTrueStereo);

    // A routing change dips the wet to silence first (routeGain) and happens at the chunk where the dip has ended; the
    // dry crossfades to its new form at once (dryStereo). A new routing has no usable state in the second engine: start
    // it silent rather than with an old tail.
    dryStereo.setTargetValue(wantTrueStereo ? 1.0f : 0.0f);
    if (wantTrueStereo == lastTrueStereo)
        routeGain.setTargetValue(1.0f); // (also: switched back during the dip)
    else if (routeGain.getTargetValue() > 0.0f)
        routeGain.setTargetValue(0.0f);
    else if (!routeGain.isSmoothing())
    {
        engines[1].reset();
        transformers[1].reset();
        engineKey = {}; // a reset may drop the settings: send them again
        lastTrueStereo = wantTrueStereo;
        routeGain.setTargetValue(1.0f);
    }
    const bool trueStereo = lastTrueStereo;

    applyParameters();

    const BlockBuffers b(scratch);
    const auto inputPeak = processInputStage(buffer.getReadPointer(0), buffer.getReadPointer(stereoIn ? 1 : 0), b, n,
                                             trueStereo);

    engines[0].process(b.engineIn0, b.wetL0, b.wetR0, n);
    if (trueStereo)
        engines[1].process(b.engineIn1, b.wetL1, b.wetR1, n);

    // OVL: a full-scale input, or the A/D out of range (the unit's REG LED; in float mode, where it would clip).
    const bool over0 = engines[0].takeOverload(), over1 = engines[1].takeOverload();
    if (inputPeak >= ovlPeak || over0 || (trueStereo && over1))
        inputClip.store(true, std::memory_order_relaxed);

    // Both read-outs through the output stage; a mono output folds them down after it (mixOutput). In true stereo the
    // two sides come from the two engines.
    processTransformerOutputs(b, n, trueStereo);

    // The input-keyed duck (its key is the undelayed input, so it leads the wet by the latency, ~1 ms: a small lookahead).
    // Width comes last (pure M/S gains: nothing moves in time).
    auto* wetR = trueStereo ? b.wetR1 : b.wetR0;
    ducker.process(b.duckKey, b.wetL0, wetR, n);

    stereoWidener.setWidth(pStereoWidth->load());
    if (stereoOut)
        stereoWidener.process(b.wetL0, wetR, n);

    mixOutput(buffer, b, trueStereo, stereoOut);
}

// Panel settings -> the engines' settings and the audio-rate ramps' targets.
void Rev246Processor::applyParameters()
{
    using namespace rv::params;

    EngineKey key;
    key.program = juce::jlimit(0, numPrograms - 1, juce::roundToInt(pProgram->load()));
    key.sixteenBit = isOn(*pSixteenBit);
    for (int i = 0; i < numJoys; ++i) // a step past the program's last (N's longer HF column) reads as the last
        key.steps[(size_t)i] = juce::jlimit(0, joyMaxStep((Program)key.program, (Joy)i),
                                            juce::roundToInt(pJoy[(size_t)i]->load()));
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
        engines[0].setSettings(s);
        if (lastTrueStereo) // engine 1 idles otherwise; a switch to true stereo resets it and sends them again
            engines[1].setSettings(s);
    }

    inGain.setTargetValue(gainLinearFromPosition(pIn->load()));
    outGain.setTargetValue(gainLinearFromPosition(pOut->load()));
    mixAmount.setTargetValue(pMix->load());
    powerOn.setTargetValue(isOn(*pPower) ? 1.0f : 0.0f);
    transformerOn.setTargetValue(isOn(*pTransformer) ? 1.0f : 0.0f);
    ducker.setEnabled(isOn(*pDuckEnabled));
    ducker.setRangeDb(pDuckRangeDb->load());
    ducker.setThresholdDb(pDuckThresholdDb->load());
    ducker.setAttackMs(pDuckAttackMs->load());
    ducker.setReleaseMs(pDuckReleaseMs->load());
    const auto fullScale = adcFullScale.load(std::memory_order_relaxed);
    for (auto& e : engines)
        e.setFullScale(fullScale);
    for (auto& t : transformers)
        t.setFullScale(fullScale);
    const auto pan = panGains(pPan->load());
    panLeft.setTargetValue(pan.left);
    panRight.setTargetValue(pan.right);
}

// The unit has one mono input: sum (monoSumGain times), apply INPUT, and read the level there. Then the optional input
// transformer; the A/D (its clip in 16-bit mode and the overload flag) is inside the engine, at the device rate, after
// its analogue input chain. The transformer always runs so its state is warm when switched on; off, its crossfade
// weight is exactly 0. The engine delays the wet path by the plugin's latency; the dry and bypassed signals get the
// same delay here (b.dryIn*).
float Rev246Processor::processInputStage(const float* inL, const float* inR, const BlockBuffers& b, int n,
                                         bool trueStereo)
{
    float peak = 0.0f, sumSq = 0.0f;
    const auto ringSize = dryDelay[0].size();
    auto ringPos = dryPos;
    for (int i = 0; i < n; ++i)
    {
        const auto g = inGain.getNextValue();
        const auto transformerMix = transformerOn.getNextValue();
        const auto l = inL[i], r = inR[i];
        b.dryInL[i] = dryDelay[0][ringPos];
        b.dryInR[i] = dryDelay[1][ringPos];
        b.dryGain[i] = dryDelay[2][ringPos];
        dryDelay[0][ringPos] = l;
        dryDelay[1][ringPos] = r;
        dryDelay[2][ringPos] = g;
        ringPos = ringPos + 1 == ringSize ? 0 : ringPos + 1;
        const auto monoInput = monoSumGain * (l + r) * g;
        peak = std::max(peak, std::abs(monoInput));
        sumSq += monoInput * monoInput;
        b.duckKey[i] = std::max(std::abs(l * g), std::abs(r * g));
        b.transformerWeight[i] = transformerMix;
        if (trueStereo)
        {
            auto left = l * g, right = r * g;
            left += transformerMix * (transformers[0].processInput(left) - left);
            right += transformerMix * (transformers[1].processInput(right) - right);
            b.engineIn0[i] = left;
            b.engineIn1[i] = right;
        }
        else
        {
            auto mono = monoInput;
            mono += transformerMix * (transformers[0].processInput(mono) - mono);
            b.engineIn0[i] = mono;
        }
    }
    dryPos = ringPos;

    const auto rms = std::sqrt(sumSq / (float)n);
    const auto a = 1.0f - std::exp(-(float)n / (float)vuTauSamples);
    vuStage1 = a * rms + (1.0f - a) * vuStage1;
    vuLevel = a * vuStage1 + (1.0f - a) * vuLevel;
    inputVu.store(vuLevel, std::memory_order_relaxed);
    return peak;
}

// The transformer stage's outputs: the captures' low band and the output transformers, crossfaded with the dry-through
// wet by the stage's weight, on both read-outs (a mono output folds them down afterwards). In true stereo each
// channel has its own transformer, and engine 0's left / engine 1's right are the two sides used.
void Rev246Processor::processTransformerOutputs(const BlockBuffers& b, int n, bool trueStereo)
{
    auto& transformerL = transformers[0];
    auto& transformerR = trueStereo ? transformers[1] : transformers[0];
    float* wetL = b.wetL0;
    float* wetR = trueStereo ? b.wetR1 : b.wetR0;
    for (int i = 0; i < n; ++i)
    {
        const auto w = b.transformerWeight[i];
        const auto xl = transformerL.processOutput(0, wetL[i]);
        const auto xr = transformerR.processOutput(1, wetR[i]);
        wetL[i] += w * (xl - wetL[i]);
        wetR[i] += w * (xr - wetR[i]);
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
        const auto route = routeGain.getNextValue();
        // The dry goes through INPUT as well as OUTPUT, so a LINKed pair (one up by what the other is down) leaves it where it was.
        const auto dg = b.dryGain[i];
        const auto dryMono = monoSumGain * (b.dryInL[i] + b.dryInR[i]) * dg;
        if (stereoOut)
        {
            // True stereo keeps the dry signal's channels; the original mode pans the mono sum as before.
            const auto st = dryStereo.getNextValue(); // 0 or 1 at rest: exactly one form, no rounding from the blend
            auto blend = [st, dryMono](float own) { return st >= 1.0f ? own : st <= 0.0f ? dryMono : dryMono + st * (own - dryMono); };
            const auto dryL = blend(b.dryInL[i] * dg) * panL;
            const auto dryR = blend(b.dryInR[i] * dg) * panR;
            const auto processedL = ((1.0f - mix) * dryL + mix * route * wetL[i]) * outputGain;
            const auto processedR = ((1.0f - mix) * dryR + mix * route * wetR[i]) * outputGain;
            // Power off: the host's signal, delayed only.
            outL[i] = power * processedL + (1.0f - power) * b.dryInL[i];
            outR[i] = power * processedR + (1.0f - power) * b.dryInR[i];
        }
        else
        {
            // A mono output cannot carry the pan, so the dry is the mono sum. The wet is the unit's two read-outs folded
            // down at equal power: they are uncorrelated, so each keeps its level (DEVELOPMENT.md, Channel modes).
            const auto wetMono = (wetL[i] + wetR[i]) * monoFoldGain;
            const auto processed = ((1.0f - mix) * dryMono + mix * wetMono) * outputGain;
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
    state.setProperty(rv::params::stateVersionProperty, rv::params::stateVersion, nullptr);
    if (auto xml = state.createXml())
        copyXmlToBinary(*xml, dest);
}

void Rev246Processor::setStateInformation(const void* data, int size)
{
    const auto xml = getXmlFromBinary(data, size);
    if (xml == nullptr || !xml->hasTagName(apvts.state.getType()))
        return; // not our state (empty, garbage or another plugin's): keep the current settings

    const auto tree = juce::ValueTree::fromXml(*xml);
    rv::params::migrateState(tree); // an older build's session: up to the current layout first
    apvts.replaceState(tree);
    const auto ui = tree.getChildWithName(UiProps::node); // may be absent (older sessions): defaults then
    applyMeterCal((int)ui.getProperty(UiProps::meterCal, rv::design::defaultMeterCal));
    presets.setCurrentName(ui.getProperty(UiProps::preset, rv::PresetManager::factoryName(0)).toString());
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Rev246Processor();
}
