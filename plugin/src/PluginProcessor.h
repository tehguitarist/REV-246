#pragma once

#include "dsp/ReverbEngine.h"
#include "dsp/StereoWidth.h"
#include "dsp/Transformer.h"
#include "params/Parameters.h"
#include "presets/PresetManager.h"

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include <array>
#include <atomic>
#include <vector>

class Rev246Processor : public juce::AudioProcessor
{
  public:
    // Properties of the "UI" child of the state: kept with the session, not parameters.
    struct UiProps
    {
        static inline const juce::Identifier node{"UI"};
        static inline const juce::Identifier width{"width"};
        static inline const juce::Identifier tooltipsOn{"tooltipsOn"};
        static inline const juce::Identifier meterCal{"meterCal"}; // index into design::meterCalDbFs
        static inline const juce::Identifier preset{"preset"};     // the current preset's name, written on save
    };
    static constexpr int defaultEditorWidth = 1398; // 40% of the 3495-wide design space

    Rev246Processor();
    ~Rev246Processor() override = default;

    void prepareToPlay(double sampleRate, int samplesPerBlock) override;
    void releaseResources() override {}
    bool isBusesLayoutSupported(const BusesLayout&) const override;
    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    using AudioProcessor::processBlock;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

    const juce::String getName() const override { return JucePlugin_Name; }
    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }
    double getTailLengthSeconds() const override { return 21.0; } // L's longest RT plus the longest DELAY
    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram(int) override {}
    const juce::String getProgramName(int) override { return {}; }
    void changeProgramName(int, const juce::String&) override {}

    void getStateInformation(juce::MemoryBlock&) override;
    void setStateInformation(const void*, int) override;

    juce::AudioProcessorValueTreeState apvts;
    rv::PresetManager presets{apvts};

    // The panel's LEVEL meter, read by the editor's timer: the input (after the INPUT fader) on a VU ballistic,
    // linear, and whether a full-scale sample occurred since the last read.
    float getInputVu() const { return inputVu.load(std::memory_order_relaxed); }
    bool takeInputClip() { return inputClip.exchange(false, std::memory_order_relaxed); }

    // The meter calibration (an index into design::meterCalDbFs), kept in the UI state. It is the unit's input
    // alignment, so it also places the A/D's full scale: 6 dB over 0 VU for a sine (dsp/AdcCalibration.h).
    int getMeterCal() const { return meterCal.load(std::memory_order_relaxed); }
    void setMeterCal(int index);

    // True only for a stereo in / stereo out bus layout: the original / true-stereo choice has no meaning otherwise.
    bool isStereoToStereo() const
    {
        return getTotalNumInputChannels() > 1 && getTotalNumOutputChannels() > 1;
    }

    juce::ValueTree uiState() { return apvts.state.getOrCreateChildWithName(UiProps::node, nullptr); }

  private:
    // Typed views of `scratch` for one chunk of a block (defined in PluginProcessor.cpp).
    struct BlockBuffers;
    // What the engines' settings are derived from; they are only sent when it changes (program -1: nothing sent yet).
    struct EngineKey
    {
        int program = -1;
        bool sixteenBit = false;
        std::array<int, rv::params::numJoys> steps{};
        bool operator==(const EngineKey& o) const
        {
            return program == o.program && sixteenBit == o.sixteenBit && steps == o.steps;
        }
    };

    // The unit's algorithm is mono in / stereo out. The second engine and transformer make the optional true-stereo
    // mode (one mono path per input channel); otherwise engine 0 carries the (possibly summed) mono input.
    std::array<rv::dsp::ReverbEngine, 2> engines;
    std::array<rv::dsp::Transformer, 2> transformers;
    rv::dsp::StereoWidener stereoWidener;
    // 20 ms ramps (set in prepareToPlay). The two "On" ones are crossfade weights, 0 or 1 at rest: power (processed
    // against bypassed) and the transformer stage.
    juce::LinearSmoothedValue<float> inGain, outGain, mixAmount, panLeft, panRight, powerOn, transformerOn;
    // The ORIGINAL / TRUE ST switch: the wet dips to silence (10 ms), the routing changes there, and it comes back
    // (10 ms); a wet crossfade would need both routings at once, i.e. a third engine. The dry needs no engine, so it
    // crossfades between its two forms (the panned mono sum, the input's own channels) over 20 ms (dryStereo).
    juce::LinearSmoothedValue<float> routeGain, dryStereo;
    juce::AudioBuffer<float> scratch; // BlockBuffers' channels; sized in prepareToPlay, so a block never allocates
    std::atomic<float> inputVu{0.0f};
    std::atomic<bool> inputClip{false};
    std::atomic<int> meterCal{0};
    std::atomic<float> adcFullScale{1.0f}; // the A/D's full scale, linear peak
    float vuLevel = 0.0f;      // the VU integrator's state (linear RMS), audio thread only
    double vuTauSamples = 1.0; // VU time constant in samples (tau * sample rate), set in prepareToPlay
    bool lastTrueStereo = false; // engine 1's state is reset when the routing turns true stereo on or off
    EngineKey engineKey;
    // The reported latency is the engine's (ReverbEngine::latency: the same for every program and both modes): the dry
    // and bypassed signals are delayed by the same amount, in a ring per input channel.
    std::array<std::vector<float>, 2> dryDelay;
    size_t dryPos = 0;

    void applyMeterCal(int index); // the atomics only (audio-thread readable), not the UI state

    void processChunk(juce::AudioBuffer<float>&);
    void applyParameters();
    // Returns the block's input peak (after INPUT), for OVL.
    float processInputStage(const float* inL, const float* inR, const BlockBuffers&, int n, bool trueStereo);
    void processTransformerOutputs(const BlockBuffers&, int n, bool trueStereo);
    void mixOutput(juce::AudioBuffer<float>&, const BlockBuffers&, bool trueStereo, bool stereoOut);

    std::atomic<float>*pProgram, *pMix, *pPan, *pIn, *pOut, *pPower, *pTransformer, *pTrueStereo, *pSixteenBit, *pStereoWidth;
    std::array<std::atomic<float>*, rv::params::numJoys> pJoy;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Rev246Processor)
};
