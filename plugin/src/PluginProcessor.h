#pragma once

#include "dsp/ReverbEngine.h"
#include "dsp/Transformer.h"
#include "params/Parameters.h"
#include "presets/PresetManager.h"

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

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
    double getTailLengthSeconds() const override { return 20.0; } // L's longest RT
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

    juce::ValueTree uiState() { return apvts.state.getOrCreateChildWithName(UiProps::node, nullptr); }

  private:
    rv::dsp::ReverbEngine engine;
    rv::dsp::Transformer transformer;
    juce::LinearSmoothedValue<float> inGain, outGain, mixAmount, panLeft, panRight, wetOn, transformerOn;
    juce::AudioBuffer<float> scratch;
    std::atomic<float> inputVu{0.0f};
    std::atomic<bool> inputClip{false};
    float vu = 0.0f;         // the VU integrator's state (linear RMS), audio thread only
    double vuTauSamples = 0.0; // VU time constant in samples (tau * sample rate), set in prepareToPlay

    std::atomic<float>*pProgram, *pMix, *pPan, *pIn, *pOut, *pPower, *pTransformer;
    std::array<std::atomic<float>*, rv::params::numJoys> pJoy;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Rev246Processor)
};
