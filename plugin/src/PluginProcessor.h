#pragma once

#include "dsp/CaptureEq.h"
#include "dsp/ReverbEngine.h"
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

    // The panel's LEVEL meter: the highest input sample (after the INPUT fader) since the last call, linear.
    float takeInputPeak() { return inputPeak.exchange(0.0f); }

    juce::ValueTree uiState() { return apvts.state.getOrCreateChildWithName(UiProps::node, nullptr); }

  private:
    rv::dsp::ReverbEngine engine;
    rv::dsp::CaptureEq captureEq;
    juce::LinearSmoothedValue<float> inGain, outGain, mixAmount, panLeft, panRight, wetOn, captureOn;
    juce::AudioBuffer<float> scratch;
    std::atomic<float> inputPeak{0.0f};

    std::atomic<float>*pProgram, *pMix, *pPan, *pIn, *pOut, *pPower, *pCaptureEq;
    std::array<std::atomic<float>*, rv::params::numJoys> pJoy;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Rev246Processor)
};
