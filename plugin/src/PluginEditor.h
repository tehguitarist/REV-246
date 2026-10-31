#pragma once

#include "PluginProcessor.h"
#include "ui/BottomBar.h"
#include "ui/Buttons.h"
#include "ui/Fader.h"
#include "ui/Joystick.h"
#include "ui/LevelMeter.h"
#include "ui/SegmentDisplay.h"
#include "ui/ValueBubble.h"

#include <array>
#include <atomic>
#include <functional>
#include <memory>
#include <vector>

// The panel: baked artwork (ui/main.png) plus image controls laid out from ui-layout.csv, and the bottom bar.
// Resizable with a locked aspect ratio; the width is kept in the session state. Nothing repaints while idle: each
// control repaints only when its value changes, and the meter only when its lit LEDs change.
class Rev246Editor : public juce::AudioProcessorEditor, private juce::Timer, private juce::AudioProcessorValueTreeState::Listener
{
  public:
    explicit Rev246Editor(Rev246Processor&);
    ~Rev246Editor() override;

    void paint(juce::Graphics&) override;
    void resized() override;

    // Test hooks.
    rv::ui::HelpButton& getHelpButton() { return helpButton; }
    bool tooltipAllowed(const juce::Component& c) const { return tooltipWindow.allows(c); }
    void refreshNow() { refreshAll(); }
    rv::ui::LevelMeter& getLevelMeter() { return levelMeter; }
    rv::ui::BarButton& getStereoButton() { return stereoButton; }
    rv::ui::BarButton& getFormatButton() { return formatButton; }
    rv::ui::WidthKnob& getWidthKnob() { return widthKnob; }
    rv::ui::WidthKnob& getPreDelayKnob() { return preDelayKnob; }
    rv::ui::BarButton& getDuckButton() { return duckButton; }
    rv::ui::PresetBar& getPresetBar() { return presetBar; }

  private:
    using Joy = rv::params::Joy;
    using Program = rv::params::Program;

    class PanelLookAndFeel : public juce::LookAndFeel_V4
    {
      public:
        PanelLookAndFeel();
    };

    void parameterChanged(const juce::String&, float) override { needsRefresh = true; } // any thread: only flags
    void timerCallback() override;
    void refreshAll();
    void refreshStateTips(bool useTrueStereo);
    void refreshJoys(Program);

    // Construction steps, in the order the constructor runs them.
    void layoutBottomBar();
    void setUpProgramControls();
    void setUpJoy(int index);
    void setUpDisplays();
    void setUpBarButtons();
    void setUpKnobsAndFaders();
    void addChildren();
    void restoreSize();
    void addParameterListeners();
    void removeParameterListeners();

    // Parameter access (message thread).
    juce::RangedAudioParameter& param(const char* id) const;
    bool isOn(const char* id) const { return param(id).getValue() >= 0.5f; }
    int joyStep(Joy) const;
    Program program() const;
    void setParamValue(juce::RangedAudioParameter&, float plainValue); // one complete host gesture
    void setToggle(const char* id, bool on) { setParamValue(param(id), on ? 1.0f : 0.0f); }
    void selectProgram(int index);

    void stepJoy(Joy, int direction, int hold);
    void setJoyStep(Joy, int step);
    void setJoy(Joy, float fraction); // 0..1 -> the nearest step
    void resetJoy(Joy);               // the default step
    void moveLinked(rv::ui::Fader& to, float oldPos, float newPos);
    void endLinkedGesture();
    juce::String calLabel() const;

    // Typed entry (double-click a display field). `onApply` returns false for text it rejects, which shows a warning
    // and then calls `onInvalid` with the same text so it can be edited.
    void promptJoyValue(Joy, const juce::String& initial = {});
    void promptWidthValue(const juce::String& initial = {});
    void promptPreDelayValue(const juce::String& initial = {});
    void showDuckPanel();
    void showValueDialog(const juce::String& title, const juce::String& message, const juce::String& shown,
                         const juce::String& hint, std::function<bool(const juce::String&)> onApply,
                         std::function<void(const juce::String&)> onInvalid);

    Rev246Processor& proc;
    PanelLookAndFeel lookAndFeel;
    rv::ui::ScaledImages images;
    rv::ui::ValueBubble bubble;

    std::vector<std::unique_ptr<rv::ui::ProgramButton>> programButtons;
    std::vector<std::unique_ptr<rv::ui::Led>> leds;
    rv::ui::PowerButton powerButton;
    std::vector<std::unique_ptr<rv::ui::Joystick>> joysticks;
    std::vector<std::unique_ptr<rv::ui::Funnel>> funnels;
    rv::ui::SegmentDisplay programDisplay, lowpassDisplay, reflectionsDisplay, decayDisplay;
    rv::ui::LevelMeter levelMeter;
    rv::ui::PanKnob panKnob;
    rv::ui::WidthKnob widthKnob, preDelayKnob;
    rv::ui::Fader mixFader, inputFader, outputFader;
    rv::ui::PresetBar presetBar;
    rv::ui::LinkButton linkButton;
    rv::ui::BarButton duckButton;
    rv::ui::BarButton transformerButton;
    rv::ui::BarButton formatButton;
    rv::ui::BarButton stereoButton;
    rv::ui::BarButton calButton;
    rv::ui::HelpButton helpButton;
    rv::ui::BarLabel titleLabel, versionLabel;
    int meterCalIndex = rv::design::defaultMeterCal; // follows the processor's, which a restored session can change

    // Tips show only while the "?" is on, except the "?"'s own.
    static constexpr int tipDelayMs = 300;
    class TipWindow : public juce::TooltipWindow
    {
      public:
        TipWindow(juce::Component* parent, const rv::ui::HelpButton& help) : juce::TooltipWindow(parent, tipDelayMs), helpButton(help) {}
        bool allows(const juce::Component& c) const { return helpButton.isOn() || &c == &helpButton; }

      private:
        juce::String getTipFor(juce::Component& c) override { return allows(c) ? juce::TooltipWindow::getTipFor(c) : juce::String(); }
        const rv::ui::HelpButton& helpButton;
    };
    TipWindow tooltipWindow{this, helpButton};

    std::vector<rv::ui::DesignComponent*> designComponents;
    juce::Component::SafePointer<juce::AlertWindow> valueDialog; // the open typed-entry dialog, closed with the editor
    std::atomic<bool> needsRefresh{true};
    bool linkedGestureOpen = false;
    int shownTipKey = -1;                        // what refreshStateTips() last built
    int shownProgram = -1;                       // what refreshJoys() last drew: it only redraws on a change
    std::array<int, rv::params::numJoys> shownSteps{};
    bool stereoAvailable = false; // a stereo in / stereo out bus layout (the stereo button's enable state)
    bool widthAvailable = true;   // a stereo output (the width halo's enable state)

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Rev246Editor)
};
