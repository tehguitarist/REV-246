#pragma once

#include "PluginProcessor.h"
#include "ui/BottomBar.h"
#include "ui/Buttons.h"
#include "ui/Fader.h"
#include "ui/Joystick.h"
#include "ui/LevelMeter.h"
#include "ui/SegmentDisplay.h"
#include "ui/ValueBubble.h"

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

  private:
    using Joy = rv::params::Joy;

    class PanelLookAndFeel : public juce::LookAndFeel_V4
    {
      public:
        PanelLookAndFeel();
    };

    void parameterChanged(const juce::String&, float) override { needsRefresh = true; }
    void timerCallback() override;
    void refreshAll();
    void stepJoy(Joy, int direction, int hold);
    juce::RangedAudioParameter& param(const char* id) const;
    rv::params::Program program() const;
    void setUpJoy(int index);
    void moveLinked(rv::ui::Fader& from, rv::ui::Fader& to, float oldPos, float newPos);

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
    rv::ui::Fader mixFader, inputFader, outputFader;
    rv::ui::PresetBar presetBar;
    rv::ui::BarButton linkButton;
    rv::ui::BarButton transformerButton;
    rv::ui::HelpButton helpButton;
    rv::ui::BarLabel titleLabel, versionLabel;

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
    std::atomic<bool> needsRefresh{true};
    bool linkedGestureOpen = false;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Rev246Editor)
};
