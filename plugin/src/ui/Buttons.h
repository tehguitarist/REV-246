#pragma once

#include "ui/DesignComponent.h"
#include "ui/Layout.h"

#include <juce_audio_processors/juce_audio_processors.h>

namespace rv::ui
{
// A programme button (S M L E N): a radio choice (the lamp beside it shows which is selected). The key is drawn raised
// at rest, whatever the selection, and pressed in only while the mouse button is held down on it.
class ProgramButton : public DesignComponent
{
  public:
    ProgramButton(ScaledImages&, layout::Slot, layout::Image out, layout::Image in);
    std::function<void()> onClick;
    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;

  private:
    ScaledImages& images;
    layout::Image outImage, inImage;
    bool held = false;
};

// An indicator lamp: one of two images.
class Led : public DesignComponent
{
  public:
    Led(ScaledImages&, layout::Slot, layout::Image on, layout::Image off);
    void setOn(bool);
    void paint(juce::Graphics&) override;

  private:
    ScaledImages& images;
    layout::Image onImage, offImage;
    bool on = false;
};

// The power button. It is pressed in while held (In Off if the plugin was off, In On if it was on) and takes its new
// place on release: Out In (lit) when on, Out Off when bypassed.
class PowerButton : public DesignComponent
{
  public:
    PowerButton(ScaledImages&);
    void setOn(bool);
    bool isOn() const { return on; }
    std::function<void(bool newState)> onToggle;
    juce::String getTooltip() override;
    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;

  private:
    ScaledImages& images;
    bool on = true, held = false;
};
} // namespace rv::ui
