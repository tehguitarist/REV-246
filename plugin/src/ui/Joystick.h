#pragma once

#include "ui/DesignComponent.h"
#include "ui/Layout.h"

namespace rv::ui
{
// A remote lever: rest image, up or down while it is pushed. Pressing the top half steps the value up, the bottom
// half down; holding repeats every half second, and after a second of holding `hold` counts up so callers can
// accelerate (the unit's levers do). The mouse wheel is one step per notch.
class Joystick : public DesignComponent, private juce::Timer
{
  public:
    static constexpr int repeatMs = 500;

    Joystick(ScaledImages&, layout::Slot);
    std::function<void(int direction, int hold)> onStep; // direction +1 up / -1 down; hold = repeats so far

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    void timerCallback() override;
    void fire();

    ScaledImages& images;
    int pushed = 0; // -1 down, 0 rest, +1 up
    int repeats = 0;
};

// The wedge above a lever: a glowing red line at the value's height, following the wedge's width.
class Funnel : public DesignComponent
{
  public:
    explicit Funnel(int column);
    void setFraction(float f); // 0 = minimum (narrow end) .. 1 = maximum
    void paint(juce::Graphics&) override;

  private:
    CachedLayer layer;
    int column;
    float fraction = 0.5f;
};
} // namespace rv::ui
