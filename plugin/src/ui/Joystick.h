#pragma once

#include "ui/DesignComponent.h"
#include "ui/Layout.h"

#include <functional>

namespace rv::ui
{
// A remote lever: rest image, up or down while it is pushed. Pressing the top half steps the value up, the bottom
// half down; holding repeats every half second, and after a second of holding `hold` counts up so callers can
// accelerate (the unit's levers do). The mouse wheel is one step per notch.
class Joystick : public DesignComponent, private juce::Timer
{
  public:
    static constexpr int repeatMs = 500;           // between repeats while held
    static constexpr int accelerateAfterMs = 1000; // held this long, onStep's `hold` starts counting

    Joystick(ScaledImages&, layout::Slot);
    std::function<void(int direction, int hold)> onStep; // direction +1 up / -1 down; hold = repeats so far
    std::function<void()> onReset;                       // Option/Alt-click: back to the default step

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    void timerCallback() override;
    void fire();

    ScaledImages& images;
    int pushed = 0;  // -1 down, 0 rest, +1 up
    int repeats = 0; // timer ticks since the press
};

// The wedge above a lever: a bar of segments (faint dividing lines always visible) with the one nearest the value lit.
// Click or drag on it to set the value at the pointer, wheel to step, Option/Alt-click for the default.
class Funnel : public DesignComponent
{
  public:
    explicit Funnel(int column);
    void setFraction(float f); // 0 = minimum (narrow end) .. 1 = maximum
    std::function<void(float fraction)> onSetFraction; // click or drag to a position (0..1)
    std::function<void(int direction)> onStep; // wheel notch
    std::function<void()> onReset;             // Option/Alt-click: back to the default step
    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    float fractionAt(const juce::MouseEvent&) const;
    CachedLayer layer;
    int column;
    int litSegment = 4; // 0 = the bottom segment
    bool dragging = false;
};
} // namespace rv::ui
