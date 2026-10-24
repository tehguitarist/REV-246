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
    std::function<void()> onReset;                       // Cmd/Alt-click: back to the default step

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

// The wedge above a lever: a bar of segments (faint dividing lines always visible) with the one nearest the value lit.
// Click or drag on it to set the value at the pointer, wheel to step, Cmd/Alt-click for the default.
class Funnel : public DesignComponent
{
  public:
    explicit Funnel(int column);
    void setFraction(float f); // 0 = minimum (narrow end) .. 1 = maximum
    int litSegment() const { return segment; }
    std::function<void(float fraction)> onUserSet; // click or drag to a position (0..1)
    std::function<void(int direction)> onUserStep; // wheel notch
    std::function<void()> onUserReset;             // Cmd/Alt-click: back to the default step
    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    float fractionAt(const juce::MouseEvent&) const;
    CachedLayer layer;
    int column;
    int segment = 4; // 0 = the bottom segment
    bool dragging = false;
};
} // namespace rv::ui
