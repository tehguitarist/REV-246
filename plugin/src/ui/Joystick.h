#pragma once

#include "ui/DesignComponent.h"
#include "ui/Layout.h"

#include <functional>

namespace rv::ui
{
// A remote lever: rest image, up or down while it is pushed. Pressing the top half steps the value up, the bottom
// half down; dragging up or down from where it was pressed pushes it that way instead, whichever half it started in,
// and back to the pressed half's direction when the pointer returns near the start. Holding repeats every half second; after a second it speeds up, and how fast scales with
// the lever's range (`range`): a 16-step lever stops at about its range per second, a 999-step one gets there in a few
// seconds. The mouse wheel is one step per notch.
class Joystick : public DesignComponent, private juce::Timer
{
  public:
    static constexpr int repeatMs = 500;           // between repeats while held
    static constexpr int accelerateAfterMs = 1000; // held this long, the repeat speeds up
    static constexpr int tickMs = 20;              // the accelerated repeat's timer

    Joystick(ScaledImages&, layout::Slot);
    std::function<void(int direction, int steps)> onStep; // direction +1 up / -1 down; steps to move (at least 1)
    std::function<int()> range;                           // the lever's largest step (its number of steps - 1)
    std::function<void()> onReset;                       // Option/Alt-click: back to the default step

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    void timerCallback() override;
    void fire(int steps = 1);
    double stepsPerSecond(double heldSeconds) const;
    void push(int direction); // a new direction: step at once and restart the repeat

    ScaledImages& images;
    int pushed = 0;  // -1 down, 0 rest, +1 up
    int pressedHalf = 0; // the direction of the half that was pressed, kept while the pointer stays near the press
    double pressMs = 0.0, lastTickMs = 0.0, owed = 0.0; // press time, last accelerated tick, steps not yet taken
    bool accelerating = false;
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
