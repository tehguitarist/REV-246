#pragma once

#include "ui/DesignComponent.h"
#include "ui/Sprite.h"
#include "ui/ValueBubble.h"

#include <juce_audio_processors/juce_audio_processors.h>

#include <functional>
#include <optional>

namespace rv::ui
{
// What the faders and knobs share: one parameter, shown through its 0..1 value; the change gesture of a drag (opened
// on mouse down, closed on mouse up or on destruction, never opened twice); the value bubble; and the tooltip.
class ParameterControl : public DesignComponent
{
  public:
    juce::String helpText; // shown, with the current value, as the tooltip while the "?" is on
    juce::String getTooltip() override;

  protected:
    // `valueText` turns the parameter's 0..1 value into the text of the bubble and the tooltip.
    ParameterControl(juce::Rectangle<float> designBounds, juce::RangedAudioParameter&, ValueBubble&,
                     std::function<juce::String(float normalised)> valueText);
    ~ParameterControl() override;

    void beginGesture();
    void endGesture(); // does nothing when no gesture is open
    bool isGestureOpen() const { return gestureOpen; }

    // Sets the 0..1 value (clamped) and shows the bubble. Inside an open gesture it is one step of it, otherwise it is
    // a complete change (a wheel notch, a double-click).
    void setNormalised(float);
    // The 0..1 value after a wheel event, which moves it by `step` of the range per unit of wheel movement (`fineStep`
    // with Shift); empty when the event moves nothing.
    std::optional<float> valueAfterWheel(const juce::MouseEvent&, const juce::MouseWheelDetails&, float step,
                                         float fineStep) const;

    juce::RangedAudioParameter& param;

  private:
    ValueBubble& bubble;
    std::function<juce::String(float)> valueText;
    juce::ParameterAttachment attachment;
    bool gestureOpen = false;
};

// A vertical fader whose cap runs between two y positions of the design space. The parameter's 0..1 value is the
// cap position along the travel (0 = bottom), so a fader's printed scale and the host's value agree. Drag (Shift =
// fine), wheel, double-click or Option/Alt-click to reset. The cap is a sprite with its drop shadow baked.
class Fader : public ParameterControl
{
  public:
    Fader(ScaledImages&, layout::Image cap, float centreX, float topY, float bottomY, juce::RangedAudioParameter&,
          ValueBubble&, std::function<juce::String(float position)> valueText);

    // The user moved the cap (not the host): the editor uses this to move a linked fader with it.
    std::function<void(float oldPosition, float newPosition)> onUserMove;
    std::function<void()> onInteractionEnd; // after a drag, a wheel notch or a reset: the user is done for now

    float position() const { return param.getValue(); }
    // For the editor's linked fader. `insideGesture`: the caller has opened this parameter's gesture, so the change is
    // a step of it; otherwise it is a complete change.
    void setPosition(float p, bool insideGesture);

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseDoubleClick(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    float capCentreY() const; // design pixels
    juce::Rectangle<float> capDesignRect() const;
    float designY(const juce::MouseEvent&) const;
    float positionAtY(float y) const; // from a design-space y; unclamped
    void userSet(float p);
    void resetToDefault();

    ScaledImages& images;
    layout::Image capImage;
    float centreX, topY, bottomY;
    ShadowedSprite sprite;
    float grabOffset = 0.0f; // design pixels from the cap's centre to where it was grabbed
    float dragStartPos = 0.0f;
    int dragStartY = 0;
};

// What the pan knob and the width halo share: a vertical drag turns the value (Shift = fine), Option/Alt-click and,
// unless a subclass says otherwise, double-click set the default.
class DragKnob : public ParameterControl
{
  public:
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseDoubleClick(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  protected:
    // `valueText` takes the parameter's own value (not 0..1). `wheelStep` and `wheelFineStep` are fractions of the range.
    DragKnob(juce::Rectangle<float> designBounds, juce::RangedAudioParameter&, ValueBubble&,
             std::function<juce::String(float value)> valueText, float wheelStep, float wheelFineStep);

    void resetToDefault();

  private:
    float wheelStep, wheelFineStep;
    float dragStartNorm = 0.5f;
    int dragStartY = 0;
};

// The pan knob: the cap image rotates over 7 to 5 o'clock.
class PanKnob : public DragKnob
{
  public:
    PanKnob(ScaledImages&, juce::RangedAudioParameter&, ValueBubble&, std::function<juce::String(float value)> valueText);

    void paint(juce::Graphics&) override;

  private:
    ScaledImages& images;
    ShadowedSprite sprite;
};

// The stereo-width halo: an arc with the value in the middle and WIDTH underneath. Drag up to widen, down to narrow,
// wheel, Option/Alt-click for 100 (the default), or double-click, which calls onDoubleClick (the editor opens its
// numeric entry dialog). Ignores the mouse while disabled.
class WidthKnob : public DragKnob
{
  public:
    WidthKnob(juce::Rectangle<float> designBounds, juce::RangedAudioParameter&, ValueBubble&,
              std::function<juce::String(float value)> valueText);

    std::function<void()> onDoubleClick;

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDoubleClick(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    void enablementChanged() override { repaint(); }
};
} // namespace rv::ui
