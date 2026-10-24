#pragma once

#include "ui/DesignComponent.h"
#include "ui/Sprite.h"
#include "ui/ValueBubble.h"

#include <juce_audio_processors/juce_audio_processors.h>

namespace rv::ui
{
// A vertical fader whose cap runs between two y positions of the design space. The parameter's 0..1 value is the
// cap position along the travel (0 = bottom), so a fader's printed scale and the host's value agree. Drag (Shift =
// fine), wheel, double-click or Cmd/Alt-click to reset. The cap is a sprite with its drop shadow baked.
class Fader : public DesignComponent
{
  public:
    Fader(ScaledImages&, layout::Image cap, float centreX, float topY, float bottomY, juce::RangedAudioParameter&,
          ValueBubble&, std::function<juce::String(float position)> valueText);

    juce::String helpText; // shown, with the current value, as the tooltip while the "?" is on
    juce::String getTooltip() override;

    // The user moved the cap (not the host): the editor uses this to move a linked fader with it.
    std::function<void(float oldPosition, float newPosition)> onUserMove;
    std::function<void()> onDragEnd;

    float position() const { return param.getValue(); }
    void setPosition(float p, bool asGesture); // asGesture: inside an open drag, otherwise a complete change

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseDoubleClick(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    float capCentreY() const; // design pixels
    juce::Rectangle<float> capDesignRect() const;
    void userSet(float p, bool gesture);

    ScaledImages& images;
    layout::Image capImage;
    float centreX, topY, bottomY;
    juce::RangedAudioParameter& param;
    ValueBubble& bubble;
    std::function<juce::String(float)> valueText;
    juce::ParameterAttachment attachment;
    ShadowedSprite sprite;
    float grabOffset = 0.0f, dragStartPos = 0.0f;
    int dragStartY = 0;
    bool dragging = false;
};

// The pan knob: the cap image rotates over 7 to 5 o'clock. Vertical drag (Shift = fine), wheel, double-click to centre.
class PanKnob : public DesignComponent
{
  public:
    PanKnob(ScaledImages&, juce::RangedAudioParameter&, ValueBubble&, std::function<juce::String(float value)> valueText);

    juce::String helpText;
    juce::String getTooltip() override;

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    void mouseDoubleClick(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;

  private:
    void setNormalised(float n, bool gesture);

    ScaledImages& images;
    juce::RangedAudioParameter& param;
    ValueBubble& bubble;
    std::function<juce::String(float)> valueText;
    juce::ParameterAttachment attachment;
    ShadowedSprite sprite;
    float dragStartNorm = 0.5f;
    int dragStartY = 0;
    bool dragging = false;
};
} // namespace rv::ui
