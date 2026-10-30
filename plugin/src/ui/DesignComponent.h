#pragma once

#include "ui/Assets.h"
#include "ui/Design.h"

namespace rv::ui
{
// A component placed in the 3495x1530 design space (Layout.h). The editor lays children out directly in its own
// pixels (not through an AffineTransform on a 3495-wide content component), so every component sits on
// whole logical pixels and its pre-scaled images stay 1:1 with the device. Components convert
// design-space rectangles to local coordinates with toLocal().
class DesignComponent : public juce::Component, public juce::SettableTooltipClient
{
  public:
    explicit DesignComponent(juce::Rectangle<float> designBoundsIn) : designRect(designBoundsIn) {}

    juce::Rectangle<float> getDesignBounds() const { return designRect; }
    void setDesignBounds(juce::Rectangle<float> r) { designRect = r; }

    // s = editor width / design width. Edges are rounded (not position and size separately), so
    // neighbouring components can't drift apart by a pixel.
    void layoutForScale(float s)
    {
        editorScale = s;
        const auto r = designRect * s;
        setBounds(juce::Rectangle<int>::leftTopRightBottom(juce::roundToInt(r.getX()), juce::roundToInt(r.getY()),
                                                           juce::roundToInt(r.getRight()),
                                                           juce::roundToInt(r.getBottom())));
    }

    float getScale() const { return editorScale; }

    // A switched-off section's controls draw at design::dimmedAlpha but stay usable.
    // Repaints only this component, and only when the state changes.
    void setDimmed(bool shouldBeDimmed)
    {
        if (shouldBeDimmed == dimmed)
            return;
        dimmed = shouldBeDimmed;
        repaint();
    }
    bool isDimmed() const { return dimmed; }

    juce::Rectangle<float> toLocal(juce::Rectangle<float> r) const
    {
        return r * editorScale - getPosition().toFloat();
    }

    juce::Point<float> toLocal(juce::Point<float> designPoint) const
    {
        return designPoint * editorScale - getPosition().toFloat();
    }

  protected:
    // 1 normally, design::dimmedAlpha while dimmed.
    float contentAlpha() const { return dimmed ? design::dimmedAlpha : 1.0f; }

  private:
    juce::Rectangle<float> designRect; // where it sits in the design space
    float editorScale = 0.5f;          // editor width / design width
    bool dimmed = false;
};
} // namespace rv::ui
