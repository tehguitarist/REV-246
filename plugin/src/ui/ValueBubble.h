#pragma once

#include "ui/Design.h"

#include <juce_gui_basics/juce_gui_basics.h>

#include <cmath>

namespace rv::ui
{
// The value readout shown beside a control while it is dragged (two decimals), whether or not the "?" tooltips are on.
// One per editor, added over the panel; controls call show() on every change and it hides itself shortly after.
class ValueBubble : public juce::Component, private juce::Timer
{
  public:
    ValueBubble()
    {
        setInterceptsMouseClicks(false, false);
        setAlwaysOnTop(true);
        setVisible(false);
    }

    // Shows `text` above `anchor` (clamped inside the parent), and keeps it up for `holdMs` after the last call.
    void show(const juce::Component& anchor, const juce::String& text, int holdMs = 900)
    {
        auto* parent = getParentComponent();
        if (parent == nullptr)
            return;
        message = text;
        const auto w = juce::GlyphArrangement::getStringWidth(font, text) + 2.0f * padX;
        const auto a = parent->getLocalArea(&anchor, anchor.getLocalBounds());
        juce::Rectangle<int> r((int)std::ceil(w), height);
        r.setCentre(a.getCentreX(), a.getY() - aboveAnchor);
        setBounds(r.constrainedWithin(parent->getLocalBounds().reduced(edgeMargin)));
        setVisible(true);
        repaint();
        startTimer(holdMs);
    }
    void hide() { stopTimer(); setVisible(false); }

    void paint(juce::Graphics& g) override
    {
        const auto b = getLocalBounds().toFloat().reduced(1.0f);
        g.setColour(juce::Colour(0xf2161615)); // near-black, 95% opaque
        g.fillRoundedRectangle(b, cornerRadius);
        g.setColour(design::gold.withAlpha(0.8f));
        g.drawRoundedRectangle(b, cornerRadius, 1.0f);
        g.setColour(design::labelSelected);
        g.setFont(font);
        g.drawText(message, getLocalBounds(), juce::Justification::centred, false);
    }

  private:
    // In logical pixels: the box's height, the text's side padding, the gap from the anchor's top edge to the box's
    // centre, the margin kept from the parent's edges, and the corner radius.
    static constexpr int height = 26, padX = 11, aboveAnchor = 18, edgeMargin = 4;
    static constexpr float cornerRadius = 5.0f;

    void timerCallback() override { hide(); }

    const juce::Font font{juce::FontOptions(15.0f, juce::Font::bold)};
    juce::String message;
};
} // namespace rv::ui
