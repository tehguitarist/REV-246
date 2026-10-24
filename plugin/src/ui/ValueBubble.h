#pragma once

#include "ui/Design.h"

#include <juce_gui_basics/juce_gui_basics.h>

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
        message = text;
        const auto font = juce::FontOptions(15.0f, juce::Font::bold);
        const auto w = juce::GlyphArrangement::getStringWidth(juce::Font(font), text) + 22.0f;
        auto* parent = getParentComponent();
        if (parent == nullptr)
            return;
        const auto a = parent->getLocalArea(&anchor, anchor.getLocalBounds());
        juce::Rectangle<int> r((int)std::ceil(w), 26);
        r.setCentre(a.getCentreX(), a.getY() - 18);
        setBounds(r.constrainedWithin(parent->getLocalBounds().reduced(4)));
        setVisible(true);
        repaint();
        startTimer(holdMs);
    }
    void hide() { stopTimer(); setVisible(false); }

    void paint(juce::Graphics& g) override
    {
        const auto b = getLocalBounds().toFloat().reduced(1.0f);
        g.setColour(juce::Colour(0xf2161615));
        g.fillRoundedRectangle(b, 5.0f);
        g.setColour(design::gold.withAlpha(0.8f));
        g.drawRoundedRectangle(b, 5.0f, 1.0f);
        g.setColour(design::labelSelected);
        g.setFont(juce::Font(juce::FontOptions(15.0f, juce::Font::bold)));
        g.drawText(message, getLocalBounds(), juce::Justification::centred, false);
    }

  private:
    void timerCallback() override { hide(); }
    juce::String message;
};
} // namespace rv::ui
