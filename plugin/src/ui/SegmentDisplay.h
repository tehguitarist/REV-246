#pragma once

#include "ui/DesignComponent.h"

#include <functional>
#include <vector>

namespace rv::ui
{
// A red seven-segment display (the plugin draws it, so there is no font): ghost segments always faintly lit, lit
// ones with a glow, a point after a digit for each '.' in the text. The digits are split into fields, each holding
// one value (the reflections box has two, the decay box three); wheel and drag over a field nudge its value,
// and double-click types one in.
class SegmentDisplay : public DesignComponent
{
  public:
    struct Field
    {
        int first = 0, count = 3;
        bool leftAlign = false;
        std::function<void(int steps)> onNudge; // may be empty (the program display)
        std::function<void()> onReset;          // may be empty; Option/Alt-click resets to the default
        std::function<void()> onDoubleClick;    // may be empty; double-click to type a value
        juce::String text;
    };

    SegmentDisplay(juce::Rectangle<float> box, int numDigits, std::vector<Field> fieldList);

    void setText(size_t field, const juce::String&);
    void setNudge(size_t field, std::function<void(int steps)> f) { fields[field].onNudge = std::move(f); }
    void setReset(size_t field, std::function<void()> f) { fields[field].onReset = std::move(f); }
    void setDoubleClick(size_t field, std::function<void()> f) { fields[field].onDoubleClick = std::move(f); }
    void setFieldTooltip(size_t field, const juce::String&);

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseDrag(const juce::MouseEvent&) override;
    void mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails&) override;
    juce::String getTooltip() override;

    // The segments lit for a character, bit 0 = a (top), then b c d e f g clockwise from the top right.
    static unsigned segmentsFor(juce::juce_wchar);

  private:
    struct Digit
    {
        unsigned segments = 0;
        bool point = false;
    };
    std::vector<Digit> layoutDigits() const;
    int fieldAt(int x) const;

    int numDigits;
    std::vector<Field> fields;
    std::vector<juce::String> tips;
    CachedLayer layer;
    int dragLastY = 0, dragRemainder = 0;
    int dragField = -1; // the field being dragged; -1 when the press started no drag (Option/Alt-click, double-click)
};
} // namespace rv::ui
