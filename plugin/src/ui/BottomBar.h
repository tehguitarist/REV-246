#pragma once

#include "presets/PresetManager.h"
#include "ui/DesignComponent.h"

namespace rv::ui
{
// A control of the bottom bar: gold-outlined text, lit (white with a glow) while `on` for the toggles.
class BarButton : public DesignComponent
{
  public:
    enum class Kind { momentary, toggle };

    BarButton(juce::Rectangle<float> slot, juce::String label, Kind = Kind::momentary);

    void setOn(bool);
    bool isOn() const { return on; }
    void setLabel(const juce::String&);
    juce::String getLabel() const { return label; }
    void setItalic(bool);
    bool isItalic() const { return italic; }
    void setAlwaysLit(bool);
    bool isAlwaysLit() const { return alwaysLit; }
    std::function<void()> onClick;
    std::function<void(bool)> onToggle; // toggles: the requested state

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    std::unique_ptr<juce::AccessibilityHandler> createAccessibilityHandler() override;

  private:
    juce::String label;
    Kind kind;
    bool on = false, held = false, italic = false, alwaysLit = false;
};

// Static text in the bottom bar (the product name, the version).
class BarLabel : public DesignComponent
{
  public:
    BarLabel(juce::Rectangle<float> slot, juce::String text, float capHeight, juce::Colour colour,
             juce::Justification justification = juce::Justification::centred);
    void paint(juce::Graphics&) override;

  private:
    juce::String text;
    float capHeight;
    juce::Colour colour;
    juce::Justification justification;
};

// The "?" that turns the tooltips on and off (off by default). Its own tooltip is always available, so the way to
// turn the others on can be found with them off.
class HelpButton : public BarButton
{
  public:
    explicit HelpButton(juce::Rectangle<float> slot) : BarButton(slot, "?", Kind::toggle) {}
    juce::String getTooltip() override { return isOn() ? "Disable tooltips" : "Enable tooltips"; }
};

// Preset name with previous / next (a click on the name opens the preset menu) and SAVE (asks for a name, then asks
// before replacing an existing user preset).
class PresetBar
{
  public:
    PresetBar(rv::PresetManager&, juce::Component& parent);
    void layoutForScale(float s);
    void refresh(); // the name, italic with a "*" while modified
    void addTo(std::vector<DesignComponent*>&);
    const BarButton& getNameButton() const { return name; }

  private:
    void showMenu();
    void promptSaveAs(juce::String initial = {});

    rv::PresetManager& presets;
    juce::Component& parent;
    BarButton prev, next, name, save;
};
} // namespace rv::ui
