#pragma once

#include "presets/PresetManager.h"
#include "ui/DesignComponent.h"

#include <functional>
#include <memory>
#include <vector>

namespace rv::ui
{
// A control of the bottom bar: gold-outlined text, lit (white with a glow) while `on` for the toggles.
class BarButton : public DesignComponent
{
  public:
    enum class Kind { momentary, toggle };

    BarButton(juce::Rectangle<float> slot, juce::String label, Kind = Kind::momentary, float capHeight = 34.0f);

    void setOn(bool);
    bool isOn() const { return on; }
    void setLabel(const juce::String&);
    juce::String getLabel() const { return label; }
    void setItalic(bool);
    bool isItalic() const { return italic; }
    void setAlwaysLit(bool);
    bool isAlwaysLit() const { return alwaysLit; }
    // A toggle that draws as disabled (greyed out) while it is off, but still takes clicks.
    void setGreyWhenOff(bool shouldBeGrey) { greyWhenOff = shouldBeGrey; repaint(); }
    std::function<void()> onClick;
    std::function<void(bool)> onToggle; // toggles: the requested state

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    std::unique_ptr<juce::AccessibilityHandler> createAccessibilityHandler() override;

  private:
    void activate(); // the button's action: a click for a momentary one, a requested change of state for a toggle
    void enablementChanged() override { repaint(); }

    juce::String label;
    Kind kind;
    float textCapHeight;
    bool on = false, held = false, italic = false, alwaysLit = false, greyWhenOff = false;
};

// The linked / unlinked pair in the bar: the artwork itself is the button, no gold frame.
class LinkButton : public DesignComponent
{
  public:
    LinkButton(ScaledImages&, juce::Rectangle<float> slot);

    void setOn(bool);
    bool isOn() const { return on; }
    std::function<void(bool)> onToggle;

    void paint(juce::Graphics&) override;
    void mouseDown(const juce::MouseEvent&) override;
    void mouseUp(const juce::MouseEvent&) override;
    std::unique_ptr<juce::AccessibilityHandler> createAccessibilityHandler() override;

  private:
    void activate();
    // The icon is small, so the clickable area extends this far (design pixels) around it.
    static constexpr float hitPad = 24.0f;
    ScaledImages& images;
    bool on = true, held = false;
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
    void saveAs(const juce::String& chosen); // asks again as needed (empty or factory name, replacing a user preset)

    rv::PresetManager& presets;
    juce::Component& parent;
    juce::Component::SafePointer<juce::Component> parentAlive; // dialogs and menus can outlive the editor
    BarButton prev, next, name, save;
};
} // namespace rv::ui
