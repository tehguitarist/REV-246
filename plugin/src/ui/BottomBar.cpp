#include "ui/BottomBar.h"

#include "ui/Design.h"

#include "params/Parameters.h"

#include <array>
#include <utility>

namespace rv::ui
{
namespace
{
// The bar buttons' look in design pixels: corner radius, outline width, and text inset.
constexpr float barCornerRadius = 14.0f, barOutline = 4.0f, barTextInset = 12.0f;
constexpr float disabledAlpha = 0.3f;
} // namespace

BarButton::BarButton(juce::Rectangle<float> slot, juce::String text, Kind k, float capHeight)
    : DesignComponent(slot), label(std::move(text)), kind(k), textCapHeight(capHeight)
{
    setMouseCursor(juce::MouseCursor::PointingHandCursor);
    setTitle(label);
}

void BarButton::setOn(bool s)
{
    if (s == on)
        return;
    on = s;
    repaint();
}

void BarButton::setLabel(const juce::String& l)
{
    if (l == label)
        return;
    label = l;
    setTitle(label);
    repaint();
}

void BarButton::setItalic(bool shouldBeItalic)
{
    if (shouldBeItalic == italic)
        return;
    italic = shouldBeItalic;
    repaint();
}

void BarButton::setAlwaysLit(bool shouldAlwaysBeLit)
{
    if (shouldAlwaysBeLit == alwaysLit)
        return;
    alwaysLit = shouldAlwaysBeLit;
    repaint();
}

void BarButton::paint(juce::Graphics& g)
{
    const auto s = getScale();
    const auto b = getLocalBounds().toFloat().reduced(2.0f);
    const auto lit = isEnabled() && (alwaysLit || (kind == Kind::toggle ? on : held));
    // Greyed out: disabled, or (for the DUCK toggle) switched off.
    const auto fade = !isEnabled() || (greyWhenOff && kind == Kind::toggle && !on) ? disabledAlpha : 1.0f;
    g.setColour(design::gold.withAlpha((lit ? 0.22f : 0.08f) * fade));
    g.fillRoundedRectangle(b, barCornerRadius * s);
    g.setColour(design::gold.withAlpha((lit ? 0.95f : 0.6f) * fade));
    g.drawRoundedRectangle(b, barCornerRadius * s, juce::jmax(1.0f, barOutline * s));
    const auto font = italic ? labelFont(textCapHeight * s).italicised() : labelFont(textCapHeight * s);
    g.setFont(font);
    const auto labelColour = !alwaysLit && kind == Kind::toggle && !on ? design::labelDimmed : design::labelSelected;
    g.setColour(labelColour.withMultipliedAlpha(fade));
    g.drawFittedText(label, getLocalBounds().reduced(juce::roundToInt(barTextInset * s), 0), juce::Justification::centred, 1, 1.0f);
}

void BarButton::mouseDown(const juce::MouseEvent&)
{
    if (!isEnabled())
        return;
    held = true;
    repaint();
}

void BarButton::mouseUp(const juce::MouseEvent& e)
{
    const auto wasHeld = std::exchange(held, false); // also when it was disabled while held
    repaint();
    if (wasHeld && isEnabled() && getLocalBounds().contains(e.getPosition()))
        activate();
}

void BarButton::activate()
{
    if (!isEnabled())
        return;
    if (kind == Kind::momentary)
    {
        if (onClick)
            onClick();
    }
    else if (onToggle)
        onToggle(!on);
    else
        setOn(!on);
}

std::unique_ptr<juce::AccessibilityHandler> BarButton::createAccessibilityHandler()
{
    return std::make_unique<juce::AccessibilityHandler>(
        *this, kind == Kind::toggle ? juce::AccessibilityRole::toggleButton : juce::AccessibilityRole::button,
        juce::AccessibilityActions().addAction(juce::AccessibilityActionType::press, [this] { activate(); }));
}

//==============================================================================
LinkButton::LinkButton(ScaledImages& i, juce::Rectangle<float> slot)
    : DesignComponent(slot.expanded(hitPad)), images(i)
{
    setMouseCursor(juce::MouseCursor::PointingHandCursor);
    setTitle("Link input/output");
}

void LinkButton::setOn(bool state)
{
    if (state == on)
        return;
    on = state;
    repaint();
}

void LinkButton::paint(juce::Graphics& g)
{
    if (held)
    {
        g.setColour(design::gold.withAlpha(0.18f));
        g.fillRoundedRectangle(getLocalBounds().toFloat(), 6.0f * getScale());
    }
    drawFitted(g, images, on ? layout::Image::linked : layout::Image::unlinked,
               getLocalBounds().toFloat().reduced(hitPad * getScale()));
}

void LinkButton::mouseDown(const juce::MouseEvent&)
{
    held = true;
    repaint();
}

void LinkButton::mouseUp(const juce::MouseEvent& e)
{
    const auto wasHeld = std::exchange(held, false);
    repaint();
    if (wasHeld && getLocalBounds().contains(e.getPosition()))
        activate();
}

void LinkButton::activate()
{
    if (onToggle)
        onToggle(!on);
    else
        setOn(!on);
}

std::unique_ptr<juce::AccessibilityHandler> LinkButton::createAccessibilityHandler()
{
    return std::make_unique<juce::AccessibilityHandler>(
        *this, juce::AccessibilityRole::toggleButton,
        juce::AccessibilityActions().addAction(juce::AccessibilityActionType::press, [this] { activate(); }));
}

BarLabel::BarLabel(juce::Rectangle<float> slot, juce::String t, float cap, juce::Colour c, juce::Justification j)
    : DesignComponent(slot), text(std::move(t)), capHeight(cap), colour(c), justification(j)
{
    setInterceptsMouseClicks(false, false);
}

void BarLabel::paint(juce::Graphics& g)
{
    g.setColour(colour);
    g.setFont(labelFont(capHeight * getScale()));
    g.drawFittedText(text, getLocalBounds(), justification, 1, 1.0f);
}

//==============================================================================
namespace
{
constexpr float barY = design::barTop + 35.0f, barH = 80.0f;
}

PresetBar::PresetBar(rv::PresetManager& p, juce::Component& parentIn)
    : presets(p), parent(parentIn), parentAlive(&parentIn), prev({110, barY, 90, barH}, juce::CharPointer_UTF8("\xe2\x97\x80")),
      next({916, barY, 90, barH}, juce::CharPointer_UTF8("\xe2\x96\xb6")), name({208, barY, 700, barH}, ""),
      save({1016, barY, 200, barH}, "SAVE")
{
    prev.setTooltip("Previous preset");
    next.setTooltip("Next preset");
    name.setTooltip("Presets: load, save, delete. The name turns italic with a * once a control differs from the saved preset.");
    save.setTooltip("Save the current settings. Factory names cannot be overwritten; a user preset with the chosen name asks first.");
    prev.onClick = [this] { presets.step(-1); };
    next.onClick = [this] { presets.step(1); };
    name.onClick = [this] { showMenu(); };
    save.onClick = [this] { promptSaveAs(); };
    for (auto* c : {&prev, &next, &name, &save})
        parent.addAndMakeVisible(*c);
    refresh();
}

void PresetBar::layoutForScale(float s)
{
    for (auto* c : {&prev, &next, &name, &save})
        c->layoutForScale(s);
}

void PresetBar::addTo(std::vector<DesignComponent*>& v)
{
    for (auto* c : {&prev, &next, &name, &save})
        v.push_back(c);
}

void PresetBar::refresh()
{
    const auto modified = presets.isModified();
    name.setItalic(modified);
    name.setLabel(presets.currentName() + (modified ? " *" : ""));
}

void PresetBar::showMenu()
{
    juce::PopupMenu menu;
    int id = 1;
    std::vector<juce::String> byId{{}};
    menu.addSectionHeader("Factory");
    // The five Init presets at the top level; the Altiverb ones in one submenu per program.
    const auto factory = presets.factoryNames();
    std::array<juce::PopupMenu, rv::params::numPrograms> altiverb;
    for (int i = 0; i < factory.size(); ++i)
    {
        const auto& n = factory[i];
        auto& m = i < rv::params::numPrograms ? menu : altiverb[(size_t)rv::PresetManager::factoryProgram(i)];
        m.addItem(id++, n, true, n == presets.currentName());
        byId.push_back(n);
    }
    for (int p = 0; p < rv::params::numPrograms; ++p)
        if (altiverb[(size_t)p].getNumItems() > 0)
            menu.addSubMenu("Altiverb " + rv::params::programLongName((rv::params::Program)p), altiverb[(size_t)p]);
    const auto users = presets.userNames();
    if (!users.isEmpty())
    {
        menu.addSectionHeader("User");
        for (const auto& n : users)
        {
            menu.addItem(id++, n, true, n == presets.currentName());
            byId.push_back(n);
        }
    }
    menu.addSeparator();
    const int saveAsId = id++, del = id++, reveal = id++;
    menu.addItem(saveAsId, "Save as...");
    menu.addItem(del, "Delete current", !presets.isFactory(presets.currentName()));
    menu.addItem(reveal, "Show presets folder");

    menu.showMenuAsync(juce::PopupMenu::Options().withTargetComponent(&name), [this, alive = parentAlive, byId, saveAsId, del, reveal](int r) {
        if (r <= 0 || alive == nullptr)
            return;
        if (r == saveAsId)
            promptSaveAs();
        else if (r == del)
            presets.remove(presets.currentName());
        else if (r == reveal)
        {
            presets.folder().createDirectory();
            presets.folder().revealToUser();
        }
        else
            presets.load(byId[(size_t)r]);
        refresh();
    });
}

void PresetBar::promptSaveAs(juce::String initial)
{
    if (initial.isEmpty())
        initial = presets.currentName();
    auto* w = new juce::AlertWindow("Save preset", "Name for the preset:", juce::MessageBoxIconType::NoIcon, &parent);
    w->addTextEditor("name", initial);
    w->addButton("Save As", 1, juce::KeyPress(juce::KeyPress::returnKey));
    w->addButton("Cancel", 0, juce::KeyPress(juce::KeyPress::escapeKey));
    w->enterModalState(true, juce::ModalCallbackFunction::create([this, alive = parentAlive, w](int r) {
                           if (r == 1 && alive != nullptr)
                               saveAs(rv::PresetManager::sanitise(w->getTextEditorContents("name")));
                       }),
                       true);
}

void PresetBar::saveAs(const juce::String& chosen)
{
    if (chosen.isEmpty())
    {
        promptSaveAs();
        return;
    }
    if (presets.isFactory(chosen))
    {
        // Factory names are read-only; return to the name field rather than failing silently.
        juce::AlertWindow::showAsync(juce::MessageBoxOptions()
                                         .withIconType(juce::MessageBoxIconType::InfoIcon)
                                         .withTitle("Factory preset")
                                         .withMessage("\"" + chosen + "\" is a factory preset name. Choose another name.")
                                         .withButton("OK")
                                         .withAssociatedComponent(&parent),
                                     [this, alive = parentAlive, chosen](int) {
                                         if (alive != nullptr)
                                             promptSaveAs(chosen);
                                     });
        return;
    }
    if (presets.userNames().contains(chosen))
    {
        // Only user presets can be replaced.
        juce::AlertWindow::showAsync(
            juce::MessageBoxOptions()
                .withIconType(juce::MessageBoxIconType::QuestionIcon)
                .withTitle("Override preset?")
                .withMessage("A user preset called \"" + chosen + "\" already exists. Do you want to override it?")
                .withButton("Yes")
                .withButton("No")
                .withAssociatedComponent(&parent),
            [this, alive = parentAlive, chosen](int choice) {
                if (alive == nullptr)
                    return;
                if (choice == 1)
                {
                    presets.save(chosen);
                    refresh();
                }
                else
                    promptSaveAs(chosen); // No: go back to changing the name.
            });
        return;
    }
    presets.save(chosen);
    refresh();
}
} // namespace rv::ui
