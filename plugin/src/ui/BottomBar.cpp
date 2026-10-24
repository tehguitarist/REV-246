#include "ui/BottomBar.h"

#include "ui/Design.h"

namespace rv::ui
{
BarButton::BarButton(juce::Rectangle<float> slot, juce::String text, Kind k)
    : DesignComponent(slot), label(std::move(text)), kind(k)
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
    repaint();
}

void BarButton::paint(juce::Graphics& g)
{
    const auto s = getScale();
    const auto b = getLocalBounds().toFloat().reduced(2.0f);
    const auto lit = kind == Kind::toggle ? on : held;
    g.setColour(design::gold.withAlpha(lit ? 0.22f : 0.08f));
    g.fillRoundedRectangle(b, 14.0f * s);
    g.setColour(design::gold.withAlpha(lit ? 0.95f : 0.6f));
    g.drawRoundedRectangle(b, 14.0f * s, juce::jmax(1.0f, 4.0f * s));
    const auto font = labelFont(34.0f * s);
    g.setFont(font);
    g.setColour(kind == Kind::toggle && !on ? design::labelDimmed : design::labelSelected);
    g.drawFittedText(label, getLocalBounds().reduced(juce::roundToInt(12.0f * s), 0), juce::Justification::centred, 1, 1.0f);
}

void BarButton::mouseDown(const juce::MouseEvent&)
{
    held = true;
    repaint();
}

void BarButton::mouseUp(const juce::MouseEvent& e)
{
    held = false;
    repaint();
    if (!getLocalBounds().contains(e.getPosition()))
        return;
    if (kind == Kind::toggle)
    {
        if (onToggle)
            onToggle(!on);
        else
            setOn(!on);
    }
    else if (onClick)
        onClick();
}

std::unique_ptr<juce::AccessibilityHandler> BarButton::createAccessibilityHandler()
{
    return std::make_unique<juce::AccessibilityHandler>(
        *this, kind == Kind::toggle ? juce::AccessibilityRole::toggleButton : juce::AccessibilityRole::button,
        juce::AccessibilityActions().addAction(juce::AccessibilityActionType::press, [this] {
            if (kind == Kind::toggle)
            {
                if (onToggle)
                    onToggle(!on);
            }
            else if (onClick)
                onClick();
        }));
}

//==============================================================================
namespace
{
constexpr float barY = design::barTop + 35.0f, barH = 80.0f;
}

PresetBar::PresetBar(rv::PresetManager& p, juce::Component& parentIn)
    : presets(p), parent(parentIn), prev({110, barY, 90, barH}, juce::CharPointer_UTF8("\xe2\x97\x80")),
      next({1000, barY, 90, barH}, juce::CharPointer_UTF8("\xe2\x96\xb6")), name({208, barY, 784, barH}, ""),
      save({1100, barY, 220, barH}, "SAVE")
{
    prev.setTooltip("Previous preset");
    next.setTooltip("Next preset");
    name.setTooltip("Presets: load, save, delete. The name is marked * once a control differs from the saved preset.");
    save.setTooltip("Save the current settings to the current preset name (a built-in preset asks for a new name).");
    prev.onClick = [this] { presets.step(-1); };
    next.onClick = [this] { presets.step(1); };
    name.onClick = [this] { showMenu(); };
    save.onClick = [this] {
        if (presets.isFactory(presets.currentName()))
            promptSaveAs();
        else
            presets.save(presets.currentName());
        refresh();
    };
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
    name.setLabel(presets.currentName() + (presets.isModified() ? " *" : ""));
}

void PresetBar::showMenu()
{
    juce::PopupMenu menu;
    int id = 1;
    std::vector<juce::String> byId{{}};
    menu.addSectionHeader("Factory");
    for (const auto& n : presets.factoryNames())
    {
        menu.addItem(id++, n, true, n == presets.currentName());
        byId.push_back(n);
    }
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
    const int saveAs = id++, del = id++, reveal = id++;
    menu.addItem(saveAs, "Save as...");
    menu.addItem(del, "Delete current", !presets.isFactory(presets.currentName()));
    menu.addItem(reveal, "Show presets folder");

    menu.showMenuAsync(juce::PopupMenu::Options().withTargetComponent(&name), [this, byId, saveAs, del, reveal](int r) {
        if (r <= 0)
            return;
        if (r == saveAs)
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

void PresetBar::promptSaveAs()
{
    auto* w = new juce::AlertWindow("Save preset", "Name for the preset:", juce::MessageBoxIconType::NoIcon, &parent);
    w->addTextEditor("name", presets.isFactory(presets.currentName()) ? "" : presets.currentName());
    w->addButton("Save", 1, juce::KeyPress(juce::KeyPress::returnKey));
    w->addButton("Cancel", 0, juce::KeyPress(juce::KeyPress::escapeKey));
    w->enterModalState(true, juce::ModalCallbackFunction::create([this, w](int r) {
                           if (r == 1)
                               presets.save(w->getTextEditorContents("name"));
                           refresh();
                       }),
                       true);
}
} // namespace rv::ui
