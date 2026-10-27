#include "ui/Buttons.h"

namespace rv::ui
{
ProgramButton::ProgramButton(ScaledImages& i, layout::Slot slot, layout::Image out, layout::Image in)
    : DesignComponent(slot.bounds()), images(i), outImage(out), inImage(in)
{
    setMouseCursor(juce::MouseCursor::PointingHandCursor);
}

void ProgramButton::paint(juce::Graphics& g)
{
    drawFitted(g, images, held ? inImage : outImage, getLocalBounds().toFloat());
}

void ProgramButton::mouseDown(const juce::MouseEvent&)
{
    held = true;
    repaint();
    if (onClick)
        onClick();
}

void ProgramButton::mouseUp(const juce::MouseEvent&)
{
    held = false;
    repaint();
}

Led::Led(ScaledImages& i, layout::Slot slot, layout::Image on_, layout::Image off)
    : DesignComponent(slot.bounds()), images(i), onImage(on_), offImage(off)
{
    setInterceptsMouseClicks(false, false);
}

void Led::setOn(bool s)
{
    if (s == on)
        return;
    on = s;
    repaint();
}

void Led::paint(juce::Graphics& g)
{
    drawFitted(g, images, on ? onImage : offImage, getLocalBounds().toFloat());
}

PowerButton::PowerButton(ScaledImages& i) : DesignComponent(layout::powerButton.bounds()), images(i)
{
    setMouseCursor(juce::MouseCursor::PointingHandCursor);
    setTitle("Power");
}

void PowerButton::setOn(bool s)
{
    if (s == on)
        return;
    on = s;
    repaint();
}

juce::String PowerButton::getTooltip()
{
    return on ? "Power: on. Press to bypass the plugin." : "Power: bypassed. Press to switch the plugin on.";
}

void PowerButton::paint(juce::Graphics& g)
{
    using layout::Image;
    const auto id = held ? (on ? Image::powerInOn : Image::powerInOff) : (on ? Image::powerOutIn : Image::powerOutOff);
    drawFitted(g, images, id, getLocalBounds().toFloat());
}

void PowerButton::mouseDown(const juce::MouseEvent&)
{
    held = true;
    repaint();
}

void PowerButton::mouseUp(const juce::MouseEvent& e)
{
    held = false;
    if (getLocalBounds().contains(e.getPosition()))
    {
        if (onToggle)
            onToggle(!on);
        else
            setOn(!on);
    }
    repaint();
}
} // namespace rv::ui
