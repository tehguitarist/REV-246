#include "ui/Fader.h"

namespace rv::ui
{
namespace
{
constexpr float dragTravelPerFull = 1.0f; // a drag of the full travel moves the fader by the full range
juce::Rectangle<float> travelBounds(float cx, float top, float bottom)
{
    return {cx - design::faderCapW * 0.5f - 20.0f, top - design::faderCapH * 0.5f - 10.0f, design::faderCapW + 40.0f,
            bottom - top + design::faderCapH + 40.0f};
}
} // namespace

Fader::Fader(ScaledImages& i, layout::Image cap, float cx, float top, float bottom, juce::RangedAudioParameter& p,
             ValueBubble& b, std::function<juce::String(float)> text)
    : DesignComponent(travelBounds(cx, top, bottom)), images(i), capImage(cap), centreX(cx), topY(top), bottomY(bottom),
      param(p), bubble(b), valueText(std::move(text)), attachment(p, [this](float) { repaint(); })
{
    setOpaque(false);
    setMouseCursor(juce::MouseCursor::UpDownResizeCursor);
    setTitle(p.getName(64));
    attachment.sendInitialUpdate();
}

float Fader::capCentreY() const
{
    return bottomY - juce::jlimit(0.0f, 1.0f, param.getValue()) * (bottomY - topY);
}

juce::Rectangle<float> Fader::capDesignRect() const
{
    return {centreX - design::faderCapW * 0.5f, capCentreY() - design::faderCapH * 0.5f, design::faderCapW,
            design::faderCapH};
}

juce::String Fader::getTooltip()
{
    return helpText + "\n" + valueText(param.getValue());
}

void Fader::paint(juce::Graphics& g)
{
    const auto scale = physicalScale(g);
    const auto placed = toLocal(capDesignRect()) * scale;
    const auto pw = juce::jmax(1, juce::roundToInt(placed.getWidth()));
    const auto ph = juce::jmax(1, juce::roundToInt(placed.getHeight()));
    const auto sp = sprite.get(images, capImage, pw, ph, design::faderShadow, (float)pw / design::faderCapW);
    g.drawImageTransformed(sp.image, juce::AffineTransform::translation(std::round(placed.getX()) - (float)sp.margin,
                                                                      std::round(placed.getY()) - (float)sp.margin)
                                         .scaled(1.0f / scale));
}

void Fader::setPosition(float p, bool asGesture)
{
    p = juce::jlimit(0.0f, 1.0f, p);
    if (!asGesture)
        param.beginChangeGesture();
    param.setValueNotifyingHost(p);
    if (!asGesture)
        param.endChangeGesture();
}

void Fader::userSet(float p, bool gesture)
{
    p = juce::jlimit(0.0f, 1.0f, p);
    const auto old = param.getValue();
    if (onUserMove)
        onUserMove(old, p);
    setPosition(p, gesture);
    bubble.show(*this, valueText(p));
}

void Fader::mouseDown(const juce::MouseEvent& e)
{
    if (e.mods.isCommandDown() || e.mods.isAltDown())
    {
        mouseDoubleClick(e);
        return;
    }
    param.beginChangeGesture();
    dragging = true;
    const auto yDesign = ((float)e.y + (float)getY()) / getScale();
    const auto onCap = capDesignRect().contains(centreX, yDesign);
    grabOffset = onCap ? yDesign - capCentreY() : 0.0f;
    dragStartPos = param.getValue();
    dragStartY = e.y;
    if (!onCap) // a click on the track brings the cap to the pointer
        userSet((bottomY - yDesign) / (bottomY - topY), true);
}

void Fader::mouseDrag(const juce::MouseEvent& e)
{
    if (!dragging)
        return;
    if (e.mods.isShiftDown())
    {
        // Fine: a tenth of the pointer's movement, relative to where it started.
        const auto dy = (float)(dragStartY - e.y) / getScale() / (bottomY - topY);
        userSet(dragStartPos + dy * 0.1f * dragTravelPerFull, true);
    }
    else
    {
        const auto yDesign = ((float)e.y + (float)getY()) / getScale() - grabOffset;
        userSet((bottomY - yDesign) / (bottomY - topY), true);
        dragStartPos = param.getValue(); // so releasing Shift mid-drag doesn't jump
        dragStartY = e.y;
    }
}

void Fader::mouseUp(const juce::MouseEvent&)
{
    if (!dragging)
        return;
    dragging = false;
    param.endChangeGesture();
    if (onDragEnd)
        onDragEnd();
}

void Fader::mouseDoubleClick(const juce::MouseEvent&)
{
    userSet(param.getDefaultValue(), false);
    if (onDragEnd)
        onDragEnd();
}

void Fader::mouseWheelMove(const juce::MouseEvent& e, const juce::MouseWheelDetails& w)
{
    const auto delta = (w.isReversed ? -w.deltaY : w.deltaY);
    if (delta == 0.0f)
        return;
    userSet(param.getValue() + juce::jlimit(-0.2f, 0.2f, delta) * (e.mods.isShiftDown() ? 0.02f : 0.1f), false);
    if (onDragEnd)
        onDragEnd();
}

//==============================================================================
PanKnob::PanKnob(ScaledImages& i, juce::RangedAudioParameter& p, ValueBubble& b, std::function<juce::String(float)> text)
    : DesignComponent(layout::panKnob.bounds().expanded(20.0f)), images(i), param(p), bubble(b), valueText(std::move(text)),
      attachment(p, [this](float) { repaint(); })
{
    setOpaque(false);
    setMouseCursor(juce::MouseCursor::UpDownResizeCursor);
    setTitle(p.getName(64));
    attachment.sendInitialUpdate();
}

juce::String PanKnob::getTooltip()
{
    return helpText + "\n" + valueText(param.convertFrom0to1(param.getValue()));
}

void PanKnob::paint(juce::Graphics& g)
{
    const auto scale = physicalScale(g);
    const auto cap = toLocal(layout::panKnob.bounds());
    const auto pw = juce::jmax(1, juce::roundToInt(cap.getWidth() * scale));
    const auto ph = juce::jmax(1, juce::roundToInt(cap.getHeight() * scale));
    const auto sp = sprite.get(images, layout::Image::knob, pw, ph, design::knobShadow, (float)pw / layout::panKnob.w,
                          design::knobBrightness, design::knobContrast);
    const auto angle = juce::degreesToRadians((param.getValue() * 2.0f - 1.0f) * design::panSweepDegrees);

    // The shadow rotates with the baked sprite; at the knob's size the difference from a fixed light is invisible.
    g.drawImageTransformed(sp.image,
                           juce::AffineTransform::translation(-0.5f * (float)sp.image.getWidth(), -0.5f * (float)sp.image.getHeight())
                               .rotated(angle)
                               .translated(cap.getCentreX() * scale, cap.getCentreY() * scale)
                               .scaled(1.0f / scale));
}

void PanKnob::setNormalised(float n, bool gesture)
{
    n = juce::jlimit(0.0f, 1.0f, n);
    if (!gesture)
        param.beginChangeGesture();
    param.setValueNotifyingHost(n);
    if (!gesture)
        param.endChangeGesture();
    bubble.show(*this, valueText(param.convertFrom0to1(n)));
}

void PanKnob::mouseDown(const juce::MouseEvent& e)
{
    if (e.mods.isCommandDown() || e.mods.isAltDown())
    {
        mouseDoubleClick(e);
        return;
    }
    param.beginChangeGesture();
    dragging = true;
    dragStartNorm = param.getValue();
    dragStartY = e.y;
}

void PanKnob::mouseDrag(const juce::MouseEvent& e)
{
    if (!dragging)
        return;
    // About 200 logical pixels (at the default size) for the full sweep, x0.1 with Shift.
    const auto full = 200.0f * getScale() / 0.4f;
    const auto dy = (float)(dragStartY - e.y) / full * (e.mods.isShiftDown() ? 0.1f : 1.0f);
    setNormalised(dragStartNorm + dy, true);
}

void PanKnob::mouseUp(const juce::MouseEvent&)
{
    if (!dragging)
        return;
    dragging = false;
    param.endChangeGesture();
}

void PanKnob::mouseDoubleClick(const juce::MouseEvent&)
{
    setNormalised(param.getDefaultValue(), false);
}

void PanKnob::mouseWheelMove(const juce::MouseEvent& e, const juce::MouseWheelDetails& w)
{
    const auto delta = w.isReversed ? -w.deltaY : w.deltaY;
    if (delta != 0.0f)
        setNormalised(param.getValue() + juce::jlimit(-0.2f, 0.2f, delta) * (e.mods.isShiftDown() ? 0.02f : 0.1f), false);
}
} // namespace rv::ui
