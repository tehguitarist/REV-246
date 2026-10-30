#include "ui/Fader.h"

#include <cmath>
#include <utility>

namespace rv::ui
{
namespace
{
// Shift-drag moves a control by this fraction of the pointer's movement.
constexpr float fineDragRatio = 0.1f;
// A wheel event moves a control by at most this fraction of its range (before the step scaling).
constexpr float maxWheelDelta = 0.2f;

// A fader's hit area: the cap's width and travel, plus a margin of design pixels around it.
constexpr float travelPadX = 20.0f, travelPadTop = 10.0f, travelPadBottom = 30.0f;

juce::Rectangle<float> travelBounds(float cx, float top, float bottom)
{
    return {cx - design::faderCapW * 0.5f - travelPadX, top - design::faderCapH * 0.5f - travelPadTop,
            design::faderCapW + 2.0f * travelPadX, bottom - top + design::faderCapH + travelPadTop + travelPadBottom};
}

// A vertical drag of this many design pixels turns a knob through its full range (200 logical pixels at the default
// editor size, which is 0.4 of the design space).
constexpr float knobDragFullRange = 500.0f;

// The WIDTH halo in design pixels, relative to its box: a 270 degree arc opening at the bottom, the value in its
// middle and the label under it.
constexpr float widthArcCentreY = 50.0f, widthArcRadius = 38.0f, widthArcThickness = 4.5f;
constexpr float widthArcStartDegrees = 225.0f, widthArcSweepDegrees = 270.0f;
constexpr float widthGlowThickness = 14.0f, widthBloomThickness = 8.0f;
constexpr float widthValueFont = 25.0f, widthValueW = 100.0f, widthValueH = 36.0f;
constexpr float widthLabelFont = 18.0f, widthLabelTop = 101.0f, widthLabelH = 24.0f;

// A knob's value from its parameter's own range for the bubble and the tooltip (the callers' text takes the real value).
std::function<juce::String(float)> fromRealValue(juce::RangedAudioParameter& p, std::function<juce::String(float)> realText)
{
    return [&p, text = std::move(realText)](float normalised) { return text(p.convertFrom0to1(normalised)); };
}
} // namespace

//==============================================================================
ParameterControl::ParameterControl(juce::Rectangle<float> designBounds, juce::RangedAudioParameter& p, ValueBubble& b,
                                   std::function<juce::String(float)> text)
    : DesignComponent(designBounds), param(p), bubble(b), valueText(std::move(text)),
      attachment(p, [this](float) { repaint(); })
{
    setOpaque(false);
    setMouseCursor(juce::MouseCursor::UpDownResizeCursor);
    setTitle(p.getName(64));
    attachment.sendInitialUpdate();
}

ParameterControl::~ParameterControl()
{
    endGesture(); // the editor closed mid-drag
}

juce::String ParameterControl::getTooltip()
{
    return helpText + "\n" + valueText(param.getValue());
}

void ParameterControl::beginGesture()
{
    if (gestureOpen)
        return;
    gestureOpen = true;
    param.beginChangeGesture();
}

void ParameterControl::endGesture()
{
    if (!gestureOpen)
        return;
    gestureOpen = false;
    param.endChangeGesture();
}

void ParameterControl::setNormalised(float n)
{
    n = juce::jlimit(0.0f, 1.0f, n);
    const auto complete = !gestureOpen;
    if (complete)
        param.beginChangeGesture();
    param.setValueNotifyingHost(n);
    if (complete)
        param.endChangeGesture();
    bubble.show(*this, valueText(n));
}

std::optional<float> ParameterControl::valueAfterWheel(const juce::MouseEvent& e, const juce::MouseWheelDetails& w,
                                                       float step, float fineStep) const
{
    const auto delta = w.isReversed ? -w.deltaY : w.deltaY;
    if (juce::exactlyEqual(delta, 0.0f))
        return std::nullopt;
    return param.getValue() + juce::jlimit(-maxWheelDelta, maxWheelDelta, delta) * (e.mods.isShiftDown() ? fineStep : step);
}

//==============================================================================
Fader::Fader(ScaledImages& i, layout::Image cap, float cx, float top, float bottom, juce::RangedAudioParameter& p,
             ValueBubble& b, std::function<juce::String(float)> text)
    : ParameterControl(travelBounds(cx, top, bottom), p, b, std::move(text)), images(i), capImage(cap), centreX(cx),
      topY(top), bottomY(bottom)
{
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

float Fader::designY(const juce::MouseEvent& e) const
{
    return ((float)e.y + (float)getY()) / getScale();
}

float Fader::positionAtY(float y) const
{
    return (bottomY - y) / (bottomY - topY);
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

void Fader::setPosition(float p, bool insideGesture)
{
    if (insideGesture)
        param.setValueNotifyingHost(juce::jlimit(0.0f, 1.0f, p));
    else
        setNormalised(p);
}

void Fader::userSet(float p)
{
    p = juce::jlimit(0.0f, 1.0f, p);
    if (onUserMove)
        onUserMove(param.getValue(), p);
    setNormalised(p);
}

void Fader::resetToDefault()
{
    userSet(param.getDefaultValue());
    if (onInteractionEnd)
        onInteractionEnd();
}

void Fader::mouseDown(const juce::MouseEvent& e)
{
    if (e.mods.isAltDown())
    {
        resetToDefault();
        return;
    }
    beginGesture();
    const auto y = designY(e);
    const auto onCap = capDesignRect().contains(centreX, y);
    grabOffset = onCap ? y - capCentreY() : 0.0f;
    if (!onCap) // a click on the track brings the cap to the pointer
        userSet(positionAtY(y));
    dragStartPos = param.getValue();
    dragStartY = e.y;
}

void Fader::mouseDrag(const juce::MouseEvent& e)
{
    if (!isGestureOpen())
        return;
    if (e.mods.isShiftDown())
    {
        // Fine: a tenth of the pointer's movement, relative to where it started.
        const auto dy = (float)(dragStartY - e.y) / getScale() / (bottomY - topY);
        userSet(dragStartPos + dy * fineDragRatio);
        grabOffset = designY(e) - capCentreY(); // so releasing Shift mid-drag doesn't jump the cap to the pointer
    }
    else
    {
        userSet(positionAtY(designY(e) - grabOffset));
        dragStartPos = param.getValue(); // so pressing Shift mid-drag starts from here
        dragStartY = e.y;
    }
}

void Fader::mouseUp(const juce::MouseEvent&)
{
    if (!isGestureOpen())
        return;
    endGesture();
    if (onInteractionEnd)
        onInteractionEnd();
}

void Fader::mouseDoubleClick(const juce::MouseEvent&)
{
    resetToDefault();
}

void Fader::mouseWheelMove(const juce::MouseEvent& e, const juce::MouseWheelDetails& w)
{
    if (const auto target = valueAfterWheel(e, w, 0.1f, 0.02f))
    {
        userSet(*target);
        if (onInteractionEnd)
            onInteractionEnd();
    }
}

//==============================================================================
DragKnob::DragKnob(juce::Rectangle<float> designBounds, juce::RangedAudioParameter& p, ValueBubble& b,
                   std::function<juce::String(float)> text, float wheel, float wheelFine)
    : ParameterControl(designBounds, p, b, fromRealValue(p, std::move(text))), wheelStep(wheel), wheelFineStep(wheelFine)
{
}

void DragKnob::resetToDefault()
{
    setNormalised(param.getDefaultValue());
}

void DragKnob::mouseDown(const juce::MouseEvent& e)
{
    if (e.mods.isAltDown())
    {
        resetToDefault();
        return;
    }
    beginGesture();
    dragStartNorm = param.getValue();
    dragStartY = e.y;
}

void DragKnob::mouseDrag(const juce::MouseEvent& e)
{
    if (!isGestureOpen())
        return;
    const auto dy = (float)(dragStartY - e.y) / (knobDragFullRange * getScale()) * (e.mods.isShiftDown() ? fineDragRatio : 1.0f);
    setNormalised(dragStartNorm + dy);
}

void DragKnob::mouseUp(const juce::MouseEvent&)
{
    endGesture();
}

void DragKnob::mouseDoubleClick(const juce::MouseEvent&)
{
    resetToDefault();
}

void DragKnob::mouseWheelMove(const juce::MouseEvent& e, const juce::MouseWheelDetails& w)
{
    if (const auto target = valueAfterWheel(e, w, wheelStep, wheelFineStep))
        setNormalised(*target);
}

//==============================================================================
PanKnob::PanKnob(ScaledImages& i, juce::RangedAudioParameter& p, ValueBubble& b, std::function<juce::String(float)> text)
    : DragKnob(layout::panKnob.bounds().expanded(20.0f), p, b, std::move(text), 0.1f, 0.02f), images(i)
{
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

//==============================================================================
WidthKnob::WidthKnob(juce::Rectangle<float> designBounds, juce::RangedAudioParameter& p, ValueBubble& b,
                     std::function<juce::String(float)> text)
    : DragKnob(designBounds, p, b, std::move(text), 0.02f, 0.002f)
{
}

void WidthKnob::paint(juce::Graphics& g)
{
    const auto s = getScale();
    const auto alpha = contentAlpha() * (isEnabled() ? 1.0f : 0.3f);
    const auto box = getDesignBounds();
    const auto centre = toLocal(juce::Point<float>{box.getCentreX(), box.getY() + widthArcCentreY});
    const auto radius = widthArcRadius * s;
    const auto thickness = widthArcThickness * s;
    const auto start = juce::degreesToRadians(widthArcStartDegrees);
    const auto sweep = juce::degreesToRadians(widthArcSweepDegrees);
    const auto value = juce::jlimit(0.0f, 1.0f, param.getValue());
    const auto stroke = [&](float width) {
        return juce::PathStrokeType(width, juce::PathStrokeType::curved, juce::PathStrokeType::rounded);
    };

    juce::Path background;
    background.addCentredArc(centre.x, centre.y, radius, radius, 0.0f, start, start + sweep, true);
    g.setColour(design::gold.withAlpha(0.30f * alpha));
    g.strokePath(background, stroke(thickness));

    if (value > 0.0001f)
    {
        juce::Path arc;
        arc.addCentredArc(centre.x, centre.y, radius, radius, 0.0f, start, start + sweep * value, true);
        g.setColour(design::gold.withAlpha(0.14f * alpha));
        g.strokePath(arc, stroke(widthGlowThickness * s));
        g.setColour(design::gold.withAlpha(0.30f * alpha));
        g.strokePath(arc, stroke(widthBloomThickness * s));
        g.setColour(design::labelSelected.withAlpha(0.95f * alpha));
        g.strokePath(arc, stroke(thickness));
    }

    g.setFont(labelFont(widthValueFont * s));
    g.setColour(design::labelSelected.withAlpha(alpha));
    const auto valueArea = juce::Rectangle<float>(widthValueW * s, widthValueH * s).withCentre(centre);
    g.drawFittedText(juce::String(juce::roundToInt(param.convertFrom0to1(param.getValue()))), valueArea.toNearestInt(),
                     juce::Justification::centred, 1);

    g.setFont(labelFont(widthLabelFont * s));
    g.setColour(design::gold.withAlpha(0.85f * alpha));
    const auto labelArea = toLocal({box.getX(), box.getY() + widthLabelTop, box.getWidth(), widthLabelH});
    g.drawFittedText("WIDTH", labelArea.toNearestInt(), juce::Justification::centred, 1);
}

void WidthKnob::mouseDown(const juce::MouseEvent& e)
{
    if (isEnabled()) // a disabled control still receives the mouse
        DragKnob::mouseDown(e);
}

void WidthKnob::mouseWheelMove(const juce::MouseEvent& e, const juce::MouseWheelDetails& w)
{
    if (isEnabled())
        DragKnob::mouseWheelMove(e, w);
}

void WidthKnob::mouseDoubleClick(const juce::MouseEvent& e)
{
    // Option/Alt-click already set the default on its mouse down.
    if (!isEnabled() || e.mods.isAltDown())
        return;
    endGesture(); // the second click of the double-click opened one
    if (onDoubleClick)
        onDoubleClick();
}
} // namespace rv::ui
