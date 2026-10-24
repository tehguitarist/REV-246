#include "ui/Joystick.h"

#include "ui/Design.h"

namespace rv::ui
{
Joystick::Joystick(ScaledImages& i, layout::Slot slot) : DesignComponent(slot.bounds()), images(i)
{
    setMouseCursor(juce::MouseCursor::UpDownResizeCursor);
}

void Joystick::paint(juce::Graphics& g)
{
    const auto id = pushed > 0 ? layout::Image::joyUp : pushed < 0 ? layout::Image::joyDown : layout::Image::joyMid;
    drawFitted(g, images, id, getLocalBounds().toFloat());
}

void Joystick::mouseDown(const juce::MouseEvent& e)
{
    pushed = e.y < getHeight() / 2 ? 1 : -1;
    repeats = 0;
    repaint();
    fire();
    startTimer(repeatMs);
}

void Joystick::mouseUp(const juce::MouseEvent&)
{
    stopTimer();
    pushed = 0;
    repaint();
}

void Joystick::mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails& w)
{
    const auto delta = w.isReversed ? -w.deltaY : w.deltaY;
    if (delta != 0.0f && onStep)
        onStep(delta > 0.0f ? 1 : -1, 0);
}

void Joystick::timerCallback()
{
    ++repeats;
    fire();
}

void Joystick::fire()
{
    if (onStep)
        onStep(pushed, repeats * repeatMs >= 1000 ? repeats : 0);
}

//==============================================================================
namespace
{
// Linear interpolation (and extrapolation at the ends) of the measured wedge edges.
void edgesAt(float y, float& left, float& right)
{
    constexpr auto n = std::size(design::funnelEdges);
    size_t i = 1;
    while (i < n - 1 && design::funnelEdges[i].y < y)
        ++i;
    const auto& p = design::funnelEdges[i - 1];
    const auto& q = design::funnelEdges[i];
    const auto t = (y - p.y) / (q.y - p.y);
    left = p.left + (q.left - p.left) * t;
    right = p.right + (q.right - p.right) * t;
}
} // namespace

Funnel::Funnel(int col)
    : DesignComponent(juce::Rectangle<float>(design::funnelEdges[0].left + design::funnelDx[col] - 30.0f, 500.0f,
                                             design::funnelEdges[0].right - design::funnelEdges[0].left + 60.0f, 560.0f)),
      column(col)
{
    setOpaque(false);
    setInterceptsMouseClicks(false, false);
}

void Funnel::setFraction(float f)
{
    f = juce::jlimit(0.0f, 1.0f, f);
    if (std::abs(f - fraction) < 1e-4f)
        return;
    fraction = f;
    layer.invalidate();
    repaint();
}

void Funnel::paint(juce::Graphics& g)
{
    layer.paint(g, getLocalBounds(), juce::String(fraction, 4) + (isDimmed() ? "d" : "n"),
                [this](juce::Graphics& lg, float scale) {
                    const auto y = design::funnelMinY + (design::funnelMaxY - design::funnelMinY) * fraction;
                    float left, right;
                    edgesAt(y, left, right);
                    const auto dx = design::funnelDx[column];
                    const auto rect = juce::Rectangle<float>(left + dx + design::funnelInset, y - design::funnelLineHeight * 0.5f,
                                                             right - left - 2.0f * design::funnelInset, design::funnelLineHeight);
                    const auto alpha = contentAlpha();
                    drawWithGlow(lg, getLocalBounds(), scale, juce::Colour(0xffff2010).withAlpha(0.8f * alpha),
                                 14.0f * getScale() * 2.5f, [&](juce::Graphics& cg) {
                                     cg.setColour(juce::Colour(0xffff3a22).withAlpha(alpha));
                                     cg.fillRoundedRectangle(toLocal(rect), 3.0f * getScale() * 2.5f);
                                     cg.setColour(juce::Colour(0xffffb0a0).withAlpha(0.7f * alpha));
                                     cg.fillRoundedRectangle(toLocal(rect.reduced(1.0f, 3.0f)), 1.5f);
                                 });
                });
}
} // namespace rv::ui
