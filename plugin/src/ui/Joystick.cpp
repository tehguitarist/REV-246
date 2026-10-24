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
                    const auto yc = design::funnelMinY + (design::funnelMaxY - design::funnelMinY) * fraction;
                    const auto y0 = yc - design::funnelBandHeight * 0.5f, y1 = yc + design::funnelBandHeight * 0.5f;
                    const auto dx = design::funnelDx[column];
                    const auto alpha = contentAlpha();

                    // The wedge's interior between y0 and y1, inset from its edges.
                    juce::Path band;
                    constexpr int steps = 10;
                    for (int i = 0; i <= steps; ++i)
                    {
                        float left, right;
                        const auto y = y0 + (y1 - y0) * (float)i / (float)steps;
                        edgesAt(y, left, right);
                        const auto p = toLocal(juce::Point<float>(left + dx + design::funnelInset, y));
                        i == 0 ? band.startNewSubPath(p) : band.lineTo(p);
                    }
                    for (int i = steps; i >= 0; --i)
                    {
                        float left, right;
                        const auto y = y0 + (y1 - y0) * (float)i / (float)steps;
                        edgesAt(y, left, right);
                        band.lineTo(toLocal(juce::Point<float>(right + dx - design::funnelInset, y)));
                    }
                    band.closeSubPath();

                    // Brightest at the centre, fading to nothing at the ends (a smooth falloff, not a line per value).
                    const auto top = toLocal(juce::Point<float>(0.0f, y0)).y, bottom = toLocal(juce::Point<float>(0.0f, y1)).y;
                    juce::ColourGradient fade(juce::Colour(0xffff2a14).withAlpha(0.0f), 0.0f, top,
                                              juce::Colour(0xffff2a14).withAlpha(0.0f), 0.0f, bottom, false);
                    const auto edge = 0.5 - 0.5 * (double)design::funnelCoreFraction; // where the core starts
                    fade.addColour(edge * 0.45, juce::Colour(0xffff2a14).withAlpha(0.22f * alpha));
                    fade.addColour(edge, juce::Colour(0xffff3a20).withAlpha(0.95f * alpha));
                    fade.addColour(1.0 - edge, juce::Colour(0xffff3a20).withAlpha(0.95f * alpha));
                    fade.addColour(1.0 - edge * 0.45, juce::Colour(0xffff2a14).withAlpha(0.22f * alpha));

                    drawWithGlow(lg, getLocalBounds(), scale, juce::Colour(0xffff1a08).withAlpha(0.35f * alpha),
                                 14.0f * getScale() * 2.5f, [&](juce::Graphics& cg) {
                                     cg.setGradientFill(fade);
                                     cg.fillPath(band);
                                 });
                });
}
} // namespace rv::ui
