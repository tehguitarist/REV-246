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
    const auto seg = juce::jlimit(0, design::funnelSegments - 1,
                                  juce::roundToInt(juce::jlimit(0.0f, 1.0f, f) * (float)(design::funnelSegments - 1)));
    if (seg == segment)
        return;
    segment = seg;
    layer.invalidate();
    repaint();
}

void Funnel::paint(juce::Graphics& g)
{
    layer.paint(g, getLocalBounds(), juce::String(segment) + (isDimmed() ? "d" : "n"), [this](juce::Graphics& lg, float scale) {
        const auto dx = design::funnelDx[column];
        const auto alpha = contentAlpha();
        const auto top = design::funnelTopY, bottom = design::funnelBottomY;
        const auto segH = (bottom - top) / (float)design::funnelSegments;

        // The wedge's interior between two heights, inset from its edges.
        auto wedge = [&](float ya, float yb) {
            juce::Path path;
            constexpr int steps = 8;
            for (int i = 0; i <= steps; ++i)
            {
                float left, right;
                const auto y = ya + (yb - ya) * (float)i / (float)steps;
                edgesAt(y, left, right);
                const auto p = toLocal(juce::Point<float>(left + dx + design::funnelInset, y));
                i == 0 ? path.startNewSubPath(p) : path.lineTo(p);
            }
            for (int i = steps; i >= 0; --i)
            {
                float left, right;
                const auto y = ya + (yb - ya) * (float)i / (float)steps;
                edgesAt(y, left, right);
                path.lineTo(toLocal(juce::Point<float>(right + dx - design::funnelInset, y)));
            }
            path.closeSubPath();
            return path;
        };
        lg.reduceClipRegion(wedge(top, bottom)); // nothing leaves the wedge

        // The lit segment (counted from the bottom): a brighter middle fading to the red at its edges, in a soft glow.
        const auto yBottom = bottom - segH * (float)segment, yTop = yBottom - segH;
        const auto lit = wedge(yTop + design::funnelGap * 0.5f, yBottom - design::funnelGap * 0.5f);
        const auto topL = toLocal(juce::Point<float>(0.0f, yTop)).y, bottomL = toLocal(juce::Point<float>(0.0f, yBottom)).y;
        juce::ColourGradient fill(juce::Colour(0xffff2e18).withAlpha(alpha), 0.0f, topL,
                                  juce::Colour(0xffff2e18).withAlpha(alpha), 0.0f, bottomL, false);
        fill.addColour(0.5, juce::Colour(0xffff5a30).withAlpha(alpha));
        drawWithGlow(lg, getLocalBounds(), scale, juce::Colour(0xffff1a08).withAlpha(0.34f * alpha),
                     12.0f * getScale() * 2.5f, [&](juce::Graphics& cg) {
                         cg.setGradientFill(fill);
                         cg.fillPath(lit);
                     });

        // Faint lines between the segments, over everything.
        lg.setColour(juce::Colours::black.withAlpha(design::funnelLineAlpha));
        for (int i = 1; i < design::funnelSegments; ++i)
        {
            const auto y = top + segH * (float)i;
            const auto line = wedge(y - design::funnelLineThickness * 0.5f, y + design::funnelLineThickness * 0.5f);
            lg.fillPath(line);
        }
    });
}
} // namespace rv::ui
