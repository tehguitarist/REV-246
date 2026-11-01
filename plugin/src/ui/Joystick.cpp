#include "ui/Joystick.h"

#include "ui/Design.h"

#include <cmath>
#include <iterator>

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
    if (e.mods.isAltDown())
    {
        if (onReset)
            onReset();
        return;
    }
    pressedHalf = e.y < getHeight() / 2 ? 1 : -1;
    pushed = 0;
    push(pressedHalf);
}

void Joystick::mouseDrag(const juce::MouseEvent& e)
{
    if (pushed == 0 || e.mods.isAltDown())
        return;
    // Dragging pushes the lever the way the pointer has moved from where it was pressed (up on screen = up), once it has
    // left a small dead zone; back inside it, the lever returns to the half that was pressed.
    const auto deadZone = juce::jmax(5.0f, 0.12f * (float)getHeight());
    const auto dy = (float)e.getDistanceFromDragStartY();
    const auto direction = std::abs(dy) <= deadZone ? pressedHalf : dy < 0.0f ? 1 : -1;
    if (direction != pushed)
        push(direction);
}

void Joystick::push(int direction)
{
    pushed = direction;
    pressMs = lastTickMs = juce::Time::getMillisecondCounterHiRes();
    owed = 0.0;
    accelerating = false;
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
    if (!juce::exactlyEqual(delta, 0.0f) && onStep)
        onStep(delta > 0.0f ? 1 : -1, 1);
}

// Held steps per second, `heldSeconds` after the press: from 4 at the end of the first second, growing by a factor
// per second that is larger for longer levers, up to the whole range per second.
double Joystick::stepsPerSecond(double heldSeconds) const
{
    const auto r = (double)juce::jmax(range ? range() : 15, 1);
    const auto growth = 2.0 + std::log10(juce::jmax(r, 16.0) / 16.0);
    const auto after = juce::jmax(0.0, heldSeconds - accelerateAfterMs * 0.001);
    return juce::jmin(r, 4.0 * std::pow(growth, after));
}

void Joystick::timerCallback()
{
    const auto now = juce::Time::getMillisecondCounterHiRes();
    if (!accelerating)
    {
        if (now - pressMs < accelerateAfterMs - 1) // still the half-second repeat
        {
            fire();
            return;
        }
        accelerating = true;
        lastTickMs = now;
        startTimer(tickMs);
        fire();
        return;
    }
    owed += stepsPerSecond(0.001 * (now - pressMs)) * 0.001 * (now - lastTickMs);
    lastTickMs = now;
    const auto n = (int)owed;
    if (n > 0)
    {
        owed -= n;
        fire(n);
    }
}

void Joystick::fire(int steps)
{
    if (onStep)
        onStep(pushed, juce::jmax(1, steps));
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

// The indicator's reds: the spill, the lit segment (edge and middle), and the wedge's ambient tint.
const juce::Colour spillOuter{0xffff1a08}, spillInner{0xffff2a14}, litEdge{0xffff2e18}, litMiddle{0xffff5a30},
    ambientTint{0xffff2410};
} // namespace

Funnel::Funnel(int col)
    : DesignComponent(juce::Rectangle<float>(design::funnelEdges[0].left + design::funnelDx[col] - 70.0f, 470.0f,
                                             design::funnelEdges[0].right - design::funnelEdges[0].left + 140.0f, 620.0f)),
      column(col)
{
    setOpaque(false);
    setMouseCursor(juce::MouseCursor::UpDownResizeCursor);
}

float Funnel::fractionAt(const juce::MouseEvent& e) const
{
    const auto yDesign = ((float)e.y + (float)getY()) / getScale();
    return juce::jlimit(0.0f, 1.0f, (design::funnelBottomY - yDesign) / (design::funnelBottomY - design::funnelTopY));
}

void Funnel::mouseDown(const juce::MouseEvent& e)
{
    if (e.mods.isAltDown())
    {
        if (onReset)
            onReset();
        return;
    }
    dragging = true;
    if (onSetFraction)
        onSetFraction(fractionAt(e));
}

void Funnel::mouseDrag(const juce::MouseEvent& e)
{
    if (!dragging)
        return;
    if (onSetFraction)
        onSetFraction(fractionAt(e));
}

void Funnel::mouseUp(const juce::MouseEvent&)
{
    dragging = false;
}

void Funnel::mouseWheelMove(const juce::MouseEvent&, const juce::MouseWheelDetails& w)
{
    const auto delta = w.isReversed ? -w.deltaY : w.deltaY;
    if (!juce::exactlyEqual(delta, 0.0f) && onStep)
        onStep(delta > 0.0f ? 1 : -1);
}

void Funnel::setFraction(float f)
{
    const auto seg = juce::jlimit(0, design::funnelSegments - 1,
                                  juce::roundToInt(juce::jlimit(0.0f, 1.0f, f) * (float)(design::funnelSegments - 1)));
    if (seg == litSegment)
        return;
    litSegment = seg;
    layer.invalidate();
    repaint();
}

void Funnel::paint(juce::Graphics& g)
{
    layer.paint(g, getLocalBounds(), juce::String(litSegment) + (isDimmed() ? "d" : "n"), [this](juce::Graphics& lg, float scale) {
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

        // The lit segment (counted from the bottom).
        const auto yBottom = bottom - segH * (float)litSegment, yTop = yBottom - segH;
        const auto lit = wedge(yTop + design::funnelGap * 0.5f, yBottom - design::funnelGap * 0.5f);
        const auto outline = wedge(top, bottom);

        // A little of the indicator's glow spills out past the wedge's edges, onto the frame and the panel.
        const auto k = getScale() * 2.5f; // 1 at the default editor size (scale 0.4): the glow radii below are in its pixels
        drawWithGlow(lg, getLocalBounds(), scale, spillOuter.withAlpha(0.7f * alpha), 34.0f * k, [&](juce::Graphics& g1) {
            drawWithGlow(g1, getLocalBounds(), scale, spillInner.withAlpha(0.9f * alpha), 14.0f * k, [&](juce::Graphics& g2) {
                g2.setColour(litEdge.withAlpha(alpha)); // drawn again below; this casts the spill
                g2.fillPath(lit);
            });
        });

        lg.reduceClipRegion(outline); // everything else stays inside the wedge

        // The glow also lights the wedge itself: a red tint, strongest around the lit segment and falling away from it.
        {
            const auto centre = (((bottom - segH * ((float)litSegment + 0.5f)) - top) / (bottom - top)); // 0..1 down the wedge
            auto tint = [&](double at) { return design::funnelAmbient * std::exp(-std::pow((at - centre) / design::funnelAmbientReach, 2.0)); };
            const auto topL = toLocal(juce::Point<float>(0.0f, top)).y, bottomL = toLocal(juce::Point<float>(0.0f, bottom)).y;
            juce::ColourGradient ambient(ambientTint.withAlpha((float)tint(0.0) * alpha), 0.0f, topL,
                                         ambientTint.withAlpha((float)tint(1.0) * alpha), 0.0f, bottomL, false);
            for (int i = 1; i < 20; ++i)
                ambient.addColour(i / 20.0, ambientTint.withAlpha((float)tint(i / 20.0) * alpha));
            lg.setGradientFill(ambient);
            lg.fillPath(outline);
        }

        // The lit segment: a brighter middle fading to the red at its edges, in a soft glow.
        const auto topL = toLocal(juce::Point<float>(0.0f, yTop)).y, bottomL = toLocal(juce::Point<float>(0.0f, yBottom)).y;
        juce::ColourGradient fill(litEdge.withAlpha(alpha), 0.0f, topL,
                                  litEdge.withAlpha(alpha), 0.0f, bottomL, false);
        fill.addColour(0.5, litMiddle.withAlpha(alpha));
        drawWithGlow(lg, getLocalBounds(), scale, spillOuter.withAlpha(0.6f * alpha),
                     17.0f * k, [&](juce::Graphics& cg) {
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
