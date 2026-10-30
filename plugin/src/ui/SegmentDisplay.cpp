#include "ui/SegmentDisplay.h"

#include <array>
#include <cmath>

namespace rv::ui
{
namespace
{
constexpr unsigned sa = 1, sb = 2, sc = 4, sd = 8, se = 16, sf = 32, sg = 64;
constexpr int numSegments = 7;

// Of a digit's pitch: how wide the digit is, how far it sits in from the pitch's left edge, and where its point sits.
constexpr float digitWidth = 0.66f, digitLeftPad = 0.06f, pointOffset = 0.15f;
constexpr float pointRadius = 0.09f; // of the digit's width
constexpr float ghostBoost = 1.3f;   // the unlit segments are drawn a little stronger than segmentGhostAlpha
// The three glow halos around the lit segments, widest and faintest first: radius in pixels at the default editor
// size, and alpha.
constexpr float glowRadii[] = {64.0f, 30.0f, 8.0f};
constexpr float glowAlphas[] = {0.5f, 0.7f, 0.8f};
constexpr int dragPixelsPerStep = 6; // vertical drag, in logical pixels

// Hexagonal segments in a digit cell of width w and height h, thickness t.
juce::Path hSegment(float x0, float x1, float yc, float t)
{
    juce::Path p;
    p.startNewSubPath(x0, yc);
    p.lineTo(x0 + t * 0.5f, yc - t * 0.5f);
    p.lineTo(x1 - t * 0.5f, yc - t * 0.5f);
    p.lineTo(x1, yc);
    p.lineTo(x1 - t * 0.5f, yc + t * 0.5f);
    p.lineTo(x0 + t * 0.5f, yc + t * 0.5f);
    p.closeSubPath();
    return p;
}

juce::Path vSegment(float xc, float y0, float y1, float t)
{
    juce::Path p;
    p.startNewSubPath(xc, y0);
    p.lineTo(xc + t * 0.5f, y0 + t * 0.5f);
    p.lineTo(xc + t * 0.5f, y1 - t * 0.5f);
    p.lineTo(xc, y1);
    p.lineTo(xc - t * 0.5f, y1 - t * 0.5f);
    p.lineTo(xc - t * 0.5f, y0 + t * 0.5f);
    p.closeSubPath();
    return p;
}

// The verticals run the full height of the digit (a "1" is as tall as an "8"), mitred at the corners, and the
// horizontals sit between them, shorter by a segment thickness on each side.
std::array<juce::Path, numSegments> segmentPaths(float w, float h)
{
    const auto t = w * 0.15f, gap = t * 0.7f, half = t * 0.5f;
    const auto mid = h * 0.5f;
    const auto hx0 = t + gap * 0.2f, hx1 = w - t - gap * 0.2f;
    return {hSegment(hx0, hx1, half, t),                              // a
            vSegment(w - half, 0.0f, mid - gap * 0.35f, t),           // b
            vSegment(w - half, mid + gap * 0.35f, h, t),              // c
            hSegment(hx0, hx1, h - half, t),                          // d
            vSegment(half, mid + gap * 0.35f, h, t),                  // e
            vSegment(half, 0.0f, mid - gap * 0.35f, t),               // f
            hSegment(hx0, hx1, mid, t)};                              // g
}
} // namespace

unsigned SegmentDisplay::segmentsFor(juce::juce_wchar ch)
{
    switch (juce::CharacterFunctions::toUpperCase(ch))
    {
    case '0': case 'O': return sa | sb | sc | sd | se | sf;
    case '1': return sb | sc;
    case '2': return sa | sb | sd | se | sg;
    case '3': return sa | sb | sc | sd | sg;
    case '4': return sb | sc | sf | sg;
    case '5': case 'S': return sa | sc | sd | sf | sg;
    case '6': return sa | sc | sd | se | sf | sg;
    case '7': return sa | sb | sc;
    case '8': return sa | sb | sc | sd | se | sf | sg;
    case '9': return sa | sb | sc | sd | sf | sg;
    case 'A': return sa | sb | sc | se | sf | sg;
    case 'B': return sc | sd | se | sf | sg;
    case 'C': return sa | sd | se | sf;
    case 'D': return sb | sc | sd | se | sg;
    case 'E': return sa | sd | se | sf | sg;
    case 'F': return sa | se | sf | sg;
    case 'G': return sa | sc | sd | se | sf;
    case 'H': return sb | sc | se | sf | sg;
    case 'I': return se | sf;
    case 'L': return sd | se | sf;
    case 'M': return sa | sb | sc | se | sf;
    case 'N': return sc | se | sg;
    case 'P': return sa | sb | se | sf | sg;
    case 'R': return se | sg;
    case 'T': return sd | se | sf | sg;
    case 'U': return sb | sc | sd | se | sf;
    case '-': return sg;
    default: return 0;
    }
}

SegmentDisplay::SegmentDisplay(juce::Rectangle<float> box, int digitCount, std::vector<Field> fieldList)
    : DesignComponent(box), numDigits(digitCount), fields(std::move(fieldList)), tips(fields.size())
{
    setOpaque(false);
}

void SegmentDisplay::setText(size_t field, const juce::String& text)
{
    if (fields[field].text == text)
        return;
    fields[field].text = text;
    layer.invalidate();
    repaint();
}

void SegmentDisplay::setFieldTooltip(size_t field, const juce::String& text)
{
    tips[field] = text;
}

std::vector<SegmentDisplay::Digit> SegmentDisplay::layoutDigits() const
{
    std::vector<Digit> out((size_t)numDigits);
    for (const auto& field : fields)
    {
        std::vector<Digit> run; // the text as digits, each '.' lighting the point of the digit before it
        for (auto ch : field.text)
        {
            if (ch == '.')
            {
                if (!run.empty())
                    run.back().point = true;
            }
            else
                run.push_back({segmentsFor(ch), false});
        }
        const auto offset = field.leftAlign ? 0 : juce::jmax(0, field.count - (int)run.size());
        for (size_t i = 0; i < run.size() && (int)i + offset < field.count; ++i)
            out[(size_t)(field.first + offset) + i] = run[i];
    }
    return out;
}

void SegmentDisplay::paint(juce::Graphics& g)
{
    // The text is not in the key: setText() invalidates the layer.
    layer.paint(g, getLocalBounds(), isDimmed() ? "d" : "n", [this](juce::Graphics& lg, float scale) {
        const auto s = getScale();
        const auto area = getLocalBounds().toFloat().reduced(design::displayInsetX * s, design::displayInsetY * s);
        const auto pitch = area.getWidth() / (float)numDigits;
        const auto w = pitch * digitWidth, h = area.getHeight();
        const auto paths = segmentPaths(w, h);
        const auto shear = std::tan(juce::degreesToRadians(design::displaySlantDegrees));
        const auto cells = layoutDigits();
        const auto pointR = w * pointRadius;
        const auto alpha = contentAlpha();

        auto place = [&](int i) {
            // Slanted about the bottom of the digit.
            return juce::AffineTransform::shear(-shear, 0.0f)
                .translated(area.getX() + pitch * (float)i + shear * h + pitch * digitLeftPad, area.getY());
        };
        auto fillPoint = [&](juce::Graphics& cg, int i) {
            const auto x = area.getX() + pitch * (float)i + w + pitch * pointOffset;
            const auto y = area.getY() + h - pointR;
            cg.fillEllipse(x - pointR, y - pointR, 2.0f * pointR, 2.0f * pointR);
        };

        // Ghost segments (everything that could light), then the lit ones with their glow.
        lg.setColour(design::segmentGlow.withAlpha(design::segmentGhostAlpha * alpha * ghostBoost)); // ghosts stay red
        for (int i = 0; i < numDigits; ++i)
        {
            for (const auto& p : paths)
                lg.fillPath(p, place(i));
            fillPoint(lg, i);
        }

        // Uniformly bright orange-red segments with a thick red glow that tapers off: three halos, wide and faint to
        // tight and strong.
        auto lit = [&](juce::Graphics& cg) {
            cg.setColour(design::segmentLit.withAlpha(alpha));
            for (int i = 0; i < numDigits; ++i)
            {
                for (int seg = 0; seg < numSegments; ++seg)
                    if (cells[(size_t)i].segments & (1u << seg))
                        cg.fillPath(paths[(size_t)seg], place(i));
                if (cells[(size_t)i].point)
                    fillPoint(cg, i);
            }
        };
        const auto k = s * 2.5f; // 1 at the default editor size (scale 0.4): the glow radii are in its pixels
        const auto halo = [&](int n) { return design::segmentGlow.withAlpha(glowAlphas[n] * alpha); };
        drawWithGlow(lg, getLocalBounds(), scale, halo(0), glowRadii[0] * k, [&](juce::Graphics& g1) {
            drawWithGlow(g1, getLocalBounds(), scale, halo(1), glowRadii[1] * k, [&](juce::Graphics& g2) {
                drawWithGlow(g2, getLocalBounds(), scale, halo(2), glowRadii[2] * k, lit);
            });
        });
    });
}

int SegmentDisplay::fieldAt(int x) const
{
    const auto area = getLocalBounds().toFloat().reduced(design::displayInsetX * getScale(), 0.0f);
    const auto digit = (int)juce::jlimit(0.0f, (float)numDigits - 0.01f,
                                         ((float)x - area.getX()) / area.getWidth() * (float)numDigits);
    for (size_t i = 0; i < fields.size(); ++i)
        if (digit >= fields[i].first && digit < fields[i].first + fields[i].count)
            return (int)i;
    return -1;
}

void SegmentDisplay::mouseDown(const juce::MouseEvent& e)
{
    dragField = -1;
    const auto i = fieldAt(e.x);

    if (e.mods.isAltDown())
    {
        if (i >= 0 && fields[(size_t)i].onReset)
            fields[(size_t)i].onReset();
        return;
    }

    // Acts on the second press rather than on mouseDoubleClick, which comes after a press that has started a drag.
    if (e.getNumberOfClicks() >= 2 && i >= 0 && fields[(size_t)i].onDoubleClick)
    {
        fields[(size_t)i].onDoubleClick();
        return;
    }

    dragLastY = e.y;
    dragRemainder = 0;
    dragField = i;
}

void SegmentDisplay::mouseDrag(const juce::MouseEvent& e)
{
    if (dragField < 0 || !fields[(size_t)dragField].onNudge)
        return;
    // Up raises the value.
    dragRemainder += dragLastY - e.y;
    dragLastY = e.y;
    const auto steps = dragRemainder / dragPixelsPerStep;
    if (steps != 0)
    {
        dragRemainder -= steps * dragPixelsPerStep;
        fields[(size_t)dragField].onNudge(steps);
    }
}

void SegmentDisplay::mouseWheelMove(const juce::MouseEvent& e, const juce::MouseWheelDetails& w)
{
    const auto i = fieldAt(e.x);
    if (i < 0 || !fields[(size_t)i].onNudge)
        return;
    const auto delta = w.isReversed ? -w.deltaY : w.deltaY;
    if (!juce::exactlyEqual(delta, 0.0f))
        fields[(size_t)i].onNudge(delta > 0.0f ? 1 : -1);
}

juce::String SegmentDisplay::getTooltip()
{
    const auto i = fieldAt(getMouseXYRelative().x);
    return i >= 0 && !tips[(size_t)i].isEmpty() ? tips[(size_t)i] : SettableTooltipClient::getTooltip();
}
} // namespace rv::ui
