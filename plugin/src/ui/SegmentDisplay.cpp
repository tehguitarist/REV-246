#include "ui/SegmentDisplay.h"

namespace rv::ui
{
namespace
{
constexpr unsigned sa = 1, sb = 2, sc = 4, sd = 8, se = 16, sf = 32, sg = 64;

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

std::array<juce::Path, 7> segmentPaths(float w, float h)
{
    const auto t = w * 0.165f, gap = t * 0.6f, half = t * 0.5f;
    const auto mid = h * 0.5f;
    return {hSegment(half + gap * 0.2f, w - half - gap * 0.2f, half, t),            // a
            vSegment(w - half, half + gap, mid - gap * 0.55f, t),                   // b
            vSegment(w - half, mid + gap * 0.55f, h - half - gap, t),               // c
            hSegment(half + gap * 0.2f, w - half - gap * 0.2f, h - half, t),        // d
            vSegment(half, mid + gap * 0.55f, h - half - gap, t),                   // e
            vSegment(half, half + gap, mid - gap * 0.55f, t),                       // f
            hSegment(half + gap * 0.2f, w - half - gap * 0.2f, mid, t)};            // g
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

SegmentDisplay::SegmentDisplay(juce::Rectangle<float> box, int numDigits, std::vector<Field> fieldList)
    : DesignComponent(box), digits(numDigits), fields(std::move(fieldList)), tips(fields.size())
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
    std::vector<Digit> out((size_t)digits);
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
    juce::String key;
    for (const auto& field : fields)
        key << field.text << '|';
    key << (isDimmed() ? "d" : "n");

    layer.paint(g, getLocalBounds(), key, [this](juce::Graphics& lg, float scale) {
        const auto s = getScale();
        const auto area = getLocalBounds().toFloat().reduced(design::displayInsetX * s, design::displayInsetY * s);
        const auto pitch = area.getWidth() / (float)digits;
        const auto w = pitch * 0.68f, h = area.getHeight();
        const auto paths = segmentPaths(w, h);
        const auto shear = std::tan(juce::degreesToRadians(design::displaySlantDegrees));
        const auto dots = std::vector<Digit>(layoutDigits());
        const auto pointR = w * 0.09f;
        const auto alpha = contentAlpha();

        auto place = [&](int i) {
            // Slanted about the bottom of the digit.
            return juce::AffineTransform::shear(-shear, 0.0f)
                .translated(area.getX() + pitch * (float)i + shear * h + pitch * 0.06f, area.getY());
        };
        auto pointAt = [&](int i) {
            return juce::Point<float>(area.getX() + pitch * (float)i + w + pitch * 0.15f, area.getY() + h - pointR);
        };

        // Ghost segments (everything that could light), then the lit ones with their glow.
        lg.setColour(design::segmentLit.withAlpha(design::segmentGhostAlpha * alpha));
        for (int i = 0; i < digits; ++i)
        {
            for (const auto& p : paths)
                lg.fillPath(p, place(i));
            lg.fillEllipse(pointAt(i).x - pointR, pointAt(i).y - pointR, 2 * pointR, 2 * pointR);
        }

        // A bright core with a thick glow that tapers off: three halos, wide and faint to tight and strong, around
        // segments that have a hot line along their middle.
        auto lit = [&](juce::Graphics& cg) {
            for (int i = 0; i < digits; ++i)
            {
                for (int sIdx = 0; sIdx < 7; ++sIdx)
                    if (dots[(size_t)i].segments & (1u << sIdx))
                    {
                        cg.setColour(design::segmentLit.withAlpha(alpha));
                        cg.fillPath(paths[(size_t)sIdx], place(i));
                        cg.setColour(design::segmentCore.withAlpha(0.75f * alpha));
                        cg.strokePath(paths[(size_t)sIdx], juce::PathStrokeType(w * 0.045f), place(i));
                    }
                if (dots[(size_t)i].point)
                {
                    cg.setColour(design::segmentLit.withAlpha(alpha));
                    cg.fillEllipse(pointAt(i).x - pointR, pointAt(i).y - pointR, 2 * pointR, 2 * pointR);
                    cg.setColour(design::segmentCore.withAlpha(0.75f * alpha));
                    cg.fillEllipse(pointAt(i).x - pointR * 0.45f, pointAt(i).y - pointR * 0.45f, pointR * 0.9f, pointR * 0.9f);
                }
            }
        };
        const auto k = s * 2.5f;
        drawWithGlow(lg, getLocalBounds(), scale, design::segmentGlow.withAlpha(0.6f * alpha), 64.0f * k, [&](juce::Graphics& g1) {
            drawWithGlow(g1, getLocalBounds(), scale, design::segmentGlow.withAlpha(0.85f * alpha), 30.0f * k, [&](juce::Graphics& g2) {
                drawWithGlow(g2, getLocalBounds(), scale, design::segmentGlow.withAlpha(0.95f * alpha), 8.0f * k, lit);
            });
        });
    });
}

int SegmentDisplay::fieldAt(int x) const
{
    const auto area = getLocalBounds().toFloat().reduced(design::displayInsetX * getScale(), 0.0f);
    const auto digit = (int)juce::jlimit(0.0f, (float)digits - 0.01f, ((float)x - area.getX()) / area.getWidth() * (float)digits);
    for (size_t i = 0; i < fields.size(); ++i)
        if (digit >= fields[i].first && digit < fields[i].first + fields[i].count)
            return (int)i;
    return -1;
}

void SegmentDisplay::mouseDown(const juce::MouseEvent& e)
{
    dragLastY = e.y;
    dragRemainder = 0;
    const auto i = fieldAt(e.x);
    dragField = i < 0 ? 0 : (size_t)i;
}

void SegmentDisplay::mouseDrag(const juce::MouseEvent& e)
{
    if (!fields[dragField].onNudge)
        return;
    // About 6 logical pixels of vertical drag per step; up raises the value.
    dragRemainder += dragLastY - e.y;
    dragLastY = e.y;
    const auto steps = dragRemainder / 6;
    if (steps != 0)
    {
        dragRemainder -= steps * 6;
        fields[dragField].onNudge(steps);
    }
}

void SegmentDisplay::mouseWheelMove(const juce::MouseEvent& e, const juce::MouseWheelDetails& w)
{
    const auto i = fieldAt(e.x);
    if (i < 0 || !fields[(size_t)i].onNudge)
        return;
    const auto delta = w.isReversed ? -w.deltaY : w.deltaY;
    if (delta != 0.0f)
        fields[(size_t)i].onNudge(delta > 0.0f ? 1 : -1);
}

juce::String SegmentDisplay::getTooltip()
{
    const auto i = fieldAt(getMouseXYRelative().x);
    return i >= 0 && !tips[(size_t)i].isEmpty() ? tips[(size_t)i] : SettableTooltipClient::getTooltip();
}
} // namespace rv::ui
