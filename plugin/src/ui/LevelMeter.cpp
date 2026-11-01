#include "ui/LevelMeter.h"

#include <juce_audio_basics/juce_audio_basics.h>

namespace rv::ui
{
namespace
{
const layout::Slot* const slots[design::numLevelLeds] = {&layout::levelOvl, &layout::levelP3, &layout::level0,
                                                         &layout::levelM3,  &layout::levelM6, &layout::levelM9,
                                                         &layout::levelM12};

juce::Rectangle<float> unionBounds()
{
    auto r = slots[0]->bounds();
    for (auto* s : slots)
        r = r.getUnion(s->bounds());
    return r.expanded(70.0f); // room for the glow
}
} // namespace

LevelMeter::LevelMeter() : DesignComponent(unionBounds())
{
    setOpaque(false);
    setInterceptsMouseClicks(false, false);
}

void LevelMeter::setLevel(float vuLinear, bool clip, float refDbFs)
{
    // Thresholds on the scale, in dB relative to the meter's 0 VU reference; the lowest LED is -12.
    const auto now = juce::Time::getMillisecondCounterHiRes();
    const auto vuDb = juce::Decibels::gainToDecibels(vuLinear, -120.0f) - refDbFs;
    if (vuDb >= heldDb)
    {
        heldDb = vuDb;
        holdUntilMs = now + holdMs;
    }
    else if (now > holdUntilMs)
        heldDb = juce::jmax(vuDb, heldDb - fallDbPerSecond * (float)(0.001 * (now - juce::jmax(lastMs, holdUntilMs))));
    lastMs = now;
    const auto db = heldDb;
    int count = 0;
    for (int i = design::numLevelLeds - 1; i >= 1; --i) // from -12 up to +3
        if (db >= design::levelStepsDb[i])
            count = design::numLevelLeds - i;
    // OVL: +6 VU, or a peak (full scale, or the unit's A/D overload), whichever comes first; held like the level.
    if (vuDb >= design::levelStepsDb[0] || clip)
        overUntilMs = now + holdMs;
    const bool ovl = now < overUntilMs;
    if (count != lit || ovl != over)
    {
        lit = count;
        over = ovl;
        layer.invalidate();
        repaint();
    }
}

juce::String LevelMeter::getTooltip()
{
    return "Input level in VU, after the INPUT fader (the 0 VU reference is set in the bottom bar). OVL lights at +6 VU, "
           "or on any peak the unit's A/D converter would clip (as its REG LED did: high frequencies and short peaks "
           "can clip before the VU reading gets there). The meter has VU ballistics (300 ms to full reading); the lights hold "
           "their peak for 0.4 s so a short hit stays readable.";
}

void LevelMeter::paint(juce::Graphics& g)
{
    layer.paint(g, getLocalBounds(), juce::String(lit) + (over ? "o" : "n"), [this](juce::Graphics& lg, float scale) {
        const auto k = getScale() * 2.5f; // 1 at the default editor size (scale 0.4): the radii below are in its pixels
        for (int i = 0; i < design::numLevelLeds; ++i)
        {
            const bool on = i == 0 ? over : (design::numLevelLeds - i) <= lit;
            if (!on)
                continue;
            const auto colour = i == 0 ? design::levelRed : i == 1 ? design::levelYellow : design::levelGreen;
            const auto box = toLocal(slots[i]->bounds().reduced(3.0f));
            // A lit lamp: a brighter middle fading to the colour at the edges, in a soft halo that tapers off.
            juce::ColourGradient lamp(colour.brighter(0.3f), box.getCentre(), colour.darker(0.25f),
                                      box.getTopLeft(), true);
            drawWithGlow(lg, getLocalBounds(), scale, colour.withAlpha(0.38f), 34.0f * k, [&](juce::Graphics& g1) {
                drawWithGlow(g1, getLocalBounds(), scale, colour.withAlpha(0.6f), 12.0f * k, [&](juce::Graphics& g2) {
                    g2.setGradientFill(lamp);
                    g2.fillRoundedRectangle(box, 4.0f * k);
                });
            });
        }
    });
}
} // namespace rv::ui
