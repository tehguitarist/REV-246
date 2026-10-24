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

void LevelMeter::setPeakLinear(float peak)
{
    // Thresholds on the scale, in dB relative to levelRefDbFs; the lowest LED is -12.
    const auto db = juce::Decibels::gainToDecibels(peak, -120.0f) - design::levelRefDbFs;
    int count = 0;
    for (int i = design::numLevelLeds - 1; i >= 1; --i) // from -12 up to +3
        if (db >= design::levelStepsDb[i])
            count = design::numLevelLeds - i;
    const bool ovl = peak >= 0.999f;
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
    return "Input level in dB, after the INPUT fader. OVL lights on a full-scale peak.";
}

void LevelMeter::paint(juce::Graphics& g)
{
    layer.paint(g, getLocalBounds(), juce::String(lit) + (over ? "o" : "n"), [this](juce::Graphics& lg, float scale) {
        const auto k = getScale() * 2.5f;
        for (int i = 0; i < design::numLevelLeds; ++i)
        {
            const bool on = i == 0 ? over : (design::numLevelLeds - i) <= lit;
            if (!on)
                continue;
            const auto colour = i == 0 ? design::levelRed : i == 1 ? design::levelYellow : design::levelGreen;
            const auto box = toLocal(slots[i]->bounds().reduced(3.0f));
            // A lit lamp: a brighter middle fading to the colour at the edges, in a soft halo that tapers off.
            juce::ColourGradient lamp(colour.brighter(0.7f), box.getCentre(), colour.darker(0.15f),
                                      box.getTopLeft(), true);
            drawWithGlow(lg, getLocalBounds(), scale, colour.withAlpha(0.5f), 34.0f * k, [&](juce::Graphics& g1) {
                drawWithGlow(g1, getLocalBounds(), scale, colour.withAlpha(0.8f), 12.0f * k, [&](juce::Graphics& g2) {
                    g2.setGradientFill(lamp);
                    g2.fillRoundedRectangle(box, 4.0f * k);
                });
            });
        }
    });
}
} // namespace rv::ui
