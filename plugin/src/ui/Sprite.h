#pragma once

#include "ui/Assets.h"
#include "ui/Design.h"

namespace rv::ui
{
// An image with its drop shadow baked in, cached per pixel size: the knob and fader caps move over the panel, so the
// shadow is composed once and the cap is blitted. `pixelsPerDesign` converts the shadow's design-pixel numbers.
class ShadowedSprite
{
  public:
    struct Result
    {
        const juce::Image& image;
        int margin; // pixels added on every side
    };

    Result get(ScaledImages& images, layout::Image id, int pw, int ph, design::Shadow shadow, float pixelsPerDesign)
    {
        if (!image.isValid() || key != std::make_tuple((int)id, pw, ph))
        {
            key = std::make_tuple((int)id, pw, ph);
            const auto blur = juce::jmax(1, juce::roundToInt(shadow.radius * pixelsPerDesign * 2.0f));
            const auto offX = shadow.offset * design::shadowDirX * pixelsPerDesign;
            const auto offY = shadow.offset * design::shadowDirY * pixelsPerDesign;
            margin = blur + (int)std::ceil(std::max(std::abs(offX), std::abs(offY))) + 2;
            image = juce::Image(juce::Image::ARGB, pw + 2 * margin, ph + 2 * margin, true);
            juce::Graphics g(image);
            const auto& cap = images.get(id, pw, ph);
            juce::Image sprite(juce::Image::ARGB, pw + 2 * margin, ph + 2 * margin, true);
            {
                juce::Graphics sg(sprite);
                sg.drawImageAt(cap, margin, margin);
            }
            juce::DropShadow(juce::Colours::black.withAlpha(shadow.alpha), blur,
                             {juce::roundToInt(offX), juce::roundToInt(offY)})
                .drawForImage(g, sprite);
            g.setOpacity(1.0f); // drawForImage leaves the shadow colour's 50% alpha, which images are drawn with
            g.drawImageAt(cap, margin, margin);
        }
        return {image, margin};
    }

  private:
    juce::Image image;
    std::tuple<int, int, int> key{-1, 0, 0};
    int margin = 0;
};
} // namespace rv::ui
