#include "ui/Assets.h"

#include <BinaryData.h>

#include <cmath>
#include <utility>

namespace rv::ui
{
namespace
{
// Embedded resources are looked up by their original file name, so Layout.h's file list is the
// only place the names live.
std::pair<const char*, int> findResource(const juce::String& fileName)
{
    for (int i = 0; i < BinaryData::namedResourceListSize; ++i)
    {
        const auto* name = BinaryData::namedResourceList[i];
        if (fileName == BinaryData::getNamedResourceOriginalFilename(name))
        {
            int size = 0;
            const auto* data = BinaryData::getNamedResource(name, size);
            return {data, size};
        }
    }
    jassertfalse; // run plugin/tools/build_assets.sh and re-run CMake
    return {nullptr, 0};
}

// Halves an ARGB image with a 2x2 box filter. Premultiplied channels average correctly.
juce::Image halve(const juce::Image& src)
{
    const auto w = juce::jmax(1, src.getWidth() / 2), h = juce::jmax(1, src.getHeight() / 2);
    juce::Image dst(juce::Image::ARGB, w, h, false, juce::SoftwareImageType());
    const juce::Image::BitmapData s(src, juce::Image::BitmapData::readOnly);
    juce::Image::BitmapData d(dst, juce::Image::BitmapData::writeOnly);
    const auto maxX = src.getWidth() - 1, maxY = src.getHeight() - 1;

    for (int y = 0; y < h; ++y)
    {
        const auto* row0 = s.getLinePointer(2 * y);
        const auto* row1 = s.getLinePointer(juce::jmin(2 * y + 1, maxY));
        auto* out = d.getLinePointer(y);
        for (int x = 0; x < w; ++x, out += d.pixelStride)
        {
            const auto x0 = 2 * x * s.pixelStride, x1 = juce::jmin(2 * x + 1, maxX) * s.pixelStride;
            for (int c = 0; c < 4; ++c)
                out[c] = (juce::uint8)((row0[x0 + c] + row0[x1 + c] + row1[x0 + c] + row1[x1 + c] + 2) / 4);
        }
    }
    return dst;
}

// Ink height of the glyph (not the font's line box).
float glyphHeight(const juce::Font& font, const juce::String& glyph)
{
    juce::GlyphArrangement ga;
    ga.addLineOfText(font, glyph, 0.0f, 0.0f);
    juce::Path path;
    ga.createPath(path);
    return path.getBounds().getHeight();
}

// The label typeface and the font height that gives a cap height of 1, found once: labelFont() is called from paint
// code, and measuring a glyph outline per call would be wasted work.
struct LabelTypeface
{
    juce::FontOptions options;
    float heightPerCapHeight;
};

const LabelTypeface& labelTypeface()
{
    static const LabelTypeface typeface = [] {
        const auto name = juce::Font::findAllTypefaceNames().contains("Arial") ? juce::String("Arial")
                                                                                : juce::String("Liberation Sans");
        const juce::FontOptions options(name, 10.0f, juce::Font::bold);
        constexpr float reference = 100.0f;
        const auto measured = glyphHeight(juce::Font(options.withHeight(reference)), "H");
        return LabelTypeface{options, measured > 0.0f ? reference / measured : 1.0f};
    }();
    return typeface;
}
} // namespace

//==============================================================================
SourceAssets::SourceAssets()
{
    for (size_t i = 0; i < images.size(); ++i)
    {
        const auto [data, size] = findResource(layout::imageFileNames[i]);
        if (data != nullptr)
            images[i] = juce::ImageFileFormat::loadFrom(data, (size_t)size);
        jassert(images[i].isValid());
    }
}

const juce::Image& ScaledImages::get(layout::Image id, int pixelWidth, int pixelHeight)
{
    const auto key = std::make_tuple((int)id, pixelWidth, pixelHeight);
    auto it = cache.find(key);
    if (it == cache.end())
        it = cache.emplace(key, resample(sourceAssets->image(id), pixelWidth, pixelHeight)).first;
    return it->second;
}

juce::Image resample(const juce::Image& source, int width, int height)
{
    if (! source.isValid() || width <= 0 || height <= 0)
        return {};

    auto current = juce::SoftwareImageType().convert(source.convertedToFormat(juce::Image::ARGB));
    while (current.getWidth() >= 2 * width && current.getHeight() >= 2 * height)
        current = halve(current);

    juce::Image result(juce::Image::ARGB, width, height, true);
    juce::Graphics g(result);
    g.setImageResamplingQuality(juce::Graphics::highResamplingQuality);
    g.drawImage(current, 0, 0, width, height, 0, 0, current.getWidth(), current.getHeight());
    return result;
}

float physicalScale(juce::Graphics& g)
{
    return g.getInternalContext().getPhysicalPixelScaleFactor();
}

juce::Rectangle<float> fittedArea(const ScaledImages& images, layout::Image id, juce::Rectangle<float> area)
{
    const auto& src = images.source().image(id);
    return juce::RectanglePlacement(juce::RectanglePlacement::centred)
        .appliedTo(juce::Rectangle<float>((float)src.getWidth(), (float)src.getHeight()), area);
}

void drawFitted(juce::Graphics& g, ScaledImages& images, layout::Image id, juce::Rectangle<float> area)
{
    const auto scale = physicalScale(g);
    const auto placed = fittedArea(images, id, area) * scale;

    // Snap to whole device pixels so the cached image maps 1:1 and stays sharp.
    const auto x = std::round(placed.getX()), y = std::round(placed.getY());
    const auto w = juce::jmax(1, juce::roundToInt(placed.getWidth()));
    const auto h = juce::jmax(1, juce::roundToInt(placed.getHeight()));
    g.drawImageTransformed(images.get(id, w, h), juce::AffineTransform::translation(x, y).scaled(1.0f / scale));
}

void drawRotated(juce::Graphics& g, ScaledImages& images, layout::Image id, juce::Rectangle<float> area, float radians)
{
    const auto scale = physicalScale(g);
    const auto placed = fittedArea(images, id, area);
    const auto w = juce::jmax(1, juce::roundToInt(placed.getWidth() * scale));
    const auto h = juce::jmax(1, juce::roundToInt(placed.getHeight() * scale));

    juce::Graphics::ScopedSaveState save(g);
    g.setImageResamplingQuality(juce::Graphics::highResamplingQuality);
    g.drawImageTransformed(images.get(id, w, h), juce::AffineTransform::translation(-0.5f * (float)w, -0.5f * (float)h)
                                                     .scaled(1.0f / scale)
                                                     .rotated(radians)
                                                     .translated(placed.getCentre()));
}

//==============================================================================
juce::Font labelFont(float capHeight)
{
    const auto& typeface = labelTypeface();
    return juce::Font(typeface.options.withHeight(capHeight * typeface.heightPerCapHeight));
}

void drawTextAt(juce::Graphics& g, const juce::Font& font, const juce::String& text, float x, float baselineY,
                juce::Justification justification)
{
    const auto width = textWidth(font, text);
    if (justification.testFlags(juce::Justification::right))
        x -= width;
    else if (justification.testFlags(juce::Justification::horizontallyCentred))
        x -= 0.5f * width;

    juce::GlyphArrangement ga;
    ga.addLineOfText(font, text, x, baselineY);
    ga.draw(g);
}

float textWidth(const juce::Font& font, const juce::String& text)
{
    return juce::GlyphArrangement::getStringWidth(font, text);
}
} // namespace rv::ui
