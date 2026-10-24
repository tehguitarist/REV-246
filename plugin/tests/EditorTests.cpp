#include "PluginEditor.h"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("the editor opens, resizes and keeps its aspect")
{
    Rev246Processor p;
    std::unique_ptr<juce::AudioProcessorEditor> base(p.createEditor());
    auto* editor = dynamic_cast<Rev246Editor*>(base.get());
    REQUIRE(editor != nullptr);
    CHECK(editor->getWidth() == Rev246Processor::defaultEditorWidth);
    CHECK(editor->getHeight() == 612);
    editor->setSize(800, 350);
    CHECK(editor->getWidth() == 800);
}

TEST_CASE("tooltips follow the ? button")
{
    Rev246Processor p;
    std::unique_ptr<juce::AudioProcessorEditor> base(p.createEditor());
    auto* editor = dynamic_cast<Rev246Editor*>(base.get());
    REQUIRE(editor != nullptr);
    auto& help = editor->getHelpButton();
    CHECK_FALSE(help.isOn());
    CHECK(editor->tooltipAllowed(help));
    CHECK_FALSE(editor->tooltipAllowed(*editor));
    help.setOn(true);
    CHECK(editor->tooltipAllowed(*editor));
}

// Renders the panel to a PNG when REV_SNAPSHOT names a file, for checking the layout against ui/main.png by eye.
TEST_CASE("snapshot")
{
    const auto path = juce::SystemStats::getEnvironmentVariable("REV_SNAPSHOT", {});
    if (path.isEmpty())
        return;
    Rev246Processor p;
    const auto program = juce::SystemStats::getEnvironmentVariable("REV_SNAPSHOT_PROGRAM", "0").getIntValue();
    p.apvts.getParameter(rv::params::id::program)->setValueNotifyingHost(
        p.apvts.getParameter(rv::params::id::program)->convertTo0to1((float)program));
    // REV_SNAPSHOT_EXTREME=min|max puts every lever at its end, to check the funnel indicators reach the wedge's ends.
    const auto extreme = juce::SystemStats::getEnvironmentVariable("REV_SNAPSHOT_EXTREME", {});
    if (extreme.isNotEmpty())
        for (int i = 0; i < rv::params::numJoys; ++i)
            p.apvts.getParameter(rv::params::joyId((rv::params::Joy)i))->setValueNotifyingHost(extreme == "max" ? 1.0f : 0.0f);
    std::unique_ptr<juce::AudioProcessorEditor> base(p.createEditor());
    dynamic_cast<Rev246Editor&>(*base).getLevelMeter().setPeakLinear(0.7f); // lights some of the LEVEL LEDs
    base->setSize(1747, 764);
    auto img = base->createComponentSnapshot(base->getLocalBounds(), true, 1.0f);
    juce::File(path).deleteFile();
    juce::FileOutputStream out{juce::File(path)};
    juce::PNGImageFormat().writeImageToStream(img, out);
}
