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
    CHECK_FALSE(help.isAlwaysLit());
    CHECK(editor->tooltipAllowed(help));
    CHECK_FALSE(editor->tooltipAllowed(*editor));
    help.setOn(true);
    CHECK(editor->tooltipAllowed(*editor));
}

TEST_CASE("the stereo button offers the original mode by default and only in stereo to stereo")
{
    Rev246Processor p;
    std::unique_ptr<juce::AudioProcessorEditor> base(p.createEditor());
    auto* editor = dynamic_cast<Rev246Editor*>(base.get());
    REQUIRE(editor != nullptr);
    auto& button = editor->getStereoButton();
    CHECK(button.isEnabled());
    CHECK(button.isAlwaysLit());
    CHECK_FALSE(button.isOn());
    CHECK(button.getLabel() == "ORIGINAL");
    CHECK(button.getTooltip().contains("Original"));

    button.onToggle(true);
    editor->refreshNow();
    CHECK(button.isOn());
    CHECK(button.getLabel() == "TRUE ST");
    CHECK(p.apvts.getRawParameterValue(rv::params::id::trueStereo)->load() >= 0.5f);

    juce::AudioProcessor::BusesLayout layout;
    layout.inputBuses.add(juce::AudioChannelSet::mono());
    layout.outputBuses.add(juce::AudioChannelSet::mono());
    REQUIRE(p.setBusesLayout(layout));
    editor->refreshNow();
    CHECK_FALSE(button.isEnabled());
}

TEST_CASE("the format button selects float or 16-bit processing")
{
    Rev246Processor p;
    std::unique_ptr<juce::AudioProcessorEditor> base(p.createEditor());
    auto* editor = dynamic_cast<Rev246Editor*>(base.get());
    REQUIRE(editor != nullptr);
    auto& button = editor->getFormatButton();
    CHECK(button.isAlwaysLit());
    CHECK_FALSE(button.isOn());
    CHECK(button.getLabel() == "FLOAT");
    CHECK(button.getTooltip().contains("16-BIT"));

    button.onToggle(true);
    editor->refreshNow();
    CHECK(button.isOn());
    CHECK(button.getLabel() == "16-BIT");
    CHECK(p.apvts.getRawParameterValue(rv::params::id::sixteenBit)->load() >= 0.5f);
}

TEST_CASE("a modified preset name is italic with a star")
{
    Rev246Processor p;
    std::unique_ptr<juce::AudioProcessorEditor> base(p.createEditor());
    auto* editor = dynamic_cast<Rev246Editor*>(base.get());
    REQUIRE(editor != nullptr);
    auto& name = editor->getPresetBar().getNameButton();
    CHECK_FALSE(name.isItalic());
    CHECK_FALSE(name.getLabel().endsWithChar('*'));

    p.apvts.getParameter(rv::params::id::mix)->setValueNotifyingHost(0.9f);
    editor->refreshNow();
    CHECK(name.isItalic());
    CHECK(name.getLabel().endsWith(" *"));

    p.presets.load("Init S");
    editor->refreshNow();
    CHECK_FALSE(name.isItalic());
    CHECK_FALSE(name.getLabel().endsWithChar('*'));
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
    dynamic_cast<Rev246Editor&>(*base).getLevelMeter().setLevel(0.3f, true, -18.0f); // lights the LEVEL LEDs incl. OVL
    base->setSize(1747, 764);
    auto img = base->createComponentSnapshot(base->getLocalBounds(), true, 1.0f);
    juce::File(path).deleteFile();
    juce::FileOutputStream out{juce::File(path)};
    juce::PNGImageFormat().writeImageToStream(img, out);
}
