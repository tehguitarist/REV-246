#include "PluginProcessor.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>

TEST_CASE("state round-trips")
{
    Rev246Processor a;
    a.apvts.getParameter(rv::params::id::program)->setValueNotifyingHost(a.apvts.getParameter(rv::params::id::program)->convertTo0to1(2.0f));
    a.apvts.getParameter(rv::params::id::delay)->setValueNotifyingHost(a.apvts.getParameter(rv::params::id::delay)->convertTo0to1(123.0f));
    a.apvts.getParameter(rv::params::id::sixteenBit)->setValueNotifyingHost(0.0f);
    a.apvts.getParameter(rv::params::id::stereoWidth)->setValueNotifyingHost(
        a.apvts.getParameter(rv::params::id::stereoWidth)->convertTo0to1(150.0f));
    juce::MemoryBlock block;
    a.getStateInformation(block);

    Rev246Processor b;
    b.setStateInformation(block.getData(), (int)block.getSize());
    CHECK(juce::roundToInt(b.apvts.getRawParameterValue(rv::params::id::program)->load()) == 2);
    CHECK(juce::roundToInt(b.apvts.getRawParameterValue(rv::params::id::delay)->load()) == 123);
    CHECK(b.apvts.getRawParameterValue(rv::params::id::sixteenBit)->load() < 0.5f);
    CHECK(std::abs(b.apvts.getRawParameterValue(rv::params::id::stereoWidth)->load() - 150.0f) < 1e-6f);
}

TEST_CASE("factory presets load and are not user-overwritable")
{
    Rev246Processor p;
    CHECK(p.presets.factoryNames().size() == rv::params::numPrograms);
    CHECK(p.presets.load("Init N"));
    CHECK(juce::roundToInt(p.apvts.getRawParameterValue(rv::params::id::program)->load()) == 4);
    CHECK_FALSE(p.presets.save("Init N"));
    CHECK_FALSE(p.presets.isModified());
    p.apvts.getParameter(rv::params::id::mix)->setValueNotifyingHost(0.9f);
    CHECK(p.presets.isModified());
}

TEST_CASE("every factory preset sits on its program's step tables; E has the Altiverb family's LF / HF")
{
    using namespace rv::params;
    Rev246Processor p;
    for (int program = 0; program < numPrograms; ++program)
    {
        REQUIRE(p.presets.load(rv::PresetManager::factoryName(program)));
        INFO("program " << program);
        CHECK(juce::roundToInt(p.apvts.getRawParameterValue(id::program)->load()) == program);
        for (int j = 0; j < numJoys; ++j)
        {
            const auto raw = p.apvts.getRawParameterValue(joyId((Joy)j))->load();
            CHECK(std::equal_to<float>{}(raw, std::round(raw)));
            CHECK((int)raw >= 0);
            CHECK((int)raw <= joyMaxStep((Program)program, (Joy)j));
        }
    }
    REQUIRE(p.presets.load("Init E"));
    auto value = [&](Joy j) { return joyNumber(Program::e, j, juce::roundToInt(p.apvts.getRawParameterValue(joyId(j))->load())); };
    CHECK(std::equal_to<double>{}(value(Joy::lf), 1.5)); // the step nearest the family's 1.49
    CHECK(std::abs(value(Joy::hf) - 0.7) < 1e-12); // 0.69
    REQUIRE(p.presets.load("Init S"));
    CHECK(std::equal_to<double>{}(joyNumber(Program::s, Joy::rt, juce::roundToInt(p.apvts.getRawParameterValue(id::rt)->load())), 1.0));
}

TEST_CASE("presets leave the stereo mode alone")
{
    Rev246Processor p;
    auto* mode = p.apvts.getParameter(rv::params::id::trueStereo);
    CHECK(p.presets.load("Init S"));
    mode->setValueNotifyingHost(1.0f);
    CHECK_FALSE(p.presets.isModified()); // a routing preference, not a setting
    CHECK(p.presets.load("Init M"));
    CHECK(mode->getValue() >= 0.5f);
}

TEST_CASE("the stereo width is part of a preset")
{
    Rev246Processor p;
    CHECK(p.presets.load("Init S"));
    CHECK_FALSE(p.presets.isModified());
    auto* width = p.apvts.getParameter(rv::params::id::stereoWidth);
    width->setValueNotifyingHost(width->convertTo0to1(150.0f));
    CHECK(p.presets.isModified()); // a sound setting, unlike the transformer
}

TEST_CASE("the 16-bit mode is part of a preset")
{
    Rev246Processor p;
    CHECK(p.presets.load("Init S"));
    CHECK_FALSE(p.presets.isModified());
    p.apvts.getParameter(rv::params::id::sixteenBit)->setValueNotifyingHost(0.0f);
    CHECK(p.presets.isModified()); // an accuracy setting, unlike the transformer
}
