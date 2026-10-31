#include "PluginProcessor.h"
#include "presets/FactoryPresets.h"

#include <catch2/catch_approx.hpp>
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
    CHECK(p.presets.factoryNames().size() == rv::params::numPrograms + (int)rv::presets::factoryPresets.size());
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
    const auto names = p.presets.factoryNames();
    for (int i = 0; i < names.size(); ++i)
    {
        const auto program = rv::PresetManager::factoryProgram(i);
        REQUIRE(p.presets.load(names[i]));
        INFO(names[i] << ", program " << program);
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

TEST_CASE("the saved state carries its layout version, and an unversioned state still loads")
{
    Rev246Processor a;
    a.apvts.getParameter(rv::params::id::delay)->setValueNotifyingHost(a.apvts.getParameter(rv::params::id::delay)->convertTo0to1(77.0f));
    juce::MemoryBlock block;
    a.getStateInformation(block);
    const auto xml = juce::AudioProcessor::getXmlFromBinary(block.getData(), (int)block.getSize());
    REQUIRE(xml != nullptr);
    CHECK(xml->getIntAttribute(rv::params::stateVersionProperty.toString()) == rv::params::stateVersion);

    // A pre-versioning session (stateVersion absent = 0) loads as the same layout.
    xml->removeAttribute(rv::params::stateVersionProperty.toString());
    juce::MemoryBlock old;
    juce::AudioProcessor::copyXmlToBinary(*xml, old);
    Rev246Processor b;
    b.setStateInformation(old.getData(), (int)old.getSize());
    CHECK(juce::roundToInt(b.apvts.getRawParameterValue(rv::params::id::delay)->load()) == 77);
    CHECK((int)b.apvts.state.getProperty(rv::params::stateVersionProperty) == rv::params::stateVersion);
}

TEST_CASE("an Altiverb factory preset sets its program and its six panel steps")
{
    using namespace rv::params;
    Rev246Processor p;
    for (const auto& f : rv::presets::factoryPresets)
    {
        REQUIRE(p.presets.load(f.name));
        INFO(f.name);
        CHECK(juce::roundToInt(p.apvts.getRawParameterValue(id::program)->load()) == f.program);
        for (int j = 0; j < numJoys; ++j)
            CHECK(juce::roundToInt(p.apvts.getRawParameterValue(joyId((Joy)j))->load()) == f.steps[(size_t)j]);
        CHECK(!p.presets.isModified());
    }
}

TEST_CASE("a session saved before pre-delay and the ducker resets them to their defaults when loaded")
{
    namespace id = rv::params::id;
    Rev246Processor fresh;
    juce::MemoryBlock block;
    fresh.getStateInformation(block);
    const auto xml = juce::AudioProcessor::getXmlFromBinary(block.getData(), (int)block.getSize());
    REQUIRE(xml != nullptr);
    for (const auto* paramId : {id::preDelayMs, id::duckEnabled, id::duckRangeDb, id::duckThresholdDb, id::duckAttackMs,
                                id::duckReleaseMs})
        for (auto* child = xml->getFirstChildElement(); child != nullptr; child = child->getNextElement())
            if (child->getStringAttribute("id") == paramId)
            {
                xml->removeChildElement(child, true);
                break;
            }
    xml->setAttribute(rv::params::stateVersionProperty.toString(), 1);
    juce::MemoryBlock old;
    juce::AudioProcessor::copyXmlToBinary(*xml, old);

    Rev246Processor p;
    p.apvts.getParameter(id::preDelayMs)->setValueNotifyingHost(p.apvts.getParameter(id::preDelayMs)->convertTo0to1(80.0f));
    p.apvts.getParameter(id::duckEnabled)->setValueNotifyingHost(1.0f);
    p.setStateInformation(old.getData(), (int)old.getSize());
    CHECK(juce::exactlyEqual(p.apvts.getRawParameterValue(id::preDelayMs)->load(), 0.0f));
    CHECK(juce::exactlyEqual(p.apvts.getRawParameterValue(id::duckEnabled)->load(), 0.0f));
    CHECK(p.apvts.getRawParameterValue(id::duckRangeDb)->load() == Catch::Approx(-10.0f));
}
