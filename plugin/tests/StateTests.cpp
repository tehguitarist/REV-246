#include "PluginProcessor.h"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("state round-trips")
{
    Rev246Processor a;
    a.apvts.getParameter(rv::params::id::program)->setValueNotifyingHost(a.apvts.getParameter(rv::params::id::program)->convertTo0to1(2.0f));
    a.apvts.getParameter(rv::params::id::delay)->setValueNotifyingHost(a.apvts.getParameter(rv::params::id::delay)->convertTo0to1(123.0f));
    juce::MemoryBlock block;
    a.getStateInformation(block);

    Rev246Processor b;
    b.setStateInformation(block.getData(), (int)block.getSize());
    CHECK(juce::roundToInt(b.apvts.getRawParameterValue(rv::params::id::program)->load()) == 2);
    CHECK(juce::roundToInt(b.apvts.getRawParameterValue(rv::params::id::delay)->load()) == 123);
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
