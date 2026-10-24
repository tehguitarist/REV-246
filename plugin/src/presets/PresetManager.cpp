#include "presets/PresetManager.h"

#include "params/Parameters.h"

namespace rv
{
namespace
{
juce::String factoryName(int program)
{
    return "Init " + juce::String::charToString("SMLEN"[program]);
}
} // namespace

PresetManager::PresetManager(juce::AudioProcessorValueTreeState& s) : apvts(s)
{
    presetFolder = juce::File::getSpecialLocation(juce::File::userApplicationDataDirectory)
#if JUCE_MAC
                       .getChildFile("Audio/Presets")
#endif
                       .getChildFile("Leigh Pierce")
                       .getChildFile("REV-246");
}

juce::StringArray PresetManager::factoryNames() const
{
    juce::StringArray n;
    for (int i = 0; i < params::numPrograms; ++i)
        n.add(factoryName(i));
    return n;
}

juce::StringArray PresetManager::userNames() const
{
    juce::StringArray n;
    for (const auto& f : presetFolder.findChildFiles(juce::File::findFiles, false, juce::String("*") + extension))
        n.add(f.getFileNameWithoutExtension());
    n.sortNatural();
    return n;
}

juce::String PresetManager::sanitise(const juce::String& name)
{
    return juce::File::createLegalFileName(name.trim()).trim();
}

juce::ValueTree PresetManager::factoryState(int program) const
{
    // Every parameter at its default, with the program set; the settings are the unit's own defaults (params).
    auto state = apvts.copyState();
    for (const auto* p : apvts.processor.getParameters())
        if (auto* r = dynamic_cast<const juce::RangedAudioParameter*>(p))
            if (r->paramID != params::id::power)
                if (auto child = state.getChildWithProperty("id", r->paramID); child.isValid())
                    child.setProperty("value", r->convertFrom0to1(r->getDefaultValue()), nullptr);
    if (auto child = state.getChildWithProperty("id", params::id::program); child.isValid())
        child.setProperty("value", program, nullptr);
    return state;
}

juce::ValueTree PresetManager::stateFor(const juce::String& name) const
{
    const auto f = factoryNames().indexOf(name);
    if (f >= 0)
        return factoryState(f);
    const auto file = presetFolder.getChildFile(name + extension);
    if (auto xml = juce::parseXML(file))
        return juce::ValueTree::fromXml(*xml);
    return {};
}

void PresetManager::apply(const juce::ValueTree& state)
{
    // Parameter by parameter, so the host sees ordinary automation-style changes and the power button is left alone.
    for (auto child : state)
    {
        const auto pid = child.getProperty("id").toString();
        if (pid == params::id::power || pid == params::id::link)
            continue;
        if (auto* p = apvts.getParameter(pid))
        {
            p->beginChangeGesture();
            p->setValueNotifyingHost(p->convertTo0to1((float)child.getProperty("value")));
            p->endChangeGesture();
        }
    }
}

bool PresetManager::load(const juce::String& name)
{
    const auto state = stateFor(name);
    if (!state.isValid())
        return false;
    apply(state);
    current = name;
    saved = state;
    if (onChange)
        onChange();
    return true;
}

bool PresetManager::save(const juce::String& rawName)
{
    const auto name = sanitise(rawName);
    if (name.isEmpty() || isFactory(name))
        return false;
    if (!presetFolder.createDirectory().wasOk())
        return false;
    const auto xml = apvts.copyState().createXml();
    if (xml == nullptr || !xml->writeTo(presetFolder.getChildFile(name + extension)))
        return false;
    current = name;
    saved = apvts.copyState();
    if (onChange)
        onChange();
    return true;
}

bool PresetManager::remove(const juce::String& name)
{
    if (isFactory(name))
        return false;
    const auto ok = presetFolder.getChildFile(name + extension).deleteFile();
    if (ok && current == name)
        current = factoryName(0);
    if (onChange)
        onChange();
    return ok;
}

void PresetManager::step(int direction)
{
    const auto names = allNames();
    if (names.isEmpty())
        return;
    const auto i = names.indexOf(current);
    load(names[(i < 0 ? 0 : i + direction + names.size()) % names.size()]);
}

bool PresetManager::isModified() const
{
    if (!saved.isValid())
        return false;
    for (auto child : saved)
    {
        const auto pid = child.getProperty("id").toString();
        if (pid == params::id::power || pid == params::id::link)
            continue;
        if (auto* p = apvts.getParameter(pid))
            if (std::abs(p->convertFrom0to1(p->getValue()) - (float)child.getProperty("value")) > 1e-4f)
                return true;
    }
    return false;
}
} // namespace rv
