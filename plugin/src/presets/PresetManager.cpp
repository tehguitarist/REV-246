#include "presets/PresetManager.h"

#include "dsp/generated/e_spec.h"
#include "params/Parameters.h"

#include <cmath>

namespace rv
{
namespace
{
// The properties of a parameter's node in the APVTS state; they are what the preset files hold.
const juce::Identifier idProperty{"id"};
const juce::Identifier valueProperty{"value"};

// The preferences (routing and monitoring choices), as opposed to the sound settings: presets neither store nor
// change them (DEVELOPMENT.md, Parameters).
bool isPreference(const juce::String& paramId)
{
    return paramId == params::id::power || paramId == params::id::link || paramId == params::id::transformer ||
           paramId == params::id::trueStereo;
}

// A preset holds parameters only: the session's own nodes (the window size, ...) have no id. Edits `copy` in place
// (a ValueTree is a shared reference), so pass a copy of the live state.
juce::ValueTree parametersOnly(juce::ValueTree copy)
{
    for (int i = copy.getNumChildren(); --i >= 0;)
        if (!copy.getChild(i).hasProperty(idProperty))
            copy.removeChild(i, nullptr);
    return copy;
}

static_assert(sizeof("SMLEN") - 1 == (size_t)params::numPrograms, "one letter per program, in Program order");
} // namespace

juce::String PresetManager::factoryName(int program)
{
    return "Init " + juce::String::charToString("SMLEN"[program]);
}

juce::String PresetManager::sanitise(const juce::String& name)
{
    return juce::File::createLegalFileName(name.trim()).trim();
}

PresetManager::PresetManager(juce::AudioProcessorValueTreeState& s) : apvts(s)
{
    presetFolder = juce::File::getSpecialLocation(juce::File::userApplicationDataDirectory)
#if JUCE_MAC
                       .getChildFile("Audio/Presets")
#endif
                       .getChildFile("Leigh Pierce")
                       .getChildFile("REV-246");

    // A fresh instance shows "Init S"; treat the default state as that preset's saved state, so the first edit marks
    // the name as modified even before a preset has been explicitly loaded.
    saved = stateFor(current);
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

juce::String PresetManager::currentName() const
{
    const juce::ScopedLock lock(stateLock);
    return current;
}

void PresetManager::setCurrent(const juce::String& name, const juce::ValueTree& state)
{
    {
        const juce::ScopedLock lock(stateLock);
        current = name;
        saved = state;
    }
    if (onChange)
        onChange();
}

void PresetManager::setCurrentName(const juce::String& name)
{
    setCurrent(name, stateFor(name));
}

juce::File PresetManager::presetFile(const juce::String& name) const
{
    // A name is a file name, never a path (it can come from a restored session): nothing outside the folder.
    if (name.isEmpty() || name.containsAnyOf("/\\:"))
        return {};
    return presetFolder.getChildFile(name + extension);
}

juce::ValueTree PresetManager::factoryState(int program) const
{
    // Every sound setting at its default, with the program set; the settings are the unit's own defaults (params).
    auto state = apvts.copyState();
    for (const auto* p : apvts.processor.getParameters())
        if (auto* r = dynamic_cast<const juce::RangedAudioParameter*>(p))
            if (!isPreference(r->paramID))
                if (auto child = state.getChildWithProperty(idProperty, r->paramID); child.isValid())
                    child.setProperty(valueProperty, r->convertFrom0to1(r->getDefaultValue()), nullptr);
    auto set = [&state](const char* paramId, int value) {
        if (auto child = state.getChildWithProperty(idProperty, paramId); child.isValid())
            child.setProperty(valueProperty, value, nullptr);
    };
    set(params::id::program, program);
    // E: the Altiverb family's LF x 1.49 / HF x 0.69 (e_plugin_spec controls.LF_HF), at the nearest panel steps (1.5,
    // 0.7). Every other value is a step index, so every preset sits on its program's table.
    if (program == (int)params::Program::e)
    {
        set(params::id::lf, params::joyNearestStep(params::Program::e, params::Joy::lf, dsp::gen::e::familyLf));
        set(params::id::hf, params::joyNearestStep(params::Program::e, params::Joy::hf, dsp::gen::e::familyHf));
    }
    return parametersOnly(state);
}

juce::ValueTree PresetManager::stateFor(const juce::String& name) const
{
    const auto f = factoryNames().indexOf(name);
    if (f >= 0)
        return factoryState(f);
    const auto file = presetFile(name);
    if (file == juce::File())
        return {};
    // Only an XML document of the state's own type is a preset: anything else is a stray or damaged file.
    if (auto xml = juce::parseXML(file); xml != nullptr && xml->hasTagName(apvts.state.getType()))
        return juce::ValueTree::fromXml(*xml);
    return {};
}

void PresetManager::apply(const juce::ValueTree& state)
{
    // Parameter by parameter, so the host sees ordinary automation-style changes and the preferences are left alone.
    for (auto child : state)
    {
        const auto pid = child.getProperty(idProperty).toString();
        if (isPreference(pid) || !child.hasProperty(valueProperty))
            continue;
        if (auto* p = apvts.getParameter(pid))
        {
            p->beginChangeGesture();
            p->setValueNotifyingHost(p->convertTo0to1((float)child.getProperty(valueProperty)));
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
    setCurrent(name, state);
    return true;
}

bool PresetManager::save(const juce::String& rawName)
{
    const auto name = sanitise(rawName);
    const auto file = presetFile(name);
    if (file == juce::File() || isFactory(name))
        return false;
    if (!presetFolder.createDirectory().wasOk())
        return false;
    const auto state = parametersOnly(apvts.copyState());
    const auto xml = state.createXml();
    if (xml == nullptr || !xml->writeTo(file))
        return false;
    setCurrent(name, state);
    return true;
}

bool PresetManager::remove(const juce::String& name)
{
    if (isFactory(name))
        return false;
    const auto file = presetFile(name);
    if (!file.existsAsFile() || !file.deleteFile())
        return false;
    if (currentName() == name) // the parameters stay as they are, now under the default name
        setCurrentName(factoryName(0));
    else if (onChange)
        onChange();
    return true;
}

void PresetManager::step(int direction)
{
    const auto names = allNames();
    const auto count = names.size();
    if (count == 0 || direction == 0)
        return;
    const auto sign = direction > 0 ? 1 : -1;
    const auto start = names.indexOf(currentName());
    // Try each preset in turn from the neighbour, so a damaged file doesn't trap the buttons on it.
    for (int offset = 1; offset <= count; ++offset)
    {
        const auto index = start < 0 ? (sign > 0 ? offset - 1 : count - offset)
                                     : ((start + sign * offset) % count + count) % count;
        if (load(names[index]))
            return;
    }
}

bool PresetManager::isModified() const
{
    juce::ValueTree reference;
    {
        const juce::ScopedLock lock(stateLock);
        reference = saved;
    }
    if (!reference.isValid())
        return false;
    for (auto child : reference)
    {
        const auto pid = child.getProperty(idProperty).toString();
        if (isPreference(pid) || !child.hasProperty(valueProperty))
            continue;
        if (auto* p = apvts.getParameter(pid))
            if (std::abs(p->convertFrom0to1(p->getValue()) - (float)child.getProperty(valueProperty)) > 1e-4f)
                return true;
    }
    return false;
}
} // namespace rv
