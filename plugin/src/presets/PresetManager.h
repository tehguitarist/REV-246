#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

namespace rv
{
// Presets are the parameter tree saved as XML files (.r246preset) in the user's application-data folder, with a few
// built-in ones (one per program) that can't be overwritten or deleted. Names are file names without the extension.
class PresetManager
{
  public:
    static constexpr const char* extension = ".r246preset";

    explicit PresetManager(juce::AudioProcessorValueTreeState&);

    juce::File folder() const { return presetFolder; }
    juce::StringArray factoryNames() const;
    juce::StringArray userNames() const;
    juce::StringArray allNames() const { auto n = factoryNames(); n.addArray(userNames()); return n; }
    bool isFactory(const juce::String& name) const { return factoryNames().contains(name); }

    bool load(const juce::String& name);
    bool save(const juce::String& name); // false for a built-in name or an unusable one
    bool remove(const juce::String& name);
    void step(int direction); // next / previous in allNames(), wrapping

    juce::String currentName() const { return current; }
    void setCurrentName(const juce::String& n) { current = n; saved = stateFor(n); if (onChange) onChange(); }
    bool isModified() const; // the parameters differ from the loaded preset (as it was loaded or saved)
    std::function<void()> onChange; // the current name or list changed (message thread)

    static juce::String sanitise(const juce::String& name);

  private:
    juce::ValueTree factoryState(int program) const;
    juce::ValueTree stateFor(const juce::String& name) const;
    void apply(const juce::ValueTree&);

    juce::AudioProcessorValueTreeState& apvts;
    juce::File presetFolder;
    juce::String current = "Init S";
    juce::ValueTree saved; // the current preset's parameters, for isModified()
};
} // namespace rv
