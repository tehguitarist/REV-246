#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

#include <functional>

namespace rv
{
// Presets are the parameter tree saved as XML files (.r246preset) in the user's application-data folder, with a few
// built-in ones (one per program) that can't be overwritten or deleted. Names are file names without the extension.
//
// Threads: everything is for the message thread, except that currentName(), setCurrentName() and isModified() may
// also be reached from a host's state calls on another thread, so the name and the reference state are locked.
class PresetManager
{
  public:
    static constexpr const char* extension = ".r246preset";

    // The built-in preset for a program ("Init S", "Init M", ...), also the initial current preset.
    static juce::String factoryName(int program);
    // A user-typed name as a file name: legal characters only, trimmed. Empty if nothing usable is left.
    static juce::String sanitise(const juce::String& name);

    explicit PresetManager(juce::AudioProcessorValueTreeState&);

    juce::File folder() const { return presetFolder; }
    // The built-in presets: the five "Init" ones, then the Altiverb EMT 246 presets as best-guess panel settings
    // (FactoryPresets.h, docs/library_notes.md section 7).
    juce::StringArray factoryNames() const;
    // The program of factoryNames()[index] (for grouping them in the menu).
    static int factoryProgram(int index);
    juce::StringArray userNames() const;
    juce::StringArray allNames() const { auto n = factoryNames(); n.addArray(userNames()); return n; }
    bool isFactory(const juce::String& name) const { return factoryNames().contains(name, true); }

    bool load(const juce::String& name); // false if there is no such preset or its file is unusable
    bool save(const juce::String& name); // false for a built-in name, an unusable one or a write error
    bool remove(const juce::String& name);
    void step(int direction); // next (> 0) / previous (< 0) in allNames(), wrapping; skips unusable files

    juce::String currentName() const;
    // Marks `name` as the current preset without applying it (a restored session has its own parameter values).
    void setCurrentName(const juce::String& name);
    bool isModified() const; // the parameters differ from the loaded preset (as it was loaded or saved)
    std::function<void()> onChange; // the current name or list changed (message thread)

  private:
    juce::ValueTree factoryState(int index) const; // index into factoryNames()
    juce::ValueTree stateFor(const juce::String& name) const;
    juce::File presetFile(const juce::String& name) const;
    void apply(const juce::ValueTree&);
    void setCurrent(const juce::String& name, const juce::ValueTree& state);

    juce::AudioProcessorValueTreeState& apvts;
    juce::File presetFolder;
    mutable juce::CriticalSection stateLock; // guards `current` and `saved`
    juce::String current = factoryName(0);
    juce::ValueTree saved; // the current preset's parameters, for isModified()
};
} // namespace rv
