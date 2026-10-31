#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

#include <memory>

namespace rv::ui
{
// The DUCK pop-up: enable/disable and the range, threshold, attack and release fields. The fields take a vertical drag
// or a typed value (two decimal places); a negative range lowers the wet signal while the input is over the threshold.
std::unique_ptr<juce::Component> createDuckPanel(juce::AudioProcessorValueTreeState&);
} // namespace rv::ui
