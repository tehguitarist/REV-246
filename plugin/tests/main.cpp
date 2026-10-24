#include <catch2/catch_session.hpp>
#include <juce_events/juce_events.h>

int main(int argc, char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui; // the message thread and a display connection, for the editor tests
    return Catch::Session().run(argc, argv);
}
