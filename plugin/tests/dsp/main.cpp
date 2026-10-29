#include <catch2/catch_session.hpp>

// The DSP test runner: no GUI, no processor, only the engine and its cores (docs/ENGINE.md, Tests).
int main(int argc, char** argv)
{
    return Catch::Session().run(argc, argv);
}
