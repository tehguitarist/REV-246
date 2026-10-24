#include "params/Parameters.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace rv::params;
using Catch::Matchers::WithinAbs;

TEST_CASE("the gain scale hits the printed ticks")
{
    CHECK_THAT(gainDbFromPosition(defaultGainPosition), WithinAbs(0.0, 0.05));
    CHECK_THAT(gainDbFromPosition(1.0f), WithinAbs(12.0, 1e-4));
    CHECK(gainDbFromPosition(0.0f) <= -99.0f);
    CHECK(gainLinear(0.0f) == 0.0f);
    for (float db : {12.f, 6.f, 0.f, -5.f, -17.f, -50.f})
        CHECK_THAT(gainDbFromPosition(gainPositionFromDb(db)), WithinAbs(db, 1e-3));
}

TEST_CASE("the gain text has two decimals")
{
    CHECK(gainText(gainPositionFromDb(-3.0f)) == "-3.00 dB");
}

TEST_CASE("the S defaults read as the reference setting")
{
    CHECK(joyValue(Program::s, Joy::lowpass, joyDefaultStep(Joy::lowpass)).text == "12");
    CHECK(joyValue(Program::s, Joy::delay, joyDefaultStep(Joy::delay)).text == "10");
    CHECK(joyValue(Program::s, Joy::ampl, joyDefaultStep(Joy::ampl)).text == "80");
    CHECK(joyValue(Program::s, Joy::lf, joyDefaultStep(Joy::lf)).text == "1.0");
    CHECK(joyValue(Program::s, Joy::hf, joyDefaultStep(Joy::hf)).text == "0.9");
}

TEST_CASE("N has no RT")
{
    CHECK_FALSE(joyValue(Program::n, Joy::rt, 3).used);
    CHECK(joyValue(Program::s, Joy::rt, 3).used);
}

TEST_CASE("step ranges")
{
    CHECK(joyMaxStep(Joy::lowpass) == 7);  // 8 LP steps
    CHECK(joyMaxStep(Joy::ampl) == 10);    // 11 AMPL steps
    CHECK(joyMaxStep(Joy::lf) == 15);      // 16 LF steps
    CHECK(joyMaxStep(Joy::hf) == 9);       // 10 HF steps
}
