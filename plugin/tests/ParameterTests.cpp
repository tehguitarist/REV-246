#include "params/Parameters.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <functional>
#include <vector>

using namespace rv::params;
using Catch::Matchers::WithinAbs;

TEST_CASE("the gain scale hits the printed ticks")
{
    CHECK_THAT(gainDbFromPosition(defaultGainPosition), WithinAbs(0.0, 0.05));
    CHECK_THAT(gainDbFromPosition(1.0f), WithinAbs(12.0, 1e-4));
    CHECK(gainDbFromPosition(0.0f) <= -99.0f);
    CHECK(std::equal_to<float>{}(gainLinearFromPosition(0.0f), 0.0f));
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
    CHECK(joyMaxStep(Joy::rt) == 15);      // 16 RT steps
    CHECK(joyMaxStep(Joy::hf) == 10);      // the parameter: N's 11 LATE values
    for (const auto p : {Program::s, Program::m, Program::l, Program::e})
        CHECK(joyMaxStep(p, Joy::hf) == 9); // 10 HF steps
    CHECK(joyMaxStep(Program::n, Joy::hf) == 10);
    CHECK(joyValue(Program::s, Joy::hf, 10).text == "1.0"); // past S's last HF step reads as it
}

TEST_CASE("the RT tables hold the library's steps and the data sheet's ranges")
{
    struct Case
    {
        Program p;
        double lo, hi;
        std::vector<double> library;
    };
    const std::vector<Case> cases{
        {Program::s, 0.2, 2.5, {0.3, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.5}},
        {Program::m, 0.5, 5.0, {0.6, 0.8, 1.2, 1.8, 2.6, 3.0, 4.0, 5.0}},
        {Program::l, 2.0, 20.0, {2.2, 2.6, 3.6, 4.5, 6.0, 8.0, 10.0, 15.0, 20.0}},
        {Program::e, 0.4, 5.0, {0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6, 4.0, 4.5}},
    };
    for (const auto& c : cases)
    {
        INFO("program " << (int)c.p);
        CHECK(std::equal_to<double>{}(joyNumber(c.p, Joy::rt, 0), c.lo));
        CHECK(std::equal_to<double>{}(joyNumber(c.p, Joy::rt, joyMaxStep(Joy::rt)), c.hi));
        for (int step = 1; step <= joyMaxStep(Joy::rt); ++step)
            CHECK(joyNumber(c.p, Joy::rt, step) > joyNumber(c.p, Joy::rt, step - 1));
        for (const auto v : c.library)
        {
            const auto step = joyStepFromValue(c.p, Joy::rt, v);
            REQUIRE(step.has_value());
            CHECK(std::equal_to<double>{}(joyNumber(c.p, Joy::rt, *step), v)); // an exact step, not a neighbour
        }
    }
    CHECK(std::equal_to<double>{}(joyNumber(Program::s, Joy::rt, joyDefaultStep(Joy::rt)), 1.0));
}

TEST_CASE("the LF and HF tables hold every library value of S, M, L and E")
{
    for (const auto v : {0.1, 0.4, 1.0, 2.0})
        CHECK(joyStepFromValue(Program::s, Joy::lf, v).has_value());
    for (const auto v : {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0})
        CHECK(joyStepFromValue(Program::m, Joy::hf, v).has_value());
}

TEST_CASE("N's columns are its spec's values")
{
    const std::vector<double> early{-3.0, -2.0, -1.0, -0.7, -0.5, -0.2, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 10.0};
    for (int step = 0; step <= joyMaxStep(Joy::lf); ++step)
        CHECK(std::equal_to<double>{}(joyNumber(Program::n, Joy::lf, step), early[(size_t)step]));
    for (int step = 0; step <= joyMaxStep(Program::n, Joy::hf); ++step)
        CHECK(std::abs(joyNumber(Program::n, Joy::hf, step) - 0.2 * step) < 1e-12);
    for (int step = 0; step <= joyMaxStep(Joy::ampl); ++step)
        CHECK(std::abs(joyNumber(Program::n, Joy::ampl, step) - 0.05 * step) < 1e-12);
    CHECK(joyValue(Program::n, Joy::lf, 3).text == "-0.7");
    CHECK(joyValue(Program::n, Joy::lf, 15).text == "10");
    CHECK(joyValue(Program::n, Joy::hf, 10).text == "2.0");
}

TEST_CASE("typed values map to exact steps where the control has a visible table")
{
    auto low = joyStepFromValue(Program::s, Joy::lowpass, 8.0);
    REQUIRE(low.has_value());
    CHECK(*low == 5);
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::lowpass, 7.0).has_value());

    auto delay = joyStepFromValue(Program::s, Joy::delay, 123.0);
    REQUIRE(delay.has_value());
    CHECK(*delay == 123);
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::delay, 12.5).has_value());
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::delay, 1000.0).has_value());

    auto ampl = joyStepFromValue(Program::s, Joy::ampl, 80.0);
    REQUIRE(ampl.has_value());
    CHECK(*ampl == 8);
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::ampl, 75.0).has_value());

    auto lf = joyStepFromValue(Program::s, Joy::lf, 1.4);
    REQUIRE(lf.has_value());
    CHECK(*lf == 12);
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::lf, 1.3).has_value());

    auto hf = joyStepFromValue(Program::s, Joy::hf, 0.9);
    REQUIRE(hf.has_value());
    CHECK(*hf == 8);
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::hf, 0.95).has_value());
}

TEST_CASE("the non-linear program uses its own typed-value ranges")
{
    auto span = joyStepFromValue(Program::n, Joy::ampl, 0.25);
    REQUIRE(span.has_value());
    CHECK(*span == 5);
    CHECK_FALSE(joyStepFromValue(Program::n, Joy::ampl, 0.26).has_value());

    auto early = joyStepFromValue(Program::n, Joy::lf, -3.0);
    REQUIRE(early.has_value());
    CHECK(*early == 0);
    CHECK_FALSE(joyStepFromValue(Program::n, Joy::lf, -3.25).has_value());
    CHECK_FALSE(joyStepFromValue(Program::n, Joy::lf, -4.0).has_value());

    CHECK_FALSE(joyStepFromValue(Program::n, Joy::rt, 1.0).has_value()); // RT is not used

    auto late = joyStepFromValue(Program::n, Joy::hf, 1.2);
    REQUIRE(late.has_value());
    CHECK(*late == 6);
    auto lateMax = joyStepFromValue(Program::n, Joy::hf, 2.0);
    REQUIRE(lateMax.has_value());
    CHECK(*lateMax == 10);
    CHECK_FALSE(joyStepFromValue(Program::n, Joy::hf, 1.3).has_value());
    CHECK_FALSE(joyStepFromValue(Program::n, Joy::hf, 2.2).has_value());
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::hf, 1.1).has_value());
}

TEST_CASE("the geometric RT accepts a range and picks the nearest step")
{
    auto rt = joyStepFromValue(Program::s, Joy::rt, 1.0);
    REQUIRE(rt.has_value());
    CHECK(*rt >= 0);
    CHECK(*rt <= joyMaxStep(Joy::rt));
    CHECK(std::abs(joyValue(Program::s, Joy::rt, *rt).value - 1.0) < 0.2);
    CHECK_FALSE(joyStepFromValue(Program::s, Joy::rt, 3.0).has_value());
}

TEST_CASE("the typed-entry hints describe the S/M/L/E and N differences")
{
    CHECK(joyInputHint(Program::s, Joy::ampl).contains("10 %"));
    CHECK(joyInputHint(Program::n, Joy::ampl).contains("0.05"));
    CHECK(joyInputHint(Program::s, Joy::lf).contains("1.4"));
    CHECK(joyInputHint(Program::n, Joy::lf).contains("-3.0"));
    CHECK(joyInputHint(Program::n, Joy::lf).contains("10 = flat"));
    CHECK(joyInputHint(Program::s, Joy::rt).contains("nearest"));
    CHECK(joyInputHint(Program::n, Joy::rt).contains("not used"));
    CHECK(joyInputHint(Program::s, Joy::hf).contains("0.1"));
    CHECK(joyInputHint(Program::n, Joy::hf).contains("0.2 s"));
    CHECK(joyInputName(Program::s, Joy::ampl).contains("amplitude"));
    CHECK(joyInputName(Program::n, Joy::ampl).contains("span"));
}
