#include "dsp/Ducker.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

TEST_CASE("Ducker lowers the wet by the range while the key is over the threshold, and not otherwise")
{
    rv::dsp::Ducker d;
    d.prepare(48000.0);
    d.setEnabled(true);
    d.setRangeDb(-10.0f);
    d.setThresholdDb(-20.0f);
    d.setAttackMs(0.0f);
    d.setReleaseMs(0.0f);

    std::vector<float> key(100, 0.5f), l(100, 1.0f), r(100, 1.0f);
    d.process(key.data(), l.data(), r.data(), 100);
    CHECK(l[50] == Catch::Approx(std::pow(10.0f, -0.5f)).margin(1e-4));

    std::fill(key.begin(), key.end(), 0.01f);
    std::fill(l.begin(), l.end(), 1.0f);
    std::fill(r.begin(), r.end(), 1.0f);
    d.process(key.data(), l.data(), r.data(), 100);
    CHECK(l[50] == Catch::Approx(1.0f).margin(1e-4));

    d.setEnabled(false);
    std::fill(key.begin(), key.end(), 0.5f);
    std::fill(l.begin(), l.end(), 1.0f);
    d.process(key.data(), l.data(), r.data(), 100);
    CHECK(l[99] == Catch::Approx(1.0f).margin(1e-4));
}
