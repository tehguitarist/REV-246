#include "dsp/Ducker.h"
#include "dsp/PreDelay.h"

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

TEST_CASE("PreDelay is transparent at 0 ms and delays by the set time")
{
    rv::dsp::PreDelay p;
    p.prepare(48000.0);
    std::vector<float> l(4800, 0.0f), r(4800, 0.0f);
    l[0] = r[0] = 1.0f;
    auto copyL = l;
    p.process(l.data(), r.data(), 4800);
    CHECK(l == copyL);

    p.reset();
    p.setDelayMs(10.0f);
    for (int i = 0; i < 40000; ++i) // let the slewed delay settle
    {
        float a = 0.0f, b = 0.0f;
        p.process(&a, &b, 1);
    }
    std::vector<float> l2(2000, 0.0f), r2(2000, 0.0f);
    l2[0] = r2[0] = 1.0f;
    p.process(l2.data(), r2.data(), 2000);
    CHECK(l2[480] == Catch::Approx(1.0f).margin(1e-3));
}
