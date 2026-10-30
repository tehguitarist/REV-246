#include "dsp/RefIr.h"

#include <catch2/catch_test_macros.hpp>

namespace
{
// A core that writes its input to both outputs, delayed by `d` samples: checks the render and null helpers.
class DelayCore final : public rv::dsp::ProgramCore
{
  public:
    explicit DelayCore(int d) : buf((size_t)d + 1, 0.0f), delay(d) {}
    void prepare() override {}
    void reset() override { std::fill(buf.begin(), buf.end(), 0.0f); }
    void setSettings(const rv::dsp::Settings&) override {}
    void process(const float* in, float* outL, float* outR, int n) override
    {
        for (int i = 0; i < n; ++i)
        {
            buf[pos] = in[i];
            pos = (pos + 1) % buf.size();
            outL[i] = outR[i] = buf[pos];
        }
    }
    int latency() const override { return delay; }

  private:
    std::vector<float> buf;
    size_t pos = 0;
    int delay;
};
} // namespace

TEST_CASE("Harness: a reference IR loads at the device rate", "[harness]")
{
    const auto ir = rv::test::loadWav("analysis/s/ref/s_1.0_1.0_1.0.wav");
    REQUIRE(ir.size() > 10000);
    double e = 0.0;
    for (auto v : ir.l)
        e += v * v;
    CHECK(e > 0.0);
}

TEST_CASE("Harness: render removes the latency and the null of an exact match is at the floor", "[harness]")
{
    DelayCore core(40);
    core.prepare();
    core.reset();
    const auto got = rv::test::renderImpulse(core, 2000);
    REQUIRE(std::abs(got.l[0] - 1.0) < 1e-12);
    rv::test::StereoIr ref{std::vector<double>(2000, 0.0), std::vector<double>(2000, 0.0)};
    ref.l[0] = ref.r[0] = 1.0;
    CHECK(rv::test::nullDb(ref.l, got.l, 0, 50) < -300.0);
    ref.l[0] = 1.001;
    CHECK(std::abs(rv::test::nullDb(ref.l, got.l, 0, 50) + 60.0) < 0.1);
}
