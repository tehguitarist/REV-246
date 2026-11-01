#pragma once

#include "ui/DesignComponent.h"

namespace rv::ui
{
// The seven LEVEL boxes (OVL, +3, 0, -3, -6, -9, -12): each lights its colour with a glow while the input VU level is
// at or above its threshold, and OVL on a full-scale peak. Repaints only when the number of lit boxes changes.
class LevelMeter : public DesignComponent
{
  public:
    LevelMeter();
    void setLevel(float vuLinear, bool clip, float refDbFs); // called from the UI timer: holds and decays the reading
    juce::String getTooltip() override;
    void paint(juce::Graphics&) override;

  private:
    CachedLayer layer;
    // The displayed level holds its highest value for holdMs, then falls at fallDbPerSecond until it meets the VU. A real VU
    // needle has no hold (a peak-programme meter does); the hold is here so a snare's hit stays readable on the LEDs.
    static constexpr double holdMs = 400.0;
    static constexpr float fallDbPerSecond = 20.0f;
    float heldDb = -120.0f;
    double holdUntilMs = 0.0, lastMs = 0.0, overUntilMs = 0.0;
    int lit = 0;      // boxes lit from the bottom: 0 = none .. 7 = all (OVL)
    bool over = false; // OVL is its own condition (a full-scale peak), not "above +3"
};
} // namespace rv::ui
