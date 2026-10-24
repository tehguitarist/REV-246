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
    void setLevel(float vuLinear, bool clip, float refDbFs);
    int litCount() const { return lit; }
    juce::String getTooltip() override;
    void paint(juce::Graphics&) override;

  private:
    CachedLayer layer;
    int lit = 0;      // boxes lit from the bottom: 0 = none .. 7 = all (OVL)
    bool over = false; // OVL is its own condition (a full-scale peak), not "above +3"
};
} // namespace rv::ui
