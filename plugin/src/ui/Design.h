#pragma once

#include "ui/Layout.h"

// Hand-measured design values that ui-layout.csv doesn't carry (measured on ui/main.png, in the same 3495x1530
// design space as Layout.h) and the behaviour the designer's notes describe in words.
namespace rv::design
{
// The reference editor is 1398 x 612 (40% of the design space); resizable from 25% to 60%.
inline constexpr int defaultWidth = 1398;
inline constexpr int defaultHeight = 612;
inline constexpr float minScale = 0.25f, maxScale = 0.60f;

// The bottom bar (presets, link, help) sits under main.png and is drawn in code.
inline constexpr float barTop = 1380.0f;
inline constexpr float barHeight = 150.0f;
inline const juce::Colour barColour{0xff1a272d}; // main.png's own bottom edge
inline const juce::Colour gold{0xfff8c55e};
inline const juce::Colour labelSelected{0xfffff3d6};
inline const juce::Colour labelDimmed{0xff6a6f66};

inline constexpr float dimmedAlpha = 0.35f;

// Knob and fader drop shadows (the sheet: 315 degrees, 50%). The offset points down and to the right, as in the
// example render; flip the sign here if the final art wants the other way.
struct Shadow
{
    float offset, radius, alpha;
};
inline constexpr Shadow knobShadow{4.7f, 0.4f, 0.5f};
inline constexpr Shadow faderShadow{5.7f, 0.8f, 0.5f};
inline constexpr float shadowDirX = 0.7071f, shadowDirY = 0.7071f;

// The pan knob runs 7 to 5 o'clock: 150 degrees either side of 12.
inline constexpr float panSweepDegrees = 150.0f;

// Fader travel: the cap centre's y at the top and bottom of the track.
inline constexpr float mixTop = 545.0f, mixBottom = 1185.0f;
inline constexpr float gainTop = 218.0f, gainBottom = 1238.0f;
inline constexpr float faderCapW = 100.0f, faderCapH = 211.0f;

// Joystick columns' centre x (matches the funnels above them).
inline constexpr float joyX[] = {897.0f, 1198.0f, 1431.0f, 1722.0f, 2004.0f, 2280.0f};

// The funnel (wedge) above each joystick. The value's red line runs from the narrow bottom (minimum) to the wide
// top (maximum), following the wedge, which is slightly curved. Edges of the wedge's dark interior measured on
// the art for the first funnel; the others are the same shape shifted by funnelDx.
inline constexpr float funnelLineHeight = 10.0f;
inline constexpr float funnelMinY = 1025.0f, funnelMaxY = 532.0f; // line centre
inline constexpr float funnelInset = 3.0f;                         // keeps the line inside the gold edge
struct FunnelPoint
{
    float y, left, right;
};
inline constexpr FunnelPoint funnelEdges[] = {{560, 817, 976}, {600, 825, 968}, {640, 831, 959}, {680, 838, 953},
                                              {720, 844, 947}, {760, 849, 941}, {800, 854, 935}, {840, 858, 930},
                                              {880, 862, 926}, {920, 867, 921}, {960, 872, 917}, {1000, 876, 912}};
inline constexpr float funnelDx[] = {0.0f, 301.0f, 534.0f, 825.0f, 1107.0f, 1383.0f};

// The seven-segment displays: digit area inset in the box, slant and colours.
inline constexpr float displayInsetX = 26.0f, displayInsetY = 24.0f;
inline constexpr float displaySlantDegrees = 4.0f;
inline const juce::Colour segmentLit{0xffff2a1a};
inline const juce::Colour segmentGlow{0xffff1500};
inline const juce::Colour segmentCore{0xffffb59a}; // the hot centre of a lit segment
inline constexpr float segmentGhostAlpha = 0.3f;

// The LEVEL LEDs, top (OVL) to bottom (-12). Thresholds in dBFS of the input after the INPUT fader: the scale reads
// 0 at -18 dBFS (H: the unit's alignment level isn't known) and OVL is a full-scale peak.
inline constexpr int numLevelLeds = 7;
inline constexpr float levelRefDbFs = -18.0f;
inline constexpr float levelStepsDb[numLevelLeds] = {0.0f, 3.0f, 0.0f, -3.0f, -6.0f, -9.0f, -12.0f}; // OVL: special
inline const juce::Colour levelRed{0xffff3a28}, levelYellow{0xffffc21f}, levelGreen{0xff39e04a};
} // namespace rv::design
