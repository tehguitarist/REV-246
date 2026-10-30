#pragma once

#include "ui/Layout.h"

// Hand-measured design values that ui-layout.csv doesn't carry (measured on ui/main.png, in the same 3495x1530
// design space as Layout.h) and the behaviour the designer's notes describe in words.
namespace rv::design
{
// The editor is resizable from 25% to 60% of the design space (the default, 40%, is 1398 x 612).
inline constexpr float minScale = 0.25f, maxScale = 0.60f;

// The bottom bar (presets, link, help) sits under main.png and is drawn in code. The version number keeps
// barEdgePad from the panel's right edge, and the "?" keeps the same distance from the version number.
inline constexpr float barEdgePad = 57.0f;
inline constexpr float barTop = 1380.0f;
// The bar's buttons share one row; the width halo is taller and sits a little higher.
inline constexpr float barRowTop = barTop + 35.0f, barRowHeight = 80.0f;
inline constexpr float barGap = 12.0f; // between neighbouring buttons
inline constexpr float linkWidth = 220.0f, transformerWidth = 180.0f, formatWidth = 210.0f, stereoModeWidth = 300.0f;
inline constexpr float helpWidth = 90.0f, calWidth = 220.0f;
inline constexpr float calLeft = 1224.0f;                          // CAL follows the preset buttons at the bar's left end
inline constexpr float widthKnobTop = barTop + 10.0f, widthKnobW = 140.0f, widthKnobH = 130.0f;
inline constexpr float titleLeft = 1597.0f, titleWidth = 300.0f;   // REV-246 centred on the panel
inline constexpr float barRuleThickness = 4.0f, barRuleAlpha = 0.55f; // the gold rule where main.png ends
inline constexpr float titleCapHeight = 42.0f, versionCapHeight = 30.0f;
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

// The pan knob's cap art is darkened and given more contrast when it is loaded (the asset itself is untouched).
inline constexpr float knobBrightness = 0.86f; // multiplies the colour
inline constexpr float knobContrast = 1.10f;   // about mid-grey

// The pan knob runs 7 to 5 o'clock: 150 degrees either side of 12.
inline constexpr float panSweepDegrees = 150.0f;

// Fader travel: the cap centre's y at the top and bottom of the track.
inline constexpr float mixTop = 545.0f, mixBottom = 1185.0f;
inline constexpr float gainTop = 218.0f, gainBottom = 1238.0f;
inline constexpr float faderCapW = 100.0f, faderCapH = 211.0f;

// The funnel (wedge) above each joystick. The value's red line runs from the narrow bottom (minimum) to the wide
// top (maximum), following the wedge, which is slightly curved. Edges of the wedge's dark interior measured on
// the art for the first funnel; the others are the same shape shifted by funnelDx.
// The indicator is a bar of funnelSegments equal segments (the real unit has 8) over the wedge's interior, which runs
// from funnelTopY to funnelBottomY (measured on the art); the segment nearest the value is lit, with faint lines
// between all of them. The segments, their lines and the glow are clipped to the wedge.
inline constexpr int funnelSegments = 10;
inline constexpr float funnelTopY = 526.0f, funnelBottomY = 1029.0f;
inline constexpr float funnelLineThickness = 3.5f;   // the faint dividing lines
inline constexpr float funnelLineAlpha = 0.45f;
// The wedge itself is tinted red around the lit segment, as if lit by it: peak alpha, and how far it reaches (as a
// fraction of the wedge's height).
inline constexpr float funnelAmbient = 0.30f;
inline constexpr float funnelAmbientReach = 0.24f;
inline constexpr float funnelGap = 3.0f;   // dark gap kept each side of a lit segment's lines
inline constexpr float funnelInset = 3.0f; // keeps the band inside the gold edge
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
inline const juce::Colour segmentLit{0xffff8a3a};
inline const juce::Colour segmentGlow{0xffff0a00};
inline constexpr float segmentGhostAlpha = 0.2f;

// The LEVEL LEDs, top (OVL) to bottom (-12). Thresholds in dB relative to the meter's 0 VU reference, of the input
// after the INPUT fader, read on a VU ballistic (not peak). OVL is +6 VU or a peak the A/D clips (or full scale),
// whichever comes first (the engine's DeviceBridge).
// What 0 VU reads in dBFS is the meter calibration, selectable in the bottom bar (default -18 dBFS: 0 = 1 kHz at
// 0.775 V, H).
inline constexpr int numLevelLeds = 7;
inline constexpr float levelStepsDb[numLevelLeds] = {6.0f, 3.0f, 0.0f, -3.0f, -6.0f, -9.0f, -12.0f}; // OVL: +6
inline constexpr int numMeterCal = 6;
inline constexpr float meterCalDbFs[numMeterCal] = {-24.0f, -20.0f, -18.0f, -16.0f, -14.0f, -12.0f};
inline constexpr int defaultMeterCal = 2; // -18 dBFS
inline const juce::Colour levelRed{0xffff3a28}, levelYellow{0xffffc21f}, levelGreen{0xff39e04a};
} // namespace rv::design
