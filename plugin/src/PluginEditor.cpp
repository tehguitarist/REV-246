#include "PluginEditor.h"

#include "ui/Design.h"

namespace
{
namespace L = rv::layout;
namespace P = rv::params;
using rv::ui::SegmentDisplay;

constexpr L::Slot joySlots[] = {L::lowpassJoy, L::delayJoy, L::amplJoy, L::lfJoy, L::rtJoy, L::hfJoy};
constexpr int joyCount = P::numJoys;

const char* joyHelp(P::Program p, P::Joy j)
{
    const bool n = p == P::Program::n;
    switch (j)
    {
    case P::Joy::lowpass: return "Lowpass: the output low-pass filter, 2 to 12 kHz. Press the lever up or down (hold to repeat), "
                                 "or use the wheel or drag on the display.";
    case P::Joy::delay: return "Reflections delay: 0 to 999 ms before the first reflections. Hold for more than a second to speed up.";
    case P::Joy::ampl: return n ? "Span: how long the early section lasts, 0 to 0.5 s." : "Reflections amplitude: the early reflections' level, 0 to 100 %.";
    case P::Joy::lf: return n ? "Early RT: the early section's decay slope (negative rises)." : "Decay LF: the multiplier on the decay time at low frequencies.";
    case P::Joy::rt: return n ? "RT is not used in the non-linear program." : "Decay RT: the reverberation time, in seconds (at 1 kHz).";
    case P::Joy::hf: return n ? "Late RT: the late section's decay time." : "Decay HF: the multiplier on the decay time at high frequencies.";
    case P::Joy::count:
    default: return "";
    }
}

juce::Rectangle<float> displayBox(const L::Slot& s) { return s.bounds(); }
} // namespace

Rev246Editor::PanelLookAndFeel::PanelLookAndFeel()
{
    setColour(juce::TooltipWindow::backgroundColourId, juce::Colour(0xf2161615));
    setColour(juce::TooltipWindow::textColourId, juce::Colour(0xffe6e6e2));
    setColour(juce::TooltipWindow::outlineColourId, juce::Colour(0xff4a4a47));
}

juce::RangedAudioParameter& Rev246Editor::param(const char* id) const
{
    return *proc.apvts.getParameter(id);
}

P::Program Rev246Editor::program() const
{
    const auto& p = param(P::id::program);
    return (P::Program)juce::jlimit(0, P::numPrograms - 1, juce::roundToInt(p.convertFrom0to1(p.getValue())));
}

Rev246Editor::Rev246Editor(Rev246Processor& p)
    : juce::AudioProcessorEditor(p), proc(p),
      powerButton(images),
      programDisplay(displayBox(L::programDisplay), 5, {{0, 5, true, {}, {}}}),
      lowpassDisplay(displayBox(L::lowpassDisplay), 3, {{0, 3, false, {}, {}}}),
      reflectionsDisplay(displayBox(L::reflectionsDisplay), 6, {{0, 3, false, {}, {}}, {3, 3, false, {}, {}}}),
      decayDisplay(displayBox(L::decayDisplay), 9, {{0, 3, false, {}, {}}, {3, 3, false, {}, {}}, {6, 3, false, {}, {}}}),
      panKnob(images, *p.apvts.getParameter(P::id::pan), bubble,
              [](float v) { return std::abs(v) < 0.005f ? juce::String("Centre") : juce::String(v < 0.f ? "L " : "R ") + juce::String(std::abs(v) * 100.f, 2) + " %"; }),
      mixFader(images, L::Image::redFader, L::mixFader.cx, rv::design::mixTop, rv::design::mixBottom,
               *p.apvts.getParameter(P::id::mix), bubble, [](float v) { return juce::String(v * 100.f, 2) + " %"; }),
      inputFader(images, L::Image::whiteFader, L::inputFader.cx, rv::design::gainTop, rv::design::gainBottom,
                 *p.apvts.getParameter(P::id::inputGain), bubble, [](float v) { return P::gainText(v); }),
      outputFader(images, L::Image::whiteFader, L::outputFader.cx, rv::design::gainTop, rv::design::gainBottom,
                  *p.apvts.getParameter(P::id::outputGain), bubble, [](float v) { return P::gainText(v); }),
      presetBar(p.presets, *this),
      linkButton({0.0f, rv::design::barTop + 35.0f, 220.0f, 80.0f}, "LINK", rv::ui::BarButton::Kind::toggle), // x: see below
      helpButton({0.0f, rv::design::barTop + 35.0f, 90.0f, 80.0f}),
      // The product name in the middle of the bar; the plugin's version (JUCE's, from the CMake project) at the right.
      titleLabel({1597.0f, rv::design::barTop + 35.0f, 300.0f, 80.0f}, "REV-246", 42.0f, rv::design::gold),
      versionLabel({0.0f, rv::design::barTop + 35.0f, 200.0f, 80.0f}, juce::String("v") + JucePlugin_VersionString, 30.0f,
                   rv::design::labelDimmed, juce::Justification::centredRight)
{
    setLookAndFeel(&lookAndFeel);

    // The right end of the bar, laid out from the version number's width (fonts differ by platform): the version sits
    // barEdgePad from the panel's right edge, the "?" the same distance to its left, and LINK just left of the "?".
    {
        const auto versionText = juce::String("v") + JucePlugin_VersionString;
        const auto textW = rv::ui::textWidth(rv::ui::labelFont(30.0f), versionText) + 2.0f;
        const auto right = rv::layout::designWidth - rv::design::barEdgePad;
        const auto y = rv::design::barTop + 35.0f;
        versionLabel.setDesignBounds({right - textW, y, textW, 80.0f});
        const auto helpRight = right - textW - rv::design::barEdgePad;
        helpButton.setDesignBounds({helpRight - 90.0f, y, 90.0f, 80.0f});
        linkButton.setDesignBounds({helpRight - 90.0f - 12.0f - 220.0f, y, 220.0f, 80.0f});
    }

    // Programme buttons and lamps.
    struct Prog { L::Slot button, led; L::Image out, in, ledOn, ledOff; };
    const Prog progs[] = {
        {L::sButton, L::sLed, L::Image::sOut, L::Image::sIn, L::Image::redLedOn, L::Image::redLedOff},
        {L::mButton, L::mLed, L::Image::mOut, L::Image::mIn, L::Image::redLedOn, L::Image::redLedOff},
        {L::lButton, L::lLed, L::Image::lOut, L::Image::lIn, L::Image::redLedOn, L::Image::redLedOff},
        {L::eButton, L::eLed, L::Image::eOut, L::Image::eIn, L::Image::redLedOn, L::Image::redLedOff},
        {L::nButton, L::nLed, L::Image::nOut, L::Image::nIn, L::Image::yellowLedOn, L::Image::yellowLedOff}};
    for (int i = 0; i < P::numPrograms; ++i)
    {
        auto b = std::make_unique<rv::ui::ProgramButton>(images, progs[i].button, progs[i].out, progs[i].in);
        b->setTooltip(P::programLongName((P::Program)i));
        b->onClick = [this, i] {
            auto& pr = param(P::id::program);
            pr.beginChangeGesture();
            pr.setValueNotifyingHost(pr.convertTo0to1((float)i));
            pr.endChangeGesture();
        };
        programButtons.push_back(std::move(b));
        leds.push_back(std::make_unique<rv::ui::Led>(images, progs[i].led, progs[i].ledOn, progs[i].ledOff));
    }

    powerButton.onToggle = [this](bool on) {
        auto& pr = param(P::id::power);
        pr.beginChangeGesture();
        pr.setValueNotifyingHost(on ? 1.0f : 0.0f);
        pr.endChangeGesture();
    };

    for (int i = 0; i < joyCount; ++i)
        setUpJoy(i);

    programDisplay.setTooltip("Programme: S short, M medium, L long, E the 250 mode, N non-linear.");
    lowpassDisplay.setTooltip("Lowpass frequency in kHz.");
    // Wheel / drag on a display field nudges its value, like the lever below it (the program display changes program).
    auto nudger = [this](Joy j) {
        return [this, j](int steps) {
            for (int k = 0; k < std::abs(steps); ++k)
                stepJoy(j, steps > 0 ? 1 : -1, 0);
        };
    };
    lowpassDisplay.setNudge(0, nudger(Joy::lowpass));
    reflectionsDisplay.setNudge(0, nudger(Joy::delay));
    reflectionsDisplay.setNudge(1, nudger(Joy::ampl));
    decayDisplay.setNudge(0, nudger(Joy::lf));
    decayDisplay.setNudge(1, nudger(Joy::rt));
    decayDisplay.setNudge(2, nudger(Joy::hf));
    programDisplay.setNudge(0, [this](int steps) {
        auto& pr = param(P::id::program);
        const auto now = juce::roundToInt(pr.convertFrom0to1(pr.getValue()));
        const auto next = juce::jlimit(0, P::numPrograms - 1, now + (steps > 0 ? 1 : -1));
        pr.beginChangeGesture();
        pr.setValueNotifyingHost(pr.convertTo0to1((float)next));
        pr.endChangeGesture();
    });

    // The "?" and the link.
    helpButton.setOn(proc.uiState().getProperty(Rev246Processor::UiProps::tooltipsOn, false));
    helpButton.onToggle = [this](bool on) {
        proc.uiState().setProperty(Rev246Processor::UiProps::tooltipsOn, on, nullptr);
        helpButton.setOn(on);
    };
    linkButton.setTooltip("Link: moving the INPUT or OUTPUT fader moves the other by the same number of dB the opposite way, keeping the overall gain.");
    linkButton.onToggle = [this](bool on) {
        auto& pr = param(P::id::link);
        pr.beginChangeGesture();
        pr.setValueNotifyingHost(on ? 1.0f : 0.0f);
        pr.endChangeGesture();
    };

    panKnob.helpText = "Pan: where the dry signal sits in the stereo field. Double-click to centre; Shift for fine.";
    mixFader.helpText = "Mix: dry to reverb, 0 to 100 %.";
    inputFader.helpText = "Input level, +12 dB to -infinity. Double-click for 0 dB; Shift for fine.";
    outputFader.helpText = "Output level, +12 dB to -infinity. Double-click for 0 dB; Shift for fine.";
    inputFader.onUserMove = [this](float o, float n) { if (param(P::id::link).getValue() >= 0.5f) moveLinked(inputFader, outputFader, o, n); };
    outputFader.onUserMove = [this](float o, float n) { if (param(P::id::link).getValue() >= 0.5f) moveLinked(outputFader, inputFader, o, n); };
    for (auto* f : {&inputFader, &outputFader})
        f->onDragEnd = [this] {
            if (linkedGestureOpen)
            {
                param(P::id::inputGain).endChangeGesture();
                param(P::id::outputGain).endChangeGesture();
                linkedGestureOpen = false;
            }
        };

    // Children, bottom to top: the funnels sit under the levers, the meter over the panel.
    for (auto& f : funnels) designComponents.push_back(f.get());
    for (auto& j : joysticks) designComponents.push_back(j.get());
    for (auto& b : programButtons) designComponents.push_back(b.get());
    for (auto& l : leds) designComponents.push_back(l.get());
    for (auto* c : std::initializer_list<rv::ui::DesignComponent*>{&powerButton, &programDisplay, &lowpassDisplay, &reflectionsDisplay,
                                                                  &decayDisplay, &levelMeter, &panKnob, &mixFader, &inputFader,
                                                                  &outputFader, &linkButton, &helpButton, &titleLabel, &versionLabel})
        designComponents.push_back(c);
    presetBar.addTo(designComponents);
    for (auto* c : designComponents)
        if (c->getParentComponent() == nullptr)
            addAndMakeVisible(*c);
    addChildComponent(bubble);

    for (const auto* id : {P::id::program, P::id::mix, P::id::pan, P::id::inputGain, P::id::outputGain, P::id::power, P::id::link})
        proc.apvts.addParameterListener(id, this);
    for (int i = 0; i < joyCount; ++i)
        proc.apvts.addParameterListener(P::joyId((Joy)i), this);
    proc.presets.onChange = [this] { needsRefresh = true; };

    // Size: remembered width, resizable at a locked aspect ratio.
    const int width = juce::jlimit(juce::roundToInt(rv::layout::designWidth * rv::design::minScale),
                                   juce::roundToInt(rv::layout::designWidth * rv::design::maxScale),
                                   (int)proc.uiState().getProperty(Rev246Processor::UiProps::width, Rev246Processor::defaultEditorWidth));
    const auto aspect = rv::layout::designWidth / rv::layout::designHeight;
    setResizable(true, true);
    setResizeLimits(juce::roundToInt(rv::layout::designWidth * rv::design::minScale),
                    juce::roundToInt(rv::layout::designHeight * rv::design::minScale),
                    juce::roundToInt(rv::layout::designWidth * rv::design::maxScale),
                    juce::roundToInt(rv::layout::designHeight * rv::design::maxScale));
    getConstrainer()->setFixedAspectRatio(aspect);
    setSize(width, juce::roundToInt((float)width / aspect));

    refreshAll();
    startTimerHz(30);
}

Rev246Editor::~Rev246Editor()
{
    stopTimer();
    proc.presets.onChange = nullptr;
    for (const auto* id : {P::id::program, P::id::mix, P::id::pan, P::id::inputGain, P::id::outputGain, P::id::power, P::id::link})
        proc.apvts.removeParameterListener(id, this);
    for (int i = 0; i < joyCount; ++i)
        proc.apvts.removeParameterListener(P::joyId((Joy)i), this);
    setLookAndFeel(nullptr);
}

void Rev246Editor::setUpJoy(int i)
{
    auto stick = std::make_unique<rv::ui::Joystick>(images, joySlots[i]);
    stick->onStep = [this, i](int dir, int hold) { stepJoy((Joy)i, dir, hold); };
    stick->setTooltip(joyHelp(P::Program::s, (Joy)i));
    joysticks.push_back(std::move(stick));
    funnels.push_back(std::make_unique<rv::ui::Funnel>(i));
}

void Rev246Editor::stepJoy(Joy j, int direction, int hold)
{
    auto& pr = param(P::joyId(j));
    const int amount = (j == Joy::delay && hold > 0) ? 10 : 1; // the unit speeds up when a lever is held
    const auto now = juce::roundToInt(pr.convertFrom0to1(pr.getValue()));
    const auto next = juce::jlimit(0, P::joyMaxStep(j), now + direction * amount);
    if (next == now)
        return;
    pr.beginChangeGesture();
    pr.setValueNotifyingHost(pr.convertTo0to1((float)next));
    pr.endChangeGesture();
}

void Rev246Editor::moveLinked(rv::ui::Fader& from, rv::ui::Fader& to, float oldPos, float newPos)
{
    const auto delta = P::gainDbFromPosition(newPos) - P::gainDbFromPosition(oldPos);
    const auto target = P::gainPositionFromDb(P::gainDbFromPosition(to.position()) - delta);
    (void)from;
    if (!linkedGestureOpen)
    {
        param(P::id::inputGain).beginChangeGesture();
        param(P::id::outputGain).beginChangeGesture();
        linkedGestureOpen = true;
    }
    to.setPosition(target, true);
}

void Rev246Editor::timerCallback()
{
    levelMeter.setPeakLinear(proc.takeInputPeak());
    if (needsRefresh.exchange(false))
        refreshAll();
}

void Rev246Editor::refreshAll()
{
    const auto prog = program();
    for (int i = 0; i < P::numPrograms; ++i)
    {
        programButtons[(size_t)i]->setSelected(i == (int)prog);
        leds[(size_t)i]->setOn(i == (int)prog);
    }
    programDisplay.setText(0, P::programName(prog));
    powerButton.setOn(param(P::id::power).getValue() >= 0.5f);
    linkButton.setOn(param(P::id::link).getValue() >= 0.5f);

    std::array<P::JoyValue, joyCount> values;
    for (int i = 0; i < joyCount; ++i)
    {
        const auto& pr = param(P::joyId((Joy)i));
        const auto step = juce::roundToInt(pr.convertFrom0to1(pr.getValue()));
        values[(size_t)i] = P::joyValue(prog, (Joy)i, step);
        funnels[(size_t)i]->setFraction((float)step / (float)P::joyMaxStep((Joy)i));
        funnels[(size_t)i]->setDimmed(!values[(size_t)i].used);
        joysticks[(size_t)i]->setDimmed(!values[(size_t)i].used);
        joysticks[(size_t)i]->setTooltip(joyHelp(prog, (Joy)i));
    }
    lowpassDisplay.setText(0, values[0].text);
    reflectionsDisplay.setText(0, values[1].text);
    reflectionsDisplay.setText(1, values[2].text);
    decayDisplay.setText(0, values[3].text);
    decayDisplay.setText(1, values[4].text);
    decayDisplay.setText(2, values[5].text);

    // Two decimals in the field tooltips, so the exact value is there whatever the display rounds to.
    const bool n = prog == P::Program::n;
    auto tip = [&](const char* name, int i) { return juce::String(name) + ": " + juce::String(values[(size_t)i].value, 2) + " " + values[(size_t)i].unit; };
    reflectionsDisplay.setFieldTooltip(0, tip("Delay", 1));
    reflectionsDisplay.setFieldTooltip(1, tip(n ? "Span" : "Amplitude", 2));
    decayDisplay.setFieldTooltip(0, tip(n ? "Early RT" : "LF x", 3));
    decayDisplay.setFieldTooltip(1, n ? juce::String("RT (not used)") : tip("RT", 4));
    decayDisplay.setFieldTooltip(2, tip(n ? "Late RT" : "HF x", 5));

    presetBar.refresh();
}

void Rev246Editor::paint(juce::Graphics& g)
{
    const auto s = (float)getWidth() / rv::layout::designWidth;
    g.fillAll(rv::design::barColour);
    rv::ui::drawFitted(g, images, rv::layout::Image::main,
                       juce::Rectangle<float>(0.0f, 0.0f, rv::layout::designWidth * s, rv::design::barTop * s));
    // The bottom bar: the panel's own dark with a gold rule where main.png ends.
    const auto barTop = rv::design::barTop * s;
    g.setColour(rv::design::gold.withAlpha(0.55f));
    g.fillRect(juce::Rectangle<float>(0.0f, barTop, (float)getWidth(), juce::jmax(1.0f, 4.0f * s)));
}

void Rev246Editor::resized()
{
    const auto s = (float)getWidth() / rv::layout::designWidth;
    images.clear();
    presetBar.layoutForScale(s);
    for (auto* c : designComponents)
        c->layoutForScale(s);
    bubble.setAlwaysOnTop(true);
    proc.uiState().setProperty(Rev246Processor::UiProps::width, getWidth(), nullptr);
}
