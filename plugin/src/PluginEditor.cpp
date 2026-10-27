#include "PluginEditor.h"

#include "ui/Design.h"
#include "ui/HelpText.h"

#include <array>
#include <cmath>
#include <iterator>
#include <optional>

namespace
{
namespace L = rv::layout;
namespace P = rv::params;
namespace D = rv::design;
using rv::ui::SegmentDisplay;

constexpr L::Slot joySlots[] = {L::lowpassJoy, L::delayJoy, L::amplJoy, L::lfJoy, L::rtJoy, L::hfJoy};
constexpr int joyCount = P::numJoys;
static_assert(std::size(joySlots) == (size_t)joyCount);

// Per programme: its button, lamp and their artwork (N, the non-linear programme, has the yellow lamp).
struct ProgramArt
{
    L::Slot button, led;
    L::Image out, in, ledOn, ledOff;
};
constexpr ProgramArt programArt[] = {
    {L::sButton, L::sLed, L::Image::sOut, L::Image::sIn, L::Image::redLedOn, L::Image::redLedOff},
    {L::mButton, L::mLed, L::Image::mOut, L::Image::mIn, L::Image::redLedOn, L::Image::redLedOff},
    {L::lButton, L::lLed, L::Image::lOut, L::Image::lIn, L::Image::redLedOn, L::Image::redLedOff},
    {L::eButton, L::eLed, L::Image::eOut, L::Image::eIn, L::Image::redLedOn, L::Image::redLedOff},
    {L::nButton, L::nLed, L::Image::nOut, L::Image::nIn, L::Image::yellowLedOn, L::Image::yellowLedOff}};
static_assert(std::size(programArt) == (size_t)P::numPrograms);

// The parameters whose changes the editor shows (the joystick steps are added by the caller).
constexpr const char* watchedParams[] = {P::id::program, P::id::mix, P::id::pan, P::id::stereoWidth, P::id::inputGain,
                                         P::id::outputGain, P::id::power, P::id::link, P::id::transformer,
                                         P::id::trueStereo, P::id::sixteenBit};

// Holding the DELAY lever for more than a second speeds it up to this many steps (ms) per press, as on the unit.
constexpr int delayHoldStep = 10;

// A bottom-bar slot of the given width and left edge, on the bar's button row.
juce::Rectangle<float> barSlot(float width, float left = 0.0f) { return {left, D::barRowTop, width, D::barRowHeight}; }

// A display field of `count` digits starting at digit `first`.
SegmentDisplay::Field displayField(int first, int count, bool leftAlign = false)
{
    SegmentDisplay::Field f;
    f.first = first;
    f.count = count;
    f.leftAlign = leftAlign;
    return f;
}

juce::String versionText() { return juce::String("v") + JucePlugin_VersionString; }

// A plain decimal number, with no unit text or other trailing characters, so a mistyped entry is rejected rather than
// silently truncated by String::getDoubleValue().
bool parseNumber(const juce::String& text, double& value)
{
    const auto s = text.trim();
    if (s.isEmpty())
        return false;
    bool digit = false, point = false;
    for (int i = 0; i < s.length(); ++i)
    {
        const auto ch = s[i];
        if (ch == '+' || ch == '-')
        {
            if (i != 0)
                return false;
        }
        else if (ch == '.')
        {
            if (point)
                return false;
            point = true;
        }
        else if (juce::CharacterFunctions::isDigit(ch))
        {
            digit = true;
        }
        else
        {
            return false;
        }
    }
    if (!digit)
        return false;
    value = s.getDoubleValue();
    return std::isfinite(value);
}
} // namespace


Rev246Editor::PanelLookAndFeel::PanelLookAndFeel()
{
    setColour(juce::TooltipWindow::backgroundColourId, juce::Colour(0xf2161615));
    setColour(juce::TooltipWindow::textColourId, juce::Colour(0xffe6e6e2));
    setColour(juce::TooltipWindow::outlineColourId, juce::Colour(0xff4a4a47));
}

//==============================================================================
Rev246Editor::Rev246Editor(Rev246Processor& p)
    : juce::AudioProcessorEditor(p), proc(p),
      powerButton(images),
      programDisplay(L::programDisplay.bounds(), 5, {displayField(0, 5, true)}),
      lowpassDisplay(L::lowpassDisplay.bounds(), 3, {displayField(0, 3)}),
      reflectionsDisplay(L::reflectionsDisplay.bounds(), 6, {displayField(0, 3), displayField(3, 3)}),
      decayDisplay(L::decayDisplay.bounds(), 9, {displayField(0, 3), displayField(3, 3), displayField(6, 3)}),
      panKnob(images, *p.apvts.getParameter(P::id::pan), bubble,
              [](float v) { return std::abs(v) < 0.005f ? juce::String("Centre") : juce::String(v < 0.f ? "L " : "R ") + juce::String(std::abs(v) * 100.f, 2) + " %"; }),
      widthKnob({0.0f, D::widthKnobTop, D::widthKnobW, D::widthKnobH}, *p.apvts.getParameter(P::id::stereoWidth), bubble,
                [](float v) { return juce::String(juce::roundToInt(v)) + " %"; }),
      mixFader(images, L::Image::redFader, L::mixFader.cx, D::mixTop, D::mixBottom,
               *p.apvts.getParameter(P::id::mix), bubble, [](float v) { return juce::String(v * 100.f, 2) + " %"; }),
      inputFader(images, L::Image::whiteFader, L::inputFader.cx, D::gainTop, D::gainBottom,
                 *p.apvts.getParameter(P::id::inputGain), bubble, [](float v) { return P::gainText(v); }),
      outputFader(images, L::Image::whiteFader, L::outputFader.cx, D::gainTop, D::gainBottom,
                  *p.apvts.getParameter(P::id::outputGain), bubble, [](float v) { return P::gainText(v); }),
      presetBar(p.presets, *this),
      // The bar's buttons are placed by layoutBottomBar(); only their widths are set here.
      linkButton(barSlot(D::linkWidth), "LINK", rv::ui::BarButton::Kind::toggle),
      transformerButton(barSlot(D::transformerWidth), "XFMR", rv::ui::BarButton::Kind::toggle),
      formatButton(barSlot(D::formatWidth), "FLOAT", rv::ui::BarButton::Kind::toggle),
      stereoButton(barSlot(D::stereoModeWidth), "ORIGINAL", rv::ui::BarButton::Kind::toggle),
      calButton(barSlot(D::calWidth, D::calLeft), "CAL -18", rv::ui::BarButton::Kind::momentary),
      helpButton(barSlot(D::helpWidth)),
      // The product name in the middle of the bar; the plugin's version (JUCE's, from the CMake project) at the right.
      titleLabel(barSlot(D::titleWidth, D::titleLeft), "REV-246", D::titleCapHeight, D::gold),
      versionLabel(barSlot(0.0f), versionText(), D::versionCapHeight, D::labelDimmed, juce::Justification::centredRight)
{
    setLookAndFeel(&lookAndFeel);
    setOpaque(true); // paint() covers the whole editor

    layoutBottomBar();
    setUpProgramControls();
    for (int i = 0; i < joyCount; ++i)
        setUpJoy(i);
    setUpDisplays();
    setUpBarButtons();
    setUpKnobsAndFaders();
    addChildren();
    addParameterListeners();
    restoreSize();

    refreshAll();
    startTimerHz(30);
}

Rev246Editor::~Rev246Editor()
{
    stopTimer();
    if (valueDialog != nullptr)
        valueDialog->exitModalState(0); // its callback runs later and finds the editor gone
    endLinkedGesture();
    proc.presets.onChange = nullptr;
    removeParameterListeners();
    setLookAndFeel(nullptr);
}

//==============================================================================
// Construction steps

// The right end of the bar is laid out from the version number's width (fonts differ by platform): the version sits
// barEdgePad from the panel's right edge, then the "?" the same distance to its left, then LINK, XFMR, FLOAT / 16-BIT
// and ORIGINAL / TRUE ST each a barGap apart, and the width halo last. CAL sits at the left end of the bar, after SAVE.
void Rev246Editor::layoutBottomBar()
{
    const auto versionWidth = rv::ui::textWidth(rv::ui::labelFont(D::versionCapHeight), versionText()) + 2.0f;
    versionLabel.setDesignBounds(barSlot(versionWidth));

    auto edge = L::designWidth - D::barEdgePad; // the right edge of the next control
    auto placeLeft = [&edge](rv::ui::DesignComponent& c, float gapAfter) {
        auto b = c.getDesignBounds();
        b.setX(edge - b.getWidth());
        c.setDesignBounds(b);
        edge = b.getX() - gapAfter;
    };
    placeLeft(versionLabel, D::barEdgePad);
    placeLeft(helpButton, D::barGap);
    placeLeft(linkButton, D::barGap);
    placeLeft(transformerButton, D::barGap);
    placeLeft(formatButton, D::barGap);
    placeLeft(stereoButton, D::barGap);
    placeLeft(widthKnob, D::barGap);
}

void Rev246Editor::setUpProgramControls()
{
    for (int i = 0; i < P::numPrograms; ++i)
    {
        const auto& art = programArt[i];
        auto button = std::make_unique<rv::ui::ProgramButton>(images, art.button, art.out, art.in);
        button->setTooltip(rv::ui::help::program((P::Program)i));
        button->onClick = [this, i] { selectProgram(i); };
        programButtons.push_back(std::move(button));
        leds.push_back(std::make_unique<rv::ui::Led>(images, art.led, art.ledOn, art.ledOff));
    }

    powerButton.onToggle = [this](bool on) { setToggle(P::id::power, on); };
}

void Rev246Editor::setUpJoy(int i)
{
    const auto joy = (Joy)i;
    // Tooltips are set by refreshJoys(), which knows the programme.
    auto stick = std::make_unique<rv::ui::Joystick>(images, joySlots[i]);
    stick->onStep = [this, joy](int dir, int hold) { stepJoy(joy, dir, hold); };
    stick->onReset = [this, joy] { resetJoy(joy); };
    joysticks.push_back(std::move(stick));

    auto funnel = std::make_unique<rv::ui::Funnel>(i);
    funnel->onSetFraction = [this, joy](float f) { setJoy(joy, f); };
    funnel->onStep = [this, joy](int dir) { stepJoy(joy, dir, 0); };
    funnel->onReset = [this, joy] { resetJoy(joy); };
    funnels.push_back(std::move(funnel));
}

void Rev246Editor::setUpDisplays()
{
    programDisplay.setTooltip("Programme display: use the wheel or drag to change the program (S short, M medium, L long, E the EMT 250 program, N non-linear), or press the buttons.");
    lowpassDisplay.setTooltip("Lowpass frequency in kHz.");

    // Wheel / drag on a display field nudges its value, like the lever below it; Option/Alt-click resets it;
    // double-click types it.
    struct Field
    {
        SegmentDisplay& display;
        size_t index;
        Joy joy;
    };
    for (const auto& f : {Field{lowpassDisplay, 0, Joy::lowpass}, Field{reflectionsDisplay, 0, Joy::delay},
                          Field{reflectionsDisplay, 1, Joy::ampl}, Field{decayDisplay, 0, Joy::lf},
                          Field{decayDisplay, 1, Joy::rt}, Field{decayDisplay, 2, Joy::hf}})
    {
        const auto joy = f.joy;
        f.display.setNudge(f.index, [this, joy](int steps) { setJoyStep(joy, joyStep(joy) + steps); });
        f.display.setReset(f.index, [this, joy] { resetJoy(joy); });
        f.display.setDoubleClick(f.index, [this, joy] { promptJoyValue(joy); });
    }

    // The programme display changes programme instead.
    programDisplay.setNudge(0, [this](int steps) { selectProgram((int)program() + (steps > 0 ? 1 : -1)); });
}

void Rev246Editor::setUpBarButtons()
{
    // ORIGINAL / TRUE ST and FLOAT / 16-BIT always show the active choice, so keep them lit. The "?" and LINK are
    // ordinary toggles.
    stereoButton.setAlwaysLit(true);
    formatButton.setAlwaysLit(true);
    helpButton.setOn(proc.uiState().getProperty(Rev246Processor::UiProps::tooltipsOn, false));
    helpButton.onToggle = [this](bool on) {
        proc.uiState().setProperty(Rev246Processor::UiProps::tooltipsOn, on, nullptr);
        helpButton.setOn(on);
    };

    linkButton.setTooltip(rv::ui::help::link());
    linkButton.onToggle = [this](bool on) { setToggle(P::id::link, on); };

    transformerButton.setTooltip("Transformer (default on): the unit with its optional input and output transformers "
                                 "(service manual). Adds the low-end colour of the commercial captures (about +3 dB at "
                                 "40-60 Hz) and a guessed transformer distortion, mostly 3rd harmonic at low frequencies "
                                 "and high levels. Off: the stock unit.");
    transformerButton.onToggle = [this](bool on) { setToggle(P::id::transformer, on); };

    formatButton.setTooltip("16-BIT (default): quantises the reverb's internal storage to 16 bits, as the original did, "
                            "for its faint signal-dependent noise floor and grit in the decays, and clips the input where "
                            "the unit's A/D converter did (set by CAL; OVL lights when it clips). FLOAT: full-precision "
                            "processing. Only affects the reverb path.");
    formatButton.onToggle = [this](bool on) { setToggle(P::id::sixteenBit, on); };

    stereoButton.onToggle = [this](bool on) { setToggle(P::id::trueStereo, on); };

    meterCalIndex = proc.getMeterCal();
    calButton.setTooltip("Meter calibration: what 0 VU reads, in dBFS. Click to cycle. -18 dBFS is the default "
                         "(the unit's 0 VU is 1 kHz at 0.775 V). It is the unit's input alignment, so it also sets where "
                         "its A/D clips: 6 dB over 0 VU for a sine below 2 kHz, up to 5 dB sooner at high frequencies. "
                         "16-BIT applies that clipping; OVL shows it in either mode.");
    calButton.setLabel(calLabel());
    calButton.onClick = [this] {
        proc.setMeterCal((proc.getMeterCal() + 1) % D::numMeterCal);
        meterCalIndex = proc.getMeterCal();
        calButton.setLabel(calLabel());
    };
}

void Rev246Editor::setUpKnobsAndFaders()
{
    widthKnob.onDoubleClick = [this] { promptWidthValue(); };
    mixFader.helpText = rv::ui::help::mix();
    inputFader.helpText = rv::ui::help::input();
    outputFader.helpText = rv::ui::help::output();

    // With LINK on, the other gain fader follows; both gains' host gestures stay open until the drag ends.
    inputFader.onUserMove = [this](float o, float n) { if (isOn(P::id::link)) moveLinked(outputFader, o, n); };
    outputFader.onUserMove = [this](float o, float n) { if (isOn(P::id::link)) moveLinked(inputFader, o, n); };
    inputFader.onInteractionEnd = outputFader.onInteractionEnd = [this] { endLinkedGesture(); };
}

// Children, bottom to top: the funnels sit under the levers, the meter over the panel.
void Rev246Editor::addChildren()
{
    for (auto& f : funnels) designComponents.push_back(f.get());
    for (auto& j : joysticks) designComponents.push_back(j.get());
    for (auto& b : programButtons) designComponents.push_back(b.get());
    for (auto& l : leds) designComponents.push_back(l.get());
    for (auto* c : std::initializer_list<rv::ui::DesignComponent*>{&powerButton, &programDisplay, &lowpassDisplay, &reflectionsDisplay,
                                                                  &decayDisplay, &levelMeter, &panKnob, &widthKnob, &mixFader, &inputFader,
                                                                  &outputFader, &linkButton, &transformerButton, &formatButton, &stereoButton, &calButton, &helpButton,
                                                                  &titleLabel, &versionLabel})
        designComponents.push_back(c);
    presetBar.addTo(designComponents);
    for (auto* c : designComponents)
        if (c->getParentComponent() == nullptr) // the preset bar's buttons are already in the editor
            addAndMakeVisible(*c);
    addChildComponent(bubble); // last, so it is on top
}

void Rev246Editor::addParameterListeners()
{
    for (const auto* id : watchedParams)
        proc.apvts.addParameterListener(id, this);
    for (int i = 0; i < joyCount; ++i)
        proc.apvts.addParameterListener(P::joyId((Joy)i), this);
    proc.presets.onChange = [this] { needsRefresh = true; };
}

void Rev246Editor::removeParameterListeners()
{
    for (const auto* id : watchedParams)
        proc.apvts.removeParameterListener(id, this);
    for (int i = 0; i < joyCount; ++i)
        proc.apvts.removeParameterListener(P::joyId((Joy)i), this);
}

// Size: the remembered width, resizable at a locked aspect ratio.
void Rev246Editor::restoreSize()
{
    const auto minW = juce::roundToInt(L::designWidth * D::minScale), maxW = juce::roundToInt(L::designWidth * D::maxScale);
    const auto minH = juce::roundToInt(L::designHeight * D::minScale), maxH = juce::roundToInt(L::designHeight * D::maxScale);
    const int width = juce::jlimit(minW, maxW, (int)proc.uiState().getProperty(Rev246Processor::UiProps::width, Rev246Processor::defaultEditorWidth));
    const auto aspect = L::designWidth / L::designHeight;
    setResizable(true, true);
    setResizeLimits(minW, minH, maxW, maxH);
    getConstrainer()->setFixedAspectRatio(aspect);
    setSize(width, juce::roundToInt((float)width / aspect));
}

//==============================================================================
// Parameter access

juce::RangedAudioParameter& Rev246Editor::param(const char* id) const
{
    return *proc.apvts.getParameter(id);
}

int Rev246Editor::joyStep(Joy j) const
{
    const auto& pr = param(P::joyId(j));
    return juce::roundToInt(pr.convertFrom0to1(pr.getValue()));
}

P::Program Rev246Editor::program() const
{
    const auto& p = param(P::id::program);
    return (P::Program)juce::jlimit(0, P::numPrograms - 1, juce::roundToInt(p.convertFrom0to1(p.getValue())));
}

void Rev246Editor::setParamValue(juce::RangedAudioParameter& pr, float plainValue)
{
    pr.beginChangeGesture();
    pr.setValueNotifyingHost(pr.convertTo0to1(plainValue));
    pr.endChangeGesture();
}

void Rev246Editor::selectProgram(int index)
{
    setParamValue(param(P::id::program), (float)juce::jlimit(0, P::numPrograms - 1, index));
}

void Rev246Editor::setJoyStep(Joy j, int step)
{
    step = juce::jlimit(0, P::joyMaxStep(j), step);
    if (step != joyStep(j))
        setParamValue(param(P::joyId(j)), (float)step);
}

void Rev246Editor::stepJoy(Joy j, int direction, int hold)
{
    const int amount = (j == Joy::delay && hold > 0) ? delayHoldStep : 1; // the unit speeds up when a lever is held
    setJoyStep(j, joyStep(j) + direction * amount);
}

void Rev246Editor::setJoy(Joy j, float fraction)
{
    setJoyStep(j, juce::roundToInt(juce::jlimit(0.0f, 1.0f, fraction) * (float)P::joyMaxStep(j)));
}

void Rev246Editor::resetJoy(Joy j)
{
    setJoyStep(j, P::joyDefaultStep(j));
}

void Rev246Editor::moveLinked(rv::ui::Fader& to, float oldPos, float newPos)
{
    const auto delta = P::gainDbFromPosition(newPos) - P::gainDbFromPosition(oldPos);
    const auto target = P::gainPositionFromDb(P::gainDbFromPosition(to.position()) - delta);
    if (!linkedGestureOpen)
    {
        param(P::id::inputGain).beginChangeGesture();
        param(P::id::outputGain).beginChangeGesture();
        linkedGestureOpen = true;
    }
    to.setPosition(target, true);
}

void Rev246Editor::endLinkedGesture()
{
    if (!linkedGestureOpen)
        return;
    param(P::id::inputGain).endChangeGesture();
    param(P::id::outputGain).endChangeGesture();
    linkedGestureOpen = false;
}

juce::String Rev246Editor::calLabel() const
{
    return "CAL " + juce::String(D::meterCalDbFs[meterCalIndex], 0);
}

//==============================================================================
// Typed entry

void Rev246Editor::showValueDialog(const juce::String& title, const juce::String& message, const juce::String& shown,
                                   const juce::String& hint, std::function<bool(const juce::String&)> onApply,
                                   std::function<void(const juce::String&)> onInvalid)
{
    auto* w = new juce::AlertWindow(title, message, juce::MessageBoxIconType::NoIcon, this);
    w->addTextEditor("value", shown);
    if (auto* editor = w->getTextEditor("value"))
        editor->selectAll();
    w->addButton("Set", 1, juce::KeyPress(juce::KeyPress::returnKey));
    w->addButton("Cancel", 0, juce::KeyPress(juce::KeyPress::escapeKey));
    valueDialog = w;

    // The dialog outlives neither the window nor the editor: the callbacks hold SafePointers.
    juce::Component::SafePointer<Rev246Editor> self(this);
    juce::Component::SafePointer<juce::AlertWindow> dialog(w);
    w->enterModalState(true, juce::ModalCallbackFunction::create([self, dialog, hint, apply = std::move(onApply), retry = std::move(onInvalid)](int result) {
                           if (result != 1 || self == nullptr || dialog == nullptr)
                               return;
                           const auto text = dialog->getTextEditorContents("value");
                           if (apply(text))
                               return;
                           juce::AlertWindow::showAsync(
                               juce::MessageBoxOptions()
                                   .withIconType(juce::MessageBoxIconType::WarningIcon)
                                   .withTitle("Invalid value")
                                   .withMessage("\"" + text.trim() + "\" is not valid for this control.\n" + hint)
                                   .withButton("OK")
                                   .withAssociatedComponent(self.getComponent()),
                               [self, retry, text](int) {
                                   if (self != nullptr)
                                       retry(text);
                               });
                       }),
                       true);
}

void Rev246Editor::promptJoyValue(Joy j, const juce::String& initial)
{
    const auto prog = program();
    if (j == Joy::rt && prog == P::Program::n)
    {
        juce::AlertWindow::showMessageBoxAsync(juce::MessageBoxIconType::InfoIcon, "RT",
                                                "RT is not used in the non-linear program.", "OK", this);
        return;
    }

    const auto current = P::joyValue(prog, j, juce::jlimit(0, P::joyMaxStep(j), joyStep(j)));
    auto currentText = current.text.trim();
    if (prog == P::Program::n && j == Joy::ampl)
        currentText = juce::String(current.value, 2); // SPAN has 0.05 steps; the display only shows one decimal
    const auto hint = P::joyInputHint(prog, j);
    showValueDialog("Enter value",
                    "Enter the " + P::joyInputName(prog, j) + ".\nCurrent: " + currentText + " " + current.unit + ".\n" + hint,
                    initial.isEmpty() ? currentText : initial, hint,
                    [this, prog, j](const juce::String& text) {
                        double value = 0.0;
                        const auto step = parseNumber(text, value) ? P::joyStepFromValue(prog, j, value) : std::nullopt;
                        if (!step)
                            return false;
                        setJoyStep(j, *step);
                        return true;
                    },
                    [this, j](const juce::String& text) { promptJoyValue(j, text); });
}

void Rev246Editor::promptWidthValue(const juce::String& initial)
{
    auto& pr = param(P::id::stereoWidth);
    const auto current = juce::String(juce::roundToInt(pr.convertFrom0to1(pr.getValue())));
    const juce::String hint("Valid: whole numbers from 0 to 200; 100 = unchanged (0 = mono, 200 = extra wide).");
    showValueDialog("Enter width", "Enter the stereo width.\nCurrent: " + current + " %.\n" + hint,
                    initial.isEmpty() ? current : initial, hint,
                    [this](const juce::String& text) {
                        double value = 0.0;
                        if (!parseNumber(text, value) || value < 0.0 || value > 200.0
                            || std::abs(value - (double)juce::roundToInt(value)) >= 1e-6)
                            return false;
                        setParamValue(param(P::id::stereoWidth), (float)juce::roundToInt(value));
                        return true;
                    },
                    [this](const juce::String& text) { promptWidthValue(text); });
}

//==============================================================================
// Refresh

void Rev246Editor::timerCallback()
{
    if (proc.getMeterCal() != meterCalIndex) // e.g. a session or state restored with the editor open
    {
        meterCalIndex = proc.getMeterCal();
        calButton.setLabel(calLabel());
    }
    levelMeter.setLevel(proc.getInputVu(), proc.takeInputClip(), D::meterCalDbFs[meterCalIndex]);

    // The host may change the bus layout while the editor is open; the stereo button and width halo follow it.
    if (proc.isStereoToStereo() != stereoAvailable || (proc.getTotalNumOutputChannels() > 1) != widthAvailable)
        needsRefresh = true;
    if (needsRefresh.exchange(false))
        refreshAll();
}

void Rev246Editor::refreshAll()
{
    const auto prog = program();
    for (int i = 0; i < P::numPrograms; ++i)
    {
        leds[(size_t)i]->setOn(i == (int)prog);
    }
    programDisplay.setText(0, P::programName(prog));
    powerButton.setOn(isOn(P::id::power));
    linkButton.setOn(isOn(P::id::link));
    transformerButton.setOn(isOn(P::id::transformer));
    const bool sixteenBit = isOn(P::id::sixteenBit);
    formatButton.setOn(sixteenBit);
    formatButton.setLabel(sixteenBit ? "16-BIT" : "FLOAT");
    stereoAvailable = proc.isStereoToStereo();
    stereoButton.setEnabled(stereoAvailable);
    widthAvailable = proc.getTotalNumOutputChannels() > 1;
    widthKnob.setEnabled(widthAvailable);
    const bool useTrueStereo = isOn(P::id::trueStereo);
    stereoButton.setOn(useTrueStereo);
    stereoButton.setLabel(useTrueStereo ? "TRUE ST" : "ORIGINAL");

    refreshStateTips(useTrueStereo);
    refreshJoys(prog);
    presetBar.refresh();
}

// The tooltips of the controls whose meaning depends on the bus layout and the stereo choice. Rebuilt only when one of
// those changes, because refreshAll() runs on every parameter change.
void Rev246Editor::refreshStateTips(bool useTrueStereo)
{
    const int key = (stereoAvailable ? 1 : 0) | (widthAvailable ? 2 : 0) | (useTrueStereo ? 4 : 0);
    if (key == shownTipKey)
        return;
    shownTipKey = key;
    const rv::ui::help::Routing routing{stereoAvailable, widthAvailable, useTrueStereo};
    stereoButton.setTooltip(rv::ui::help::stereoMode(routing));
    panKnob.helpText = rv::ui::help::pan(routing);
    widthKnob.helpText = rv::ui::help::width(routing);
}

// The levers, funnels and displays. Any parameter change (a MIX drag, say) refreshes the editor, so this skips the
// work, which builds strings, unless the programme or a lever step changed.
void Rev246Editor::refreshJoys(P::Program prog)
{
    std::array<int, joyCount> steps;
    for (int i = 0; i < joyCount; ++i)
        steps[(size_t)i] = joyStep((Joy)i);
    if ((int)prog == shownProgram && steps == shownSteps)
        return;
    shownProgram = (int)prog;
    shownSteps = steps;

    std::array<P::JoyValue, joyCount> values;
    for (int i = 0; i < joyCount; ++i)
    {
        const auto step = steps[(size_t)i];
        values[(size_t)i] = P::joyValue(prog, (Joy)i, step);
        funnels[(size_t)i]->setFraction((float)step / (float)P::joyMaxStep((Joy)i));
        funnels[(size_t)i]->setDimmed(!values[(size_t)i].used);
        joysticks[(size_t)i]->setDimmed(!values[(size_t)i].used);
        joysticks[(size_t)i]->setTooltip(rv::ui::help::joy(prog, (Joy)i));
        funnels[(size_t)i]->setTooltip(rv::ui::help::joy(prog, (Joy)i));
    }
    lowpassDisplay.setText(0, values[0].text);
    reflectionsDisplay.setText(0, values[1].text);
    reflectionsDisplay.setText(1, values[2].text);
    decayDisplay.setText(0, values[3].text);
    decayDisplay.setText(1, values[4].text);
    decayDisplay.setText(2, values[5].text);

    // Two decimals in the field tooltips, so the exact value is there whatever the display rounds to. Double-click
    // opens the typed-entry dialog.
    const bool n = prog == P::Program::n;
    auto tip = [&](const char* name, int i) {
        return juce::String(name) + ": " + juce::String(values[(size_t)i].value, 2) + " " + values[(size_t)i].unit
               + ". Double-click to type.";
    };
    lowpassDisplay.setFieldTooltip(0, tip("Lowpass", 0));
    reflectionsDisplay.setFieldTooltip(0, tip("Delay", 1));
    reflectionsDisplay.setFieldTooltip(1, tip(n ? "Span" : "Amplitude", 2));
    decayDisplay.setFieldTooltip(0, tip(n ? "Early RT" : "LF x", 3));
    decayDisplay.setFieldTooltip(1, n ? juce::String("RT (not used).") : tip("RT", 4));
    decayDisplay.setFieldTooltip(2, tip(n ? "Late RT" : "HF x", 5));
}

//==============================================================================
void Rev246Editor::paint(juce::Graphics& g)
{
    const auto s = (float)getWidth() / L::designWidth;
    const auto barTop = D::barTop * s;
    g.fillAll(D::barColour);
    rv::ui::drawFitted(g, images, L::Image::main, juce::Rectangle<float>(0.0f, 0.0f, L::designWidth * s, barTop));
    // The bottom bar: the panel's own dark with a gold rule where main.png ends.
    g.setColour(D::gold.withAlpha(D::barRuleAlpha));
    g.fillRect(juce::Rectangle<float>(0.0f, barTop, (float)getWidth(), juce::jmax(1.0f, D::barRuleThickness * s)));
}

void Rev246Editor::resized()
{
    const auto s = (float)getWidth() / L::designWidth;
    images.clear();
    presetBar.layoutForScale(s);
    for (auto* c : designComponents)
        c->layoutForScale(s);
    proc.uiState().setProperty(Rev246Processor::UiProps::width, getWidth(), nullptr);
}
