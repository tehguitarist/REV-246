#include "ui/DuckPanel.h"

#include "params/Parameters.h"
#include "ui/Assets.h"
#include "ui/Design.h"

#include <cmath>

namespace rv::ui
{
namespace
{
namespace P = rv::params;

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
            digit = true;
        else
            return false;
    }
    if (!digit)
        return false;
    value = s.getDoubleValue();
    return std::isfinite(value);
}

// One boxed value: drag it up/down or double-click to type. The label is only for display; the component handles the
// mouse so a single click starts a drag and a double-click opens the text editor.
class ValueField final : public juce::Component
{
  public:
    ValueField(juce::RangedAudioParameter& p, bool negativeOnlyIn)
        : param(p), negativeOnly(negativeOnlyIn),
          attachment(p, [this](float) { updateText(); })
    {
        value.setJustificationType(juce::Justification::centred);
        value.setFont(labelFont(18.0f));
        value.setColour(juce::Label::textColourId, design::labelSelected);
        value.setEditable(false, false, true);
        value.setInterceptsMouseClicks(false, false);
        value.onTextChange = [this] { applyText(); };
        addAndMakeVisible(value);
        setMouseCursor(juce::MouseCursor::UpDownResizeCursor);
        updateText();
        attachment.sendInitialUpdate();
    }

    ~ValueField() override { endDrag(); }

    void paint(juce::Graphics& g) override
    {
        const auto b = getLocalBounds().toFloat().reduced(1.0f);
        g.setColour(design::gold.withAlpha(0.10f));
        g.fillRoundedRectangle(b, 6.0f);
        g.setColour(design::gold.withAlpha(0.65f));
        g.drawRoundedRectangle(b, 6.0f, 1.5f);
    }

    void resized() override { value.setBounds(getLocalBounds().reduced(6, 2)); }

    void mouseDown(const juce::MouseEvent& e) override
    {
        if (value.isBeingEdited())
            return;
        dragStartNorm = param.getValue();
        dragStartY = e.y;
        param.beginChangeGesture();
        dragging = true;
    }

    void mouseDrag(const juce::MouseEvent& e) override
    {
        if (!dragging)
            return;
        auto delta = (float)(dragStartY - e.y) / 260.0f;
        if (e.mods.isShiftDown())
            delta *= 0.1f;
        param.setValueNotifyingHost(juce::jlimit(0.0f, 1.0f, dragStartNorm + delta));
    }

    void mouseUp(const juce::MouseEvent&) override { endDrag(); }

    void mouseDoubleClick(const juce::MouseEvent&) override
    {
        endDrag();
        value.showEditor();
        if (auto* editor = value.getCurrentTextEditor())
            editor->selectAll();
    }

  private:
    void endDrag()
    {
        if (!dragging)
            return;
        dragging = false;
        param.endChangeGesture();
    }

    void updateText()
    {
        updating = true;
        value.setText(param.getText(param.getValue(), 32), juce::dontSendNotification);
        updating = false;
    }

    void applyText()
    {
        if (updating)
            return;
        double real = 0.0;
        const auto& range = param.getNormalisableRange();
        if (!parseNumber(value.getText(), real) || real < (double)range.start || real > (double)range.end
            || (negativeOnly && real > 0.0))
        {
            updateText();
            return;
        }
        real = std::round(real * 100.0) / 100.0;
        const auto normalised = param.convertTo0to1((float)real);
        param.beginChangeGesture();
        param.setValueNotifyingHost(normalised);
        param.endChangeGesture();
        updateText();
    }

    juce::RangedAudioParameter& param;
    bool negativeOnly = false;
    juce::ParameterAttachment attachment;
    juce::Label value;
    bool dragging = false, updating = false;
    float dragStartNorm = 0.0f;
    int dragStartY = 0;
};

class DuckPanel final : public juce::Component
{
  public:
    explicit DuckPanel(juce::AudioProcessorValueTreeState& apvts)
        : enabledParam(*apvts.getParameter(P::id::duckEnabled)),
          enabledAttachment(enabledParam, [this](float) { refreshEnableButton(); }),
          rangeField(*apvts.getParameter(P::id::duckRangeDb), true),
          thresholdField(*apvts.getParameter(P::id::duckThresholdDb), true),
          attackField(*apvts.getParameter(P::id::duckAttackMs), false),
          releaseField(*apvts.getParameter(P::id::duckReleaseMs), false)
    {
        title.setText("DUCK", juce::dontSendNotification);
        title.setFont(labelFont(22.0f));
        title.setColour(juce::Label::textColourId, design::gold);
        title.setJustificationType(juce::Justification::centredLeft);
        addAndMakeVisible(title);

        enableButton.setButtonText("Enable");
        enableButton.setColour(juce::TextButton::buttonColourId, design::barColour.brighter(0.08f));
        enableButton.setColour(juce::TextButton::textColourOffId, design::gold);
        enableButton.setColour(juce::TextButton::textColourOnId, design::labelSelected);
        enableButton.onClick = [this] { setDuckEnabled(enabledParam.getValue() < 0.5f); };
        addAndMakeVisible(enableButton);

        addName(rangeName, "RANGE");
        addName(thresholdName, "THRESHOLD");
        addName(attackName, "ATTACK");
        addName(releaseName, "RELEASE");
        addAndMakeVisible(rangeField);
        addAndMakeVisible(thresholdField);
        addAndMakeVisible(attackField);
        addAndMakeVisible(releaseField);

        hint.setText("Keys off the input signal; lowers the wet while the input is over the threshold.",
                     juce::dontSendNotification);
        hint.setFont(labelFont(13.0f));
        hint.setColour(juce::Label::textColourId, design::labelDimmed);
        hint.setJustificationType(juce::Justification::centredLeft);
        addAndMakeVisible(hint);

        refreshEnableButton();
        setSize(370, 268);
    }

    void paint(juce::Graphics& g) override
    {
        const auto b = getLocalBounds().toFloat().reduced(0.5f);
        g.setColour(design::barColour);
        g.fillRoundedRectangle(b, 8.0f);
        g.setColour(design::gold.withAlpha(0.55f));
        g.drawRoundedRectangle(b, 8.0f, 1.5f);
    }

    void resized() override
    {
        auto r = getLocalBounds().reduced(16);
        title.setBounds(r.removeFromTop(26));
        r.removeFromTop(6);
        enableButton.setBounds(r.removeFromTop(34));
        r.removeFromTop(10);
        constexpr int rowH = 30, gap = 7, nameW = 108;
        auto row = [&](juce::Label& name, ValueField& field) {
            auto line = r.removeFromTop(rowH);
            name.setBounds(line.removeFromLeft(nameW));
            field.setBounds(line);
            r.removeFromTop(gap);
        };
        row(rangeName, rangeField);
        row(thresholdName, thresholdField);
        row(attackName, attackField);
        row(releaseName, releaseField);
        hint.setBounds(r.removeFromTop(30));
    }

  private:
    void addName(juce::Label& name, const char* text)
    {
        name.setText(text, juce::dontSendNotification);
        name.setFont(labelFont(15.0f));
        name.setColour(juce::Label::textColourId, design::labelDimmed);
        name.setJustificationType(juce::Justification::centredLeft);
        addAndMakeVisible(name);
    }

    void setDuckEnabled(bool on)
    {
        enabledParam.beginChangeGesture();
        enabledParam.setValueNotifyingHost(on ? 1.0f : 0.0f);
        enabledParam.endChangeGesture();
    }

    void refreshEnableButton()
    {
        const auto on = enabledParam.getValue() >= 0.5f;
        enableButton.setButtonText(on ? "Disable" : "Enable");
        enableButton.setColour(juce::TextButton::buttonOnColourId, design::gold.withAlpha(0.22f));
        enableButton.setToggleState(on, juce::dontSendNotification);
    }

    juce::RangedAudioParameter& enabledParam;
    juce::ParameterAttachment enabledAttachment;
    juce::TextButton enableButton;
    juce::Label title, rangeName, thresholdName, attackName, releaseName, hint;
    ValueField rangeField, thresholdField, attackField, releaseField;
};
} // namespace

std::unique_ptr<juce::Component> createDuckPanel(juce::AudioProcessorValueTreeState& apvts)
{
    return std::make_unique<DuckPanel>(apvts);
}
} // namespace rv::ui
