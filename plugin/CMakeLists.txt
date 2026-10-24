cmake_minimum_required(VERSION 3.15)

project(REV246 VERSION 0.1.0)   # release.yml greps this line for the version; keep the format

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum macOS deployment version")
# macOS ships two separate builds, not a universal binary: arm64, the default here and locally, and x86_64
# ("Intel"), which release.yml cross-compiles with -DCMAKE_OSX_ARCHITECTURES=x86_64.

option(RV_BUILD_TESTS "Build the unit tests (tests/, Catch2 v3)" OFF)
option(RV_STANDALONE "Also build a Standalone app, for checking the UI without a host" OFF)

# -Wextra is a GCC/Clang flag; MSVC's cl.exe misparses it as the /W switch (error D8021).
if (MSVC)
    set(RV_WARNING_FLAGS /W4)
else()
    set(RV_WARNING_FLAGS -Wall -Wextra)
endif()

add_subdirectory(libs/JUCE)

# UI artwork (built from ../ui by tools/build_assets.sh). Resources are looked up by file name
# (src/ui/Assets.cpp), so new files only need to land in assets/images.
file(GLOB RV_ASSET_FILES CONFIGURE_DEPENDS assets/images/*.png)
juce_add_binary_data(REV246Assets SOURCES ${RV_ASSET_FILES})

# Plugin sources, shared by the plugin and the test runner.
set(RV_SOURCES
    src/PluginProcessor.cpp
    src/PluginEditor.cpp
    src/params/Parameters.cpp
    src/presets/PresetManager.cpp
    src/ui/Assets.cpp
    src/ui/BottomBar.cpp
    src/ui/Buttons.cpp
    src/ui/Fader.cpp
    src/ui/Joystick.cpp
    src/ui/LevelMeter.cpp
    src/ui/SegmentDisplay.cpp
)

set(RV_FORMATS AU VST3)
if (RV_STANDALONE)
    list(APPEND RV_FORMATS Standalone)
endif()

juce_add_plugin(REV246
    COMPANY_NAME "Leigh Pierce"
    BUNDLE_ID com.leighpierce.rev246
    PLUGIN_MANUFACTURER_CODE Lprc   # 4 chars, one uppercase
    PLUGIN_CODE Rv24                # 4 chars, one uppercase, unique per plugin
    FORMATS ${RV_FORMATS}
    PRODUCT_NAME "REV-246"
    IS_SYNTH FALSE
    NEEDS_MIDI_INPUT FALSE
    NEEDS_MIDI_OUTPUT FALSE
    VST3_CATEGORIES Fx Reverb
    AU_MAIN_TYPE kAudioUnitType_Effect
)

target_sources(REV246 PRIVATE ${RV_SOURCES})
target_include_directories(REV246 PRIVATE src)

target_compile_definitions(REV246
    PUBLIC
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
)

target_link_libraries(REV246
    PRIVATE
        REV246Assets
        juce::juce_audio_utils
        juce::juce_dsp
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags
)

target_compile_options(REV246 PRIVATE ${RV_WARNING_FLAGS})

# --- Tests ---------------------------------------------------------------------------------------
# One Catch2 runner, built as a JUCE console app from the same sources as the plugin (the state and editor tests
# need the real processor and editor). It compiles its own copy of the JUCE modules, so it roughly doubles a clean
# build: off by default locally, on in CI. Run with ctest.
if (RV_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.16.0
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(Catch2)

    enable_testing()

    juce_add_console_app(REV246Tests PRODUCT_NAME "REV246Tests")
    target_sources(REV246Tests PRIVATE
        ${RV_SOURCES}
        tests/main.cpp
        tests/ParameterTests.cpp
        tests/StateTests.cpp
        tests/ProcessorTests.cpp
        tests/CaptureEqTests.cpp
        tests/EditorTests.cpp
    )
    target_include_directories(REV246Tests PRIVATE src tests)
    target_compile_definitions(REV246Tests PRIVATE
        _USE_MATH_DEFINES
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_STANDALONE_APPLICATION=1
        JUCE_MODAL_LOOPS_PERMITTED=1
        JucePlugin_Name="REV-246"
        JucePlugin_VersionString="${PROJECT_VERSION}"
    )
    target_link_libraries(REV246Tests PRIVATE
        REV246Assets
        juce::juce_audio_utils
        juce::juce_dsp
        juce::juce_recommended_config_flags
        juce::juce_recommended_warning_flags
        Catch2::Catch2
    )
    target_compile_options(REV246Tests PRIVATE ${RV_WARNING_FLAGS})
    add_test(NAME REV246Tests COMMAND REV246Tests)
endif()
