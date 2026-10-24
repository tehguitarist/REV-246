; REV-246 VST3 installer (NSIS). Windows has no AU format, so there's no plugin-type choice
; here — VST3 only, installed to the shared system VST3 folder.
;
; Build with (from repo root, after building REV246_VST3):
;   bash installer/stage_docs.sh build/docs 0.4.0
;   makensis /DVERSION=0.4.0 /DARTEFACTS_DIR=build\REV246_artefacts\Release\VST3 /DDOCS_DIR=build\docs installer\windows\REV246.nsi
;
; ARTEFACTS_DIR should point at the directory CONTAINING "REV-246.vst3" (i.e. the VST3 release
; output folder), not the bundle itself. DOCS_DIR is the readme and licences staged by stage_docs.sh,
; installed to Program Files\REV-246 (not into the signed bundle).

!ifndef VERSION
  !define VERSION "0.0.0"
!endif
!ifndef ARTEFACTS_DIR
  !define ARTEFACTS_DIR "..\..\build\REV246_artefacts\Release\VST3"
!endif
!ifndef DOCS_DIR
  !define DOCS_DIR "..\..\build\docs"
!endif

Name "REV-246"
OutFile "REV246-Windows-v${VERSION}-Installer.exe"
InstallDir "$COMMONFILES64\VST3"
RequestExecutionLevel admin
SetCompressor /SOLID lzma

Page directory
Page instfiles
UninstPage uninstConfirm
UninstPage instfiles

Section "REV-246 VST3 Plugin" SecVST3
    SetOutPath "$INSTDIR\REV-246.vst3"
    File /r "${ARTEFACTS_DIR}\REV-246.vst3\*.*"

    WriteUninstaller "$INSTDIR\REV-246.vst3\Uninstall-REV246.exe"

    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\REV-246" \
        "DisplayName" "REV-246 VST3"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\REV-246" \
        "UninstallString" "$INSTDIR\REV-246.vst3\Uninstall-REV246.exe"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\REV-246" \
        "DisplayVersion" "${VERSION}"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\REV-246" \
        "Publisher" "Leigh Pierce"
SectionEnd

Section "Read me and licences" SecDocs
    SectionIn RO
    SetOutPath "$PROGRAMFILES64\REV-246"
    File /r "${DOCS_DIR}\*.*"
SectionEnd

Section "Uninstall"
    RMDir /r "$INSTDIR\REV-246.vst3"
    RMDir /r "$PROGRAMFILES64\REV-246"
    DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\REV-246"
SectionEnd
