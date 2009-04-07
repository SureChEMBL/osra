!include Sections.nsh

; The name of the installer
Name "osra"

; The file to write
OutFile "osra-install.exe"

; The default installation directory
InstallDir $PROGRAMFILES\osra

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\osra" "Install_Dir"

;LicenseData "license.txt"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

;--------------------------------

; Pages

;Page license
Page components
Page directory
Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

;--------------------------------

; The stuff to install
Section "osra (required)"

  SectionIn RO
  
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  
  ; Put file there
  File "osra.exe"
  File "delegates.xml"
  File "README.txt"
  File "libopenbabel-3.dll"
  
  strcpy $3 "GPL Ghostscript"
  call checkSoftVersion
  strcmp $2 "" no_gs
  call getGhostscriptInstPath
  strcmp $1 "" no_gs +2
  no_gs:
  MessageBox MB_OK "Ghostscript interpreter not found" IDOK 0 
  call createOSRAbat
  
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\osra "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "DisplayName" "osra"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
  
SectionEnd

Section /o "Symyx Draw plugin" symyx_draw
 call MakeSureIGotGFL
 strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
 call CheckSoftVersion
 strcmp $2 "" no_symyx
 call getSymyxPath
 strcmp $1 "" no_symyx +2
 no_symyx:
  MessageBox MB_OK "Symyx Draw not found" IDOK done
 SetOutPath "$1\AddIns"
 File "plugins\symyx_draw\README.txt"
 File "plugins\symyx_draw\OSRAAction.xml"
 SetOutPath "$1\AddIns\OSRAAction"
 File "plugins\symyx_draw\OSRAAction\OSRAAction.dll"
 File "plugins\symyx_draw\OSRAAction\Interop.GflAx.dll"
 call createXMLConfig
 done:
SectionEnd

; Uninstaller

Section "Uninstall"
  ReadRegStr $0 HKLM SOFTWARE\osra "Install_Dir"
  strcpy $INSTDIR $0
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra"
  DeleteRegKey HKLM SOFTWARE\osra

  ; Remove files and uninstaller
  Delete $INSTDIR\osra.exe
  Delete $INSTDIR\delegates.xml
  Delete $INSTDIR\README.txt
  Delete $INSTDIR\osra.bat
  Delete $INSTDIR\libopenbabel-3.dll
  Delete $INSTDIR\uninstall.exe
  RMDir "$INSTDIR"
  strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
  call un.CheckSoftVersion
  strcmp $2 "" no_symyx
  call un.getSymyxPath
  strcmp $1 "" no_symyx 
  Delete "$1\AddIns\README.txt"
  Delete "$1\AddIns\OSRAAction.xml"
  Delete "$1\AddIns\OSRAAction\OSRAAction.dll"
  Delete "$1\AddIns\OSRAAction\Interop.GflAx.dll"
  Delete "$1\AddIns\OSRAAction\OSRAAction.dll.config"
  RMDir "$1\AddIns\OSRAAction"
  no_symyx:
SectionEnd

Function CheckSoftVersion
StrCpy $0 0
StrCpy $2 ""
loop:
  EnumRegKey $1 HKLM "Software\$3" $0
  StrCmp $1 "" done
  StrCpy $2 $1
  IntOp $0 $0 + 1
  Goto loop
done:
; $2 contains the version of Soft now or empty
FunctionEnd

Function un.CheckSoftVersion
StrCpy $0 0
StrCpy $2 ""
loop:
  EnumRegKey $1 HKLM "Software\$3" $0
  StrCmp $1 "" done
  StrCpy $2 $1
  IntOp $0 $0 + 1
  Goto loop
done:
; $2 contains the version of Soft now or empty
FunctionEnd

Function getSoftInstPath

  Push $0
  Push $1
  Push $2
  ReadRegStr $0 HKLM \
     "Software\Microsoft\Windows\CurrentVersion\Uninstall\$3" \ 
     "UninstallString"
  StrCmp $0 "" fin
  
    StrCpy $1 $0 1 0 ; get firstchar
    StrCmp $1 '"' "" getparent 
      ; if first char is ", let's remove "'s first.
      StrCpy $0 $0 "" 1
      StrCpy $1 0
      rqloop:
        StrCpy $2 $0 1 $1
        StrCmp $2 '"' rqdone
        StrCmp $2 "" rqdone
        IntOp $1 $1 + 1
        Goto rqloop
      rqdone:
      StrCpy $0 $0 $1
    getparent:
    ; the uninstall string goes to an EXE, let's get the directory.
    StrCpy $1 -1
    gploop:
      StrCpy $2 $0 1 $1
      StrCmp $2 "" gpexit
      StrCmp $2 "\" gpexit
      IntOp $1 $1 - 1
      Goto gploop
    gpexit:
    StrCpy $0 $0 $1

    StrCmp $0 "" fin
    IfFileExists $0\$4 fin
      StrCpy $0 ""
  fin:
  Pop $2
  Pop $1
  Exch $0
  
FunctionEnd

Function getGhostscriptInstPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\GPL Ghostscript\$2" \ 
     "GS_DLL"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 -16
  IfFileExists $1\bin\gswin32c.exe fin
  StrCpy $1 ""
  fin:
  ;$1 contains the folder of Ghostscript or empty
FunctionEnd

Function getSymyxPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\Symyx Technologies, Inc.\Symyx Draw\Client\$2" \ 
     "Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 
  IfFileExists $1\SymyxDraw.exe fin
  StrCpy $1 ""
  fin:
  ;$1 contains the folder of Symyx Draw or empty
FunctionEnd

Function un.getSymyxPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\Symyx Technologies, Inc.\Symyx Draw\Client\$2" \ 
     "Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 
  IfFileExists $1\SymyxDraw.exe fin
  StrCpy $1 ""
  fin:
  ;$1 contains the folder of Symyx Draw or empty
FunctionEnd

Function createOSRAbat
fileOpen $0 "$INSTDIR\osra.bat" w
  fileWrite $0 "\
@echo off$\r$\n\
setlocal$\r$\n\
set OMP_NUM_THREADS=1$\r$\n\
set PATH=%~dp0;$1\bin;$1\lib;%PATH%$\r$\n\
set MAGICK_CONFIGURE_PATH=%~dp0%$\r$\n\
osra.exe %*$\r$\n\
endlocal$\r$\n\
"
fileClose $0
FunctionEnd

Function MakeSureIGotGFL
  strcpy $3 "GflAx_is1"
  strcpy $4 "GflAx\lib\GflAx.dll"
  Call getSoftInstPath
  Pop $0
  StrCmp $0 "" getgfl
  Return
    
  getgfl:
   DetailPrint "need to download and install GflAxSetup.exe"
   Call ConnectInternet ;Make an internet connection (if no connection available)
   StrCpy $2 "$TEMP\GflAxSetup.exe"
   NSISdl::download http://download.xnview.com/GflAxSetup.exe $2
   Pop $0
   StrCmp $0 success success
    SetDetailsView show
    DetailPrint "download failed: $0"
    Abort
  success:
    ExecWait "$2"
    Delete $2
    
FunctionEnd

Function ConnectInternet

  Push $R0
    
    ClearErrors
    Dialer::AttemptConnect
    IfErrors noie3
    
    Pop $R0
    StrCmp $R0 "online" connected
      MessageBox MB_OK|MB_ICONSTOP "Cannot connect to the internet."
      Quit
    
    noie3:
  
    ; IE3 not installed
    MessageBox MB_OK|MB_ICONINFORMATION "Please connect to the internet now."
    
    connected:
  
  Pop $R0
  
FunctionEnd

Function createXMLConfig
fileOpen $0 "$1\AddIns\OSRAAction\OSRAAction.dll.config" w
 fileWrite $0 '<?xml version="1.0"?>$\r$\n'
 fileWrite $0 '<configuration>$\r$\n'
 fileWrite $0 '<configSections>$\r$\n'
 fileWrite $0 '<sectionGroup name="userSettings" type="System.Configuration.UserSettingsGroup, System, Version=2.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089" >$\r$\n'
 fileWrite $0 '<section name="Symyx.Draw.Addins.OSRAAddinAction.OSRA" type="System.Configuration.ClientSettingsSection, System, Version=2.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089" allowExeDefinition="MachineToLocalUser" requirePermission="false" />$\r$\n'
 fileWrite $0 '<section name="Symyx.Draw.Addins.OSRAAddinAction.Settings1" type="System.Configuration.ClientSettingsSection, System, Version=2.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089" allowExeDefinition="MachineToLocalUser" requirePermission="false" />$\r$\n'
 fileWrite $0 '</sectionGroup>$\r$\n'
 fileWrite $0 '</configSections>$\r$\n'
 fileWrite $0 '<startup><supportedRuntime version="v2.0.50727"/></startup><userSettings>$\r$\n'
 fileWrite $0 '<Symyx.Draw.Addins.OSRAAddinAction.OSRA>$\r$\n'
 fileWrite $0 '<setting name="OSRAPath" serializeAs="String">$\r$\n'
 fileWrite $0 '<value>$INSTDIR\osra.bat</value>$\r$\n'
 fileWrite $0 '</setting>$\r$\n'
 fileWrite $0 '<setting name="OSRACommandLine" serializeAs="String">$\r$\n'
 fileWrite $0 '<value> -f sdf </value>$\r$\n'
 fileWrite $0 '</setting>$\r$\n'
 fileWrite $0 '</Symyx.Draw.Addins.OSRAAddinAction.OSRA>$\r$\n'
 fileWrite $0 '<Symyx.Draw.Addins.OSRAAddinAction.Settings1>$\r$\n'
 fileWrite $0 '<setting name="OSRAPath" serializeAs="String">$\r$\n'
 fileWrite $0 '<value />$\r$\n'
 fileWrite $0 '</setting>$\r$\n'
 fileWrite $0 '<setting name="Resolution" serializeAs="String">$\r$\n'
 fileWrite $0 '<value />$\r$\n'
 fileWrite $0 '</setting>$\r$\n'
 fileWrite $0 '</Symyx.Draw.Addins.OSRAAddinAction.Settings1>$\r$\n'
 fileWrite $0 '</userSettings>$\r$\n'
 fileWrite $0 '</configuration>$\r$\n'
fileClose $0
FunctionEnd

Function .onInit
 strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
 call CheckSoftVersion
 strcmp $2 "" no_symyx
 call getSymyxPath
 strcmp $1 "" no_symyx
  SectionGetFlags "${symyx_draw}" $0
  IntOp $0 $0 | ${SF_SELECTED}
  SectionSetFlags "${symyx_draw}" $0
 no_symyx:
FunctionEnd