!include Sections.nsh

; The name of the installer
Name "Optical Structure Recognition Application"

; The file to write
OutFile "osra-setup-1-3-0.exe"

; The default installation directory
InstallDir $PROGRAMFILES\osra\1.3.0

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\osra\1.3.0" "Install_Dir"

LicenseData "license.txt"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

;--------------------------------

; Pages

Page license
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
  File "spelling.txt"
  File "superatom.txt"
  
  strcpy $3 "GPL Ghostscript"
  call checkSoftVersion
  call getGhostscriptInstPath
  strcmp $1 "" no_gs +2
  no_gs:
  MessageBox MB_OK "Ghostscript interpreter not found" IDOK 0 
  call createOSRAbat
  
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\osra\1.3.0 "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "DisplayName" "OSRA 1.3.0"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
  
SectionEnd

Section /o "Symyx Draw plugin" symyx_draw
 strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
 call CheckSoftVersion
 strcmp $2 "" no_symyx
 call getSymyxPath
 strcmp $1 "" no_symyx
 SetOutPath "$1\AddIns"
 File "plugins\symyx_draw\OSRAAction.xml"
 SetOutPath "$1\AddIns\OSRAAction"
 File "plugins\symyx_draw\README.txt"
 File "plugins\symyx_draw\OSRAAction\OSRAAction.dll"
 File "plugins\symyx_draw\OSRAAction\OSRAAction.dll.config"
 Goto done
 no_symyx:
  MessageBox MB_OK "Symyx Draw not found" IDOK done
 done:
SectionEnd

Section /o "ChemBioOffice 12 plugin" chemoffice
 strcpy $3 "CambridgeSoft\ChemScript"
 call CheckSoftVersion
 strcmp $2 "12.0" +1 no_chemoffice
 call getChemScriptPath
 strcmp $1 "" no_chemoffice
 ReadRegStr $2 HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\PIL-py2.5" "DisplayName"
 strcmp $2 "" +1 pil_exists
 call downloadPIL
 pil_exists:
 SetOutPath "$1\Scripts"
 File "Import Structures with OSRA.py"
 Goto done
 no_chemoffice:
  MessageBox MB_OK "ChemScript 12.0 not found" IDOK done
 done:
SectionEnd

; Uninstaller

Section "Uninstall"
  ReadRegStr $0 HKLM SOFTWARE\osra\1.3.0 "Install_Dir"
  strcpy $INSTDIR $0
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra"
  DeleteRegKey HKLM SOFTWARE\osra\1.3.0

  ; Remove files and uninstaller
  Delete $INSTDIR\osra.exe
  Delete $INSTDIR\delegates.xml
  Delete $INSTDIR\README.txt
  Delete $INSTDIR\osra.bat
  Delete $INSTDIR\libopenbabel-3.dll
  Delete $INSTDIR\superatom.txt
  Delete $INSTDIR\spelling.txt
  Delete $INSTDIR\uninstall.exe
  RMDir "$INSTDIR"
  strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
  call un.CheckSoftVersion
  strcmp $2 "" no_symyx
  call un.getSymyxPath
  strcmp $1 "" no_symyx 
  Delete "$1\AddIns\OSRAAction.xml"
  Delete "$1\AddIns\OSRAAction\README.txt"
  Delete "$1\AddIns\OSRAAction\OSRAAction.dll"
  Delete "$1\AddIns\OSRAAction\OSRAAction.dll.config"
  RMDir "$1\AddIns\OSRAAction"
  no_symyx:
  strcpy $3 "CambridgeSoft\ChemScript"
  call un.CheckSoftVersion
  strcmp $2 "12.0" +1 no_chemoffice
  call un.getChemScriptPath
  strcmp $1 "" no_chemoffice
  Delete "$1\Scripts\Import Structures with OSRA.py"
  no_chemoffice:
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

Function downloadPIL
   DetailPrint "need to download and install Python Imaging Library"
   Call ConnectInternet ;Make an internet connection (if no connection available)
   StrCpy $2 "$TEMP\PIL-1.1.6.win32-py2.5.exe"
   NSISdl::download http://effbot.org/media/downloads/PIL-1.1.6.win32-py2.5.exe $2
   Pop $0
   StrCmp $0 success success
    SetDetailsView show
    DetailPrint "download failed: $0"
    Abort
   success:
    ExecWait "$2"
    Delete $2
FunctionEnd

Function getGhostscriptInstPath
 strcmp $2 "" download_gs get_path
 download_gs:
   DetailPrint "need to download and install Ghostscript"
   Call ConnectInternet ;Make an internet connection (if no connection available)
   StrCpy $2 "$TEMP\gs870w32.exe"
   NSISdl::download http://voxel.dl.sourceforge.net/sourceforge/ghostscript/gs870w32.exe $2
   Pop $0
   StrCmp $0 success success
    SetDetailsView show
    DetailPrint "download failed: $0"
    Abort
   success:
    ExecWait "$2"
    Delete $2
    strcpy $2 "8.70"
	
 get_path:
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

Function getChemScriptPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\CambridgeSoft\ChemDraw\12.0\General" \ 
     "ChemDraw Items Default Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 -15
  IfFileExists "$1\Scripts\Get 3D Structure.py" fin
  StrCpy $1 ""
  fin:
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

Function un.getChemScriptPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\CambridgeSoft\ChemDraw\12.0\General" \ 
     "ChemDraw Items Default Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 -15
  IfFileExists "$1\Scripts\Get 3D Structure.py" fin
  StrCpy $1 ""
  fin:
FunctionEnd

Function createOSRAbat
fileOpen $0 "$INSTDIR\osra.bat" w
  fileWrite $0 '\
@echo off$\r$\n\
setlocal$\r$\n\
set exec_dir=%~dp0%$\r$\n\
set OMP_NUM_THREADS=1$\r$\n\
set PATH=%exec_dir%;$1\bin;$1\lib;%PATH%$\r$\n\
set MAGICK_CONFIGURE_PATH=%exec_dir%$\r$\n\
"%exec_dir%osra.exe" %*$\r$\n\
endlocal$\r$\n\
'
fileClose $0
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


Function .onInit
 strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
 call CheckSoftVersion
 strcmp $2 "" no_symyx
 call getSymyxPath
 strcmp $1 "" no_symyx
 SectionGetFlags "${symyx_draw}" $0
 IntOp $0 $0 | ${SF_SELECTED}
 SectionSetFlags "${symyx_draw}" $0
 strcpy $3 "CambridgeSoft\ChemScript"
 call CheckSoftVersion
 strcmp $2 "12.0" +1 no_chemoffice
 call getChemScriptPath
 strcmp $1 "" no_chemoffice
 SectionGetFlags "${chemoffice}" $0
 IntOp $0 $0 | ${SF_SELECTED}
 SectionSetFlags "${chemoffice}" $0
 no_symyx:
 no_chemoffice:
FunctionEnd
