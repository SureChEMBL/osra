@echo off
setlocal
set PATH=%~dp0;%ProgramFiles(x86)%\gs\gs8.61\bin;%ProgramFiles(x86)%\gs\gs8.61\lib;%ProgramFiles%\gs\gs8.61\bin;%ProgramFiles%\gs\gs8.61\lib;%PATH%
set MAGICK_CONFIGURE_PATH=%~dp0%
osra.exe %*
endlocal