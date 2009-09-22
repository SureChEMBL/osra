@echo off
setlocal
set current=%CD%
cd %~dp0%
set OMP_NUM_THREADS=1
set PATH=%~dp0;%ProgramFiles(x86)%\gs\gs8.63\bin;%ProgramFiles(x86)%\gs\gs8.63\lib;%ProgramFiles%\gs\gs8.63\bin;%ProgramFiles%\gs\gs8.63\lib;%PATH%
set MAGICK_CONFIGURE_PATH=%~dp0%
osra.exe "%current%\\%*"
cd %current%
endlocal
