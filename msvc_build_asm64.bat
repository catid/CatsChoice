@ECHO OFF

ECHO ------------------------------------------------------------------
ECHO Assumes FASM is located at C:\FASM\FASM.EXE or in the system path.
ECHO Download FASM from http://flatassembler.net/

PATH=%PATH%;C:\FASM

ECHO Assembling 64-bit Pseudo-Mersenne library...
ECHO ------------------------------------------------------------------

FASM.EXE "big_x64_mscoff.asm" "big_x64.obj"

ECHO ------------------------------------------------------------------

ECHO Now you need to build the .lib file.  Open up the MSVC x64 Command
ECHO Prompt, change to this directory and type:
ECHO    lib /out:big_x64.lib big_x64.obj
ECHO This will produce the big_x64.lib library that MSVC can link to.

PAUSE
