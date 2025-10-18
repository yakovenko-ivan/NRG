#!/bin/bash
del dir.txt
pushd %~dp0
dir /b *.plt /s 2> nul | find "" /v /c > tmp 
set /p count=<tmp 
del tmp
cd ..
echo %count% >> dir.txt
