@echo off
cls

.paket\paket.bootstrapper.exe prerelease
if errorlevel 1 (
  exit /b %errorlevel%
)

.paket\paket.exe  install --hard
if errorlevel 1 (
  exit /b %errorlevel%
)

MSBuild SequenceAlignment.sln /p:Configuration=Release