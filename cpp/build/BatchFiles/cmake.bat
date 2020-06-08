echo off

cd /d %~dp0
cd ../

for %%f in ( * ) do call :rmfiles "%%f"
for /D %%f in ( * ) do call :rmfiles "%%f" d

cmake ../src

msbuild Test_NumCpp.sln /p:Configuration=Release

exit /b

:rmfiles
if %1=="BatchFiles" goto :EOF
if "%2"=="" del %1
if "%2"=="d" rd /S /Q %1
goto :EOF
