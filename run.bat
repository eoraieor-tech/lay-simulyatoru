@echo off
cd /d "%~dp0"

if exist "..\venv\Scripts\activate.bat" (
    call "..\venv\Scripts\activate.bat"
) else if exist "venv\Scripts\activate.bat" (
    call "venv\Scripts\activate.bat"
) else (
    echo.
    echo XETA: venv tapilmadi.
    echo Gozlenilen yerler:
    echo   %~dp0..\venv
    echo   %~dp0venv
    echo.
    pause
    exit /b 1
)

python app.py
if errorlevel 1 pause