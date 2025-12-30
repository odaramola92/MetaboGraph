@echo off
REM MetaboGraph Launcher for Windows
REM This batch file launches the MetaboGraph GUI application

echo Starting MetaboGraph...
echo.

REM Run the GUI
python "%~dp0main_script\run_gui.py"

REM Pause if there was an error
if errorlevel 1 (
    echo.
    echo An error occurred. Press any key to exit...
    pause > nul
)
