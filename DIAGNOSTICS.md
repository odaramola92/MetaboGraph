# MetaboGraph Diagnostics Guide

## Overview

MetaboGraph has been enhanced with comprehensive diagnostic logging to help identify startup issues, especially in the frozen PyInstaller executable.

## Diagnostic Features

### 1. **Startup Script Diagnostics** (`metabograph.py`)

The launcher script now:
- ✓ Logs all startup steps to a file in `%TEMP%\metabograph_startup_YYYYMMDD_HHMMSS.log`
- ✓ Prints messages to console with timestamps
- ✓ Tests critical imports (tkinter, pandas, numpy) before loading GUI
- ✓ Provides full traceback with file path to log in error dialog
- ✓ Works in both normal and frozen (PyInstaller) modes

**What it logs:**
```
[timestamp] Python version: 3.12.0 (tags/v3.12.0:...)
[timestamp] Working directory: C:\path\to\MetaboGraph
[timestamp] Script directory: C:\path\to\MetaboGraph
[timestamp] Testing Python standard library imports...
[timestamp] ✓ tkinter imported successfully
[timestamp] ✓ pandas imported successfully
[timestamp] ✓ numpy imported successfully
[timestamp] Importing MetaboGraph GUI main module...
[timestamp] ✓ gui.main imported successfully
[timestamp] Launching main application...
```

### 2. **GUI Main Module Diagnostics** (`gui/main.py`)

Enhanced logging tracks:
- ✓ File-based logging to `%TEMP%\metabograph_gui_YYYYMMDD_HHMMSS.log`
- ✓ Root window creation
- ✓ ttk styles setup
- ✓ DataManager initialization
- ✓ Title bar creation
- ✓ Notebook (tab container) setup
- ✓ Each individual tab loading with detailed error reporting

**What it logs:**
```
[HH:MM:SS] ================================================================================
[HH:MM:SS] MetaboGraph Main Application Entry Point
[HH:MM:SS] ================================================================================
[HH:MM:SS] Creating MetaboGraphGUI instance...
[HH:MM:SS] ================================================================================
[HH:MM:SS] Initializing MetaboGraph GUI
[HH:MM:SS] ================================================================================
[HH:MM:SS] Creating root Tk window...
[HH:MM:SS] ✓ Root window created successfully
[HH:MM:SS] Root window configured
[HH:MM:SS] Setting up ttk styles...
[HH:MM:SS] ✓ ttk styles configured
[HH:MM:SS] Initializing DataManager...
[HH:MM:SS] ✓ DataManager initialized successfully
[HH:MM:SS] Setting up tabs...
[HH:MM:SS] Loading tab: 📊 Data Cleaning
[HH:MM:SS]   → Instantiating DataCleaningTab...
[HH:MM:SS]   ✓ DataCleaningTab instantiated
[HH:MM:SS]   → Getting frame from tab...
[HH:MM:SS]   ✓ Frame retrieved
[HH:MM:SS]   → Adding frame to notebook...
[HH:MM:SS]   ✓ Frame added to notebook
[HH:MM:SS] ✓ Successfully loaded tab: 📊 Data Cleaning
```

### 3. **Per-Tab Error Reporting**

If a tab fails to load:
```
[HH:MM:SS] Loading tab: 🧬 Pathway Analysis
[HH:MM:SS]   → Instantiating PathwayAnalysisParentTab...
[HH:MM:SS] ✗ Error loading tab '🧬 Pathway Analysis': [specific error message]
[HH:MM:SS] Traceback (most recent call last):
[HH:MM:SS]   File "gui/main.py", line XXX, in _add_tab
[HH:MM:SS]     tab_instance = TabClass(self.notebook, self.data_manager)
[HH:MM:SS] [Full traceback...]
[HH:MM:SS]   ✓ Error tab placeholder added
```

## How to Use Diagnostics

### When Application Closes Immediately

1. **Run the application from command line:**
   ```cmd
   cd C:\path\to\MetaboGraph
   python metabograph.py
   ```
   
   Watch the console output for error messages.

2. **Check log files:**
   - Startup log: `%TEMP%\metabograph_startup_*.log`
   - GUI log: `%TEMP%\metabograph_gui_*.log`
   
   These files contain full diagnostic output including timestamps and tracebacks.

3. **Error Dialog**
   If an error dialog appears, it will show:
   - Brief error message
   - Path to the diagnostic log file
   
   Example:
   ```
   MetaboGraph Failed to Start
   
   ImportError: cannot import name 'xyz' from 'module' (...)
   
   Diagnostic log:
   C:\Users\username\AppData\Local\Temp\metabograph_startup_20240115_143022.log
   ```

### Common Issues and Solutions

#### Console Closes Too Fast to Read

- Run from Command Prompt/PowerShell instead of double-clicking the EXE
- The diagnostic log will still be saved even if console closes
- Check the log file: `%TEMP%\metabograph_startup_*.log`

#### Error Dialog Shows "Check Console Output"

- Log file path will be displayed in error dialog
- Open the diagnostic log file with any text editor
- Search for `EXCEPTION` or `Error` keywords for the root cause

#### "Tab Failed to Load"

- Check the GUI log file
- Find the failing tab name (e.g., "🧬 Pathway Analysis")
- The log will show which specific import or initialization failed
- Report the full traceback section to developers

## Interpreting Log Files

### Log File Format

Each line has this format:
```
[HH:MM:SS.mmm] LEVEL: MESSAGE
```

- `HH:MM:SS.mmm` = timestamp (hour:minute:second.milliseconds)
- `LEVEL` = DEBUG, INFO, WARNING, ERROR, CRITICAL
- `MESSAGE` = diagnostic message

### Key Log Markers

- ✓ = Successful step
- ✗ = Failed step  
- → = Starting a substep
- `="*80` = Section boundary (major phase)

### Search Terms for Troubleshooting

| Issue | Search For |
|-------|-----------|
| Import errors | "ImportError", "ModuleNotFoundError" |
| Missing files | "FileNotFoundError", "not found" |
| GUI widget creation | "_setup_", "Frame", "Notebook" |
| Specific tab failure | tab name (e.g., "Pathway Analysis") |
| Full error trace | "Traceback", "EXCEPTION" |

## Development Notes

### Adding Diagnostics to New Code

Use the standard logger:

```python
import logging

logger = logging.getLogger(__name__)

# In __init__ or startup function:
logger.info("Starting critical operation...")
logger.debug("  → Substep: loading data...")
logger.debug("  ✓ Data loaded (5,000 rows)")

try:
    # Do something
except Exception as e:
    logger.error(f"Operation failed: {e}", exc_info=True)
    raise
```

### Log Levels

- **DEBUG**: Low-level details (frame created, variable values)
- **INFO**: Major progress checkpoints (tab loaded, operation complete)
- **WARNING**: Non-fatal issues (deprecated function, fallback used)
- **ERROR**: Failed operations (tab load failed, file not found)
- **CRITICAL**: Fatal errors (will crash application)

## PyInstaller Executable

When building with `pyinstaller metabograph.spec`:

1. The `console=True` setting ensures a console window appears
2. Diagnostics are written to both console and log file
3. If the executable crashes silently:
   - Still check the log file in `%TEMP%`
   - May indicate very early import error before GUI creation
   - Try running from command line to see console output

### Building with Enhanced Diagnostics

```cmd
# Standard build
pyinstaller metabograph.spec

# Debug build (adds more logging)
pyinstaller metabograph.spec --debug all
```

## Support

When reporting a crash:

1. **Provide the full log file contents** from:
   - `%TEMP%\metabograph_startup_*.log` (most recent)
   - `%TEMP%\metabograph_gui_*.log` (most recent)

2. **Include any error dialog text** that was shown

3. **Note the exact sequence**: 
   - Click "start app" → app closes (how many seconds?)
   - Does error dialog appear?
   - Is there a console window?

4. **Python and system info**:
   - Python version: `python --version`
   - OS: Windows 10/11?
   - Antivirus software might be blocking executable

This information helps developers quickly identify and fix startup issues.
