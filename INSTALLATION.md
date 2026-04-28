# MetaboGraph Installation & Build Guide

## Quick Start (Running from Source)

```bash
# 1. Navigate to MetaboGraph directory
cd MetaboGraph

# 2. Install all dependencies
pip install -r requirements.txt

# 3. Run the application
python metabograph.py
```

## Building the Windows Executable

### Prerequisites

1. **Python 3.10 or higher** (tested on Python 3.12+)
   ```bash
   python --version  # Should show 3.10+
   ```

2. **Visual C++ Build Tools** (required for building Cython extensions)
   - Download: https://visualstudio.microsoft.com/visual-cpp-build-tools/
   - Or install Visual Studio Community with C++ components

3. **PyInstaller** (included in requirements.txt)
   ```bash
   pip install PyInstaller>=5.0.0
   ```

### Step 1: Install All Dependencies

This is **CRITICAL** - all dependencies must be installed before building:

```bash
cd MetaboGraph
pip install -r requirements.txt
```

**Key packages that must be installed:**
- requests ✓ (for downloading databases)
- pandas ✓
- numpy ✓
- scipy ✓
- scikit-learn ✓
- matplotlib ✓
- seaborn ✓
- networkx ✓
- And all others in requirements.txt

### Step 2: Verify Installation

```bash
# Test that the application runs from source
python metabograph.py

# Accept/close the GUI to verify it launches
```

If you see an error about missing modules, **go back to Step 1** and ensure pip install completes without errors.

### Step 3: Build the Executable

```bash
# Run PyInstaller with the spec file
pyinstaller metabograph.spec --noconfirm --clean

# Building will take 2-5 minutes...
```

**What PyInstaller does:**
- Bundles Python runtime + all installed packages
- Collects required data files
- Creates executable: `dist\MetaboGraph\MetaboliteAnnotationTool.exe`

### Step 4: Distribute the Application

The complete application is in the `dist\MetaboGraph\` folder. You can:

```bash
# Zip the entire folder for distribution
Compress-Archive -Path "dist\MetaboGraph" -DestinationPath "MetaboGraph_installer.zip"

# Or create an installer using Inno Setup, NSIS, etc.
# (Not included in this package)
```

### Step 5: Users Install the Application

1. Extract `MetaboGraph_installer.zip` to desired location
   - Example: `C:\Program Files\MetaboGraph`

2. Run the application:
   - Double-click: `MetaboliteAnnotationTool.exe`
   - Or from command line: `MetaboliteAnnotationTool.exe`

3. On first run, copy database files:
   - Create folder: `Databases\` (if it doesn't exist)
   - Copy `.feather` files to the `Databases\` folder

## Troubleshooting Build Issues

### Error: "ModuleNotFoundError: No module named 'X'"

**Solution:** The package isn't installed
```bash
# Reinstall all requirements
pip install -r requirements.txt --upgrade --force-reinstall

# Or install specific missing package
pip install requests beautifulsoup4 lxml tqdm
```

### Error: "Cannot find Visual C++ build tools"

**Solution:** Install Visual Studio Build Tools
- Download: https://visualstudio.microsoft.com/visual-cpp-build-tools/
- Run installer, select "C++ build tools"
- Restart Python/terminal after installation

### Executable closes immediately

**Solution:** Check the diagnostic log
```
1. The executable creates a console window
2. Look for the error dialog with log file location
3. Log file location: C:\Users\<username>\AppData\Local\Temp\metabograph_startup_*.log
4. Share this log file for debugging
```

### Executable is very large (>500 MB)

**Normal behavior** - PyInstaller bundles the entire Python runtime plus all scientific packages (numpy, scipy, sklearn, matplotlib, etc.)

## Technical Details

### PyInstaller Configuration (metabograph.spec)

The spec file includes:

1. **Hidden Imports** - Packages that need explicit inclusion:
   - All sklearn submodules and Cython extensions
   - scipy.stats, scipy.special, scipy.sparse
   - matplotlib.backends.backend_tkagg
   - plotly, dash, dash_cytoscape
   - And 50+ other scientific computing packages

2. **Data Files** - Non-Python files needed at runtime:
   - `gui/` module files
   - `main_script/` analysis modules
   - Configuration JSON files
   - Database template files

3. **Runtime Hooks** - Pre-initialization scripts:
   - `main_script/runtime_license_hook.py` - License validation (currently disabled)

4. **Excluded Packages** - Not needed:
   - pytest, jupyter, IPython
   - PyQt5/PyQt6, PySide, wx (GUI frameworks we don't use)

### Build Paths

When building on Windows with a path containing spaces:
```
✓ Works: C:\Program Files\MetaboGraph
✓ Works: C:\Users\Name\Documents\MyApps\MetaboGraph
✗ May fail: Paths with non-ASCII characters in username
```

If your username has non-ASCII characters, build to a temporary directory without spaces and move after building.

## Advanced: Custom Builds

### For Development Testing

```bash
# Build with console window visible (for debugging)
# (Already set in metabograph.spec: console=True)
pyinstaller metabograph.spec --console
```

### For Release Distribution

```bash
# Build without console window
# Edit metabograph.spec: change console=True to console=False
# Then rebuild
pyinstaller metabograph.spec --noconfirm --clean
```

### For Minimal Size

```bash
# Remove logging/debug files to reduce size
pyinstaller metabograph.spec --noconfirm --clean --onefile
# Note: --onefile creates single executable but slower startup
```

## Verifying the Build

After build completes, test the executable:

```bash
# Navigate to built executable
cd dist\MetaboGraph

# Run it
MetaboliteAnnotationTool.exe

# Check for startup errors in dialog
# Log file: C:\Users\<username>\AppData\Local\Temp\metabograph_startup_*.log
```

## Diagnosing Installation Issues

If the installed executable closes immediately:

1. **Run from command prompt** (not by double-clicking):
   ```cmd
   cd "C:\Program Files\MetaboGraph"
   MetaboliteAnnotationTool.exe
   ```
   Watch the console for error messages

2. **Check diagnostic log**:
   ```
   Folder: C:\Users\<username>\AppData\Local\Temp\
   File: metabograph_startup_*.log (most recent)
   ```
   Open with any text editor and look for errors

3. **Common issues**:
   - Missing database files in `Databases\` folder
   - Corrupt installation - try reinstalling
   - Antivirus blocking execution - whitelist the application
   - Missing Visual C++ runtime - run Windows Update

## Support

For build/installation issues:

1. Share the **diagnostic log file** from `%TEMP%\metabograph_startup_*.log`
2. Include Python version: `python --version`
3. Include OS: Windows 10/11/Server version
4. Steps to reproduce the issue
5. Full error message from console or error dialog

See [DIAGNOSTICS.md](DIAGNOSTICS.md) for more logging information.
