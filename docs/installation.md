# Installation Guide

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation Methods](#installation-methods)
  - [Option 1: Standalone Executable (Recommended)](#option-1-standalone-executable-recommended)
  - [Option 2: From Source](#option-2-from-source)
- [Database Setup](#database-setup)
- [Verification](#verification)
- [Troubleshooting Installation](#troubleshooting-installation)

---

## System Requirements

### Minimum Requirements
- **Operating System:** Windows 10/11 (64-bit)
- **RAM:** 8 GB minimum, 16 GB recommended
- **Storage:** 10 GB free space (for application + databases)
- **Display:** 1280x720 minimum resolution

### For Source Installation
- **Python:** 3.8 or higher (developed and tested on Python 3.13)
- **pip:** Latest version recommended

---

## Installation Methods

### Option 1: Standalone Executable (Recommended)

**Best for:** Users who want to use MetaboGraph without installing Python.

#### Step 1: Download
Download the latest release from:
- **Zenodo:** [Insert DOI link]
- **OneDrive:** [Insert download link]
- **GitHub Releases:** [Insert GitHub release link]

#### Step 2: Extract
1. Extract `MetaboGraph_v1.0.0.zip` to your desired location
2. Recommended: `C:\Program Files\MetaboGraph\` or `C:\Users\YourName\MetaboGraph\`

#### Step 3: Run
1. Navigate to the extracted folder
2. Double-click `MetaboliteAnnotationTool.exe`
3. First launch may take 10-15 seconds

#### Step 4: Database Setup
See [Database Setup](#database-setup) section below.

---

### Option 2: From Source

**Best for:** Developers or users who want to modify the code.

#### Step 1: Clone Repository
```bash
git clone https://github.com/yourusername/MetaboGraph.git
cd MetaboGraph
```

#### Step 2: Create Virtual Environment (Recommended)
```bash
# Windows
python -m venv venv
venv\Scripts\activate

# macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

#### Step 3: Install Dependencies
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

#### Step 4: Verify Installation
```bash
python metabograph.py
```

#### Step 5: Database Setup
See [Database Setup](#database-setup) section below.

---

## Database Setup

MetaboGraph requires several databases for annotation and pathway analysis. Databases are **not included** in the installation due to their size (~2-3 GB total).

### Automatic Setup (Recommended)

1. **Open MetaboGraph**
2. **Go to "Database Setup" tab**
3. **Click "Download All Databases"**
4. **Wait for processing** (15-30 minutes depending on internet speed)

The application will:
- Download databases from official sources
- Process and optimize for fast lookup
- Save to `Databases/` folder
- Verify integrity

### Manual Setup

If automatic download fails, you can download databases manually:

#### Required Databases

| Database | Source | Purpose |
|----------|--------|---------|
| **HMDB** | [hmdb.ca](https://hmdb.ca/downloads) | Metabolite identification |
| **LipidMaps** | [lipidmaps.org](https://www.lipidmaps.org/resources/downloads) | Lipid identification |
| **PathBank** | [pathbank.org](https://pathbank.org/downloads) | Pathway annotation |
| **SMPDB** | [smpdb.ca](https://smpdb.ca/downloads) | Small molecule pathways |
| **WikiPathways** | [wikipathways.org](https://www.wikipathways.org/index.php/Download_Pathways) | Curated pathways |

#### Manual Processing Steps

1. Download database files from sources above
2. Place raw files in `Databases/raw/` folder
3. Open MetaboGraph → Database Setup tab
4. Click "Process Local Databases"
5. Select organism (Human/Mouse/Rat) for WikiPathways
6. Wait for processing to complete

### Verify Database Installation

In the Database Setup tab, you should see:
- ✅ Green checkmarks for all databases
- File sizes displayed
- "Last Updated" dates

---

## Verification

### Test Your Installation

1. **Launch MetaboGraph**
   - Executable: Double-click `MetaboliteAnnotationTool.exe`
   - Source: Run `python metabograph.py`

2. **Check All Tabs Open:**
   - Data Cleaning
   - ID Annotation
   - Pathway Annotation
   - Pathway Network
   - Comparative Analysis
   - Multi-Omics Integration
   - Database Setup
   - Help

3. **Test Database Connection:**
   - Go to Database Setup tab
   - Verify all databases show green checkmarks

4. **Run Example Workflow (Optional):**
   - Load sample data from `examples/` folder
   - Try basic annotation in ID Annotation tab
   - Verify results appear

---

## Troubleshooting Installation

### Common Issues

#### "Python not found" (Source Installation)
**Solution:**
- Install Python 3.8 or higher from [python.org](https://www.python.org/downloads/)
- During installation, check "Add Python to PATH"
- Restart terminal/command prompt

#### "No module named tkinter"
**Solution:**
```bash
# Windows: Reinstall Python with tkinter
# Linux (Ubuntu/Debian):
sudo apt-get install python3-tk

# macOS:
brew install python-tk
```

#### "Application won't start" (Executable)
**Solution:**
- Run as Administrator (right-click → Run as administrator)
- Check Windows Defender/Antivirus isn't blocking
- Verify all files extracted properly
- Check `logs/` folder for error messages

#### "Database download fails"
**Solution:**
- Check internet connection
- Try manual database download
- Disable VPN/proxy temporarily
- Check firewall settings
- See [Troubleshooting](troubleshooting.md) for detailed steps

#### "Out of memory" during database processing
**Solution:**
- Close other applications
- Increase available RAM
- Process databases one at a time
- Use manual download + processing

#### "Missing DLL" errors (Windows)
**Solution:**
- Install [Microsoft Visual C++ Redistributable](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist)
- Update Windows to latest version

---

## Post-Installation

### Recommended Next Steps

1. **Read Documentation:**
   - [User Guide](user-guides/) - How to use each feature
   - [FAQ](faq.md) - Common questions
   - [Troubleshooting](troubleshooting.md) - Fix common issues

2. **Try Example Data:**
   - Load files from `examples/` folder
   - Follow example workflows
   - Verify outputs match expected results

3. **Configure Settings:**
   - Set default output directories
   - Configure database preferences
   - Adjust analysis parameters

4. **Join Community:**
   - Report issues on GitHub
   - Share feedback
   - Request features

---

## Updating MetaboGraph

### Update Standalone Executable
1. Download latest version
2. Extract to new folder (keep old version as backup)
3. Copy `Databases/` folder to new location (avoid re-downloading)
4. Test new version before deleting old

### Update Source Installation
```bash
cd MetaboGraph
git pull origin main
pip install -r requirements.txt --upgrade
```

---

## Uninstallation

### Remove Executable
1. Delete MetaboGraph folder
2. Delete any desktop shortcuts
3. Optional: Remove `Databases/` if no longer needed

### Remove Source Installation
```bash
# Deactivate virtual environment
deactivate

# Remove directory
cd ..
rm -rf MetaboGraph  # Linux/macOS
# Or delete folder manually on Windows
```

---

## Getting Help

- **Documentation:** [docs/](../docs/)
- **FAQ:** [faq.md](faq.md)
- **Troubleshooting:** [troubleshooting.md](troubleshooting.md)
- **GitHub Issues:** [Report a bug](https://github.com/yourusername/MetaboGraph/issues)
- **Email Support:** support@metabograph.org

---

## License

MetaboGraph is licensed under the MIT License. See [LICENSE](../LICENSE) for details.

Databases used by MetaboGraph are subject to their respective licenses. See [Database Setup](user-guides/04-database-setup.md) for details.
