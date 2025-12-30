# Troubleshooting Guide

## Table of Contents
- [General Issues](#general-issues)
- [Installation Problems](#installation-problems)
- [Database Issues](#database-issues)
- [Data Import/Export Issues](#data-importexport-issues)
- [Analysis Errors](#analysis-errors)
- [Performance Issues](#performance-issues)
- [GUI Problems](#gui-problems)
- [Getting Additional Help](#getting-additional-help)

---

## General Issues

### Application Won't Start

**Symptoms:** Double-clicking executable does nothing, or window closes immediately.

**Solutions:**

1. **Run as Administrator:**
   - Right-click `MetaboliteAnnotationTool.exe`
   - Select "Run as administrator"

2. **Check Logs:**
   - Navigate to `logs/` folder
   - Open latest log file
   - Look for error messages at the end
   - Share log contents when reporting issues

3. **Verify Complete Installation:**
   - Ensure all files extracted properly
   - Check `_internal/` folder exists (for PyInstaller version)
   - Total folder size should be ~200-250 MB

4. **Antivirus/Firewall:**
   - Some antivirus software blocks unsigned executables
   - Add MetaboGraph folder to exclusion list
   - Temporarily disable antivirus to test

5. **Missing Dependencies (Windows):**
   - Install [Microsoft Visual C++ Redistributable](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist)
   - Install [.NET Framework 4.8](https://dotnet.microsoft.com/download/dotnet-framework/net48)

---

### Application Crashes During Use

**Symptoms:** Application closes unexpectedly, freezes, or shows error dialog.

**Solutions:**

1. **Check Available Memory:**
   - Open Task Manager (Ctrl+Shift+Esc)
   - Check RAM usage
   - Close other applications
   - Minimum 8 GB RAM recommended

2. **Check Disk Space:**
   - Ensure at least 5 GB free space
   - Large datasets require temporary storage

3. **Review Log Files:**
   - Check `logs/` folder for stack traces
   - Look for "MemoryError", "PermissionError", etc.

4. **Try Smaller Dataset:**
   - Test with example data first
   - If works, issue may be dataset-specific

5. **Update Application:**
   - Download latest version
   - Bug may already be fixed

---

## Installation Problems

### "Python not found" (Source Installation)

**Solution:**
```bash
# Verify Python installation
python --version

# Should show Python 3.8+
# If not found:
# 1. Download from python.org
# 2. During installation, check "Add Python to PATH"
# 3. Restart terminal
```

### "No module named tkinter"

**Symptoms:** Error when launching from source.

**Solutions:**

**Windows:**
- Tkinter should be included with Python
- Reinstall Python, ensure "tcl/tk and IDLE" is checked

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install python3-tk
```

**macOS:**
```bash
brew install python-tk
```

### "pip install" Fails

**Common Errors:**

1. **"Could not find a version that satisfies the requirement..."**
   ```bash
   # Update pip
   python -m pip install --upgrade pip
   
   # Try again
   pip install -r requirements.txt
   ```

2. **"Permission denied"**
   ```bash
   # Use --user flag
   pip install -r requirements.txt --user
   
   # Or use virtual environment (recommended)
   python -m venv venv
   venv\Scripts\activate
   pip install -r requirements.txt
   ```

3. **"Microsoft Visual C++ 14.0 is required"**
   - Install [Microsoft C++ Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
   - Or install pre-compiled wheels

---

## Database Issues

### Databases Won't Download

**Symptoms:** Download button does nothing, or download fails partway.

**Solutions:**

1. **Check Internet Connection:**
   - Ensure stable connection
   - Try different network (some institutions block large downloads)
   - Disable VPN temporarily

2. **Manual Download:**
   - Download databases manually from source websites
   - Place in `Databases/raw/` folder
   - Use "Process Local Databases" button

3. **Firewall/Proxy:**
   - Configure proxy settings if behind corporate firewall
   - Add database URLs to firewall whitelist

4. **Disk Space:**
   - Ensure 5+ GB free space
   - Databases require ~3 GB total

### "Database file not found" Error

**Symptoms:** Application says database missing even after download.

**Solutions:**

1. **Verify Database Files Exist:**
   ```
   Databases/
   ├── hmdb_database.feather
   ├── lipidmap.feather
   ├── pathbank_selected.feather
   ├── merged_SMP_metabolites.feather
   ├── wikipathways_homo_sapiens.feather
   └── ... (other files)
   ```

2. **Check File Permissions:**
   - Right-click database files → Properties
   - Ensure not "Read-only"
   - Ensure your user has read permissions

3. **Reprocess Databases:**
   - Delete corrupted `.feather` files
   - Re-run database processing
   - Check logs for processing errors

### Database Processing Stuck/Slow

**Symptoms:** Processing bar stuck at certain percentage, or takes hours.

**Solutions:**

1. **Be Patient:**
   - First-time processing takes 15-30 minutes
   - Large databases (HMDB, KEGG) take longest

2. **Check Activity:**
   - Open Task Manager
   - Verify MetaboGraph is using CPU (processing active)
   - If CPU idle for >5 minutes, likely hung

3. **Increase Memory:**
   - Close all other applications
   - Restart computer to free up RAM

4. **Process One at a Time:**
   - Instead of "Download All", download individually
   - Process each database separately

---

## Data Import/Export Issues

### "Could not read file" Error

**Common Causes:**

1. **File Format Issues:**
   - Ensure file is `.xlsx` (not `.xls` or `.csv`)
   - Open in Excel and re-save as `.xlsx`
   - Check file isn't corrupted

2. **File Path Issues:**
   - Avoid special characters in filename: `#`, `%`, `&`
   - Avoid very long file paths (>256 characters)
   - Move file to simpler location (e.g., `C:\Data\`)

3. **File Locked:**
   - Close file in Excel
   - Check if another program is using it
   - Restart computer if necessary

4. **Encoding Issues:**
   - Ensure Excel file uses UTF-8 encoding
   - Re-save file with "Save As" → `.xlsx`

### Export Fails or Creates Empty File

**Solutions:**

1. **Check Permissions:**
   - Ensure you have write permissions to output folder
   - Try exporting to Desktop first
   - Run as Administrator if needed

2. **Check Disk Space:**
   - Large exports require significant space
   - Ensure 1+ GB free on target drive

3. **Close Output File:**
   - If file already exists and open in Excel, close it
   - Application can't overwrite open files

4. **Verify Data:**
   - Check that analysis produced results
   - Empty export may indicate no data to export

---

## Analysis Errors

### ID Annotation Returns No Results

**Possible Causes:**

1. **Database Not Loaded:**
   - Go to Database Setup tab
   - Verify all databases have green checkmarks
   - Re-download if needed

2. **Column Names Incorrect:**
   - Check input file has required columns:
     - `Compounds` or `Name`
     - `RT [min]` or `Retention Time`
     - `m/z` or `Mass`
   - Column names must match exactly

3. **API Rate Limiting:**
   - PubChem/KEGG may limit requests
   - Wait 5-10 minutes and try again
   - Use smaller dataset for testing

4. **Poor Data Quality:**
   - Check mass values are reasonable (50-2000 m/z)
   - Check retention times are positive
   - Remove duplicate rows

### Pathway Annotation Shows "No Pathways Found"

**Solutions:**

1. **Verify IDs Present:**
   - Check input file has HMDB, KEGG, or PubChem IDs
   - At least one ID column required
   - IDs must be in correct format (e.g., `HMDB0000001`)

2. **Check Pathway Database:**
   - Go to Database Setup
   - Verify PathBank or SMPDB loaded
   - Reprocess if needed

3. **Metabolite Not in Pathways:**
   - Some metabolites have IDs but no pathway annotations
   - This is expected for some compounds
   - Check with known metabolite (e.g., glucose)

### Network Generation Fails

**Symptoms:** "Could not generate network" error.

**Solutions:**

1. **Insufficient Data:**
   - Network requires at least 3-5 pathways
   - Check pathway annotation produced results

2. **Missing P-values:**
   - Network tab requires enrichment analysis results
   - Ensure Fisher's exact test completed
   - Check "Selected_Pathways" sheet exists

3. **Memory Issues:**
   - Large networks (>100 pathways) require significant RAM
   - Filter to top 50 pathways
   - Close other applications

---

## Performance Issues

### Application Runs Slowly

**Solutions:**

1. **Reduce Dataset Size:**
   - Start with subset of data
   - Filter before analysis
   - Process in batches

2. **Close Unused Tabs:**
   - Each open plot/visualization uses memory
   - Close completed analyses

3. **Optimize Database:**
   - Databases are already optimized (`.feather` format)
   - If slow, regenerate databases

4. **System Resources:**
   - Check Task Manager for memory usage
   - Close background applications
   - Restart computer

### Large File Export Takes Forever

**Solutions:**

1. **Use Filtered Data:**
   - Export only significant results
   - Filter before exporting

2. **Choose Efficient Formats:**
   - CSV faster than Excel for large datasets
   - Use "Export Summary" instead of full data

3. **Disable Formatting:**
   - Skip conditional formatting
   - Use plain export

---

## GUI Problems

### Buttons Don't Work

**Solutions:**

1. **Wait for Completion:**
   - Some operations take time
   - Check if progress bar is active
   - Don't click multiple times

2. **Check Error Messages:**
   - Look at bottom status bar
   - Red text indicates errors

3. **Restart Application:**
   - Close and reopen
   - Check if issue persists

### Text Too Small/Large

**Solutions:**

1. **Adjust Windows Display Settings:**
   - Right-click Desktop → Display Settings
   - Change scaling (100%, 125%, 150%)
   - Restart application

2. **Change Resolution:**
   - Higher resolution = more space but smaller text
   - Minimum 1280x720 recommended

### Plots Don't Display

**Solutions:**

1. **Wait for Generation:**
   - Large plots take time to render
   - Check progress in status bar

2. **Check Data:**
   - Verify input data has values to plot
   - Check for missing or NaN values

3. **Update Graphics Drivers:**
   - Especially for network visualizations
   - Update from manufacturer website

---

## Getting Additional Help

### Before Reporting Issues

**Gather This Information:**

1. **Version Information:**
   - MetaboGraph version
   - Python version (if source installation)
   - Operating system version

2. **Log Files:**
   - Latest log from `logs/` folder
   - Include full error message and stack trace

3. **Steps to Reproduce:**
   - What you did step-by-step
   - What you expected to happen
   - What actually happened

4. **Sample Data:**
   - If possible, provide minimal dataset that reproduces issue
   - Remove sensitive/proprietary data

### Report Issues

**GitHub Issues:** [https://github.com/yourusername/MetaboGraph/issues](https://github.com/yourusername/MetaboGraph/issues)

**Email Support:** support@metabograph.org

**Include:**
- Problem description
- Error messages/logs
- System information
- Steps to reproduce

### Community Resources

- **Documentation:** [docs/](../docs/)
- **FAQ:** [faq.md](faq.md)
- **Installation Guide:** [installation.md](installation.md)
- **User Guides:** [user-guides/](user-guides/)

---

## Common Error Messages

### "ImportError: No module named 'X'"

**Meaning:** Required Python package not installed.

**Solution:**
```bash
pip install X
# Or reinstall all dependencies:
pip install -r requirements.txt
```

### "PermissionError: [Errno 13] Permission denied"

**Meaning:** Application can't write to file/folder.

**Solutions:**
- Run as Administrator
- Check file isn't open in another program
- Check folder permissions
- Change output location

### "MemoryError"

**Meaning:** Not enough RAM for operation.

**Solutions:**
- Close other applications
- Use smaller dataset
- Add more RAM to computer
- Process in batches

### "KeyError: 'column_name'"

**Meaning:** Expected column not found in data.

**Solutions:**
- Check column names in input file
- Ensure exact spelling and capitalization
- Check for extra spaces
- See user guide for required columns

### "ValueError: could not convert string to float"

**Meaning:** Non-numeric data in numeric column.

**Solutions:**
- Check for text in numeric columns
- Remove empty rows
- Check for special characters
- Validate data in Excel first

---

## Still Having Issues?

If you've tried these solutions and still have problems:

1. **Try Example Data:**
   - Use provided example files
   - If examples work, issue is likely your data format

2. **Fresh Installation:**
   - Completely remove application
   - Re-download and reinstall
   - Don't copy old settings/databases

3. **Check System Requirements:**
   - Windows 10/11 64-bit
   - 8+ GB RAM
   - 10+ GB free disk space
   - Modern CPU (2015 or newer)

4. **Contact Support:**
   - Provide detailed information
   - Include log files
   - Be patient - we'll help you resolve it!

---

**Last Updated:** December 2024
