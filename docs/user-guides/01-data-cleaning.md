# Data Cleaning Tab - User Guide

## Overview
The Data Cleaning tab is designed to clean and standardize raw mass spectrometry output files from **metabolomics** and **lipidomics** experiments. It handles data from instruments like Thermo Fisher Compound Discoverer and LipidSearch, processing positive and negative ionization modes separately or combined.

## Key Purpose
This tab cleans raw MS data by:
- **Filtering by Reference Ions/Adducts**: Select which ion adducts to keep
- **Removing Duplicates**: Uses PPM threshold and RT deviation to identify duplicates
- **Cleaning Column Names**: Removes prefixes/suffixes from sample column names
- **Merging MZCloud Data**: Integrates additional compound information (optional)
- **Splitting by Polarity**: Handles combined or separate positive/negative mode files

## Features

### Two Sub-Tabs
1. **🧪 Metabolites**: Clean metabolomics data
2. **🧬 Lipids**: Clean lipidomics data

---

## Metabolite Cleaning

### Processing Modes
1. **Combined File Mode**: Single Excel file containing both positive and negative mode data (split by polarity internally)
2. **Separate Files Mode**: Individual positive and negative mode Excel files

### File Inputs
**Required**:
- Metabolite Excel file(s) from Compound Discoverer

**Optional**:
- MZCloud Excel file(s) for additional compound annotations

### Basic Settings Tab

#### File Selection
1. Choose **Processing Mode**: Combined or Separate
2. Browse and select your input file(s)
3. *(Optional)* Add MZCloud file(s) for enriched annotations
4. Set **Output File Path** for cleaned data

### Advanced Settings Tab

#### 📊 Filtering Thresholds
| Parameter | Default | Description |
|-----------|---------|-------------|
| **PPM Threshold** | 10.0 | Maximum allowed ppm difference for duplicate detection |
| **RT Threshold** | 4.0 | Retention time deviation tolerance for deduplication |
| **MS2 Filter** | None | Optional MS2 filtering (DDA preferred ion, non-preferred, or none) |

#### ⚡ Reference Ion Filtering
Select which reference ions to include in analysis:

**Positive Mode Ions**:
- `[M+H]+1` - Protonated molecule
- `[M+2H]+2` - Doubly protonated
- `[M+H-H2O]+1` - Protonated with water loss
- `[M+H+MeOH]+1` - Methanol adduct
- `[M+FA+H]+1` - Formic acid adduct

**Negative Mode Ions**:
- `[M-H]-1` - Deprotonated molecule
- `[M-2H]-2` - Doubly deprotonated
- `[M-H-H2O]-1` - Deprotonated with water loss
- `[M-H-MeOH]-1` - Methanol adduct loss
- `[M+FA-H]-1` - Formate adduct

**Controls**:
- ✓ Select All
- ✗ Deselect All
- ↻ Reset to Default (all checked)

#### 🧹 Column Name Cleaning
Remove specified text patterns from sample column names.

**Default Patterns**:
```
Group Area:, _pos, _neg, .raw, Area[, ], [, ], pos, pos_, neg, neg_, positive, negative
```

**Examples**:
- `Group Area: Sample1` → `Sample1`
- `Control_pos` → `Control`
- `Treated.raw` → `Treated`

### Workflow
1. Select **Processing Mode** (Combined or Separate)
2. **Browse** and select input file(s)
3. *(Optional)* Add MZCloud file(s)
4. Set **Output Path**
5. Click **🔍 Verify Columns** to validate file structure
6. *(Optional)* Adjust **Advanced Settings** (thresholds, ions, patterns)
7. Click **🚀 Start Data Cleaning**
8. Monitor progress in the **Progress Log**
9. Click **📂 Open Output Folder** when complete

---

## Lipid Cleaning

### Processing Modes
1. **Combined File Mode**: Single LipidSearch export with both polarities
2. **Separate Files Mode**: Individual positive and negative mode files

### File Inputs
**Required**:
- Lipid Excel file(s) from LipidSearch

**Optional**:
- File ID Mapping file(s) for sample name cleanup

### Adduct Ion Filtering
Select which lipid adducts to include:

**Positive Mode Adducts**:
- M+H, M+2H, M+NH4, M+H-H2O
- M+HCOO, M+CH3COO
- M+Na, M+K, M+Li

**Negative Mode Adducts**:
- M-H, M-2H, M-CH3

### Column Name Cleaning
**Default Patterns**:
```
Area[, ], [, ], _pos, _neg, .raw, _Area, Area_
```

### Workflow
1. Switch to **🧬 Lipids** sub-tab
2. Select **Processing Mode**
3. **Browse** and select lipid input file(s)
4. *(Optional)* Add File ID mapping file(s)
5. Set **Output Path**
6. Click **🔍 Verify Columns**
7. Adjust **Adduct Filtering** as needed
8. Click **🚀 Start Lipid Cleaning**
9. Monitor progress and open output folder when complete

---

## Output Files

### Metabolite Cleaned Output
**Filename**: `cleaned_metabolites.xlsx` (or user-specified)

**Sheets**:
- **Positive**: Cleaned positive mode data
- **Negative**: Cleaned negative mode data  
- **Combined**: Merged positive and negative data

### Lipid Cleaned Output
**Filename**: `cleaned_lipids.xlsx` (or user-specified)

**Sheets**:
- **Positive**: Cleaned positive mode lipids
- **Negative**: Cleaned negative mode lipids
- **Class_Table**: Lipid class summary

---

## Tips & Best Practices

### Before Cleaning
- ✅ Export data from Compound Discoverer / LipidSearch in Excel format
- ✅ Keep original raw files as backup
- ✅ Use consistent file naming conventions
- ✅ Verify sample names are present in column headers

### Ion Selection
- ✅ Use default (all ions) for comprehensive analysis
- ✅ Deselect specific ions only if you have reason to exclude them
- ✅ Consider your ionization source and mobile phases

### Thresholds
- ✅ PPM 10.0 is appropriate for most high-resolution instruments
- ✅ Adjust RT threshold based on your chromatography reproducibility
- ✅ Tighter thresholds = fewer duplicates removed, looser = more aggressive deduplication

### Column Cleaning
- ✅ Add instrument-specific prefixes to removal patterns
- ✅ Ensure cleaned column names will match your experimental metadata
- ✅ Review changes in output to verify correct pattern removal

---

## Troubleshooting

### "Verify Columns" shows errors
- **Cause**: Required columns not found in input file
- **Solution**: Check that file is proper Compound Discoverer/LipidSearch export

### No data after cleaning
- **Cause**: Too restrictive ion filtering or thresholds
- **Solution**: Select more reference ions, loosen PPM/RT thresholds

### Column names not cleaned properly
- **Cause**: Patterns don't match your file's naming convention
- **Solution**: Add your specific prefixes/suffixes to pattern list

### Processing takes too long
- **Cause**: Very large dataset (>50,000 features)
- **Solution**: Process positive and negative modes separately

---

## What's Next?

After data cleaning, proceed to:
1. **ID Annotation Tab** → Match metabolites/lipids to database identifiers
2. **Pathway Annotation Tab** → Map to biological pathways

---

## Support
For questions or issues, refer to the Help tab or submit an issue on GitHub.
- Check the [FAQ](../faq.md)
- Open an issue on [GitHub](github-link)
- Contact: your-email@example.com
