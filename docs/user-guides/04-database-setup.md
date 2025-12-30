# Database Setup Tab - User Guide

## Overview
The Database Setup tab manages the local database files required for MetaboGraph's annotation and pathway analysis features. The tab is divided into **two columns**:
1. **Left Column**: ID Annotation Databases (HMDB, LipidMaps)
2. **Right Column**: Pathway Annotation Databases (SMPDB, WikiPathways, PathBank)

## Database Location
All database files are stored in the **Databases/** folder:
- **For Executable**: `MetaboGraph/Databases/`
- **For Python Script**: Same location relative to `metabograph.py`

---

## Required Database Files

### ID Annotation Databases (Left Column)

| File | Description | Size |
|------|-------------|------|
| `hmdb_database.feather` | HMDB metabolite database | ~150 MB |
| `All_metabolites_synonyms_hmdb.feather` | HMDB synonym lookup | ~50 MB |
| `lipidmap.feather` | LipidMaps lipid database | ~50 MB |

### Pathway Annotation Databases (Right Column)

| File | Description | Size |
|------|-------------|------|
| `merged_SMP_metabolites.feather` | SMPDB pathway-metabolite mapping | ~20 MB |
| `pathbank_selected.feather` | PathBank pathways (Human/Mouse/Rat) | ~30 MB |
| `wikipathways_homo_sapiens.feather` | WikiPathways (Human) | ~20 MB |
| `wikipathways_mus_musculus.feather` | WikiPathways (Mouse) | ~15 MB |
| `wikipathways_rattus_norvegicus.feather` | WikiPathways (Rat) | ~10 MB |

---

## Features

### Database Status Display
Each database shows:
- ✅ **Found**: Database file exists and is valid
- ❌ **Missing**: File not found
- 📊 **Size**: File size on disk
- 📁 **Path**: Location of database file

### Download Links
Click **"Download"** button for each database to open the official download page:

| Database | Download URL |
|----------|--------------|
| HMDB | https://hmdb.ca/downloads |
| LipidMaps | https://www.lipidmaps.org/databases/lmsd/download |
| SMPDB | https://smpdb.ca/downloads |
| PathBank | https://pathbank.org/downloads |
| WikiPathways | https://www.wikipathways.org/ |

### Process Downloaded Files
After downloading raw files, click **"Process"** to convert to `.feather` format:

1. **Browse** to select downloaded file (XML, SDF, or CSV)
2. **Processing** starts automatically
3. **Progress** displayed in log
4. **Output** saved to Databases/ folder
5. **Status** updates when complete

---

## Processing Details

### HMDB Processing
- **Input**: `hmdb_metabolites.xml` (~1.5 GB)
- **Output**: 
  - `hmdb_database.feather` (~150 MB)
  - `All_metabolites_synonyms_hmdb.feather` (~50 MB)
- **Time**: 5-10 minutes
- **Extracts**: HMDB ID, name, formula, mass, SMILES, InChI, synonyms

### LipidMaps Processing
- **Input**: `structures.sdf` (~150 MB)
- **Output**: `lipidmap.feather` (~50 MB)
- **Time**: 2-5 minutes
- **Extracts**: LipidMaps ID, name, class, subclass, SMILES

### SMPDB Processing
- **Input**: SMPDB CSV files
- **Output**: `merged_SMP_metabolites.feather`
- **Time**: 2-3 minutes
- **Extracts**: Pathway names, metabolite mappings

### PathBank Processing
- **Input**: PathBank CSV files
- **Output**: `pathbank_selected.feather`
- **Time**: 3-5 minutes
- **Filters**: Human, Mouse, Rat pathways only

### WikiPathways Processing
- **Input**: GMT files per species
- **Output**: Species-specific `.feather` files
- **Time**: 1-2 minutes per species
- **Extracts**: Pathway-gene/metabolite mappings

---

## Installation Directory Structure

```
MetaboGraph/
├── metabograph.exe (or metabograph.py)
├── Databases/                          ← All databases here
│   ├── hmdb_database.feather
│   ├── All_metabolites_synonyms_hmdb.feather
│   ├── lipidmap.feather
│   ├── merged_SMP_metabolites.feather
│   ├── pathbank_selected.feather
│   ├── wikipathways_homo_sapiens.feather
│   ├── wikipathways_mus_musculus.feather
│   └── wikipathways_rattus_norvegicus.feather
├── gui/
├── main_script/
└── (other files)
```

---

## Workflow

### First-Time Setup
1. **Open** Database Setup tab
2. **Check** status of each database (red = missing)
3. **Download** raw files from official sources
4. **Process** each downloaded file
5. **Verify** all databases show green checkmark
6. **Proceed** to Data Cleaning tab

### Updating Databases
1. **Download** newer version from source
2. **Process** to regenerate feather file
3. **Old version** is overwritten
4. **Restart** application to use new database

---

## Tips & Best Practices

### Download Tips
- ✅ Use stable internet connection for large files
- ✅ HMDB XML file is ~1.5 GB - allow time
- ✅ Save downloaded files before processing
- ✅ Keep backup of raw downloaded files

### Processing Tips
- ✅ Processing is CPU-intensive - be patient
- ✅ Don't interrupt during processing
- ✅ Watch progress log for errors
- ✅ Sufficient disk space required (~500 MB)

### Maintenance
- ✅ Check for database updates periodically
- ✅ HMDB updates quarterly
- ✅ WikiPathways updates monthly
- ✅ Re-process after major database updates

---

## Troubleshooting

### "Database not found"
- **Cause**: Feather file missing from Databases/ folder
- **Solution**: Download and process the database

### "Processing failed"
- **Cause**: Corrupted download or wrong file format
- **Solution**: Re-download raw file, verify file integrity

### "File too large to process"
- **Cause**: Insufficient RAM for XML parsing
- **Solution**: Close other applications, use 64-bit Python

### "Permission denied"
- **Cause**: Databases/ folder is read-only
- **Solution**: Check folder permissions, run as administrator

### Database files not detected after processing
- **Cause**: Output saved to wrong location
- **Solution**: Verify Databases/ folder path, move files manually if needed

---

## Database Sources

| Database | Organization | URL |
|----------|--------------|-----|
| HMDB | Wishart Lab, U of Alberta | https://hmdb.ca |
| LipidMaps | UCSD | https://www.lipidmaps.org |
| SMPDB | Wishart Lab | https://smpdb.ca |
| PathBank | Wishart Lab | https://pathbank.org |
| WikiPathways | Community | https://www.wikipathways.org |

---

## Related Tabs
- **ID Annotation Tab**: Uses HMDB, LipidMaps databases
- **Pathway Annotation Tab**: Uses SMPDB, PathBank, WikiPathways
- **Pathway Network Tab**: Uses all pathway databases

---

## Support
For database issues, check that all files are present in the Databases/ folder and have valid file sizes. Refer to the Help tab or submit an issue on GitHub.

### Database Not Found Error
**Symptoms**: "Database file not found" message when annotating

**Solutions**:
1. Check Databases/ folder exists
2. Verify .feather files are present
3. Ensure correct filenames (case-sensitive)
4. Check file permissions (not locked)
5. Refresh database status

### Processing Failed
**Common Causes**:
- Corrupted download (incomplete file)
- Insufficient disk space
- Wrong file format selected
- Outdated source file format

**Solutions**:
1. Re-download source file
2. Check available disk space (need 2-3 GB)
3. Verify file integrity (check file size)
4. Use latest database version
5. Check processing logs for specific error

### Processing Takes Too Long
**Normal Times**:
- HMDB: 5-10 minutes (large file)
- Others: 2-5 minutes each

**If Exceeds Normal**:
- Check CPU usage (should be high during processing)
- Ensure no antivirus interference
- Close other heavy applications
- Consider processing one at a time

### Invalid Format Error
**Causes**:
- Downloaded HTML instead of actual file
- Wrong file selected
- Corrupted file

**Solutions**:
1. Re-download from official source
2. Verify file extension (.xml, .sdf, .csv)
3. Open file in text editor to check content
4. Download from direct link (not through browser redirect)

### Out of Memory Error
**Causes**:
- Large database file
- Insufficient RAM
- Multiple processing tasks

**Solutions**:
1. Close other applications
2. Process one database at a time
3. Restart computer to free memory
4. Consider 64-bit Python if using 32-bit

## Best Practices

### Database Management
- ✅ Download all databases during initial setup
- ✅ Process databases immediately after download
- ✅ Keep source files as backup
- ✅ Check for updates quarterly
- ✅ Document database versions used in publications

### Storage
- Keep source files separate from processed files
- Databases/ folder for .feather files only
- Create backup of processed databases
- Version control: name files with date if updating

### Updates
- Subscribe to database update notifications
- Check before starting new projects
- Update all databases together (consistency)
- Test after updating before critical analysis

### Documentation
- Record database versions in analysis logs
- Note processing date
- Keep changelog of updates
- Include in methods section of papers

## Database Information

### HMDB
- **Source**: Human Metabolome Database
- **URL**: https://hmdb.ca
- **Update Frequency**: Quarterly
- **Current Version**: 5.0 (as of 2024)
- **License**: Free for academic use

### LipidMaps
- **Source**: LIPID MAPS Consortium
- **URL**: https://www.lipidmaps.org
- **Update Frequency**: Annually
- **Current Version**: 2024
- **License**: Free for academic use

### PathBank
- **Source**: PathBank Database
- **URL**: https://pathbank.org
- **Update Frequency**: Bi-annually
- **Coverage**: 110,000+ pathways
- **License**: Open access

### SMPDB
- **Source**: Small Molecule Pathway Database
- **URL**: https://smpdb.ca
- **Update Frequency**: Annually
- **Coverage**: 30,000+ pathway associations
- **License**: Free for academic use

### WikiPathways
- **Source**: WikiPathways Community
- **URL**: https://www.wikipathways.org
- **Update Frequency**: Monthly (community-driven)
- **Coverage**: 2,800+ pathways
- **License**: CC BY 3.0

## FAQs

**Q: Do I need all databases?**
A: HMDB is essential for ID annotation. Others enhance pathway analysis. Start with HMDB, add others as needed.

**Q: How much disk space needed?**
A: ~500 MB for all processed databases, plus ~2 GB for source files during processing.

**Q: Can I use older database versions?**
A: Yes, but newer versions have more metabolites and better annotations. Update for best results.

**Q: Processing failed, what now?**
A: Check the error log, verify source file integrity, ensure sufficient disk space, try reprocessing.

**Q: How often should I update?**
A: Check quarterly. Update before major projects or publications.

**Q: Can I use custom databases?**
A: Advanced users can add custom databases. Contact support for format specifications.

## Next Steps

After Database Setup:
1. ✅ Verify all databases show as found
2. ✅ Test with sample data (ID annotation)
3. ✅ Proceed to ID Annotation tab
4. ✅ Document database versions for your project

## Additional Resources
- [Database Format Specifications](./database-formats.md)
- [Custom Database Integration](./custom-databases.md)
- [Database Update Guide](./database-updates.md)
- [Troubleshooting Common Issues](../troubleshooting.md)

## Support
For database setup questions:
- [FAQ - Database Setup](../faq.md#database-setup)
- [Video Tutorial - Database Setup](link-to-video)
- GitHub Issues: [github-link]
- Email: support@metabograph.org
