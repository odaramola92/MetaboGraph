# ID Annotation Tab - User Guide

## Overview
The ID Annotation tab maps metabolite and lipid names to standard database identifiers using a multi-database approach. It queries online APIs (PubChem, KEGG) and local databases (HMDB, LipidMaps) to provide comprehensive compound annotations with cross-reference IDs.

## Key Purpose
This tab annotates cleaned metabolite/lipid data with:
- **Database IDs**: HMDB, PubChem CID, KEGG, ChEBI, CAS, LipidMaps
- **Chemical Properties**: Formula, molecular weight, SMILES, InChI
- **Compound Classification**: Super class, class, sub class
- **Synonyms**: Alternative names for each compound

## Two Sub-Tabs
1. **🧪 Metabolites**: Annotate metabolomics data
2. **🧬 Lipids**: Annotate lipidomics data

---

## Metabolite ID Annotation

### Input Data
- **Auto-loaded** from Data Cleaning tab (if available)
- OR **Browse** to select a custom Excel file
- Expected sheets: `Positive`, `Negative`, or `Combined`

### Annotation Modes

#### Standard Mode (Default)
- Requires **Positive** and/or **Negative** sheets in input file
- Processes both polarities automatically
- Full annotation with all database cross-references

#### Custom ID Search Mode
- Works with **any sheet** structure
- ID lookup only (fewer annotations)
- Use when file doesn't have standard Pos/Neg sheets
- Requires column verification before running

### Database Sources

**Online APIs (Primary)**:
- **PubChem**: Primary compound lookup
- **KEGG**: Cross-reference retrieval

**Offline Databases (Fallback)**:
- **HMDB**: Human Metabolome Database (local feather file)
- **Synonyms Database**: Extended name matching

### Annotation Output

| Column | Description |
|--------|-------------|
| Name | Original compound name |
| HMDB_ID | Human Metabolome Database ID |
| PubChem_CID | PubChem Compound ID |
| KEGG_ID | KEGG Compound ID |
| ChEBI_ID | ChEBI database ID |
| CAS_Number | CAS Registry Number |
| LipidMaps_ID | LipidMaps ID (if applicable) |
| Formula | Molecular formula |
| Exact_Mass | Exact molecular mass |
| SMILES | Chemical structure notation |
| InChI | International Chemical Identifier |
| InChIKey | Hashed InChI |
| Super_Class | Chemical taxonomy level 1 |
| Class | Chemical taxonomy level 2 |
| Sub_Class | Chemical taxonomy level 3 |
| Synonyms | Alternative compound names |

### Workflow

1. **Select Input File**
   - Auto-loaded from Data Cleaning OR browse manually
   
2. **Choose Annotation Mode**
   - Standard: For files with Pos/Neg sheets
   - Custom: For any file structure
   
3. **Set Output File Path**
   - Default: `metabolite_ids_annotated.xlsx`
   
4. **Configure Workers**
   - Parallel processing threads (default: 2-6)
   - More workers = faster but higher API load
   
5. **Click "🚀 Start ID Annotation"**
   - Monitor progress bar and log
   - ETA displayed during processing
   
6. **Review Results**
   - Open output folder when complete
   - Check annotation coverage statistics

### Re-Filter Feature
Apply new filters to previously annotated data without re-running annotation:
1. Upload existing annotated Excel file
2. Apply new filter criteria
3. Export filtered results

---

## Lipid ID Annotation

### Input Data
- **Auto-loaded** from Data Cleaning tab (lipids)
- OR **Browse** to select lipid Excel file

### Lipid-Specific Features
- **Class Assignment**: Lipid class from LipidMaps taxonomy
- **Adduct Handling**: Processes all detected adducts
- **Chain Composition**: Parses fatty acid chains where available

### Supported Lipid Classes
| Class Code | Full Name |
|------------|-----------|
| FA | Fatty Acids |
| GL | Glycerolipids (DG, TG) |
| GP | Glycerophospholipids (PC, PE, PI, PS, PG, PA) |
| SP | Sphingolipids (Cer, SM, GlcCer) |
| ST | Sterol Lipids (CE, cholesterol) |
| PR | Prenol Lipids |

### Lipid Output
Same columns as metabolite annotation plus:
- **Class_Name**: Lipid class (from LipidSearch)
- **Lipid_Species**: Parsed species name

---

## Processing Settings

### Parallel Workers
- **Range**: 1 to CPU core count
- **Recommended**: 2-6 workers
- **Impact**: More workers = faster processing, but may hit API rate limits

### API Rate Limiting
- Automatic delays between requests
- Retry logic for failed requests
- Timeout handling

---

## Output Files

### Annotated Output
**Filename**: `metabolite_ids_annotated.xlsx` or `lipid_ids_annotated.xlsx`

**Sheets**:
- **Positive**: Annotated positive mode data
- **Negative**: Annotated negative mode data
- **Combined**: Merged positive and negative (if both present)
- **Statistics**: Annotation coverage summary

---

## Tips & Best Practices

### Before Annotation
- ✅ Clean data first using Data Cleaning tab
- ✅ Ensure compound names are accurate
- ✅ Check internet connection for API access
- ✅ Start with fewer workers if experiencing API errors

### During Annotation
- ✅ Monitor progress log for errors
- ✅ Don't interrupt during batch processing
- ✅ Watch for rate limit warnings

### After Annotation
- ✅ Review annotation coverage percentage
- ✅ Check unmatched compounds manually
- ✅ Proceed to Pathway Annotation tab

### Improving Match Rate
- Use standard nomenclature (IUPAC names)
- Remove adduct information from names (e.g., "[M+H]+")
- Check for typos in compound names
- Try synonyms for unmatched compounds

---

## Troubleshooting

### Low annotation rate (<50%)
- **Cause**: Non-standard compound names
- **Solution**: Check naming conventions, use synonyms

### "API rate limit exceeded"
- **Cause**: Too many parallel workers
- **Solution**: Reduce workers to 2, wait and retry

### "File not found" errors
- **Cause**: Missing database files
- **Solution**: Check Database Setup tab, re-download databases

### "No Positive/Negative sheet found"
- **Cause**: Wrong annotation mode selected
- **Solution**: Use Custom ID Search mode for non-standard files

### Slow processing
- **Cause**: Large dataset or slow internet
- **Solution**: Increase workers, check network connection

---

## What's Next?

After ID annotation, proceed to:
1. **Pathway Annotation Tab** → Map metabolites to biological pathways

---

## Support
For questions or issues, refer to the Help tab or submit an issue on GitHub.
- ✅ Review statistics: Expect 60-80% annotation rate
- ✅ Manually curate high-priority/abundant features
- ✅ Cross-validate with MS/MS data if available
- ✅ Document annotation confidence in reports

### Performance Optimization
- Process large datasets in batches (<1000 features)
- Enable caching to avoid redundant queries
- Use offline mode first, online for gaps
- Schedule long annotations during off-hours

## Troubleshooting

### Low Annotation Rate (<40%)
**Possible Causes**:
- Incorrect m/z column (might be neutral mass instead of m/z)
- Too strict tolerance settings
- Species-specific metabolites not in human databases
- Xenobiotics or novel compounds

**Solutions**:
- Verify m/z calculation (consider adducts: [M+H]+, [M-H]-)
- Increase m/z tolerance slightly
- Try alternative databases
- Consider manual annotation for key features

### Slow Performance
**Causes**:
- Large dataset (>2000 features)
- Online API queries with rate limiting
- Network connectivity issues

**Solutions**:
- Annotate in batches
- Disable online APIs for first pass
- Check internet connection
- Run overnight for large datasets

### Multiple Matches for Single Feature
**Causes**:
- Common mass (e.g., 180.06 Da)
- Insufficient specificity
- Isomers with same formula

**Solutions**:
- Use RT information if available
- Consider biological context
- Review all candidates, choose most likely
- Flag for MS/MS confirmation

### Missing Database Files
**Error**: "Database file not found"
**Solution**:
1. Go to Database Setup tab
2. Download required databases (HMDB, LipidMaps)
3. Process into .feather format
4. Verify files in Databases/ folder

## Understanding Output Columns

### Added Annotation Columns
| Column | Description | Example |
|--------|-------------|---------|
| HMDB_ID | HMDB identifier | HMDB0000001 |
| PubChem_CID | PubChem compound ID | 5950 |
| KEGG_ID | KEGG compound ID | C00031 |
| ChEBI_ID | ChEBI identifier | CHEBI:15377 |
| InChI | International Chemical Identifier | InChI=1S/C6H12O6/... |
| InChIKey | Hashed InChI | WQZGKKKJIJFFOK-GASJEMHNSA-N |
| SMILES | Simplified molecular input | C(C1C(C(C(C(O1)O)O)O)O)O |
| Formula | Molecular formula | C6H12O6 |
| MonoMass | Monoisotopic mass | 180.0634 |
| Synonyms | Alternative names | glucose; dextrose; D-glucose |
| SuperClass | Top-level classification | Organic compounds |
| Class | Mid-level classification | Carbohydrates |
| SubClass | Specific classification | Monosaccharides |
| Annotation_Confidence | Confidence score | High, Medium, Low |
| Match_Method | How it was matched | Name+m/z, m/z_only, Formula |
| Database_Sources | Which DBs had hits | HMDB,KEGG,PubChem |

## Next Steps

After ID Annotation:
1. **Quality Check**: Review annotation statistics and confidence distribution
2. **Manual Curation**: Address unmatched or low-confidence features
3. **Pathway Analysis**: Proceed to map metabolites to pathways
4. **Export Data**: Save annotated dataset for external tools (MetaboAnalyst, XCMS, etc.)

## Additional Resources
- [ID Annotation Best Practices](./annotation-best-practices.md)
- [Database Setup Guide](./database-setup.md)
- [Troubleshooting Common Issues](../troubleshooting.md)
- [API Rate Limits and Policies](./api-guidelines.md)

## Support
For annotation questions:
- Review [FAQ](../faq.md)
- Check [Common Annotation Issues](../troubleshooting.md#annotation-issues)
- Open GitHub issue with example data
- Email: support@metabograph.org
