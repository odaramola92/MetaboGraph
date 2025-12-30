# Multi-Omics Integration Tab - User Guide

## Overview
The Multi-Omics Integration tab is designed to **merge multiple pathway-annotated datasets** into a single unified file for downstream network analysis. This is specifically for combining related datasets (e.g., metabolite data + lipid data from the same samples) before pathway network visualization.

## Key Purpose
This tab does **NOT**:
- Perform correlation analysis between omics layers
- Do transcriptomics/proteomics integration
- Run statistical comparisons

This tab **DOES**:
- Merge multiple pathway-annotated Excel files into one
- Validate that each file has the required columns
- Consolidate pathway annotations from all sources
- Push the merged dataset directly to the Pathway Network tab

## When to Use

Use this tab when you have:
- **Metabolite + Lipid data** from the same experiment
- **Multiple annotated datasets** that need to be analyzed together
- **Separate positive/negative mode** results to combine
- Any scenario where you need a **unified dataset for network analysis**

## Features

### 1. Upload Multiple Datasets
- Upload **unlimited files** for merging
- Each file should be a pathway-annotated dataset (from Pathway Annotation tab)
- Name each dataset for tracking

### 2. Required Columns
Each uploaded file should contain:
- **Feature ID**: Metabolite/compound identifier
- **P-Value**: Statistical significance (adjusted p-value preferred)
- **Log2 Fold Change**: Direction and magnitude of change
- **Pathway columns** (optional but recommended): HMDB_Pathways, KEGG_Pathways, etc.

### 3. Column Validation
- Auto-detects required columns
- Shows confirmation dialog for each file
- Allows manual column assignment if auto-detection fails
- Supports adjusted p-values (adj.P, FDR, q-value) with priority

### 4. Smart Merging
When merging multiple datasets:
- **P-values**: Keeps the most significant (minimum) value
- **Log2FC**: Keeps the value with largest absolute magnitude
- **Pathways**: Concatenates all pathway annotations, removing duplicates
- **Source**: Tracks which dataset each feature came from

## Workflow

### Step 1: Upload Datasets
1. Click **"Add Files"**
2. Select multiple pathway-annotated Excel files
3. Give each dataset a descriptive name
4. Files appear in the list showing rows and columns

### Step 2: Process & Validate
1. Click **"Process & Validate Data"**
2. System auto-detects columns for each file
3. **Column Assignment Dialog** appears for each file:
   - Verify Feature ID column
   - Verify P-Value column (system prefers adjusted p-values)
   - Verify Log2FC column
4. Confirm or adjust column assignments

### Step 3: Merge & Export
1. Once all files validated (green checkmarks)
2. Click **"Merge & Export to Network Tab"**
3. System merges all datasets:
   - Combines on Feature_ID
   - Consolidates statistics
   - Merges pathway columns
   - Preserves upstream regulators and disease associations
4. Creates unified Excel file
5. Automatically pushes to Pathway Network tab

## Example Use Cases

### Example 1: Metabolites + Lipids
```
Uploads:
  - pathway_annotated_metabolites.xlsx (500 metabolites)
  - pathway_annotated_lipids.xlsx (300 lipids)

Result:
  - merged_multiomics.xlsx (800 features)
  - Both metabolite and lipid pathways included
  - Ready for unified network analysis
```

### Example 2: Multiple Experiments
```
Uploads:
  - treatment_A_annotated.xlsx
  - treatment_B_annotated.xlsx
  - treatment_C_annotated.xlsx

Result:
  - Merged file with all unique features
  - Source column tracks origin
  - Network shows combined metabolic changes
```

### Example 3: Positive + Negative Mode MS
```
Uploads:
  - positive_mode_annotated.xlsx
  - negative_mode_annotated.xlsx

Result:
  - Unified dataset covering both ionization modes
  - No duplicate features (consolidated by Feature_ID)
```

## Log Output

```
📂 Loading 2 file(s)...
✅ Loaded Metabolites (metabolites.xlsx): 500 rows, 45 columns
✅ Loaded Lipids (lipids.xlsx): 300 rows, 42 columns

🔄 Starting data validation...

📊 Validating Metabolites...
  ✓ Feature ID: Name
  ✓ P-Value: adj.P (adjusted)
  ✓ Log2FC: log2FoldChange
  ✓ Found 3 pathway column(s)

📊 Validating Lipids...
  ✓ Feature ID: Lipid_Name
  ✓ P-Value: FDR_pvalue (adjusted)
  ✓ Log2FC: Log2_FC
  ✓ Found 2 pathway column(s)

✅ All 2 file(s) validated successfully!

🔄 Merging datasets...
  ✓ Prepared Metabolites: 500 rows, 3 pathway cols
  ✓ Prepared Lipids: 300 rows, 2 pathway cols (Class_name ✓)

🔧 Consolidating columns...
✅ Merge complete: 800 total features
✅ Exported to: merged_multiomics.xlsx
✅ Pushed to Network Tab
```

## Merged Output Structure

The merged file contains:

| Column | Description |
|--------|-------------|
| Feature_ID | Unique identifier for each feature |
| pvalue | Most significant p-value across sources |
| log2FC | Fold change with largest magnitude |
| Source | Origin dataset(s) |
| HMDB_Pathways | Combined HMDB pathway annotations |
| KEGG_Pathways | Combined KEGG pathway annotations |
| SMPDB_Pathways | Combined SMPDB pathway annotations |
| Upstream_Enzymes | Combined enzyme regulators |
| Disease_Associations | Combined disease links |
| Class_name | Lipid class (if from lipid data) |

## Tips & Best Practices

### Before Merging
- ✅ Complete pathway annotation for all datasets first
- ✅ Ensure datasets are from same experiment/samples
- ✅ Use consistent statistical thresholds
- ✅ Verify Feature_ID format is compatible

### Column Naming
- ✅ Use standard column names (Name, pvalue, log2FC)
- ✅ Include "pathway" in pathway column names
- ✅ Use adj.P or FDR for adjusted p-values (auto-prioritized)

### After Merging
- ✅ Review merged file statistics
- ✅ Proceed to Network tab for visualization
- ✅ Source column helps track data origin

## Troubleshooting

### "Multi-omics merging requires at least 2 datasets"
- Upload at least 2 files before processing

### Column auto-detection fails
- Manually select correct columns in dialog
- Check column names match expected patterns

### Duplicate features after merge
- Features are matched by Feature_ID
- Different Feature_IDs = separate rows
- Same Feature_ID = values consolidated

### Missing pathway columns
- Pathway columns are optional
- Merging still works but network analysis limited
- Ensure pathway annotation was completed before exporting

## Next Steps

After Multi-Omics merge:
1. **Review merged data** in the export Excel file
2. **Go to Pathway Network tab** (auto-loaded)
3. **Run network analysis** on unified dataset
4. **Generate visualizations** with combined metabolite+lipid pathways

## Related Tabs
- **Pathway Annotation Tab**: Annotate datasets before merging
- **Pathway Network Tab**: Visualize merged data
- **Comparative Analysis Tab**: Compare pathway results (different purpose)
